import json
import sqlite3
from functools import lru_cache

import pandas as pd
from gypsum_client import (
    cache_directory,
    fetch_metadata_database,
)

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


@lru_cache
def list_datasets(
    cache_dir: str = cache_directory(), overwrite: bool = False, latest: bool = True
) -> pd.DataFrame:
    """List all available datasets.

    Example:

        .. code-block:: python

            datasets = list_datasets()

    Args:
        cache_dir:
            Path to cache directory.

        overwrite:
            Whether to overwrite the database in cache.
            Defaults to False.

        latest:
            Whether to only fetch the latest version of each dataset.
            Defaults to True.

    Returns:
        A pandas DataFrame where each row corresponds to a dataset.
        Each row contains title and description for each dataset,
        the number of rows and columns, the organisms and genome builds involved,
        whether the dataset has any pre-computed reduced dimensions, and so on.
        More details can be found in the
        `Bioconductor metadata schema <https://github.com/ArtifactDB/bioconductor-metadata-index>`_.
    """
    db_path = fetch_metadata_database(cache_dir=cache_dir, overwrite=overwrite)
    conn = sqlite3.connect(db_path, check_same_thread=False)

    stmt = "SELECT json_extract(metadata, '$') AS meta, versions.asset AS asset, versions.version AS version, path"
    if latest is not True:
        stmt = f"{stmt} versions.latest AS latest"

    stmt = f"{stmt} FROM paths LEFT JOIN versions ON paths.vid = versions.vid WHERE versions.project = 'scRNAseq'"
    if latest is True:
        stmt = f"{stmt} AND versions.latest = 1"

    results = pd.read_sql_query(stmt, conn)
    conn.close()

    return _sanitize_query_to_output(results, latest)


def _sanitize_query_to_output(results, latest, meta_name="meta"):
    results["path"] = results["path"].apply(
        lambda p: None if "/" not in p else p.rsplit("/", 1)[0]
    )
    df = pd.DataFrame(
        {
            "name": results["asset"],
            "version": results["version"],
            "path": results["path"],
        }
    )
    if not latest:
        df["latest"] = results["latest"] == 1

    all_meta = results[meta_name].apply(json.loads)

    df["object"] = _extract_atomic_from_json(
        all_meta, lambda x: x.get("applications", {}).get("takane", {}).get("type")
    )
    df["title"] = _extract_atomic_from_json(all_meta, lambda x: x.get("title"))
    df["description"] = _extract_atomic_from_json(all_meta, lambda x: x.get("title"))
    df["taxonomy_id"] = _extract_charlist_from_json(
        all_meta, lambda x: x.get("taxonomy_id")
    )
    df["genome"] = _extract_charlist_from_json(all_meta, lambda x: x.get("genome"))

    df["rows"] = _extract_atomic_from_json(
        all_meta,
        lambda x: x.get("applications", {})
        .get("takane", {})
        .get("summarized_experiment", {})
        .get("rows"),
    )

    df["columns"] = _extract_atomic_from_json(
        all_meta,
        lambda x: x.get("applications", {})
        .get("takane", {})
        .get("summarized_experiment", {})
        .get("columns"),
    )

    df["assays"] = _extract_charlist_from_json(
        all_meta,
        lambda x: x.get("applications", {})
        .get("takane", {})
        .get("summarized_experiment", {})
        .get("assays"),
    )
    df["column_annotations"] = _extract_charlist_from_json(
        all_meta,
        lambda x: x.get("applications", {})
        .get("takane", {})
        .get("summarized_experiment", {})
        .get("column_annotations"),
    )
    df["reduced_dimensions"] = _extract_charlist_from_json(
        all_meta,
        lambda x: x.get("applications", {})
        .get("takane", {})
        .get("single_cell_experiment", {})
        .get("reduced_dimensions"),
    )
    df["alternative_experiments"] = _extract_charlist_from_json(
        all_meta,
        lambda x: x.get("applications", {})
        .get("takane", {})
        .get("single_cell_experiment", {})
        .get("alternative_experiments"),
    )

    df["bioconductor_version"] = _extract_atomic_from_json(
        all_meta, lambda x: x.get("bioconductor_version")
    )
    df["maintainer_name"] = _extract_atomic_from_json(
        all_meta, lambda x: x.get("maintainer_name")
    )
    df["maintainer_email"] = _extract_atomic_from_json(
        all_meta, lambda x: x.get("maintainer_email")
    )

    sources = []
    for meta in all_meta:
        cursources = meta.get("sources")
        if cursources is None:
            sources.append(pd.DataFrame(columns=["provider", "id", "version"]))
        else:
            sources.append(
                pd.DataFrame(
                    {
                        "provider": [s.get("provider") for s in cursources],
                        "id": [s.get("id") for s in cursources],
                        "version": [s.get("version") for s in cursources],
                    }
                )
            )
    df["sources"] = sources

    return df


def _extract_atomic_from_json(metadata, extract):
    return [
        extract(_meta) if extract(_meta) is not None else None for _meta in metadata
    ]


def _extract_charlist_from_json(metadata, extract):
    return [extract(_meta) if extract(_meta) is not None else [] for _meta in metadata]
