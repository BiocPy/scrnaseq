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
def list_datasets(cache_dir: str = cache_directory(), overwrite: bool = False, latest: bool = True) -> pd.DataFrame:
    """List all available datasets.

    Example:

        .. code-block:: python

            datasets = (
                list_datasets()
            )

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
        A :py:class:`~pandas.DataFrame` where each row corresponds to a dataset.
        Each row contains title and description for each dataset,
        the number of rows and columns, the organisms and genome builds involved,
        whether the dataset has any pre-computed reduced dimensions, and so on.
        More details can be found in the
        `Bioconductor metadata schema <https://github.com/ArtifactDB/bioconductor-metadata-index>`_.
    """
    db_path = fetch_metadata_database(cache_dir=cache_dir, overwrite=overwrite)
    conn = sqlite3.connect(db_path, check_same_thread=False)

    stmt = "SELECT json_extract(metadata, '$') AS meta, versions.asset AS asset, versions.version AS version, path"
    key_names = ["meta", "asset", "version", "path"]
    if latest is not True:
        stmt = f"{stmt} versions.latest AS latest"
        key_names.append("latest")

    stmt = f"{stmt} FROM paths LEFT JOIN versions ON paths.vid = versions.vid WHERE versions.project = 'scRNAseq'"
    if latest is True:
        stmt = f"{stmt} AND versions.latest = 1"

    _qresults = conn.execute(stmt).fetchall()
    conn.close()

    results = _format_query_results(_qresults, key_names)

    return _sanitize_query_to_output(results, latest)


def _format_query_results(results: list, key_names: list):
    """Format the results from sqlite as a pandas dataframe.

    Key names must be in the exact same order as the query.
    """
    _out = {}
    for k in key_names:
        _out[k] = []

    for r in results:
        for idx, k in enumerate(key_names):
            _out[k].append(r[idx])

    return _out


def _sanitize_query_to_output(results: list, latest: bool, meta_name: str = "meta"):
    _all_paths = [None if "/" not in p else p.rsplit("/", 1)[0] for p in results["path"]]

    df = pd.DataFrame(
        {
            "name": results["asset"],
            "version": results["version"],
            "path": _all_paths,
        }
    )
    if not latest:
        _all_latest = [s == 1 for s in results["latest"]]
        df["latest"] = _all_latest

    _all_metas = [json.loads(s) for s in results[meta_name]]

    df["object"] = _extract_atomic_from_json(
        _all_metas, lambda x: x.get("applications", {}).get("takane", {}).get("type")
    )
    df["title"] = _extract_atomic_from_json(_all_metas, lambda x: x.get("title"))
    df["description"] = _extract_atomic_from_json(_all_metas, lambda x: x.get("title"))
    df["taxonomy_id"] = _extract_charlist_from_json(_all_metas, lambda x: x.get("taxonomy_id"))
    df["genome"] = _extract_charlist_from_json(_all_metas, lambda x: x.get("genome"))

    df["rows"] = _extract_atomic_from_json(
        _all_metas,
        lambda x: x.get("applications", {}).get("takane", {}).get("summarized_experiment", {}).get("rows"),
    )

    df["columns"] = _extract_atomic_from_json(
        _all_metas,
        lambda x: x.get("applications", {}).get("takane", {}).get("summarized_experiment", {}).get("columns"),
    )

    df["assays"] = _extract_charlist_from_json(
        _all_metas,
        lambda x: x.get("applications", {}).get("takane", {}).get("summarized_experiment", {}).get("assays"),
    )
    df["column_annotations"] = _extract_charlist_from_json(
        _all_metas,
        lambda x: x.get("applications", {})
        .get("takane", {})
        .get("summarized_experiment", {})
        .get("column_annotations"),
    )
    df["reduced_dimensions"] = _extract_charlist_from_json(
        _all_metas,
        lambda x: x.get("applications", {})
        .get("takane", {})
        .get("single_cell_experiment", {})
        .get("reduced_dimensions"),
    )
    df["alternative_experiments"] = _extract_charlist_from_json(
        _all_metas,
        lambda x: x.get("applications", {})
        .get("takane", {})
        .get("single_cell_experiment", {})
        .get("alternative_experiments"),
    )

    df["bioconductor_version"] = _extract_atomic_from_json(_all_metas, lambda x: x.get("bioconductor_version"))
    df["maintainer_name"] = _extract_atomic_from_json(_all_metas, lambda x: x.get("maintainer_name"))
    df["maintainer_email"] = _extract_atomic_from_json(_all_metas, lambda x: x.get("maintainer_email"))

    sources = []
    for meta in _all_metas:
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
    return [extract(_meta) if extract(_meta) is not None else None for _meta in metadata]


def _extract_charlist_from_json(metadata, extract):
    return [extract(_meta) if extract(_meta) is not None else [] for _meta in metadata]
