import sqlite3
from functools import lru_cache
from typing import Union

import pandas as pd
from gypsum_client import cache_directory, fetch_metadata_database
from gypsum_client.search_metadata import (
    GypsumSearchClause,
    search_metadata_text_filter,
)

from .list_datasets import _format_query_results, _sanitize_query_to_output

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


@lru_cache
def search_datasets(
    query: Union[str, GypsumSearchClause],
    cache_dir: str = cache_directory(),
    overwrite: bool = False,
    latest: bool = True,
) -> pd.DataFrame:
    """Search for datasets of interest based on matching text in the associated metadata.

    This is a wrapper around
    :py:func:`~gypsum_client.search_metadata.search_metadata_text`.

    The returned DataFrame contains the usual suspects like the title
    and description for each dataset, the number of rows and columns,
    the organisms and genome builds involved, whether the dataset has
    any pre-computed reduced dimensions, and so on.

    More details can be found in the Bioconductor
    `metadata index <https://github.com/ArtifactDB/bioconductor-metadata-index>`_.

    See Also:
        :py:func:`~scrnaseq.list_datasets.list_datasets`, to list all
        available datasets.

        :py:func:`~gypsum_client.search_metadata.search_metadata_text`,
        to search metadata.

    Examples:

    .. code-block:: python

        res = search_datasets("brain")

        res = search_datasets(define_text_query("Neuro%", partial=True")

        res = search_datasets(define_text_query("10090", field="taxonomy_id")

        res = search_datasets(
            define_text_query("GRCm38", field="genome") &
            (define_text_query("neuro%", partial=True) |
                define_text_query("pancrea%", partial=True))
        )

    Args:
        query:
            The search query string or a gypsum.search.object for
            more complex queries.

        cache_directory:
            Path to cache directory.

        overwrite:
            Whether to overwrite the existing cache.
            Defaults to False.

        latest:
            Whether to fetch only the latest versions of datasets.
            Defaults to True.

    Returns:
        A :py:class:`~pandas.DataFrame` where each row corresponds to
        a dataset, containing various columns of metadata.
        Some columns may be lists to capture 1:many mappings.
    """

    bpath = fetch_metadata_database(cache_dir=cache_dir, overwrite=overwrite)

    where = search_metadata_text_filter(query)
    cond = where["where"]
    params = where["parameters"]

    conn = sqlite3.connect(bpath, check_same_thread=False)
    stmt = "SELECT json_extract(metadata, '$') AS meta, versions.asset AS asset, versions.version AS version, path"
    key_names = ["meta", "asset", "version", "path"]

    if not latest:
        stmt += ", versions.latest AS latest"
        key_names.append("latest")

    stmt += " FROM paths LEFT JOIN versions ON paths.vid = versions.vid WHERE versions.project = 'scRNAseq'"

    if latest:
        stmt += " AND versions.latest = 1"

    if cond:
        stmt += " AND " + " AND ".join(cond)
        cursor = conn.execute(stmt, params)
    else:
        cursor = conn.execute(stmt)

    _qresults = cursor.fetchall()
    conn.close()

    results = _format_query_results(_qresults, key_names)
    return _sanitize_query_to_output(results, latest)
