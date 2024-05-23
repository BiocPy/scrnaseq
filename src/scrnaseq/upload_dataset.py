import json
import os
import shutil
from functools import singledispatch
from typing import Any

import dolomite_base as dl
from gypsum_client import cache_directory, rest_url, upload_directory
from singlecellexperiment import SingleCellExperiment

from .utils import format_object_metadata

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def upload_dataset(
    name: str,
    version: str,
    directory: str = None,
    package: str = "scRNAseq",
    cache_dir: str = cache_directory(),
    deduplicate: bool = True,
    probation: bool = False,
    url: str = rest_url(),
    token: str = None,
    concurrent: int = 1,
    abort_failed: bool = True,
):
    return upload_directory(
        directory,
        package,
        name,
        version,
        cache_dir=cache_dir,
        deduplicate=deduplicate,
        probation=probation,
        url=url,
        token=token,
        concurrent=concurrent,
        abort_failed=abort_failed
    )
