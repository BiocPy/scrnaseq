import sys

if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional) when `python_requires = >= 3.8`
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

from .fetch_dataset import fetch_dataset, fetch_metadata
from .list_datasets import list_datasets
from .list_versions import fetch_latest_version, list_versions
from .polish_dataset import polish_dataset
from .save_dataset import save_dataset
from .search_datasets import search_datasets
from .upload_dataset import upload_dataset
