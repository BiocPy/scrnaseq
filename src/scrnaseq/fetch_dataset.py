import atexit
import json
import os

from dolomite_base import alt_read_object, alt_read_object_function
from gypsum_client import cache_directory, save_file, save_version
from summarizedexperiment import SummarizedExperiment

from .utils import single_cell_load_object

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def fetch_dataset(
    name: str,
    version: str,
    path: str = None,
    package: str = "scRNAseq",
    cache_dir: str = cache_directory(),
    overwrite: bool = False,
    realize_assays: bool = False,
    realize_reduced_dims: bool = True,
    **kwargs,
) -> SummarizedExperiment:
    """Fetch a single-cell dataset from the gypsum backend.

    See Also:
        `metadata index <https://github.com/ArtifactDB/bioconductor-metadata-index>`_,
        on the expected schema for the metadata.

        :py:func:`~scrnaseq.save_dataset.save_dataset` and
        :py:func:`~gypsum_client.upload_file_operations.upload_directory`,
        to save and upload a dataset.

        :py:func:`~scrnaseq.list_datasets.list_datasets` and :py:func:`~scrnaseq.list_versions.list_versions`,
        to get possible values for `name` and `version`.

    Example:

        .. code-block:: python

            sce = fetch_dataset(
                "zeisel-brain-2015",
                "2023-12-14",
            )

    Args:
        name:
            Name of the dataset.

        version:
            Version of the dataset.

        path:
            Path to a subdataset, if name contains multiple datasets.
            Defaults to None.

        package:
            Name of the package.
            Defaults to "scRNAseq".

        cache_dir:
            Path to cache directory.

        overwrite:
            Whether to overwrite existing files.
            Defaults to False.

        realize_assays:
            Whether to realize assays into memory.
            Defaults to False.

        realize_reduced_dims:
            Whether to realize reduced dimensions into memory.
            Defaults to True.

        **kwargs:
            Further arguments to pass to
            :py:func:`~dolomite_base.read_object.read_object`.

    Returns:
        The dataset as a
        :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`
        or one of its subclasses.
    """

    version_path = save_version(package, name, version, cache_dir=cache_dir, overwrite=overwrite)
    obj_path = version_path if path is None else os.path.join(version_path, path.rstrip("/"))

    old = alt_read_object_function(single_cell_load_object)

    def reset_alt_read_func():
        alt_read_object_function(old)

    atexit.register(reset_alt_read_func)
    return alt_read_object(
        obj_path,
        scrnaseq_realize_assays=realize_assays,
        scrnaseq_realize_reduced_dims=realize_reduced_dims,
        **kwargs,
    )


def fetch_metadata(
    name: str,
    version: str,
    path: str = None,
    package: str = "scRNAseq",
    cache_dir: str = cache_directory(),
    overwrite: bool = False,
):
    """Fetch metadata for a dataset from the gypsum backend.

    See Also:
        :py:func:`~.fetch_dataset`,
        to fetch a dataset.

    Example:

    .. code-block:: python

        meta = fetch_metadata(
            "zeisel-brain-2015",
            "2023-12-14",
        )

    Args:
        name:
            Name of the dataset.

        version:
            Version of the dataset.

        path:
            Path to a subdataset, if name contains multiple datasets.
            Defaults to None.

        package:
            Name of the package.
            Defaults to "scRNAseq".

        cache:
            Path to the cache directory.

        overwrite:
            Whether to overwrite existing files.
            Defaults to False.

    Returns:
        Dictionary containing metadata for the specified dataset.
    """
    remote_path = "_bioconductor.json" if path is None else f"{path}/_bioconductor.json"
    local_path = save_file(package, name, version, remote_path, cache_dir=cache_dir, overwrite=overwrite)

    with open(local_path, "r") as f:
        metadata = json.load(f)

    return metadata
