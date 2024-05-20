import atexit
import json
import os

from delayedarray import is_sparse, to_dense_array, to_scipy_sparse_matrix
from dolomite_base import alt_read_object, alt_read_object_function, read_object
from gypsum_client import cache_directory, save_file, save_version
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment

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
    """Fetch a dataset from the gypsum backend.

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

    cache_dir = cache_directory(cache_dir)

    version_path = save_version(
        package, name, version, cache_dir=cache_dir, overwrite=overwrite
    )
    obj_path = (
        version_path if path is None else os.path.join(version_path, path.rstrip("/"))
    )

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
    local_path = save_file(
        package, name, version, remote_path, cache_dir=cache_dir, overwrite=overwrite
    )

    with open(local_path, "r") as f:
        metadata = json.load(f)

    return metadata


def single_cell_load_object(
    path: str,
    metadata: dict = None,
    scrnaseq_realize_assays: bool = False,
    scrnaseq_realize_reduced_dims: bool = True,
    **kwargs,
):
    """Load a ``SummarizedExperiment`` or ``SingleCellExperiment`` object from a file.

    Args:
        path:
            Path to the dataset.

        metadata:
            Metadata for the dataset.
            Defaults to None.

        scrnaseq_realize_assays:
            Whether to realize assays into memory.
            Defaults to False.

        scrnaseq_realize_reduced_dims:
            Whether to realize reduced dimensions into memory.
            Defaults to True.

        **kwargs:
            Further arguments to pass to
            :py:func:`~dolomite_base.read_object.read_object`.

    Returns:
        A `SummarizedExperiment` of the object.
    """
    obj = read_object(
        path,
        metadata=metadata,
        scrnaseq_realize_assays=scrnaseq_realize_assays,
        scrnaseq_realize_reduced_dims=scrnaseq_realize_reduced_dims,
        **kwargs,
    )

    if isinstance(obj, SummarizedExperiment):
        if scrnaseq_realize_assays:
            _assays = {}
            for y in obj.get_assay_names():
                _assays[y] = realize_array(obj.assay(y))

            obj = obj.set_assays(_assays)

        if isinstance(obj, SingleCellExperiment):
            if scrnaseq_realize_reduced_dims:
                _red_dims = {}
                for z in obj.get_reduced_dim_names():
                    _red_dims[z] = realize_array(obj.reduced_dim(z))

                obj = obj.set_reduced_dims(_red_dims)

    return obj


def realize_array(x):
    """
    Realize a `ReloadedArray` into a dense array or sparse matrix.

    Args:
        x:
            `ReloadedArray` object.

    Returns:

        Realized array or matrix.
    """
    from dolomite_matrix import ReloadedArray

    if isinstance(x, ReloadedArray):
        if is_sparse(x):
            x = to_scipy_sparse_matrix(x, "csr")
        else:
            x = to_dense_array(x)

    return x
