import json
import os
import shutil
from functools import singledispatch
from typing import Any

import dolomite_base as dl
from gypsum_client import fetch_metadata_schema, validate_metadata
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment

from .utils import format_object_metadata

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


@singledispatch
def save_dataset(x: Any, path, metadata):
    """Save a dataset to disk.

    Args:
        x:
            An object containing single-cell data.
            May be a derivative of
            :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`
            or :py:class:`~anndata.AnnData`.

        path:
            Path to a new directory to save the dataset.

        metadata:
            Dictionary containing the metadata for this dataset.
            see the schema returned by
            :py:func:`~gypsum_client.fetch_metadata_schema.fetch_metadata_schema`.

            Note that the ``applications.takane`` property will be automatically
            added by this function and does not have to be supplied.

    See Also:
        `metadata index <https://github.com/ArtifactDB/bioconductor-metadata-index>`_,
        on the expected schema for the metadata.

        :py:func:`~scrnaseq.polish_dataset.polish_dataset`,
        to polish ``x`` before saving it.

        :py:func:`~scrnaseq.upload_dataset.upload_dataset`, to upload the saved contents.

    Example:

        .. code-block:: python

            # Fetch an existing dataset
            # or create your own ``SingleCellExperiment``
            # or ``AnnData`` object.
            sce = scrnaseq.fetch_dataset(
                "zeisel-brain-2015",
                "2023-12-14",
            )

            # Provide dataset level metadata for search and findability
            meta = {
                "title": "My dataset made from ziesel brain",
                "description": "This is a copy of the ziesel",
                "taxonomy_id": [
                    "10090"
                ],  # NCBI ID
                "genome": [
                    "GRCh38"
                ],  # genome build
                "sources": [
                    {
                        "provider": "GEO",
                        "id": "GSE12345",
                    }
                ],
                "maintainer_name": "Shizuka Mogami",
                "maintainer_email": "mogami.shizuka@765pro.com",
            }

            import shutil
            import tempfile

            cache_dir = tempfile.mkdtemp()

            # Make sure the directory is clean
            shutil.rmtree(
                cache_dir
            )

            # Save the dataset
            scrnaseq.save_dataset(
                sce,
                cache_dir,
                meta,
            )
    """
    raise NotImplementedError(f"'save_dataset' is not supported for objects of class: {type(x)}")


def _save_se(x, path, metadata):
    schema = fetch_metadata_schema()

    if "bioconductor_version" not in metadata:
        metadata["bioconductor_version"] = "3.19"  # current release

    validate_metadata(metadata, schema)

    if os.path.exists(path):
        shutil.rmtree(path)

    dl.save_object(x, path, reloaded_array_reuse_mode="symlink")

    takane = format_object_metadata(x)
    takane["type"] = dl.read_object_file(path)["type"]

    if "applications" not in metadata:
        metadata["applications"] = {}

    metadata["applications"]["takane"] = takane

    # Second validation with the takane metadata.
    contents = json.dumps(metadata, indent=4)
    validate_metadata(json.loads(contents), schema=schema)
    with open(os.path.join(path, "_bioconductor.json"), "w") as f:
        f.write(contents)


@save_dataset.register
def save_dataset_sce(x: SingleCellExperiment, path: str, metadata: dict):
    """Save :py:class:`~singlecellexperiment.SingleCellExperiment.SingleCellExperiment` to disk."""
    return _save_se(x, path, metadata)


@save_dataset.register
def save_dataset_se(x: SummarizedExperiment, path: str, metadata: dict):
    """Save :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment` to disk."""
    return _save_se(x, path, metadata)


has_anndata = False
try:
    import anndata

    has_anndata = True
except Exception:
    pass

if has_anndata:

    @save_dataset.register
    def save_dataset_anndata(x: anndata.AnnData, path: str, metadata: dict):
        """Save :py:class:`~anndata.AnnData` to disk."""
        _sce = SingleCellExperiment.from_anndata(x)
        return save_dataset(_sce, path, metadata)
