import json
import os
import shutil
import tempfile

import anndata as ad
import dolomite_base as dl
import numpy as np
import pandas as pd
import pytest
from dolomite_matrix import ReloadedArray
from scrnaseq import fetch_dataset, save_dataset
from delayedarray import to_dense_array

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def test_save_dataset_sce():
    sce = fetch_dataset("zeisel-brain-2015", "2023-12-14")

    meta = {
        "title": "My dataset forked from ziesel brain",
        "description": "This is a copy of the ziesel",
        "taxonomy_id": ["10090"],  # NCBI ID
        "genome": ["GRCh38"],  # genome build
        "sources": [{"provider": "GEO", "id": "GSE12345"}],
        "maintainer_name": "Shizuka Mogami",
        "maintainer_email": "mogami.shizuka@765pro.com",
    }

    tmp = tempfile.mkdtemp()
    save_dataset(sce, tmp, meta)

    # Load the saved AnnData object
    roundtrip = dl.read_object(tmp)

    assert len(roundtrip.get_column_data()) == len(sce.get_column_data())
    assert isinstance(roundtrip.get_assays()["counts"], ReloadedArray)
    assert isinstance(sce.get_assays()["counts"], ReloadedArray)

    # Load and check the metadata
    with open(os.path.join(tmp, "_bioconductor.json")) as f:
        saved_meta = json.load(f)

    assert saved_meta["bioconductor_version"] == "3.19"

    # Test validation failure
    meta["title"] = 1234
    with pytest.raises(Exception):
        save_dataset(sce, tmp, meta)

    shutil.rmtree(tmp)


def test_save_dataset_anndata():
    data = np.random.poisson(1, (10, 100))
    adata = ad.AnnData(data)
    adata.obs["foo"] = np.random.choice(list("abcdefghijklmnopqrstuvwxyz"), 10)
    adata.var_names = [f"GENE_{i+1}" for i in range(adata.n_vars)]
    adata.obs_names = list("ABCDEFGHIJ")
    adata.layers["counts"] = data

    meta = {
        "title": "My dataset forked from ziesel brain",
        "description": "This is a copy of the ziesel",
        "taxonomy_id": ["10090"],  # NCBI ID
        "genome": ["GRCh38"],  # genome build
        "sources": [{"provider": "GEO", "id": "GSE12345"}],
        "maintainer_name": "Shizuka Mogami",
        "maintainer_email": "mogami.shizuka@765pro.com",
    }

    tmp = tempfile.mkdtemp()
    save_dataset(adata, tmp, meta)

    # Load the saved AnnData object
    roundtrip = dl.read_object(tmp)

    assert len(roundtrip.get_column_data()) == 10
    assert isinstance(roundtrip.get_assays()["counts"], ReloadedArray)
    assert isinstance(adata.layers["counts"], np.ndarray)
    assert np.array_equal(
        to_dense_array(roundtrip.get_assays()["counts"]).transpose(),
        adata.layers["counts"],
    )

    # Load and check the metadata
    with open(os.path.join(tmp, "_bioconductor.json")) as f:
        saved_meta = json.load(f)

    assert saved_meta["bioconductor_version"] == "3.19"

    # Test validation failure
    meta["title"] = 1234
    with pytest.raises(Exception):
        save_dataset(adata, tmp, meta)

    shutil.rmtree(tmp)
