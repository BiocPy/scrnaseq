import json
import os
import shutil
import tempfile

import anndata as ad
import dolomite_base as dl
import dolomite_matrix as dlm
import numpy as np
import pandas as pd
import pytest
from biocframe import BiocFrame
from gypsum_client import prepare_directory_upload
from scipy import sparse as sp
from scrnaseq import fetch_dataset, save_dataset, upload_dataset
from singlecellexperiment import SingleCellExperiment

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def generate_sce():
    mat = np.random.poisson(1, (100, 10))
    row_names = [f"GENE_{i}" for i in range(mat.shape[0])]
    col_names = list("ABCDEFGHIJ")
    sce = SingleCellExperiment(
        assays={"counts": mat},
        row_data=BiocFrame(row_names=row_names),
        column_data=BiocFrame(row_names=col_names),
    )
    return sce


def test_file_listing_without_special_elements():
    sce = generate_sce()

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

    listing = prepare_directory_upload(tmp, links="always")
    assert len(listing["links"]) == 0
    assert sorted(listing["files"]) == sorted(
        [
            os.path.relpath(os.path.join(dp, f), tmp)
            for dp, dn, filenames in os.walk(tmp)
            for f in filenames
        ]
    )

    shutil.rmtree(tmp)


def test_file_listing_with_reloaded_arrays():
    sce = generate_sce()

    cache = tempfile.mkdtemp()
    tmp0 = os.path.join(cache, "bucket", "scRNAseq", "test", "foo")
    os.makedirs(os.path.dirname(tmp0), exist_ok=True)

    meta = {
        "title": "My dataset forked from ziesel brain",
        "description": "This is a copy of the ziesel",
        "taxonomy_id": ["10090", "9606"],  # NCBI ID
        "genome": ["GRCh38", "GRCm38"],
        "sources": [
            {"provider": "GEO", "id": "GSE12345"},
            {"provider": "DOI", "id": "123asd/231.123"},
        ],
        "maintainer_name": "Kaori Sakuramori",
        "maintainer_email": "sakuramori.kaori@765pro.com",
    }

    save_dataset(sce, tmp0, meta)

    sce2 = sce.copy()
    sce2.assays["reloaded_assay"] = dlm.ReloadedArray(
        path=f"{tmp0}/assays/0", seed=np.random.poisson(1, (100, 10))
    )

    tmp = tempfile.mkdtemp()
    save_dataset(sce2, tmp, meta)
    with pytest.raises(Exception):
        listing = prepare_directory_upload(tmp, links="always")

    listing = prepare_directory_upload(tmp, links="always", cache_dir=cache)
    assert sorted(x["to.path"] for x in listing["links"]) == sorted(
        ["assays/0/array.h5", "assays/0/OBJECT"]
    )
    assert sorted(x["from.path"] for x in listing["links"]) == sorted(
        ["assays/1/array.h5", "assays/1/OBJECT"]
    )
    assert sorted(
        listing["files"] + [x["from.path"] for x in listing["links"]]
    ) == sorted(
        [
            os.path.relpath(os.path.join(dp, f), tmp)
            for dp, dn, filenames in os.walk(tmp)
            for f in filenames
        ]
    )

    shutil.rmtree(tmp)
    shutil.rmtree(cache)


@pytest.mark.skipif(
    "gh_token" not in os.environ, reason="GitHub token not in environment"
)
def test_actual_upload_works_correctly():
    sce = generate_sce()

    gh_token = os.environ.get("gh_token", None)
    if gh_token is None:
        raise ValueError("GitHub token not in environment")

    meta = {
        "title": "My dataset",
        "description": "This is my dataset",
        "taxonomy_id": ["10090"],
        "genome": ["GRCh38"],
        "sources": [{"provider": "GEO", "id": "GSE12345"}],
        "maintainer_name": "Shizuka Mogami",
        "maintainer_email": "mogami.shizuka@765pro.com",
    }

    tmp = tempfile.mkdtemp()
    save_dataset(sce, tmp, meta)

    app_url = "https://gypsum.artifactdb.com"

    version = str(pd.Timestamp.today().date())
    upload_dataset(tmp, "test", version, probation=True, url=app_url, token=gh_token)
    fetch_dataset.cache_clear()  # Clear cache before fetching

    cache = tempfile.mkdtemp()
    roundtrip = fetch_dataset("test", version, cache=cache, url=app_url)

    assert roundtrip.get_column_data().equals(sce.get_column_data())
    assert np.array_equal(roundtrip.assays["counts"], sce.assays["counts"])
    assert isinstance(roundtrip.assay("counts"), type(sce.assay("counts")))

    shutil.rmtree(tmp)
    shutil.rmtree(cache)
