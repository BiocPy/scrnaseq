import numpy as np
import pytest
import scipy.sparse as sp
from delayedarray import is_sparse
from dolomite_matrix import ReloadedArray
from scrnaseq import fetch_dataset, fetch_metadata
from singlecellexperiment import SingleCellExperiment

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def test_fetch_dataset():
    sce = fetch_dataset("zeisel-brain-2015", "2023-12-14")
    assert isinstance(sce, SingleCellExperiment)

    # Correctly creates ReloadedMatrix objects.
    ass = sce.get_assays()
    assert all(isinstance(a, ReloadedArray) for _, a in ass.items())
    assert all(is_sparse(a) for _, a in ass.items())

    ass_0 = ass["counts"]
    assert "zeisel-brain-2015" in ass_0.seed.path
    assert "2023-12-14" in ass_0.seed.path

    # Works with realization options.
    sce = fetch_dataset("zeisel-brain-2015", "2023-12-14", realize_assays=True)
    ass = sce.get_assays()
    assert all(isinstance(a, (sp.csc_matrix, sp.csr_matrix)) for _, a in ass.items())

    alt_exps = sce.get_alternative_experiments()
    for altname, alt in alt_exps.items():
        alt_exp_ass = alt.get_assays()
        assert all(isinstance(a, (np.ndarray)) for _, a in alt_exp_ass.items())


@pytest.mark.skip("takes too long")
def test_fetch_dataset_realizes_reduced_dimensions():
    sce = fetch_dataset("aztekin-tail-2019", "2023-12-14", realize_reduced_dims=False)
    red_dim = sce.get_reduced_dims()
    assert all(isinstance(a, ReloadedArray) for _, a in red_dim.items())

    sce = fetch_dataset("aztekin-tail-2019", "2023-12-14", realize_reduced_dims=True)
    red_dim = sce.get_reduced_dims()
    assert all(isinstance(a, np.ndarray) for _, a in red_dim.items())


def test_fetch_metadata():
    meta = fetch_metadata("zeisel-brain-2015", "2023-12-14")
    assert "Brain structure" in meta["title"]
    assert meta["taxonomy_id"][0] == "10090"
