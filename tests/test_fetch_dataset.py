import pytest
from scrnaseq import fetch_dataset
from singlecellexperiment import SingleCellExperiment

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def test_fetch_dataset():
    sce = fetch_dataset("zeisel-brain-2015", "2023-12-14")

    print(sce)
    assert isinstance(sce, SingleCellExperiment)
