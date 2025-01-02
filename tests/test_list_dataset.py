import tempfile

from biocframe import BiocFrame
from scrnaseq import list_datasets

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def test_list_dataset():
    datasets = list_datasets(cache_dir=tempfile.mkdtemp())

    assert isinstance(datasets, BiocFrame)
    assert len(datasets) > 80
