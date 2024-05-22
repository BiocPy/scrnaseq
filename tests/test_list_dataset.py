import tempfile

import pandas as pd
from scrnaseq import list_datasets

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def test_list_dataset():
    datasets = list_datasets(cache_dir=tempfile.mkdtemp())

    assert isinstance(datasets, pd.DataFrame)
    assert len(datasets) > 80
