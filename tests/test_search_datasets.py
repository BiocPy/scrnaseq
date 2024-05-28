import pandas as pd
from gypsum_client import define_text_query
from scrnaseq import search_datasets

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def test_search_datasets():
    res = search_datasets("brain")
    assert len(res) > 10
    assert isinstance(res, pd.DataFrame)

    res = search_datasets(define_text_query("Neuro%", partial=True))
    assert isinstance(res, pd.DataFrame)
    assert len(res) > 0

    res = search_datasets(define_text_query("10090", field="taxonomy_id"))
    assert isinstance(res, pd.DataFrame)
    assert len(res) > 0

    res = search_datasets(
        define_text_query("GRCm38", field="genome")
        & (
            define_text_query("neuro%", partial=True)
            | define_text_query("pancrea%", partial=True)
        )
    )
    assert isinstance(res, pd.DataFrame)
    assert len(res) > 0
