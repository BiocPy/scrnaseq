import numpy as np
import pytest
from biocframe import BiocFrame
from scipy import sparse
from scrnaseq import fetch_dataset, polish_dataset
from singlecellexperiment import SingleCellExperiment

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def generate_matrix():
    mat = np.random.poisson(1, (100, 10))
    return mat


def test_polish_dataset_strips_assay_dimnames():
    mat = generate_matrix()
    row_names = [f"GENE_{i}" for i in range(mat.shape[0])]
    col_names = list("ABCDEFGHIJ")
    sce = SingleCellExperiment(
        assays={"counts": mat},
        row_data=BiocFrame(row_names=row_names),
        column_data=BiocFrame(row_names=col_names),
    )

    y = polish_dataset(sce)
    assert sce.shape == y.shape
    assert sce.get_assay_names() == y.get_assay_names()
    assert y.get_row_names() is not None
    assert y.get_column_names() is not None

    sce.row_data = sce.row_data.set_row_names([name.lower() for name in row_names])
    sce.column_names = sce.column_data.set_row_names(
        [name.lower() for name in col_names]
    )
    y = polish_dataset(sce)
    assert sce.shape == y.shape
    assert sce.get_assay_names() == y.get_assay_names()
    assert y.get_row_names() == sce.get_row_names()
    assert y.get_column_names() == sce.get_column_names()


def test_polish_dataset_strips_reduced_dimension_names():
    mat = generate_matrix()
    row_names = [f"GENE_{i}" for i in range(mat.shape[0])]
    col_names = list("ABCDEFGHIJ")
    pca = np.random.randn(10, 5)
    sce = SingleCellExperiment(
        assays={"counts": mat},
        row_data=BiocFrame(row_names=row_names),
        column_data=BiocFrame(row_names=col_names),
        reduced_dims={"pca": pca},
    )

    y = polish_dataset(sce)
    assert sce.shape == y.shape
    assert sce.get_assay_names() == y.get_assay_names()
    assert y.get_row_names() is not None
    assert y.get_column_names() is not None
    assert sce.get_reduced_dim_names() == y.get_reduced_dim_names()
    assert y.reduced_dims["pca"].shape == pca.shape


def test_polish_dataset_strips_alternative_experiment_names():
    mat = generate_matrix()
    row_names = [f"GENE_{i}" for i in range(mat.shape[0])]
    col_names = list("ABCDEFGHIJ")
    pca = np.random.randn(10, 5)
    altexp = SingleCellExperiment(
        assays={"counts": mat},
        row_data=BiocFrame(row_names=row_names),
        column_data=BiocFrame(row_names=col_names),
    )
    sce = SingleCellExperiment(
        assays={"counts": mat},
        row_data=BiocFrame(row_names=row_names),
        column_data=BiocFrame(row_names=col_names),
        reduced_dims={"pca": pca},
        alternative_experiments={"rando": altexp},
    )

    y = polish_dataset(sce)
    assert y.alternative_experiments["rando"].shape == altexp.shape
    assert y.alternative_experiments["rando"].get_row_names() == altexp.get_row_names()


def test_polish_dataset_converts_dense_to_sparse():
    mat = np.random.poisson(0.2, (100, 10))
    sce = SingleCellExperiment(assays={"counts": mat})

    y = polish_dataset(sce)
    assert sparse.issparse(y.assays["counts"])

    y = polish_dataset(sce, reformat_assay_by_density=0)
    assert not sparse.issparse(y.assays["counts"])
    assert y.shape == mat.shape

    y = polish_dataset(sce, reformat_assay_by_density=None)
    assert not sparse.issparse(y.assays["counts"])


def test_polish_dataset_converts_sparse_to_dense():
    mat = np.random.poisson(3, (100, 10))
    sce = SingleCellExperiment(assays={"counts": sparse.csr_matrix(mat)})

    y = polish_dataset(sce)
    assert not sparse.issparse(y.assays["counts"])

    y = polish_dataset(sce, reformat_assay_by_density=1)
    assert sparse.issparse(y.assays["counts"])

    y = polish_dataset(sce, reformat_assay_by_density=None)
    assert sparse.issparse(y.assays["counts"])


def test_polish_dataset_attempts_integer_conversions():
    mat = np.random.poisson(3, (100, 10)).astype(float)
    sce = SingleCellExperiment(assays={"counts": mat})

    y = polish_dataset(sce)
    assert np.issubdtype(y.assays["counts"].dtype, np.integer)

    mat = np.random.poisson(0.1, (100, 10)).astype(float)
    sce = SingleCellExperiment(assays={"counts": sparse.csr_matrix(mat)})

    y = polish_dataset(sce)
    assert sparse.issparse(y.assays["counts"])
    assert np.issubdtype(y.assays["counts"].dtype, np.integer)

    mat = np.random.poisson(3, (100, 10)) * 1.5
    sce = SingleCellExperiment(assays={"counts": mat})

    y = polish_dataset(sce)
    assert np.issubdtype(y.assays["counts"].dtype, np.floating)


def test_polish_dataset_works_with_na_values():
    mat = np.random.poisson(0.1, (100, 10))
    mat = mat.astype(np.float64)
    mat.ravel()[np.random.choice(mat.size, 10, replace=False)] = np.nan
    sce = SingleCellExperiment(assays={"counts": sparse.csr_matrix(mat)})

    y = polish_dataset(sce)
    assert sparse.issparse(y.assays["counts"])
    assert np.issubdtype(y.assays["counts"].dtype, np.floating)


def test_polish_dataset_forbids_highly_nested_altexps():
    mat = np.random.poisson(1, (100, 10))
    sce2 = SingleCellExperiment(assays={"counts": mat[:5] + 2})

    sce1 = SingleCellExperiment(
        assays={"counts": mat[:10] + 1}, alternative_experiments={"rando": sce2}
    )

    sce0 = SingleCellExperiment(
        assays={"counts": mat}, alternative_experiments={"rando": sce1}
    )

    with pytest.raises(Exception):
        polish_dataset(sce0)

    y = polish_dataset(sce0, forbid_nested_altexp=False)
    assert (
        y.get_alternative_experiment_names() == sce0.get_alternative_experiment_names()
    )


@pytest.mark.skip("takes too long")
def test_polish_existing_dataset():
    sce = fetch_dataset("zeisel-brain-2015", "2023-12-14")

    y = polish_dataset(sce)

    assert y.shape == sce.shape
    assert type(y.assays["counts"]) != type(sce.assays["counts"])
