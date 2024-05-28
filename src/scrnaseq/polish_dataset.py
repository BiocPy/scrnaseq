from typing import Type

import numpy as np
from scipy import sparse as sp
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def polish_dataset(
    x: Type[SummarizedExperiment],
    reformat_assay_by_density: float = 0.3,
    attempt_integer_conversion: bool = True,
    remove_altexp_coldata: bool = True,
    forbid_nested_altexp: bool = True,
) -> Type[SummarizedExperiment]:
    """Optimize dataset for saving.

    Prepare a
    :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment` or
    :py:class:`~singlecellexperiment.SingleCellExperiment.SingleCellExperiment`
    to be saved with :py:func:`scrnaseq.save_dataset.save_dataset`.

    This performs minor changes to improve storage efficiency, especially
    with matrices.

    Args:
        x:
            A :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`
            or one of its derivative.

        reformat_assay_by_density:
            Whether to optimize assay formats based on the density of non-zero values.
            Assays with densities above this number are converted to ordinary dense
            arrays (if they are not already), while those with lower densities are
            converted to sparse matrices.

            This can be disabled by setting it to `None`.

        attempt_integer_conversion:
            Whether to convert double-precision assays containing integer values
            to actually have the integer type.

            This can improve efficiency of downstream applications by avoiding
            the need to operate in double precision.

        remove_altexp_coldata:
            Whether column data for alternative experiments should be removed.
            Defaults to `True` as the alternative experiment column data is
            usually redundant compared to the main experiment.

        forbid_nested_altexp:
            Whether nested alternative experiments (i.e., alternative experiments of
            alternative experiments) should be forbidden.

    Returns:
        A modifed object with the same type as ``x``.
    """
    return _polish_dataset(
        x,
        reformat_assay_by_density,
        attempt_integer_conversion,
        remove_altexp_coldata,
        forbid_nested_altexp,
    )


def _polish_dataset(
    x: Type[SummarizedExperiment],
    reformat_assay_by_density: float,
    attempt_integer_conversion: bool,
    remove_altexp_coldata: bool,
    forbid_nested_altexp: bool,
    level: int = 0,
):
    new_assays = {}
    for asyname, asy in x.assays.items():
        if reformat_assay_by_density is not None:
            density = min(np.mean(asy != 0), np.mean(asy != np.nan))
            if density < reformat_assay_by_density:
                if not sp.issparse(asy):
                    asy = sp.csr_matrix(asy)
            else:
                if sp.issparse(asy):
                    asy = asy.toarray()

        if attempt_integer_conversion:
            if np.issubdtype(asy.dtype, np.floating):
                _cast = False
                if sp.issparse(asy):
                    if not np.any(asy.data % 1 != 0):
                        _cast = True
                elif not np.any(asy % 1 != 0):
                    _cast = True

                if _cast is True:
                    asy = asy.astype(np.int_)

        new_assays[asyname] = asy

    x = x.set_assays(new_assays)

    if isinstance(x, SingleCellExperiment):
        if len(x.get_alternative_experiment_names()) > 0:
            if forbid_nested_altexp and level > 0:
                raise ValueError("Nested alternative experiments are forbidden.")

        new_alts = {}
        for altname, altexp in x.alternative_experiments.items():
            if remove_altexp_coldata:
                altexp = altexp.set_column_data(None)

            altexp = _polish_dataset(
                altexp,
                reformat_assay_by_density=reformat_assay_by_density,
                attempt_integer_conversion=attempt_integer_conversion,
                remove_altexp_coldata=remove_altexp_coldata,
                forbid_nested_altexp=forbid_nested_altexp,
                level=level + 1,
            )

            new_alts[altname] = altexp

        x = x.set_alternative_experiments(new_alts)

    return x
