from typing import List

import gypsum_client as gypc

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def list_versions(name: str) -> List[str]:
    """List all available versions for a dataset.

    Example:

        .. code-block:: python

            versions = list_versions(
                "romanov-brain-2017"
            )

    Args:
        name:
            Name of the dataset.

    Returns:
        A list of version names.
    """
    return gypc.list_versions("scRNAseq", name)


def fetch_latest_version(name: str) -> str:
    """Fetch latest version for a dataset.

    Example:

        .. code-block:: python

            version = fetch_latest_version(
                "romanov-brain-2017"
            )

    Args:
        name:
            Name of the dataset.

    Returns:
        Latest version name.
    """
    return gypc.fetch_latest("scRNAseq", name)
