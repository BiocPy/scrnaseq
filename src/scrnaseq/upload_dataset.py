from gypsum_client import cache_directory, rest_url, upload_directory

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def upload_dataset(
    directory: str,
    name: str,
    version: str,
    package: str = "scRNAseq",
    cache_dir: str = cache_directory(),
    deduplicate: bool = True,
    probation: bool = False,
    url: str = rest_url(),
    token: str = None,
    concurrent: int = 1,
    abort_failed: bool = True,
):
    """Upload the dataset to the gypsum bucket.

    This is a wrapper around
    :py:func:`~gypsum_client.upload_file_actions.upload_directory`
    specific to the `scRNAseq` package.

    See Also:
        :py:func:`~gypsum_client.upload_file_actions.upload_directory`,
        to upload a directory to the gypsum backend.

    Args:
        Name:
            Dataset name.

        version:
            Version name.

        directory:
            Path to a directory containing the ``files`` to be uploaded.
            This directory is assumed to correspond to a version of an asset.

        cache_dir:
            Path to the cache for saving files, e.g., in
            :py:func:`~gypsum_client.save_operations.save_version`.

            Used to convert symbolic links to upload links,see
            :py:func:`~gypsum_client.prepare_directory_for_upload.prepare_directory_upload`.

        deduplicate:
            Whether the backend should attempt deduplication of ``files``
            in the immediately previous version.
            Defaults to True.

        probation:
            Whether to perform a probational upload.
            Defaults to False.

        url:
            URL of the gypsum REST API.

        token:
            GitHub access token to authenticate to the gypsum REST API.

        concurrent:
            Number of concurrent downloads.
            Defaults to 1.

        abort_failed:
            Whether to abort the upload on any failure.

            Setting this to `False` can be helpful for diagnosing upload problems.

    Returns:
        `True` if successfull, otherwise `False`.
    """
    return upload_directory(
        directory,
        package,
        name,
        version,
        cache_dir=cache_dir,
        deduplicate=deduplicate,
        probation=probation,
        url=url,
        token=token,
        concurrent=concurrent,
        abort_failed=abort_failed,
    )
