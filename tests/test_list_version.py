from scrnaseq import fetch_latest_version, list_versions

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def test_list_versions():
    versions = list_versions("romanov-brain-2017")

    assert isinstance(versions, list)
    assert "2023-12-19" in versions


def test_latest_version():
    version = fetch_latest_version("romanov-brain-2017")

    assert isinstance(version, str)
    assert "2023-12-19" == version
