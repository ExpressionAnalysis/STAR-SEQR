---------------
Developer Notes
---------------

General Procedures
------------------

Versioning:
 * Follows semantic versioning (major.minor.patch)
 * Pre-release candidates (major.minor.path.dev0)


Contributing:
 * Always use PRs and add tests when appropriate
 * Merges will be considered after travis-ci passes
 * Please add updates to **CHANGES.rst** which will be used for release notes


Local Tests:
 * python setup.py build
 * python setup.py install
 * nosetests
 * check py2 and py3 environments if possible

Deployment
----------

PypI Release:
 * First make sure your PyPI config is properly setup::

    [distutils]
    index-servers =
      pypi
      pypitest

    [pypi]
    repository=https://upload.pypi.org/legacy/
    username=<NAME>
    password=<PASS>

    [pypitest]
    repository=https://test.pypi.org/legacy/
    username=<NAME>
    password=<PASS>

 * All releases are uploaded to PyPI using twine::

        # test
        python setup.py sdist
        twine upload dist/* -r pypitest
        pip install -U --pre -i https://test.pypi.org/legacy/  starseqr

        # final
        python setup.py sdist
        twine upload dist/* -r pypi
        pip install -U starseqr


Docker:
 * Update release in docker/Dockerfile
 * The following scripts get the version tag from the package and tag the docker image::

        sudo bash build_docker.sh
        sudo bash push_docker.sh





Bioconda:
 * Notes

