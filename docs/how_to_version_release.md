# PySEP Version Release Checklist

Checklist for any version release relating to semantic version number
incrementation. Please include in any PR to master and make sure all points
are checked off when incrementing any of the version numbers (major, minor, 
patch), or provide reason in PR for why certain points are not checked.

## Prior to PR merge:
- [] Merge `devel` -> `master`
- [] Ensure `devel` is up to date with `master` ($ git merge master)
- [] Bump version number `pyproject.toml`
- [] Bump version number `docs/conf.py`
- [] Bump version number `pysep/__init__.py`
- [] Ensure all tests still pass, fix broken tests
- [] Update `CHANGELOG` to include all major changes since last version

## Following PR merge:
- [] Create GitHub version release
- [] Publish latest version on PyPi
- [] Update PyPi badge in `README.md`
- [] Post on adjTomo Discussion for major and minor version releases
