# Contributing

--8<-- "CONTRIBUTING.md:4:63"

## Developer guide

### Setting up a development environment

```doscon
git clone https://github.com/kcibv/pyflange.git
cd pyflange
py -m pip install -e .
py -m pip install pytest -r docsrc/requirements.txt
```

The `-e` option installs PyFlange in *editable* mode: changes to the source
code in `src/pyflange` are used immediately, without reinstalling.

### Running the tests and the validation

```doscon
py -m pytest
py validations/validate.py validations/references/*.json --failed
```

The tests also run the example scripts of the tutorials
(`docsrc/examples/*.py`). See [Validation](guide/validation.md).

### Writing documentation

The documentation is written in Markdown in the `docsrc` folder and built
with [MkDocs](https://www.mkdocs.org/) and the
[Material](https://squidfunk.github.io/mkdocs-material/) theme. The API
reference is generated from the docstrings by
[mkdocstrings](https://mkdocstrings.github.io/).

```doscon
mkdocs serve            # live preview at http://127.0.0.1:8000
mkdocs build --strict   # checks and builds the site into the site folder
```

(If the `mkdocs` command is not found, use `py -m mkdocs`.)

Every change to the code should come with the corresponding documentation
change, in the same pull request:

- **Docstrings** follow the
  [Google style](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings):
  a summary line, a description, then `Args:`, `Returns:` (or `Yields:`),
  `Raises:` and, where useful, `Example:` sections. Attributes of dataclasses
  are documented in an `Attributes:` section of the class docstring.
- Docstrings cite the documents of the [References](references.md) page as
  `ref. [n]`, with the same numbers in every module, followed by the section
  or equation, e.g. `ref. [2], eq. G.36`. Each module docstring lists the
  references it cites in a `References:` section. A new reference gets the
  next free number and is added to `docsrc/references.md`.
- New public functions or classes appear automatically in the API
  reference. New modules need a page in `docsrc/apis/` and an entry in the
  `nav` section of `mkdocs.yml`.
- Changes to the command line must be reflected in `docsrc/cli.md`, changes
  to the case files in `docsrc/json-format.md`, and user-visible changes in
  `docsrc/changelog.md`.
- Tutorial code lives in `docsrc/examples/*.py` and is included in the pages
  with `--8<-- "examples/file.py:section"`. Update the outputs shown on the
  pages when the results change, and regenerate the figures with
  `py docsrc/assets/img/make_figures.py`. The drawings of the flange
  geometry are inline SVG files in `docsrc/assets/img/*-geometry.svg`,
  styled by `docsrc/assets/css/extra.css` so that they follow the light and
  dark themes.

### Publishing and versions

The site is published on GitHub Pages from the `docs` folder of the `main`
branch, with one version of the documentation per release (`docs/0.13/`,
`docs/0.14/`, ...), a `latest` alias pointing to the newest release, and a
`dev` version for the development branch. The reader chooses the version from
the selector next to the site name.

Publishing is automatic: the *Documentation* GitHub workflow uses
[mike](https://github.com/jimporter/mike) to build the site and commit it to
`main`:

- at every push to `main`, as version *major.minor* of the version in
  `pyproject.toml` (e.g. `0.13` for 0.13.2), moving `latest` to it;
- at every push to `dev`, as version `dev`.

So, to publish the documentation of a new release, bump the version in
`pyproject.toml` and merge into `main`. Never edit the `docs` folder by hand,
and do not run `mike` on a local `main` branch. To remove an obsolete
version, run `mike delete --branch main --deploy-prefix docs --push VERSION`
from a clean clone.

The *Tests* workflow runs the test suite on Linux and Windows, with the
oldest and newest supported Python versions, at every push and pull request
to `main` and `dev`.

### Releasing on PyPI

The *Publish to PyPI* workflow publishes the package on
[PyPI](https://pypi.org/project/pyflange/) when a tag named
`pyflange-v<version>` is pushed. To release, say, version 0.14.0:

1. set `version = "0.14.0"` in `pyproject.toml`, date the release in
   `docsrc/changelog.md` and merge into `main`;
2. on `main`, run the tests and build the package locally:

    ```doscon
    py -m pytest
    py -m build
    py -m twine check --strict dist/pyflange-0.14.0*
    ```

3. push `main`, then tag the merged commit and push the tag:

    ```doscon
    git push origin main
    git tag -a pyflange-v0.14.0 -m "PyFlange 0.14.0"
    git push origin pyflange-v0.14.0
    ```

The workflow does not run the tests and does not check the branch: it only
checks that the tag and `pyproject.toml` versions are equal, then builds the
package, checks it with `twine check --strict` and uploads it. The upload
uses PyPI
[trusted publishing](https://docs.pypi.org/trusted-publishers/), so no API
token is stored in GitHub. After the upload, the workflow creates the GitHub
release of the tag, with the packages attached and the changelog section of
the version (`## 0.14.0 (...)`) as release notes; the version shown next to
the repository link at the top of this site is the one of the latest GitHub
release. A version cannot be uploaded to PyPI twice: to fix
a release, publish a new patch version.
