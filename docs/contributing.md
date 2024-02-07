<!--
SPDX-FileCopyrightText: 2024 Helmholtz-Zentrum hereon GmbH

SPDX-License-Identifier: CC-BY-4.0
-->

(contributing)=
# Contribution and development hints

```{warning}

This page has been automatically generated as has not yet been reviewed
by the authors of generalized-aquatic-ecosystem-model!
```

The _generalized-aquatic-ecosystem-model_ project is developed by the
Helmholtz-Zentrum Hereon[institution-url]. It is open-source
as we believe that this analysis can be helpful for reproducibility and
collaboration, and we are looking forward for your feedback,
questions and especially for your contributions.

- If you want to ask a question, are missing a feature or have
  comments on the docs, please [open an issue at the source code
  repository][issues]
- If you have suggestions for improvement, please let us know in an
  issue, or fork the repository and create a merge request. See also
  {ref}`development`.

[institution-url]: https://www.hereon.de
[issues]: https://codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model

(development)=
## Contributing in the development

```{info}
We use automated formatters to ensure a high quality and maintanability of
our source code. Getting familiar with these techniques can take quite some
time and you might get error messages that are hard to understand.

We not slow down your development and we do our best to support you with
these techniques. If you have any troubles, just commit with
`git commit --no-verify` (see below) and the maintainers will take care
of the tests and continuous integration.
```

Thanks for your wish to contribute to this project!! The source code of
the _generalized-aquatic-ecosystem-model_ package is hosted at https://codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model.


This is an open gitlab where you can register via the Helmholtz AAI. If your
home institution is not listed in the Helmholtz AAI, please use one of the
social login providers, such as Google, GitHub or OrcID.


Once you created an account in this gitlab, you can [fork][fork] this
repository to your own user account and implement the changes.

Afterwards, please make a merge request into the main repository. If you
have any questions, please do not hesitate to create an issue on gitlab
and contact the maintainers of this package.

Once you created you fork, you can clone it via

```bash
git clone https://codebase.helmholtz.cloud/<your-user>/generalized-aquatic-ecosystem-model.git
```

we recommend that you change into the directory and create a virtual
environment via:

```bash
cd generalized-aquatic-ecosystem-model
python -m venv venv
source venv/bin/activate # (or venv/Scripts/Activate.bat on windows)
```

and install it in development mode with the `[dev]` option via:

```bash
pip install -e ./generalized-aquatic-ecosystem-model/[dev]
```

[fork]: https://codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model/-/forks/new

## Helpers

### Shortcuts with make

There are several shortcuts available with the `Makefile` in the root of the
repository. On Linux, you can execute `make help` to get an overview.

### Annotating licenses

If you want to create new files, you need to set license and copyright
statements correctly. We use `reuse` to check that the licenses are
correctly encoded. As a helper script, you can use the script at
`.reuse/add_license.py` that provides several shortcuts from
`.reuse/shortcuts.yaml`. Please select the correct shortcut, namely

- If you create a new python file, you should run

  ```bash
  python .reuse/add_license.py code <file-you-created>.py
  ```
- If you created a new file for the docs, you should run

  ```bash
  python .reuse/add_license.py docs <file-you-created>.py
  ```
- If you created any other non-code file, you should run

  ```bash
  python .reuse/add_license.py supp <file-you-created>.py
  ```

If you have any questions on how licenses are handled, please do not hesitate
to contact the maintainers of _generalized-aquatic-ecosystem-model_.


## Fixing the docs

The documentation for this package is written in restructured Text and
built with [sphinx](https://www.sphinx-doc.org) and deployed on
[readthedocs](https://readthedocs.org).

If you found something in the docs that you want to fix, head over to
the `docs` folder, install the necessary requirements via
`pip install -r requirements.txt ../[docs]` and build the docs with `make html` (or
`make.bat` on windows). The docs are then available in
`docs/_build/html/index.html` that you can open with your local browser.

Implement your fixes in the corresponding `.rst`-file and push them to
your fork on gitlab.

## Contributing to the code

We use automated formatters (see their config in `pyproject.toml`), namely

- [Black](https://black.readthedocs.io/en/stable/) for standardized
  code formatting
- [blackdoc](https://blackdoc.readthedocs.io/en/latest/) for
  standardized code formatting in documentation
- [Flake8](http://flake8.pycqa.org/en/latest/) for general code
  quality
- [isort](https://github.com/PyCQA/isort) for standardized order in
  imports.
- [mypy](http://mypy-lang.org/) for static type checking on [type
  hints](https://docs.python.org/3/library/typing.html)
- [reuse](https://reuse.readthedocs.io/) for handling of licenses
- [cffconvert](https://github.com/citation-file-format/cff-converter-python)
  for validating the `CITATION.cff` file.

We highly recommend that you setup [pre-commit hooks](https://pre-commit.com/)
to automatically run all the above tools every time you make a git commit. This
can be done by running:

```bash
pre-commit install
```

from the root of the repository. You can skip the pre-commit checks with
`git commit --no-verify` but note that the CI will fail if it encounters
any formatting errors.

You can also run the `pre-commit` step manually by invoking:

```bash
pre-commit run --all-files
```

## Updating the skeleton for this package

This package has been generated from the template
https://codebase.helmholtz.cloud/hcdc/software-templates/python-package-template.git.

See the template repository for instructions on how to update the skeleton for
this package.
