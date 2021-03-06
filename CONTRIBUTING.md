# Contributing Guide

### Creating useful git commit messages
For most additions to this project, we expect users to provide useful/detailed commit messages.
Please refer to the commit messages in the current `git log` for examples.
[This is a great resource for creating detailed commit messages.](https://chris.beams.io/posts/git-commit/)

### Using `pre-commit` to auto-format commits

**NOTE** We cannot set up `pre-commit-ci` to auto format this repository at the moment since it is a private repo.
`pre-commit-ci` is free for **public** repositories, but not private ones like this project.

Users are expected to set up `pre-commit` on their local machines for the time being.

```bash
conda activate mosdef-study38
mamba install -c conda-forge pre-commit
pre-commit install --install-hooks
pre-commit run --all-files
```
