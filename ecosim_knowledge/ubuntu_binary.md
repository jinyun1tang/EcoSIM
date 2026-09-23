# EcoSIM binary for Ubuntu 24.04 x86_64

This workflow targets the supplied machine: Ubuntu 24.04.3 LTS, x86_64,
glibc 2.39 (reported package revision `2.39-0ubuntu8.6`). It builds on the
Ubuntu 24.04 runner with GCC/GFortran 13 and generic x86-64 CPU flags.
It does not pin Ubuntu's security-update package revision or promise a
byte-for-byte reproducible build. The artifact records the actual versions.

## Build and download

1. Commit and push `.github/workflows/ecosim-ubuntu-binary.yml`,
   `.github/scripts/package-ubuntu-binary.py`, and this document together.
2. Once the workflow is on the repository's default branch, open GitHub
   **Actions → EcoSIM Ubuntu 24.04 binary → Run workflow** and select the
   branch to build. The selected branch must contain these files too.
   Pushing a `v*` tag containing the workflow also triggers a build.
3. After success, download the `ecosim-ubuntu24.04-x86_64-<commit>` artifact
   from the run. It is retained for 30 days, subject to repository limits.
4. Unzip the artifact, then verify and extract the enclosed archive:

```bash
sha256sum -c ecosim-ubuntu24.04-x86_64.tar.gz.sha256
tar -xzf ecosim-ubuntu24.04-x86_64.tar.gz
cd ecosim-ubuntu24.04-x86_64
sha256sum -c SHA256SUMS
```

The tar archive preserves executable permissions. No GitHub release is
published automatically. Builds use committed source and pinned submodule
commits; uncommitted workstation changes are not included.

## Install runtime dependencies and run

The package contains `bin/ecosim.f90.x`, the standalone double-precision
Release executable. EcoSIM and its scientific dependencies use the default
static-library build, but the executable still depends on Ubuntu runtime
libraries, including the Fortran runtime. This is not a fully static binary.

The packaging step derives the runtime package list from the executable's
resolved dependencies. On the target Ubuntu 24.04 machine, from the
extracted directory, run:

```bash
sudo apt-get update
xargs -r -a runtime-packages.txt sudo apt-get install --no-install-recommends -y
ldd bin/ecosim.f90.x
```

Run EcoSIM from the intended case working directory, using the extracted
binary's absolute path and the case's namelist as its first argument:

```bash
cd /path/to/your/case
/path/to/ecosim-ubuntu24.04-x86_64/bin/ecosim.f90.x your_case.nml
```

Input datasets, forcing, namelists, and restart files are not in the binary
archive. Preserve their relative paths and use a scratch case when testing
to avoid overwriting existing outputs. The main executable does not implement
a `--help` or `--version` option; those strings would be treated as input paths.

## What the workflow verifies

- The runner is x86_64 with glibc 2.39, and the installed program is ELF64 x86-64.
- Required GLIBC symbol versions in the binary are no newer than 2.39.
- Build-directory RPATHs are removed, and shared dependencies resolve without
  `LD_LIBRARY_PATH`; libraries without an Ubuntu package owner are rejected.
- A fresh Ubuntu 24.04 container installs only the recorded runtime packages
  and runs `ldd -r` to check library resolution and relocation symbols.

The container check uses currently available Ubuntu 24.04 packages; it is not
an exact recreation of the target's original security-update state. It is a
loader check, not an EcoSIM simulation or proof of scientific correctness.
No claim of support is made for an older Ubuntu/glibc release or ARM hardware.

`build-info.txt` records source and submodule commits, OS/compiler versions,
and resolved runtime package versions. `runtime-libraries.txt` and
`elf-version-info.txt` give dependency details. Failed builds upload diagnostics.

GitHub references: [runner labels](https://docs.github.com/en/actions/reference/runners/github-hosted-runners),
[manual workflow runs](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow),
and [artifact uploads](https://github.com/actions/upload-artifact).
