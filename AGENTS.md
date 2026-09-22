# EcoSIM repository guidance

## Scope and purpose

This repository contains the EcoSIM biogeochemical model and a biocrust simulation
workspace. Use these instructions throughout the repository; follow more specific
`AGENTS.md` files in subdirectories where applicable. Explicit user instructions
take precedence over repository guidance.

## Repository map

- `f90src/`: Fortran model source. `Main/` and `Ecosim_mods/` coordinate model
  execution; `Modelforc/` handles forcing; `IOutils/` handles inputs, history,
  and restarts.
- `f90src/Plant_bgc/`, `Microbial_bgc/`, `Geochem/`, `HydroTherm/`, and
  `Transport/`: biological, chemical, hydrological, thermal, and transport processes.
- `f90src/Ecosim_datatype/`, `APIData/`, `Modelconfig/`, and `Modelpars/`:
  state/flux definitions, API data, identifiers, and parameters.
- `f90src/Balances/` and `ModelDiags/`: accounting and diagnostics.
- `drivers/`, `f90src/APIs/`, and `f90src/ATSUtils/`: drivers and coupling interfaces.
- `input_data/`: shared model inputs and PFT parameter NetCDF/CDL files.
- `examples/inputs/`: case-specific inputs. `examples/run_dir/`: case namelists
  and simulation working directories, including `biocrust/`.
- `python_tools/`: a separate Git submodule for preprocessing and analysis.
  Read its `AGENTS.md` before working there. Its canonical skills are under
  `python_tools/.agents/skills/`; generated preprocessing outputs belong under
  `python_tools/result/` unless the user specifies another destination.
- `tests/` and `regression-tests/`: validation resources and regression tooling.
- `build/` and `local/`: generated build/install products.
- `3rd-partylibs/`: third-party dependency submodules.

## Knowledge documents and user preferences

- Save generated EcoSIM process explanations, code walkthroughs, and similar
  knowledge documents (including HTML guides) under `ecosim_knowledge/` by
  default, unless the user specifies another destination. This is a persistent
  user preference for future requests.
- Keep relative source links valid from that directory when creating or moving
  a document.

## Working with model code

- Inspect `git status` before editing. Preserve existing user changes, untracked
  inputs, simulation outputs, and restart files; do not clean or reset them.
- Trace callers and CMake source lists to establish the active implementation.
  Similar routines occur in initialization, process, and coupling modules; do not
  assume one is obsolete or change every copy without checking its role.
- Follow nearby Fortran conventions, precision kinds, array ordering, and module
  interfaces. Update allocation, initialization, API mapping, history, and restart
  handling as needed when adding persistent state or flux variables.
- Keep edits focused. Do not hand-edit generated build files or modify dependency
  submodules unless that is part of the requested work.
- Treat `python_tools` changes as submodule changes, separate from root-repository
  changes. Do not automatically commit, update submodule pointers, or publish work.

## Scientific and input conventions

- When asked to evaluate a plant-trait `.desc` file (including
  `*.plant_trait.desc.*`), read `python_tools/AGENTS.md` and use the relevant
  skills under `python_tools/.agents/skills/`, starting with
  `ecosim-plant-trait-sanity-check/SKILL.md`. Treat the `.desc` file as read-only.
- Track units explicitly: distinguish carbon mass from dry mass, per-individual
  quantities from per-area quantities, and per-area values from grid-cell totals.
  Check the actual equations when comments and metadata disagree.
- Preserve C/N/P, water, and energy accounting for affected processes. Explain
  newly introduced sources, sinks, transfers, or numerical clipping.
- Read enum values from `f90src/Modelconfig/ElmIDMod.F90`; do not rely solely on
  NetCDF flag descriptions. Validate PFT codes against the parameter file selected
  by the case, and check whether grid Koppen settings replace the PFT suffix.
- EcoSIM treats lichen and moss through nonvascular plant pathways. Interpret
  seed/planting terminology according to the implemented pools and equations;
  do not equate establishment reserves with an initialized mature community.
- Planting strings use `DDMMYYYY`, population density in units per square meter,
  and depth in meters. Honor user-provided values without silently substituting
  biological defaults. For a single grid cell, use
  `NH1 = NV1 = NH2 = NV2 = 1`; `NZ` is the active PFT count.
- Use the relevant Python-tools skill for management or other preprocessing.
  Validate NetCDF dimensions, strings, units, and active entries by reading the
  generated file back. Keep editable source inputs alongside generated artifacts
  when practical. Keep paired CDL/NetCDF representations consistent when changing
  an input dataset, and identify which file the model actually reads.

## Build and validation

- Run build commands from the repository root. Inspect supported options with
  `bash build_EcoSIM.sh --help`; the standard build is
  `bash build_EcoSIM.sh`. Debug/regression configuration is available through
  `bash build_EcoSIM.sh --debug --regression_test`.
- For an existing configured build, use
  `cmake --build build/<configuration> --parallel <jobs>`, selecting an actual
  configuration directory and checking its compiler/cache first. Locate the
  resulting executable rather than assuming a fixed path from older documentation.
- `build_EcoSIM.sh --clean` deletes both `build/` and `local/`; do not use it as
  routine validation.
- Select validation appropriate to the change. Fortran changes generally need a
  build and a focused test or short simulation exercising the modified process.
  Input-only changes need schema/value checks; documentation-only changes do not
  require compiling the model.
- To screen the example cases after a Fortran change, use the
  `ecosim-smoke-run` skill (`.agents/skills/ecosim-smoke-run/`, also exposed at
  `.claude/skills/` and `.codex/skills/`). It runs each example for a few
  simulated days in a scratch mirror and checks exit status, steps completed,
  and output for NaN/Inf. A pass means the model ran and produced finite
  numbers; it is not evidence of scientific correctness.
- Inspect `regression-tests/Makefile` before running its targets. It provides
  `mtest` and `rtest`; verify the executable path and required data first. Do not
  update regression baselines merely to make a changed result pass.
- Run simulations in a scratch/copy of the case when existing outputs or restart
  pointers could be overwritten. Inspect the namelist, input paths, run duration,
  and restart mode before launching.
- Report what changed, validation actually performed, and any unresolved limits.
  Distinguish a successful build from demonstrated scientific correctness.
