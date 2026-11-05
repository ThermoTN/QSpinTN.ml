# Agent method interface

This directory provides reusable QSpinTN method guidance, source-level interface limitations, and a thermal-unit contract. It is maintained with the numerical package so consumers can read documentation that belongs to the same Git revision. The initial skill extraction and documentation were developed with GPT assistance.

Start with [method selection](qspintn-method-selection/SKILL.md) and [units](qspintn-units/SKILL.md), then load only the selected method skill. Source paths written in the skills are relative to the **QSpinTN.ml repository root**, not to a skill directory. Markdown links between skill documents are relative links.

For new interactions or lattice definitions, read [model building](qspintn-model-building/SKILL.md). It explains the native model/map contract, constrained Hamiltonian presentation, energy scaling, finite-bond counting and validation before fitting, and distinguishes native models from a consumer's declarative adapter.

The [units skill](qspintn-units/SKILL.md) covers Cm, susceptibility and M–H data: natural versus physical solver outputs, J/mJ, tesla/oersted, SI/CGS susceptibility, and μB/spin versus literal μB/mol, emu/mol and A·m²/mol. It includes conversions in both directions, molar and mass normalization, and numerical examples. The existing thermal-unit JSON contract remains unchanged; documenting magnetization conversion does not add an M–H fitting adapter. An application's default experimental normalization must be recorded separately from the package's output convention.

## Discovery contract

[manifest.json](manifest.json) uses schema `qspintn.skills.v1`. Each `skills` entry supplies `name`, `description`, `path`, and `method`; paths resolve relative to this directory. `method` is the package solver identifier (`ED`, `iLTRG`, `XTRG`, `tanTRG`, `DMRG`, or `TDVP`), or null for shared guidance. `contracts.thermal_units` points to [the thermal-unit contract](contracts/thermal_units.json), also relative to this directory.

Consumers should validate the schema and referenced files, record the package Git revision and content hashes, and load only explicitly supported skills. Discovery is documentation access, not permission to execute code, a tool registration mechanism, or a claim that the consumer implements every solver. Consumers own orchestration, fitting optimizers, data handling, resource limits, and execution adapters.

## Maintaining the interface

When changing a MATLAB entry point, parameter, normalization, output grid, or known limitation, update the corresponding skill in the same branch. Keep each skill's YAML name/description consistent with its manifest entry. Revise the unit contract identifier when its physical meaning changes; change schema versions for incompatible structural changes. Preserve the explicit distinction between direct natural-unit results and converted physical outputs.

The present skills document the checked-in interfaces and known restrictions; they do not certify numerical execution on every platform. Examples may require prepared inputs, MATLAB toolboxes, platform-specific binaries, or local saving hooks. Keep numerical validation evidence separate from documentation/schema validation.
