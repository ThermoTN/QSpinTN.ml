---
name: qspintn-units
description: Interpret and convert QSpinTN physical thermal outputs, especially molar SI susceptibility versus CGS susceptibility.
---

# Thermal units

All source paths below are relative to the QSpinTN.ml repository root. The machine-readable contract is [thermal_units.json](../contracts/thermal_units.json), identifier `qspintn_chi_si_cm3_v1`.

## Scope

This contract covers the physical thermal outputs `T_l`, `Cm_l`, and `Chi_l` produced by `ManyBodySolver/GetResult.m` using `UtilityFunc/GetUnitCon.m` and `UtilityFunc/GetField.m`. A direct solver result is normally in natural units and needs its own conversion. In particular the standalone iLTRG example, DMRG, and real-time TDVP do not acquire this contract merely by sharing the repository.

| Output/input | Convention |
| --- | --- |
| `Model.ES` | Positive energy scale expressed in kelvin |
| Natural `beta` | `T_l = ES/beta`, in K |
| Applied `Field.B` | Magnetic flux density vector in tesla |
| `Cm_l` | J/(mol K), per mole of magnetic spins |
| `Chi_l` | **cm^3/mol (SI)**, per mole of magnetic spins |

`qspintn_native_chi` is a consumer identifier for this `Chi_l` convention, not a MATLAB field name. It is distinct from the unconverted ratio `M/h`.

## Why susceptibility is SI cm^3/mol

With `kB = Const.kb*Const.MJ` in joules per kelvin, the code uses

```text
h_con = ES*kB/muB
h = |B|/h_con
Chi_con = NA*mu0*muB^2/(ES*kB) * 10^6
Chi_l = (M/h)*Chi_con = 10^6*mu0*NA*muB*M/|B|
```

`M` is the per-spin moment coefficient along the field; the model's moment operator already contains its g factor. Do not multiply a physical `Chi_l` curve by g squared again.

Conversions of the numerical `Chi_l` values are:

- SI m^3/mol: multiply by `1e-6`.
- CGS cm^3/mol (emu/mol/Oe): divide by `4*pi`.
- To compare per mole of formula units, apply the independently established number of magnetic spins per formula unit to both Cm and Chi. Mass-normalized data also need the material's molar mass and normalization convention.

The factor `4*pi` is essential: the familiar CGS Curie constant for S=1/2 and g=2 is approximately 0.375 cm^3 K/mol, whereas this output convention gives approximately 4.71 cm^3 K/mol. Package constants are rounded, so use tolerances consistent with those constants.

## Physical limits

`GetResult(..., 'Chi')` computes finite-field `M/H`, with applied `H=|B|/mu0`. It does not calculate the zero-field differential susceptibility or a demagnetization correction. Check field linearity if comparing to a low-field susceptibility measurement. The shared `GetField` divides by `norm(B)` to construct a direction; a literal zero field needs route-specific handling and is not a valid finite-field susceptibility request.

Report the SI/CGS label, per-spin/formula-unit normalization, applied field, and finite-field definition together with curves. A dimensional conversion cannot correct phonon backgrounds, impurity terms, offsets, or model mismatch.
