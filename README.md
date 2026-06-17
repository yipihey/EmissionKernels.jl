# EmissionKernels.jl

Per-channel and per-line **radiative emissivity** for a primordial + metal astrophysical
plasma — the foundation layer beneath [`ChemistryKernels.jl`](https://github.com/yipihey/ChemistryKernels.jl).
One source of truth for two uses:

- **Cooling** — `cooling_rate_total` sums every radiative channel (H/He collisional
  excitation & ionisation, recombination, bremsstrahlung, H₂, HD, CMB-Compton, and metal
  fine-structure) into the total volumetric cooling rate `[erg s⁻¹ cm⁻³]` that a chemistry
  network needs. `ChemistryKernels.cooling_edot = -cooling_rate_total`.
- **Synthetic emission** — per-channel (`emiss_*`) and per-line (`lya_emissivity`,
  `metal_line_emissivities`) volumetric emissivities for mock spectra and line maps.

Every coefficient is a pure, `@inline`, precision-generic (`Float64`/`Float32`),
allocation-free scalar function, runnable on CPU and on GPU via KernelAbstractions
(Metal/CUDA package extensions). Physics: Abel/Anninos et al. (1997) primordial cooling;
Galli & Palla (2008) H₂; Glover & Jappsen (2007) metal fine-structure.

> Most users want **ChemistryKernels.jl** (chemistry + cooling); pull this in directly only
> for the emission/cooling‑by‑channel surface. It carries the radiative physics that
> ChemistryKernels depends on.

## Install (custom registry — not in Julia's General registry)

```julia
pkg> registry add https://github.com/yipihey/VespaRegistry
pkg> add EmissionKernels
```

## Quick look

```julia
using EmissionKernels
# total radiative cooling [erg/s/cm³] at a thermodynamic state (densities in cm⁻³):
Λ = cooling_rate_total(nHI, nHII, nHeI, nde, nH2, nHD, T, z)
# per-channel / per-line for synthetic observations:
ch    = radiative_channels(nHI, nHII, nHeI, nde, nH2, nHD, T, z; nH=nH, metals=metal_abund(solar=1.0))
lines = metal_line_emissivities(T, z, nHI, nHII, nde, nH2, nH, metal_abund(solar=1.0))  # [C II] 158µm, [O I] 63µm, …
```

See `ChemistryKernels.jl` for the integration guide and the public-API reference.

## License

University of Illinois/NCSA Open Source License (the Enzo Public License) — see `LICENSE`.
Extracted from the Vespa/EnzoNG project (`github.com/yipihey/enzo-dev`).
