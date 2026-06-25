# cooling_dust.jl — grain-mediated thermal channels of the gas energy budget.
#
# Three volumetric rates that cross the gas<>dust interface:
#   Gamma_PE   — photoelectric heating of GAS by UV-illuminated PAHs/grains [erg/s/H]
#   Lambda_gr  — gas-grain collisional coupling [erg/cm3/s]; positive = gas cools
#   Lambda_dust — dust continuum emission [erg/cm3/s] (diagnostic; NOT direct gas heat)
#
# Pattern mirrors cooling_metal.jl: pure @inline, precision-generic R = typeof(first_arg),
# allocation-free scalars; CGS in/out.  Imported and re-exported by ChemistryKernels so
# its cooling_edot / evolve_cell can call them without an explicit using EmissionKernels.
#
# References:
#   Bakes & Tielens   (1994)  ApJ  427  822    — photoelectric heating efficiency
#   Hollenbach & McKee (1989)  ApJ  342  306    — gas-grain collisional coupling
#   Hollenbach & McKee (1979)  ApJS  41  555    — dust modified-blackbody emission

export Gamma_PE, Lambda_gr, Lambda_dust

# Photoelectric heating rate per H nucleus [erg/s]
# Bakes & Tielens (1994), Eq. 2.  UV photons eject electrons from PAHs/small
# grains; the heating efficiency e depends on the grain charging parameter
# psi = G0*sqrt(T_gas) / n_e.  Multiply by n_H to get volumetric rate [erg/cm3/s].
# Z_rel scales the grain abundance (proportional to metallicity); Gamma_PE -> 0 in pristine gas.
@inline function Gamma_PE(T_gas::Real, G0::Real, Z_rel::Real, n_e::Real)
    R = typeof(T_gas)
    psi = G0 * sqrt(T_gas) / max(n_e, R(1e-20))
    e1  = R(4.87e-2) / (one(R) + R(4e-3) * psi^R(0.73))
    e2  = R(3.65e-2) * (T_gas / R(1e4))^R(0.7) / (one(R) + R(2e-4) * psi)
    return R(1.3e-24) * (e1 + e2) * G0 * Z_rel
end

# Gas-grain collisional coupling [erg/cm3/s]
# Hollenbach & McKee (1989), Lambda_gr per unit volume.  Positive when gas is hotter
# than dust (gas cools); negative when T_gas < T_dust (gas gains energy from
# contact with warmer grains, e.g. in strongly irradiated PDR surfaces).
@inline function Lambda_gr(T_gas::Real, T_dust::Real, n_H::Real, Z_rel::Real)
    R = typeof(T_gas)
    return R(2e-33) * sqrt(T_gas) * (T_gas - T_dust) * n_H^2 * Z_rel
end

# Dust continuum emission [erg/cm3/s]
# Modified blackbody with kappa ~ nu^2: Lambda_dust ~ n_H * Z_rel * T_dust^6
# (Hollenbach & McKee 1979, Krumholz et al. 2011).  Energy radiated by dust grains
# to the radiation field; does NOT enter the GAS energy equation directly (gas
# couples to dust only through Lambda_gr).  Exported for diagnostics and post-processing.
@inline function Lambda_dust(T_dust::Real, Z_rel::Real, n_H::Real)
    R = typeof(T_dust)
    return R(2.0e-27) * n_H * Z_rel * T_dust^R(6)
end
