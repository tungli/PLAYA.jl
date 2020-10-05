using PLAYA
using Statistics
using UnicodeFun # to_latex
using DifferentialEquations
using DiffEqSensitivity
using ODEInterfaceDiffEq # radau
using Catalyst

# parse BOLOS output
bolos_interpolations = parse_bolos_output("assets/air-kinetics-bolos.dat")

EN = 200. # Townsend
Tgas = 300.
Te = bolos_interpolations["Mean energy"](EN) * 1.1604505e4  # Temperature of electrons in K
dTion  = 2.0e0 / ( 3.0e0 * 1.3807e-16 ) * 1.6605e-24 * ( 1.0e-17 * EN )^2
TionN  = Tgas + dTion * 14.0e0 * 8.0e19^2
TionN2 = Tgas + dTion * 28.0e0 * 4.1e19^2
TionN3 = Tgas + dTion * 42.0e0 * 6.1e19^2
TionN4 = Tgas + dTion * 56.0e0 * 7.0e19^2
TeffN  = ( TionN  + 0.5e0 * Tgas ) / ( 1.0e0 + 0.5e0 )
TeffN2 = ( TionN2 + 1.0e0 * Tgas ) / ( 1.0e0 + 1.0e0 )
TeffN3 = ( TionN3 + 1.5e0 * Tgas ) / ( 1.0e0 + 1.5e0 )
TeffN4 = ( TionN4 + 2.0e0 * Tgas ) / ( 1.0e0 + 2.0e0 )
Teff3Q = 46.501347e-27*((1.87e-4*(1.0e5/(273.0*1.38064e-23))*( 1.0e-21 * EN ))^2)/(3.0e0*1.38064e-23) + Tgas


globals = Dict{Symbol, Any}([
                :E => EN,
                :Te => Te,
                :Tgas => Tgas,
                :dTion   => dTion,
                :TionN   => TionN,
                :TionN2  => TionN2,
                :TionN3  => TionN3,
                :TionN4  => TionN4,
                :TeffN   => TeffN,
                :TeffN2  => TeffN2,
                :TeffN3  => TeffN3,
                :TeffN4  => TeffN4,
                :Teff3Q  => Teff3Q,
               ])

# create functions that will get evaluated to rates
globals[:BOLOS] = (id, env) -> bolos_interpolations[id](env[:E])

println(globals)

# replace stuff so it is Catalyst-compliant
replacements = [
                "^+" => to_latex("^+"),
                "^-" => to_latex("^-"),
                "(" => to_latex("\\llcorner"),
                ")" => to_latex("\\lrcorner"),
                "4.5eV" => to_latex("4°5eV"),
                "`" => to_latex("\\tilde"),
                "=>" => "-->",
    ]

# parse the reaction file
(rs, ps) = parse_reactions("assets/air-kinetics.reactions", globals; to_replace=replacements)

# initial conditions
n_gas = 2.5e19
conc = Dict(["e" => 1e12, "N2" => 0.8*n_gas, "O2" => 0.2*n_gas])
u0 = [ get(conc, string(s), 0.) for s in species(rs) ]
tspan = (0., 1e-4)

print("Converting to ODEs....")
odesys = convert(ODESystem, rs)

# Keep electrons and neutrals fixed
keep_fixed!(odesys, "e")
keep_fixed!(odesys, "ANY_NEUTRAL")
