module MassSpecChemicals

using Combinatorics, TypedTables, MLStyle
using UnitfulMoles: parse_compound, ustrip, @u_str
using SentinelArrays: ChainedVector
import Base: show, length, +, -, *, /, isless, isequal, in, union, intersect

export Ion, Chemical, IonCluster,
    PosAdduct, NegAdduct,
    LossElectron,
    Protonation,
    ProtonationNLH2O,
    ProtonationNL2H2O,
    ProtonationNL3H2O,
    DiProtonation,
    TriProtonation,
    AddNH4,
    AddHNH4,
    Add2NH4,
    Sodization,
    SodizationProtonation,
    DiSodization,

    AddElectron,
    Deprotonation,
    DeprotonationNLH2O,
    DiDeprotonation,
    TriDeprotonation,
    AddOAc,
    AddHCOO,

    LossCH2O,
    AddO,
    AddC2H2O,
    LossCH8NO,
    LossC2H8NO,
    AddC3H5NO,
    AddC2H5NO,
    LossCH3,
    DeprotonationLossSerineAddH2O,

    chemicalformula, chemicalname, ioncore, ionadduct, kmer, charge, rt,
    infonames, getinfo, infopairs,
    transform_chemicals, transform_ions,
    mw, mz, abundantion,
    abundance, isotope_abundance, getisotopes,
    isobar_rt_mz1, getisobars,
    acrit, rcrit, crit, @ri_str

"""
    AbstractChemical

Abstract type for chemicals.
"""
abstract type AbstractChemical end

include(joinpath("type", "adduct.jl"))
include(joinpath("type", "chemical.jl"))
include(joinpath("type", "ion.jl"))
include(joinpath("type", "utils.jl"))
include(joinpath("interface", "adduct.jl"))
include(joinpath("interface", "chemical.jl"))
include(joinpath("interface", "ion.jl"))
include(joinpath("interface", "utils.jl"))
include("elements.jl")
include("formula.jl")
include("mw.jl")
include("isotopes.jl")
include("rt.jl")
include("io.jl")
include("utils.jl")
include(joinpath("biochemical", "BioChemicals.jl"))

end