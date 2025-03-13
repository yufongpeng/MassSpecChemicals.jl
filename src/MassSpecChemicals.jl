module MassSpecChemicals

using Combinatorics, TypedTables, MLStyle, Statistics, StatsBase, Dictionaries
using UnitfulMoles: parse_compound, ustrip, @u_str
using SentinelArrays: ChainedVector
import Base: show, length, +, -, *, /, isless, isequal, in, union, intersect, iterate, Broadcast.broadcastable

export AdductIon, Chemical, Isobars, Isotopomers, 
    AbstractPosAdduct, AbstractNegAdduct,
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
    SodizationAddNH4, 
    Potassiation, 
    PotassiationProtonation, 
    DiPotassiation, 
    Lithiation, 
    LithiationProtonation, 
    DiLithiation, 
    Silveration, 

    AddElectron,
    Deprotonation,
    DeprotonationNLH2O,
    DiDeprotonation,
    TriDeprotonation,
    AddOAc,
    AddOFo,

    # LossCH2O,
    # AddO,
    # AddC2H2O,
    # LossCH8NO,
    # LossC2H8NO,
    # AddC3H5NO,
    # AddC2H5NO,
    # LossCH3,
    # DeprotonationLossSerineAddH2O,
    Fluoridation,
    Chloridation,
    Demethylation,

    parse_chemical, parse_adduct, chemicalformula, chemicalelements, chemicalname, chemicalabbr, chemicalsmiles, ioncore, ionadduct, kmer, charge, ncharge, rt,
    getchemicalattr, ischemicalequal, isadductequal, 

    set_adduct_abbr!, set_adduct_name!, ADDUCTS, ELEMENTS, set_elements!, 

    mw, mz, abundantchemical,

    isotopicabundance, isotopologues_table, isotopologues,

    isobars_rt, isobars_rt_table,

    acrit, rcrit, crit, @ri_str, real_interval

"""
    AbstractChemical

Abstract type for chemicals. 
    
The following attributes are required for a concrete type of `AbstractChemical`
* `name`: a unique chemical name (`String`).

At least one of the following attributes are required. They are interchangable.
* `formula`: chemical formula (`String`).
* `elements`: chemical elements (`Vector{Pair{String, Int}}`).

The following attributes are optional, but user interfaces are implemented 
* `abbreviation`: common abbreviation (`String`); defaults to `:name`. 
* `SMILES`: SMILES (`String`); defaults to `""`.
* `charge`: charge of `chemical` (positive or negative). 
* `abundant_chemical`: the most abundant chemical from a chemical (itself) or isobars. 
* `rt`: retention time. 

These attributes are not neccessary to be type fields, but `getchemicalattr(chemical, Val(attr); kwargs...)` should return the correct value. See `getchemicalattr` for details.
"""
abstract type AbstractChemical end
# mt, ccs 
include(joinpath("type", "adduct.jl"))
include(joinpath("type", "chemical.jl"))
include(joinpath("type", "adduction.jl"))
include(joinpath("type", "utils.jl"))
include(joinpath("interface", "adduct.jl"))
include(joinpath("interface", "chemical.jl"))
include(joinpath("interface", "adduction.jl"))
include(joinpath("interface", "utils.jl"))
include("attr.jl")
include("elements.jl")
include("formula.jl")
include("mw.jl")
include("isotopologues.jl")
include("rt.jl")
include("io.jl")
include("utils.jl")

end