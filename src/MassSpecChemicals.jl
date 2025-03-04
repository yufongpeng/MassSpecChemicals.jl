module MassSpecChemicals

using Combinatorics, TypedTables, MLStyle, Statistics, StatsBase
using UnitfulMoles: parse_compound, ustrip, @u_str
using SentinelArrays: ChainedVector
import Base: show, length, +, -, *, /, isless, isequal, in, union, intersect, Broadcast.broadcastable

export Ion, Chemical, Isobars,
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

    parse_chemical, parse_adduct, chemicalformula, chemicalname, chemicalabbr, chemicalsmiles, ioncore, ionadduct, kmer, charge, ncharge, rt,
    getchemicalattr, isionequal, ischemicalequal, isadductequal, 

    push_adduct_abbr!, push_adduct_name!, ADDUCT_ABBR, ADDUCT_NAME,

    chemicalvariants, ionvariants,

    mw, mz, abundantchemical,

    isotopicabundance, isotopetable, getisotopes,

    isobartable_rt, getisobars_rt,

    acrit, rcrit, crit, @ri_str

"""
    AbstractChemical

Abstract type for chemicals. 
    
The following attributes are required for a concrete type of `AbstractChemical`
* `name`: a unique chemical name (`String`).
* `formula`: chemical formula (`String`).

The following attributes are optional, but has been implemented for `AbstractChemical`
* `abbreviation`: common abbreviation (`String`); defaults to `:name`. 
* `SMILES`: SMILES (`String`)

These attributes are not neccessary be type fields, but `getchemicalattr(chemical, Val(attr)` should return the correct value. See `getchemicalattr` for details.
"""
abstract type AbstractChemical end
# mt, ccs 
include(joinpath("type", "adduct.jl"))
include(joinpath("type", "chemical.jl"))
include(joinpath("type", "ion.jl"))
include(joinpath("type", "utils.jl"))
include(joinpath("interface", "adduct.jl"))
include(joinpath("interface", "chemical.jl"))
include(joinpath("interface", "ion.jl"))
include(joinpath("interface", "utils.jl"))
include("attr.jl")
include("elements.jl")
include("formula.jl")
include("mw.jl")
include("isotopes.jl")
include("rt.jl")
include("io.jl")
include("utils.jl")

end