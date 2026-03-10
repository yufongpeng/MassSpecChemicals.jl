module MassSpecChemicals

using Combinatorics, TypedTables, MLStyle, Statistics, StatsBase, Dictionaries, Intervals, SplitApplyCombine, Plots
using UnitfulMoles: parse_compound, ustrip, @u_str
using SentinelArrays: ChainedVector
using Bessels: besseli, besseli!
import Base: show, length, +, -, *, /, isless, isequal, in, union, intersect, iterate, Broadcast.broadcastable

export 
    # Types
    AbstractChemical, Chemical, FormulaChemical, ChemicalSeries, Isobars, Isotopomers, ChemicalLoss, ChemicalPair,
    AbstractAdductIon, AdductIon,
    AbstractAdduct, Adduct, ComposedAdduct, 
    LossElectron,
    Protonation,
    ProtonationNLH2O,
    ProtonationNL2H2O,
    ProtonationNL3H2O,
    DiProtonation,
    TriProtonation,
    AddNH4,
    AddNH4Protonation,
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

    # Setter 
    set_adduct!, set_addabbr!, set_element!, 

    # Parser
    parse_chemical, parse_adduct, 

    # Attributes
    chemicalformula, chemicalelements, chemicalname, chemicalabbr, chemicalsmiles, ioncore, ionadduct, kmer, charge, ncharge, retentiontime,
    mmi, molarmass, mz, 
    chemicalentity, chemicalspecies, chemicalpair, 
    chemicalparent, isotopomersisotopes,
    analyzedchemical, detectedchemical, detectedelements, detectedcharge, detectedisotopes,
    inputchemical, outputchemical,
    analyzedprecursor, detectedproduct, detectedproductelements, detectedproductcharge, detectedproductisotopes, 
    inputprecursor, outputproduct, 
    msstage, chemicaltransitions, 
    getchemicalproperty, 

    # Derived attributes
    isotopicabundance, Isotopologues, TandemIsotopologues, 
    isobars_rt, isobars_rt_table, CoelutingIsobars, 
    
    # Spectrum 
    Spectrum, MSScan, AllIons, TargetIon, Fragmentation, peak_table, 
    MSAnalyzer, Quadrupole, QuadrupoleIonTrap, QIT, LinearIonTrap, LIT, TimeOfFlight, TOF, Orbitrap, FourierTransformIonCyclotronResonance, FTICR, 
    resolving_power, 
    
    # Visualization
    plot_spectrum, plot_spectrum!, 
    plot_resolving_power, plot_resolving_power!,
    plot_window, plot_window!,

    # Utils
    match_chemical, 
    ischemicalequal, isadductequal, 
    acrit, rcrit, crit, @ri_str, makecrit_delta, makecrit_value

"""
    AbstractChemical

Abstract type for chemicals. 
    
The attribute `chemicalname -> String` (unique chemical name) is required for a concrete type of `AbstractChemical`. 
It defaults to property `name`.

At least one of the following attributes are required. They are interchangable.
* `chemicalformula -> String`: chemical formula. It defaults to  property `:formula`.
* `chemicalelements -> Vector{Pair{String, Int}}`: chemical elements. It defaults to property `:elements`.

The following attributes are optional, but generic functions are defined.
* `chemicalabbr -> String`: abbreviation. It defaults to property `:abbreviation` and `chemicalname`. 
* `chemicalsmiles -> String`: SMILES. It defaults to property `:SMILES` and `""`.
* `charge -> Int`: charge state; positive for cation, negative for anion. It defaults to property `:charge` annd `0`.
* `chemicalentity -> AbstractChemical`: chemical entity (representative). 
* `chemicalspecies -> Vector{<: AbstractChemical}`: chemical species. 
* `chemicalpair -> Pair{<: AbstractChemical, <: AbstractChemical}`: chemical pair. 
* `rt -> AbstractFloat`: retention time. It defaults to property `:rt` and `NaN`
* `ncharge -> Int`: number of charges. 
* `mmi -> AbstractFloat`: monoisotopic mass.
* `molarmass -> AbstractFloat`: molar mass.
* `mz -> AbstractFloat`: mass to charge ratio.
"""
abstract type AbstractChemical end
# mt, ccs 
include(joinpath("type", "adduct.jl"))
include(joinpath("type", "chemical.jl"))
include(joinpath("type", "adduction.jl"))
include(joinpath("type", "utils.jl"))
include(joinpath("type", "spectrum.jl"))
include("attr.jl")
include(joinpath("interface", "adduct.jl"))
include(joinpath("interface", "chemical.jl"))
include(joinpath("interface", "adduction.jl"))
include(joinpath("interface", "utils.jl"))
include(joinpath("interface", "spectrum.jl"))
include("elements.jl")
include("formula.jl")
include("mass.jl")
include("isotopologues.jl")
include("spectrum.jl")
include("rt.jl")
include("isobars.jl")
include("io.jl")
include("utils.jl")

end