module MassSpecChemicals

using Combinatorics, TypedTables, MLStyle, Statistics, StatsBase, Dictionaries, Intervals, SplitApplyCombine, Plots
using UnitfulMoles: parse_compound, ustrip, @u_str
using SentinelArrays: ChainedVector
using Bessels: besseli, besseli!
import Base: show, length, +, -, *, /, isless, isequal, in, union, intersect, iterate, Broadcast.broadcastable

export 
    # Types
    AbstractChemical, Chemical, FormulaChemical, ChemicalTransition, ChemicalSeries, Isobars, Isotopomers, ChemicalLoss, ChemicalGain, 
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
    seriesanalyzedchemical, seriesanalyzedisotopes, seriesanalyzedcharge, seriesanalyzedelements, 
    msstage, chemicaltransition, isotopomerstate, 
    getchemicalproperty, 

    # Derived attributes
    isotopicabundance, Isotopologues, TandemIsotopologues, 
    CoelutingIsobars, isobar_table, group_isotopologues, 
    
    # Spectrum 
    Spectrum, MSScan, AllIons, Isolation, SelectedIonMonitor, Fragmentation, peak_table, 
    MSAnalyzer, Quadrupole, QuadrupoleIonTrap, QIT, LinearIonTrap, LIT, TimeOfFlight, TOF, Orbitrap, FourierTransformIonCyclotronResonance, FTICR, 
    resolving_power, 
    
    # Visualization
    plot_spectrum, plot_spectrum!, 
    plot_resolving_power, plot_resolving_power!,
    plot_window, plot_window!,

    # Utils
    match_chemical, 
    ischemicalequal, isadductequal, 
    acrit, rcrit, crit, @ri_str, makecrit_delta, makecrit_value, 
    value_error, relative_error, percentage_error, ppm_error, relative_error_mean, percentage_error_mean, ppm_error_mean

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
* `ncharge -> Int`: number of charges
* `retentiontime -> AbstractFloat`: retention time. It defaults to property `:rt` and `NaN`.
* `mmi -> AbstractFloat`: monoisotopic mass
* `molarmass -> AbstractFloat`: molar mass
* `mz -> AbstractFloat`: mass to charge ratio
* `chemicalparent -> AbstractChemical`: parent chemical without delocalized isotopes replacement
* `isotopomersisotopes -> Vector{Pair{String, Int}}`: delocalized isotopes replacement of isotopomers
* `isotopomerstate -> Int`: isotopomers state, i.e. equivalent number of isotope
* `chemicalentity -> AbstractChemical`: a single chemical entity representing the chemical
* `chemicalspecies -> Vector{<: AbstractChemical}`: multiple chemical entities having shared properties. 
* `chemicaltransition -> Vector{<: AbstractChemical}`: chemical entities analyzed in each stage of instrumental analysis.  
* `inputchemical -> AbstractChemical`: a single chemical entity that is the input of the very ending of instrumental analysis. 
* `outputchemical -> AbstractChemical`: a single chemical entity that is the output of the very ending of instrumental analysis. 
* `analyzedchemical -> AbstractChemical`: a single chemical entity directly detected in the very begining of instrumental analysis. 
* `detectedchemical -> AbstractChemical`: a single chemical entity directly detected in the very ending of instrumental analysis. 
* `seriesanalyzedchemical -> Vector{<: AbstractChemical}`: chemical entities directly analyzed in each stage of instrumental analysis.
* `detectedcharge -> Int`: charge of detected chemical
* `detectedelements -> Vector{Pair{String, Int}}`: elements of detected chemical
* `detectedisotopes -> Vector{Pair{String, Int}}`: isotopes replacement of detected chemical
* `seriesanalyzedcharge -> Vector{Int}`: charge of sereially analyzed chemicals
* `seriesanalyzedelements -> Vector{Vector{Pair{String, Int}}}`: elements of sereially analyzed chemicals
* `seriesanalyzedisotopes -> Vector{Vector{Pair{String, Int}}}`: isotopes replacements of sereially analyzed chemicals
* `msstage -> Int`: number of stages of MS the chemical has been through.

Specific Methods for the attributes are defined for other intrinsic chemical type on different chemical level
* Entity Level: attribute of the corresponding chemical entity
* Species Level: attribute of the corresponding chemical species
* Transition Level: attribute of the corresponding chemical transition
"""
abstract type AbstractChemical end
# mt, ccs 
include(joinpath("type", "adduct.jl"))
include(joinpath("type", "chemical.jl"))
include(joinpath("type", "adduction.jl"))
include(joinpath("type", "utils.jl"))
include(joinpath("type", "spectrum.jl"))
include("interface.jl")
include("attr.jl")
include("attr_adduct.jl")
include("attr_measure.jl")
include("chemical.jl")
include("elements.jl")
include("adduct.jl")
include("formula_elements.jl")
include("measure.jl")
include("isotopologues.jl")
include("spectrum.jl")
include("isobars.jl")
include("plot.jl")
include("io.jl")
include("utils.jl")

end