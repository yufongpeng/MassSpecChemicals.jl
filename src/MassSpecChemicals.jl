module MassSpecChemicals

using Combinatorics, TypedTables, MLStyle, Statistics, StatsBase, Dictionaries, Intervals, SplitApplyCombine, Plots
using UnitfulMoles: parse_compound
using SentinelArrays: ChainedVector
import Base: show, length, +, -, *, /, isless, isequal, in, union, intersect, iterate, Broadcast.broadcastable, ==, hash, copy, axes, isempty, eltype

export 
    # Types
    AbstractChemical, AttributedChemical, ChemicalAttribute, Chemical, FormulaChemical, ChemicalTransition, ChemicalSeries, Isobars, Isotopomers,
    AbstractAdductIon, AdductIon,
    AbstractScheme, AbstractElementalScheme, AbstractStructuralScheme, StructuralChemicalScheme, AbstractCompleteScheme, 
    ElementalScheme, ChemicalGain, ChemicalLoss, 
    ChemicalSchemes, IsotopomerizedSchemes, StructuralElementalScheme, 
    CompleteSchemes, StructuralSchemes, ElementalSchemes, 

    # Default chemical for scheme 
    StandardChemical,
    Electron, 
    Proton, 
    Water, 
    Ammonia, 
    Ammonium, 
    Sodium, 
    Potassium, 
    Lithium, 
    Silver, 
    Acetate, 
    Formate, 
    AceticAcid, 
    FormicAcid, 
    MethylAcetate, 
    MethylFormate, 
    Fluoride, 
    Chloride, 
    Methenium, 
    
    # Setter 
    set_scheme!, set_schabbr!, set_element!, 

    # Elements 
    parent_element, major_isotope, minor_isotope, iselement, isisotope, 

    # Parser
    parse_chemical, parse_adduct, 
    ChemicalParser, FormulaChemicalParser, AdductParser, ChemicalExpressionParser, ChemicalTransitionParser, ChemicalEntityParser, ChemicalSchemeParser, 

    # Attributes manipulation
    getchemicalattr, getchemicalattrname, pushchemicalattr!, pushchemicalattr_new!,

    # Attributes
    chemicalformula, chemicalelements, chemicalname, chemicalabbr, chemicalsmiles, ioncore, ionadduct, ncore, charge, ncharge, retentiontime,
    mmi, molarmass, mz, 
    chemicalentity, chemicalspecies, chemicalpair, 
    chemicalparent, isotopomersisotopes, groupedisotopomersisotopes, groupedisotopomersabundance,
    analyzedchemical, detectedchemical, detectedcharge, detectedisotopes,
    inputchemical, outputchemical,
    seriesanalyzedchemical, seriesanalyzedisotopes, seriesanalyzedcharge, 
    msstage, chemicaltransition, mass_shift_index, 
    elementalscheme, structuralscheme, 

    # Scheme 
    ionize, isotopomerize, groupedisotopomerize, 

    # Isotopologues
    isotopicabundance, Isotopologues, TandemIsotopologues, group_isotopologues, 
    
    # MS analysis 
    Ionization, Spectrum, MSScan, AllIons, Isolation, SelectedIonMonitor, Fragmentation, peak_table, 
    MSAnalyzer, Quadrupole, QuadrupoleIonTrap, QIT, LinearIonTrap, LIT, TimeOfFlight, TOF, Orbitrap, FourierTransformIonCyclotronResonance, FTICR, 
    resolving_power, 

    # Co-elution
    CoelutingIsobars, isobar_table, 
    
    # Visualization
    plot_spectrum, plot_spectrum!, 
    plot_resolving_power, plot_resolving_power!,
    plot_window, plot_window!,

    # Utils
    match_chemical, ischemicalequal, 
    acrit, rcrit, crit, @ri_str, 
    value_error, relative_error, percentage_error, ppm_error, relative_error_mean, percentage_error_mean, ppm_error_mean


abstract type AbstractChemicalScheme end
@deprecate AbstractChemicalsSchema AbstractChemicalScheme # del
"""
    AbstractChemical <: AbstractChemicalScheme

Abstract type for chemicals. 
    
The attribute function [`chemicalname`](@ref) `-> String` (unique chemical name) is required for a concrete type of `AbstractChemical`. 
It defaults to attribute `name`.

At least one of the following attributes is required. They are interchangeable.
* [`chemicalformula`](@ref) `-> String`: chemical formula. It defaults to attribute `:formula`.
* [`chemicalelements`](@ref) `-> Vector{Pair{String, Int}}`: chemical elements. It defaults to attribute `:elements`.

The following attributes are optional, but generic functions are defined.
* [`chemicalabbr`](@ref) `-> String`: abbreviation. It defaults to attribute `:abbreviation` and `chemicalname`. 
* [`chemicalsmiles`](@ref) `-> String`: SMILES. It defaults to attribute `:SMILES` and `""`.
* [`charge`](@ref) `-> Int`: charge state; positive for cation, negative for anion. It defaults to attribute `:charge` and `0`.
* [`ncharge`](@ref) `-> Int`: number of charges.
* [`retentiontime`](@ref) `-> AbstractFloat`: retention time. It defaults to attribute `:retentiontime` and `NaN`.
* [`mmi`](@ref) `-> AbstractFloat`: monoisotopic mass.
* [`molarmass`](@ref) `-> AbstractFloat`: molar mass.
* [`mz`](@ref) `-> AbstractFloat`: mass to charge ratio (m/z).
* [`chemicalparent`](@ref) `-> AbstractChemical`: parent chemical without delocalized isotopes replacement.
* [`isotopomersisotopes`](@ref) `-> Vector{Pair{String, Int}}`: delocalized isotopes replacement of isotopomers.
* [`mass_shift_index`](@ref) `-> Int`: nominal index of mass shift between exact mass and monoisotopic mass using mass difference of an isotope and its parent.
* [`groupedisotopomersisotopes`](@ref) `-> Vector{Vector{Pair{String, Int}}}`: delocalized isotopes replacements of each isotopomers.
* [`groupedisotopomersabundance`](@ref) `-> AbstractFloat`: abundance of each isotopomers.
* [`chemicalentity`](@ref) `-> AbstractChemical`: a single chemical entity representing the chemical.
* [`chemicalspecies`](@ref) `-> Vector{<: AbstractChemical}`: multiple chemical entities having shared properties. 
* [`chemicaltransition`](@ref) `-> Vector{<: AbstractChemical}`: chemical entities analyzed in each stage of instrumental analysis.  
* [`inputchemical`](@ref) `-> AbstractChemical`: a single chemical entity that is the input at the very beginning of instrumental analysis. 
* [`outputchemical`](@ref) `-> AbstractChemical`: a single chemical entity that is the output at the very end of instrumental analysis. 
* [`analyzedchemical`](@ref) `-> AbstractChemical`: a single chemical entity directly detected at the very beginning of instrumental analysis. 
* [`detectedchemical`](@ref) `-> AbstractChemical`: a single chemical entity directly detected at the very end of instrumental analysis. 
* [`seriesanalyzedchemical`](@ref) `-> Vector{<: AbstractChemical}`: chemical entities directly analyzed in each stage of instrumental analysis.
* [`detectedcharge`](@ref) `-> Int`: charge of detected chemical.
* [`detectedisotopes`](@ref) `-> Vector{Pair{String, Int}}`: isotopes replacement of detected chemical.
* [`seriesanalyzedcharge`](@ref) `-> Vector{Int}`: charge of serially analyzed chemicals.
* [`seriesanalyzedisotopes`](@ref) `-> Vector{Vector{Pair{String, Int}}}`: isotope replacements of serially analyzed chemicals.
* [`msstage`](@ref) `-> Int`: number of stages of MS the chemical has been through.

Specific methods for these attributes are defined for other intrinsic chemical types at different chemical levels.
* Entity Level: attribute of the corresponding chemical entity.
* Species Level: attribute of the corresponding chemical species.
* Transition Level: attribute of the corresponding chemical transition.

Base functions, `copy`, `hash`, and `==`, should be overloaded for any new concrete types.
"""
abstract type AbstractChemical <: AbstractChemicalScheme end

"""
    AbstractScheme <: AbstractChemicalScheme

Abstract type for all kinds of chemical schemes.

The following attributes are implemented.
* [`elementalscheme`](@ref) `-> AbstractScheme`.
* [`structuralscheme`](@ref) `-> AbstractScheme`.
* [`chemicalname`](@ref) `-> String`.
* [`chemicalformula`](@ref) `-> String`.
* [`chemicalelements`](@ref) `-> Vector{Pair{String, Int}}`.
* [`chemicalabbr`](@ref) `-> String`: abbreviation.
* [`charge`](@ref) `-> Int`: charge state; positive for cation, negative for anion.
* [`ncharge`](@ref) `-> Int`: number of charges.
* [`mmi`](@ref) `-> AbstractFloat`: monoisotopic mass.
* [`molarmass`](@ref) `-> AbstractFloat`: molar mass.
* [`chemicalparent`](@ref) `-> AbstractChemical`: parent scheme without delocalized isotope replacement.
* [`isotopomersisotopes`](@ref) `-> Vector{Pair{String, Int}}`: delocalized isotopes replacement of isotopomers.
* [`mass_shift_index`](@ref) `-> Int`: nominal index of mass shift between exact mass and monoisotopic mass using mass difference of an isotope and its parent.
* [`groupedisotopomersisotopes`](@ref) `-> Vector{Vector{Pair{String, Int}}}`: delocalized isotopes replacements of each isotopomers.
* [`groupedisotopomersabundance`](@ref) `-> AbstractFloat`: abundance of each isotopomers.

Specific methods for these attributes are defined for other intrinsic scheme types at different scheme levels.
* Entity Level: attribute of the elemental scheme.
* Species Level: attribute of the scheme itself.
* Transition Level: only apply to schemes in `ChemicalTransition`; `Species Level` for each scheme.

Base functions, `copy`, `hash`, and `==`, should be overloaded for any new concrete types.
"""
abstract type AbstractScheme <: AbstractChemicalScheme end
# mt, ccs 
include(joinpath("type", "chemical.jl"))
include(joinpath("type", "scheme.jl"))
include(joinpath("type", "parser.jl"))
include(joinpath("type", "utils.jl"))
include(joinpath("type", "msanalysis.jl"))
include("interface.jl")
include("attr.jl")
include("chemical.jl")
include("scheme.jl")
include("elements.jl")
include("measure.jl")
include("isotopologues.jl")
include("msanalyzer.jl")
include("msanalysis.jl")
include("coelution.jl")
include("plot.jl")
include("input.jl")
include("output.jl")
include("utils.jl")

end