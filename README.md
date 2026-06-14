# MassSpecChemicals

[![Build Status](https://github.com/yufongpeng/MassSpecChemicals.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/yufongpeng/MassSpecChemicals.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/yufongpeng/MassSpecChemicals.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yufongpeng/MassSpecChemicals.jl)

A package for representing molecules or ions formed in mass spectrometers (MS). 

All chemicals are instances of abstract type `AbstractChemical`.
Charged chemicals with specific adduct or molecule loss that formed in MS (adduct ion) are instances of abstract type `AbstractAdductIon`.

# Built-in chemical types
1. `Chemical`: unstructured chemicals, storing name, elements, and other attributes

    ```julia
    Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}; property...)

    Chemical(name::AbstractString, formula::String; property...)

    Chemical(name::String, elements::Vector{Pair{String, Int}}, property::Vector{Pair{Symbol, Any}})
    ```

2. `FormulaChemical`: unstructured chemicals using formula as name

    ```julia
    FormulaChemical(elements::Vector{Pair{String, Int}}; property...)

    FormulaChemical(formula::String; property...)

    FormulaChemical(elements::Vector{Pair{String, Int}}, property::Vector{Pair{Symbol, Any}})
    ```

2. `AdductIon`: charged chemicals with specific adduct or molecule loss.

    ```julia
    AdductIon(core::AbstractChemical, adduct::AbstractScheme, ncore::Int)

    AdductIon(core::AbstractChemical, adduct::AbstractString)
    ```
3. `ChemicalTransition`: MS/MS transition.

    ```julia
    ChemicalTransition(transition::Vector)

    ChemicalTransition(precursor::AbstractChemical, products...)
    ```

4. `Isobars`: multiple chemicals with similar m/z.
        
    ```julia
    Isobars(chemical::Vector{<: AbstractChemical}, abundance::VecOrMat)
    ```

5. `Isotopomers`: multiple chemicals differed from isotopic replacement location
        
    ```julia
    Isotopomers(parent::AbstractChemical, isotopes::Vector{Pair{String, Int}})

    Isotopomers(parent::AbstractChemical, fullformula::String)

    Isotopomers(parent::AbstractChemical, fullelements::Dictionary)

    Isotopomers(parent::AbstractChemical, fullelements::Vector{Pair{String, Int}})
    ```
6. `Groupedisotopomers`: isotopomers grouped by isotopomer state
    ```julia
    Groupedisotopomers(parent::AbstractChemicalsSchema, state::Int, isotope::String, isotopes::Vector{Vector{Pair{String, Int}}}
    , abundance::Vector)
    ```

Users can parse strings and pairs by `parse_chemical`.

```julia
parse_chemical("[H3O]+") # H3O with 1 positive charge
parse_chemical("[H2O+H]+") # Protonated H2O 
parse_chemical("[C6H12O6+H]+ -> [H2O+H]+") # Protonated C6H12O6 framented to Protonated H2O 
parse_chemical("[C6H12O6+H]+" => "-H2O") # Protonated C6H12O6 and neutral loss H2O
parse_chemical("[C6H12O6+2H]2+" => "-[H2O+H]+") # DiProtonated C6H12O6 and loss Protonated H2O 
parse_chemical(ChemicalTransitionParser(ChemicalExpressionParser(; charge = 2, loss = 1)), "C6H14O6" => "-H3O") # C6H14O6 (charge = 2) and loss H3O (charge = 1)
parse_chemical("[C6H12O6+2H]2+" => ChemicalLoss(Water())) # mixing with other chemical type 
```

# Elements
## Parent Elements and Major isotopes
|Symbol|Major isotopes|Atomic number|Mass number|
|--------|-------------|-----------|------------------|
|C|[12C]|6|12|
|H|[1H]|1|1|
|O|[16H]|8|16|
|N|[14N]|7|14|
|P|[31P]|15|31|
|S|[32S]|16|32|
|Li|[7Li]|3|7|
|Na|[23Na]|11|23|
|K|[39K]|19|39|
|F|[19F]|9|19|
|Cl|[35Cl]|17|35|
|Ag|[108Ag]|47|108|
|Se|[80Se]|34|80|

## Minor Isotopes
|Minor isotopes|Atomic number|Mass number|Alternative symbol|
|--------|-------------|-----------|------------------|
|[13C]|6|13||
|D|1|2|[2H]|
|[17O]|8|17||
|[18O]|8|18||
|[15N]|7|15||
|[33S]|16|33||
|[34S]|16|34||
|[36S]|16|36||
|[6Li]|3|6||
|[40K]|19|40||
|[41K]|19|41||
|[37Cl]|17|37||
|[109Ag]|47|109||
|[74Se]|34|74||
|[76Se]|34|76||
|[77Se]|34|77||
|[78Se]|34|78||
|[82Se]|34|82||


By default, parent elements are considered as major isotopes possibly replaced by minor isotopes. For instance,
* CO2 has a carbon-12 and two oxygen-16, but any minor isotopes replacements are possible.
* [13C][16O]O has a carbon-13, an oxygen-16, and an oxygen-16 possibly replaced by other minor isotopes.

One exception is that in `parent` chemical of `Isotopomers`, parent elements are major isotopes, and the number of replacement is restricted by the field `isotopes`. 

User can use 
```julia
    set_elements!(element, mass, abundance; minor_name = nothing)
```
to add new elements with mass and natural abundance of all isotopes.
Custumized minor element names (`minor_name`) is optional.  

# AbstractSchema
Any chemical gain and loss is an instance of `AbstractSchema`. This type has three abstract sub types:
1. `AbstractElementalScheme`: scheme containing real elemental information, including isotopic replacement.
2. `AbstractStructuralScheme`: scheme containing only structural information. This is useful for defining rule-based fragmentation.
3. `AbstractCompleteScheme`: scheme containg both elemental and structural scheme. It is the final scheme storing in `AdductIon`. 

In addition to single scheme, multiple schema are wrapped in `ChemicalSchema`. 

Predefined chemicals used in adduct:
|Chemical|abbreviation|
|-------------|-----------------|
|`Electron`|`"e"`|
|`Proton`|`"H"`|
|`Water`|`"H2O"`|
|`Ammonia`|`"NH3"`|
|`Ammonium`|"`[NH4]+`"|
|`Sodium`|`"[Na]+"`|
|`Potassium`|`"[K]+"`|
|`Silver`|`"[Ag]+"`|
|`OAc`|`"[OAc]-"`|
|`OFo`|`"[OFo]-"`|
|`Fluiride`|`"[F]-"`|
|`Chloride`|`"[Cl]-"`|
|`Methinium`|`"[Me]+"`|
|`AceticAcid`|`"HOAc"`|
|`FormicaAcid`|`"HOFo"`|
|`Methylacetate`|`"MeOAc"`|
|`Methylformate`|`"MeOFo"`|

# API
## Attributes of `AbstractChemical`
Attributes are interfaces directly accessing properties and fields through `getchemicalproperty` or deriving values from other attributes.

|Attribute|Return type|Description|
|----|----|-----------|
|`chemicalname`|`String`|unique chemical name|
|`chemicalformula`|`String`|chemical formula|
|`chemicalelements`|`Vector{Pair{String, Int}}`|chemical elements| 
|`chemicalabbr`|`String`|common abbreviation; defaults to `chemicalname`| 
|`chemicalsmiles`|`String`|SMILES; defaults to `""`| 
|`charge`|`Int`|charges (positive or negative); defaults to 0|
|`ncharge`|`Int`|number of charges|
|`retentiontime`|`Float64`|retention time; defaults to `NaN`| 
|`chemicalparent`|`AbstractChemical`|parent chemical without delocalized isotopes replacement| 
|`isotopomersisotopes`|`Vector{Pair{String, Int}}`|delocalized isotopes replacement of isotopomers|
|`isotopomerstate`|`Int`|isotopomers state, i.e. equivalent number of isotope|
|`groupedisotopomersisotopes`|`Vector{Vector{Pair{String, Int}}}`|delocalized isotopes replacements of each isotopomers|
|`groupedisotopomersabundance`|`AbstractFloat`|abundance of each isotopomers|
|`chemicalentity`|`AbstractChemical`|a single chemical entity representing the chemical| 
|`chemicalspecies`|`Vector{<: AbstractChemical}`|multiple chemical entities having shared properties| 
|`chemicaltransitions`|`Vector{<: AbstractChemical}`|chemical entities analyzed in each stage of instrumental analysis|
|`inputchemical`|`AbstractChemical`|a single chemical entity that is the input of the very begining of instrumental analysis| 
|`outputchemical`|`AbstractChemical`|a single chemical entity that is the output of the very ending of instrumental analysis| 
|`analyzedchemical`|`AbstractChemical`|a single chemical entity that is directly analyzed in the very begining of instrumental analysis| 
|`detectedchemical`|`AbstractChemical`|a single chemical entity that is directly detected in the very ending of instrumental analysis| 
|`detectedisotopes`|`Vector{Pair{String, Int}}`|delocalized isotopes replacement of detected chemical| 
|`detectedcharge`|`Int`|charge state of detected chemical| 
|`detectedelements`|`Vector{Pair{String, Int}}`|elements of detected chemical| 
|`seriesanalyzedchemical`|`Vector{<: AbstractChemical}`|chemical entites that are directly analyzed in each stage of instrumental analysis| 
|`seriesanalyzedisotopes`|`Vector{Vector{Pair{String, Int}}}`|delocalized isotopes replacements of serially analyzed chemical| 
|`seriesanalyzedcharge`|`Vector{Int}`|charge state of serially analyzed chemical| 
|`seriesanalyzedelements`|`Vector{Vector{Pair{String, Int}}}`|elements of serially analyzed chemical| 
|`msstage`|`Int`|number of stages of MS the chemical has been through| 
|`mmi`|`Float64`|monoisotopic mass|
|`molarmass`|`Float64`|molar mass|
|`mz`|`Float64`|m/z; mass-to-charge ratio|

Specific Methods for the attributes are defined for other intrinsic chemical type on different chemical level
* Entity Level: attribute of the corresponding chemical entity
* Species Level: attribute of the corresponding chemical species
* Transition Level: attribute of the corresponding chemical transitions

See documentation of each attribute for the exact level.

### Type-specific attributes 
Users can define new attributes or directly overload existing attribute functions for any chemical types
```julia
abstract type Lipid <: AbstractChemical end 
struct Fattyacid <: Lipid 
    ncarbon::Int
    ndoublebond::Int 
end
struct Acylglycerol <: Lipid 
    carbonchains::Vector{Fattyacid}
end
# Different attribute methods for different chemical types 
ncarbon(chemical::Fattyacid) = chemical.ncarbon 
ncarbon(chemical::Acylglycerol) = sum(ncarbon, chemical.carbonchains) + 3
fa1 = Fattyacid(18, 0)
fa2 = Fattyacid(18, 1)
dg = Acylglycerol([fa1, fa2])
ncarbon(fa1) == 18 
ncarbon(dg) == 39 
```
### Instance-specific attributes 
For attributes defined in only some instances, users can use `getchemicalproperty` to handle the missing attributes. 
```julia
# Access the property through :carbonchains getchemicalproperty, default to [chemical]
carbonchains(chemical::Lipid) = getchemicalproperty(chemical, :carbonchains, [chemical])
carbonchains(dg) == [fa1, fa2]
carbonchains(fa1) == [fa1]
```
The type `Chemical` stores any non-default attributes in the field `property`, user can create object with these attributes as keyword arguments or add attributes by directly pushing the attr_name-attr_value pair to `chemical.property`, and define the attrubte function to access it. 
```julia
    chemical = Chemical("18:0 PC-d9", "	C44H79NO8PD9"; lipidclass = "PC")
    
    lipidclass(chemical::Chemical) = getchemicalproperty(chemical, :lipidclass) # new attribute
    lipidclass(chemical) == "PC"

    push!(chemical.property, :retentiontime => 10)
    retentiontime(chemical) == 10 # Defined as accessing :retentiontime through getchemicalproperty
```

### Additional Atttributes of `AbstractAdductIon`
|Attribute|Return type|Description|
|----|----|-----------|
|`ioncore`|`AbstractChemical`|core chemical|
|`ionadduct`|`AbstractAdduct`|adduct originated from ionization|
|`ncore`|`Int`|number of core chemical|

When isotopes are involved in addut ion formation for an object `adduct_ion` which `chemical = ioncore(adduct_ion)::ChemicalType` and `adduct = ionadduct(adduct_ion)::Existing_Scheme`, there are two solutions.
1. If `ChemicalType` is a customized chemical type, define type-specific `completescheme`
    ```julia
    completescheme(adduct_ion::ChemicalType, ::Affected_Scheme)
    completescheme(adduct_ion::AbstractAdductIon{ChemicalType, Existing_Scheme}, ::Affected_Scheme)
    ```
    
    `completescheme` returns `CompleteSchema` which wraps `Affected_Scheme` and corresponding elemental scheme. For instance, [M-Me]- of Deuterium-labeled phosphatidylcholine (`DL_PC`) may turn out to be [M-CD3]- (`ElementalScheme{false, DL_Methenium}`) rather than [M-CH3]- (`ElementalScheme{false, Methenium}`) if Deuteriums are labeled on the methyl group of choline. In this case, extend `completescheme(::DL_PC, ::ElementalScheme{false, Methenium})` such that
    ```julia 
    # dlmcpc: [M-Me]- of Deuterium-labeled phosphatidylcholine on the methyl group of choline
    # dlpc: [M-Me]- of Deuterium-labeled phosphatidylcholine on other part
    completescheme(dlmcpc, x::ElementalScheme{false, Methenium}) == StructuralElementalScheme(x, ElementalScheme(false, DL_Methenium()))
    completescheme(dlpc, x::ElementalScheme{false, Methenium}) == StructuralElementalScheme(x, x)
    ```
2. If `ChemicalType` is `Chemical`, define an attribute `:structure` for the `chemical`. The attribute should be ionadduct-(scheme-scheme pairs) pairs. `:structure` finds this attribute, and extracts the value of key `adduct`. 
    ```julia
    chemical = Chemical("18:0 PC-d9", "	C44H79NO8PD9")
    push!(chemical.property, :structure => [nothing => [ElementalScheme(false, Methenium()) => StructuralElementalScheme(x, ElementalScheme(false, DL_Methenium()))]])
    completescheme(chemical, ChemicalGain(Proton())) == StructuralElementalScheme(ChemicalGain(Proton()), ChemicalGain(Proton())) # No key Protonation()
    completescheme(chemical, ChemicalLoss(Methenium())) == StructuralElementalScheme(ChemicalLoss(Methenium()), ChemicalLoss(Methenium()))
    ```
## Isotopic abundance and Isotopologues
There are three related functions
* `isotopicabundance`
    
    This function computes isotopic abundance of the input elements composition (`Vector{Pair{Int}}` or `Dictionary`), formula, and chemical (converted to elements by `chemicalelements`). Parent elements are viewed as major isotopes, and isotopic abundances of all elements are considered in computation. To compute isotopic abundance of chemicals with all isotopes labeled intentionally and not following natural distribution, set keyword argument `ignore_isotopes` true, and only parent elements are considered.

* `Isotopologues`

    This function computes isotopologues of formula, single chemical, MS/MS precursor-product pairs (formula pairs or `ChemimcalTransition`) or multiple chemicals. Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. Isotopologues can be filtered by abundance threshold. MS/MS product can be neutral loss/gain or chemical loss/gain. 

* `TandemIsotopologues`

    This function is similar to `Isotopologues`; it computes isotopologues of given precursor(s) and additionally computes the abundance of fragments with given fragmentation patterns. The key difference is that this function is recursive and abundance is calculated from the begining. It generally performs slightly slower for multiple MS stages and abundance is normalized in the first stage and filtered in all stages.

Isotopologues table can be aggregated using `group_isotopologues`.
## Mass Spectrometer 
There are five functions to simulate ions in mass spectrometer. 
1. `Isolation`: isolating target ion(s) with specific m/z values and resolution to enter the next MS stage. 
2. `AllIons`: allow all Ions within m/z range entering the next MS stage. 
3. `Fragmentation`: create a table of fragments with given fragmentation patterns. 
4. `MSScan`: perform MS scan using given mass analyzer. This function creates `Spectrum` objects, which can be visualized by `plot_spectrum` and `plot_spectrum!`. Peak lists can be extracted with function `peak_table`. 
5. `SelectedIonMonitor`: perform selected ion monitoring for target transition(s). Peak lists can be extracted with function `peak_table`. 

The following common mass analyzer are defined.
1. `Quadrupole`
2. `QuadrupoleIon` or `QIT`
3. `LinearIonTrap` or `LIT`
4. `TimeOfFlight` or `TOF`
5. `Orbitrap`
6. `FourierTransformIonCyclotronResonance` or `FTICR`

Default settings related to resolution, and isolation window are also defined for each analyzer. To create generic mass analyzer, use `MSAnalyzer`.
## Co-eluting isobars
The function `CoelutingIsobars` creates an object `CoelutingIsobars` with a vector of elution function-criteria pair, a vector of ms analyzer-criteria pair, and a target chemical table.

This object can be further aggregated using `isobar_table`.

## Other Functions
* `parent_element`: return parent element.
* `major_isotope`: return major isotope.
* `minor_isotope`: return minor isotope.
* `iselement`: whether it is an element.
* `isisotope`: whether it is an isotope (including element).
* `ischemicalequal`: whether two chemicals are chemically equivalent.
* `ischemicalequaltransform`: return an object for comparison with other chemicals by `istransformedadduct`.
* `istransformedchemicalequal`: whether two chemicals are chemically equivalent after applying `istransformedchemicalequal`.
* `isadductequal`: whether two adducts chemically equivalent.
* `isadductequaltransform`: return an object for comparison with other adducts by `istransformedadduct`.
* `istransformedadduct`: whether two adducts are chemically equivalent after applying `isadductequaltransform`.
* `match_chemical`: match detected chemicals with reference library.
* `acrit`: create absolute criterion.
* `rcrit`: create relative criterion.
* `crit`: create both absolute and relative criterion.
* `@ri_str`: real number interval.
* `plot_resolving_power` and `plot_resolving_power!`: plot the function of m/z to resolving_power.
* `plot_window` and `plot_window!`: plot the window function.