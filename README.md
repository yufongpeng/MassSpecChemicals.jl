# MassSpecChemicals

[![Build Status](https://github.com/yufongpeng/MassSpecChemicals.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/yufongpeng/MassSpecChemicals.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/yufongpeng/MassSpecChemicals.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yufongpeng/MassSpecChemicals.jl)

A package for representing molecules or ions formed in mass spectrometers (MS). 

All chemicals are instances of abstract type `AbstractChemical`.
Charged chemicals with specific adduct or molecule loss that formed in MS (adduct ion) are instances of abstract type `AbstractAdductIon`.

# Built-in chemical types
1. `Chemical`: unstructured chemicals, storing name, elements, and other attributes

    ```julia
        Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}; attr...)

        Chemical(name::AbstractString, formula::String; attr...)

        Chemical(name::String, elements::Vector{Pair{String, Int}}, attr::Vector{Pair{Symbol, Any}})
    ```

2. `AdductIon`: charged chemicals with specific adduct or molecule loss.

    ```julia
        AdductIon(core::AbstractChemical, adduct::AbstractAdduct)
    ```

3. `Isobars`: multiple chemicals with similar m/z.
        
    ```julia
        Isobars(chemical::Vector{<: AbstractChemical}, abundance::Vector{Float64})
    ```

4. `Isotopomers`: multiple chemicals differed from isotopic replacement location
        
    ```julia
        Isotopomers(parent::AbstractChemical, isotopes::Vector{Pair{String, Int}})
    ```
        
5. `ChemicalLoss`: chemical loss from a precursor.

    ```julia
        ChemicalLoss(chemical::AbstractChemical)
    ```

6. `ChemicalPair`: a pair of precursor and product in MS/MS.

    ```julia
        ChemicalPair(precursor::AbstractChemical, product::AbstractChemical)
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

# Adduct
All adducts are instances of abstract type `AbstractAdduct` and `AbstractPosAdduct` or `AbstractNegAdduct`.

Prefefined non-generic adducts:
|Adduct Object|Adduct expression|
|-------------|-----------------|
|`LossElectron`|[M]+|
|`Protonation`|[M+H]+|
|`ProtonationNLH2O`|[M+H-H2O]+|
|`ProtonationNL2H2O`|[M+H-2H2O]+|
|`ProtonationNL3H2O`|[M+H-3H2O]+|
|`DiProtonation`|[M+2H]2+|
|`TriProtonation`|[M+3H]3+|
|`AddNH4`|[M+NH4]+|
|`AddHNH4`|[M+H+NH4]2+|
|`Add2NH4`|[M+2NH4]2+|
|`Sodization`|[M+Na]+|
|`SodizationProtonation`|[M+Na+H]2+|
|`DiSodization`|[M+2Na]2+|
|`AddElectron`|[M]-|
|`Deprotonation`|[M-H]-|
|`DeprotonationNLH2O`|[M-H-H2O]-|
|`DiDeprotonation`|[M-2H]2-|
|`TriDeprotonation`|[M-3H]3-|
|`AddOAc`|[M+OAc]-|
|`AddOFo`|[M+OFo]-|
|`Fluoridation`|[M+F]-|
|`Chloridation`|[M+Cl]-|
|`Demethylation`|[M-Me]-|

User can customize adduct by generic type `PosAdduct` and `NegAdduct` or define subtypes of `AbstractPosAdduct` or `AbstractNegAdduct` for more customized usage.

For example, [2M+H]+, and [M-2H-2H2O]2- can be created by 
```julia
    # kmer, adductformula, ncharge
    PosAdduct(2, "+H", 1)
    NegAdduct(1, "-2H-2H2O", 2) 
```

# API
## Attributes of `AbstractChemical`
`getchemicalattr(chemical::ChemicalType, attr::Symbol; kwargs...)` returns the attribute `attr` of `chemical` with additional settings by `kwargs`. 

To define new attributes or overwrite defaults for `ChemicalType`, define the following method 
```julia
    getchemicalattr(chemical::ChemicalType, attr::Val{attr}; kwargs...)
``` 

For attributes exclusive for the instance, the attibutes have to be stored in the instance and `getchemicalattr(chemical::ChemicalType, attr::Val{attr}; kwargs...)` is defined as extracting the value from `chemical`. For example, the default `Chemical` stores any non-default attributes in the field `attr`, and `getchemicalattr` extracts value from the field. User can also add attributes by directly pushing the attr_name-attr_value pair to `chemical.attr`.
```julia
    chemical = Chemical("18:0 PC-d9", "	C44H79NO8PD9")
    push!(chemical.attr, :rt => 10)
    rt(chemical) == 10
```

The following getter functions are user friendly interfaces for `getchemicalattr`.
|Function|Attr|Return type|Description|
|--------|----|----|-----------|
|`chemicalname`|`:name`|`String`|unique chemical name|
|`chemicalformula`|`:formula`|`String`|chemical formula|
|`chemicalelements`|`:elements`|`Vector{Pair{String, Int}}`|chemical elements| 
|`chemicalabbr`|`:abreviation`|`String`|common abbreviation; defaults to `chemicalname`| 
|`chemicalsmiles`|`:SMILES`|`String`|SMILES; defaults to `""`| 
|`charge`|`:charge`|`Int`|charges (positive or negative); defaults to 0|
|`abundantchemical`|`:abundant_chemical`|`AbstractChemical`|the most abundant chemical from a chemical (itself) or isobars|
|`rt`|`:rt`|`Float64`|retention time; defaults to `NaN`|  

### Additional Atttributes of `AbstractAdductIon`
|Function|Attr|Return type|Description|
|--------|----|----|-----------|
|`ioncore`|`:core`|`AbstractChemical`|core chemical|
|`ionadduct`|`:adduct`|`AbstractAdduct`|adduct originated from ionization|
|`kmer`|`:kmer`|`Int`|number of core chemical|

When isotopes are involved in addut ion formation for an object `adduct_ion` which `chemical = ioncore(adduct_ion)::ChemicalTypeT` and `adduct = ionadduct(adduct_ion)::Affected_Adduct`, there are two solutions.
1. If `ChemicalType` is a customized chemical type, define type-specific `adductisotopes`
    ```julia
        adductisotopes(adduct_ion::AbstractAdductIon{Affected_Adduct, T})
    ```
    
    `adductisotopes` returns element-number pairs which are the elements changed when the core chemical has isotopic labeling that is lost during adduct formation. For instance, [M-Me]- (`Demethyl`) of Deuterium-labeled phosphatidylcholine (`DL_PC`) may turn out to be [M-CD3]- rather than [M-CH3]- if Deuteriums are labeled on the methyl group of choline. In this case, extend `adductisotopes(::AbstractAdductIon{Demethyl, DL_PC})` such that
    ```julia 
        # dlmcpc: [M-Me]- of Deuterium-labeled phosphatidylcholine on the methyl group of choline
        # dlpc: [M-Me]- of Deuterium-labeled phosphatidylcholine on other part
        adductisotopes(dlmcpc) == ["H" => 3, "D" => -3]
        adductisotopes(dlpc) == []
    ```
2. If `ChemicalType` is `Chemical`, define an attribute `:adductisotopes` for the `chemical`. The attribute should be ionadduct-(element-number pairs) pairs. `adductisotopes` finds this attribute, and extracts the value of key `adduct`. If the attribute or the key does not exist, empty vector is returned. 
    ```julia
        chemical = Chemical("18:0 PC-d9", "	C44H79NO8PD9")
        push!(chemical.attr, :adductisotopes => [Demethyl() => ["H" => 3, "D" => -3]])
        adductisotopes(Adduction(chemical, "[M+H]+")) == []
        adductisotopes(Adduction(chemical, "[M-Me]-")) == ["H" => 3, "D" => -3]
    ```

### Additional Atttributes of `Isobars`
|Function|Attr|Return type|Description|
|--------|----|----|-----------|
||`:chemicals`|`Vector{<: AbstractChemical}`|composition of isabars|
||`:abundance`|`Vector{Float64}`|abundance of each isabars|

### Additional Atttributes of `Isotopomers`
|Function|Attr|Return type|Description|
|--------|----|----|-----------|
||`:parent`|`AbstractChemical`|shared chemical structure of isotopomers prior to isotopic replacement|
||`:isotopes`|`Vector{Pair{String, Int}}`|isotopes-number pairs of isotopic replacement|

### Additional Atttributes of `ChemicalPair`
|Function|Attr|Return type|Description|
|--------|----|----|-----------|
||`:precursor`|`AbstractChemical`|chemical prior to fragmentation|
||`:product`|`AbstractChemical`|one of the fragments|

### Additional Atttributes of `ChemicalLoss`
|Function|Attr|Return type|Description|
|--------|----|----|-----------|
||`:chemical`|`AbstractChemical`|chemical loss from a precursor|


### Functions for Attributes derived from other attributes
|Function|Required attr|Return type|Description|
|--------|----|----|-----------|
|`mmi`|`:elements` or `:formula`|`Float64`|monoisotopic mass|
|`molarmass`|`:elements` or `:formula`|`Float64`|molar mass|
|`mz`|`:elements` or `:formula` and `:charge`|`Float64`|m/z|
|`ncharge`|`:charge`|`Int`|number of charges|
|`isotopicabundance`|`:elements` or `:formula`|`Float64`|isotopic abundance of a single monoisotopic chemical or formula|
|`isotopologues`|`:elements` or `:formula`|`Vector{<:　ＡｂｓｔｒａｃｔＣｈｅｍｉｃａｌ}`|isotopologues of a single chemical or formula, and MS/MS transition|
|`isotopologues_table`|`:elements` or `:formula`|`Table`|a table containing isotopologues information|
|`isobars_rt`|`:elements` or `:formula`, and `rt`|`Isobars`|isobars under specific m/z, mass, and fwhm|
|`isobars_rt_table`|`:elements` or `:formula`, and `rt`|`Table`|a table containing co-eluting isobars information|

## Functions for Attributes of `AbstractAdduct`
|Function|Return type|Description|
|--------|----|-----------|
|`kmer`|`Int`|number of core chemical|
|`charge`|`Int`|charges (positive or negative); defaults to 0|
|`ncharge`|`Int`|number of charges|
|`adductformula`|`String`|the formula for adduct. For instance, `"-H"` for [M-H]-, `"+OAc"` for [M+OAc]-.|
|`adductelements`|`Vector{Pair{String, Int}}`|elements changed with adduct formation. For generic adduct, it uses `adductformula` to calculate elements. If non-element strings are used in `adductformula`, defining customized `adductelements` is required.| 

## Isotopic abundance and Isotopologues
There are three related functions
* `isotopicabundance`
    
    This function computes isotopic abundance of the input elements composition (`Vector{Pair{Int}}` or `Dictionary`), formula, and chemical (converted to elements by `chemicalelements`). Parent elements are viewed as major isotopes, and isotopic abundances of all elements are considered in computation. To compute isotopic abundance of chemicals with all isotopes labeled intentionally and not following natural distribution, set keyword argument `ignore_isotopes` true, and only parent elements are considered.

* `isotopologues`

    This function computes isotopologues of formula, single chemical, or MS/MS precursor-product pairs (formula pairs or `ChemicalPair`). Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. Isotopologues can be filtered by abundance threshold, and grouped by m/z or monoisotopic mass. MS/MS product can be neutral loss or chemical loss. 

* `isotopologues_table`

    It is similar to `isotopologues` but stores all information in a `Table`, including mass or m/z, and abundance. Multiple chemicals are supported.

## Co-eluting isobars
There are two related functions
* `isobars_rt`

    This function finds all potential isobars of a detected chemical from a library (`AbstractVector` or `Table`). It returns only isobars without any additional information.

* `isotopologues_table`

    This function works very similar, but stores information including difference of retention time, m/z, ID, etc, in a `Table`. Another difference is not restricted to one chemical (`AbstractChemical`, `AbstractVector` or `Table`). 

## Other Utility Functions
* `ischemicalequal`: whether two chemicals are chemically equivalent.
* `ischemicalequaltransform`: return an object for comparison with other chemicals by `istransformedadduct`.
* `istransformedchemicalequal`: whether two chemicals are chemically equivalent after applying `istransformedchemicalequal`.
* `isadductequal`: whether two adducts chemically equivalent.
* `isadductequaltransform`: return an object for comparison with other adducts by `istransformedadduct`.
* `istransformedadduct`: whether two adducts are chemically equivalent after applying `isadductequaltransform`.
* `acrit`: create absolute criterion.
* `rcrit`: create relative criterion.
* `crit`: create both absolute and relative criterion.
* `@ri_str`: real number interval.
