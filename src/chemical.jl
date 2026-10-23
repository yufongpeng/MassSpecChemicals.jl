"""
    match_chemical(exp, lib; colexp = :Chemical, collib = :Chemical, fnexp = detectedchemical) -> Table

Match chemicals in `exp` (a `Table` or vector) against those in `lib` (a `Table` or vector).
The returned table contains the matching rows from `exp`, the matching library indices in the `LibID` column, the matched chemicals in the `Match` column, and any other columns from `lib`.

Before matching, each chemical in `exp` is transformed by `fnexp`.
"""
function match_chemical(exp, lib; colexp = :Chemical, collib = :Chemical, fnexp = detectedchemical)
    del = Int[]
    libid = Int[]
    exp = hasproperty(exp, colexp) ? exp : Table(; Chemical = exp)
    chemical_exp = fnexp.(getproperty(exp, colexp))
    chemical_lib = hasproperty(lib, collib) ? getproperty(lib, collib) : lib
    for i in eachindex(exp)
        j = findfirst(x -> ischemicalequal(x, chemical_exp[i]), chemical_lib)
        isnothing(j) ? push!(del, i) : push!(libid, j)
    end
    id = setdiff(eachindex(exp), del)
    ps = filter(!=(collib), propertynames(lib))
    Table(exp[id]; LibID = libid, Match = chemical_lib[libid], [p => [getproperty(lib, p)[i] for i in libid] for p in ps]...)
end

"""
    ischemicalequal(x::AbstractChemicalScheme, y::AbstractChemicalScheme) -> Bool

Determine whether two chemicals are chemically equivalent.
By default, both chemicals are transformed with [`ischemicalequaltransform`](@ref) and compared with [`istransformedchemicalequal`](@ref).

For chemicals, their names, elemental composition, and other attributes for wrapper types  are compared.
For scehemes, internal chemicals and other attributes for wrapper types are compared.
"""
ischemicalequal(x::AbstractChemical, y::AbstractChemical) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))
ischemicalequal(x::AbstractScheme, y::AbstractScheme) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))
ischemicalequal(x::Isobars, y::Isobars) = istransformedchemicalequal(x, y)
ischemicalequal(x::Isotopomers, y::Isotopomers) = istransformedchemicalequal(x, y)
ischemicalequal(x::Groupedisotopomers, y::Groupedisotopomers) = istransformedchemicalequal(x, y)
ischemicalequal(x::ChemicalTransition, y::ChemicalTransition) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))
ischemicalequal(x::ChemicalSchemes, y::ChemicalSchemes) = istransformedchemicalequal(x, y)
ischemicalequal(x::StructuralElementalScheme, y::StructuralElementalScheme) = istransformedchemicalequal(x, y)
ischemicalequal(x::ElementalScheme{true}, y::ElementalScheme{true}) = istransformedchemicalequal(x, y)
ischemicalequal(x::ElementalScheme{false}, y::ElementalScheme{false}) = istransformedchemicalequal(x, y)
ischemicalequal(x::IsotopomerizedSchemes, y::IsotopomerizedSchemes) = istransformedchemicalequal(x, y)
ischemicalequal(x::Groupedisotopomerizedschemes, y::Groupedisotopomerizedschemes) = istransformedchemicalequal(x, y)
ischemicalequal(x::AbstractChemicalScheme, y::AbstractChemicalScheme) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))

"""
    ischemicalequaltransform(x::AbstractChemicalScheme) -> AbstractChemicalScheme

Return the transformed object used when comparing chemicals with [`istransformedchemicalequal`](@ref).
"""
ischemicalequaltransform(x::AbstractChemical) = x 
ischemicalequaltransform(x::AbstractScheme) = x 
ischemicalequaltransform(x::Isobars) = length(x) == 1 ? ischemicalequaltransform(chemicalentity(x)) : x
ischemicalequaltransform(x::Isotopomers) = isempty(unique_elements(x.isotopes)) ? x.parent : x 
ischemicalequaltransform(x::Groupedisotopomers) = length(x.isotopes) > 1 ? x : isempty(unique_elements(x.isotopes[begin])) ? x.parent : Isotopomers(x.parent, x.isotopes[begin]) 
# ischemicalequaltransform(x::IsotopomerizedSchemes) = isempty(unique_elements(x.isotopes)) ? x.parent : x 
# ischemicalequaltransform(x::Groupedisotopomerizedschemes) = length(x.isotopes) > 1 ? x : isempty(unique_elements(x.isotopes[begin])) ? x.parent : IsotopomerizedSchemes(x.parent, x.isotopes[begin]) 
ischemicalequaltransform(x::ElementalScheme{T}) where T = ElementalScheme(T, ischemicalequaltransform(x.chemical))
ischemicalequaltransform(x::StructuralElementalScheme) = StructuralElementalScheme(ischemicalequaltransform(x.structuralscheme), ischemicalequaltransform(x.elementalscheme))
ischemicalequaltransform(x::ChemicalTransition) = ChemicalTransition([ischemicalequaltransform(c) for c in chemicaltransition(x)])
ischemicalequaltransform(x::ChemicalSchemes) = (length(x.schemes) == 1 && first(x.number) == 1) ? first(x.schemes) : x
ischemicalequaltransform(x::IsotopomerizedSchemes) = isempty(unique_elements(x.isotopes)) ? ischemicalequaltransform(x.parent) : (length(x.parent.schemes) == 1 && first(x.parent.number) == 1) ? isotopomerize(first(x.parent.schemes), x.isotopes) : x
ischemicalequaltransform(x::Groupedisotopomerizedschemes) = length(x.isotopes) > 1 ? x : isempty(unique_elements(x.isotopes[begin])) ? ischemicalequaltransform(x.parent) : isotopomerize(ischemicalequaltransform(x.parent), x.isotopes[begin])
"""
    istransformedchemicalequal(x::AbstractChemicalScheme, y::AbstractChemicalScheme) -> Bool

Determine whether two chemicals are chemically equivalent after applying [`ischemicalequaltransform`](@ref).

For chemicals, their names, elemental compositions, and other attributes for wrapper types are compared.
For scehemes, internal chemicals and other attributes for wrapper types are compared.
"""
istransformedchemicalequal(x::AbstractChemicalScheme, y::AbstractChemicalScheme) = false
istransformedchemicalequal(x::AbstractChemical, y::AbstractChemical) = chemicalname(x) == chemicalname(y) && _chemicalelements_equal(chemicalelements(x), chemicalelements(y))
istransformedchemicalequal(x::AbstractScheme, y::AbstractScheme) = isequal(x, y)
istransformedchemicalequal(x::AbstractAdductIon, y::AbstractAdductIon) = ischemicalequal(ionadduct(x), ionadduct(y)) && ischemicalequal(ioncore(x), ioncore(y))
istransformedchemicalequal(x::Chemical, y::Chemical) = x.name == y.name && _chemicalelements_equal(x.elements, y.elements)
istransformedchemicalequal(x::FormulaChemical, y::FormulaChemical) = _chemicalelements_equal(x.elements, y.elements)
istransformedchemicalequal(x::AbstractChemicalWrapper, y::AbstractChemical) = false
istransformedchemicalequal(x::AbstractChemical, y::AbstractChemicalWrapper) = false
istransformedchemicalequal(x::AbstractChemicalWrapper, y::AbstractChemicalWrapper) = x == y
istransformedchemicalequal(x::Isobars, y::Isobars) = all(ischemicalequal(a, b) for (a, b) in zip(x.chemicals, y.chemicals)) && all(isapprox(a, b) for (a, b) in zip(x.abundance, y.abundance))
istransformedchemicalequal(x::Isotopomers, y::Isotopomers) = ischemicalequal(x.parent, y.parent) && x.isotopes == y.isotopes
istransformedchemicalequal(x::Groupedisotopomers, y::Groupedisotopomers) = ischemicalequal(x.parent, y.parent) && x.index == y.index && x.isotope == y.isotope && all(splat(==), zip(x.isotopes, y.isotopes)) && all(splat(isapprox), zip(x.abundance, y.abundance))
istransformedchemicalequal(x::ChemicalTransition, y::ChemicalTransition) = all(ischemicalequal.(x.transition, y.transition))
istransformedchemicalequal(x::AbstractSchemeWrapper, y::AbstractScheme) = false
istransformedchemicalequal(x::AbstractScheme, y::AbstractSchemeWrapper) = false
istransformedchemicalequal(x::AbstractSchemeWrapper, y::AbstractSchemeWrapper) = x == y
istransformedchemicalequal(x::IsotopomerizedSchemes, y::IsotopomerizedSchemes) = istransformedchemicalequal(x.parent, y.parent) && x.isotopes == y.isotopes
istransformedchemicalequal(x::Groupedisotopomerizedschemes, y::Groupedisotopomerizedschemes) = ischemicalequal(x.parent, y.parent) && x.index == y.index && x.isotope == y.isotope && all(splat(==), zip(x.isotopes, y.isotopes)) && all(splat(isapprox), zip(x.abundance, y.abundance))
function istransformedchemicalequal(x::ChemicalSchemes, y::ChemicalSchemes) 
    uk = [false for _ in eachindex(y.schemes)]
    for (kx, vx) in zip(x.schemes, x.number)
        pass = false
        for (i, ky) in enumerate(y.schemes)
            uk[i] && continue 
            if ischemicalequal(kx, ky) && vx == y.number[i]
                pass = true
                uk[i] = true
                break
            end
        end
        pass || return false
    end
    all(uk)
end
istransformedchemicalequal(x::ElementalScheme, y::ElementalScheme) = ischemicalequal(x.chemical, y.chemical)
istransformedchemicalequal(x::StructuralElementalScheme, y::StructuralElementalScheme) = ischemicalequal(structuralscheme(x), structuralscheme(y)) && ischemicalequal(elementalscheme(x), elementalscheme(y))

"""
    ionize([constructor,] chemical; kwargs...) -> constructor
    ionize([constructor,] chemical, args...; kwargs...) -> constructor

Ionize `chemical` and wrap it with `constructor`. The default constructor is [`AdductIon`](@ref).
"""
ionize(chemical::AbstractChemical; kwargs...) = ionize(AdductIon, chemical; kwargs...)
ionize(chemical::AbstractChemical, adduct, ncore = 1; kwargs...) = ionize(AdductIon, chemical, adduct, ncore; kwargs...)

"""
    ionize(::Type{AdductIon}, chemical::AbstractChemical; adduct, ncore = 1, kwargs...) -> AdductIon
    ionize(::Type{AdductIon}, chemical::AbstractChemical, adduct, ncore = 1; kwargs...) -> AdductIon

Ionize `chemical` and wrap with [`AdductIon`](@ref).

# Arguments
* `adduct::AbstractScheme`: adduct scheme.
* `ncore::Int`: number of core chemicals.
"""
ionize(::Type{AdductIon}, chemical::AbstractChemical; adduct, ncore = 1, kwargs...) = AdductIon(chemical, adduct, ncore)
ionize(::Type{AdductIon}, chemical::AbstractAdductIon; adduct, ncore = 1, kwargs...) = AdductIon(ioncore(chemical), adduct, ncore)
ionize(::Type{AdductIon}, chemical::Isobars; adduct, ncore = 1, kwargs...) = Isobars([ionize(AdductIon, x, adduct, ncore; kwargs...) for x in chemicalspecies(chemical)], chemical.abundance)
ionize(::Type{AdductIon}, chemical::Isotopomers; adduct, ncore = 1, kwargs...) = Isotopomers(ionize(AdductIon, chemicalparent(chemical), adduct, ncore; kwargs...), chemical.isotopes)
ionize(::Type{AdductIon}, chemical::Groupedisotopomers; adduct, ncore = 1, kwargs...) = Groupedisotopomers(ionize(AdductIon, chemicalparent(chemical), adduct, ncore; kwargs...), chemical.index, chemical.isotope, chemical.isotopes, chemical.abundance)

ionize(::Type{AdductIon}, chemical::AbstractChemical, adduct, ncore = 1; kwargs...) = AdductIon(chemical, adduct, ncore)
ionize(::Type{AdductIon}, chemical::AbstractAdductIon, adduct, ncore = 1; kwargs...) = AdductIon(ioncore(chemical), adduct, ncore)
ionize(::Type{AdductIon}, chemical::Isobars, adduct, ncore = 1; kwargs...) = Isobars([ionize(AdductIon, x, adduct, ncore; kwargs...) for x in chemicalspecies(chemical)], chemical.abundance)
ionize(::Type{AdductIon}, chemical::Isotopomers, adduct, ncore = 1; kwargs...) = Isotopomers(ionize(AdductIon, chemicalparent(chemical), adduct, ncore; kwargs...), chemical.isotopes)
ionize(::Type{AdductIon}, chemical::Groupedisotopomers, adduct, ncore = 1; kwargs...) = Groupedisotopomers(ionize(AdductIon, chemicalparent(chemical), adduct, ncore; kwargs...), chemical.index, chemical.isotope, chemical.isotopes, chemical.abundance)

"""
    isotopomerize(chemical::AbstractChemicalScheme, isotopes) -> AbstractChemicalScheme

Add the delocalized isotopic replacements in `isotopes` to `chemical`.
"""
isotopomerize(chemical::AbstractChemical, isotopes) = Isotopomers(chemical, isotopes)
isotopomerize(chemical::Isotopomers, isotopes) = Isotopomers(chemical.parent, gain_elements(chemical.isotopes, isotopes))
isotopomerize(sch::StructuralElementalScheme, isotopes) = StructuralElementalScheme(structuralscheme(sch), isotopomerize(elementalscheme(sch), isotopes))
# isotopomerize(sch::AbstractCompleteScheme{T,<:AbstractChemical}, isotopes) where T = isotopomerize(elementalscheme(sch), isotopes)
# isotopomerize(sch::StructuralElementalScheme{T,<:AbstractChemical}, isotopes) where T = isotopomerize(elementalscheme(sch), isotopes)
isotopomerize(sch::ElementalScheme{true}, isotopes) = ElementalScheme(true, isotopomerize(sch.chemical, isotopes))
isotopomerize(sch::ElementalScheme{false}, isotopes) = ElementalScheme(false, isotopomerize(sch.chemical, reverse_elements(isotopes, true)))
# isotopomerize(sch::ElementalScheme{false}, isotopes::Tuple) = ElementalScheme(false, isotopomerize(sch.chemical, loss_elements(isotopes...)))
isotopomerize(sch::ChemicalSchemes, isotopes) = IsotopomerizedSchemes(sch, isotopes)
# isotopomerize(sch::ChemicalSchemes, isotopes::Tuple) = IsotopomerizedSchemes(sch, loss_elements(last(isotopes), first(isotopes)))
isotopomerize(sch::IsotopomerizedSchemes, isotopes) = IsotopomerizedSchemes(sch.parent, gain_elements(sch.isotopes, isotopes))
# isotopomerize(sch::IsotopomerizedSchemes, isotopes::Tuple) = IsotopomerizedSchemes(sch.parent, loss_elements!(gain_elements(sch.isotopes, last(isotopes), first(isotopes))))
isotopomerize(sch::T, isotopes) where {T<:AbstractScheme} = throw(ArgumentError("Cannot add isotopes information to $T."))

"""
    groupedisotopomerize(chemical::AbstractChemicalScheme, index, isotope, isotopes, abundance) -> AbstractChemicalScheme

Add grouped isotopomer information to `chemical`, using `index`, reference `isotope`, replacement sets in `isotopes`, and their `abundance` values.
"""
groupedisotopomerize(chemical::AbstractChemical, index, isotope, isotopes, abundance) = Groupedisotopomers(chemical, index, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(chemical::Isotopomers, index, isotope, isotopes, abundance) = Groupedisotopomers(chemical.parent, index, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(chemical::Groupedisotopomers, index, isotope, isotopes, abundance) = Groupedisotopomers(chemical.parent, index, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(sch::StructuralElementalScheme, index, isotope, isotopes, abundance) = StructuralElementalScheme(structuralscheme(sch), groupedisotopomerize(elementalscheme(sch), index, isotope, isotopes, abundance))
# groupedisotopomerize(sch::CompleteSchemeChemical, index, isotope, isotopes, abundance) = groupedisotopomerize(elementalscheme(sch), index, isotope, isotopes, abundance)
# groupedisotopomerize(sch::StructuralElementalScheme{T,<:AbstractChemical}, index, isotope, isotopes, abundance) where T = groupedisotopomerize(elementalscheme(sch), index, isotope, isotopes, abundance)
groupedisotopomerize(sch::ElementalScheme{true}, index, isotope, isotopes, abundance) = ElementalScheme(true, groupedisotopomerize(sch.chemical, index, isotope, isotopes, abundance))
groupedisotopomerize(sch::ElementalScheme{false}, index, isotope, isotopes, abundance) = ElementalScheme(false, groupedisotopomerize(sch.chemical, index, isotope, reverse_elements.(isotopes, true), abundance))
groupedisotopomerize(sch::ChemicalSchemes, index, isotope, isotopes, abundance) = Groupedisotopomerizedschemes(sch, index, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(sch::IsotopomerizedSchemes, index, isotope, isotopes, abundance) = Groupedisotopomerizedschemes(sch.parent, index, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(sch::Groupedisotopomerizedschemes, index, isotope, isotopes, abundance) = Groupedisotopomerizedschemes(sch.parent, index, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(sch::T, index, isotope, isotopes, abundance) where {T<:AbstractScheme} = throw(ArgumentError("Cannot add isotopes information to $T."))

groupedisotopomersisotopes(isotopes::Vector{ElementsVector}) = isotopes
groupedisotopomersisotopes(isotopes::Vector{Vector{Pair{String, Int64}}}) = [ElementsVector(first.(x), last.(x)) for x in isotopes]
groupedisotopomersisotopes(isotopes::Vector{Dict{String, Int64}}) = [ElementsVector(collect(keys(x)), collect(values((x)))) for x in isotopes]
groupedisotopomersisotopes(isotopes::Vector{Dictionary{String, Int64}}) = [ElementsVector(collect(keys(x)), collect(values((x)))) for x in isotopes]

"""
    isgainscheme(sch::AbstractScheme) -> Bool

Return `true` if `sch` contains only chemical gains; otherwise, return `false`.
"""
isgainscheme(sch) = false
isgainscheme(sch::ElementalScheme{true}) = true
isgainscheme(sch::AbstractCompleteScheme) = isgainscheme(elementalscheme(sch))
isgainscheme(sch::ChemicalSchemes) = all(isgainscheme, sch.schemes)
isgainscheme(sch::IsotopomerizedSchemes) = isgainscheme(sch.parent)
isgainscheme(sch::Groupedisotopomerizedschemes) = isgainscheme(sch.parent)

"""
    islossscheme(sch::AbstractScheme) -> Bool

Return `true` if `sch` contains only chemical losses; otherwise, return `false`.
"""
islossscheme(sch) = false
islossscheme(sch::ElementalScheme{false}) = true
islossscheme(sch::AbstractCompleteScheme) = islossscheme(elementalscheme(sch))
islossscheme(sch::ChemicalSchemes) = all(islossscheme, sch.schemes)
islossscheme(sch::IsotopomerizedSchemes) = islossscheme(sch.parent) 
islossscheme(sch::Groupedisotopomerizedschemes) = islossscheme(sch.parent) 

"""
    completescheme(precursor::AbstractChemical, product::AbstractChemical) -> StructuralElementalScheme
    completescheme(precursor::AbstractChemical, product::AbstractScheme) -> StructuralElementalScheme    
    completescheme(precursor::AbstractChemical, product::GenericChemical) -> StructuralElementalScheme
    completescheme(precursor::AbstractChemical, product::AdductIon{<:GenericChemical}) -> StructuralElementalScheme
    completescheme(precursor::AbstractChemical, product::CompleteSchemes) -> CompleteSchemes
    completescheme(precursor::AbstractChemical, product::ChemicalSchemes) -> ChemicalSchemes 
    completescheme(precursor::AbstractChemical, product::IsotopomerizedSchemes) -> IsotopomerizedSchemes
    completescheme(precursor::AbstractChemical, product::Groupedisotopomerizedschemes) -> Groupedisotopomerizedschemes
    completescheme(precursor::Nothing, product::AbstractChemical) -> StructuralElementalScheme

Complete `product` using information from `precursor`.

For generic precursors, this function calls `schemetransform_complete`, which first looks for `:_schemestructure` attribute, then looks for the entry matching precursor scheme (`nothing` representing no scheme) in the `:transformer` property, and lasly looks for the entry matching product scheme.

For other chemical types, it calls `elementalscheme(precursor, product)` to produce an elemental scheme that incorporates information from `precursor`.

# Examples
## Generic chemical
```julia-repl
julia> chemical = Chemical("18:0 PC-d9", "C44H79NO8PD9") # Deuterium-labeled PC on the methyl group
18:0 PC-d9

julia> loss_me = ElementalScheme(false, Chemical("Me", "CH3"; charge = 1))
Loss_Me

julia> loss_cd3 = ElementalScheme(false, Chemical("Me[D3]", "CD3"; charge = 1))
Loss_Me[D3]

julia> pushschemetransformer!(chemical, nothing, loss_me, loss_cd3);

julia> ionize(chemical, loss_me)
[(18:0 PC-d9)-Me[D3]]-
```
The following, more complex example shows how to define structure and scheme mappings for labeled chemicals.
```julia-repl
julia> ps1 = Chemical("PS[D3,13C3] 18:0/20:4", "C41[13C]3H75D3NO10P") # Deuterium- and carbon-13-labeled PS on the serine moiety
PS[D3,13C3] 18:0/20:4

julia> ps2 = Chemical("PS 18:0[D5]/20:4(5Z,8Z,11Z,14Z)", "C44H73D5NO10P") # Deuterium-labeled PS on the sn-1 fatty acid
PS 18:0[D5]/20:4(5Z,8Z,11Z,14Z)

julia> begin
       cserine = Chemical("Serine", "C3H5NO2"; abbreviation = "Ser") # Unlabeled serine
       cserinei = Chemical("Serine[D3,13C3]", "[13C]3H2D3NO2"; abbreviation = "Ser[D3,13C3]") # Labeled serine
       lossserine = ChemicalLoss(cserine)
       lossserinei = ChemicalLoss(cserinei)
       losshserine = ChemicalLoss(AdductIon(cserine, "[M+H]+"))
       losshserinei = ChemicalLoss(AdductIon(cserinei, "[M+H]+"))
       fa1 = Chemical("FA 18:0", "C18H36O2")
       fa1i = Chemical("FA 18:0[D5]", "C13D5H36O2") # Labeled fatty acid
       ifa1 = AdductIon(fa1, ChemicalLoss(Proton))
       ifa1i = AdductIon(fa1i, ChemicalLoss(Proton))
       end;

julia> pushschemetransformer!(ps1, nothing, [lossserine => lossserinei, losshserine => losshserinei]);

julia> pushschemetransformer!(ps1, ChemicalLoss(Proton), lossserine, lossserinei);

julia> pushschememerger!(ps1, ChemicalLoss(Proton), lossserine, losshserine);

julia> detectedchemical(ChemicalSeries(ionize(ps1, ChemicalLoss(Proton)) => lossserine))
[(PS[D3,13C3] 18:0/20:4)-Ser[D3,13C3]-H]-

julia> pushschemetransformer!(ps2, ChemicalLoss(Proton), ifa1, ifa1i);

julia> pushschemetransformer!(ps2, losshserine, ifa1 => ifa1i);
    
julia> pushschememerger!(ps2, ChemicalLoss(Proton), lossserine, losshserine);

julia> detectedchemical(ChemicalSeries(ionize(ps2, ChemicalLoss(Proton)) => ifa1))
[FA 18:0[D5]-H]-
```
## Custom chemical type
```julia-repl
julia> abstract type AbstractPC <: AbstractChemical end

julia> struct PC <: AbstractPC end # Unlabeled PC

julia> struct DLPC <: AbstractPC 
           location::Symbol
       end # Deuterium-labeled PC; location is :Me for the methyl group or another value for a different part

julia> struct DLMe <: AbstractChemical end; charge(::DLMe; loss = false, kwargs...) = loss ? -1 : 1; # Deuterium-labeled Methinium

julia> elementalscheme(pc::DLPC, ::CLType(Methenium)) = pc.location == :Me ? ChemicalLoss(DLMe) : ChemicalLoss(Methenium)
elementalscheme (generic function with 12 methods)

julia> ion1 = ionize(DLPC(:Me); adduct = "[M-Me]-")
[DLPC-DLMe]-

julia> detectedscheme(pc::AIType(AbstractPC, SESType(CGType(Acetate))), ::CLType(MethylAcetate)) = completescheme(ioncore(pc), ChemicalLoss(Methenium))
detectedscheme (generic function with 5 methods)

julia> ion2 = detectedchemical(ChemicalSeries(ionize(DLPC(:Me), "[M+OAc]-") => ChemicalLoss(MethylAcetate)))
[DLPC-DLMe]-

julia> ion1 == ion2
true
```
"""
completescheme(precursor::AbstractChemical, product::AbstractChemicalScheme) = completescheme_abstract(precursor, product)
completescheme(precursor::Nothing, product::AbstractChemical) = StructuralElementalScheme(RandomProductScheme(), product)
completescheme(precursor::GenericChemical, product::AbstractChemicalScheme) = completescheme_generic(precursor, product)
completescheme(precursor::AdductIon{<:GenericChemical}, product::AbstractChemicalScheme) = completescheme_generic(precursor, product)

completescheme_abstract(precursor::AbstractChemical, product::AbstractChemical) = StructuralElementalScheme(RandomProductScheme(), product)
completescheme_abstract(precursor::AbstractChemical, product::AbstractScheme) = StructuralElementalScheme(product, elementalscheme(precursor, product))
completescheme_abstract(precursor::AbstractChemical, product::CompleteSchemes) = product
# Direct construct; ensure schemes are complete
completescheme_abstract(precursor::AbstractChemical, product::ChemicalSchemes) = ChemicalSchemes(completescheme.(Ref(precursor), product.schemes), product.number)
completescheme_abstract(precursor::AbstractChemical, product::IsotopomerizedSchemes) = IsotopomerizedSchemes(completescheme(precursor, product.parent), product.isotopes)
completescheme_abstract(precursor::AbstractChemical, product::Groupedisotopomerizedschemes) = Groupedisotopomerizedschemes(completescheme(precursor, product.parent), product.index, product.isotope, product.isotopes, product.abundance)

# Chemical Wrappers are complete
completescheme_generic(precursor::GenericChemical, product::AbstractChemical) = schemetransform_complete(precursor, nothing, product)
completescheme_generic(precursor::GenericChemical, product::AbstractChemicalWrapper) = completescheme_abstract(precursor, product)
# Go for elemental scheme
completescheme_generic(precursor::GenericChemical, product::AbstractScheme) = completescheme_abstract(precursor, product)
completescheme_generic(precursor::GenericChemical, product::CompleteSchemes) = completescheme_abstract(precursor, product)
# Direct search; ensure schemes are complete
completescheme_generic(precursor::GenericChemical, product::AbstractSchemeWrapper) = schemetransform_complete(precursor, nothing, product)

completescheme_generic(precursor::AdductIon{<:GenericChemical}, product::AbstractChemical) = schemetransform_complete(ioncore(precursor), ionadduct(precursor), product)
completescheme_generic(precursor::AdductIon{<:GenericChemical}, product::AbstractChemicalWrapper) = completescheme_abstract(precursor, product)
completescheme_generic(precursor::AdductIon{<:GenericChemical}, product::AbstractScheme) = completescheme_abstract(precursor, product)
completescheme_generic(precursor::AdductIon{<:GenericChemical}, product::CompleteSchemes) = completescheme_abstract(precursor, product)
completescheme_generic(precursor::AdductIon{<:GenericChemical}, product::AbstractSchemeWrapper) = schemetransform_complete(ioncore(precursor), ionadduct(precursor), product)

"""
    completeschemechemical(precursor::AbstractChemicalScheme, product::AbstractChemicalScheme) -> AbstractChemicalScheme
    completeschemechemical(sch::CompleteSchemeChemical) -> AbstractChemical
    completeschemechemical(sch::AbstractScheme) -> AbstractScheme

First, call [`completescheme`](@ref) to generate `sch` from `precursor` and `product`. If `sch` is a `CompleteSchemeChemical`, return `elementalscheme(sch)`; otherwise, return `sch` unchanged.
"""
completeschemechemical(precursor, product) = completeschemechemical(completescheme(precursor, product))
completeschemechemical(sch::CompleteSchemeChemical) = elementalscheme(sch)
completeschemechemical(sch::AbstractScheme) = sch

"""
    elementalscheme(precursor::AbstractChemical, product::AbstractChemical) -> AbstractChemical
    elementalscheme(precursor::AbstractChemical, product::AbstractScheme) -> StructuralElementalScheme    
    elementalscheme(precursor::AbstractChemical, product::CompleteSchemes) -> CompleteSchemes
    elementalscheme(precursor::AbstractChemical, product::ChemicalSchemes) -> ChemicalSchemes 
    elementalscheme(precursor::AbstractChemical, product::IsotopomerizedSchemes) -> IsotopomerizedSchemes
    elementalscheme(precursor::AbstractChemical, product::Groupedisotopomerizedschemes) -> Groupedisotopomerizedschemes
    elementalscheme(precursor::AbstractChemical, product::AbstractElementalScheme) -> AbstractElementalScheme
    elementalscheme(precursor::Nothing, product::AbstractChemical) -> AbstractChemical

Return the elemental scheme for `product`, using `precursor` when needed.

For generic precursors, this function calls `schemetransform_elemental`, which first looks for `:_schemestructure` attribute, then looks for the entry matching precursor scheme (`nothing` representing no scheme) in the `:transformer` property, and lasly looks for the entry matching product scheme.

For a new structural scheme and chemical type, define the applicable method:
* `elementalscheme(::new_chemical_type, ::new_structural_type)` when the chemical type represents an ion analyzed by MS.
* `elementalscheme(::AdductIon{new_chemical_type, StructuralElementalScheme{structural_type}}, ::new_structural_type)` when the chemical type forms an adduct ion with a single scheme.
* `elementalscheme(::AdductIon{new_chemical_type, ChemicalSchemes}, ::new_structural_type)` when the chemical type forms an adduct ion with multiple schemes.

See [`completescheme`](@ref) for examples of defining new methods.
```
"""
elementalscheme(precursor::AbstractChemical, product::AbstractChemicalScheme) = elementalscheme_abstract(precursor, product)
elementalscheme(precursor::Nothing, product::AbstractChemical) = product
elementalscheme(precursor::GenericChemical, product::AbstractChemicalScheme) = elementalscheme_generic(precursor, product)
elementalscheme(precursor::AdductIon{<:GenericChemical}, product::AbstractChemicalScheme) = elementalscheme_generic(precursor, product)

elementalscheme_abstract(precursor::AbstractChemical, product::AbstractChemical) = product
elementalscheme_abstract(precursor::T, product::S) where {T<:AbstractChemical, S<:AbstractScheme} = throw(ArgumentError("Specific `elementalscheme(precursor::$T, scheme::$S)` method has to be implemented."))
elementalscheme_abstract(precursor::AbstractChemical, product::CompleteSchemes) = elementalscheme(product)
elementalscheme_abstract(precursor::AbstractChemical, product::ChemicalSchemes) = ChemicalSchemes(elementalscheme.(Ref(precursor), product.schemes), product.number)
elementalscheme_abstract(precursor::AbstractChemical, product::IsotopomerizedSchemes) = IsotopomerizedSchemes(elementalscheme(precursor, product.parent), product.isotopes)
elementalscheme_abstract(precursor::AbstractChemical, product::Groupedisotopomerizedschemes) = Groupedisotopomerizedschemes(elementalscheme(precursor, product.parent), product.index, product.isotope, product.isotopes, product.abundance)
elementalscheme_abstract(precursor::AbstractChemical, product::AbstractElementalScheme) = copy(product)
# elementalscheme(precursor::AbstractChemical, product::ElementalScheme{T}) where T = StructuralElementalScheme(product, ElementalScheme(T, product))

# structure search for generic types
elementalscheme_generic(precursor::GenericChemical, product::AbstractChemical) = schemetransform_elemental(precursor, nothing, product)
elementalscheme_generic(precursor::GenericChemical, product::AbstractChemicalWrapper) = product
elementalscheme_generic(precursor::GenericChemical, product::AbstractScheme) = schemetransform_elemental(precursor, nothing, product)
elementalscheme_generic(precursor::GenericChemical, product::CompleteSchemes) = elementalscheme(product)

elementalscheme_generic(precursor::AdductIon{<:GenericChemical}, product::AbstractChemical) = schemetransform_elemental(ioncore(precursor), ionadduct(precursor), product)
elementalscheme_generic(precursor::AdductIon{<:GenericChemical}, product::AbstractChemicalWrapper) = product
elementalscheme_generic(precursor::AdductIon{<:GenericChemical}, product::AbstractScheme) = schemetransform_elemental(ioncore(precursor), ionadduct(precursor), product)
elementalscheme_generic(precursor::AdductIon{<:GenericChemical}, product::CompleteSchemes) = elementalscheme(product)

"""
    detectedscheme(precursor::AdductIon, product::AbstractScheme) -> CompleteSchemes

Return a new scheme that merges `ionadduct(precursor)` with `product`.

For generic chemical types, this function calls `schememerge`, which first looks for `:_schemestructure` attribute, then looks for the entry matching precursor scheme (`nothing` representing no scheme) in the `:merger` property, and lasly looks for the entry matching product scheme.
It then calls [`schemetransform_complete`](@ref) and searches `nothing` and the matched scheme.

For other chemical types, it returns a complete scheme without incorporating information from `precursor`.

Defining a new method is optional for new structural schemes and chemical types unless their schemes need to be merged.
* `detectedscheme(::AdductIon{new_chemical_type, StructuralElementalScheme{structural_type}}, ::new_structural_type)`.
* `detectedscheme(::AdductIon{new_chemical_type, ChemicalSchemes}, ::new_structural_type)`.
* `detectedscheme(::AdductIon{new_chemical_type, StructuralElementalScheme{structural_type}}, ::ChemicalSchemes)`.
* `detectedscheme(::AdductIon{new_chemical_type, ChemicalSchemes}, ::ChemicalSchemes)`.

See [`completescheme`](@ref) for examples of defining new methods.
```
"""
detectedscheme(precursor::AdductIon, product::CompleteSchemes) = detectedscheme(precursor, structuralscheme(product))
detectedscheme(precursor::AdductIon, product::AbstractScheme) = ChemicalSchemes(ionadduct(precursor), completescheme(precursor, product))
# detectedscheme(precursor::AdductIon, product::CompleteSchemes) = StructuralElementalScheme(ChemicalSchemes(structuralscheme(ionadduct(precursor)), structuralscheme(product)), ChemicalSchemes(elementalscheme(ionadduct(precursor)), elementalscheme(product)))
# schemes search for generic adduction
function detectedscheme(precursor::AdductIon{<:GenericChemical}, product::CompleteSchemes) 
    sch = schememerge(ioncore(precursor), ionadduct(precursor), product)
    isnothing(sch) ? ChemicalSchemes(ionadduct(precursor), product) : schemetransform_complete(ioncore(precursor), nothing, sch)
end
function detectedscheme(precursor::AdductIon{<:GenericChemical}, product::AbstractScheme) 
    sch = schememerge(ioncore(precursor), ionadduct(precursor), product)
    isnothing(sch) ? ChemicalSchemes(ionadduct(precursor), completescheme(precursor, product)) : schemetransform_complete(ioncore(precursor), nothing, sch)
end

# For customized adduction type of specific chemical type, implement 
# detectedchemical(::chemicaltype, ::CompleteSchemes) -> adductiontype
# detectedchemical(::adductiontype, ::CompleteSchemes) -> adductiontype
# detectedscheme(::adductiontype, ::CompleteSchemes) -> CompleteSchemes
# completescheme(::chemicaltype, ::AbstractScheme) -> CompleteSchemes
# completescheme(::adductiontype, ::AbstractScheme) -> CompleteSchemes

# Internal interfaces for detectedchemical
"""
    detectedchemical(precursor::AbstractChemical, product::AbstractChemical) -> AbstractChemical 
    detectedchemical(precursor::AbstractChemical, sch::AbstractScheme) -> AbstractChemical 
    detectedchemical(precursor::AbstractChemical, sch::CompleteSchemeChemical) -> AbstractChemical 
    detectedchemical(precursor::AbstractChemical, sch::CompleteSchemes) -> AbstractChemical 
    detectedchemical(precursor::AbstractChemical, sch::StructuralChemicalScheme) -> AbstractChemical 
    detectedchemical(precursor::Nothing, product::AbstractChemicalScheme) -> AbstractChemical 

Return the chemical directly detected by MS.

By default, [`completescheme`](@ref) is applied first when `product` is a scheme. The completed scheme is then converted to the detected chemical, which may be an [`AdductIon`](@ref) or an elemental chemical, depending on its type.

When `precursor` is an [`Isotopomers`](@ref) or [`Groupedisotopomers`](@ref), the detected chemical is wrapped in the corresponding type.

For [`AdductIon`](@ref), [`detectedscheme`](@ref) combines the precursor adduct scheme with the product scheme.

You only need to define additional methods when using another subtype of [`AbstractAdductIon`](@ref):
* `detectedchemical(::new_adduction_type, ::CompleteSchemes)`.
* `detectedchemical(::new_adduction_type, ::AbstractScheme)`.
* Do not define `detectedchemical(::new_adduction_type, ::StructuralChemicalScheme)` specifically. Instead, define `elementalscheme(::new_adduction_type, ::structural_type)` for each `structural_type <: StructuralChemicalScheme`.
* Do not define `detectedchemical(::new_adduction_type, ::CompleteSchemeChemical)` specifically; its behavior depends only on the `elementalscheme` method for `sch`.
"""
detectedchemical(precursor::AbstractChemical, product::AbstractChemical) = detectedchemical(precursor, completescheme(precursor, product))
# detectedchemical(precursor::AbstractChemical, product::AbstractScheme) = detectedchemical(precursor, completescheme(precursor, product))
detectedchemical(precursor::AbstractChemical, product::CompleteSchemeChemical) = elementalscheme(product)
detectedchemical(precursor::AbstractChemical, product::CompleteSchemes) = AdductIon(precursor, product, 1) 
detectedchemical(precursor::AbstractChemical, product::AbstractScheme) = AdductIon(precursor, completescheme(precursor, product), 1) 
detectedchemical(precursor::AbstractChemical, product::StructuralChemicalScheme) = elementalscheme(precursor, product)

detectedchemical(precursor::Nothing, product::AbstractChemical) = product
detectedchemical(precursor::Nothing, product::CompleteSchemeChemical) = elementalscheme(product)
detectedchemical(precursor::Nothing, product::CompleteSchemes) = throw(ArgumentError("`detectedchemical` requires precursor for scheme product."))
detectedchemical(precursor::Nothing, product::AbstractScheme) = throw(ArgumentError("`detectedchemical` requires precursor for scheme product."))
detectedchemical(precursor::Nothing, product::StructuralChemicalScheme) = throw(ArgumentError("`detectedchemical` requires precursor for scheme product."))

detectedchemical(precursor::AbstractAdductIon, product::CompleteSchemeChemical) = elementalscheme(product)
detectedchemical(precursor::T, product::CompleteSchemes) where {T<:AbstractAdductIon} = throw(ArgumentError("Specific `detectedchemical(precursor::$T, scheme::CompleteSchemes)` method has to be implemented."))
detectedchemical(precursor::T, product::AbstractScheme) where {T<:AbstractAdductIon} = throw(ArgumentError("Specific `detectedchemical(precursor::$T, scheme::AbstractScheme)` method has to be implemented."))
detectedchemical(precursor::AbstractAdductIon, product::StructuralChemicalScheme) = elementalscheme(precursor, product)

detectedchemical(precursor::AdductIon, product::CompleteSchemeChemical) = elementalscheme(product)
detectedchemical(precursor::AdductIon, product::CompleteSchemes) = AdductIon(ioncore(precursor), detectedscheme(precursor, product), ncore(precursor))
detectedchemical(precursor::AdductIon, product::AbstractScheme) = AdductIon(ioncore(precursor), detectedscheme(precursor, product), ncore(precursor))
detectedchemical(precursor::AdductIon, product::StructuralChemicalScheme) = elementalscheme(precursor, product)

detectedchemical(precursor::Isotopomers, product::AbstractChemical) = Isotopomers(detectedchemical(chemicalparent(precursor), product), Pair{String, Int}[])
detectedchemical(precursor::Isotopomers, product::Isotopomers) = product
detectedchemical(precursor::Isotopomers, product::Groupedisotopomers) = Isotopomers(chemicalparent(product), isotopomersisotopes(product))
detectedchemical(precursor::Isotopomers, product::CompleteSchemeChemical) = detectedchemical(precursor, elementalscheme(product))
function detectedchemical(precursor::Isotopomers, product::CompleteSchemes)
    chemical = detectedchemical(chemicalparent(precursor), chemicalparent(product))
    isotopes = gain_elements(isotopomersisotopes(precursor), isotopomersisotopes(product))
    Isotopomers(chemical, isotopes)
end
function detectedchemical(precursor::Isotopomers, product::AbstractScheme)
    product = completescheme(precursor, product)
    chemical = detectedchemical(chemicalparent(precursor), chemicalparent(product))
    isotopes = gain_elements(isotopomersisotopes(precursor), isotopomersisotopes(product))
    Isotopomers(chemical, isotopes)
end
detectedchemical(precursor::Isotopomers, product::StructuralChemicalScheme) = detectedchemical(precursor, elementalscheme(chemicalparent(precursor), product))

detectedchemical(precursor::Groupedisotopomers, product::AbstractChemical) = Groupedisotopomers(detectedchemical(chemicalparent(precursor), product), 0, precursor.isotope, groupedisotopomersisotopes(product), groupedisotopomersabundance(product))
detectedchemical(precursor::Groupedisotopomers, product::Isotopomers) = Groupedisotopomers(chemicalparent(product), mass_shift_index(product; isotope = precursor.isotope), precursor.isotope, groupedisotopomersisotopes(product), groupedisotopomersabundance(product))
detectedchemical(precursor::Groupedisotopomers, product::Groupedisotopomers) = product
detectedchemical(precursor::Groupedisotopomers, product::CompleteSchemeChemical) = detectedchemical(precursor, elementalscheme(product))
function detectedchemical(precursor::Groupedisotopomers, product::CompleteSchemes)
    chemical = detectedchemical(chemicalparent(precursor), chemicalparent(product))
    isotopes = map((x, y) -> gain_elements(x, y), groupedisotopomersisotopes(precursor), groupedisotopomersisotopes(product))
    index = _mass_shift_index(first(isotopes), elements_mass()[precursor.isotope] - elements_mass()[elements_parents()[precursor.isotope]])
    Groupedisotopomers(chemical, index, precursor.isotope, groupedisotopomersisotopes(isotopes), groupedisotopomersabundance(product))
end
function detectedchemical(precursor::Groupedisotopomers, product::AbstractScheme)
    product = completescheme(precursor, product)
    chemical = detectedchemical(chemicalparent(precursor), chemicalparent(product))
    isotopes = map((x, y) -> gain_elements(x, y), groupedisotopomersisotopes(precursor), groupedisotopomersisotopes(product))
    index = _mass_shift_index(first(isotopes), elements_mass()[precursor.isotope] - elements_mass()[elements_parents()[precursor.isotope]])
    Groupedisotopomers(chemical, index, precursor.isotope, groupedisotopomersisotopes(isotopes), groupedisotopomersabundance(product))
end
detectedchemical(precursor::Groupedisotopomers, product::StructuralChemicalScheme) = detectedchemical(precursor, elementalscheme(chemicalparent(precursor), product))

chemicalentity(isobars::Isobars; kwargs...) = chemicalentity(first(chemicalspecies(isobars)))
chemicalentity(isotopomers::Groupedisotopomers; kwargs...) = Isotopomers(chemicalparent(isotopomers), isotopomersisotopes(isotopomers))
chemicalentity(ct::ChemicalTransition; kwargs...) = chemicalentity(first(chemicaltransition(ct)))

elementalscheme(sch::Groupedisotopomerizedschemes; kwargs...) = Groupedisotopomerizedschemes(elementalscheme(sch.parent; kwargs...), sch.index, sch.isotope, sch.isotopes, sch.abundance)
elementalscheme(sch::IsotopomerizedSchemes; kwargs...) = IsotopomerizedSchemes(elementalscheme(sch.parent; kwargs...), sch.isotopes)
elementalscheme(sch::ChemicalSchemes; kwargs...) = ChemicalSchemes(elementalscheme.(sch.schemes; kwargs...), sch.number)
structuralscheme(sch::Groupedisotopomerizedschemes; kwargs...) = Groupedisotopomerizedschemes(structuralscheme(sch.parent; kwargs...), sch.index, sch.isotope, sch.isotopes, sch.abundance)
structuralscheme(sch::IsotopomerizedSchemes; kwargs...) = IsotopomerizedSchemes(structuralscheme(sch.parent; kwargs...), sch.isotopes)
structuralscheme(sch::ChemicalSchemes; kwargs...) = ChemicalSchemes(structuralscheme.(sch.schemes; kwargs...), sch.number)
structuralscheme(::Nothing; kwargs...) = nothing 
structuralscheme(x::AbstractChemical; kwargs...) = x 
elementalscheme(::Nothing; kwargs...) = nothing 
elementalscheme(x::AbstractChemical; kwargs...) = x 

chemicalspecies(isobars::Isobars; kwargs...) = isobars.chemicals

function chemicaltransition(isobars::Isobars{<:ChemicalTransition}; kwargs...) 
    ct = chemicaltransition.(chemicalspecies(isobars))
    [Isobars(getindex.(ct, i), abundance) for (i, abundance) in enumerate(eachcol(isobars.abundance))]
end
chemicaltransition(ct::ChemicalTransition; kwargs...) = ct.transition

chemicalparent(isobars::Isobars; kwargs...) = chemicalparent(chemicalentity(isobars); kwargs...)
chemicalparent(isotopomers::Isotopomers; kwargs...) = isotopomers.parent 
chemicalparent(isotopomers::Groupedisotopomers; kwargs...) = isotopomers.parent 
chemicalparent(ct::ChemicalTransition; kwargs...) = ChemicalTransition(chemicalparent.(chemicaltransition(ct); kwargs...))

chemicalparent(sch::StructuralElementalScheme; kwargs...) = StructuralElementalScheme(structuralscheme(sch), chemicalparent(elementalscheme(sch); kwargs...))
chemicalparent(sch::ElementalScheme{T}; kwargs...) where T = ElementalScheme(T, chemicalparent(sch.chemical; kwargs...))
chemicalparent(sch::IsotopomerizedSchemes; kwargs...) = sch.parent
chemicalparent(sch::ChemicalSchemes; kwargs...) = ChemicalSchemes(chemicalparent.(sch.schemes; kwargs...), sch.number)
chemicalparent(sch::Groupedisotopomerizedschemes; kwargs...) = sch.parent 

inputchemical(isobars::Isobars; kwargs...) = Isobars([inputchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
inputchemical(ct::ChemicalTransition; kwargs...) = first(chemicaltransition(ct))

outputchemical(isobars::Isobars; kwargs...) = Isobars([outputchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
outputchemical(ct::ChemicalTransition; kwargs...) = last(chemicaltransition(ct))

analyzedchemical(isobars::Isobars; kwargs...) = detectedchemical(isobars; kwargs...)
analyzedchemical(isobars::Isobars{<:ChemicalTransition}; kwargs...) = Isobars([analyzedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
analyzedchemical(ct::ChemicalTransition; kwargs...) = detectedchemical(inputchemical(ct); kwargs...)

function seriesanalyzedchemical(isobars::Isobars{<:ChemicalTransition}; kwargs...) 
    ct = chemicaltransition.(chemicalspecies(isobars))
    [analyzedchemical(Isobars(getindex.(ct, i), abundance); kwargs...) for (i, abundance) in enumerate(eachcol(isobars.abundance))]
end
function seriesanalyzedchemical(ct::ChemicalTransition; kwargs...) 
    v = AbstractChemical[]
    precursor = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedchemical(c; precursor))
        precursor = last(v)
    end
    v
end
function seriesanalyzedisotopes(ct::ChemicalTransition; kwargs...)
    v = Vector{Pair{String, Int}}[]
    precursorisotopes = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedisotopes(c; precursorisotopes))
        precursorisotopes = last(v)
    end
    v
end
function seriesanalyzedcharge(ct::ChemicalTransition; kwargs...)
    v = Int[]
    precursorcharge = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedcharge(c; precursorcharge))
        precursorcharge = last(v)
    end
    v
end

detectedchemical(isobars::Isobars; kwargs...) = Isobars([detectedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance)
detectedchemical(isobars::Isobars{<:ChemicalTransition}; kwargs...) = Isobars([detectedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, end])

detectedisotopes(sch::CompleteSchemeChemical; precursor = nothing, precursorisotopes = nothing, kwargs...) = isotopomersisotopes(sch; kwargs...)
detectedcharge(sch::CompleteSchemeChemical; precursor = nothing, precursorcharge = nothing, kwargs...) = charge(sch; kwargs...)
detectedelements(sch::CompleteSchemeChemical; precursor = nothing, precursorelements = nothing, kwargs...) = chemicalelements(sch; kwargs...)

function detectedchemical(ct::ChemicalTransition; kwargs...) 
    precursor = nothing
    for c in chemicaltransition(ct) 
        precursor = detectedchemical(precursor, c)
    end
    precursor
end
detectedisotopes(ct::ChemicalTransition; kwargs...) = 
    detectedisotopes(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...)::Vector{Pair{String, Int}}
detectedisotopes(ct::AbstractVector; kwargs...) = 
    (length(ct) > 1 ? detectedisotopes(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) : isotopomersisotopes(first(ct); kwargs...))::Vector{Pair{String, Int}}
detectedcharge(ct::ChemicalTransition; kwargs...) = 
    detectedcharge(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...)::Int
detectedcharge(ct::AbstractVector; kwargs...) = 
    (length(ct) > 1 ? detectedcharge(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) : charge(first(ct); kwargs...))::Int

# GenericChemical attribute search
# structural -> complete, elemental -> complete
schemetransform_complete(chemical, precursor_schemes, product_schemes::CompleteSchemes) = product_schemes
schemetransform_elemental(chemical, precursor_schemes, product_schemes::IsotopomerizedSchemes) = IsotopomerizedSchemes(schemetransform_elemental(chemical, precursor_schemes, product_schemes.parent), product_schemes.isotopes)
schemetransform_complete(chemical, precursor_schemes, product_schemes::IsotopomerizedSchemes) = IsotopomerizedSchemes(schemetransform_complete(chemical, precursor_schemes, product_schemes.parent), product_schemes.isotopes)
schemetransform_elemental(chemical, precursor_schemes, product_schemes::ChemicalSchemes) = ChemicalSchemes([schemetransform_elemental(chemical, precursor_schemes, k) for k in product_schemes.schemes], product_schemes.number)
schemetransform_complete(chemical, precursor_schemes, product_schemes::ChemicalSchemes) = ChemicalSchemes([schemetransform_complete(chemical, precursor_schemes, k) for k in product_schemes.schemes], product_schemes.number)
function schemetransform_complete(chemical, precursor_schemes, product::AbstractChemical) 
    sch = schemetransform(chemical, precursor_schemes, product)
    isnothing(sch) ? StructuralElementalScheme(RandomProductScheme(), product) : StructuralElementalScheme(product, sch)
end
function schemetransform_elemental(chemical, precursor_schemes, product::AbstractChemical) 
    sch = schemetransform(chemical, precursor_schemes, product)
    isnothing(sch) ? product : sch
end
function schemetransform_complete(chemical, precursor_schemes, product_schemes::AbstractElementalScheme) 
    StructuralElementalScheme(product_schemes, schemetransform_elemental(chemical, precursor_schemes, product_schemes))
end

function schemetransform_elemental(chemical, precursor_schemes, product_schemes::AbstractElementalScheme) 
    sch = schemetransform(chemical, precursor_schemes, product_schemes)
    isnothing(sch) ? copy(product_schemes) : sch
end

"""
    schemestructure(chemical::AbstractChemicalScheme) -> Union{SchemeStructure, Nothing}

Get scheme structure.
"""
schemestructure(chemical::AbstractChemicalScheme) = getchemicalattr(chemical, :_schemestructure, nothing)

"""
    schemetransformer(chemical::AbstractChemicalScheme) -> Union{SchemeVector, Nothing}
    schemetransformer(chemical::AbstractChemicalScheme, precursor_schemes) -> Union{Pair, Nothing}

Show how schemes will be transformed for given `chemical` and `precursor_schemes`.
"""
function schemetransformer(chemical::AbstractChemicalScheme, precursor_schemes = missing)
    schemes = schemestructure(chemical)
    ismissing(precursor_schemes) && return isnothing(schemes) ? nothing : schemes.transformer
    isnothing(schemes) && return nothing 
    transformer = schemes.transformer
    isnothing(transformer) && return nothing 
    i = findfirst(x -> first(x) == structuralscheme(precursor_schemes), transformer.vector)
    isnothing(i) ? nothing : last(transformer.vector[i])
end

"""
    schemetransform(chemical::AbstractChemicalScheme, precursor_schemes, product_schemes) -> AbstractChemicalScheme

Transform `product_schemes` with `chemical` and `precursor_schemes`.
"""
function schemetransform(chemical::AbstractChemicalScheme, precursor_schemes, product_schemes)
    scheme = schemetransformer(chemical, precursor_schemes)
    isnothing(scheme) && return nothing
    i = findfirst(x -> first(x) == structuralscheme(product_schemes), scheme)
    isnothing(i) ? nothing : last(scheme[i])
end

# scheme -> complete

"""
    schememerger(chemical::AbstractChemicalScheme) -> Union{SchemeVector, Nothing}
    schememerger(chemical::AbstractChemicalScheme, precursor_schemes) -> Union{Pair, Nothing}

Show how schemes will be merged for given `chemical` and `precursor_schemes`.
"""
function schememerger(chemical::AbstractChemicalScheme, precursor_schemes = missing)
    schemes = schemestructure(chemical)
    ismissing(precursor_schemes) && return isnothing(schemes) ? nothing : schemes.merger
    isnothing(schemes) && return nothing 
    merger = schemes.merger
    isnothing(merger) && return nothing 
    i = findfirst(x -> first(x) == structuralscheme(precursor_schemes), merger.vector)
    isnothing(i) ? nothing : last(merger.vector[i])
end

"""
    schememerge(chemical::AbstractChemicalScheme, precursor_schemes, product_schemes) -> AbstractChemicalScheme

Merge `product_schemes` with `chemical` and `precursor_schemes`.
"""
function schememerge(chemical::AbstractChemicalScheme, precursor_schemes, product_schemes)
    scheme = schememerger(chemical, precursor_schemes)
    isnothing(scheme) && return nothing
    i = findfirst(x -> first(x) == structuralscheme(product_schemes), scheme)
    isnothing(i) ? nothing : last(scheme[i])
end

schememerge(chemical::AbstractChemicalScheme, precursor_schemes::Nothing, product_schemes::CompleteSchemes) = product_schemes
schememerge(chemical::AbstractChemicalScheme, precursor_schemes::Nothing, product_schemes) = completescheme(chemical, product_schemes)

"""
    pushschemestructure!(chemical::AbstractChemicalScheme, value::SchemeStructure)

Push `value` into `chemical` as a new attribute `:_schemestructure`. If `:_schemestructure` already exists, the `value` replaces the old one.
"""
pushschemestructure!(chemical::AbstractChemicalScheme, value::SchemeStructure) = 
    pushchemicalattr!(chemical, :_schemestructure, value)

"""
    pushschemetransformer!(chemical::AbstractChemicalScheme, value...)

Push `value` into `chemical` as transformer of the attribute `:_schemestructure`. 

`value` can be:
* A schemevector (`SchemeVetor`), or anything that can be converted into a [`SchemeVector`](@ref).
* A precursor scheme (`Union{AbstractChemicalScheme, Nothing}`) and product schemes pair(s). The pair can be any iteratable object.
* A precursor scheme (`Union{AbstractChemicalScheme, Nothing}`), input product scheme, and output product scheme.

If any values already exist, they are updated.
"""
pushschemetransformer!(chemical::AbstractChemicalScheme, value...) = pushschemestructure!(chemical, :transformer, value...)

"""
    pushschememerger!(chemical::AbstractChemicalScheme, value...)

Push `value` into `chemical` as merger of the attribute `:_schemestructure`. 

`value` can be:
* A schemevector (`SchemeVetor`), or anything that can be converted into a [`SchemeVector`](@ref).
* A with-chemical scheme (`Union{AbstractChemicalScheme, Nothing}`) and input-merged schemes pair(s). The pair can be any iteratable object.
* A with-chemical scheme (`Union{AbstractChemicalScheme, Nothing}`), input scheme, and merged scheme.

If any values already exist, they are updated.
"""
pushschememerger!(chemical::AbstractChemicalScheme, value...) = pushschemestructure!(chemical, :merger, value...)

pushschemestructure!(chemical::AbstractChemicalScheme, attr::AbstractString, value...) = pushschemestructure!(chemical, Symbol(attr), value...)
function pushschemestructure!(chemical::T, attr::Symbol, value...) where {T <: AbstractChemicalScheme}
    if hasproperty(T, :_schemestructure) 
        if isnothing(getproperty(chemical, :_schemestructure))
            return setproperty!(chemical, :_schemestructure, SchemeStructure(; (attr => SchemeVector(value...), )...))
        else
            return pushschemevector!(getproperty(chemical, :_schemestructure), attr, value...)
        end
    end
    attribute = getchemicalattr(chemical)
    isnothing(attribute) && return pushchemicalattr_new!(chemical, :_schemestructure, SchemeStructure(; (attr => SchemeVector(value...), )...))
    if isnothing(attribute.attribute)
        attribute.attribute = Pair{Symbol, Any}[:_schemestructure => SchemeStructure(; (attr => SchemeVector(value...), )...)]
    else
        id = findfirst(==(:_schemestructure) ∘ first, attribute.attribute)
        if isnothing(id)
            push!(attribute.attribute, :_schemestructure => SchemeStructure(; (attr => SchemeVector(value...), )...))
        else
            pushschemevector!(last(attribute.attribute[id]), attr, value...) 
        end
    end
end

function pushschemevector!(schemeattr::SchemeStructure, attr::Symbol, value...)
    schemevec = getproperty(schemeattr, attr)
    if isnothing(schemevec)
        setproperty!(schemeattr, attr, SchemeVector(value...))
    else
        pushschemevector!(schemevec, value...)
    end
end

function pushschemevector!(schemevec::SchemeVector, pre, value...)
    id = findfirst(==(pre) ∘ first, schemevec.vector)
    if isnothing(id)
        topush = pre => _scheme_vec_pair(value...)
        # if topush isa eltype(schemevec.vector)
            push!(schemevec.vector, topush)
            sort!(schemevec.vector; by = first)
        # else
            # schemevec.vector = vcat(schemevec.vector, topush)
        # end
    else
        schemes = last(schemevec.vector[id])
        topush = _scheme_pair(value...)
        id2 = findfirst(==(first(topush)) ∘ first, schemes)
        if isnothing(id2)
            # if topush isa eltype(schemes)
                push!(schemes, topush)
                sort!(schemes; by = first)
            # else
                # topush = pre => vcat(schemes, topush)
                # schemevec.vector = vcat(deleteat!(schemevec.vector, id), topush)
            # end
        # elseif topush isa eltype(schemes)
        else
            schemes[id2] = topush 
        # else
            # topush = pre => vcat(deleteat!(schemes, id2), topush)
            # schemevec.vector = vcat(deleteat!(schemevec.vector, id), topush)
        end
    end
end