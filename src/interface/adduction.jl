"""
    istransformedchemicalequal(x, y)

Determine whether two adduct ions or isobars are chemically equivalent. The equality relies on both `isadductequal` of adducts and `ischemicalequal` of core chemicals.
"""
istransformedchemicalequal(x::AbstractAdductIon, y::AbstractAdductIon) = isadductequal(ionadduct(x), ionadduct(y)) && ischemicalequal(ioncore(x), ioncore(y))

"""
    adductelements(adduct_ion)

The elements changed with adduct of `adduct_ion`. It contains the elements of adduct itself and isotopic labeling related to the adduct for the `adduct_ion` (`adductisotopes(adduct_ion)`). 
"""
adductelements(adduct_ion::AbstractAdductIon) = vcat(adductelements(ionadduct(adduct_ion)), adductisotopes(adduct_ion))

"""
    adductisotopes(adduct_ion)
    adductisotopes(adduct_ion::AbstractAdductIon{Chemical})

The elements changed when the core chemical has isotopic labeling that is lost in adduct formation. The returned vector is element-number pairs.

For instance, [M-Me]- (`Demethylation`) of Deuterium-labeled phosphatidylcholine (PC) may turn out to be [M-CD3]- rather than [M-CH3]- if Deuteriums are labeled on the methyl group of choline (`DLMC_PC`). In this case, `adductisotopes(::AbstractAdductIon{Demethylation, DLMC_PC})` should return `["H" => 3, "D" => -3]`.

For `adduct_ion::AbstractAdductIon{Chemical}`, user can define an attribute `:adductisotopes` for `ioncore(chemical)`. The attribute should be ionadduct-(elements-number pairs) pairs. 
This function finds this attribute, and extracts the value of key `ionadduct(adduct_ion)`. If the attribute or the key does not exist, empty vector is returned. 
"""
adductisotopes(adduct_ion::AbstractAdductIon) = Pair{String, Int}[] # ex 1D: ["H" => 1, "D" => -1]
function adductisotopes(adduct_ion::AbstractAdductIon{Chemical})
    a = ionadduct(adduct_ion)
    v = getchemicalattr(ioncore(adduct_ion), :adductisotopes)
    isnothing(v) && return Pair{String, Int}[]
    i = findfirst(x -> first(x) == a, v)
    isnothing(i) ? Pair{String, Int}[] : convert(Vector{Pair{String, Int}}, last(v[i])) # Force convert to pair as attr does not restrict input type
end

"""
    getchemicalattr(adduct_ion::AbstractAdductIon, attr::Symbol; kwargs...)
    getchemicalattr(adduct_ion::AbstractAdductIon, val_attr::Val{T}; kwargs...)
    getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:elements}; kwargs...)
    getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:name}; kwargs...)
    getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:formula}; kwargs...)
    getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:kmer}; kwargs...)
    getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:charge}; kwargs...)
    getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:abundant_chemical}; kwargs...)

Get attribute (`attr`) from `adduct_ion`. By default, it return attributes of core chemical without specialized methods.
"""
getchemicalattr(adduct_ion::AbstractAdductIon, attr::Symbol; kwargs...) = getchemicalattr(adduct_ion, Val(attr); kwargs...)
getchemicalattr(adduct_ion::AbstractAdductIon, val_attr::Val; kwargs...) = getchemicalattr(ioncore(adduct_ion), val_attr; kwargs...)
function getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:elements}; kwargs...) 
    el = deepcopy(chemicalelements(ioncore(adduct_ion); kwargs...))
    nm = kmer(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    vcat(el, adductelements(adduct_ion))
end
function getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:name}; kwargs...) 
    r = chemicalname(ioncore(adduct_ion); kwargs...)
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    s = replace(string(ionadduct(adduct_ion)), "M" => r; count = 1)
    s = replace(s, r"\d*[+-]$" => ""; count = 1)
    c = charge(adduct_ion)
    c == 0 ? s : string(s, abs(c) > 1 ? abs(c) : "", c > 0 ? "+" : "-") 
end
function getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:formula}; kwargs...)
    el = deepcopy(chemicalelements(ioncore(adduct_ion); kwargs...))
    nm = kmer(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    chemicalformula(add_elements!(unique_elements(el), adductelements(adduct_ion)))
end
function getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:abbreviation}; kwargs...) 
    r = chemicalabbr(ioncore(adduct_ion); kwargs...)
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    s = replace(string(ionadduct(adduct_ion)), "M" => r, r"\d*[+-]$" => ""; count = 2)
    c = charge(adduct_ion)
    c == 0 ? s : string(s, abs(c) > 1 ? abs(c) : "", c > 0 ? "+" : "-") 
end
getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:kmer}; kwargs...) = kmer(ionadduct(adduct_ion)) 
getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:charge}; kwargs...) = kmer(adduct_ion) * charge(ioncore(adduct_ion); kwargs...) + charge(ionadduct(adduct_ion))
getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:abundant_chemical}; kwargs...) = adduct_ion

"""
    getchemicalattr(adduct_ion::AdductIon, ::Val{:core}; kwargs...)
    getchemicalattr(adduct_ion::AdductIon, ::Val{:adduct}; kwargs...)

Get attribute (`attr`) from `adduct_ion`. 
"""
getchemicalattr(adduct_ion::AdductIon, ::Val{:core}; kwargs...) = adduct_ion.core 
getchemicalattr(adduct_ion::AdductIon, ::Val{:adduct}; kwargs...) = adduct_ion.adduct 
