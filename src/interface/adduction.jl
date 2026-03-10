"""
    istransformedchemicalequal(x::AbstractAdductIon, y::AbstractAdductIon)

Determine whether two adduct ions or isobars are chemically equivalent. The equality relies on both `isadductequal` of adducts and `ischemicalequal` of core chemicals.
"""
istransformedchemicalequal(x::AbstractAdductIon, y::AbstractAdductIon) = isadductequal(ionadduct(x), ionadduct(y)) && ischemicalequal(ioncore(x), ioncore(y))

"""
    adductelements(adduct_ion::AbstractAdductIon)

The elements changed with adduct of `adduct_ion`. It contains the elements of adduct itself and isotopic labeling related to the adduct for the `adduct_ion` (`adductisotopes(adduct_ion)`). 
"""
adductelements(adduct_ion::AbstractAdductIon) = vcat(adductelements(ionadduct(adduct_ion)), adductisotopes(adduct_ion))

"""
    adductisotopes(adduct_ion::AbstractAdductIon)
    adductisotopes(adduct_ion::AbstractAdductIon{Chemical})

The elements changed when the core chemical has isotopic labeling that is lost in adduct formation. The returned vector is element-number pairs.

For instance, [M-Me]- (`Demethylation`) of Deuterium-labeled phosphatidylcholine (PC) may turn out to be [M-CD3]- rather than [M-CH3]- if Deuteriums are labeled on the methyl group of choline (`DLMC_PC`). In this case, `adductisotopes(::AbstractAdductIon{Demethylation, DLMC_PC})` should return `["H" => 3, "D" => -3]`.

For `adduct_ion::AbstractAdductIon{Chemical}`, user can define an additional property `:adductisotopes` for `ioncore(chemical)`. The property should be ionadduct-(elements-number pairs) pairs. 
This function finds this property, and extracts the value of key `ionadduct(adduct_ion)`. If the property or the key does not exist, empty vector is returned. 
"""
adductisotopes(adduct_ion::AbstractAdductIon) = Pair{String, Int}[] # ex 1D: ["H" => 1, "D" => -1]
function adductisotopes(adduct_ion::AbstractAdductIon{Chemical})
    a = ionadduct(adduct_ion)
    v = getchemicalproperty(ioncore(adduct_ion), :adductisotopes, Pair{String, Int}[])
    isempty(v) && return v
    i = findfirst(x -> first(x) == a, v)
    isnothing(i) ? Pair{String, Int}[] : convert(Vector{Pair{String, Int}}, last(v[i])) # Force convert to pair as property does not restrict input type
end

getchemicalproperty(adduct_ion::AbstractAdductIon, property::Symbol, default) = hasproperty(adduct_ion, property) ? getproperty(adduct_ion, property) : getchemicalproperty(ioncore(adduct_ion), property, default)

function chemicalelements(adduct_ion::AbstractAdductIon; kwargs...) 
    el = deepcopy(chemicalelements(ioncore(adduct_ion); kwargs...))
    nm = kmer(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    vcat(el, adductelements(adduct_ion))
end
function chemicalname(adduct_ion::AbstractAdductIon; corename = nothing, kwargs...) 
    r = isnothing(corename) ? chemicalname(ioncore(adduct_ion); kwargs...) : corename
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    if isnothing(corename)
        s = replace(string(ionadduct(adduct_ion)), "M" => r; count = 1)
        s = replace(s, r"\d*[+-]$" => ""; count = 1)
    else
        s = replace(string(ionadduct(adduct_ion)), r"\d*M" => r; count = 1)
        s = replace(s, r"\d*[+-]$" => ""; count = 1)
    end
    c = charge(adduct_ion)
    c == 0 ? s : string(s, abs(c) > 1 ? abs(c) : "", c > 0 ? "+" : "-") 
end
function chemicalformula(adduct_ion::AbstractAdductIon; kwargs...)
    el = deepcopy(chemicalelements(ioncore(adduct_ion); kwargs...))
    nm = kmer(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    chemicalformula(add_elements!(unique_elements(el), adductelements(adduct_ion)))
end
function chemicalabbr(adduct_ion::AbstractAdductIon; kwargs...) 
    r = chemicalabbr(ioncore(adduct_ion); kwargs...)
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    s = replace(string(ionadduct(adduct_ion)), "M" => r, r"\d*[+-]$" => ""; count = 2)
    c = charge(adduct_ion)
    c == 0 ? s : string(s, abs(c) > 1 ? abs(c) : "", c > 0 ? "+" : "-") 
end
chemicalsmiles(adduct_ion::AbstractAdductIon; kwargs...) = chemicalsmiles(ioncore(adduct_ion); kwargs...) 