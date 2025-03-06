"""
    ischemicalequal(x, y)

Determine whether two adduct ions or isobars are chemically equivalent. The equality relies on both `isadductequal` of adducts and `ischemicalequal` of core chemicals.
"""
ischemicalequal(x::AbstractAdductIon, y::AbstractAdductIon) = isadductequal(ionadduct(x), ionadduct(y)) && ischemicalequal(ioncore(x), ioncore(y))

"""
    ioncore(adduct_ion)

Core chemical of `adduct_ion`.
"""
ioncore(adduct_ion::AdductIon) = adduct_ion.core

"""
    ionadduct(adduct_ion)

Adduct of `adduct_ion`.
"""
ionadduct(adduct_ion::AdductIon) = adduct_ion.adduct

"""
    kmer(adduct_ion)

The number of core chemical. For instance, 2 for "[2M+H]+".
"""
kmer(adduct_ion::AbstractAdductIon) = kmer(ionadduct(adduct_ion))

"""
    charge(adduct_ion)

The charge of `adduct_ion` (positive or negative). For instance, -1 for "[M-H]-", 2 for "[M+2H]2+". The default value for positive/negative ion is 1/-1.
"""
charge(adduct_ion::AbstractAdductIon) = kmer(adduct_ion) * charge(ioncore(adduct_ion)) + charge(ionadduct(adduct_ion))

"""
    ncharge(adduct_ion)

The number of charges of `adduct_ion`. For instance, 1 for "[M-H]-", 2 for "[M+2H]2+".
"""
ncharge(adduct_ion::AbstractAdductIon) = abs(charge(adduct_ion))

"""
    adductelements(adduct_ion)

The elements changed with adduct of `adduct_ion`. It contains the elements of adduct itself and isotopic labeling related to the adduct for the `adduct_ion` (`adductisotopes(adduct_ion)`). 
"""
adductelements(adduct_ion::AbstractAdductIon) = vcat(adductelements(ionadduct(adduct_ion)), adductisotopes(adduct_ion))

"""
    adductisotopes(adduct_ion)

The elements changed with adduct when the core chemical has isotopic labeling related to the adduct. 

For instance, [M-Me]- of D-labeled PC may turn out to be [M-CD3]- rather than [M-CH3]- if Ds are on the methyl group. In this case, `adductisotopes` of [M-Me]- of PC should be `["H" => 3, "D" => -3]`.
"""
adductisotopes(adduct_ion::AbstractAdductIon) = Pair{String, Int}[] # ex 1D: ["H" => 1, "D" => -1]

"""
    getchemicalattr(adduct_ion::AbstractAdductIon, attr::Symbol; kwargs...)
    getchemicalattr(adduct_ion::AbstractAdductIon, val_attr::Val{T}; kwargs...)

Get attribute (`attr`) from core chemical of `adduct_ion`. 
"""
getchemicalattr(adduct_ion::AbstractAdductIon, attr::Symbol; kwargs...) = getchemicalattr(adduct_ion, Val(attr); kwargs...)
getchemicalattr(adduct_ion::AbstractAdductIon, attr::Val; kwargs...) = getchemicalattr(ioncore(adduct_ion), attr; kwargs...)
getchemicalattr(adduct_ion::AbstractAdductIon, attr::Val{:elements}; kwargs...) = vcat(getchemicalattr(ioncore(adduct_ion), attr; kwargs...), adductelements(adduct_ion))
function getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:name}) 
    r = chemicalname(ioncore(adduct_ion))
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    s = replace(string(ionadduct(adduct_ion)), "M" => r, r"\d*[+-]$" => ""; count = 2)
    c = charge(adduct_ion)
    c == 0 ? s : string(s, abs(c) > 1 ? abs(c) : "", c > 0 ? "+" : "-") 
end
function getchemicalattr(adduct_ion::AbstractAdductIon, ::Val{:formula}; kwargs...)
    el = chemicalelements(chemicalformula(ioncore(adduct_ion)))
    nm = kmer(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    chemicalformula(add_elements!(unique_elements(el), adductelements(adduct_ion)))
end
