"""
    ischemicalequal(x, y)

Determine whether two ions are chemically equivalent. The equality relies on both `isadductequal` of adducts and `ischemicalequal` of core chemicals.
"""
isionequal(x::AbstractIon, y::AbstractIon) = isadductequal(ionadduct(x), ionadduct(y)) && ischemicalequal(ioncore(x), ioncore(y))
"""
    ioncore(ion)

Core chemical of `ion`.
"""
ioncore(ion::Ion) = ion.core
"""
    ionadduct(ion)

Adduct of `ion`.
"""
ionadduct(ion::Ion) = ion.adduct
"""
    kmer(ion)

The number of core chemical. For instance, 2 for "[2M+H]+".
"""
kmer(ion::AbstractIon) = kmer(ionadduct(ion))
"""
    charge(ion)

The charge of `ion`. For instance, -1 for "[M-H]-", 2 for "[M+2H]2+". The default value for positive/negative ion is 1/-1.
"""
charge(ion::AbstractIon) = charge(ionadduct(ion))
"""
    ncharge(ion)

The number of charges of `ion`. For instance, 1 for "[M-H]-", 2 for "[M+2H]2+".
"""
ncharge(ion::AbstractIon) = ncharge(ionadduct(ion))
"""
    adductelement(ion)

The elements changed with adduct of `ion`. It contains the elements of adduct itself and isotopic labeling related to the adduct for the `ion` (`adductisotope(ion)`). 
"""
adductelement(ion::AbstractIon) = vcat(adductelement(ionadduct(ion)), adductisotope(ion))
"""
    adductisotope(ion)

The elements changed with adduct when the core chemical has isotopic labeling related to the adduct. 

For instance, [M-Me]- of D labeled PC may turn out to be [M-CD3]- rather than [M-CH3]- if Ds are on the methyl group. In this case, `adductisotope` of [M-Me]- of PC should be `["H" => 3, "D" => -3]`.
"""
adductisotope(ion::AbstractIon) = Pair{String, Int}[] # ex 1D: ["H" => 1, "D" => -1]
"""
    getchemicalattr(ion::AbstractIon, attr::Symbol; kwargs...)
    getchemicalattr(ion::AbstractIon, val_attr::Val{T}; kwargs...)

Get attribute (`attr`) from core chemical of `ion`. 
"""
getchemicalattr(ion::AbstractIon, attr::Symbol; kwargs...) = getchemicalattr(ion, Val(attr); kwargs...)
getchemicalattr(ion::AbstractIon, attr::Val; kwargs...) = getchemicalattr(ioncore(ion), attr; kwargs...)
getchemicalattr(ion::AbstractIon, ::Val{:name}) = string(chemicalname(ioncore(ion)), ionadduct(ion))
function getchemicalattr(ion::AbstractIon, ::Val{:formula}; kwargs...)
    el = parse_compound(transform_isotope_repr(chemicalformula(ioncore(ion))))
    nm = kmer(ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    ael = adductelement(ion)
    for e in ael
        id = findfirst(x -> first(x) == first(e), el)
        if isnothing(id)
            push!(el, e)
        else
            el[id] = first(el[id]) => (last(el[id]) + last(e))
        end
    end
    mapreduce(*, el) do e
        last(e) == 1 ? first(e) : string(first(e), last(e))
    end
end

"""
    ionvariants(ion::AbstractIon, formulas::Vector{String}; names = :a)
    ionvariants(ion::AbstractIon, elements::Vector{<: Union{<: Vector{<: Pair}, <: Dict}}; names = :a)

Create variants of `ion` using new formulas or elements (only core chemical part) as a vector of `Ion`. By default, it uses `Chemical` as core chemical. The attributes for new chemicals are determined by `attrpairs`.

The keyword argument `names` controls how chemical variants are nammed.
* `nothing`: no change to name.
* `:a`: append the difference of elements at the end of name.
* `:n`: enumerate name at the end.
* `::Vector`: use the vector as new names.
"""
ionvariants(ion::AbstractIon, x; names = :a) = Ion.(chemicalvariants(ioncore(ion), x; names), ionadduct(ion))