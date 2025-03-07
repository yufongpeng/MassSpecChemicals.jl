"""
    parse_chemical(name, formula; kwargs...)
    parse_chemical(::Type{<: Chemical}, name::AbstractString, formula::AbstractString; kwargs...)

Parse chemical name and construct a chemical object. The default type is `Chemical`. `kwargs` are attributes.
"""
parse_chemical(name; kwargs...) = parse_chemical(Chemical, name; kwargs...)
function parse_chemical(::Type{Chemical}, name::AbstractString; kwargs...) 
    ks = Dict(kwargs)
    f = get!(ks, :formula, "")
    delete!(ks, :formula)
    Chemical(name, f; ks...)
end

"""
    ischemicalequal(x::AbstractChemical, y::AbstractChemical)

Determine whether two chemicals are chemically equivalent. It defaults to `isequal`.
"""
ischemicalequal(x::AbstractChemical, y::AbstractChemical) = isequal(x, y)
ischemicalequal(x::Isobars, y::Isobars) = all(ischemicalequal(a, b) for (a, b) in zip(x.chemicals, y.chemicals)) && all(isapprox(a, b) for (a, b) in zip(x.abundance, y.abundance))
ischemicalequal(x::Isobars, y::AbstractChemical) = length(x) == 1 && ischemicalequal(abundantchemical(x), y)
ischemicalequal(x::AbstractChemical, y::Isobars) = length(y) == 1 && ischemicalequal(x, abundantchemical(y))
ischemicalequal(x::Isotopomers, y::Isotopomers) = ischemicalequal(x.parent, y.parent) && isequal(sort!(collect(unique_elements(x.isotopes))), sort!(collect(unique_elements(y.isotopes))))
ischemicalequal(x::Isotopomers, y::AbstractChemical) = isempty(unique_elements(x.isotopes)) && ischemicalequal(x.parent, y)
ischemicalequal(x::AbstractChemical, y::Isotopomers) = isempty(unique_elements(y.isotopes)) && ischemicalequal(x, y.parent)
ischemicalequal(x::Isobars, y::Isotopomers) = length(x) == 1 && isempty(unique_elements(y.isotopes)) && ischemicalequal(abundantchemical(x), y.parent)
ischemicalequal(x::Isotopomers, y::Isobars) = length(y) == 1 && isempty(unique_elements(x.isotopes)) &&ischemicalequal(x.parent, abundantchemical(y))

"""
    ischemicalequal(x::Chemical, y::Chemical)

Determine whether two chemicals are chemically equivalent. It first test if the names are equal than the elements composition.
"""
ischemicalequal(x::Chemical, y::Chemical) = 
    isequal(chemicalname(x), chemicalname(y)) && isequal(sort!(collect(unique_elements(chemicalelements(x)))), sort!(collect(unique_elements(chemicalelements(y)))))

"""
    getchemicalattr(chemical::AbstractChemical, attr::Symbol; kwargs...)
    getchemicalattr(chemical::AbstractChemical, ::Val{T}; kwargs...)
    getchemicalattr(cc::AbstractChemical, ::Val{:formula}; shallow = false, kwargs...)
    getchemicalattr(cc::AbstractChemical, ::Val{:elements}; shallow = false, kwargs...)
    getchemicalattr(cc::AbstractChemical, ::Val{:charge}; kwargs...)
    getchemicalattr(cc::AbstractChemical, ::Val{:abundant_chemical}; kwargs...)

Get attribute (`attr`) from `chemical`. This function defaults to take `attr` as a property name, and return the property. If `attr` is not available, it returns `nothing`.

To define specific method for a concrete type `C <: AbstractChemical`, and an attribute (`attr`), use the following function signature:

`getchemicalattr(::C, ::Val{attr}; kwargs...)`
"""
getchemicalattr(cc::AbstractChemical, attr::Symbol; kwargs...) = getchemicalattr(cc, Val(attr); kwargs...)
getchemicalattr(cc::AbstractChemical, ::Val{T}; kwargs...) where T = hasproperty(cc, T) ? getproperty(cc, T) : nothing 
function getchemicalattr(cc::AbstractChemical, val_attr::Val{:formula}; shallow = false, kwargs...)
    if hasproperty(cc, :formula) 
        getproperty(cc, :formula)
    elseif shallow
        nothing
    else 
        chemicalformula(getchemicalattr(cc, Val(:elements); shallow = true))
    end
end
function getchemicalattr(cc::AbstractChemical, ::Val{:elements}; shallow = false, kwargs...)
    if hasproperty(cc, :elements) 
        getproperty(cc, :elements)
    elseif shallow
        nothing
    else 
        chemicalelements(getchemicalattr(cc, Val(:formula); shallow = true))
    end
end
getchemicalattr(cc::AbstractChemical, ::Val{:charge}; kwargs...) = 0
getchemicalattr(cc::AbstractChemical, ::Val{:abundant_chemical}; kwargs...) = cc

"""
    getchemicalattr(chemical::Chemical, ::Val{T}; kwargs...)
    getchemicalattr(chemical::Chemical, ::Val{:name}; kwargs...)
    getchemicalattr(chemical::Chemical, ::Val{:formula}; kwargs...) 
    getchemicalattr(chemical::Chemical, ::Val{:elements}; kwargs...)
    getchemicalattr(cc::Chemical, ::Val{:charge}; kwargs...)
    getchemicalattr(chemical::Chemical, ::Val{:abundant_chemical}; kwargs...)

Get attribute (`attr`) from `chemical`. For attributes other than `:name` and `:formula`, it iterates through `chemical.attr`. If no matched attribute name is found, it returns `nothing`.
"""
function getchemicalattr(cc::Chemical, ::Val{T}; kwargs...) where T
    hasproperty(cc, T) && return getproperty(cc, T)
    for (p, v) in cc.attr
        p == T && return v
    end
    return nothing
end
getchemicalattr(cc::Chemical, ::Val{:name}; kwargs...) = cc.name
getchemicalattr(cc::Chemical, ::Val{:formula}; kwargs...) = cc.formula
getchemicalattr(cc::Chemical, ::Val{:elements}; kwargs...) = chemicalelements(cc.formula; kwargs...)
getchemicalattr(cc::Chemical, ::Val{:charge}; kwargs...) = 0
getchemicalattr(cc::Chemical, ::Val{:abundant_chemical}; kwargs...) = cc

"""
    getchemicalattr(isobars::Isobars, ::Val{:name}; kwargs...)
    getchemicalattr(isobars::Isobars, ::Val{:formula}; kwargs...) 
    getchemicalattr(isobars::Isobars, ::Val{:elements}; kwargs...) 
    getchemicalattr(isobars::Isobars, ::Val{:chemicals}; kwargs...) 
    getchemicalattr(isobars::Isobars, ::Val{:abundance}; kwargs...) 
    getchemicalattr(isobars::Isobars, ::Val{:rt}; kwargs...)
    getchemicalattr(isobars::Isobars, ::Val{:abbreviation}; kwargs...) 
    getchemicalattr(isobars::Isobars, ::Val{:SMILES}; kwargs...) 
    getchemicalattr(isobars::Isobars, ::Val{:charge}; kwargs...) 
    getchemicalattr(isobars::Isobars, ::Val{:abundant_chemical}; kwargs...) 


Get attribute (`attr`) from `isobars`. 
"""
getchemicalattr(isobars::Isobars, ::Val{:name}; verbose = true, kwargs...) = verbose ? string("Isobars[", join(chemicalname.(isobars.chemicals; kwargs...), ", "), "]") : string("Isobars[", chemicalname(first(isobars.chemicals; kwargs...)), ", …]")
getchemicalattr(isobars::Isobars, ::Val{:formula}; kwargs...) = chemicalformula.(isobars.chemicals; kwargs...)
getchemicalattr(isobars::Isobars, ::Val{:elements}; kwargs...) = chemicalelements.(isobars.chemicals; kwargs...)
getchemicalattr(isobars::Isobars, ::Val{:chemicals}; kwargs...) = isobars.chemicals
getchemicalattr(isobars::Isobars, ::Val{:abundance}; kwargs...) = isobars.abundance
getchemicalattr(isobars::Isobars, ::Val{:rt}; kwargs...) = mean(rt.(isobars.chemicals; kwargs...), weights(isobars.abundance))
getchemicalattr(isobars::Isobars, ::Val{:abbreviation}; verbose = true, kwargs...) = verbose ? string("Isobars[", join(chemicalabbr.(isobars.chemicals; kwargs...), ", "), "]") : string("Isobars[", chemicalabbr(first(isobars.chemicals; kwargs...)), ", …]")
getchemicalattr(isobars::Isobars, ::Val{:SMILES}; verbose = true, kwargs...) = verbose ? string("Isobars[", join(chemicalsmiles.(isobars.chemicals; kwargs...), ", "), "]") : string("Isobars[", chemicalsmiles(first(isobars.chemicals; kwargs...)), ", …]")
getchemicalattr(isobars::Isobars, ::Val{:charge}; kwargs...) = mean(charge.(isobars.chemicala; kwargs...), weights(isobars.abundance))
getchemicalattr(isobars::Isobars, ::Val{:abundant_chemical}; kwargs...) = first(isobars.chemicals)

"""
    getchemicalattr(isotopomers::Isotopomers, ::Val{:name}; kwargs...)
    getchemicalattr(isotopomers::Isotopomers, ::Val{:formula}; kwargs...) 
    getchemicalattr(isotopomers::Isotopomers, ::Val{:elements}; kwargs...) 
    getchemicalattr(isotopomers::Isotopomers, ::Val{:parent}; kwargs...) 
    getchemicalattr(isotopomers::Isotopomers, ::Val{:isotopes}; kwargs...) 
    getchemicalattr(isotopomers::Isotopomers, ::Val{:rt}; kwargs...)
    getchemicalattr(isotopomers::Isotopomers, ::Val{:abbreviation}; kwargs...) 
    getchemicalattr(isotopomers::Isotopomers, ::Val{:SMILES}; kwargs...) 
    getchemicalattr(isotopomers::Isotopomers, ::Val{:charge}; kwargs...) 
    getchemicalattr(isotopomers::Isotopomers, ::Val{:abundant_chemical}; kwargs...) 

Get attribute (`attr`) from `isotopomers`. 
"""
getchemicalattr(isotopomers::Isotopomers, ::Val{:name}; kwargs...) = string(chemicalname(isotopomers.parent; kwargs...), isempty(unique_elements(isotopomers.isotopes)) ? "" : string("[", replace(chemicalformula(isotopomers.isotopes; delim = ","), "[" => "", "]" => ""), "]"))
function getchemicalattr(isotopomers::Isotopomers, ::Val{:formula}; kwargs...) 
    d = unique_elements(chemicalelements(isotopomers.parent; kwargs...))
    for (k, v) in isotopomers.isotopes
        e = get(ELEMENTS[:PARENTS], k, k) 
        k == e && continue 
        d[e] -= v 
        get!(d, k, 0)
        d[k] += v 
    end
    chemicalformula(d)
end
getchemicalattr(isotopomers::Isotopomers, ::Val{:elements}; kwargs...) = chemicalelements.(isotopomers.chemicals; kwargs...)
getchemicalattr(isotopomers::Isotopomers, ::Val{:parent}; kwargs...) = isotopomers.parent
getchemicalattr(isotopomers::Isotopomers, ::Val{:isotopes}; kwargs...) = isotopomers.isotopes
getchemicalattr(isotopomers::Isotopomers, ::Val{:rt}; kwargs...) = rt(isotopomers.parent; kwargs...)
getchemicalattr(isotopomers::Isotopomers, ::Val{:abbreviation}; kwargs...) = string(chemicalabbr(isotopomers.parent; kwargs...), "[", chemicalformula(isotopomers.isotopes; delim = ","), "]")
getchemicalattr(isotopomers::Isotopomers, ::Val{:SMILES}; kwargs...) = chemicalsmiles(isotopomers.parent; kwargs...)
getchemicalattr(isotopomers::Isotopomers, ::Val{:charge}; kwargs...) = charge(isotopomers.parent; kwargs...)
getchemicalattr(isotopomers::Isotopomers, ::Val{:abundant_chemical}; kwargs...) = isotopomers

in(cc::AbstractChemical, isobars::Isobars) = any(i -> ischemicalequal(i, cc), isobars)
length(isobars::Isobars) = length(isobars.chemicals)
length(cc::AbstractChemical) = 1
Broadcast.broadcastable(cc::AbstractChemical) = Ref(cc)
