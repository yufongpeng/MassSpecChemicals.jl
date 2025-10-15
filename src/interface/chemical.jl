"""
    parse_chemical(name; kwargs...)
    parse_chemical(::Type{<: Chemical}, name::AbstractString; kwargs...)

Parse chemical name and construct a chemical object. The default type is `Chemical`. `kwargs` are attributes.
"""
parse_chemical(name; kwargs...) = parse_chemical(Chemical, name; kwargs...)
function parse_chemical(::Type{Chemical}, name::AbstractString; kwargs...) 
    ks = Dict{Symbol, Any}(kwargs)
    f = get!(ks, :eleements, [])
    delete!(ks, :eleements)
    if isempty(f) 
        f = get!(ks, :formula, "")
        delete!(ks, :formula)
        f = chemicalelements(f)
    end
    Chemical(name, f; ks...)
end

"""
    ischemicalequal(x::AbstractChemical, y::AbstractChemical)

Determine whether two chemicals are chemically equivalent. By default, it transforms both chemicals by `ischemicalequaltransform` and compares them by `istransformedchemicalequal`.
"""
ischemicalequal(x::AbstractChemical, y::AbstractChemical) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))
ischemicalequal(x::Isobars, y::Isobars) = all(ischemicalequal(a, b) for (a, b) in zip(x.chemicals, y.chemicals)) && all(isapprox(a, b) for (a, b) in zip(x.abundance, y.abundance))
ischemicalequal(x::Isotopomers, y::Isotopomers) = ischemicalequal(x.parent, y.parent) && isequal(sort!(collect(unique_elements(x.isotopes))), sort!(collect(unique_elements(y.isotopes))))
ischemicalequal(x::ChemicalLoss, y::ChemicalLoss) = ischemicalequal(x.chemical, y.chemical)
ischemicalequal(x::ChemicalPair, y::ChemicalPair) = ischemicalequal(x.precursor, y.precursor) && ischemicalequal(x.product, y.product)

"""
    ischemicalequaltransform(x::AbstractChemical) 

Return an object for comparison with other chemicals by `istransformedchemicalequal`. 
"""
ischemicalequaltransform(x::AbstractChemical) = x 
ischemicalequaltransform(x::Isobars) = length(x) == 1 ? abundantchemical(x) : x
ischemicalequaltransform(x::Isotopomers) = isempty(unique_elements(x.isotopes)) ? x.parent : x 
ischemicalequaltransform(x::ChemicalLoss) = x 
ischemicalequaltransform(x::ChemicalPair) = x 

"""
    istransformedchemicalequal(x::AbstractChemical, y::AbstractChemical)
    istransformedchemicalequal(x::Chemical, y::Chemical)

Determine whether two chemicals are chemically equivalent after applying `ischemicalequaltransform`. For `Chemical`, It first test if the names are equal than the elements composition. For others, it defaults to `isequal`.
"""
istransformedchemicalequal(x::AbstractChemical, y::AbstractChemical) = isequal(x, y)
istransformedchemicalequal(x::Chemical, y::Chemical) = 
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
        chemicalformula(getchemicalattr(cc, Val(:elements); shallow = true); kwargs...)
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
    getchemicalattr(chemical::Chemical, ::Val{:formula}; unique = false, kwargs...) 
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
getchemicalattr(cc::Chemical, ::Val{:elements}; kwargs...) = cc.elements
getchemicalattr(cc::Chemical, ::Val{:formula}; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
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
getchemicalattr(isobars::Isobars, ::Val{:name}; verbose = true, kwargs...) = (length(isobars.chemicals) == 1 || verbose) ? string("Isobars[", join(chemicalname.(isobars.chemicals; kwargs...), ", "), "]") : string("Isobars[", chemicalname(first(isobars.chemicals; kwargs...)), ", …]")
getchemicalattr(isobars::Isobars, ::Val{:formula}; kwargs...) = chemicalformula.(isobars.chemicals; kwargs...)
getchemicalattr(isobars::Isobars, ::Val{:elements}; kwargs...) = chemicalelements.(isobars.chemicals; kwargs...)
getchemicalattr(isobars::Isobars, ::Val{:chemicals}; kwargs...) = isobars.chemicals
getchemicalattr(isobars::Isobars, ::Val{:abundance}; kwargs...) = isobars.abundance
getchemicalattr(isobars::Isobars, ::Val{:rt}; kwargs...) = mean(rt.(isobars.chemicals; kwargs...), weights(isobars.abundance))
getchemicalattr(isobars::Isobars, ::Val{:abbreviation}; verbose = true, kwargs...) = (length(isobars.chemicals) == 1 || verbose) ? string("Isobars[", join(chemicalabbr.(isobars.chemicals; kwargs...), ", "), "]") : string("Isobars[", chemicalabbr(first(isobars.chemicals; kwargs...)), ", …]")
getchemicalattr(isobars::Isobars, ::Val{:SMILES}; verbose = true, kwargs...) = (length(isobars.chemicals) == 1 || verbose) ? string("Isobars[", join(chemicalsmiles.(isobars.chemicals; kwargs...), ", "), "]") : string("Isobars[", chemicalsmiles(first(isobars.chemicals; kwargs...)), ", …]")
getchemicalattr(isobars::Isobars, ::Val{:charge}; kwargs...) = mean(charge.(isobars.chemicals; kwargs...), weights(isobars.abundance))
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
function getchemicalattr(isotopomers::Isotopomers, ::Val{:elements}; kwargs...) 
    d = unique_elements(chemicalelements(isotopomers.parent; kwargs...))
    for (k, v) in isotopomers.isotopes
        e = get(ELEMENTS[:PARENTS], k, k) 
        k == e && continue 
        d[e] -= v 
        get!(d, k, 0)
        d[k] += v 
    end
    [v for v in pairs(d)]
end
getchemicalattr(isotopomers::Isotopomers, ::Val{:parent}; kwargs...) = isotopomers.parent
getchemicalattr(isotopomers::Isotopomers, ::Val{:isotopes}; kwargs...) = isotopomers.isotopes
getchemicalattr(isotopomers::Isotopomers, ::Val{:rt}; kwargs...) = rt(isotopomers.parent; kwargs...)
getchemicalattr(isotopomers::Isotopomers, ::Val{:abbreviation}; kwargs...) = string(chemicalabbr(isotopomers.parent; kwargs...), isempty(unique_elements(isotopomers.isotopes)) ? "" : string("[", chemicalformula(isotopomers.isotopes; delim = ","), "]"))
getchemicalattr(isotopomers::Isotopomers, ::Val{:SMILES}; kwargs...) = chemicalsmiles(isotopomers.parent; kwargs...)
getchemicalattr(isotopomers::Isotopomers, ::Val{:charge}; kwargs...) = charge(isotopomers.parent; kwargs...)
getchemicalattr(isotopomers::Isotopomers, ::Val{:abundant_chemical}; kwargs...) = isotopomers

"""
    getchemicalattr(chemicalloss::ChemicalLoss, ::Val{:name}; kwargs...)
    getchemicalattr(chemicalloss::ChemicalLoss, ::Val{:formula}; kwargs...)
    getchemicalattr(chemicalloss::ChemicalLoss, ::Val{:elements}; kwargs...)
    getchemicalattr(chemicalloss::ChemicalLoss, ::Val{:chemical}; kwargs...)
    getchemicalattr(chemicalloss::ChemicalLoss, ::Val{:rt}; kwargs...)
    getchemicalattr(chemicalloss::ChemicalLoss, ::Val{:abbreviation}; kwargs...)
    getchemicalattr(chemicalloss::ChemicalLoss, ::Val{:SMILES}; kwargs...)
    getchemicalattr(chemicalloss::ChemicalLoss, ::Val{:charge}; kwargs...)


Get attribute (`attr`) from `chemicalpair`.
"""
getchemicalattr(loss::ChemicalLoss, ::Val{:name}; kwargs...) = string("Loss_", chemicalname(loss.chemical; kwargs...))
getchemicalattr(loss::ChemicalLoss, ::Val{:formula}; kwargs...) = chemicalformula(loss.chemical; kwargs...)
getchemicalattr(loss::ChemicalLoss, ::Val{:elements}; kwargs...) = chemicalelements(loss.chemical; kwargs...) 
getchemicalattr(loss::ChemicalLoss, ::Val{:chemical}; kwargs...) = loss.chemical
getchemicalattr(loss::ChemicalLoss, ::Val{:rt}; kwargs...) = rt(loss.chemical; kwargs...)
getchemicalattr(loss::ChemicalLoss, ::Val{:abbreviation}; kwargs...) = string("Loss_", chemicalabbr(loss.chemical; kwargs...))
getchemicalattr(loss::ChemicalLoss, ::Val{:SMILES}; kwargs...) = chemicalsmiles(loss.chemical; kwargs...) 
getchemicalattr(loss::ChemicalLoss, ::Val{:charge}; kwargs...) = charge(loss.chemical; kwargs...) 

"""
    getchemicalattr(chemicalpair::ChemicalPair, ::Val{:name}; kwargs...)
    getchemicalattr(chemicalpair::ChemicalPair, ::Val{:formula}; kwargs...)
    getchemicalattr(chemicalpair::ChemicalPair, ::Val{:elements}; kwargs...)
    getchemicalattr(chemicalpair::ChemicalPair, ::Val{:precursor}; kwargs...)
    getchemicalattr(chemicalpair::ChemicalPair, ::Val{:product}; kwargs...)
    getchemicalattr(chemicalpair::ChemicalPair, ::Val{:rt}; kwargs...)
    getchemicalattr(chemicalpair::ChemicalPair, ::Val{:abbreviation}; kwargs...)
    getchemicalattr(chemicalpair::ChemicalPair, ::Val{:SMILES}; kwargs...)
    getchemicalattr(chemicalpair::ChemicalPair, ::Val{:charge}; kwargs...)


Get attribute (`attr`) from `chemicalpair`.
"""
getchemicalattr(cp::ChemicalPair, ::Val{:name}; kwargs...) = string(chemicalname(cp.precursor; kwargs...), " -> ", chemicalname(cp.product; kwargs...))
getchemicalattr(cp::ChemicalPair, ::Val{:formula}; kwargs...) = chemicalformula(cp.precursor; kwargs...) => chemicalformula(cp.product; kwargs...)
getchemicalattr(cp::ChemicalPair, ::Val{:elements}; kwargs...) = chemicalelements(cp.precursor; kwargs...) => chemicalelements(cp.product; kwargs...)
getchemicalattr(cp::ChemicalPair, ::Val{:precursor}; kwargs...) = cp.precursor
getchemicalattr(cp::ChemicalPair, ::Val{:product}; kwargs...) = cp.product
getchemicalattr(cp::ChemicalPair, ::Val{:rt}; kwargs...) = rt(cp.product; kwargs...)
getchemicalattr(cp::ChemicalPair, ::Val{:abbreviation}; kwargs...) = string(chemicalabbr(cp.precursor; kwargs...), " -> ", chemicalabbr(cp.product; kwargs...))
getchemicalattr(cp::ChemicalPair, ::Val{:SMILES}; kwargs...) = chemicalsmiles(cp.precursor; kwargs...) => chemicalsmiles(cp.product; kwargs...)
getchemicalattr(cp::ChemicalPair, ::Val{:charge}; kwargs...) = charge(cp.precursor; kwargs...) => charge(cp.product; kwargs...)


in(cc::AbstractChemical, isobars::Isobars) = any(i -> ischemicalequal(i, cc), isobars)
length(isobars::Isobars) = length(isobars.chemicals)
length(cc::AbstractChemical) = 1
Broadcast.broadcastable(cc::AbstractChemical) = Ref(cc)
