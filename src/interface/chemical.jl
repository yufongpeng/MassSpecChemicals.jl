"""
    parse_chemical(name, formula; kwargs...)
    parse_chemical(::Type{<: Chemical}, name::AbstractString, formula::AbstractString; kwargs...)

Parse chemical name and formula. If name is not given, it defaults to formula. The default type is `Chemical`.
"""
parse_chemical(name, formula; kwargs...) = parse_chemical(Chemical, name, formula; kwargs...)
parse_chemical(::Type{<: Chemical}, name::AbstractString, formula::AbstractString; kwargs...) = Chemical(name, formula, collect(kwargs))

"""
    ischemicalequal(x::AbstractChemical, y::AbstractChemical)

Determine whether two chemicals or isobars are chemically equivalent. It defaults to `isequal`.
"""
ischemicalequal(x::AbstractChemical, y::AbstractChemical) = isequal(x, y)
ischemicalequal(x::Isobars, y::Isobars) = all(ischemicalequal(a, b) for (a, b) in zip(x.chemicals, y.chemicals)) && all(isapprox(a, b) for (a, b) in zip(x.abundance, y.abundance))
ischemicalequal(x::Isobars, y::AbstractIon) = length(x) == 1 && ischemicalequal(abundantchemical(x), y)
ischemicalequal(x::AbstractIon, y::Isobars) = length(y) == 1 && ischemicalequal(x, abundantchemical(y))
ischemicalequal(x::AbstractIon, y::AbstractIon) = isionequal(x, y)
"""
    ischemicalequal(x::Chemical, y::Chemical)

Determine whether two chemicals are chemically equivalent. It first test if the names are equal than the elements composition.
"""
ischemicalequal(x::Chemical, y::Chemical) = isequal(chemicalname(x), chemicalname(y)) || isequal(unique_elements(parse_compound(transform_isotope_repr(chemicalformula(x)))), unique_elements(parse_compound(transform_isotope_repr(chemicalformula(y)))))

"""
    abundantchemical(m::AbstractChemical)
    abundantchemical(m::Isobars)

The most abundant chemical from a chemical (itself) or isobars. 
"""
abundantion(m::AbstractChemical) = m
abundantion(m::Isobars) = first(m.chemicals)

"""
    getchemicalattr(m::AbstractChemical, attr::Symbol; kwargs...)
    getchemicalattr(m::AbstractChemical, val_attr::Val{T}; kwargs...)

Get attribute (`attr`) from `m`. This function defaults to take `attr` as a property name, and return the property. If `attr` is not a property name, it returns `nothing`.

To define specific method for a concrete type `C <: AbstractChemical`, and an attribute (`attr`), use the following function signature:

`getchemicalattr(::C, ::Val{attr}; kwargs...)`
"""
getchemicalattr(m::AbstractChemical, attr::Symbol; kwargs...) = getchemicalattr(m, Val(attr); kwargs...)
getchemicalattr(m::AbstractChemical, val_attr::Val{T}; kwargs...) where T = hasproperty(m, T) ? getproperty(m, T) : nothing 

"""
    getchemicalattr(m::Chemical, ::Val{T}; kwargs...)
    getchemicalattr(m::Chemical, ::Val{:name}; kwargs...)
    getchemicalattr(m::Chemical, ::Val{:formula}; kwargs...) 

Get attribute (`attr`) from a `m`. For attributes other than `:name` and `:formula`, it iterates through `m.attr`. If no matched attribute name is found, it returns `nothing`.
"""
function getchemicalattr(m::Chemical, ::Val{T}; kwargs...) where T
    hasproperty(m, T) && return getproperty(m, T)
    for (p, v) in m.attr
        p == T && return v
    end
    return nothing
end
getchemicalattr(m::Chemical, ::Val{:name}; kwargs...) = m.name
getchemicalattr(m::Chemical, ::Val{:formula}; kwargs...) = m.formula

"""
    getchemicalattr(m::Isobars, ::Val{:name}; kwargs...)
    getchemicalattr(m::Isobars, ::Val{:formula}; kwargs...) 
    getchemicalattr(m::Isobars, ::Val{:rt}; kwargs...)
    getchemicalattr(m::Isobars, ::Val{:chemicals}; kwargs...) 
    getchemicalattr(m::Isobars, ::Val{:abundance}; kwargs...) 

Get attribute (`attr`) from a `m`. 
"""
getchemicalattr(m::Isobars, ::Val{:name}; verbose = true, kwargs...) = verbose ? string("Isobars[", join(chemicalname.(m.chemicals), ", "), "]") : string("Isobars[", chemicalname(first(m.chemicals)), ", …]")
getchemicalattr(m::Isobars, ::Val{:formula}; kwargs...) = chemicalformula.(m.chemicals)
getchemicalattr(m::Isobars, ::Val{:chemicals}; kwargs...) = m.chemicals
getchemicalattr(m::Isobars, ::Val{:abundance}; kwargs...) = m.abundance
getchemicalattr(m::Isobars, ::Val{:rt}; kwargs...) = mean(rt.(m.chemicals; kwargs...), weights(m.abundance))
getchemicalattr(m::Isobars, ::Val{:abbreviation}; verbose = true, kwargs...) = verbose ? string("Isobars[", join(chemicalabbr.(m.chemicals), ", "), "]") : string("Isobars[", chemicalabbr(first(m.chemicals)), ", …]")
getchemicalattr(m::Isobars, ::Val{:SMILES}; verbose = true, kwargs...) = verbose ? string("Isobars[", join(chemicalsmiles.(m.chemicals), ", "), "]") : string("Isobars[", chemicalsmiles(first(m.chemicals)), ", …]")

"""
    attrpairs(m::AbstractChemical)
    attrpairs(m::Chemical)

Optional attribute names and values pairs that are transferred when transforming `m` into a `Chemical`. Default attibutes include `:rt``, `:abbreviation`, `:SMILES`.
"""
attrpairs(m::AbstractChemical) = [getchemicalattr(m, attr) for attr in [:rt, :abbreviation, :SMILES]]
attrpairs(m::Isobars) = [getchemicalattr(m, attr) for attr in [:chemicals, :abundance, :rt, :abbreviation, :SMILES]]
attrpairs(m::Chemical) = m.attr
# nothing, nochange
# :a, append
# :n, number
# vector
"""
    chemicalvariants(chemical::AbstractChemical, formulas::Vector{String}; names = :a)
    chemicalvariants(chemical::AbstractChemical, elements::Vector{<: Union{<: Vector{<: Pair}, <: Dict}}; names = :a)

Create variants of `chemical` using new formulas or elements as a vector of chemicals. By default, it uses `Chemical` as chemical. The attributes for new chemicals are determined by `attrpairs`.

The keyword argument `name` controls how chemical variants are nammed.
* `nothing`: no change to name.
* `:a`: append the difference of elements at the end of name.
* `:n`: enumerate name at the end.
* `::Vector`: use the vector as new names.
"""
function chemicalvariants(chemical::AbstractChemical, formulas::Vector{String}; names = :a)
    if isnothing(names)
        [Chemical(chemicalname(chemical), f, collect(attrpairs(chemical))) for f in formulas]
    elseif names == :n
        el = unique_elements(parse_compound(transform_isotope_repr(chemicalformula(chemical))))
        i = findfirst(x -> ==(unique_elements(parse_compound(transform_isotope_repr(x))), el), formulas)
        if isnothing(i)
            newname = chemicalname(chemical) .* string.(1:length(formulas))
            [Chemical(n, f, collect(attrpairs(chemical))) for (n, f) in zip(newname, formulas)]
        else
            newname = chemicalname(chemical) .* string.(1:length(formulas) - 1)
            insert!(newname, i, chemicalname(chemical))
            [Chemical(n, f, collect(attrpairs(chemical))) for (n, f) in zip(newname, formulas)]
        end
    elseif names == :a
        n = 1
        oldname = chemicalname(chemical)
        oldformula = chemicalformula(chemical)
        el = unique!(first.(parse_compound(transform_isotope_repr(oldformula))))
        newname = map(formulas) do f
            eln = unique_elements(parse_compound(transform_isotope_repr(f)))
            filter!(x -> !in(first(x), el), eln)
            chemicalformula(eln; delim = ",")
        end
        for (i, x) in enumerate(newname)
            if isempty(x)
                if unique_elements(parse_compound(transform_isotope_repr(formulas[i]))) == el
                    newname[i] = oldname
                else
                    newname[i] = string(oldname, "[", n, "]")
                    n += 1
                end
            else
                newname[i] = string(oldname, "[", x, "]")
            end
        end
        [Chemical(n, f, collect(attrpairs(chemical))) for (n, f) in zip(newname, formulas)]
    else
        [Chemical(n, f, collect(attrpairs(chemical))) for (n, f) in zip(names, formulas)]
    end
end

function chemicalvariants(chemical::AbstractChemical, elements::Vector{<: Union{<: Vector{<: Pair}, <: Dict}}; names = :a)
    if isnothing(names)
        [Chemical(chemicalname(chemical), chemicalformula(e), collect(attrpairs(chemical))) for e in elements]
    elseif names == :n
        el = unique_elements(parse_compound(transform_isotope_repr(chemicalformula(chemical))))
        i = findfirst(e -> ==(unique_elements(e), el), elements)
        if isnothing(i)
            newname = chemicalname(chemical) .* string.(1:length(elements))
            [Chemical(n, chemicalformula(e), collect(attrpairs(chemical))) for (n, e) in zip(newname, elements)]
        else
            newname = chemicalname(chemical) .* string.(1:length(elements) - 1)
            insert!(newname, i, chemicalname(chemical))
            [Chemical(n, chemicalformula(e), collect(attrpairs(chemical))) for (n, e) in zip(newname, elements)]
        end
    elseif names == :a
        n = 1
        oldname = chemicalname(chemical)
        oldformula = chemicalformula(chemical)
        el = unique!(first.(parse_compound(transform_isotope_repr(oldformula))))
        newname = map(elements) do f
            eln = unique_elements(f)
            filter!(x -> !in(first(x), el), eln)
            chemicalformula(eln; delim = ",")
        end
        for (i, x) in enumerate(newname)
            if isempty(x)
                if unique_elements(elements[i]) == el
                    newname[i] = oldname
                else
                    newname[i] = string(oldname, "[",  n, "]")
                    n += 1
                end
            else
                newname[i] = string(oldname, "[", x, "]")
            end
        end
        [Chemical(n, chemicalformula(e), collect(attrpairs(chemical))) for (n, e) in zip(newname, elements)]
    else
        [Chemical(n, chemicalformula(e), collect(attrpairs(chemical))) for (n, e) in zip(names, elements)]
    end
end

in(x::AbstractChemical, y::Isobars) = any(i -> ischemicalequal(i, x), y)
length(m::Isobars) = length(m.chemicals)
length(x::AbstractChemical) = 1
Broadcast.broadcastable(x::AbstractChemical) = Ref(x)
