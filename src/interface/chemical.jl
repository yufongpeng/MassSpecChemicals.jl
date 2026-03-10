"""
    parse_chemical(name; kwargs...)
    parse_chemical(::Type{Chemical}, name::AbstractString; kwargs...)

Parse chemical name and construct a chemical object. The default type is `Chemical`. `kwargs` are properties.
"""
parse_chemical(name; kwargs...) = parse_chemical(Chemical, name; kwargs...)
function parse_chemical(::Type{Chemical}, name::AbstractString; kwargs...) 
    ks = Dict{Symbol, Any}(kwargs)
    f = get!(ks, :elements, [])
    delete!(ks, :elements)
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
ischemicalequaltransform(x::Isobars) = length(x) == 1 ? chemicalentity(x) : x
ischemicalequaltransform(x::Isotopomers) = isempty(unique_elements(x.isotopes)) ? x.parent : x 
ischemicalequaltransform(x::ChemicalLoss) = x 
ischemicalequaltransform(x::ChemicalPair) = x 

"""
    istransformedchemicalequal(x::AbstractChemical, y::AbstractChemical)
    istransformedchemicalequal(x::Chemical, y::Chemical)

Determine whether two chemicals are chemically equivalent after applying `ischemicalequaltransform`. For `Chemical` and `FormulaChemical`, It tests the name and the elements composition. For others, it defaults to `isequal`.
"""
istransformedchemicalequal(x::AbstractChemical, y::AbstractChemical) = isequal(x, y)
istransformedchemicalequal(x::Chemical, y::Chemical) = 
    isequal(chemicalname(x), chemicalname(y)) && isequal(sort!(collect(unique_elements(chemicalelements(x)))), sort!(collect(unique_elements(chemicalelements(y)))))
istransformedchemicalequal(x::FormulaChemical, y::FormulaChemical) = 
    isequal(chemicalname(x), chemicalname(y)) 

"""
    getchemicalproperty(chemical::AbstractChemical, property::Symbol, default = nothing)

Get property from `chemical`. 

This function defaults to finds the property, and returns `default` If it is not available.

For type `Chemical` and properties other than `name`, `elements` and `formula`, it iterates through `chemical.property`. If no matched property name is found, it returns `default`.

For type `AbstractAdductIon`, it searches for properties of itself and then the properties of core chemical without specialized methods.
"""
getchemicalproperty(chemical::AbstractChemical, property::Symbol, default = nothing) = hasproperty(chemical, property) ? getproperty(chemical, property) : default 

function getchemicalproperty(chemical::Union{Chemical, FormulaChemical}, property::Symbol, default)
    hasproperty(chemical, property) && return getproperty(chemical, property)
    for (p, v) in chemical.property
        p == property && return v
    end
    return default
end

# Chemical
chemicalformula(cc::Chemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalname(cc::Chemical; kwargs...) = cc.name
chemicalelements(cc::Chemical; kwargs...) = cc.elements

# FormulaChemical
chemicalformula(cc::FormulaChemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalname(cc::FormulaChemical; kwargs...) = chemicalformula(cc; kwargs...)
chemicalelements(cc::FormulaChemical; kwargs...) = cc.elements

# Isobars
chemicalname(isobars::Isobars; verbose = true, kwargs...) = (length(isobars.chemicals) == 1 || verbose) ? string("Isobars[", join(chemicalname.(isobars.chemicals; kwargs...), ", "), "]") : string("Isobars[", chemicalname(first(isobars.chemicals; kwargs...)), ", …]")
chemicalformula(isobars::Isobars; kwargs...) = chemicalformula(chemicalentity(isobars); kwargs...)
chemicalelements(isobars::Isobars; kwargs...) = chemicalelements(chemicalentity(isobars); kwargs...)
chemicalsmiles(isobars::Isobars; kwargs...) = chemicalsmiles(chemicalentity(isobars); kwargs...)
chemicalabbr(isobars::Isobars; verbose = true, kwargs...) = (length(isobars.chemicals) == 1 || verbose) ? string("Isobars[", join(chemicalabbr.(isobars.chemicals; kwargs...), ", "), "]") : string("Isobars[", chemicalabbr(first(isobars.chemicals; kwargs...)), ", …]")
charge(isobars::Isobars; kwargs...) = mean(charge.(isobars.chemicals; kwargs...), weights(isobars.abundance))
retentiontime(isobars::Isobars; kwargs...) = mean(retentiontime.(isobars.chemicals; kwargs...), weights(isobars.abundance))
chemicalentity(isobars::Isobars; kwargs...) = first(isobars.chemicals)
chemicalspecies(isobars::Isobars; kwargs...) = isobars.chemicals

# Isotopomers
chemicalname(isotopomers::Isotopomers; kwargs...) = string(chemicalname(isotopomers.parent; kwargs...), isempty(unique_elements(isotopomers.isotopes)) ? "" : string("[", replace(chemicalformula(isotopomers.isotopes; delim = ","), "[" => "", "]" => ""), "]"))
chemicalname(isotopomers::Isotopomers{FormulaChemical}; kwargs...) = string(chemicalformula(isotopomers))
function chemicalname(isotopomers::Isotopomers{<: AbstractAdductIon{<: FormulaChemical}}; kwargs...) 
    if isempty(adductelements(ionadduct(isotopomers.parent)))
        chemicalname(isotopomers.parent; corename = chemicalformula(isotopomers))
    else
        string(chemicalname(isotopomers.parent; kwargs...), isempty(unique_elements(isotopomers.isotopes)) ? "" : string("[", replace(chemicalformula(isotopomers.isotopes; delim = ","), "[" => "", "]" => ""), "]"))
    end
end
function chemicalformula(isotopomers::Isotopomers; kwargs...) 
    d = unique_elements(chemicalelements(isotopomers.parent; kwargs...))
    for (k, v) in isotopomers.isotopes
        e = get(elements_parents(), k, k) 
        k == e && continue 
        v > 0 || continue
        d[e] -= v 
        get!(d, k, 0)
        d[k] += v 
    end
    chemicalformula(d)
end
function chemicalelements(isotopomers::Isotopomers; kwargs...) 
    d = unique_elements(chemicalelements(isotopomers.parent; kwargs...))
    for (k, v) in isotopomers.isotopes
        e = get(elements_parents(), k, k) 
        k == e && continue 
        v > 0 || continue
        d[e] -= v 
        get!(d, k, 0)
        d[k] += v 
    end
    [v for v in pairs(d)]
end
chemicalsmiles(isotopomers::Isotopomers; kwargs...) = chemicalsmiles(isotopomers.parent; kwargs...)
retentiontime(isotopomers::Isotopomers; kwargs...) = retentiontime(isotopomers.parent; kwargs...)
chemicalabbr(isotopomers::Isotopomers; kwargs...) = string(chemicalabbr(isotopomers.parent; kwargs...), isempty(unique_elements(isotopomers.isotopes)) ? "" : string("[", chemicalformula(isotopomers.isotopes; delim = ","), "]"))
charge(isotopomers::Isotopomers; kwargs...) = charge(isotopomers.parent; kwargs...)
isotopomersisotopes(isotopomers::Isotopomers; kwargs...) = isotopomers.isotopes
chemicalparent(isotopomers::Isotopomers; kwargs...) = isotopomers.parent 

# ChemicalLoss
chemicalname(loss::ChemicalLoss; kwargs...) = string("Loss_", chemicalname(loss.chemical; kwargs...))
chemicalname(loss::ChemicalLoss{FormulaChemical}; kwargs...) = string("-", chemicalname(loss.chemical; kwargs...))
chemicalname(loss::ChemicalLoss{<: Isotopomers{FormulaChemical}}; kwargs...) = string("-", chemicalname(loss.chemical; kwargs...))
chemicalname(loss::ChemicalLoss{<: Isotopomers{<: AbstractAdductIon{<: FormulaChemical}}}; kwargs...) = string("-", chemicalname(loss.chemical; kwargs...))
chemicalformula(loss::ChemicalLoss; kwargs...) = chemicalformula(loss.chemical; kwargs...)
chemicalelements(loss::ChemicalLoss; kwargs...) = chemicalelements(loss.chemical; kwargs...) 
chemicalabbr(loss::ChemicalLoss; kwargs...) = string("Loss_", chemicalabbr(loss.chemical; kwargs...))
chemicalsmiles(loss::ChemicalLoss; kwargs...) = chemicalsmiles(loss.chemical; kwargs...) 
retentiontime(loss::ChemicalLoss; kwargs...) = retentiontime(loss.chemical; kwargs...)
charge(loss::ChemicalLoss; kwargs...) = charge(loss.chemical; kwargs...)
isotopomersisotopes(loss::ChemicalLoss; kwargs...) = isotopomersisotopes(loss.chemical; kwargs...) 
chemicalparent(loss::ChemicalLoss; kwargs...) = ChemicalLoss(chemicalparent(loss.chemical; kwargs...))
chemicalentity(loss::ChemicalLoss; kwargs...) = loss.chemical
analyzedchemical(loss::ChemicalLoss; kwargs...) = throw(ArgumentError("Chemical loss cannot be directly analyzed."))
function detectedchemical(loss::ChemicalLoss; precursor = nothing, kwargs...) 
    isnothing(precursor) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor."))
    precursor
end
function detectedisotopes(loss::ChemicalLoss; precursorisotopes = nothing, kwargs...) 
    isnothing(precursorisotopes) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor."))
    precursorisotopes
end
function detectedcharge(loss::ChemicalLoss; precursorcharge = nothing, kwargs...) 
    isnothing(precursorcharge) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor."))
    precursorcharge
end
function detectedelements(loss::ChemicalLoss; precursorelements = nothing, kwargs...) 
    isnothing(precursorelements) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor."))
    precursorelements
end
function detectedproduct(loss::ChemicalLoss; precursor = nothing, kwargs...) 
    isnothing(precursor) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor."))
    precursor = detectedchemical(precursor)
    _detectedproduct(precursor, loss)
end
function detectedproductisotopes(loss::ChemicalLoss; precursor = nothing, precursorisotopes = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorisotopes) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor information."))
    pre = isnothing(precursorisotopes) ? detectedisotopes(precursor; kwargs...) : precursorisotopes
    loss_elements(pre, isotopomersisotopes(loss; kwargs...))
end
function detectedproductcharge(loss::ChemicalLoss; precursor = nothing, precursorcharge = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorcharge) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor information."))
    (isnothing(precursorcharge) ? detectedcharge(precursor; kwargs...) : precursorcharge) - charge(loss; kwargs...)

end
function detectedproductelements(loss::ChemicalLoss; precursor = nothing, precursorelements = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorelements) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor information."))
    collect(pairs(loss_elements(isnothing(precursorelements) ? detectedelements(precursor; kwargs...) : precursorelements, chemicalelements(loss; kwargs...))))
end

chemicalname(cp::ChemicalPair; kwargs...) = string(chemicalname(cp.precursor; kwargs...), " -> ", chemicalname(cp.product; kwargs...))
chemicalabbr(cp::ChemicalPair; kwargs...) = string(chemicalabbr(cp.precursor; kwargs...), " -> ", chemicalabbr(cp.product; kwargs...))
chemicalformula(cp::ChemicalPair; kwargs...) = chemicalformula(cp.precursor; kwargs...)
chemicalelements(cp::ChemicalPair; kwargs...) = chemicalelements(cp.precursor; kwargs...)
chemicalsmiles(cp::ChemicalPair; kwargs...) = chemicalsmiles(cp.precursor; kwargs...)
retentiontime(cp::ChemicalPair; kwargs...) = retentiontime(cp.precursor; kwargs...)
charge(cp::ChemicalPair; kwargs...) = charge(cp.precursor; kwargs...) 
isotopomersisotopes(cp::ChemicalPair; kwargs...) = isotopomersisotopes(cp.precursor; kwargs...)

chemicalpair(cp::ChemicalPair; kwargs...) = cp.precursor => cp.product 
chemicalparent(cp::ChemicalPair; kwargs...) = ChemicalPair(chemicalparent(cp.precursor; kwargs...), chemicalparent(cp.product; kwargs...))
chemicalentity(cp::ChemicalPair; kwargs...) = chemicalentity(cp.precursor)
analyzedchemical(cp::ChemicalPair; kwargs...) = analyzedchemical(cp.precursor)
detectedchemical(cp::ChemicalPair; precursor = nothing, output = nothing, kwargs...) = 
    _detectedchemical(isnothing(output) ? outputchemical(cp.product) : output, cp; precursor, kwargs...) 
_detectedchemical(output, cp::ChemicalPair; kwargs...) = output 
_detectedchemical(output::ChemicalLoss, cp::ChemicalPair; kwargs...) = __detectedchemical(inputchemical(cp.product), cp; output, kwargs...)
__detectedchemical(post, cp::ChemicalPair; precursor = nothing, output = nothing, kwargs...) = detectedchemical(cp.product; output, kwargs...)
function __detectedchemical(post::ChemicalLoss, cp::ChemicalPair; precursor = nothing, output = nothing, kwargs...) 
    precursor = detectedchemical(cp.precursor; precursor, kwargs...)
    precursor = _detectedproduct(precursor, post)
    detectedchemical(cp.product; precursor, output, kwargs...)
end

_detectedproduct(precursor, product) = product 
function _detectedproduct(precursor, product::ChemicalLoss)
    prename = chemicalname(precursor)
    preelements = chemicalelements(precursor)
    precharge = charge(precursor)
    preisotopomersisotopes = isotopomersisotopes(precursor)
    postname = chemicalname(product)
    postelements = chemicalelements(product)
    postcharge = charge(product)
    postisotopomersisotopes = isotopomersisotopes(product)
    Chemical(string(prename, postname), collect(pairs(loss_elements(preelements, postelements))); charge = precharge - postcharge, isotopomersisotopes = loss_elements(preisotopomersisotopes, postisotopomersisotopes))
end

detectedisotopes(cp::ChemicalPair; precursorisotopes = nothing, output = nothing, kwargs...) = 
    _detectedisotopes(isnothing(output) ? outputchemical(cp.product) : output, cp; precursorisotopes, output, kwargs...)
_detectedisotopes(output, cp::ChemicalPair; precursorisotopes = nothing, kwargs...) = isotopomersisotopes(output; kwargs...)
_detectedisotopes(output::ChemicalLoss, cp::ChemicalPair; kwargs...) = __detectedisotopes(inputchemical(cp.product), cp; kwargs...)
__detectedisotopes(product, cp::ChemicalPair; precursorisotopes = nothing, output = nothing, kwargs...) = detectedisotopes(cp.product; output, kwargs...)
function __detectedisotopes(product::ChemicalLoss, cp::ChemicalPair; precursorisotopes = nothing, output = nothing, kwargs...)
    post = isotopomersisotopes(product; kwargs...)
    pre = detectedisotopes(cp.precursor; precursorisotopes, kwargs...) 
    precursorisotopes = loss_elements(pre, post)
    detectedisotopes(cp.product; precursorisotopes, output, kwargs...)
end

detectedcharge(cp::ChemicalPair; precursorcharge = nothing, output = nothing, kwargs...) = 
    _detectedcharge(isnothing(output) ? outputchemical(cp.product) : output, cp; precursorcharge, kwargs...)  
_detectedcharge(output, cp::ChemicalPair; precursorcharge = nothing, kwargs...) = charge(output; kwargs...)
_detectedcharge(output::ChemicalLoss, cp::ChemicalPair; kwargs...) = __detectedcharge(inputchemical(cp.product), cp; kwargs...)
__detectedcharge(product, cp::ChemicalPair; precursorcharge = nothing, output = nothing, kwargs...) = detectedcharge(cp.product; output, kwargs...)
function __detectedcharge(product::ChemicalLoss, cp::ChemicalPair; precursorcharge = nothing, output = nothing, kwargs...)
    post = charge(product; kwargs...)
    pre = detectedcharge(cp.precursor; precursorcharge, kwargs...)
    precursorcharge = pre - post
    detectedcharge(cp.product; precursorcharge, output, kwargs...)
end

detectedelements(cp::ChemicalPair; precursorelements = nothing, output = nothing, kwargs...) = 
    _detectedelements(isnothing(output) ? outputchemical(cp.product) : output, cp; precursorelements, output, kwargs...)
_detectedelements(output, cp::ChemicalPair; precursorelements = nothing, kwargs...) = chemicalelements(output; kwargs...)
_detectedelements(output::ChemicalLoss, cp::ChemicalPair; kwargs...) = __detectedelements(inputchemical(cp.product), cp; kwargs...)
__detectedelements(product, cp::ChemicalPair; precursorelements = nothing, output = nothing, kwargs...) = detectedelements(cp.product; output, kwargs...)
function __detectedelements(product::ChemicalLoss, cp::ChemicalPair; precursorelements = nothing, output = nothing, kwargs...)
    post = chemicalelements(product; kwargs...)
    pre = detectedelements(cp.precursor; precursorelements, kwargs...)
    precursorelements = collect(pairs(loss_elements(pre, post)))
    detectedelements(cp.product; precursorelements, output, kwargs...)
end

inputchemical(cp::ChemicalPair; kwargs...) = inputchemical(cp.precursor)
outputchemical(cp::ChemicalPair; kwargs...) = outputchemical(cp.product)

# truly formed chemical in a tandem MS
analyzedprecursor(cp::ChemicalPair; kwargs...) = detectedchemical(cp.precursor)
detectedproduct(cp::ChemicalPair; precursor = nothing, kwargs...) = 
    detectedproduct(inputchemical(cp.product); precursor = isnothing(precursor) ? cp.precursor : precursor, kwargs...)

detectedproductisotopes(cp::ChemicalPair; precursorisotopes = nothing, kwargs...) = 
    detectedproductisotopes(inputchemical(cp.product); precursor = cp.precursor, precursorisotopes, kwargs...)

detectedproductcharge(cp::ChemicalPair; precursorcharge = nothing, kwargs...) = 
    detectedproductcharge(inputchemical(cp.product); precursor = cp.precursor, precursorcharge, kwargs...)

detectedproductelements(cp::ChemicalPair; precursorelements = nothing, kwargs...) = 
    detectedproductelements(inputchemical(cp.product); precursor = cp.precursor, precursorelements, kwargs...)

# MS representation of chemical
inputprecursor(cp::ChemicalPair; kwargs...) = outputchemical(cp.precursor)
outputproduct(cp::ChemicalPair; kwargs...) = inputchemical(cp.product)

msstage(cp::ChemicalPair; n = 0) = msstage(cp.precursor; n = n + 1)

chemicaltransitions(cp::ChemicalPair) = _chemicaltransitions(cp)
function _chemicaltransitions(cp::ChemicalPair, products = AbstractChemical[])
    push!(products, cp.product)
    _chemicaltransitions(cp.precursor, products)
end

function _chemicaltransitions(chemical::AbstractChemical, products = [])
    precursors = [chemical]
    __chemicaltransitions(precursors, products)
end

function __chemicaltransitions(precursors, products)
    push!(precursors, _detectedproduct(last(precursors), pop!(products)))
    isempty(products) ? precursors : __chemicaltransitions(precursors, products)
end

in(cc::AbstractChemical, isobars::Isobars) = any(i -> ischemicalequal(i, cc), isobars)
length(isobars::Isobars) = length(isobars.chemicals)
length(cc::AbstractChemical) = 1
Broadcast.broadcastable(cc::AbstractChemical) = Ref(cc)
