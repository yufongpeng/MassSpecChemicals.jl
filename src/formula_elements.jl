chemicalformula(cc::Chemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalformula(cc::FormulaChemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalformula(isobars::Isobars; kwargs...) = chemicalformula(chemicalentity(isobars); kwargs...)
function chemicalformula(isotopomers::Isotopomers; kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(isotopomers); kwargs...))
    isotopeformula(elements, isotopomers.isotopes)
end

function isotopeformula(elements, isotopes; kwargs...)
    for (k, v) in isotopes
        e = get(elements_parents(), k, k) 
        k == e && continue 
        elements[e] -= v 
        get!(elements, k, 0)
        elements[k] += v 
    end
    chemicalformula(elements; kwargs...)
end

chemicalformula(isotopomers::Groupedisotopomers; kwargs...) = chemicalformula(chemicalentity(isotopomers); kwargs...)
chemicalformula(ct::ChemicalTransition; kwargs...) = chemicalformula(chemicalentity(ct); kwargs...)

chemicalformula(sch::AbstractCompleteScheme; kwargs...) = chemicalformula(elementalscheme(sch); kwargs...)
chemicalformula(sch::ElementalScheme{false}; loss = true, kwargs...) = chemicalformula(sch.chemical; loss, kwargs..., ischemical = false)
chemicalformula(sch::ElementalScheme{true}; loss = true, kwargs...) = chemicalformula(sch.chemical; loss = !loss, kwargs..., ischemical = false)
chemicalformula(x::ChemicalSchema; loss = true, kwargs...) = chemicalformula(chemicalelements(x); loss, kwargs..., ischemical = false)
function chemicalformula(x::IsotopomerizedSchema; loss = true, kwargs...)
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x)))
    isotopeformula(elements, x.isotopes; loss, kwargs..., ischemical = false)
end
chemicalformula(x::Groupedisotopomerizedschema; kwargs...) = chemicalformula(IsotopomerizedSchema(chemicalparent(x), isotopomersisotopes(x)); kwargs...)

chemicalelements(cc::Chemical; kwargs...) = cc.elements
chemicalelements(cc::FormulaChemical; kwargs...) = cc.elements
chemicalelements(isobars::Isobars; kwargs...) = chemicalelements(chemicalentity(isobars); kwargs...)
function chemicalelements(isotopomers::Isotopomers; kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(isotopomers); kwargs...))
    isotopeelements(elements, isotopomers.isotopes)
end

function isotopeelements(elements, isotopes)
    for (k, v) in isotopes
        e = get(elements_parents(), k, k) 
        k == e && continue 
        elements[e] -= v 
        get!(elements, k, 0)
        elements[k] += v 
    end
    collect(pairs(elements))
end

chemicalelements(isotopomers::Groupedisotopomers; kwargs...) = chemicalelements(chemicalentity(isotopomers); kwargs...)
chemicalelements(ct::ChemicalTransition; kwargs...) = chemicalelements(chemicalentity(ct); kwargs...)

chemicalelements(sch::AbstractCompleteScheme; kwargs...) = chemicalelements(elementalscheme(sch); kwargs...)
chemicalelements(sch::ElementalScheme{false}; loss = true, kwargs...) = loss ? chemicalelements(sch.chemical; kwargs...) : [k => -v for (k, v) in chemicalelements(sch.chemical; kwargs...)]
chemicalelements(sch::ElementalScheme{true}; loss = true, kwargs...) = loss ? [k => -v for (k, v) in chemicalelements(sch.chemical; kwargs...)] : chemicalelements(sch.chemical; kwargs...)
function chemicalelements(x::IsotopomerizedSchema; kwargs...)
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs...))
    isotopeelements(elements, x.isotopes)
end
chemicalelements(x::ChemicalSchema; kwargs...) = vcat((repeat(chemicalelements(k; kwargs...), v) for (k, v) in zip(x.schema, x.number))...)
chemicalelements(x::Groupedisotopomerizedschema; kwargs...) = chemicalelements(IsotopomerizedSchema(chemicalparent(x), isotopomersisotopes(x)); kwargs...)

isotopomersisotopes(isobars::Isobars; kwargs...) = isotopomersisotopes(chemicalentity(isobars); kwargs...)
isotopomersisotopes(isotopomers::Isotopomers; kwargs...) = isotopomers.isotopes
isotopomersisotopes(isotopomers::Groupedisotopomers; kwargs...) = first(isotopomers.isotopes)
isotopomersisotopes(ct::ChemicalTransition; kwargs...) = isotopomersisotopes(chemicalentity(ct); kwargs...)

isotopomersisotopes(sch::ElementalScheme{true}; loss = true, kwargs...) = loss ? [k => -v for (k, v) in isotopomersisotopes(sch.chemical; kwargs...)] : isotopomersisotopes(sch.chemical; kwargs...)
isotopomersisotopes(sch::ElementalScheme{false}; loss = true, kwargs...) = loss ? isotopomersisotopes(sch.chemical; kwargs...) : [k => -v for (k, v) in isotopomersisotopes(sch.chemical; kwargs...)]
isotopomersisotopes(x::IsotopomerizedSchema; loss = true, kwargs...) = loss ? x.isotopes : [k => -v for (k, v) in x.isotopes]
isotopomersisotopes(x::ChemicalSchema; loss = true, kwargs...) = vcat((repeat(isotopomersisotopes(k; kwargs..., loss), v) for (k, v) in zip(x.schema, x.number))...)

isotopomerstate(sch::ElementalScheme{true}; isotope_unit = nothing, isotope = "[13C]", loss = true, kwargs...) = _isotopomerstate(isotopomersisotopes(sch), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss = !loss, ischemical = false, kwargs...)
isotopomerstate(sch::ElementalScheme{false}; isotope_unit = nothing, isotope = "[13C]", loss = true, kwargs...) = _isotopomerstate(isotopomersisotopes(sch), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss, ischemical = false, kwargs...)

function _isotopomerstate(isotopes::Vector, isotope_unit; ischemical = true, loss = false)
    ds = 0
    if ischemical || loss
        for (e, n) in isotopes
            ds += (elements_mass()[e] - elements_mass()[elements_parents()[e]]) * n
        end
    else
        for (e, n) in isotopes
            ds += (elements_mass()[e] - elements_mass()[elements_parents()[e]]) * (-n)
        end
    end
    round(Int, ds / isotope_unit)
end

"""
    chemicalformula(elements::Dict; delim = "", unique = true, ischemical = true, loss = false) -> String
    chemicalformula(elements::Dictionary; delim = "", unique = true, ischemical = true, loss = false) -> String
    chemicalformula(elements::Vector{<: Pair}; delim = "", unique = true, ischemical = true, loss = false) -> String

Create chemical formula using given element-number pairs. 

* `delim` assigns the delimiter between each element.
* `unique` determines whether combines the elements to become unique or not.
* `ischemical` determines whether the chemical is a chemical or a scheme. 
* `loss` determines whether the chemical is part of chemical loss, and signs are factored out from elements. 
"""
function chemicalformula(elements::Vector{<: Pair}; delim = "", unique = true, loss = false, ischemical = true)
    if unique 
        chemicalformula(dictionary_elements(Dictionary, elements); unique = false, loss, ischemical, delim)
    elseif all(x -> last(x) >= 0, elements)
        string(ischemical ? "" : loss ? "-" : "+", join((v == 1 ? k : string(k, v) for (k, v) in elements if v != 0), delim))
    elseif all(x -> last(x) <= 0, elements)
        string(loss ? "+" : "-", join((v == -1 ? k : string(k, abs(v)) for (k, v) in elements if v != 0), delim))
    elseif ischemical
        chemicalformula(dictionary_elements(Dictionary, elements); unique = false, loss, ischemical, delim)
    else
        elements_pos = filter(x -> last(x) > 0, elements)
        elements_neg = filter(x -> last(x) < 0, elements)
        string(loss ? "+" : "-", join((v == -1 ? k : string(k, abs(v)) for (k, v) in elements_neg), delim), loss ? "-" : "+", join((v == 1 ? k : string(k, v) for (k, v) in elements_pos), delim))
    end
end

chemicalformula(elements::Dict; delim = "", unique = true, loss = false, ischemical = true) = _chemicalformula(elements; delim, unique, loss, ischemical)
chemicalformula(elements::Dictionary; delim = "", unique = true, loss = false, ischemical = true) = _chemicalformula(pairs(elements); delim, unique, loss, ischemical)

function _chemicalformula(elements; delim = "", unique = false, loss = false, ischemical = true)
    if all(x -> last(x) >= 0, elements)
        string(ischemical ? "" : loss ? "-" : "+", join((v == 1 ? k : string(k, v) for (k, v) in elements if v != 0), delim))
    elseif all(x -> last(x) <= 0, elements)
        string(loss ? "+" : "-", join((v == -1 ? k : string(k, abs(v)) for (k, v) in elements if v != 0), delim))
    else
        elements_pos = filter(x -> last(x) > 0, elements)
        elements_neg = filter(x -> last(x) < 0, elements)
        string(loss ? "+" : "-", join((v == -1 ? k : string(k, abs(v)) for (k, v) in elements_neg), delim), loss ? "-" : "+", join((v == 1 ? k : string(k, v) for (k, v) in elements_pos), delim))
    end
end

"""
    chemicalelements(formula::AbstractString; loss = false) -> Vector{Pair{String, Int}}

Create element-number pairs from chemical formula. 

* `loss` determines whether the chemical is part of chemical loss, and sign flips are propagated into elements.
"""
function chemicalelements(formula::AbstractString; loss = false, kwargs...)
    fs = split(formula, "+")
    v = Vector{Pair{String, Int}}[]
    for f in fs 
        fns = split(f, "-")
        fp = popfirst!(fns)
        isempty(fp) || push!(v, [elements_decodes()[k] => v * (loss ? -1 : 1) for (k, v) in parse_compound(encode_isotopes(fp))])
        for fn in fns
            isempty(fn) || push!(v, [elements_decodes()[k] => v * (loss ? 1 : -1) for (k, v) in parse_compound(encode_isotopes(fn))])
        end
    end
    vcat(v...)
end

"""
    unique_elements(elements::Vector{<: Pair}) -> Vector{<: Pair}
    unique_elements(elements::Dict) -> Dict
    unique_elements(elements::Dictionary) -> Dictionary

Elements container with no duplicated element keys.
"""
function unique_elements(elements::Vector{<: Pair})
    collect(pairs(filter!(!=(0), groupsum(first, last, elements))))
end
unique_elements(elements::Dict) = elements
unique_elements(elements::Dictionary) = elements

"""
    dictionary_elements([Dicttype = Dict], elements::Vector{<: Pair}) -> Dicttype
    dictionary_elements([Dicttype = Dict], elements::Dict) -> Dicttype
    dictionary_elements([Dicttype = Dict], elements::Dictionary) -> Dicttype

Create a dictionary from element-number pairs. As elements can be duplicated in the original vector, the new dictionary is convenient for updating elements number.
"""
dictionary_elements(elements) = dictionary_elements(Dict, elements)
dictionary_elements(::Type{Dict}, elements::Vector{<: Pair}) = Dict(pairs(dictionary_elements(Dictionary, elements)))
dictionary_elements(::Type{Dictionary}, elements::Vector{<: Pair}) = filter!(!=(0), groupsum(first, last, elements))
dictionary_elements(::Type{Dict}, elements::Dict) = elements
dictionary_elements(::Type{Dictionary}, elements::Dict) = Dictionary(keys(elements), values(elements))
dictionary_elements(::Type{Dict}, elements::Dictionary) = Dict(pairs(elements))
dictionary_elements(::Type{Dictionary}, elements::Dictionary) = elements

"""
    gain_elements(elements::Vector{<: Pair}, y...) -> Vector{<: Pair}
    gain_elements(elements::Dict, y...) -> Dict
    gain_elements(elements::Dictionary, y...) -> Dictionary

Add elements in `y` to copied `elements`.
"""
gain_elements(elements, y...) = gain_elements!(copy(elements), y...)

"""
    gain_elements!(elements::Vector{<: Pair}, y...) -> Vector{<: Pair}
    gain_elements!(elements::Dict, y...) -> Dict
    gain_elements!(elements::Dictionary, y...) -> Dictionary

Add elements in `y` to `elements`.
"""
function gain_elements!(elements::Dict, y...) 
    for d in y 
        _gain_elements!(elements, d)
    end
    filter!(!=(0), elements)
end

function gain_elements!(elements::Dictionary, y...) 
    for d in y 
        _gain_elements!(elements, d)
    end
    filter!(!=(0), elements)
end

function gain_elements!(elements::Vector{<: Pair}, y...)
    for d in y 
        _gain_elements!(elements, d)
    end
    elements
end

_gain_elements!(elements, y::Dictionary) = __gain_elements!(elements, pairs(y))
_gain_elements!(elements, y) = __gain_elements!(elements, y)

function __gain_elements!(elements::Dict, y)
    for (k, v) in y
        get!(elements, k, 0)
        elements[k] += v
    end
end

function __gain_elements!(elements::Dictionary, y)
    for (k, v) in y
        get!(elements, k, 0)
        elements[k] += v
    end
end

function __gain_elements!(elements::Vector{<: Pair}, y)
    for k in y
        last(k) != 0 && push!(elements, k)
    end
end

"""
    loss_elements(elements::Vector{<: Pair}, y...) -> Vector{<: Pair}
    loss_elements(elements::Dict, y...) -> Dict
    loss_elements(elements::Dictionary, y...) -> Dictionary

Substract elements in `y` from copied `elements`.
"""
loss_elements(elements, y...) = loss_elements!(copy(elements), y...)

"""
    loss_elements!(elements::Vector{<: Pair}, y...) -> Vector{<: Pair}
    loss_elements!(elements::Dict, y...) -> Dict
    loss_elements!(elements::Dictionary, y...) -> Dictionary

Substract elements in `y` from `elements`.
"""
function loss_elements!(elements::Dict, y...) 
    for d in y 
        _loss_elements!(elements, d)
    end
    filter!(!=(0), elements)
end

function loss_elements!(elements::Dictionary, y...) 
    for d in y 
        _loss_elements!(elements, d)
    end
    filter!(!=(0), elements)
end

function loss_elements!(elements::Vector{<: Pair}, y...)
    for d in y 
        _loss_elements!(elements, d)
    end
    elements
end

_loss_elements!(elements, y::Dictionary) = __loss_elements!(elements, pairs(y))
_loss_elements!(elements, y) = __loss_elements!(elements, y)

function __loss_elements!(elements::Dict, y)
    for (k, v) in y
        get!(elements, k, 0)
        elements[k] -= v
    end
end

function __loss_elements!(elements::Dictionary, y)
    for (k, v) in y
        get!(elements, k, 0)
        elements[k] -= v
    end
end

function __loss_elements!(elements::Vector{<: Pair}, y)
    for (k, v) in y
        v != 0 && push!(elements, k => -v)
    end
end

function encode_isotopes(formula::AbstractString)
    f = string(formula)
    f2 = f
    for i in eachmatch(r"\[(\d*)([^\]]*)\]", f)
        m, e = i
        delta = isempty(m) ? 0 : (parse(Int, m) - round(Int, elements_mass()[e]))
        e = delta > 0 ? string(e, "it") * "n" ^ delta :
            delta < 0 ? string(e, "it") * "p" ^ abs(delta) : string(e, "itz")
        f2 = replace(f2, i.match => e)
    end
    f2
end
