chemicalformula(cc::Chemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalformula(cc::FormulaChemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalformula(isobars::Isobars; kwargs...) = chemicalformula(chemicalentity(isobars); kwargs...)
function chemicalformula(isotopomers::Isotopomers; kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(isotopomers); kwargs...))
    isotopeformula(elements, isotopomers.isotopes)
end

function isotopeformula(elements, isotopes)
    for (k, v) in isotopes
        e = get(elements_parents(), k, k) 
        k == e && continue 
        v > 0 || continue
        elements[e] -= v 
        get!(elements, k, 0)
        elements[k] += v 
    end
    chemicalformula(elements)
end
chemicalformula(isotopomers::Groupedisotopomers; kwargs...) = chemicalformula(chemicalentity(isotopomers); kwargs...)
chemicalformula(ct::ChemicalTransition; kwargs...) = chemicalformula(chemicalentity(ct); kwargs...)

chemicalformula(sch::StructuralElementalScheme; kwargs...) = chemicalformula(chemicalentity(sch); kwargs...)
chemicalformula(loss::ElementalScheme{false}; kwargs...) = string("-", chemicalformula(chemicalentity(loss); kwargs...))
chemicalformula(gain::ElementalScheme{true}; kwargs...) = string("+", chemicalformula(chemicalentity(gain); kwargs...))
function chemicalformula(x::IsotopomerizedSchema; kwargs...)
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs...))
    isotopeformula(elements, x.isotopes)
end
chemicalformula(x::ChemicalSchema) = join((join(repeat(chemicalformula(k), v), "") for (k, v) in zip(x.schema, x.number)), "")

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
        v > 0 || continue
        elements[e] -= v 
        get!(elements, k, 0)
        elements[k] += v 
    end
    collect(pairs(elements))
end
chemicalelements(isotopomers::Groupedisotopomers; kwargs...) = chemicalelements(chemicalentity(isotopomers); kwargs...)
chemicalelements(ct::ChemicalTransition; kwargs...) = chemicalelements(chemicalentity(ct); kwargs...)

chemicalelements(sch::AbstractCompleteScheme; kwargs...) = chemicalelements(elementalscheme(sch); kwargs...)
chemicalelements(sch::ElementalScheme{false}; loss = true, kwargs...) = loss ? chemicalelements(chemicalentity(sch); kwargs...) : [k => -v for (k, v) in chemicalelements(chemicalentity(sch); kwargs...)]
chemicalelements(sch::ElementalScheme{true}; loss = true, kwargs...) = loss ? [k => -v for (k, v) in chemicalelements(chemicalentity(sch); kwargs...)] : chemicalelements(chemicalentity(sch); kwargs...)
function chemicalelements(x::IsotopomerizedSchema; kwargs...)
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs...))
    isotopeelements(elements, x.isotopes)
end
chemicalelements(x::ChemicalSchema; kwargs...) = vcat((repeat(chemicalelements(k; kwargs...), v) for (k, v) in zip(x.schema, x.number))...)

"""
    chemicalformula(elements::Dict; delim = "", unique = true) -> String
    chemicalformula(elements::Dictionary; delim = "", unique = true) -> String
    chemicalformula(elements::Vector{<: Pair}; delim = "", unique = true) -> String

Create chemical formula using given element-number pairs. 

* `delim` assigns the delimiter between each element.
* `unique` determines whether combines the elements to become unique or not.
"""
function chemicalformula(elements::Vector{<: Pair}; delim = "", unique = true, loss = false)
    f = if unique 
        chemicalformula(dictionary_elements(Dictionary, elements); unique = false, delim)
    elseif all(x -> last(x) < 0, elements)
        string("-", join((v == -1 ? k : string(k, abs(v)) for (k, v) in elements if v != 0), delim))
    elseif any(x -> last(x) < 0, elements)
        chemicalformula(dictionary_elements(Dictionary, elements); unique = false, delim)
    else
        join((v == 1 ? k : string(k, v) for (k, v) in elements if v != 0), delim)
    end
end

function chemicalformula(elements::Dict; delim = "", unique = true)
    join((v == 1 ? k : string(k, v) for (k, v) in elements if v != 0), delim)
end

function chemicalformula(elements::Dictionary; delim = "", unique = true)
    join((v == 1 ? k : string(k, v) for (k, v) in pairs(elements) if v != 0), delim)
end

"""
    chemicalelements(formula::AbstractString) -> Vector{Pair{String, Int}}

Create element-number pairs from chemical formula. 
"""
function chemicalelements(formula::AbstractString; kwargs...)
    f = startswith(formula, r"[+-]") ? formula[begin + 1:end] : formula
    startswith(formula, "-") ? [elements_decodes()[k] => -v for (k, v) in parse_compound(encode_isotopes(f))] : 
    [elements_decodes()[k] => v for (k, v) in parse_compound(encode_isotopes(f))]
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
