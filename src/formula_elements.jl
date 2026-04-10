chemicalformula(cc::Chemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalformula(cc::FormulaChemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalformula(isobars::Isobars; kwargs...) = chemicalformula(chemicalentity(isobars); kwargs...)
function chemicalformula(isotopomers::Isotopomers; kwargs...) 
    d = dictionary_elements(chemicalelements(chemicalparent(isotopomers); kwargs...))
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
chemicalformula(isotopomers::Groupedisotopomers; kwargs...) = chemicalformula(chemicalentity(isotopomers); kwargs...)
chemicalformula(loss::ChemicalLoss; kwargs...) = chemicalformula(chemicalentity(loss); kwargs...)
chemicalformula(gain::ChemicalGain; kwargs...) = chemicalformula(chemicalentity(gain); kwargs...)
chemicalformula(ct::ChemicalTransition; kwargs...) = chemicalformula(chemicalentity(ct); kwargs...)

chemicalelements(cc::Chemical; kwargs...) = cc.elements
chemicalelements(cc::FormulaChemical; kwargs...) = cc.elements
chemicalelements(isobars::Isobars; kwargs...) = chemicalelements(chemicalentity(isobars); kwargs...)
function chemicalelements(isotopomers::Isotopomers; kwargs...) 
    d = dictionary_elements(chemicalelements(chemicalparent(isotopomers); kwargs...))
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
chemicalelements(isotopomers::Groupedisotopomers; kwargs...) = chemicalelements(chemicalentity(isotopomers); kwargs...)
chemicalelements(loss::ChemicalLoss; kwargs...) = chemicalelements(chemicalentity(loss); kwargs...) 
chemicalelements(gain::ChemicalGain; kwargs...) = chemicalelements(chemicalentity(gain); kwargs...) 
chemicalelements(ct::ChemicalTransition; kwargs...) = chemicalelements(chemicalentity(ct); kwargs...)

"""
    chemicalformula(elements::Dictionary; delim = "", unique = true) -> String
    chemicalformula(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}; delim = "", unique = true) -> String

Create chemical formula using given element-number pairs. 

* `delim` assigns the delimiter between each element.
* `unique` determines whether combines the elements to become unique or not.
"""
function chemicalformula(elements::Vector{<: Pair}; delim = "", unique = true)
    if unique || any(x -> last(x) < 0, elements)
        chemicalformula(dictionary_elements(elements); unique = false, delim)
    else
        join((v == 1 ? k : string(k, v) for (k, v) in elements if v != 0), delim)
    end
end
function chemicalformula(elements::Dictionary; delim = "", unique = true)
    join((v == 1 ? k : string(k, v) for (k, v) in pairs(elements) if v != 0), delim)
end

"""
    chemicalelements(formula::AbstractString) -> Vector{Pair{String, Int}}

Create element-number pairs from chemical formula. 
"""
function chemicalelements(formula::AbstractString; kwargs...)
    formula = startswith(formula, r"[+-]") ? formula[begin + 1:end] : formula
    [elements_decodes()[k] => v for (k, v) in parse_compound(encode_isotopes(formula))]
end

"""
    unique_elements(elements::Vector{<: Pair}) -> Vector{<: Pair}
    unique_elements(elements::Dictionary) -> Dictionary

Elements container with no duplicated element keys.
"""
function unique_elements(elements::Vector{<: Pair})
    collect(pairs(filter!(!=(0), groupsum(first, last, elements))))
end
unique_elements(elements::Dictionary) = elements

"""
    unique_elements(elements::Vector{<: Pair}) -> Dictionary
    unique_elements(elements::Dictionary) -> Dictionary

Create a `Dictionary` from element-number pairs. As elements can be duplicated in the original vector, the new dictionary is convenient for updating elements number.
"""
dictionary_elements(elements::Vector{<: Pair}) = filter!(!=(0), groupsum(first, last, elements))
dictionary_elements(elements::Dictionary) = elements

"""
    gain_elements(elements::Vector{<: Pair}, y...) -> Vector{<: Pair}
    gain_elements(elements::Dictionary, y...) -> Dictionary

Add elements in `y` to deepcopied `elements`.
"""
gain_elements(elements, y...) = gain_elements!(deepcopy(elements), y...)
gain_elements(elements::Vector, y...) = vcat(elements, y...)

"""
    gain_elements!(elements::Vector{<: Pair}, y...) -> Vector{<: Pair}
    gain_elements!(elements::Dictionary, y...) -> Dictionary

Add elements in `y` to `elements`.
"""
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
    loss_elements(elements::Dictionary, y...) -> Dictionary

Substract elements in `y` from deepcopied `elements`.
"""
loss_elements(elements, y...) = loss_elements!(deepcopy(elements), y...)

"""
    loss_elements!(elements::Vector{<: Pair}, y...) -> Vector{<: Pair}
    loss_elements!(elements::Dictionary, y...) -> Dictionary

Substract elements in `y` from `elements`.
"""
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
        delta = isempty(m) ? 0 : (parse(Int, m) - round(Int, ustrip(elements_mass()[e])))
        e = delta > 0 ? string(e, "it") * "n" ^ delta :
            delta < 0 ? string(e, "it") * "p" ^ abs(delta) : string(e, "itz")
        f2 = replace(f2, i.match => e)
    end
    f2
end
