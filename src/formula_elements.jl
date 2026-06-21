chemicalformula(cc::Chemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalformula(cc::FormulaChemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalformula(isobars::Isobars; kwargs...) = chemicalformula(chemicalentity(isobars); kwargs...)
function chemicalformula(x::Isotopomers; kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); loss = false))
    chemicalformula(isotopeelements(elements, x.isotopes); kwargs...)
end

function isotopeelements(elements, isotopes)
    for (k, v) in isotopes
        e = get(elements_parents(), k, k) 
        k == e && continue 
        v == 0 && continue
        elements[e] -= v 
        get!(elements, k, 0)
        elements[k] += v 
    end
    elements
end

function chemicalformula(x::Groupedisotopomers; kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); loss = false))
    chemicalformula(isotopeelements(elements, first(x.isotopes)); kwargs...)
end
chemicalformula(ct::ChemicalTransition; kwargs...) = chemicalformula(chemicalentity(ct); kwargs...)

chemicalformula(sch::AbstractCompleteScheme; kwargs...) = chemicalformula(elementalscheme(sch); kwargs...)
chemicalformula(sch::ElementalScheme{false}; loss = false, kwargs...) = chemicalformula(sch.chemical; loss = !loss, kwargs..., ischemical = false)
chemicalformula(sch::ElementalScheme{true}; loss = false, kwargs...) = chemicalformula(sch.chemical; loss, kwargs..., ischemical = false)
chemicalformula(x::ChemicalSchema; kwargs...) = chemicalformula(chemicalelements(x; loss = false); kwargs..., ischemical = false)
function chemicalformula(x::IsotopomerizedSchema; kwargs...)
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); loss = false))
    chemicalformula(isotopeelements(elements, x.isotopes); kwargs..., ischemical = false)
end
function chemicalformula(x::Groupedisotopomerizedschema; kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); loss = false))
    chemicalformula(isotopeelements(elements, first(x.isotopes)); kwargs..., ischemical = false)
end

function reverse_formula(x, ischemical, loss) 
    if isempty(x)
        x
    elseif ischemical 
        x 
    elseif starswith(x, r"[^+-]")
        loss ? string("-", replace(x, "+" => "-", "-" => "+")) : string("+", x)
    else
        loss ? replace(x, "+" => "-", "-" => "+") : x
    end
end
    
reverse_elements(x::Vector{<: Pair}, loss) = loss ? [k => -v for (k, v) in x] : x
reverse_elements(x::Dict, loss) = loss ? Dict(k => -v for (k, v) in x) : x
reverse_elements(x::Dictionary, loss) = loss ? Dictionary(keys(x), [-v for v in x]) : x

chemicalelements(cc::Chemical; loss = false, kwargs...) = reverse_elements(cc.elements, loss)
chemicalelements(cc::FormulaChemical; loss = false, kwargs...) = reverse_elements(cc.elements, loss)
chemicalelements(isobars::Isobars; kwargs...) = chemicalelements(chemicalentity(isobars); kwargs...)
function chemicalelements(x::Isotopomers; loss = false, kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs..., loss = false))
    reverse_elements(collect(pairs(isotopeelements(elements, x.isotopes))), loss)
end

function chemicalelements(x::Groupedisotopomers; loss = false, kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs..., loss = false))
    reverse_elements(collect(pairs(isotopeelements(elements, first(x.isotopes)))), loss)
end
chemicalelements(ct::ChemicalTransition; loss = false, kwargs...) = chemicalelements(chemicalentity(ct); loss, kwargs...)

chemicalelements(sch::AbstractCompleteScheme; kwargs...) = chemicalelements(elementalscheme(sch); kwargs...)
chemicalelements(sch::ElementalScheme{false}; loss = false, kwargs...) = chemicalelements(sch.chemical; loss = !loss, kwargs...) 
chemicalelements(sch::ElementalScheme{true}; loss = false, kwargs...) = chemicalelements(sch.chemical; loss, kwargs...) 
function chemicalelements(x::IsotopomerizedSchema; loss = false, kwargs...)
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs..., loss = false))
    reverse_elements(collect(pairs(isotopeelements(elements, x.isotopes))), loss)
end
chemicalelements(x::ChemicalSchema; kwargs...) = vcat((repeat(chemicalelements(k; kwargs...), v) for (k, v) in zip(x.schema, x.number))...)
function chemicalelements(x::Groupedisotopomerizedschema; loss = false, kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs..., loss = false))
    reverse_elements(collect(pairs(isotopeelements(elements, first(x.isotopes)))), loss)
end

isotopomersisotopes(isobars::Isobars; kwargs...) = isotopomersisotopes(chemicalentity(isobars); kwargs...)
isotopomersisotopes(isotopomers::Isotopomers; loss = false, kwargs...) = reverse_elements(isotopomers.isotopes, loss)
isotopomersisotopes(isotopomers::Groupedisotopomers; loss = false, kwargs...) = reverse_elements(first(isotopomers.isotopes), loss)
isotopomersisotopes(ct::ChemicalTransition; kwargs...) = isotopomersisotopes(chemicalentity(ct); kwargs...)

isotopomersisotopes(sch::ElementalScheme{true}; loss = false, kwargs...) = isotopomersisotopes(sch.chemical; loss, kwargs...)
isotopomersisotopes(sch::ElementalScheme{false}; loss = false, kwargs...) = isotopomersisotopes(sch.chemical; loss = !loss, kwargs...)
isotopomersisotopes(x::IsotopomerizedSchema; loss = false, kwargs...) = reverse_elements(x.isotopes, loss)
isotopomersisotopes(x::Groupedisotopomerizedschema; loss = false, kwargs...) = reverse_elements(first(x.isotopes), loss)

isotopomerstate(sch::ElementalScheme{true}; isotope_unit = nothing, isotope = "[13C]", loss = false, kwargs...) = _isotopomerstate(isotopomersisotopes(sch; loss = false), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss, kwargs..., ischemical = false)
isotopomerstate(sch::ElementalScheme{false}; isotope_unit = nothing, isotope = "[13C]", loss = false, kwargs...) = _isotopomerstate(isotopomersisotopes(sch; loss = false), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss = !loss, kwargs..., ischemical = false)

function _isotopomerstate(isotopes::Vector, isotope_unit; ischemical = true, loss = false)
    ds = 0
    if ischemical || !loss
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

groupedisotopomersisotopes(x::ElementalScheme{true}; loss = false, kwargs...) = groupedisotopomersisotopes(x.chemical; loss, kwargs...)
groupedisotopomersisotopes(x::ElementalScheme{false}; loss = false, kwargs...) = groupedisotopomersisotopes(x.chemical; loss = !loss, kwargs...)
groupedisotopomersisotopes(x::ChemicalSchema; loss = false, kwargs...) = Pair{String, Int}[]
groupedisotopomersisotopes(x::Groupedisotopomers; loss = false, kwargs...) = reverse_elements.(x.isotopes, loss)
groupedisotopomersisotopes(x::Groupedisotopomerizedschema; loss = false, kwargs...) = reverse_elements.(x.isotopes, loss)

groupedisotopomersabundance(x::ElementalScheme; kwargs...) = groupedisotopomersabundance(x.chemical; kwargs...)
groupedisotopomersabundance(x::ChemicalSchema; kwargs...) = [1.0]
groupedisotopomersabundance(x::Groupedisotopomers; kwargs...) = x.abundance
groupedisotopomersabundance(x::Groupedisotopomerizedschema; kwargs...) = x.abundance

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

Elements container with no duplicated element keys or zeros.
"""
unique_elements(elements::T) where T = unique_elements(T, elements)
unique_elements(::Type{Dictionary}, elements) = filter!(!=(0), dictionary_elements(Dictionary, elements))
unique_elements(::Type{Dict}, elements) = Dict(pairs(unique_elements(Dictionary, elements)))
unique_elements(::Type{<: Vector{<: Pair}}, elements) = collect(pairs(unique_elements(Dictionary, elements)))
unique_elements(::Type{Dict}, elements::Dict) = filter(x -> last(x) != 0, elements)
unique_elements(::Type{<: Vector{<: Pair}}, elements::Dict) = filter!(x -> last(x) != 0, collect(elements))
unique_elements(::Type{Dictionary}, elements::Dictionary) = filter(!=(0), elements)

sort_unique_elements(x) = sort!(unique_elements(x))

"""
    dictionary_elements([Dicttype = Dict], elements::Vector{<: Pair}) -> Dicttype
    dictionary_elements([Dicttype = Dict], elements::Dict) -> Dicttype
    dictionary_elements([Dicttype = Dict], elements::Dictionary) -> Dicttype

Create a dictionary from element-number pairs. As elements can be duplicated in the original vector, the new dictionary is convenient for updating elements number.
"""
dictionary_elements(elements) = dictionary_elements(Dict, elements)
dictionary_elements(::Type{Dict}, elements::Vector{<: Pair}) = Dict(pairs(dictionary_elements(Dictionary, elements)))
dictionary_elements(::Type{Dictionary}, elements::Vector{<: Pair}) = groupsum(first, last, elements)
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
