"""
    Chemical <: AbstractChemical

Unstructured chemical type with its name, elements (formula), and additional properties.

# Fields 
* `name::String`: a unique chemical name.
* `elements::Vector{Pair{String, Int}}`: chemical elements.
* `property::Vector{Pair{Symbol, Any}}`: additional properties; the pairs repressent names and values.

# Constructors
* `Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}; kwargs_as_property_pairs...)`
* `Chemical(name::AbstractString, formula::AbstractString; kwargs_as_property_pairs...)`
"""
struct Chemical <: AbstractChemical
    name::String
    elements::Vector{Pair{String, Int}}
    property::Vector{Pair{Symbol, Any}}
end

Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}; kwargs...) = Chemical(name, elements, collect(kwargs))
Chemical(name::AbstractString, formula::AbstractString; kwargs...) = Chemical(name, chemicalelements(formula), collect(kwargs))

"""
    FormulaChemical <: AbstractChemical

Unstructured chemical type with elements (formula), and additional properties. Chemical name will be formula.

# Fields 
* `elements::Vector{Pair{String, Int}}`: chemical elements.
* `property::Vector{Pair{Symbol, Any}}`: additional properties; the pairs repressent names and values.

# Constructors
* `FormulaChemical(elements::Vector{Pair{String, Int}}; kwargs_as_property_pairs...)`
* `FormulaChemical(formula::AbstractString; kwargs_as_property_pairs...)`
"""
struct FormulaChemical <: AbstractChemical
    elements::Vector{Pair{String, Int}}
    property::Vector{Pair{Symbol, Any}}
end

FormulaChemical(elements::Vector{Pair{String, Int}}; kwargs...) = FormulaChemical(elements, collect(kwargs))
FormulaChemical(formula::AbstractString; kwargs...) = FormulaChemical(chemicalelements(formula), collect(kwargs))

"""
    Isobars{T <: AbstractChemical} <: AbstractChemical

Chemicals with similar m/z.

# Fields 
* `chemicals::Vector{T}`: a vector of chemicals.
* `abundnace::Vector{Float64}`: the abundance of each chemical.
"""
struct Isobars{T <: AbstractChemical, N} <: AbstractChemical
    chemicals::Vector{T}
    abundance::Vector{N}
end

"""
    Isotopomers{T <: AbstractChemical} <: AbstractChemical

Chemicals differed from isotopic replacement location.

# Fields 
* `parent::T`: shared chemical structure of isotopomers prior to isotopic replacement. 
* `isotopes::Vector{Pair{String, Int}}`: Isotopes-number pairs of isotopic replacement.

# Constructors
* `Isotopomers(parent::AbstractChemical, fullformula::String)`
* `Isotopomers(parent::AbstractChemical, fullelements::Dictionary)`
* `Isotopomers(parent::AbstractChemical, fullelements::Vector{Pair{String, Int}})`
All minor isotopes are regarded as isotopic replacement in `fullformula` and `fullelements`.
"""
struct Isotopomers{T <: AbstractChemical} <: AbstractChemical
    parent::T 
    isotopes::Vector{Pair{String, Int}}
end

function Isotopomers(parent::AbstractChemical, fullformula::String)
    Isotopomers(parent, unique_elements(chemicalelements(fullformula)))
end

function Isotopomers(parent::AbstractChemical, fullelements::Dictionary)
    dp = unique_elements(chemicalelements(parent))
    dr = deepcopy(fullelements)
    for k in keys(fullelements)
        haskey(elements_isotopes(), k) && (delete!(dr, k); continue)
        dr[k] -= get(dp, k, 0) 
    end
    Isotopomers(parent, [k => v for (k, v) in pairs(dr)])
end

"""
    ChemicalLoss{T <: AbstractChemical} <: AbstractChemical

Chemical loss from a precursor. This product is not detected in MS; the other part of precursor is detected instead.

# Fields 
* `chemical::T`
"""
struct ChemicalLoss{T <: AbstractChemical} <: AbstractChemical
    chemical::T 
end

"""
    ChemicalPair{T <: AbstractChemical, S <: AbstractChemical} <: AbstractChemical

A pair of precursor and product in MS/MS. Products can also be a `ChemicalLoss`.

# Fields 
* `precursor::T`
* `product::S`
"""
struct ChemicalPair{T <: AbstractChemical, S <: AbstractChemical} <: AbstractChemical
    precursor::T
    product::S
end

"""
    ChemicalSeries(chemical) -> AbstractChemical

Convert formula into `FormulaChemical`, and reconstruct the chemical pairs into `((c1 => c2) => c3)...` structure. This function do nothing for other valid chemical types.

For formula with adduct ion like structure, i.e. `"[formula+adduct]x+"` or `"[formula+adduct]x-"`, `FormulaChemical` will be wrapped by `AdductIon` with adduct (optional).
"""
ChemicalSeries(cc::AbstractChemical; kwargs...) = cc
ChemicalSeries(cp::ChemicalPair; kwargs...) = ChemicalPair(cp.precursor, cp.product)
ChemicalSeries(formula; charge = 0, loss = 0, kwargs...) = Formula2Chemical(formula, charge; loss)

ChemicalPair(cp1::ChemicalPair, cp2::ChemicalPair) = PostChemicalPair(ChemicalPair(cp1.precursor, cp1.product), cp2)
ChemicalPair(precursor::AbstractChemical, cp2::ChemicalPair) = PostChemicalPair(precursor, cp2)
PostChemicalPair(cc::AbstractChemical, cp::ChemicalPair) = PostChemicalPair(PostChemicalPair(cc, cp.precursor), cp.product)
PostChemicalPair(cc::AbstractChemical, cp::AbstractChemical) = ChemicalPair(cc, cp)

function Formula2Chemical(formula::Pair, net_charge = 0; default_charge = net_charge, loss = 0) 
    pre = Formula2Chemical(first(formula), net_charge; default_charge, loss)
    PostFormula2Chemical(pre, last(formula), detectedcharge(pre); default_charge, loss)
end
function parse_formula(formula::AbstractString; charge = 0, loss = 0) 
    cl, core = startswith(formula, "-") ? (true, formula[begin + 1:end]) : (false, formula)
    m = match(r"^\[(\d*)(.*)\](\d*[+-])*$", core)
    if isnothing(m) 
        n, f = match(r"^(\d*)(.*)$", core) 
        c = nothing
    else
        n, f, c = m 
    end
    n = isempty(n) ? 0 : parse(Int, n)
    if isnothing(c) && cl
        if abs(loss) == 0 
            ChemicalLoss(FormulaChemical(f; charge = loss)) 
        else 
            nm = n > 1 ? n : ""
            lm = abs(loss) > 1 ? abs(loss) : ""
            a = loss > 0 ? string("[", nm, "M]", lm , "+") : string("[", nm, "M]", lm, "-")
            ChemicalLoss(AdductIon(FormulaChemical(f), a))
        end
    elseif isnothing(c)
        if abs(charge) == 0 
            FormulaChemical(f; charge)
        else 
            nm = n > 1 ? n : ""
            lm = abs(loss) > 1 ? abs(charge) : ""
            a = charge > 0 ? string("[", nm, "M]", lm , "+") : string("[", nm, "M]", lm, "-")
            AdductIon(FormulaChemical(f), a)
        end
    else
        plusa = split(f, "+"; limit = 2)
        minusa = split(popfirst!(plusa), "-"; limit = 2)
        m = popfirst!(minusa)
        a = string("[", n > 1 ? n : "", "M", isempty(minusa) ? "" : string("-", first(minusa)), isempty(plusa) ? "" : string("+", first(plusa)), "]", c)
        ai = AdductIon(FormulaChemical(m), a)
        cl ? ChemicalLoss(ai) : ai
    end
end
Formula2Chemical(chemical::AbstractChemical, net_charge = 0; default_charge = net_charge, loss = 0) = chemical
Formula2Chemical(chemical::ChemicalPair, net_charge = 0; default_charge = net_charge, loss = 0) = ChemicalSeries(chemical)
function Formula2Chemical(formula::AbstractString, net_charge = 0; default_charge = net_charge, loss = 0) 
    # net_charge = abs(net_charge) < abs(default_charge) ? net_charge : default_charge
    net_charge = net_charge == 0 ? default_charge : net_charge
    if net_charge == 0
        loss = 0 
    else
        while (net_charge - loss) * net_charge <= 0 
            loss -= sign(net_charge)
        end
    end
    # check charge
    parse_formula(formula; charge = net_charge, loss)
end
PostFormula2Chemical(pre, chemical::AbstractChemical, net_charge = 0; default_charge = net_charge, loss = 0) = ChemicalPair(pre, Formula2Chemical(chemical, net_charge; default_charge, loss))
function PostFormula2Chemical(pre, chemical::ChemicalPair, net_charge = 0; default_charge = net_charge, loss = 0)
    pre = PostFormula2Chemical(pre, chemical.precursor, net_charge; default_charge, loss)
    PostFormula2Chemical(pre, chemical.product, detectedcharge(pre); default_charge, loss)
end
function PostFormula2Chemical(pre, formula::Pair, net_charge = 0; default_charge = net_charge, loss = 0) 
    pre = PostFormula2Chemical(pre, first(formula), net_charge; default_charge, loss)
    PostFormula2Chemical(pre, last(formula), detectedcharge(pre); default_charge, loss)
end
PostFormula2Chemical(pre, formula, net_charge = 0; default_charge = net_charge, loss = 0) = ChemicalPair(pre, Formula2Chemical(formula, net_charge; default_charge, loss))