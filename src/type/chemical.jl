"""
    Chemical <: AbstractChemical

Unstructured chemical type with its name, elements (formula), and additional properties.

# Fields 
* `name::String`: a unique chemical name.
* `elements::Vector{Pair{String, Int}}`: chemical elements.
* `property::Vector{Pair{Symbol, Any}}`: additional properties; the pairs repressent names and values.

# Constructors
* `Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}, property::Vector{Pair{Symbol, Any}})`
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
* `FormulaChemical(elements::Vector{Pair{String, Int}}, property::Vector{Pair{Symbol, Any}})`
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
    ChemicalTransition{T} <: AbstractChemical

Chemical transition in MSⁿ. Products can be a `ChemicalLoss` or `ChemicalGain`.

# Fields 
* `transition::Vector{T}`

# Constructors
* `ChemicalTransition(transition::Vector)`
* `ChemicalTransition(precursor, products...)`: push `products` into `precursor`.
"""
struct ChemicalTransition{T} <: AbstractChemical
    transition::Vector{T}
end

_ChemicalTransition(ct::ChemicalTransition, product::AbstractChemical) = push!(copy(ct.transition), product)
_ChemicalTransition(ct::ChemicalTransition, product::ChemicalTransition) = vcat(ct.transition, product.transition)
_ChemicalTransition(cc::AbstractChemical, product::AbstractChemical) = [cc, product]
_ChemicalTransition(cc::AbstractChemical, product::ChemicalTransition) = vcat(cc, product.transition)

function ChemicalTransition(ct::AbstractChemical, product...) 
    for p in product
        ct = _ChemicalTransition(ct, p)
    end
    ChemicalTransition(ct)
end

"""
    Isobars{T <: AbstractChemical, N} <: AbstractChemical

Chemicals with similar m/z.

# Fields 
* `chemicals::Vector{T}`: a vector of chemicals.
* `abundnace::VecOrMat{Float64}`: the abundance of each chemical. If chemicals are trasitions, this should be a matrix, and each column is the abundance of each ms stage.

# Constructors 
* `Isobars(chemicals::Vector, abundance::Vector)`
* `Isobars(chemicals::Vector{<: ChemicalTransition}, abundance::Vector)`
* `Isobars(chemicals::Vector{<: ChemicalTransition}, abundance::Matrix)`
"""
struct Isobars{T <: AbstractChemical, N} <: AbstractChemical
    chemicals::Vector{T}
    abundance::VecOrMat{N}
    function Isobars(chemicals::Vector{T}, abundance::Vector{N}) where {T, N}
        id = sortperm(abundance; rev = true)
        new{T, N}(chemicals[id], abundance[id])
    end
    function Isobars(chemicals::Vector{T}, abundance::Vector{N}) where {T <: ChemicalTransition, N}
        allequal(msstage, chemicals) || throw(ArgumentError("All chemicals should have the same `msstage`."))
        id = sortperm(abundance; rev = true)
        ab = hcat([abundance[id] for _ in eachindex(chemicals)]...)
        new{T, N}(chemicals[id], ab)
    end
    function Isobars(chemicals::Vector{T}, abundance::Matrix{N}) where {T <: ChemicalTransition, N}
        allequal(msstage, chemicals) || throw(ArgumentError("All chemicals should have the same `msstage`."))
        id = sortperm(abundance[:, end]; rev = true)
        new{T, N}(chemicals[id], abundance[id, :])
    end
end

Isobars(chemicals::AbstractVector, abundance::AbstractArray) = _Isobars(collect(chemicals), abundance)
_Isobars(chemicals::AbstractVector, abundance::AbstractArray) = Isobars(chemicals, collect(abundance))
_Isobars(chemicals::AbstractVector, abundance::VecOrMat) = Isobars(chemicals, abundance)

"""
    Isotopomers{T <: AbstractChemical} <: AbstractChemical

Chemicals differed from isotopic replacement location.

# Fields 
* `parent::T`: shared chemical structure of isotopomers prior to isotopic replacement. 
* `isotopes::Vector{Pair{String, Int}}`: Isotopes-number pairs of isotopic replacement.

# Constructors
* `Isotopomers(parent::AbstractChemical, isotopes::Vector{Pair{String, Int}})`
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
    Isotopomers(parent, dictionary_elements(chemicalelements(fullformula)))
end

function Isotopomers(parent::AbstractChemical, fullelements::Dictionary)
    dp = dictionary_elements(chemicalelements(parent))
    dr = deepcopy(fullelements)
    for k in keys(fullelements)
        haskey(elements_isotopes(), k) && (delete!(dr, k); continue)
        dr[k] -= get(dp, k, 0) 
    end
    Isotopomers(parent, [k => v for (k, v) in pairs(dr)])
end

"""
    Groupedisotopomers{T <: AbstractChemical, N} <: AbstractChemical

Isotopomers grouped by isotopomer state.

# Fields 
* `parent::T`: shared chemical structure of isotopomers prior to isotopic replacement. 
* `state::Int`: isotopomer state.
* `isotope::String`: isotope for computing isotopomer state.
* `isotopes::Vector{Vector{Pair{String, Int}}}`: Isotopes-number pairs of isotopic replacements of each isotopomers.
* `abundance::Vector{N}`: abundance of each isotopomers.
"""
struct Groupedisotopomers{T <: AbstractChemical, N} <: AbstractChemical
    parent::T 
    state::Int
    isotope::String
    isotopes::Vector{Vector{Pair{String, Int}}}
    abundance::Vector{N}
    function Groupedisotopomers(parent::T, state::Int, isotope::String, isotopes::Vector{Vector{Pair{String, Int}}}, abundance::Vector{N}) where {T, N}
        id = sortperm(abundance)
        new{T, N}(parent, state, isotope, isotopes[id], abundance[id])
    end
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
    ChemicalGain{T <: AbstractChemical} <: AbstractChemical

Chemical gain to a precursor. This product is not detected in MS; the merged chemical is detected instead.

# Fields 
* `chemical::T`
"""
struct ChemicalGain{T <: AbstractChemical} <: AbstractChemical
    chemical::T 
end

"""
    ChemicalSeries(chemical::AbstractChemical)
    ChemicalSeries(formula::AbstractString; charge = 0, gain = 0, loss = 0)
    ChemicalSeries(pair::Pair; charge = 0, gain = 0, loss = 0)

Transform input into valid chemical structure. 

* `formula` is parsed into `FormulaChemical`. 
If it has adduct ion like structure, i.e. `"[formula+adduct]x+"` or `"[formula+adduct]x-"`, `FormulaChemical` will be wrapped by `AdductIon` with adduct (optional) if the charge is nonzero.
It is further wrapped as chemical gain or loss if starting with + or -. 
* `pair` can be a series of precursor-product pair. Individual chemical is splated before parsing. 
"""
ChemicalSeries(cc::AbstractChemical) = cc
ChemicalSeries(cc::ChemicalTransition) = ChemicalTransition(_ChemicalSeries(cc))
_ChemicalSeries(cc::AbstractChemical, v = AbstractChemical[]) = push!(v, cc)
function _ChemicalSeries(cc::ChemicalTransition, v = AbstractChemical[]) 
    for c in chemicaltransition(cc)
        _ChemicalSeries(c, v)
    end
    v 
end
ChemicalSeries(ct, product...) = ChemicalTransition(ct, product...) 
ChemicalSeries(v::AbstractVector; charge = 0, gain = 0, loss = 0) = length(v) < 2 ? ChemicalSeries(first(v)) : ChemicalTransition(v)
ChemicalSeries(formula; charge = 0, gain = 0, loss = 0) = Formula2Chemical(formula, charge; gain, loss)

function Formula2Chemical(formula::Pair, net_charge = 0; default_charge = net_charge, gain = 0, loss = 0) 
    pre = Formula2Chemical(first(formula), net_charge; default_charge, gain, loss)
    ChemicalSeries(PostFormula2Chemical(pre, last(formula), detectedcharge(pre); default_charge, gain, loss))
end

Formula2Chemical(chemical::AbstractChemical, net_charge = 0; default_charge = net_charge, gain = 0, loss = 0) = chemical
Formula2Chemical(chemical::ChemicalTransition, net_charge = 0; default_charge = net_charge, gain = 0, loss = 0) = chemical.transition
function Formula2Chemical(formula::AbstractString, net_charge = 0; default_charge = net_charge, gain = 0, loss = 0) 
    if startswith(formula, "-") 
        if net_charge == 0 
            loss = 0
        else
            while (net_charge - loss) * net_charge <= 0 
                loss -= sign(net_charge)
            end
        end
        ChemicalLoss(parse_chemical(FormulaChemical, formula[begin + 1:end]; charge = loss))
    elseif startswith(formula, "+")
        if net_charge == 0 
            gain = 0
        else
            while (net_charge + gain) * net_charge <= 0 
                gain += sign(net_charge)
            end
        end
        ChemicalGain(parse_chemical(FormulaChemical, formula[begin + 1:end]; charge = gain))
    else
        charge = min(abs(net_charge), abs(default_charge)) * sign(net_charge)
        parse_chemical(FormulaChemical, formula; charge)
    end
end
PostFormula2Chemical(pre, chemical::AbstractChemical, net_charge = 0; default_charge = net_charge, gain = 0, loss = 0) = vcat(vectorize(pre), Formula2Chemical(chemical, net_charge; default_charge, gain, loss))

function PostFormula2Chemical(pre, formula::Pair, net_charge = 0; default_charge = net_charge, gain = 0, loss = 0) 
    pre = PostFormula2Chemical(pre, first(formula), net_charge; default_charge, gain, loss)
    PostFormula2Chemical(pre, last(formula), detectedcharge(pre); default_charge, gain, loss)
end
PostFormula2Chemical(pre, formula, net_charge = 0; default_charge = net_charge, gain = 0, loss = 0) = vcat(vectorize(pre), Formula2Chemical(formula, net_charge; default_charge, gain, loss))