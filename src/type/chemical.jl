"""
    AttributedChemical <: AbstractChemical

Abstract type for chemicals with additional attributes. The storage can be accessed by [`getchemicalattr`](@ref), and the property name can be accessed by [`getchemicalattrname`](@ref).

For Base functions, `copy`, `hash`, and `==`, generic methods are defined.
"""
abstract type AttributedChemical <: AbstractChemical end

"""
    ChemicalAttribute

Mutable object carrying attributes.

# Fields
* `attribute::Union{Nothing, Vector{Pair{Symbol, Any}}}`: the pairs represent names and values.

# Constructors
    ChemicalAttribute(attribute::ChemicalAttribute)
    ChemicalAttribute(attribute::AbstractVector{<:Pair{<:AbstractString}})
    ChemicalAttribute(attribute::AbstractVector{<:Tuple})
"""
mutable struct ChemicalAttribute
    attribute::Union{Nothing, Vector{Pair{Symbol, Any}}}
end

ChemicalAttribute(attribute::ChemicalAttribute) = attribute
ChemicalAttribute(attribute::AbstractVector{<:Pair{<:AbstractString}}) = 
    ChemicalAttribute(Pair{Symbol, Any}[Symbol(first(pair)) => last(pair) for pair in attribute])
ChemicalAttribute(attribute::AbstractVector{<:Tuple}) = 
    ChemicalAttribute(Pair{Symbol, Any}[Symbol(first(pair)) => last(pair) for pair in attribute])

"""
    buildattribute(::Nothing) -> Nothing
    buildattribute(attribute) -> ChemicalAttribute
    buildattribute(new_attr::Bool, kwargs) -> Union{Nothing, ChemicalAttribute}

Constructs a `ChemicalAttribute` object or returns `nothing` if `new_attr` is `false` and no keyword arguments are provided.
"""
buildattribute(::Nothing) = nothing
buildattribute(attribute) = ChemicalAttribute(attribute)
buildattribute(new_attr::Bool, kwargs) =
    isempty(kwargs) ? (new_attr ? ChemicalAttribute(nothing) : nothing) :
    ChemicalAttribute(Pair{Symbol, Any}[key => value for (key, value) in kwargs])

"""
    buildelements(elements) -> Vector{Pair{String, Int}}
    buildelements(formula::AbstractString) -> Vector{Pair{String, Int}}

Constructs a vector of element-count pairs.
"""
buildelements(elements::Vector{Pair{String, Int}}) = elements
buildelements(elements) = collect(Pair{String, Int}, elements)
buildelements(formula::AbstractString) = chemicalelements(formula)

const ChemicalSchemeNull = Union{Nothing, AbstractChemicalScheme}

"""
    SchemeVector

# Fields
* `vector::Vector{Pair{ChemicalSchemeNull, Vector{Pair{AbstractChemicalScheme, AbstractChemicalScheme}}}}`

The constructor will try to convert any input into the correct vector type. Non vector input is vectorized.     
Then the vector is converted to vector of pairs by fetching the first and the last element of each element in the vector. 
Each of the last element is additionally converted to vector of pairs by the same method.

For two-argument method, the aruments are regarded as a pair. For three-argument method, the argument becomes `[arg1 => [arg2 => arg3]]`.
"""
struct SchemeVector
    vector::Vector{Pair{ChemicalSchemeNull, Vector{Pair{AbstractChemicalScheme, AbstractChemicalScheme}}}}
    function SchemeVector(vec::Vector{Pair{ChemicalSchemeNull, Vector{Pair{AbstractChemicalScheme, AbstractChemicalScheme}}}}) 
        new(_schemesort!(vec))
    end
    function SchemeVector(vec::AbstractVector) 
        new(_schemesort!([first(x) => _scheme_vec_pair(last(x)) for x in vec]))
    end
end
SchemeVector(x::SchemeVector) = x
SchemeVector() = SchemeVector(Vector{Pair{ChemicalSchemeNull, Vector{Pair{AbstractChemicalScheme, AbstractChemicalScheme}}}}[])
SchemeVector(x) = SchemeVector(vectorize(x))
function _schemesort!(vec)
    sort!(vec; by = first)
    for x in vec 
        sort!(last(x); by = first)
    end
    vec
end
_scheme_vec_pair(x) = _scheme_vec_pair(vectorize(x))
_scheme_vec_pair(x::AbstractVector{<:Pair}) = x
_scheme_vec_pair(x::AbstractVector) = [first(y) => last(y) for y in x]
_scheme_vec_pair(x, y) = [x => y]
_scheme_pair(x::Pair) = x
_scheme_pair(x) = first(x) => last(x)
_scheme_pair(x, y) = x => y

SchemeVector(pre, postpair) = SchemeVector([pre => _scheme_vec_pair(postpair)])
SchemeVector(pre, postold, postnew) = SchemeVector([pre => [postold => postnew]])

"""
    SchemeStructure

# Fields
* `transformer::Union{Nothing, SchemeVector}`
* `merger::Union{Nothing, SchemeVector}`

# Constructors
    SchemeStructure(; transformer = nothing, merger = nothing)
"""
mutable struct SchemeStructure
    transformer::Union{Nothing, SchemeVector}
    merger::Union{Nothing, SchemeVector}
end

SchemeStructure(; transformer = nothing, merger = nothing) = SchemeStructure(transformer, merger)

"""
    Chemical <: AttributedChemical

Unstructured chemical type with its name, elements (formula), and additional attributes.

# Fields 
* `name::String`: a unique chemical name.
* `elements::Vector{Pair{String, Int}}`: chemical elements.
* `attribute::Union{Nothing, ChemicalAttribute}`: additional attributes.

# Constructors
    Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}, attribute)
    Chemical(name::AbstractString, elements::Dict{String, Int}, attribute)
    Chemical(name::AbstractString, formula::AbstractString, attribute)
    Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}; new_attr = true, kwargs...)
    Chemical(name::AbstractString, elements::Dict{String, Int}; new_attr = true, kwargs...)
    Chemical(name::AbstractString, formula::AbstractString; new_attr = true, kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.
"""
struct Chemical <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    attribute::Union{Nothing, ChemicalAttribute}
end

Chemical(name::AbstractString, elements, attribute) =
    Chemical(String(name), buildelements(elements), buildattribute(attribute))
Chemical(name::AbstractString, elements; new_attr = true, kwargs...) =
    Chemical(String(name), buildelements(elements), buildattribute(new_attr, kwargs))

"""
    FormulaChemical <: AttributedChemical

Unstructured chemical type with elements (formula), and additional attributes. Chemical name will be formula.

# Fields 
* `elements::Vector{Pair{String, Int}}`: chemical elements.
* `attribute::Union{Nothing, ChemicalAttribute}`: additional attributes.

# Constructors
    FormulaChemical(elements::Vector{Pair{String, Int}}, attribute)
    FormulaChemical(elements::Dict{String, Int}, attribute)
    FormulaChemical(formula::AbstractString, attribute)
    FormulaChemical(elements::Vector{Pair{String, Int}}; new_attr = true, kwargs...)
    FormulaChemical(elements::Dict{String, Int}; new_attr = true, kwargs...)
    FormulaChemical(formula::AbstractString; new_attr = true, kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.
"""
struct FormulaChemical <: AttributedChemical
    elements::Vector{Pair{String, Int}}
    attribute::Union{Nothing, ChemicalAttribute}
    FormulaChemical(elements, attribute) = new(buildelements(elements), buildattribute(attribute))
end

FormulaChemical(elements; new_attr = true, kwargs...) =
    FormulaChemical(buildelements(elements), buildattribute(new_attr, kwargs))

abstract type AbstractChemicalWrapper <: AbstractChemical end
"""
    ChemicalTransition{T<:AbstractChemicalScheme} <: AbstractChemicalWrapper

Chemical transition in MSⁿ. Products can be any subtype of `AbstractChemicalScheme` representing chemical entity or species.

# Fields 
* `transition::Vector{T}`.

# Constructors
    ChemicalTransition(transition::Vector)
    ChemicalTransition(precursor, products...)
    
`products` are pushed into `precursor` to construct `transition`.
"""
struct ChemicalTransition{T<:AbstractChemicalScheme} <: AbstractChemicalWrapper
    transition::Vector{T}
end

function ChemicalTransition(ct...) 
    ChemicalTransition(mapreduce(_transition, vcat, ct))
end

_transition(x::ChemicalTransition) = chemicaltransition(x)
_transition(x) = x

"""
    Isobars{T<:AbstractChemical, N} <: AbstractChemicalWrapper

Chemicals with similar m/z.

# Fields 
* `chemicals::Vector{T}`: a vector of chemicals.
* `abundnace::VecOrMat{N}`: the abundance of each chemical. If chemicals are trasitions, this should be a matrix, and each column is the abundance of each ms stage.

# Constructors 
    Isobars(chemicals::Vector, abundance::Vector)
    Isobars(chemicals::Vector{<:ChemicalTransition}, abundance::Vector)
    Isobars(chemicals::Vector{<:ChemicalTransition}, abundance::Matrix)
"""
struct Isobars{T<:AbstractChemical, N} <: AbstractChemicalWrapper
    chemicals::Vector{T}
    abundance::VecOrMat{N}
    function Isobars(chemicals::Vector{T}, abundance::Vector{N}) where {T, N}
        id = sortperm(abundance; rev = true)
        new{T, N}(chemicals[id], abundance[id])
    end
    function Isobars(chemicals::Vector{T}, abundance::Vector{N}) where {T<:ChemicalTransition, N}
        allequal(msstage, chemicals) || throw(ArgumentError("All chemicals should have the same `msstage`."))
        length(chemicals) == length(abundance) || throw(ArgumentError("The length of `chemicals` and `abundance` are not equaled."))
        id = sortperm(abundance; rev = true)
        ab = hcat([abundance[id] for _ in 1:msstage(chemicals[begin])]...)
        new{T, N}(chemicals[id], ab)
    end
    function Isobars(chemicals::Vector{T}, abundance::Vector{Vector{N}}) where {T<:ChemicalTransition, N}
        allequal(msstage, chemicals) || throw(ArgumentError("All chemicals should have the same `msstage`."))
        length(abundance) == msstage(chemicals[begin]) || throw(ArgumentError("The length of `abundance` should equal to `msstage`."))
        for ab in abundance
            length(chemicals) == length(ab) || throw(ArgumentError("The length of chemicals and each abundance are not equaled."))
        end
        id = sortperm(abundance[end]; rev = true)
        ab = hcat([ab[id] for ab in abundance]...)
        new{T, N}(chemicals[id], ab)
    end
    function Isobars(chemicals::Vector{T}, abundance::Matrix{N}) where {T<:ChemicalTransition, N}
        allequal(msstage, chemicals) || throw(ArgumentError("All chemicals should have the same `msstage`."))
        size(abundance, 1) == length(chemicals) || throw(ArgumentError("The column size of `abundance` should equal to the length of `chemicals`."))
        size(abundance, 2) == msstage(chemicals[begin]) || throw(ArgumentError("The row size of `abundance` should equal to `msstage`."))
        id = sortperm(abundance[:, end]; rev = true)
        new{T, N}(chemicals[id], abundance[id, :])
    end
end

Isobars(chemicals::AbstractVector, abundance::AbstractArray) = _Isobars(collect(chemicals), abundance)
_Isobars(chemicals::AbstractVector, abundance::AbstractArray) = Isobars(chemicals, collect(abundance))
_Isobars(chemicals::AbstractVector, abundance::VecOrMat) = Isobars(chemicals, abundance)

struct ElementsVector
    elements::Vector{String}
    numbers::Vector{Int}
end

ElementsVector() = ElementsVector(String[], Int[])

"""
    Isotopomers{T<:AbstractChemical} <: AbstractChemicalWrapper

Chemicals differed from isotopic replacement location.

# Fields 
* `parent::T`: shared chemical structure of isotopomers prior to isotopic replacement. 
* `isotopes::ElementsVector`: isotopes-number pairs of isotopic replacement.

# Constructors
    Isotopomers(parent::AbstractChemical, isotopes::ElementsVector)
    Isotopomers(parent::AbstractChemical, fullformula::String)
    Isotopomers(parent::AbstractChemical, fullelements::Dict)
    Isotopomers(parent::AbstractChemical, fullelements::Vector{Pair{String, Int}})

All minor isotopes are regarded as isotopic replacement in `fullformula` and `fullelements`.
"""
struct Isotopomers{T<:AbstractChemical} <: AbstractChemicalWrapper
    parent::T 
    isotopes::ElementsVector
end

function Isotopomers(chemical::AbstractChemical, fullformula::String)
    Isotopomers(chemicalparent(chemical), dictionary_elements(chemicalelements(fullformula)))
end

function Isotopomers(chemical::AbstractChemical, fullelements::Vector{Pair{String, Int}})
    Isotopomers(chemicalparent(chemical), dictionary_elements(fullelements))
end

function Isotopomers(chemical::AbstractChemical, fullelements::Dict)
    parent = chemicalparent(chemical)
    dp = dictionary_elements(chemicalelements(parent))
    ev = ElementsVector(collect(keys(fullelements)), collect(values(fullelements)))
    del = Int[]
    for (i, (k, n)) in enumerate(ev)
        iselement(k) && (push!(del, i); continue)
        ev.numbers[i] = n - get(dp, k, 0) 
    end
    deleteat!(ev.elements, del)
    deleteat!(ev.numbers, del)
    Isotopomers(parent, ev)
end

"""
    Groupedisotopomers{T<:AbstractChemical, N} <: AbstractChemicalWrapper

Isotopomerized chemicals grouped by mass-shift index.

# Fields 
* `parent::T`: shared chemical structure prior to isotopic replacement. 
* `index::Int`: mass-shift index.
* `isotope::String`: reference isotope for [`mass_shift_index`](@ref).
* `isotopes::Vector{ElementsVector}`: Isotopes-number pairs of isotopic replacements of each isotopomers.
* `abundance::Vector{N}`: abundance of each isotopomers.
"""
struct Groupedisotopomers{T<:AbstractChemical, N} <: AbstractChemicalWrapper
    parent::T 
    index::Int
    isotope::String
    isotopes::Vector{ElementsVector}
    abundance::Vector{N}
    function Groupedisotopomers(parent::T, index::Int, isotope::String, isotopes::Vector{ElementsVector}, abundance::Vector{N}) where {T, N}
        id = sortperm(abundance)
        new{T, N}(parent, index, isotope, isotopes[id], abundance[id])
    end
end

"""
    ChemicalSeries(chemical::AbstractChemicalScheme)
    ChemicalSeries(pair::Pair)
    ChemicalSeries(chemicals::AbstractVector)

Transform chemical into valid chemical structure. Multiple chemicals are converted into `ChemicalTransition`.
"""
ChemicalSeries(cc::AbstractChemicalScheme) = cc
ChemicalSeries(cc::ChemicalTransition) = cc
ChemicalSeries(ct...) = ChemicalTransition(ct...) 
ChemicalSeries(v::AbstractVector) = length(v) < 2 ? ChemicalSeries(first(v)) : ChemicalTransition(v...)
ChemicalSeries(v::Pair) = ChemicalTransition(_ChemicalSeries(v)...)
_ChemicalSeries(v::Pair) = (_ChemicalSeries(first(v))..., _ChemicalSeries(last(v))...)
_ChemicalSeries(v::ChemicalTransition) = (chemicaltransition(v)..., )
_ChemicalSeries(v::AbstractChemicalScheme) = (v, ) 

"""
    Electron <: AttributedChemical

Electron.

# Fields 
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Electron(; new_attr = true, name = "Electron", elements = Pair{String, Int}[], abbreviation = "e", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `-1`.
"""
struct Electron <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Proton <: AttributedChemical

Proton.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Proton(; new_attr = true, name = "Proton", elements = ["H" => 1], abbreviation = "H", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `1`.
"""
struct Proton <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Water <: AttributedChemical

Water.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Water(; new_attr = true, name = "Water", elements = ["H" => 2, "O" => 1], abbreviation = "H2O", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `0`.
"""
struct Water <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Ammonia <: AttributedChemical

Ammonia.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Ammonia(; new_attr = true, name = "Ammonia", elements = ["N" => 1, "H" => 3], abbreviation = "NH3", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `0`.
"""
struct Ammonia <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Ammonium <: AttributedChemical

Ammonium.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Ammonium(; new_attr = true, name = "Ammonium", elements = ["N" => 1, "H" => 4], abbreviation = "NH4", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `1`.
"""
struct Ammonium <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Sodium <: AttributedChemical

Sodium.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Sodium(; new_attr = true, name = "Sodium", elements = ["Na" => 1], abbreviation = "Na", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `1`.
"""
struct Sodium <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Potassium <: AttributedChemical

Potassium.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Potassium(; new_attr = true, name = "Potassium", elements = ["K" => 1], abbreviation = "K", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `1`.
"""
struct Potassium <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Lithium <: AttributedChemical

Lithium.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Lithium(; new_attr = true, name = "Lithium", elements = ["Li" => 1], abbreviation = "Li", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `1`.
"""
struct Lithium <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Silver <: AttributedChemical

Silver.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Silver(; new_attr = true, name = "Silver", elements = ["Ag" => 1], abbreviation = "Ag", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `1`.
"""
struct Silver <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Acetate <: AttributedChemical

Acetate.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Acetate(; new_attr = true, name = "Acetate", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1], abbreviation = "OAc", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `-1`.
"""
struct Acetate <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Formate <: AttributedChemical

Formate.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Formate(; new_attr = true, name = "Formate", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1], abbreviation = "OFo", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `-1`.
"""
struct Formate <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    AceticAcid <: AttributedChemical

Acetic acid.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    AceticAcid(; new_attr = true, name = "Acetic Acid", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1], abbreviation = "HOAc", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `0`.
"""
struct AceticAcid <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    FormicAcid <: AttributedChemical

Formic acid.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    FormicAcid(; new_attr = true, name = "Formic Acid", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1], abbreviation = "HOFo", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `0`.
"""
struct FormicAcid <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    MethylAcetate <: AttributedChemical

Methyl acetate.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    MethylAcetate(; new_attr = true, name = "Methyl Acetate", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3], abbreviation = "MeOAc", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `0`.
"""
struct MethylAcetate <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    MethylFormate <: AttributedChemical

Methyl formate.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    MethylFormate(; new_attr = true, name = "Methyl Formate", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3], abbreviation = "MeOFo", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `0`.
"""
struct MethylFormate <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Fluoride <: AttributedChemical

Fluoride.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Fluoride(; new_attr = true, name = "Fluoride", elements = ["F" => 1], abbreviation = "F", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `-1`.
"""
struct Fluoride <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Chloride <: AttributedChemical

Chloride.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Chloride(; new_attr = true, name = "Chloride", elements = ["Cl" => 1], abbreviation = "Cl", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `-1`.
"""
struct Chloride <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 
"""
    Methenium <: AttributedChemical

Methenium.

# Fields
* `name::String`
* `elements::Vector{Pair{String, Int}}`
* `abbreviation::String`
* `attribute::Union{Nothing, ChemicalAttribute}`

# Constructors
    Methenium(; new_attr = true, name = "Methenium", elements = ["C" => 1, "H" => 3], abbreviation = "Me", kwargs...)

## Arguments
* `new_attr::Bool`: whether allowing adding new attributes afterwards. If `kwargs` is not empty, it is automatically set `true`.
* `kwargs` are collected into field `attribute`.

# Additional Attributes
* `charge::Int`: `1`.
"""
struct Methenium <: AttributedChemical
    name::String
    elements::Vector{Pair{String, Int}}
    abbreviation::String
    attribute::Union{Nothing, ChemicalAttribute}
end 

"""
    const StandardChemical::Dict{Type, AbstractChemical}

Standard or default chemicals of chemical types. The types are keys, and chemical instances are values. This is mainly used for schemes.

Default keys:
* `Electron`
* `Proton`
* `Water`
* `Ammonia`
* `Ammonium`
* `Sodium`
* `Potassium`
* `Lithium`
* `Silver`
* `Acetate`
* `Formate`
* `AceticAcid`
* `FormicAcid`
* `MethylAcetate`
* `MethylFormate`
* `Fluoride`
* `Chloride`
* `Methenium`
"""
const StandardChemical = Dict{Type, AbstractChemical}()

"""
    defaultchemicalkwargs(::Type{T}) -> NamedTuple

Default `name`, `elements` and `abbreviation` of standard chemical type `T`.
"""
defaultchemicalkwargs(::Type{<:AbstractChemical}) = (; )
defaultchemicalkwargs(::Type{Electron}) = (; name = "Electron", elements = Pair{String, Int}[], abbreviation = "e")
defaultchemicalkwargs(::Type{Proton}) = (; name = "Proton", elements = ["H" => 1], abbreviation = "H")
defaultchemicalkwargs(::Type{Water}) = (; name = "Water", elements = ["H" => 2, "O" => 1], abbreviation = "H2O")
defaultchemicalkwargs(::Type{Ammonia}) = (; name = "Ammonia", elements = ["N" => 1, "H" => 3], abbreviation = "NH3")
defaultchemicalkwargs(::Type{Ammonium}) = (; name = "Ammonium", elements = ["N" => 1, "H" => 4], abbreviation = "NH4")
defaultchemicalkwargs(::Type{Sodium}) = (; name = "Sodium", elements = ["Na" => 1], abbreviation = "Na")
defaultchemicalkwargs(::Type{Potassium}) = (; name = "Potassium", elements = ["K" => 1], abbreviation = "K")
defaultchemicalkwargs(::Type{Lithium}) = (; name = "Lithium", elements = ["Li" => 1], abbreviation = "Li")
defaultchemicalkwargs(::Type{Silver}) = (; name = "Silver", elements = ["Ag" => 1], abbreviation = "Ag")
defaultchemicalkwargs(::Type{Acetate}) = (; name = "Acetate", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1], abbreviation = "OAc")
defaultchemicalkwargs(::Type{Formate}) = (; name = "Formate", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1], abbreviation = "OFo")
defaultchemicalkwargs(::Type{AceticAcid}) = (; name = "Acetic Acid", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1], abbreviation = "HOAc")
defaultchemicalkwargs(::Type{FormicAcid}) = (; name = "Formic Acid", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1], abbreviation = "HOFo")
defaultchemicalkwargs(::Type{MethylAcetate}) = (; name = "Methyl Acetate", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3], abbreviation = "MeOAc")
defaultchemicalkwargs(::Type{MethylFormate}) = (; name = "Methyl Formate", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3], abbreviation = "MeOFo")
defaultchemicalkwargs(::Type{Fluoride}) = (; name = "Fluoride", elements = ["F" => 1], abbreviation = "F")
defaultchemicalkwargs(::Type{Chloride}) = (; name = "Chloride", elements = ["Cl" => 1], abbreviation = "Cl")
defaultchemicalkwargs(::Type{Methenium}) = (; name = "Methenium", elements = ["C" => 1, "H" => 3], abbreviation = "Me")

# User-provided kwargs override defaults; name, elements and abbreviation are fields, the rest become attributes
"""
    make_attributedchemical(::Type{T}; new_attr = true, kwargs...) -> AttributedChemical

Create an attributed chemical of type `T` with default and user-provided keyword arguments. Default values are taken from [`defaultchemicalkwargs`](@ref).

# Arguments
* `::Type{T}`: the type of the attributed chemical to create.
* `new_attr`: whether to build a new attribute.
* `kwargs`: additional keyword arguments to override defaults.
"""
function make_attributedchemical(::Type{T}; new_attr = true, kwargs...) where {T<:AttributedChemical}
    kw = (; pairs(defaultchemicalkwargs(T))..., kwargs...)
    attr = Base.structdiff(kw, NamedTuple{(:name, :elements, :abbreviation)})
    T(String(kw.name), collect(Pair{String, Int}, kw.elements), String(kw.abbreviation), buildattribute(new_attr, pairs(attr)))
end

Electron(; kwargs...) = make_attributedchemical(Electron; kwargs...)
Proton(; kwargs...) = make_attributedchemical(Proton; kwargs...)
Water(; kwargs...) = make_attributedchemical(Water; kwargs...)
Ammonia(; kwargs...) = make_attributedchemical(Ammonia; kwargs...)
Ammonium(; kwargs...) = make_attributedchemical(Ammonium; kwargs...)
Sodium(; kwargs...) = make_attributedchemical(Sodium; kwargs...)
Potassium(; kwargs...) = make_attributedchemical(Potassium; kwargs...)
Lithium(; kwargs...) = make_attributedchemical(Lithium; kwargs...)
Silver(; kwargs...) = make_attributedchemical(Silver; kwargs...)
Acetate(; kwargs...) = make_attributedchemical(Acetate; kwargs...)
Formate(; kwargs...) = make_attributedchemical(Formate; kwargs...)
AceticAcid(; kwargs...) = make_attributedchemical(AceticAcid; kwargs...)
FormicAcid(; kwargs...) = make_attributedchemical(FormicAcid; kwargs...)
MethylAcetate(; kwargs...) = make_attributedchemical(MethylAcetate; kwargs...)
MethylFormate(; kwargs...) = make_attributedchemical(MethylFormate; kwargs...)
Fluoride(; kwargs...) = make_attributedchemical(Fluoride; kwargs...)
Chloride(; kwargs...) = make_attributedchemical(Chloride; kwargs...)
Methenium(; kwargs...) = make_attributedchemical(Methenium; kwargs...)

push!(StandardChemical, Electron => Electron(; new_attr = false))
push!(StandardChemical, Proton => Proton(; new_attr = false))
push!(StandardChemical, Water => Water(; new_attr = false))
push!(StandardChemical, Ammonia => Ammonia(; new_attr = false))
push!(StandardChemical, Ammonium => Ammonium(; new_attr = false))
push!(StandardChemical, Sodium => Sodium(; new_attr = false))
push!(StandardChemical, Potassium => Potassium(; new_attr = false))
push!(StandardChemical, Lithium => Lithium(; new_attr = false))
push!(StandardChemical, Silver => Silver(; new_attr = false))
push!(StandardChemical, Acetate => Acetate(; new_attr = false))
push!(StandardChemical, Formate => Formate(; new_attr = false))
push!(StandardChemical, AceticAcid => AceticAcid(; new_attr = false))
push!(StandardChemical, FormicAcid => FormicAcid(; new_attr = false))
push!(StandardChemical, MethylAcetate => MethylAcetate(; new_attr = false))
push!(StandardChemical, MethylFormate => MethylFormate(; new_attr = false))
push!(StandardChemical, Fluoride => Fluoride(; new_attr = false))
push!(StandardChemical, Chloride => Chloride(; new_attr = false))
push!(StandardChemical, Methenium => Methenium(; new_attr = false))

"""
    const GenericChemical = Union{Chemical, FormulaChemical}

Generic chemical types.
"""
const GenericChemical = Union{Chemical, FormulaChemical}

"""
    AbstractAdductIon{S, T} <: AbstractChemical

Abstract type for adduct ions with core chemical type `S` and adduct type `T`.

# Special attributes
* [`ncore`](@ref) `-> Int`: number of core chemical "M" in adduct ion representation "[M+X]n+". 
* [`ioncore`](@ref) `-> S`: the core chemical undergoing ionization. 
* [`ionadduct`](@ref) `-> T`: the adduct formed during ionization. 
"""
abstract type AbstractAdductIon{S, T} <: AbstractChemical end

"""
    AdductIon{S<:AbstractChemical, T<:AbstractScheme} <: AbstractAdductIon{S, T}

Adduct ions formed in mass spectrometry.

# Fields
* `core`: the core chemical undergoing ionization. 
* `adduct`: the adduct formed during ionization.
* `ncore`: number of core chemical undergoing ionization. 

# Constructors
    AdductIon(core::AbstractChemical, adduct::AbstractScheme, ncore = 1)
    AdductIon(core::AbstractChemical, adduct_string::AbstractString)

`adduct_string` is parsed by `parse_adduct(adduct_string; args = true)`.
"""
struct AdductIon{S<:AbstractChemical, T<:AbstractScheme} <: AbstractAdductIon{S, T}
    core::S
    adduct::T
    ncore::Int
    function AdductIon(core::S, adduct::AbstractScheme, ncore::Int) where {S<:AbstractChemical}
        adduct = completescheme(core, adduct)
        new{S, typeof(adduct)}(core, adduct, ncore)
    end
end

AdductIon(cc::AbstractChemical, a::AbstractScheme) = AdductIon(cc, a, 1)
AdductIon(cc::AbstractChemical, a::AbstractString, n = 1) = AdductIon(cc, parse_adduct(a; args = true)...)