abstract type AbstractChemicalParser end
abstract type AbstractAdductParser end

"""
    AdductParser{T} <: AbstractAdductParser
    AdductParser(T::Bool = false, parser::AbstractChemicalParser = ChemicalSchemeParser())

Default adduct parser. `T` determines whether to return a tuple (`true`) or a named tuple (`false`). A named tuple is useful for user input from a data table and for ionization via `ionize`; a tuple is useful for the fixed adduct-ion constructor interface.

# Fields 
* `chemicalparser::T`: chemical parser for each chemical scheme. The default parser is [`ChemicalSchemeParser`](@ref).
"""
struct AdductParser{T} <: AbstractAdductParser
    chemicalparser::AbstractChemicalParser
end
AdductParser() = AdductParser{false}(ChemicalSchemeParser())
AdductParser(x::Bool) = AdductParser{x}(ChemicalSchemeParser())
AdductParser(parser::AbstractChemicalParser) = AdductParser{false}(parser)
AdductParser(x::Bool, parser::AbstractChemicalParser) = AdductParser{x}(parser)

"""
    FormulaChemicalParser <: AbstractChemicalParser
    FormulaChemicalParser(; kwargs...)

Default chemical parser that parses an input string (chemical formula) into [`FormulaChemical`](@ref).

# Fields 
* `attribute::Vector{Pair{Symbol, Any}}`: additional attributes. These can be mutated, added, or deleted with keyword arguments of [`parse_chemical`](@ref). Values of `nothing` are deleted.

In function [`parse_chemical`](@ref), keyword argument `new_attr` determines whether allowing adding new attributes afterwards. The default is `true`.
If the final attribute is not empty, it is automatically set `true`.

`kwargs` are collected into field `attribute`. 
"""
struct FormulaChemicalParser <: AbstractChemicalParser
    attribute::Vector{Pair{Symbol, Any}}
end
FormulaChemicalParser(; kwargs...) = FormulaChemicalParser(collect(kwargs))

"""
    ChemicalParser <: AbstractChemicalParser
    ChemicalParser(; kwargs...)

Chemical parser that parses an input string (chemical name) into [`Chemical`](@ref).

# Fields 
* `attribute::Vector{Pair{Symbol, Any}}`: additional attributes. It must include `:formula` or `:elements`. These can be mutated, added, or deleted with keyword arguments of [`parse_chemical`](@ref). Values of `nothing` are deleted.

In function [`parse_chemical`](@ref), keyword argument `new_attr` determines whether allowing adding new attributes afterwards. The default is `true`.
If the final attribute is not empty, it is automatically set `true`.

`kwargs` are collected into field `attribute`.
"""
struct ChemicalParser <: AbstractChemicalParser
    attribute::Vector{Pair{Symbol, Any}}
end
ChemicalParser(; kwargs...) = ChemicalParser(collect(kwargs))

"""
    ChemicalTransitionParser{T} <: AbstractChemicalParser
    ChemicalTransitionParser(chemicalparser::AbstractChemicalParser = ChemicalExpressionParser)

Default chemical parser that parses an input string into [`ChemicalTransition`](@ref).

The input string is treated as a series of chemicals separated by " -> ". Pairs and vectors can also be parsed into [`ChemicalTransition`](@ref).

If the input can only be parsed into a single chemical entity, it is returned directly without being wrapped in [`ChemicalTransition`](@ref).

# Fields 
* `chemicalparser::T`: chemical parser for each chemical in the transition. Default parser is [`ChemicalExpressionParser`](@ref).
"""
struct ChemicalTransitionParser{T<:AbstractChemicalParser} <: AbstractChemicalParser
    chemicalparser::T
end 
ChemicalTransitionParser() = ChemicalTransitionParser(ChemicalExpressionParser())

"""
    ChemicalExpressionParser{T} <: AbstractChemicalParser
    ChemicalExpressionParser(chemicalparser::AbstractChemicalParser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0)

Default chemical parser that parses an input string into [`AbstractChemical`](@ref) or [`AbstractScheme`](@ref).

# Fields
* `chemicalparser::T`: chemical parser for the core chemical. The default parser is [`FormulaChemicalParser`](@ref).
* `charge::Int`: default chemical entity charge.
* `gain::Int`: default chemical gain charge.
* `loss::Int`: default chemical loss charge.
* `entity::Bool`: whether parse input into a chemical entity.
* `scheme::Bool`: whether parse input into a chemical scheme.

The input is parsed recursively.
First check `entity` and `scheme`. If both `scheme` and `entity` are true, run steps 1–8; if only `scheme` is true, run steps 1–6; if only `entity` is true, run steps 7–8.
1. Match keys of `scheme_name()` -> `AbstractChemicalScheme`
2. `"+core` -> `"+[core]\$(gain)"` -> step 4
3. `"-core` -> `"-[core]\$(loss)"` -> step 5
4. `"+[core+scheme...]n+"`, `"+[core+scheme...]n-"`, `"+[core+scheme...]"` -> step 8 -> `ChemicalGain` of a chemical entity.
5. `"-[core+scheme...]n+"`, `"-[core+scheme...]n-"`, `"-[core+scheme...]"` -> step 8 -> `ChemicalLoss` of a chemical entity.
6. `"[scheme...]n+"`, `"[scheme...]n-"`, `"[scheme...]"` -> parse each scheme
7. `core` -> `"[core]\$(charge)` -> step 8
8. `"[core+scheme...]n+"`, `"[core+scheme...]n-"`, `"[core+scheme...]"` -> parse `core` with `chemicalparser`; parse the rest part recursively -> a chemical entity
"""
struct ChemicalExpressionParser{T<:AbstractChemicalParser} <: AbstractChemicalParser
    chemicalparser::T 
    charge::Int 
    gain::Int 
    loss::Int
    entity::Bool 
    scheme::Bool
end
ChemicalExpressionParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0) = ChemicalExpressionParser(chemicalparser, charge, gain, loss, true, true)
"""
    ChemicalEntityParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0) 

A `ChemicalExpressionParser` parsing only chemical entity.
"""
ChemicalEntityParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0) = ChemicalExpressionParser(chemicalparser, charge, gain, loss, true, false)
"""
    ChemicalSchemeParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0) 

A `ChemicalExpressionParser` parsing only chemical scheme.
"""
ChemicalSchemeParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0) = ChemicalExpressionParser(chemicalparser, charge, gain, loss, false, true)