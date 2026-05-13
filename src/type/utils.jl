"""
    Criteria{A, R}

A type for creating criteria of numeric value. Criterion may be a cutoff value, allowed tolerance relative to an input value, or acceptable range (interval(s)).

# Fields
* `aval`: criterion for the absolute value (unnormalized).
* `rval`: criterion for the relative value (normalized by an input)

# Examples
```julia
julia> peak_crit = Criteria(10, 0.2);

julia> qualified_peak(x, basepeak, crit) = all(<=(x), makecrit_value(crit, basepeak)); # filter out small peaks based on base peak. 

julia> makecrit_value(peak_crit, 100) # create criteria for peak value 100
(10, 20)

julia> makecrit_value(peak_crit, 10) # create criteria for peak value 10
(10, 2)

julia> qualified_peak(10, 100, peak_crit)
false

julia> qualified_peak(20, 100, peak_crit)
true

julia> qualified_peak(10, 10, peak_crit)
true

julia> qualified_peak(2, 10, peak_crit)
false

julia> qc_peak(x, refpeak, crit) = all(in(x), makecrit_delta(crit, refpeak)); # filter out deviated peaks.

julia> makecrit_delta(peak_crit, 100) # create criteria for reference peak value 100
(Intervals.Interval{Float64, Intervals.Closed, Intervals.Closed}(90.0, 110.0), Intervals.Interval{Float64, Intervals.Closed, Intervals.Closed}(80.0, 120.0))

julia> qc_peak(121, 100, peak_crit)
false

"""
struct Criteria{A, R}
    aval::A
    rval::R
end

"""
    value_error(x, y)

Error function of difference; `y - x`.
"""
value_error(x, y) = y - x

"""
    relative_error(x, y)

Error function of relative difference (relative to true value x); `(y - x) / x`.
"""
relative_error(x, y) = (y - x) / x

"""
    relative_error_mean(x, y)

Error function of relative difference (relative to mean); `(y - x) / ((x + y) / 2)`.
"""
relative_error_mean(x, y) = 2 * (y - x) / (x + y) 

"""
    percentage_error(x, y)

Error function of relative difference in percentage (relative to true value x); `(y - x) / x * 100`.
"""
percentage_error(x, y) = (y - x) / x * 100

"""
    percentage_error_mean(x, y)

Error function of relative difference in percentage (relative to mean); `(y - x) / ((x + y) / 2) * 100`.
"""
percentage_error_mean(x, y) = 2 * (y - x) / (x + y) * 100

"""
    ppm_error(x, y)

Error function of relative difference in ppm (relative to true value x); `(y - x) / x * 1e6`.
"""
ppm_error(x, y) = (y - x) / x * 1e6

"""
    ppm_error_mean(x, y)

Error function of relative difference in ppm (relative to mean); `(y - x) / ((x + y) / 2) * 1e6`.
"""
ppm_error_mean(x, y) = 2 * (y - x) / (x + y) * 1e6

abstract type AbstractChemicalParser end
abstract type AbstractAdductParser end

"""
    AdductParser <: AbstractAdductParser

Default adduct parser.
"""
struct AdductParser <: AbstractAdductParser end

"""
    FormulaChemicalParser <: AbstractChemicalParser
    FormulaChemicalParser(; kwargs...)

Default chemical parser which parses input string (chemical formula) into `FormulaChemical`. 

# Fields 
* `property::Vector{Pair{Symbol, Any}}`: additional attributes. These can be muatated, added, or deleted with keyword arguments of `parse_chemical`. Values of `nothing` are deleted.
"""
struct FormulaChemicalParser <: AbstractChemicalParser
    property::Vector{Pair{Symbol, Any}}
end
FormulaChemicalParser(; kwargs...) = FormulaChemicalParser(collect(kwargs))

"""
    ChemicalParser <: AbstractChemicalParser
    ChemicalParser(; kwargs...)

Chemical parser which parses input string (chemical name) into `Chemical`. 

# Fields 
* `property::Vector{Pair{Symbol, Any}}`: additional attributes. It must include `:formula` or `:elements`. These can be muatated, added, or deleted with keyword arguments of `parse_chemical`. Values of `nothing` are deleted.
"""
struct ChemicalParser <: AbstractChemicalParser
    property::Vector{Pair{Symbol, Any}}
end
ChemicalParser(; kwargs...) = ChemicalParser(collect(kwargs))

"""
    AdductIonParser{T, S} <: AbstractChemicalParser
    AdductIonParser([chemicalparser = FormulaChemicalParser()]; adductparser = AdductParser(), charge = 0)

Default chemical parser which parses input string into `AdductIon` if it is charged; otherwise, core chemical is returned.

The input can contain structure like 
* `"[core+adduct]x+"`, `"[core+adduct]x-"`, `"[core]x+"`, `"[core]x-"` for typical `AdductIon`
* `"[core]"` for chemical with no charge
* `"core"` for `AdductIon` with charge determined by field `charge`.

# Fields
* `chemicalparser::T`: chemical parser for core chemical. 
* `adductparser::S: adduct parser for adduct. 
"""
struct AdductIonParser{T <: AbstractChemicalParser, S <: AbstractAdductParser} <: AbstractChemicalParser
    chemicalparser::T
    adductparser::S
    charge::Int
end
AdductIonParser(chemicalparser = FormulaChemicalParser(); adductparser = AdductParser(), charge = 0) = AdductIonParser(chemicalparser, adductparser, charge)

"""
    ChemicalTransitionParser{T} <: AbstractChemicalParser
    ChemicalTransitionParser([chemicalparser = ChemicalGainLossParser()])

Default chemical parser which parses input string into `ChemicalTransition`. 

The input string is regarded as series of chemicals separated by " -> ". Pairs and vectors can also be parsed into `ChemicalTransition`. 

If the input can only be parsed into single chemical entity, it is returned directly without wrapping into `ChemicalTransition`.

# Fields 
* `chemicalparser::T`: chemical parser for each chemical in the transition.
"""
struct ChemicalTransitionParser{T <: AbstractChemicalParser} <: AbstractChemicalParser
    chemicalparser::T
end 
ChemicalTransitionParser() = ChemicalTransitionParser(ChemicalGainLossParser())

"""
    ChemicalGainLossParser{T} <: AbstractChemicalParser
    ChemicalGainLossParser([chemicalparser = AdductIonParser()]; charge = 0, gain = 0, loss = 0)

Default chemical parser which parses input string into `ChemicalGain` or `ChemicalLoss`. 

The input string started with "+" and "-"  is parsed into `ChemicalGain` and `ChemicalLoss`, respectively. For other input, the parsed chemical entity is returned.

# Fields
* `chemicalparser::T`: chemical parser for chemical entity. 
* `charge::Int`: default chemical entity charge.
* `gain::Int`: default chemical gain charge.
* `loss::Int`: default chemical loss charge.
"""
struct ChemicalGainLossParser{T <: AbstractChemicalParser} <: AbstractChemicalParser
    chemicalparser::T
    charge::Int 
    gain::Int 
    loss::Int
end
ChemicalGainLossParser(chemicalparser = AdductIonParser(); charge = 0, gain = 0, loss = 0) = ChemicalGainLossParser(chemicalparser, charge, gain, loss)