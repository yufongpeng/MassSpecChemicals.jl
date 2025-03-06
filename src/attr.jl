"""
    chemicalname(chemical::T; kwargs...) 

Get name of `chemical`. It is equivalent to `getchemicalattr(chemical, :name; kwargs...)` except that it returns `string("::", T)` when name is not available.
"""
function chemicalname(cc::T; kwargs...) where {T <: AbstractChemical}
    result = getchemicalattr(cc, :name; kwargs...)
    isnothing(result) ? string("::", T) : result
end

"""
    chemicalformula(chemical::AbstractChemical; kwargs...)

Get formula of `chemical`. It is equivalent to `getchemicalattr(chemical, :formula; kwargs...)` except that it returns `""` when formula is not available.
"""
function chemicalformula(cc::AbstractChemical; kwargs...) 
    result = getchemicalattr(cc, :formula; kwargs...)
    isnothing(result) ? "" : result
end

"""
    chemicalelements(chemical::AbstractChemical; kwargs...)

Get elements of `chemical`. It is equivalent to `getchemicalattr(chemical, :elements; kwargs...)` except that it returns `Pair{String, Int}[]` when elements are not available.
"""
function chemicalelements(cc::AbstractChemical; kwargs...) 
    result = getchemicalattr(cc, :elements; kwargs...)
    isnothing(result) ? Pair{String, Int}[] : result
end

"""
    chemicalabbr(chemical::AbstractChemical; kwargs...)

Get abbreviation of `chemical`. It is equivalent to `getchemicalattr(chemical, :abbreviation; kwargs...)` except that it returns `chemicalname(chemical; kwargs...)` when abbreviation is not available.
"""
function chemicalabbr(cc::AbstractChemical; kwargs...) 
    result = getchemicalattr(cc, :abbreviation; kwargs...)
    isnothing(result) ? chemicalname(m; kwargs...) : result
end

"""
    chemicalsmiles(chemical::AbstractChemical; kwargs...)

Get SMILES of `chemical`. It is equivalent to `getchemicalattr(chemical, :SMILES; kwargs...)` except that it returns `""` when SMILES is not available.
"""
function chemicalsmiles(cc::AbstractChemical; kwargs...) 
    result = getchemicalattr(cc, :SMILES; kwargs...)
    isnothing(result) ? "" : result
end

"""
    rt(chemical::AbstractChemical; kwargs...)

Get retention time of `chemical`. It is equivalent to `getchemicalattr(chemical, :rt; kwargs...)` except that it returns `NaN` when rt is not available.
"""
function rt(cc::AbstractChemical; kwargs...) 
    result = getchemicalattr(cc, :rt; kwargs...)
    isnothing(result) ? NaN : result 
end
