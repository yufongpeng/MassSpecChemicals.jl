"""
    chemicalname(m::T; kwargs...) 

Get name of `m`. It is equivalent to `getchemicalattr(m, :name; kwargs...)` except that it returns `string("::", T)` when name is not available.
"""
function chemicalname(m::T; kwargs...) where {T <: AbstractChemical}
    result = getchemicalattr(m, :name; kwargs...)
    isnothing(result) ? string("::", T) : result
end

"""
    chemicalformula(m::AbstractChemical; kwargs...)

Get formula of `m`. It is equivalent to `getchemicalattr(m, :formula; kwargs...)` except that it returns `""` when formula is not available.
"""
function chemicalformula(m::AbstractChemical; kwargs...) 
    result = getchemicalattr(m, :formula; kwargs...)
    isnothing(result) ? "" : result
end
"""
    chemicalabbr(m::AbstractChemical; kwargs...)

Get abbreviation of `m`. It is equivalent to `getchemicalattr(m, :abbreviation; kwargs...)` except that it returns `chemicalname(m; kwargs...)` when abbreviation is not available.
"""
function chemicalabbr(m::AbstractChemical; kwargs...) 
    result = getchemicalattr(m, :abbreviation; kwargs...)
    isnothing(result) ? chemicalname(m; kwargs...) : result
end
"""
    chemicalsmiles(m::AbstractChemical; kwargs...)

Get SMILES of `m`. It is equivalent to `getchemicalattr(m, :SMILES; kwargs...)` except that it returns `""` when SMILES is not available.
"""
function chemicalsmiles(m::AbstractChemical; kwargs...) 
    result = getchemicalattr(m, :SMILES; kwargs...)
    isnothing(result) ? "" : result
end
"""
    rt(m::AbstractChemical; kwargs...)

Get retention time of `m`. It is equivalent to `getchemicalattr(m, :rt; kwargs...)` except that it returns `NaN` when rt is not available.
"""
function rt(m::AbstractChemical; kwargs...) 
    result = getchemicalattr(m, :rt; kwargs...)
    isnothing(result) ? NaN : result 
end

