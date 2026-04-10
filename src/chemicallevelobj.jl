"""
    ischemicalequal(x::AbstractChemical, y::AbstractChemical)

Determine whether two chemicals are chemically equivalent. By default, it transforms both chemicals by `ischemicalequaltransform` and compares them by `istransformedchemicalequal`.
"""
ischemicalequal(x::AbstractChemical, y::AbstractChemical) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))
ischemicalequal(x::Isobars, y::Isobars) = all(ischemicalequal(a, b) for (a, b) in zip(x.chemicals, y.chemicals)) && all(isapprox(a, b) for (a, b) in zip(x.abundance, y.abundance))
ischemicalequal(x::Isotopomers, y::Isotopomers) = ischemicalequal(x.parent, y.parent) && isequal(sort!(unique_elements(x.isotopes)), sort!(unique_elements(y.isotopes)))
ischemicalequal(x::ChemicalLoss, y::ChemicalLoss) = ischemicalequal(x.chemical, y.chemical)
ischemicalequal(x::ChemicalGain, y::ChemicalGain) = ischemicalequal(x.chemical, y.chemical)
ischemicalequal(x::ChemicalTransition, y::ChemicalTransition) = all(ischemicalequal.(x.transition, y.transition))

"""
    ischemicalequaltransform(x::AbstractChemical) 

Return an object for comparison with other chemicals by `istransformedchemicalequal`. 
"""
ischemicalequaltransform(x::AbstractChemical) = x 
ischemicalequaltransform(x::Isobars) = length(x) == 1 ? chemicalentity(x) : x
ischemicalequaltransform(x::Isotopomers) = isempty(unique_elements(x.isotopes)) ? x.parent : x 
ischemicalequaltransform(x::ChemicalLoss) = x 
ischemicalequaltransform(x::ChemicalGain) = x 
ischemicalequaltransform(x::ChemicalTransition) = x 

"""
    istransformedchemicalequal(x::AbstractChemical, y::AbstractChemical)
    istransformedchemicalequal(x::Chemical, y::Chemical)

Determine whether two chemicals are chemically equivalent after applying `ischemicalequaltransform`. For `Chemical` and `FormulaChemical`, It tests the name and the elements composition. For others, it defaults to `isequal`.
"""
istransformedchemicalequal(x::AbstractChemical, y::AbstractChemical) = isequal(x, y)
istransformedchemicalequal(x::Chemical, y::Chemical) = 
    isequal(chemicalname(x), chemicalname(y)) && isequal(sort!(unique_elements(chemicalelements(x))), sort!(unique_elements(chemicalelements(y))))
istransformedchemicalequal(x::FormulaChemical, y::FormulaChemical) = 
    isequal(chemicalname(x), chemicalname(y)) 

"""
    istransformedchemicalequal(x::AbstractAdductIon, y::AbstractAdductIon)

Determine whether two adduct ions or isobars are chemically equivalent. The equality relies on both `isadductequal` of adducts and `ischemicalequal` of core chemicals.
"""
istransformedchemicalequal(x::AbstractAdductIon, y::AbstractAdductIon) = isadductequal(ionadduct(x), ionadduct(y)) && ischemicalequal(ioncore(x), ioncore(y))

chemicalentity(isobars::Isobars; kwargs...) = chemicalentity(first(chemicalspecies(isobars)))
chemicalentity(isotopomers::Groupedisotopomers; kwargs...) = Isotopomers(chemicalparent(isotopomers), isotopomersisotopes(isotopomers))
chemicalentity(loss::ChemicalLoss; kwargs...) = loss.chemical
chemicalentity(gain::ChemicalGain; kwargs...) = gain.chemical
chemicalentity(ct::ChemicalTransition; kwargs...) = first(chemicaltransition(ct))

chemicalspecies(isobars::Isobars; kwargs...) = isobars.chemicals

function chemicaltransition(isobars::Isobars{<: ChemicalTransition}; kwargs...) 
    ct = chemicaltransition.(chemicalspecies(isobars))
    [Isobars(getindex.(ct, i), abundance) for (i, abundance) in enumerate(eachcol(isobars.abundance))]
end
chemicaltransition(ct::ChemicalTransition; kwargs...) = ct.transition

chemicalparent(isobars::Isobars; kwargs...) = chemicalparent(chemicalentity(isobars); kwargs...)
chemicalparent(isotopomers::Isotopomers; kwargs...) = isotopomers.parent 
chemicalparent(isotopomers::Groupedisotopomers; kwargs...) = isotopomers.parent 
chemicalparent(loss::ChemicalLoss; kwargs...) = ChemicalLoss(chemicalparent(chemicalentity(loss); kwargs...))
chemicalparent(gain::ChemicalGain; kwargs...) = ChemicalGain(chemicalparent(chemicalentity(gain); kwargs...))
chemicalparent(ct::ChemicalTransition; kwargs...) = ChemicalTransition(chemicalparent.(chemicaltransition(ct); kwargs...))

isotopomersisotopes(isobars::Isobars; kwargs...) = isotopomersisotopes(chemicalentity(isobars); kwargs...)
isotopomersisotopes(isotopomers::Isotopomers; kwargs...) = isotopomers.isotopes
isotopomersisotopes(isotopomers::Groupedisotopomers; kwargs...) = first(isotopomers.isotopes)
isotopomersisotopes(loss::ChemicalLoss; kwargs...) = isotopomersisotopes(chemicalentity(loss); kwargs...) 
isotopomersisotopes(gain::ChemicalGain; kwargs...) = isotopomersisotopes(chemicalentity(gain); kwargs...) 
isotopomersisotopes(ct::ChemicalTransition; kwargs...) = isotopomersisotopes(chemicalentity(ct); kwargs...)

inputchemical(isobars::Isobars; kwargs...) = Isobars([inputchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
inputchemical(ct::ChemicalTransition; kwargs...) = first(chemicaltransition(ct))

outputchemical(isobars::Isobars; kwargs...) = Isobars([outputchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
outputchemical(ct::ChemicalTransition; kwargs...) = last(chemicaltransition(ct))

analyzedchemical(isobars::Isobars; kwargs...) = detectedchemical(isobars; kwargs...)
analyzedchemical(isobars::Isobars{<: ChemicalTransition}; kwargs...) = Isobars([analyzedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
analyzedchemical(loss::ChemicalLoss; kwargs...) = throw(ArgumentError("Chemical loss cannot be directly analyzed."))
analyzedchemical(gain::ChemicalGain; kwargs...) = throw(ArgumentError("Chemical gain cannot be directly analyzed."))
analyzedchemical(ct::ChemicalTransition; kwargs...) = detectedchemical(inputchemical(ct); kwargs...)

function seriesanalyzedchemical(isobars::Isobars{<: ChemicalTransition}; kwargs...) 
    ct = chemicaltransition.(chemicalspecies(isobars))
    [analyzedchemical(Isobars(getindex.(ct, i), abundance); kwargs...) for (i, abundance) in enumerate(eachcol(isobars.abundance))]
end
function seriesanalyzedchemical(ct::ChemicalTransition; kwargs...) 
    v = AbstractChemical[]
    precursor = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedchemical(c; precursor))
        precursor = last(v)
    end
    v
end
function seriesanalyzedisotopes(ct::ChemicalTransition; kwargs...)
    v = Vector{Pair{String, Int}}[]
    precursorisotopes = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedisotopes(c; precursorisotopes))
        precursorisotopes = last(v)
    end
    v
end
function seriesanalyzedcharge(ct::ChemicalTransition; kwargs...)
    v = Int[]
    precursorcharge = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedcharge(c; precursorcharge))
        precursorcharge = last(v)
    end
    v
end
function seriesanalyzedelements(ct::ChemicalTransition; kwargs...)
    v = Vector{Pair{String, Int}}[]
    precursorelements = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedelements(c; precursorelements))
        precursorelements = last(v)
    end
    v
end

detectedchemical(isobars::Isobars; kwargs...) = Isobars([detectedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance)
detectedchemical(isobars::Isobars{<: ChemicalTransition}; kwargs...) = Isobars([detectedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, end])

function detectedchemical(loss::ChemicalLoss; precursor = nothing, kwargs...) 
    isnothing(precursor) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor."))
    precursor = detectedchemical(precursor; kwargs...)
    _detectedproduct(precursor, loss)
end
function detectedisotopes(loss::ChemicalLoss; precursor = nothing, precursorisotopes = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorisotopes) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor information."))
    pre = isnothing(precursorisotopes) ? detectedisotopes(precursor; kwargs...) : precursorisotopes
    loss_elements(pre, isotopomersisotopes(loss; kwargs...))
end
function detectedcharge(loss::ChemicalLoss; precursor = nothing, precursorcharge = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorcharge) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor information."))
    (isnothing(precursorcharge) ? detectedcharge(precursor; kwargs...) : precursorcharge) - charge(loss; kwargs...)

end
function detectedelements(loss::ChemicalLoss; precursor = nothing, precursorelements = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorelements) && throw(ArgumentError("Chemical loss cannot be directly detected without precursor information."))
    loss_elements(isnothing(precursorelements) ? detectedelements(precursor; kwargs...) : precursorelements, chemicalelements(loss; kwargs...))
end

function detectedchemical(gain::ChemicalGain; precursor = nothing, kwargs...) 
    isnothing(precursor) && throw(ArgumentError("Chemical gain cannot be directly detected without precursor."))
    precursor = detectedchemical(precursor; kwargs...)
    _detectedproduct(precursor, gain)
end
function detectedisotopes(gain::ChemicalGain; precursor = nothing, precursorisotopes = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorisotopes) && throw(ArgumentError("Chemical gain cannot be directly detected without precursor information."))
    pre = isnothing(precursorisotopes) ? detectedisotopes(precursor; kwargs...) : precursorisotopes
    gain_elements(pre, isotopomersisotopes(gain; kwargs...))
end
function detectedcharge(gain::ChemicalGain; precursor = nothing, precursorcharge = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorcharge) && throw(ArgumentError("Chemical gain cannot be directly detected without precursor information."))
    (isnothing(precursorcharge) ? detectedcharge(precursor; kwargs...) : precursorcharge) + charge(gain; kwargs...)

end
function detectedelements(gain::ChemicalGain; precursor = nothing, precursorelements = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorelements) && throw(ArgumentError("Chemical gain cannot be directly detected without precursor information."))
    gain_elements(isnothing(precursorelements) ? detectedelements(precursor; kwargs...) : precursorelements, chemicalelements(gain; kwargs...))
end

detectedchemical(ct::ChemicalTransition; kwargs...) = 
    detectedchemical(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...) 
detectedchemical(ct::AbstractVector; kwargs...) = 
    detectedchemical(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) 
detectedisotopes(ct::ChemicalTransition; kwargs...) = 
    detectedisotopes(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...) 
detectedisotopes(ct::AbstractVector; kwargs...) = 
    detectedisotopes(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) 
detectedcharge(ct::ChemicalTransition; kwargs...) = 
    detectedcharge(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...) 
detectedcharge(ct::AbstractVector; kwargs...) = 
    detectedcharge(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) 
detectedelements(ct::ChemicalTransition; kwargs...) = 
    detectedelements(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...) 
detectedelements(ct::AbstractVector; kwargs...) = 
    detectedelements(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) 

    
_detectedproduct(precursor, product) = product 
function _detectedproduct(precursor, product::ChemicalLoss)
    preparent = chemicalparent(precursor)
    postarent = chemicalparent(product)
    prename = chemicalname(preparent)
    preelements = detectedelements(preparent)
    precharge = detectedcharge(precursor)
    preisotopomersisotopes = detectedisotopes(precursor)
    postname = chemicalname(postarent)
    postelements = chemicalelements(postarent)
    postcharge = charge(product)
    postisotopomersisotopes = isotopomersisotopes(product)
    isotopes = unique_elements(loss_elements(preisotopomersisotopes, postisotopomersisotopes))
    chemical = Chemical(string(prename, postname), unique_elements(loss_elements(preelements, postelements)); charge = precharge - postcharge)
    isempty(isotopes) ? chemical : Isotopomers(chemical, isotopes)
end

function _detectedproduct(precursor, product::ChemicalGain)
    preparent = chemicalparent(precursor)
    postarent = chemicalparent(product)
    prename = chemicalname(preparent)
    preelements = detectedelements(preparent)
    precharge = detectedcharge(precursor)
    preisotopomersisotopes = detectedisotopes(precursor)
    postname = chemicalname(postarent)
    postelements = chemicalelements(postarent)
    postcharge = charge(product)
    postisotopomersisotopes = isotopomersisotopes(product)
    isotopes = gain_elements(preisotopomersisotopes, postisotopomersisotopes)
    chemical = Chemical(string(prename, postname), unique_elements(gain_elements(preelements, postelements)); charge = precharge + postcharge)
    isempty(isotopes) ? chemical : Isotopomers(chemical, isotopes)
end

msstage(isobars::Isobars{<: ChemicalTransition}; kwargs...) = only(unique(msstage.(chemicalspecies(isobars); kwargs...)))
msstage(ct::ChemicalTransition; kwargs...) = length(ct.transition)

