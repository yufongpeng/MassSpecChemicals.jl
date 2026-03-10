"""
    Isotopologues(chemical; abundance = 1, abtype = :max, threshold = rcrit(1e-4))
    Isotopologues(formula; charge = 1, abundance = 1, abtype = :max, threshold = rcrit(1e-4))
    Isotopologues(chemicalpair; abundance = 1, abtype = :max, threshold = rcrit(1e-4), proportion = 1) 
    Isotopologues(formulapair; charge = 1, loss = 0, abundance = 1, abtype = :max, threshold = rcrit(1e-4), proportion = 1) 
    Isotopologues(tbl; threshold = rcrit(1e-4), kwargs...)
    Isotopologues(chemicals; kwargs...)

A `Table` of isotopologues of 
* `chemical::AbstractChemical`: a single chemical entity.
* `formula::AbstractString`: a chemical formula.
* `chemicalpair::ChemicalPair`: MS/MS precursor-product pair.
* `formulapair::Pair`: MS/MS precursor-product pair of `formula`s, 
* `tbl::Table`: multiple chemicals in column `Chemical` with abundance in column `Abundance` (optional). 
* `chemicals::Vector`: multiple chemicals.

Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. 

* `abundance` sets the abundance of the isotope specified by `abtype`. 
* `abtype`
    * `:max`: the most abundant isotopologue.
    * `:input`: the input isotopologue.
    * `:list`: sum of listed isotopologues.
    * `:total`: sum of total isotopologues.
* `threshold` can be a number or criteria, representing the lower limit of abundance (absolute and/or relative to maximal value of each spectrum). 
* `proportion`: proportion of fragmentation relative to precursor . 
* `charge` specifys charge for chemicals without intrinsic charge (pure formula). 
* `loss` specifys charge for chemical losses without intrinsic charge (pure formula). 

For MS/MS precursor-product pairs, product can be neutral loss or ion loss (`ChemicalLoss` or formula starting with `-`).

!!! Special precaution for applying to MS/MS precursor-product pairs
    Product must come from a single part or mutiple non-overlapping parts of precursor. Isobaric or isomeric products are not considered. For instance, 
    * PC 18:0/18:0 and fatty acyl 18:0 fragment is valid because two fatty acids are independent and identical. 
    * PC 18:0[D5]/18:0 and fatty acyl 18:0[D5] fragment is valid but the contribution of another fatty acid 18:0 is not considered and addional computation of this pair and summation with knowledge of fragmentation efficiency are required for the correct abundances. 
    * PC 18:1/18:0 and fatty acyl 18:0 fragment is also valid but requires additional computation of isobaric contribution of another fatty acid 18:1. 

"""
function Isotopologues(input_chemical::AbstractChemical; 
        id = (1, ), 
        abundance = 1, 
        abtype = :max, 
        threshold = rcrit(1e-4), 
        proportion = [1]
    ) 
    net_charge = charge(input_chemical)
    it = _isotopologues_elements(chemicalelements(input_chemical), abundance, abtype, threshold, net_charge)
    mass = it.Mass ./ max(1, abs(net_charge))
    chemical = Isotopomers.(input_chemical, it.Element)
    if net_charge == 0 
        Table(; ID = [id for _ in eachindex(chemical)], Chemical = chemical, Mass1 = mass, Abundance1 = it.Abundance) 
    else
        Table(; ID = [id for _ in eachindex(chemical)], Chemical = chemical, MZ1 = mass, Abundance1 = it.Abundance) 
    end
end

Isotopologues(input_chemical::AbstractString; 
        id = (1, ), 
        charge = 1, 
        abundance = 1, 
        abtype = :max, 
        threshold = rcrit(1e-4),
        proportion = [1]) = 
    Isotopologues(ChemicalSeries(input_chemical; charge); id, abundance, abtype, threshold, proportion)

# Isotopologues(formula::Tuple{<: AbstractString, Int}; 
#         id = (1, ), 
#         abundance = 1, 
#         abtype = :max, 
#         threshold = rcrit(1e-4),
#         proportion = 1) = 
#     Isotopologues(ChemicalSeries(formula); id, abundance, abtype, threshold, proportion)

function Isotopologues(cp::ChemicalPair; 
        precursor_table = nothing, 
        id = ntuple(i -> 1, msstage(cp)), 
        precursor = nothing, 
        elements_precursor = nothing, 
        element_precursor = nothing, 
        isotope_precursor = nothing, 
        abundance = 1, 
        abtype = :max, 
        threshold = rcrit(1e-4),
        proportion = [1 for i in 1:msstage(cp)]
    ) 
    if isnothing(precursor)
        precursor = analyzedprecursor(cp)
    end
    if isnothing(elements_precursor)
        elements_precursor = chemicalelements(precursor)
    end
    if isnothing(element_precursor) || isnothing(isotope_precursor)
        element_precursor = filter(x -> haskey(elements_isotopes(), first(x)) && last(x) != 0, elements_precursor)
        isotope_precursor = filter(x -> !haskey(elements_isotopes(), first(x)) && last(x) != 0, elements_precursor)
    end
    product = detectedchemical(cp; precursor)
    elements_product = chemicalelements(product)
    net_charge = charge(product)
    element_product = unique_elements(filter(x -> haskey(elements_isotopes(), first(x)) && last(x) != 0, elements_product))
    isotope_product = unique_elements(filter(x -> !haskey(elements_isotopes(), first(x)) && last(x) != 0, elements_product))
    element_residual = loss_elements(element_precursor, element_product)
    proportion = vectorize(proportion)
    prop = isempty(proportion) ? 1 : pop!(proportion)
    if isnothing(precursor_table)
        precursor_table = Isotopologues(cp.precursor; id = id[begin:end - 1], abundance = abundance / prop, abtype, threshold, proportion)
    end
    nms = length(id)
    colab = Symbol(string("Abundance", nms - 1))
    itp = Table(; Element = [detectedisotopes(x) for x in precursor_table.Chemical], Abundance = getproperty(precursor_table, colab) .* prop)
    it = _isotopologues_elements_ms2(itp, isotope_precursor, element_product, isotope_product, element_residual, abundance, abtype, threshold, net_charge)
    mass = it.Mass ./ max(1, abs(net_charge))
    raw_product = inputchemical(cp.product)
    if raw_product isa ChemicalLoss
        precursors_element = detectedelements.(precursor_table.Chemical)
        chemical = [ChemicalPair(precursor_table.Chemical[r.ID], ChemicalLoss(Isotopomers(raw_product.chemical, loss_elements(precursors_element[r.ID], r.Element)))) for r in it]
    else
        chemical = [ChemicalPair(precursor_table.Chemical[r.ID], Isotopomers(raw_product, r.Element)) for r in it]
    end
    abpre = map(1:nms - 1) do i 
        s = Symbol(string("Abundance", i))
        s => [getproperty(precursor_table, s)[r.ID] for r in it]
    end
    if net_charge == 0 
        mspre = map(1:nms - 1) do i 
            s = Symbol(string("Mass", i))
            s => [getproperty(precursor_table, s)[r.ID] for r in it]
        end
        Table(; ID = [id for _ in eachindex(chemical)], Chemical = chemical, mspre..., [Symbol(string("Mass", nms)) => mass]..., abpre..., [Symbol(string("Abundance", nms)) => it.Abundance]...) 
    else
        mspre = map(1:nms - 1) do i 
            s = Symbol(string("MZ", i))
            s => [getproperty(precursor_table, s)[r.ID] for r in it]
        end
        Table(; ID = [id for _ in eachindex(chemical)], Chemical = chemical, mspre..., [Symbol(string("MZ", nms)) => mass]..., abpre..., [Symbol(string("Abundance", nms)) => it.Abundance]...) 
    end
end

function Isotopologues(input_chemical::Pair; 
        precursor_table = nothing, 
        id = nothing, 
        charge = 1, 
        loss = 0, 
        precursor = nothing, 
        elements_precursor = nothing, 
        element_precursor = nothing, 
        isotope_precursor = nothing, 
        abundance = 1, 
        abtype = :max, 
        threshold = rcrit(1e-4),
        proportion = nothing
    ) 
    cp = ChemicalSeries(input_chemical; charge, loss)
    Isotopologues(cp; precursor_table, id = isnothing(id) ? ntuple(i -> 1, msstage(cp)) : id, precursor, elements_precursor, element_precursor, isotope_precursor, abundance, abtype, threshold, proportion = isnothing(proportion) ? [1 for i in 1:msstage(cp)] : proportion)
end

function Isotopologues(mztable::Table; threshold = rcrit(1e-4), kwargs...)
    :Chemical in propertynames(mztable) || throw(ArgumentError("No column `Chemical` in input table."))
    mztable = Table(mztable; Chemical = ChemicalSeries.(mztable.Chemical; kwargs...))
    allequal(msstage, mztable.Chemical) || throw(ArgumentError("Chemicals have to be in the same MS stage."))
    kwargs = Dict(kwargs...)
    vec_key = Symbol[]
    if :charge in keys(kwargs)
        delete!(kwargs, :charge)
    end
    if :loss in keys(kwargs)
        delete!(kwargs, :loss)
    end
    ab_icol = findlast(x -> startswith(x, "Abundance"), string.(propertynames(mztable)))
    if !isnothing(ab_icol)
        mztable = Table(mztable; abundance = getproperty(mztable, propertynames(mztable)[ab_icol]))
        push!(vec_key, :abundance)
        delete!(kwargs, :abundance)
    elseif !isnothing(get(kwargs, :abundance, nothing))
        mztable = Table(mztable; abundance = vectorize(get(kwargs, :abundance, nothing), length(mztable)))
        push!(vec_key, :abundance)
        delete!(kwargs, :abundance)
    end
    if !in(:ID, propertynames(mztable))
        msn = msstage(first(mztable.Chemical)) - 1
        mztable = Table(mztable; ID = [(i, ntuple(j -> 1, msn)...) for i in eachindex(mztable)])
    end
    # if haskey(kwargs, :threshold)
    #     kwargs[:threshold] = acrit(minimum(makecrit_value(crit(kwargs[:threshold]), maximum(mztable.Abundance))))
    # end

    # if colcharge in propertynames(mztable)
    #     mztable = Table(mztable; net_charge = getproperty(mztable, colcharge))
    #     push!(vec_key, :net_charge)
    #     delete!(kwargs, :net_charge)
    # elseif !isnothing(get(kwargs, :net_charge, nothing))
    #     mztable = Table(mztable; net_charge => vectorize(get(kwargs, :net_charge, nothing), length(mztable)))
    #     push!(vec_key, :net_charge)
    #     delete!(kwargs, :net_charge)
    # end
    # del = Int[]
    # for (i, p) in enumerate(colpreserve)
    #     if !in(p, propertynames(mztable))
    #         @warn "Column `$p` does not exist. Ignore this column."
    #         push!(del, i)
    #     elseif p == :Chemical || p == :MZ || p == :Abundance
    #         @warn "Column `$p` is preserved. Ignore this column."
    #         push!(del, i)
    #     end
    # end
    # colpreserve = deleteat!(collect(colpreserve), del)
    if length(mztable) < Threads.nthreads()
        tbl = Table(vcat((Isotopologues(r.Chemical; id = r.ID, [k => getproperty(r, k) for k in vec_key]..., threshold, kwargs...) for r in mztable)...))
        # mapreduce(vcat, mztable) do r
        #     x = Isotopologues(getproperty(r, colchemical); id = getproperty(r, colid), [k => getproperty(r, k) for k in vec_key]..., kwargs...)
        #     Table(x; (map(colpreserve) do p
        #         p => vectorize(getproperty(r, p), length(x))
        #     end)...)
        # end |> Table
    else
        t = Vector{Table}(undef, length(mztable))
        Threads.@threads for i in eachindex(t)
            t[i] = Isotopologues(mztable.Chemical[i]; id = mztable.ID[i], [k => getproperty(mztable, k)[i] for k in vec_key]..., threshold, kwargs...)
            # x = Isotopologues(getproperty(mztable, colchemical)[i]; [k => getproperty(mztable, k)[i] for k in vec_key]..., kwargs...)
            # t[i] = Table(x; (map(colpreserve) do p
            #     p => vectorize(getproperty(mztable, p)[i], length(x))
            # end)...)
        end
        tbl = Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
    # spectrum specific threshold ?
    colab = propertynames(tbl)[findlast(x -> startswith(x, "Abundance"), string.(propertynames(tbl)))]
    ab = getproperty(tbl, colab)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>=(abundance_cutoff), ab)
    tbl[id]
end

Isotopologues(v::Vector; kwargs...) = Isotopologues(Table(; Chemical = v); kwargs...) 

# function Isotopologues(v::Vector{<: AbstractString}; charge = nothing, kwargs...) 
#     isnothing(charge) ? Isotopologues(Table(; Chemical = v); kwargs...) : 
#         Isotopologues(Table(; Chemical = tuple.(v, charge)); kwargs...)
# end

Isotopologues(::Isobars; kwargs...) = throw(ArgumentError("`Isobars` is not supported by `Isotopologues.`"))
Isotopologues(::Isotopomers; kwargs...) = throw(ArgumentError("`Isotopomers` is not supported by `Isotopologues.`"))
Isotopologues(::ChemicalLoss; kwargs...) = throw(ArgumentError("`ChemicalLoss` is not supported by `Isotopologues.`"))

"""
    TandemIsotopologues(chemical; product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4))
    TandemIsotopologues(formula; charge = 1, product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4))
    TandemIsotopologues(chemicalpair; product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4)) 
    TandemIsotopologues(formulapair; charge = 1, loss = 0, product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4)) 
    TandemIsotopologues(tbl; threshold = rcrit(1e-4), kwargs...)
    TandemIsotopologues(chemicals; kwargs...)

A `Table` of isotopologues of the following chemicals and their products. 
* `chemical::AbstractChemical`: a single chemical entity.
* `formula::AbstractString`: a chemical formula.
* `chemicalpair::ChemicalPair`: MS/MS precursor-product pair.
* `formulapair::Pair`: MS/MS precursor-product pair of `formula`s, 
* `tbl::Table`: multiple chemicals in column `Chemical` with abundance in column `Abundance` (optional). 
* `chemicals::Vector`: multiple chemicals.

Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. 

* `abundance` sets the abundance of the isotope of the product with the largest proportion of fragmentation. The exact isotope is specified by `abtype`. 
* `abtype`
    * `:max`: the most abundant isotopologue.
    * `:input`: the input isotopologue.
    * `:list`: sum of listed isotopologues.
    * `:total`: sum of total isotopologues.
* `threshold` can be a number or criteria, representing the lower limit of abundance (absolute and/or relative to maximal value of each spectrum). 
* `product::Vector`: product chemicals. It can also be column `Product` in `tbl`.
* `proportion::Vector`: proportion of fragmentation relative to precursor. It can also be column `Proportion` in `tbl`.
* `charge` specifys charge for chemicals without intrinsic charge (pure formula). 
* `loss` specifys charge for chemical losses without intrinsic charge (pure formula). 

For MS/MS precursor-product pairs, product can be neutral loss or ion loss (`ChemicalLoss` or formula starting with `-`).

!!! Special precaution for applying to MS/MS precursor-product pairs
    Product must come from a single part or mutiple non-overlapping parts of precursor. Isobaric or isomeric products are not considered. For instance, 
    * PC 18:0/18:0 and fatty acyl 18:0 fragment is valid because two fatty acids are independent and identical. 
    * PC 18:0[D5]/18:0 and fatty acyl 18:0[D5] fragment is valid but the contribution of another fatty acid 18:0 is not considered and addional computation of this pair and summation with knowledge of fragmentation efficiency are required for the correct abundances. 
    * PC 18:1/18:0 and fatty acyl 18:0 fragment is also valid but requires additional computation of isobaric contribution of another fatty acid 18:1. 

"""
function TandemIsotopologues(input_chemical::AbstractChemical; 
            charge = 1, 
            loss = 0, 
            abundance = 1, 
            abtype = :max, 
            threshold = rcrit(1e-4), 
            id = nothing, 
            precursor_table = nothing, 
            product = nothing, 
            proportion = nothing
        )
    isnothing(product) && throw(ArgumentError("Require keyword argument `product`."))
    product = ChemicalSeries.(product; charge, loss)
    all(x -> msstage(x) < 2, product) || throw(ArgumentError("Products should not MS/MS pairs."))
    proportion = if isnothing(proportion) 
        repeat(vcat([1 for i in 1:msstage(input_chemical)], 1/length(product)), length(product))
    else
        [vectorize(v, msstage(input_chemical)) for v in vectorize(proportion, length(product))]
    end
    abundance = vectorize(abundance, length(product)) .* last.(proportion) ./ maximum(last.(proportion))
    precursor = detectedchemical(input_chemical)
    elements_precursor = chemicalelements(input_chemical)
    element_precursor = filter(x -> haskey(elements_isotopes(), first(x)) && last(x) != 0, elements_precursor)
    isotope_precursor = filter(x -> !haskey(elements_isotopes(), first(x)) && last(x) != 0, elements_precursor)
    id = if isnothing(id) 
        isnothing(precursor_table) ? msstage(input_chemical) : first(precursor_table.ID)
    else
        id 
    end
    if isnothing(precursor_table)
        tbl = vcat((Isotopologues(ChemicalPair(input_chemical, product[i]); 
            id = (id..., i),
            precursor_table, 
            abundance = abundance[i],
            proportion = proportion[i], 
            abtype, 
            threshold, 
            precursor, 
            elements_precursor,
            element_precursor,
            isotope_precursor) for i in eachindex(product))...) 
    else
        tbl = vcat((Isotopologues(ChemicalPair(input_chemical, product[i]); 
            id = (id..., i),
            precursor_table = Table(precursor_table), 
            proportion = proportion[i], 
            abtype = :total, 
            threshold, 
            precursor, 
            elements_precursor,
            element_precursor,
            isotope_precursor) for i in eachindex(product))...) 
    end
    # spectrum specific threshold ?
    ab = getproperty(tbl, Symbol(string("Abundance", length(id) + 1)))
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>=(abundance_cutoff), ab)
    tbl[id]
end

TandemIsotopologues(input_chemical::Pair; 
            charge = 1, 
            loss = 0, 
            abundance = 1, 
            abtype = :max, 
            threshold = rcrit(1e-4), 
            id = nothing, 
            precursor_table = nothing, 
            product = nothing, 
            proportion = nothing
        ) = TandemIsotopologues(Formula2Chemical(input_chemical; charge, loss); charge, loss, abundance, abtype, threshold, id, precursor_table, product, proportion)

# TandemIsotopologues(input_chemical::Tuple{<: AbstractString, Int}; 
#             charge = 1, 
#             loss = 0, 
#             abundance = 1, 
#             abtype = :max, 
#             threshold = rcrit(1e-4), 
#             precursor_table = nothing, 
#             product = nothing, 
#             proportion = nothing
#         ) = TandemIsotopologues(Formula2Chemical(input_chemical); charge, loss, abundance, abtype, threshold, precursor_table, product, proportion)

TandemIsotopologues(input_chemical::AbstractString; 
            charge = 1, 
            loss = 0, 
            abundance = 1, 
            abtype = :max, 
            threshold = rcrit(1e-4), 
            id = nothing, 
            precursor_table = nothing, 
            product = nothing, 
            proportion = nothing
        ) = TandemIsotopologues(ChemicalSeries(input_chemical; charge, loss); charge, loss, abundance, abtype, threshold, id, precursor_table, product, proportion)

function TandemIsotopologues(mztable::Table; threshold = rcrit(1e-4), kwargs...)
    :Chemical in propertynames(mztable) || throw(ArgumentError("No column `Chemical in input table.`"))
    mztable = Table(mztable; Chemical = ChemicalSeries.(mztable.Chemical; kwargs...))
    allequal(msstage, mztable.Chemical) || throw(ArgumentError("Chemicals have to be in the same MS stage."))
    kwargs = Dict(kwargs...)
    vec_key = Symbol[]
    # if :charge in keys(kwargs)
    #     delete!(kwargs, :charge)
    # end
    # if :loss in keys(kwargs)
    #     delete!(kwargs, :loss)
    # end
    if :Product in propertynames(mztable)
        mztable = Table(mztable; product = mztable.Product)
        push!(vec_key, :product)
        delete!(kwargs, :product)
    elseif !isnothing(get(kwargs, :product, nothing))
        mztable = Table(mztable; product = vectorize(get(kwargs, :product, nothing), length(mztable)))
        push!(vec_key, :product)
        delete!(kwargs, :product)
    else
        throw(ArgumentError("Require products information. Use column `Product` of the table or keyword arguments `product`."))
    end
    # all(x -> all(y -> msstage(y) < 2, x), mztable.Product) || throw(ArgumentError("Products should not MS/MS pairs."))
    ab_icol = findlast(x -> startswith(x, "Abundance"), string.(propertynames(mztable)))
    if !isnothing(ab_icol)
        mztable = Table(mztable; abundance = getproperty(mztable, propertynames(mztable)[ab_icol]))
        push!(vec_key, :abundance)
        delete!(kwargs, :abundance)
    elseif !isnothing(get(kwargs, :abundance, nothing))
        mztable = Table(mztable; abundance = vectorize(get(kwargs, :abundance, nothing), length(mztable)))
        push!(vec_key, :abundance)
        delete!(kwargs, :abundance)
    end
    if :Proportion in propertynames(mztable)
        mztable = Table(mztable; proportion = mztable.Proportion)
        push!(vec_key, :proportion)
        delete!(kwargs, :proportion)
    elseif !isnothing(get(kwargs, :proportion, nothing))
        mztable = Table(mztable; proportion = vectorize(get(kwargs, :proportion, nothing), length(mztable)))
        push!(vec_key, :proportion)
        delete!(kwargs, :proportion)
    end
    if !in(:ID, propertynames(mztable))
        msn = msstage(first(mztable.Chemical)) - 1
        mztable = Table(mztable; ID = [(i, ntuple(j -> 1, msn)...) for i in eachindex(mztable)])
    end
    if length(mztable) < Threads.nthreads()
        tbl = Table(vcat((TandemIsotopologues(r.Chemical; id = r.ID, [k => getproperty(r, k) for k in vec_key]..., threshold, kwargs...) for r in mztable)...))
    else
        t = Vector{Table}(undef, length(mztable))
        Threads.@threads for i in eachindex(t)
            t[i] = TandemIsotopologues(mztable.Chemical[i]; id = mztable.ID[i], [k => getproperty(mztable, k)[i] for k in vec_key]..., threshold, kwargs...)
        end
        tbl = Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
    colab = propertynames(tbl)[findlast(x -> startswith(x, "Abundance"), string.(propertynames(tbl)))]
    ab = getproperty(tbl, colab)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>=(abundance_cutoff), ab)
    tbl[id]
end

TandemIsotopologues(v::Vector; kwargs...) = TandemIsotopologues(Table(; Chemical = v); kwargs...) 

# function TandemIsotopologues(v::Vector{<: AbstractString}; charge = nothing, kwargs...) 
#     isnothing(charge) ? TandemIsotopologues(Table(; Chemical = v); kwargs...) : 
#         TandemIsotopologues(Table(; Chemical = tuple.(v, charge)); kwargs...)
# end

TandemIsotopologues(::Isobars; kwargs...) = throw(ArgumentError("`Isobars` is not supported by `TandemIsotopologues`."))
TandemIsotopologues(::Isotopomers; kwargs...) = throw(ArgumentError("`Isotopomers` is not supported by `TandemIsotopologues`."))
TandemIsotopologues(::ChemicalLoss; kwargs...) = throw(ArgumentError("`ChemicalLoss` is not supported by `TandemIsotopologues`."))

"""
    isotopicabundance(chemical::AbstractChemical; ignore_isotopes = false)
    isotopicabundance(formula::AbstractString; ignore_isotopes = false)
    isotopicabundance(elements::Union{<: Vector, Dictionary}; ignore_isotopes = false)

Compute isotopic abundance of `chemical`, `formula`, vector of element-number pairs or dictionary mapping element to number. 

Parent elements are viewed as major isotopes, and isotopic abundances of all elements are considered in computation. 
To compute isotopic abundance of chemicals with all isotopes labeled intentionally and not following natural distribution, set keyword argument `ignore_isotopes` true, and only parent elements are considered.  
"""
isotopicabundance(cc::AbstractChemical; ignore_isotopes = false) = isotopicabundance(chemicalformula(cc); ignore_isotopes)
isotopicabundance(chemicalpair::ChemicalPair; ignore_isotopes = false) = isotopicabundance(chemicalformula(chemicalpair.precursor) => chemicalformula(chemicalpair.product); ignore_isotopes)
isotopicabundance(formula::AbstractString; ignore_isotopes = false) = isotopicabundance(chemicalelements(formula); ignore_isotopes)
isotopicabundance(elements::Dictionary; ignore_isotopes = false) = isotopicabundance(pairs(elements); ignore_isotopes)
function isotopicabundance(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}; ignore_isotopes = false)
    element_dict = Dictionary{String, Vector{Int}}()
    # use [12C] for fix 
    for (e, n) in elements
        haskey(elements_parents(), e) || continue
        k = elements_parents()[e]
        haskey(elements_isotopes(), k) || continue
        get!(element_dict, k, zeros(Int, length(elements_isotopes()[k])))
        haskey(elements_isotopes(), e) ? (element_dict[e][begin] += n) : (ignore_isotopes || (element_dict[k][findfirst(==(e), elements_isotopes()[k])] += n))
    end # use two vector?
    abundance_sum = 1
    for (e, ns) in pairs(element_dict)
        abundance = get.(Ref(elements_abundunce()), elements_isotopes()[e], 1)
        any(==(1), abundance) && continue
        n = sum(ns)
        id = sortperm(ns; rev = true)
        ns = ns[id]
        abundance = abundance[id]
        nspop1 = popfirst!(ns)
        f = try 
            factorial(n, nspop1)
        catch 
            factorial(big(n), nspop1)
        end
        abundance_sum *= f * popfirst!(abundance) ^ nspop1
        isempty(abundance) && continue
        # take care if big number
        abundance_sum *= mapreduce(*, ns, abundance) do x, y
            y ^ x / factorial(x)
        end
    end
    Float64(abundance_sum)
end

# ==========================================================================================================================
# Internal
_isotopologues_elements(x::AbstractString, abundance, abtype, threshold, net_charge; table = true, normalize = true) = 
    _isotopologues_elements(chemicalelements(x), abundance, abtype, threshold, net_charge; table, normalize)
function _isotopologues_elements(input_element::Vector, abundance, abtype, threshold, net_charge; table = true, normalize = true)
    # record elements change
    element_dictionary = Dictionary{String, Int}()
    isotope_dictionary = Dictionary{String, Int}()
    fix_dictionary = Dictionary{String, Int}()
    first_element_dictionary = Dictionary{String, Int}()
    for (e, n) in input_element
        if haskey(elements_isotopes(), e)
            get!(element_dictionary, e, 0)
            element_dictionary[e] += n
            get!(first_element_dictionary, e, 0)
            first_element_dictionary[e] += n
        elseif first(elements_isotopes()[get(elements_parents(), e, e)]) == e 
            p = get(elements_parents(), e, e)
            get!(fix_dictionary, p, 0)
            fix_dictionary[p] += n
            get!(element_dictionary, p, 0)
            element_dictionary[p] += n
            get!(first_element_dictionary, p, 0)
            first_element_dictionary[p] += n
        else 
            p = get(elements_parents(), e, e)
            get!(fix_dictionary, e, 0)
            fix_dictionary[e] += n
            get!(isotope_dictionary, e, 0)
            isotope_dictionary[e] += n
            get!(first_element_dictionary, e, 0)
            first_element_dictionary[e] += n
        end
    end
    # element => isoptope pairs
    # remove first
    element_isotope_pair = mapreduce(vcat, collect(keys(element_dictionary))) do e
        v = map(get(elements_isotopes(), e, e)) do x
            (e, x)
        end
        deleteat!(v, 1)
    end
    sort!(element_isotope_pair; by = x -> elements_abundunce()[last(x)], rev = true)
    first_proportion = isotopicabundance(input_element; ignore_isotopes = true)
    # abundance_cutoff_normalize = @match abtype begin
    #     :input => abundance_cutoff * first_proportion
    #     :max  => abundance_cutoff * first_proportion
    #     _     => abundance_cutoff
    # end
    abundance_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * first_proportion
    # serve abundance as sums
    element_chemical, abundance_chemical = rec_addisotopes!([first_element_dictionary], [abundance * first_proportion], element_dictionary, isotope_dictionary, isotope_dictionary, element_isotope_pair, 1, abundance * first_proportion, abundance_cutoff)
    # Normalize after resolution check?
    abundance_chemical_normalize = abtype == :max ? abundance_chemical ./ maximum(abundance_chemical) .* abundance : 
                                    abtype == :input ? abundance_chemical ./ first(abundance_chemical) .* abundance : 
                                    abtype == :list ? abundance_chemical ./ sum(abundance_chemical) .* abundance : abundance_chemical
    mass_chemical = map(mmi, element_chemical, repeat([net_charge], length(element_chemical)))
    id = sortperm(mass_chemical)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(abundance_chemical_normalize)))
    filter!(x -> >=(abundance_chemical_normalize[x], abundance_cutoff), id)
    element_chemical = element_chemical[id]
    abundance_chemical = normalize ? abundance_chemical_normalize[id] : abundance_chemical[id]
    mass_chemical = mass_chemical[id]
    table ? Table(; Element = element_chemical, Mass = mass_chemical, Abundance = abundance_chemical) : element_chemical
end

function _isotopologues_elements_ms2(it1, isotope_precursor, element_product, isotope_product, element_residual, abundance, abtype, threshold, net_charge; normalize = true, table = true)
    msfix = mmi(isotope_product, net_charge)
    data = map(it1) do r 
        element_product_is, mass_product_is, proportion_product_is = isotopes_proportion(loss_elements(r.Element, isotope_precursor), element_product, element_residual, msfix, net_charge)
        id = sortperm(mass_product_is)
        (; Element = add_elements.(Ref(isotope_product), element_product_is[id]), 
           Mass = mass_product_is[id], 
           Abundance = r.Abundance .* proportion_product_is[id])
    end
    element_product = vcat((getproperty(x, :Element) for x in data)...)
    mass_product = vcat((getproperty(x, :Mass) for x in data)...)
    abundance_pair = vcat((getproperty(x, :Abundance) for x in data)...)
    id_pair = vcat(([i for _ in eachindex(getproperty(x, :Abundance))] for (i, x) in enumerate(data))...)
    abundance_pair_normalize = abtype == :max ? abundance_pair ./ maximum(abundance_pair) .* abundance : 
                                abtype == :input ? abundance_pair ./ first(abundance_pair) .* abundance : 
                                abtype == :list ? abundance_pair ./ sum(abundance_pair) .* abundance : abundance_pair
    # spectrum specific threshold ?
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(abundance_pair_normalize)))
    id = findall(>=(abundance_cutoff), abundance_pair_normalize)
    id_pair = id_pair[id]
    element_product = element_product[id]
    mass_product = mass_product[id]
    abundance_pair = normalize ? abundance_pair_normalize[id] : abundance_pair[id]
    table ? Table(; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = abundance_pair) : element_product
end

function rec_addisotopes!(element_vec::Vector, abundance_vec::Vector, element_dictionary::Dictionary, isotope_dictionary::Dictionary, fix_dictionary::Dictionary, element_isotope_pair::Vector, isotope_position::Int, prev_abundance, threshold)
    new_isotope_position = -1
    for (e, i) in @views element_isotope_pair[isotope_position:end]
        ne = get(element_dictionary, e, 0)
        new_isotope_position += 1
        if ne > get(fix_dictionary, e, 0)
            new_element_dictionary = deepcopy(element_dictionary)
            new_isotope_dictionary = deepcopy(isotope_dictionary)
            new_element_dictionary[e] -= 1
            get!(new_isotope_dictionary, i, 0)
            abundance = update_abundance(prev_abundance, e, i, element_dictionary[e] - get(fix_dictionary, e, 0), new_isotope_dictionary[i] - get(fix_dictionary, i, 0), 1)
            new_isotope_dictionary[i] += 1
            abundance < threshold && break
            push!(element_vec, add_elements(new_element_dictionary, new_isotope_dictionary))
            push!(abundance_vec, abundance)
            rec_addisotopes!(element_vec, abundance_vec, new_element_dictionary, new_isotope_dictionary, fix_dictionary, element_isotope_pair, isotope_position + new_isotope_position, abundance, threshold)
        end
    end
    element_vec, abundance_vec
end

function update_abundance(prev_abundance, old_element, new_element, nold, nnew, delta)
    x = get(elements_abundunce(), old_element, 1)
    y = get(elements_abundunce(), new_element, 1)
    (x == 1 || y == 1) && return prev_abundance
    # take care of large numbers
    prev_abundance * factorial(nold, nold - delta) / factorial(nnew + delta, nnew) * (y / x) ^ delta
end

# Distribute isotopes: setdiff without deleting overlapped/calculate 1 fragment * n
function isotopes_proportion(element_precursor_dictionary::Dictionary, element_product_dictionary::Dictionary, element_residual_dictionary::Dictionary, msfix, net_charge)
    first_isotope_product_dictionary = Dictionary{String, Int}()
    first_isotope_residual_dictionary = Dictionary{String, Int}()
    first_element_product_dictionary = deepcopy(element_product_dictionary)
    first_element_residual_dictionary = deepcopy(element_residual_dictionary)
    # element => isoptope pairs
    # remove first
    element_isotope_pair = mapreduce(vcat, collect(keys(element_product_dictionary))) do e
        v = map(get(elements_isotopes(), e, e)) do x
            (e, x)
        end
        deleteat!(v, 1)
    end
    sort!(element_isotope_pair; by = x -> elements_abundunce()[last(x)], rev = true)
    element_product_vec, isotope_product_vec, element_residual_vec, isotope_residual_vec = distribute_isotopes!(element_precursor_dictionary, first_element_product_dictionary, first_isotope_product_dictionary, first_element_residual_dictionary, first_isotope_residual_dictionary)
    # initial proportion
    element_pair = Dictionary[]
    proportion_pair = Float64[]
    for (element_product, isotope_product, element_residual, isotope_residual) in zip(element_product_vec, isotope_product_vec, element_residual_vec, isotope_residual_vec)
        proportion = initial_proportion(element_product, isotope_product, element_residual, isotope_residual)
        element_pair_vec, proportion_pair_vec = rec_moveisotopes!([add_elements(element_product, isotope_product)], [proportion], element_product, isotope_product, element_residual, isotope_residual, element_isotope_pair, 1, proportion)
        for (element_pairi, proportion_pairi) in zip(element_pair_vec, proportion_pair_vec) 
            any(==(element_pairi), element_pair) && continue
            push!(element_pair, element_pairi)
            push!(proportion_pair, proportion_pairi)
        end
    end
    mass_pair = map(mmi, element_pair) .+ msfix
    id = sortperm(mass_pair)
    element_pair = element_pair[id]
    proportion_pair = proportion_pair[id]
    mass_pair = mass_pair[id]
    element_pair, mass_pair, proportion_pair
end

function distribute_isotopes!(element_precursor_dictionary::Dictionary, element_product_dictionary::Dictionary, isotope_product_dictionary::Dictionary, element_residual_dictionary::Dictionary, isotope_residual_dictionary::Dictionary)
    elements = unique(map(e -> get(elements_parents(), e, e), collect(keys(element_precursor_dictionary))))
    element_product_vec = [element_product_dictionary]
    isotope_product_vec = [isotope_product_dictionary]
    element_residual_vec = [element_residual_dictionary]
    isotope_residual_vec = [isotope_residual_dictionary]
    for e in elements
        isotopes = filter(i -> i != e && get(element_precursor_dictionary, i, 0) > 0, elements_isotopes()[e])
        nisotopes_vec = map(x -> element_precursor_dictionary[x], isotopes)
        nisotopes_total = sum(nisotopes_vec)
        if get(last(element_residual_vec), e, 0) >= nisotopes_total 
            for element_residual in element_residual_vec
                get!(element_residual, e, 0)
                element_residual[e] -= nisotopes_total 
            end
            for isotope in isotopes
                for isotope_residual in isotope_residual_vec
                    get!(isotope_residual, isotope, 0)
                    isotope_residual[isotope] += element_precursor_dictionary[isotope]
                end
            end
            continue
        end
        new_element_product_vec = empty(element_product_vec)
        new_isotope_product_vec = empty(isotope_product_vec)
        new_element_residual_vec = empty(element_residual_vec)
        new_isotope_residual_vec = empty(isotope_residual_vec)
        for nisotopes_residual_vec in multiexponents(length(isotopes), get(last(element_residual_vec), e, 0)) 
            any(x -> <(x...), zip(nisotopes_vec, nisotopes_residual_vec)) && continue 
            delta = nisotopes_vec .- nisotopes_residual_vec
            for isotope_residual in isotope_residual_vec
                push!(new_isotope_residual_vec, deepcopy(isotope_residual))
                for (x, isotope) in zip(nisotopes_residual_vec, isotopes)
                    get!(last(new_isotope_residual_vec), isotope, 0)
                    last(new_isotope_residual_vec)[isotope] += x
                end
            end
            for isotope_product in isotope_product_vec
                push!(new_isotope_product_vec, deepcopy(isotope_product))
                for (d, isotope) in zip(delta, isotopes)
                    get!(last(new_isotope_product_vec), isotope, 0)
                    last(new_isotope_product_vec)[isotope] += d
                end
            end
        end
        for element_product in element_product_vec
            for _ in eachindex(new_isotope_product_vec)
                push!(new_element_product_vec, deepcopy(element_product))
                get!(last(new_element_product_vec), e, 0)
                last(new_element_product_vec)[e] -= nisotopes_total - get(last(element_residual_vec), e, 0) - get(last(isotope_residual_vec), e, 0)
            end
        end
        for element_residual in element_residual_vec
            for _ in eachindex(new_isotope_product_vec)
                push!(new_element_residual_vec, deepcopy(element_residual))
                get!(last(new_element_residual_vec), e, 0)
                last(new_element_residual_vec)[e] = 0
            end
        end
        element_product_vec = new_element_product_vec
        isotope_product_vec = new_isotope_product_vec 
        element_residual_vec = new_element_residual_vec
        isotope_residual_vec = new_isotope_residual_vec
    end
    element_product_vec, isotope_product_vec, element_residual_vec, isotope_residual_vec
end

function initial_proportion(element_product::Dictionary, isotope_product::Dictionary, element_residual::Dictionary, isotope_residual::Dictionary)
    p = 1
    element_total = add_elements(element_product, element_residual)
    isotope_total = add_elements(isotope_product, isotope_residual)
    # Take care of large numbers
    for e in keys(element_total)
        n = 0
        for i in ISOTOPES[e]
            m = get(isotope_total, i, 0)
            n += m
            p *= factorial(m)
        end
        p /= factorial(n + get(element_total, e, 0), get(element_total, e, 0))
    end
    for e in keys(element_product)
        n = 0
        for i in ISOTOPES[e]
            m = get(isotope_product, i, 0)
            n += m
            p /= factorial(m)
        end
        p *= factorial(n + get(element_product, e, 0), get(element_product, e, 0))
    end 
    for e in keys(element_residual)
        n = 0
        for i in ISOTOPES[e]
            m = get(isotope_residual, i, 0)
            n += m
            p /= factorial(m)

        end
        p *= factorial(n + get(element_residual, e, 0), get(element_residual, e, 0))
    end 
    p
end

# ab as proportion
function rec_moveisotopes!(element_vec::Vector, abundance_vec::Vector, element_product_dictionary::Dictionary, isotope_product_dictionary::Dictionary, element_residual_dictionary::Dictionary, isotope_residual_dictionary::Dictionary, element_isotope_pair::Vector, isotope_position::Int, prev_proportion)
    new_isotope_position = -1
    for (e, i) in @views element_isotope_pair[isotope_position:end]
        np = get(element_product_dictionary, e, 0)
        nr = get(isotope_residual_dictionary, i, 0)
        new_isotope_position += 1
        if np > 0 && nr > 0
            new_element_product_dictionary = deepcopy(element_product_dictionary)
            new_isotope_product_dictionary = deepcopy(isotope_product_dictionary)
            new_element_residual_dictionary = deepcopy(element_residual_dictionary)
            new_isotope_residual_dictionary = deepcopy(isotope_residual_dictionary)
            new_element_product_dictionary[e] -= 1
            new_isotope_residual_dictionary[i] -= 1
            get!(new_isotope_product_dictionary, i, 0)
            get!(new_element_residual_dictionary, e, 0)
            proportion = update_proportion(prev_proportion, element_product_dictionary[e], new_isotope_product_dictionary[i], 1)
            proportion = update_proportion(proportion, isotope_residual_dictionary[i], new_element_residual_dictionary[e], 1)
            new_isotope_product_dictionary[i] += 1
            new_element_residual_dictionary[e] += 1
            push!(element_vec, add_elements(new_element_product_dictionary, new_isotope_product_dictionary))
            push!(abundance_vec, proportion)
            rec_moveisotopes!(element_vec, abundance_vec, new_element_product_dictionary, new_isotope_product_dictionary, new_element_residual_dictionary, new_isotope_residual_dictionary, element_isotope_pair, isotope_position + new_isotope_position, proportion)
        end
    end
    element_vec, abundance_vec
end

update_proportion(prev_proportion, nold, nnew, delta) = prev_proportion * factorial(nold, nold - delta) / factorial(nnew + delta, nnew)