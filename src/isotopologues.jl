"""
    Isotopologues(chemical; abundance = 1, abtype = :max, threshold = rcrit(1e-4))
    Isotopologues(formula_name; chemicalparser, abundance = 1, abtype = :max, threshold = rcrit(1e-4))
    Isotopologues(chemicaltransition; abundance = 1, abtype = :max, threshold = rcrit(1e-4)) 
    Isotopologues(formula_name_pair; chemicalparser, abundance = 1, abtype = :max, threshold = rcrit(1e-4)) 
    Isotopologues(tbl; threading = nothing, chemicalparser, threshold = rcrit(1e-4), kwargs...)
    Isotopologues(chemicals; kwargs...)

A `Table` of isotopologues of 
* `chemical::AbstractChemical`: a single chemical entity.
* `formula_name::AbstractString`: a chemical formula or name.
* `chemicaltransition::ChemicalTransition`: MS/MS transition.
* `formula_name_pair::Pair`: MS/MS transition of formulas or names.
* `tbl::Table`: multiple chemicals in column `Chemical` with abundance in column `Abundance1`, `Abundance2`, ... (optional). 
* `chemicals::Vector`: multiple chemicals.

This function is similar to `TandemIsotopologues`; the key difference is that it is iterative and abundance is calculated in the last stage. It performs faster for multiple MS stages and abundance is normalized and filtered at the end.
Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. 

# Keyword Arguments
* `chemicalparser::AbstractChemicalParser`: parser for `formula_name` or `formula_name_pair`. The default parser is `ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0))`.
* `abundance` sets the abundance of the isotope specified by `abtype`. When the input is MS/MS transition, this sets the abundanc of detected chemical. It can also be column `Abundance` in `tbl`.
* `abtype`
    * `:max`: the most abundant isotopologue.
    * `:input`: the input isotopologue.
    * `:list`: sum of listed isotopologues.
    * `:total`: sum of total isotopologues.
* `threshold` can be a number or criteria, representing the lower limit of abundance (absolute and/or relative to maximal value of each spectrum). 
* `threading`: force to use multiple threads (`true`) or single thread (`false`); `nothing` lets the program determine. 

For MS/MS transition, product can be any scheme, including `<:AbstractStructuralScheme`, `<:AbstractStructuralScheme` or formula starting with `-` or `+` (See `ChemicalExpressionParser` for valid string). Gain scheme can only used in the last MS stage.

!!! note "Special precaution for applying to MS/MS precursor-product pairs"
    Product must come from a single part or mutiple non-overlapping parts of precursor. Isobaric or isomeric products are not considered. For instance, 
    * PC 18:0/18:0 and fatty acyl 18:0 fragment is valid because two fatty acids are independent and identical. 
    * PC 18:0[D5]/18:0 and fatty acyl 18:0[D5] fragment is valid but the contribution of another fatty acid 18:0 is not considered and addional computation of this pair and summation with knowledge of fragmentation efficiency are required for the correct abundances. 
    * PC 18:1/18:0 and fatty acyl 18:0 fragment is also valid but requires additional computation of isobaric contribution of another fatty acid 18:1. 

"""
function Isotopologues(input_chemical::AbstractChemical; 
        chemicalparser = ChemicalTransitionParser(),
        id = (1, ), 
        abundance = 1, 
        abtype = Max(), 
        threshold = rcrit(1e-4), 
    ) 
    net_charge = charge(input_chemical)
    it = isotopologues_elements(chemicalelements(input_chemical), first(abundance), abtype, threshold)
    abs_charge = max(1, abs(net_charge))
    net_charge == 0 ? Table(; 
        ID = [id for _ in eachindex(it)], 
        Chemical = [Isotopomers(input_chemical, x) for x in it.Element], 
        Mass1 = it.Mass, 
        Abundance1 = it.Abundance
    ) : 
    Table(; 
        ID = [id for _ in eachindex(it.Element)], 
        Chemical = [Isotopomers(input_chemical, x) for x in it.Element], 
        MZ1 = [m / abs_charge + (net_charge < 0) * ME for m in it.Mass],  
        Abundance1 = it.Abundance
    ) 
end

Isotopologues(input_chemical::AbstractString; 
        chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
        id = (1, ), 
        abundance = 1, 
        abtype = Max(), 
        threshold = rcrit(1e-4)) = 
    Isotopologues(parse_chemical(chemicalparser, input_chemical); id, abundance, abtype, threshold)

function Isotopologues(ct::ChemicalTransition; 
        chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
        id = ntuple(i -> 1, msstage(ct)), 
        abundance = 1, 
        abtype = Max(), 
        threshold = rcrit(1e-4)
    ) 
    abtype = abtyped(abtype)
    msstage(ct) == 1 && return Isotopologues(analyzedchemical(ct); id, abundance, abtype, threshold)
    abundance = float(first(abundance))
    trans = chemicaltransition(ct)
    for c in @view trans[begin:end - 1]
        isgainscheme(c) && throw(ArgumentError("Chemical gain can only be in the last stage for `Isotopologues`; use `TandemIsotopologues` instead."))
    end
    precursor_info = serieschemicaldata(ct)
    if isgainscheme(last(trans))
        if length(trans) > 2
            @inbounds for (pre, post) in @views zip(precursor_info[begin:end - 2], precursor_info[begin + 1:end - 1])
                loss_elements!(last(pre), last(post))
            end
        end
        precursor_info[end] = (first(precursor_info[end]), dictionary_elements(chemicalelements(last(trans))))
    else
        @inbounds for (pre, post) in @views zip(precursor_info[begin:end - 1], precursor_info[begin + 1:end])
            loss_elements!(last(pre), last(post))
        end
    end
    element_dictionary_vec = Vector{Dict}(undef, length(precursor_info))
    msfix_vec = Vector{float(Int)}(undef, length(precursor_info))
    @inbounds for (i, info) in enumerate(precursor_info)
        element_dictionary_vec[i], msfix_vec[i] = get_element_dictinonary_fixmass(last(info))
    end
    it = isotopologues_elements_msn(element_dictionary_vec, msfix_vec, abundance, abtype, threshold) 
    chemical = seriesisotopomerize(trans, it.Element)
    net_charge = [charge(first(p)) for p in precursor_info]
    abs_charge = max.(1, net_charge)
    colab = Symbol(string("Abundance", length(precursor_info))) 
    colmz = all(==(0), net_charge) ? [Symbol(string("Mass", i)) for i in eachindex(precursor_info)] : [Symbol(string("MZ", i)) for i in eachindex(precursor_info)]
    if isgainscheme(last(trans)) 
        for m in it.Mass
            m[end] += m[end - 1]
        end
    end
    mass = [colmz[i] => net_charge[i] == 0 ? [m[i] for m in it.Mass] : [m[i] / abs_charge[i] + (net_charge[i] < 0) * ME for m in it.Mass] for i in eachindex(precursor_info)]
    Table(; ID = [id for _ in eachindex(chemical)], Chemical = chemical, mass..., (colab => it.Abundance, )...) 
end

function Isotopologues(input_chemical::Pair; 
        chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
        id = nothing, 
        abundance = 1, 
        abtype = Max(), 
        threshold = rcrit(1e-4),
    ) 
    ct = parse_chemical(chemicalparser, input_chemical)
    Isotopologues(ct; id = isnothing(id) ? ntuple(i -> 1, msstage(ct)) : id, abundance, abtype, threshold)
end

function Isotopologues(mztable::Table; threading = nothing, chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)), threshold = rcrit(1e-4), kwargs...)
    :Chemical in propertynames(mztable) || throw(ArgumentError("No column `Chemical` in input table."))
    mztable = Table(mztable; Chemical = parse_chemical.(Ref(chemicalparser), mztable.Chemical; kwargs...))
    allequal(msstage, mztable.Chemical) || throw(ArgumentError("Chemicals have to be in the same MS stage."))
    kwargs = Dict(kwargs...)
    vec_key = Symbol[]
    colab = allcolnum(propertynames(mztable), "Abundance"; error = false)
    if !isempty(colab)
        mztable = Table(mztable; abundance = collect.(getproperties(mztable, colab)))
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
    rn = min(length(mztable), Threads.nthreads())
    if isnothing(threading)
        threading = mean(mmi, mztable.Chemical) ^ (msstage(first(mztable.Chemical))) * sqrt(-log2(min(1, minimum(makecrit_value(crit(threshold), 1))))) * (rn - 1) > 100000 
    end
    if threading
        t = Vector{Table}(undef, length(mztable))
        Threads.@threads for i in eachindex(t)
            t[i] = Isotopologues(mztable.Chemical[i]; id = mztable.ID[i], [k => getproperty(mztable, k)[i] for k in vec_key]..., threshold, kwargs...)
        end
        tbl = Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    else
        tbl = Table(vcat((Isotopologues(r.Chemical; id = r.ID, [k => getproperty(r, k) for k in vec_key]..., threshold, kwargs...) for r in mztable)...))
    end
    # spectrum specific threshold ?
    colab = lastcolnum(propertynames(tbl), "Abundance"; error = false)
    ab = getproperty(tbl, colab)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>=(abundance_cutoff), ab)
    tbl[id]
end

Isotopologues(v::Vector; kwargs...) = Isotopologues(Table(; Chemical = v); kwargs...) 
Isotopologues(::Isobars; kwargs...) = throw(ArgumentError("`Isobars` is not supported by `Isotopologues.`"))
Isotopologues(::Isotopomers; kwargs...) = throw(ArgumentError("`Isotopomers` is not supported by `Isotopologues.`"))

"""
    TandemIsotopologues(chemical; chemicalparser, product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4))
    TandemIsotopologues(formula_name; chemicalparser, product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4))
    TandemIsotopologues(chemicaltransition; chemicalparser, product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4)) 
    TandemIsotopologues(formula_name_pair; chemicalparser, product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4)) 
    TandemIsotopologues(tbl; threading = nothing, chemicalparser, threshold = rcrit(1e-4), kwargs...)
    TandemIsotopologues(chemicals; kwargs...)

A `Table` of isotopologues of the following chemicals and their products. 
* `chemical::AbstractChemical`: a single chemical entity.
* `formula_name::AbstractString`: a chemical formula or name.
* `chemicaltransition::ChemicalTransition`: MS/MS transition.
* `formula_name_pair::Pair`: MS/MS transition of formulas or names. 
* `tbl::Table`: multiple chemicals in column `Chemical` with abundance in column `Abundance1`, `Abundance2`, ... (optional). 
* `chemicals::Vector`: multiple chemicals.

This function is similar to `Isotopologues`; the key difference is that it is recursive and abundance is calculated from the begining. It generally performs slightly slower for multiple MS stages and abundance is normalized in the first stage and filtered in all stages.
Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. 

# Keyword Arguments
* `chemicalparser::AbstractChemicalParser`: parser for `formula_name`, `formula_name_pair` or `product`. The default parser is `ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0))`.
* `abundance` sets the abundance of the precursor isotope specified by `abtype`. It can be a vector when the input is MS/MS transition, and abundances are matched to MS stages from the end (the last element matches to the last MS stage). If the length of abundance is smaller, the remaining elemets are filled using `transmission`. 
Notice that setting abundance does not guarantee the equality abundance of the particular isotopes in each MS stage, but the sum of isotopes derived from the particular input precursor isotopes. 
It can also be column `Abundance` in `tbl`.
* `abtype`
    * `:max`: the most abundant isotopologue.
    * `:input`: the input isotopologue.
    * `:list`: sum of listed isotopologues.
    * `:total`: sum of total isotopologues.
* `threshold` can be a number or criteria, representing the lower limit of abundance (absolute and/or relative to maximal value of each spectrum). 
* `product::Vector`: product chemicals. It can also be column `Product` in `tbl`.
* `transmission`: transmission rate between precursors (MS/MS transition). It is utlized when all elements of `abundance` are used out. It can also be column `Transmission` in `tbl`.
* `proportion::Vector`: proportion of fragmentation relative to precursor. It can also be column `Proportion` in `tbl`.
* `threading`: force to use multiple threads (`true`) or single thread (`false`); `nothing` lets the program determine. 

For MS/MS precursor-product pairs, product can be neutral loss/gain or ion loss/gain, including `AbstractElementalScheme`, `AbstractStructuralScheme` and formula starting with `-` or `+` (See `ChemicalExpressionParser` for valid string). 

!!! note "Special precaution for applying to MS/MS precursor-product pairs"
    Product must come from a single part or mutiple non-overlapping parts of precursor. Isobaric or isomeric products are not considered. For instance, 
    * PC 18:0/18:0 and fatty acyl 18:0 fragment is valid because two fatty acids are independent and identical. 
    * PC 18:0[D5]/18:0 and fatty acyl 18:0[D5] fragment is valid but the contribution of another fatty acid 18:0 is not considered and addional computation of this pair and summation with knowledge of fragmentation efficiency are required for the correct abundances. 
    * PC 18:1/18:0 and fatty acyl 18:0 fragment is also valid but requires additional computation of isobaric contribution of another fatty acid 18:1. 

!!! note "Special precaution for applying to MS/MS precursor-product pairs with chemical gain"
    After any chemical gain, the subsequent products are considered randomly fragmented from the gained precursor without considering any structure introduced by chemical gain.
"""
function TandemIsotopologues(input_chemical::AbstractChemical; 
            chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
            abundance = 1, 
            transmission = 1, 
            abtype = Max(), 
            threshold = rcrit(1e-4), 
            id = nothing, 
            precursor_table = nothing, 
            precursor_info = nothing, 
            product = nothing, 
            product_info = nothing, 
            proportion = nothing
        ) 
    id = if isnothing(id) 
        isnothing(precursor_table) ? ntuple(x -> 1, msstage(input_chemical)) : first(precursor_table.ID)
    else
        id 
    end
    (isnothing(product) || isempty(product)) && length(id) == 1 && return Isotopologues(input_chemical; chemicalparser, abundance, abtype, threshold, id)
    precursor_info = isnothing(precursor_info) ? serieschemicaldata(input_chemical) : precursor_info
    precursor_info = vectorize(precursor_info)
    ct = chemicaltransition(input_chemical)
    abundance = vectorize(abundance)
    if length(abundance) > length(ct)
        abundance = abundance[begin:begin + length(ct) - 1]
    elseif length(abundance) < length(ct)
        abundance = vcat(reverse([first(abundance) / transmission ^ i for i in 1:(length(ct) - length(abundance))]), abundance)
    end
    _TandemIsotopologues(precursor_table, ct, precursor_info, id, abundance, lastindex(ct); chemicalparser, abtype, threshold, product, product_info, proportion)
end

function _TandemIsotopologues(precursor_table, transition, precursor_info, id, abundance, ip::Int; 
            chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
            abtype = Max(), 
            threshold = rcrit(1e-4), 
            product = nothing, 
            product_info = nothing, 
            proportion = nothing
        )
    precursor, element_precursor = precursor_info[ip]
    if !isnothing(product) && !isempty(product)
        product = parse_chemical.(Ref(chemicalparser), product)
        all(x -> msstage(x) < 2, product) || throw(ArgumentError("Products should not be MS/MS pairs."))
        proportion = if isnothing(proportion) 
            [1/length(product) for x in eachindex(product)]
        else
            proportion = vectorize(proportion)
        end
        length(proportion) == length(product) || throw(ArgumentError("The length of `proportion` does not mactch the length of `product`"))
        product_info = if isnothing(product_info) || isempty(product_info)
            map(raw_product -> detectedchemicaldata(raw_product; precursor), product) 
        else
            vectorize(product_info)
        end
        length(product_info) == length(product) || throw(ArgumentError("The length of `product_info` does not mactch the length of `product`"))
    end

    if isnothing(precursor_table)
        if ip == firstindex(transition)
            precursor_table = Isotopologues(precursor; chemicalparser, id = id[begin:begin], abundance = abundance[begin], abtype, threshold)
        else
            precursor_table = _TandemIsotopologues(precursor_table, transition, precursor_info, id, abundance, ip - 1; 
                                chemicalparser, abtype, threshold = threshold, 
                                product = nothing, 
                                product_info = nothing, 
                                proportion = nothing)
            colab = Symbol(string("Abundance", ip - 1))
            element_precursor_dictionary = get_element_dictinonary(last(precursor_info[ip - 1]))
            el = map(precursor_table.Chemical) do x 
                first_element_precursor_dictionary = copy(element_precursor_dictionary)
                for (k, v) in detectedisotopes(x)
                    p = parent_element(k)
                    if p != k
                        get!(first_element_precursor_dictionary, k, 0)
                        first_element_precursor_dictionary[k] += v
                        get!(first_element_precursor_dictionary, p, 0)
                        first_element_precursor_dictionary[p] -= v
                    end
                end
                first_element_precursor_dictionary
            end
            itp = (; Element = el, Abundance = getproperty(precursor_table, colab))
            tbl = __TandemIsotopologues(precursor_table, itp, id[begin:ip], precursor_info[ip - 1], precursor_info[ip], transition[ip], abundance[ip] / abundance[ip - 1], abundance[ip], Total(), threshold)
            colab = lastcolnum(propertynames(tbl), "Abundance")
            ab = getproperty(tbl, colab)
            abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
            id = findall(>=(abundance_cutoff), ab)
            precursor_table = tbl[id]
        end
    end
    (isnothing(product) || isempty(product)) && return precursor_table
    colab = Symbol(string("Abundance", length(id)))
    element_precursor_dictionary = get_element_dictinonary(element_precursor)
    el = map(precursor_table.Chemical) do x 
        first_element_precursor_dictionary = copy(element_precursor_dictionary)
        for (k, v) in detectedisotopes(x)
            p = parent_element(k)
            if p != k
                get!(first_element_precursor_dictionary, k, 0)
                first_element_precursor_dictionary[k] += v
                get!(first_element_precursor_dictionary, p, 0)
                first_element_precursor_dictionary[p] -= v
            end
        end
        first_element_precursor_dictionary
    end
    itp = (; Element = el, Abundance = getproperty(precursor_table, colab))
    tbl = vcat([__TandemIsotopologues(precursor_table, itp, (id..., i), precursor_info[ip], prod_info, raw_product, prop, abundance[ip], Total(), threshold; check_product = true) for (i, raw_product, prop, prod_info) in zip(eachindex(product), product, proportion, product_info)]...)
    colab = lastcolnum(propertynames(tbl), "Abundance")
    ab = getproperty(tbl, colab)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>=(abundance_cutoff), ab)
    tbl[id]
end

function __TandemIsotopologues(precursor_table, itp, id, precursor_info, product_info, raw_product, proportion, abundance, abtype, threshold; check_product = false)
    precursor, element_precursor = precursor_info
    product, element_product = product_info
    nms = length(id)
    if check_product && !isgainscheme(raw_product) 
        for (k, v) in element_product
            get(element_precursor, k, 0) < v && throw(ArgumentError("Product can only contain elements restricted by precursor."))
        end
    end
    net_charge = charge(product)
    if isgainscheme(raw_product) 
        element_precursor = chemicalelements(raw_product)
    end
    it = isotopologues_elements_ms2(itp, element_precursor, element_product, abundance, abtype, proportion, threshold; gain = isgainscheme(raw_product), loss = islossscheme(raw_product))
    abs_charge = max(1, net_charge)
    mass = net_charge == 0 ? it.Mass : [m / abs_charge + (net_charge < 0) * ME for m in it.Mass]
    chemical = [ChemicalTransition(precursor_table.Chemical[id], isotopomerize(raw_product, element)) for (id, element) in zip(it.ID, it.Element)]
    abpre = map(1:nms - 1) do i 
        s = Symbol(string("Abundance", i))
        s => [getproperty(precursor_table, s)[id] for id in it.ID]
    end
    if net_charge == 0 
        mspre = map(1:nms - 1) do i 
            s = Symbol(string("Mass", i))
            s => [getproperty(precursor_table, s)[id] for id in it.ID]
        end
        Table(; ID = [id for _ in eachindex(chemical)], Chemical = chemical, mspre..., [Symbol(string("Mass", nms)) => mass]..., abpre..., [Symbol(string("Abundance", nms)) => it.Abundance]...) 
    else
        mspre = map(1:nms - 1) do i 
            s = Symbol(string("MZ", i))
            s => [getproperty(precursor_table, s)[id] for id in it.ID]
        end
        Table(; ID = [id for _ in eachindex(chemical)], Chemical = chemical, mspre..., [Symbol(string("MZ", nms)) => mass]..., abpre..., [Symbol(string("Abundance", nms)) => it.Abundance]...) 
    end
end

TandemIsotopologues(input_chemical::Pair; 
            chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
            abundance = 1, 
            transmission = 1, 
            abtype = Max(), 
            threshold = rcrit(1e-4), 
            id = nothing, 
            precursor_table = nothing, 
            precursor_info = nothing, 
            product = nothing, 
            product_info = nothing, 
            proportion = nothing
        ) = TandemIsotopologues(parse_chemical(chemicalparser, input_chemical); chemicalparser, abundance, transmission, abtype, threshold, id, precursor_table, precursor_info, product, product_info, proportion)

TandemIsotopologues(input_chemical::AbstractString; 
            chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
            abundance = 1, 
            transmission = 1, 
            abtype = Max(), 
            threshold = rcrit(1e-4), 
            id = nothing, 
            precursor_table = nothing, 
            precursor_info = nothing, 
            product = nothing, 
            product_info = nothing, 
            proportion = nothing
        ) = TandemIsotopologues(parse_chemical(chemicalparser, input_chemical); chemicalparser, abundance, transmission, abtype, threshold, id, precursor_table, precursor_info, product, product_info, proportion)

function TandemIsotopologues(mztable::Table; threading = nothing, chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)), threshold = rcrit(1e-4), kwargs...)
    :Chemical in propertynames(mztable) || throw(ArgumentError("No column `Chemical in input table.`"))
    mztable = Table(mztable; Chemical = parse_chemical.(Ref(chemicalparser), mztable.Chemical; kwargs...))
    allequal(msstage, mztable.Chemical) || throw(ArgumentError("Chemicals have to be in the same MS stage."))
    kwargs = Dict(kwargs...)
    vec_key = Symbol[]
    if :Product in propertynames(mztable)
        mztable = Table(mztable; product = mztable.Product)
        push!(vec_key, :product)
        delete!(kwargs, :product)
    elseif !isnothing(get(kwargs, :product, nothing))
        mztable = Table(mztable; product = vectorize(get(kwargs, :product, nothing), length(mztable)))
        push!(vec_key, :product)
        delete!(kwargs, :product)
    end
    # all(x -> all(y -> msstage(y) < 2, x), mztable.Product) || throw(ArgumentError("Products should not MS/MS pairs."))
    colab = allcolnum(propertynames(mztable), "Abundance"; error = false)
    if !isempty(colab)
        mztable = Table(mztable; abundance = collect.(getproperties(mztable, colab)))
        push!(vec_key, :abundance)
        delete!(kwargs, :abundance)
    elseif !isnothing(get(kwargs, :abundance, nothing))
        mztable = Table(mztable; abundance = vectorize(get(kwargs, :abundance, nothing), length(mztable)))
        push!(vec_key, :abundance)
        delete!(kwargs, :abundance)
    end
    if :Transmission in propertynames(mztable)
        mztable = Table(mztable; transmission = mztable.Transmission)
        push!(vec_key, :transmission)
        delete!(kwargs, :transmission)
    elseif !isnothing(get(kwargs, :transmission, nothing))
        mztable = Table(mztable; transmission = vectorize(get(kwargs, :transmission, nothing), length(mztable)))
        push!(vec_key, :transmission)
        delete!(kwargs, :transmission)
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
    rn = min(length(mztable), Threads.nthreads())
    if isnothing(threading)
        threading = mean(mmi, mztable.Chemical) ^ (msstage(first(mztable.Chemical)) + 1) * sqrt(-log2(min(1, minimum(makecrit_value(crit(threshold), 1))))) * (rn - 1) > 100000
    end
    if threading
        t = Vector{Table}(undef, length(mztable))
        Threads.@threads for i in eachindex(t)
            t[i] = TandemIsotopologues(mztable.Chemical[i]; id = mztable.ID[i], [k => getproperty(mztable, k)[i] for k in vec_key]..., threshold, kwargs...)
        end
        tbl = Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    else
        tbl = Table(vcat((TandemIsotopologues(r.Chemical; id = r.ID, [k => getproperty(r, k) for k in vec_key]..., threshold, kwargs...) for r in mztable)...))
    end
    colab = propertynames(tbl)[findlast(x -> startswith(string(x), "Abundance"), propertynames(tbl))]
    ab = getproperty(tbl, colab)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>=(abundance_cutoff), ab)
    tbl[id]
end

TandemIsotopologues(v::Vector; kwargs...) = TandemIsotopologues(Table(; Chemical = v); kwargs...) 
TandemIsotopologues(::Isobars; kwargs...) = throw(ArgumentError("`Isobars` is not supported by `TandemIsotopologues`."))
TandemIsotopologues(::Isotopomers; kwargs...) = throw(ArgumentError("`Isotopomers` is not supported by `TandemIsotopologues`."))

"""
    group_isotopologues(mztable::Table; isotope = "[13C]")

Group isotopologues by isotopomer state based on `isotope`.

* `isotope::String`: minor isotope.
"""
function group_isotopologues(mztable::Table; isotope = "[13C]")
    sp = string.(propertynames(mztable))
    colmz = allcolnum(sp, "MZ")
    colab = allcolnum(sp, "Abundance")
    gf = gf_parent_isotope(isotope)
    transitions = chemicaltransition.(mztable.Chemical)
    gid, gt = unique_group_mz_ab(mztable, transitions, colmz, colab, gf)
    # gt = map(gid) do v 
    #     (; [c => mean(getproperty(mztable, c)[v], weights(getproperty(mztable, d)[v])) for (c, d) in zip(colmz, colab)]..., [c => sum(getproperty(mztable, c)[v]) for c in colab]...)
    # end
    chemcial_parent_state = collect(keys(gid))
    chemical_isotopes = [collect(zip([isotopomersisotopes.(x) for x in transitions[id]]...)) for id in gid]
    chemical_abundance = [[getproperty(mztable, a)[id] for a in colab] for id in gid]
    chemical = [ChemicalSeries([groupedisotopomerize(p..., isotope, collect(i), a) for (p, i, a) in zip(pa, iso, ab)]) for (pa, iso, ab) in zip(chemcial_parent_state, chemical_isotopes, chemical_abundance)]
    Table(Table(; Chemical = chemical), Table(collect(NamedTuple, gt)))
end

"""
    isotopicabundance(chemical::AbstractChemical; ignore_isotopes = false, precise = false)
    isotopicabundance(formula::AbstractString; ignore_isotopes = false, precise = false)
    isotopicabundance(elements::Union{<: Vector, Dict}; ignore_isotopes = false, precise = false)

Compute isotopic abundance of `chemical`, `formula`, vector of element-number pairs or dictionary mapping element to number. 

Parent elements are viewed as major isotopes, and isotopic abundances of all elements are considered in computation. 
To compute isotopic abundance of chemicals with all isotopes labeled intentionally and not following natural distribution, set keyword argument `ignore_isotopes` true, and only parent elements are considered.  
"""
isotopicabundance(cc::AbstractChemical; ignore_isotopes = false, precise = false) = isotopicabundance(chemicalformula(cc); ignore_isotopes, precise)
isotopicabundance(formula::AbstractString; ignore_isotopes = false, precise = false) = isotopicabundance(chemicalelements(formula); ignore_isotopes, precise)
isotopicabundance(elements::Vector{<: Pair}; ignore_isotopes = false, precise = false) = _isotopicabundance(unique_elements(elements); ignore_isotopes, precise)
isotopicabundance(elements::Dict; ignore_isotopes = false, precise = false) = _isotopicabundance(collect(elements); ignore_isotopes, precise)
function _isotopicabundance(elements::Vector{<: Pair}; ignore_isotopes = false, precise = false)
    elements = ignore_isotopes ? filter(iselement ∘ first, elements) : elements
    element_dict = groupfind(parent_element ∘ first, elements)
    abundance_sum = 1
    @inbounds for id in element_dict
        abundance = [get(elements_abundance(), first(elements[i]), 1) for i in id]
        any(==(1), abundance) && continue
        ns = [last(elements[i]) for i in id]
        f = precise ? multinomial(big.(ns)...) : try 
            multinomial(ns...)
        catch 
            multinomial(big.(ns)...)
        end
        abundance_sum *= f * prod(y ^ x for (x, y) in zip(ns, abundance))
    end
    precise ? float(abundance_sum) : convert(float(Int), abundance_sum)
end