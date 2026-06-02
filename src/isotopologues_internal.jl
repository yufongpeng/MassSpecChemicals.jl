function detectedchemicaldata(raw_product; precursor)
    prod = detectedchemical(raw_product; precursor)
    prod, dictionary_elements(filter(x -> last(x) != 0, chemicalelements(prod)))
end

function serieschemicaldata(input_chemical)
    precursor = seriesanalyzedchemical(input_chemical)
    v = map(precursor) do p 
        elements_precursor = chemicalelements(p)
        p, dictionary_elements(filter(x -> last(x) != 0, elements_precursor))
    end
    for (p, d) in v 
        all(>=(0), d) || throw(ArgumentError("Product can only contain elements restricted by precursor."))
    end
    v
end

function chain_chemicals_isotopes(transitions::Vector{<: AbstractChemicalsSchema}, els::Vector{<: Vector{<: Dictionary}})
    @inbounds map(els) do el
        ChemicalTransition(map(enumerate(transitions)) do (i, trans)
            (islossscheme(trans) || isgainscheme(trans)) && i == firstindex(transitions) && throw(ArgumentError("$(tyoeof(trans)) cannot be input chemical."))
            islossscheme(trans) ? isotopomerize(trans, collect(pairs(loss_elements(el[i - 1], el[i])))) : isotopomerize(trans, el[i])
        end
        )
    end
end

function maximal_elements(elements)
    max_dictionary = copy(elements)
    for (e, m) in pairs(elements) 
        x = minor_isotope(e)
        n = floor(Int, m * (1 - elements_abundunce()[e]) / (elements_abundunce()[e] + elements_abundunce()[x]))
        if (max_dictionary[e] - n) / (n + 1) * elements_abundunce()[x] / elements_abundunce()[e] > 1
            n += 1
        end
        max_dictionary[e] -= n 
        set!(max_dictionary, minor_isotope(e), n)
    end
    max_dictionary
end

function initial_proportion(element_product::Dictionary, element_residual::Dictionary; precise = false)
    (any(<(0), element_product) || any(<(0), element_residual)) && return 0
    p = 1
    keproduct = [k for k in keys(element_product) if iselement(k)]
    keresidual = [k for k in keys(element_residual) if iselement(k)]
    # Take care of large numbers
    for e in union(keproduct, keresidual)
        v = [get(element_product, i, 0) + get(element_residual, i, 0) for i in elements_isotopes()[e]]
        lv = get(element_product, e, 0) + get(element_residual, e, 0)
        f = precise ? multinomial(big(lv), v...) : try 
            multinomial(lv, v...)
        catch
            multinomial(big(lv), v...)
        end
        p /= f
    end
    for e in keproduct
        v = [get(element_product, i, 0) for i in elements_isotopes()[e]]
        lv = get(element_product, e, 0)
        f = precise ? multinomial(big(lv), v...) : try 
            multinomial(lv, v...)
        catch
            multinomial(big(lv), v...)
        end
        p *= f
    end 
    for e in keresidual
        v = [get(element_residual, i, 0) for i in elements_isotopes()[e]]
        lv = get(element_residual, e, 0)
        f = precise ? multinomial(big(lv), v...) : try 
            multinomial(lv, v...)
        catch
            multinomial(big(lv), v...)
        end
        p *= f
    end 
    precise ? float(p) : convert(float(Int), p)
end

function maximal_proportion(isotope_precursor, element_product, element_residual; precise = false)
    el = group(parent_element, keys(isotope_precursor))
    p = 1
    for (e, v) in pairs(el) 
        pn = get(element_product, e, 0)
        rn = get(element_residual, e, 0)
        en = pn + rn
        pro = [round(Int, isotope_precursor[x] * pn / en) for x in v]
        pre = [isotope_precursor[x] for x in v]
        if abs(pn + sum(pro) - rn - sum(pre) + sum(pro)) > en - sum(pre)
            f = precise ? multinomial(big(pn), rn) : try 
                multinomial(pn, rn)
            catch 
                multinomial(big(pn), rn)
            end
            p *= f
            for (i, j) in zip(pre, pro)
                f = precise ? multinomial(big(i - j), j) : try 
                    multinomial(i - j, j)
                catch
                    multinomial(big(i - j), j)
                end
                p *= f 
            end
        else
            res = pre .- pro
            f = precise ? multinomial(big(pn), pro...) * multinomial(big(rn), res...) / multinomial(big(en), pre...) : try 
                multinomial(pn - sum(pro), pro...) * multinomial(rn - sum(res), res...) / multinomial(en - sum(pre), pre...)
            catch
                multinomial(big(pn), pro...) * multinomial(big(rn), res...) / multinomial(big(en), pre...)
            end
            p *= f
        end
    end
    precise ? float(p) : convert(float(Int), p)
end

function maximal_proportion_elements(isotope_precursor, element_product, element_residual; precise = false)
    el = group(parent_element, keys(isotope_precursor))
    p = 1
    first_element_product_dictionary = copy(element_product)
    first_element_residual_dictionary = copy(element_residual)
    for (e, v) in pairs(el) 
        pn = get(element_product, e, 0)
        rn = get(element_residual, e, 0)
        en = pn + rn
        pro = [round(Int, isotope_precursor[x] * pn / en) for x in v]
        pre = [isotope_precursor[x] for x in v]
        for (i, x) in enumerate(v) 
            get!(first_element_product_dictionary, x, 0)
            first_element_product_dictionary[x] += pro[i]
            get!(first_element_residual_dictionary, x, 0)
            first_element_residual_dictionary[x] += pre[i] - pro[i]
            get!(first_element_product_dictionary, e, 0)
            first_element_product_dictionary[e] -= pro[i]
            get!(first_element_residual_dictionary, e, 0)
            first_element_residual_dictionary[e] -= pre[i] - pro[i]
        end
        if abs(pn + sum(pro) - rn - sum(pre) + sum(pro)) > en - sum(pre)
            f = precise ? multinomial(big(pn), rn) : try 
                multinomial(pn, rn)
            catch 
                multinomial(big(pn), rn)
            end
            p *= f
            for (i, j) in zip(pre, pro)
                f = precise ? multinomial(big(i - j), j) : try 
                    multinomial(i - j, j)
                catch
                    multinomial(big(i - j), j)
                end
                p *= f 
            end
        else
            res = pre .- pro
            f = precise ? multinomial(big(pn), pro...) * multinomial(big(rn), res...) / multinomial(big(en), pre...) : try 
                multinomial(pn - sum(pro), pro...) * multinomial(rn - sum(res), res...) / multinomial(en - sum(pre), pre...)
            catch
                multinomial(big(pn), pro...) * multinomial(big(rn), res...) / multinomial(big(en), pre...)
            end
            p *= f
        end
    end
    precise ? float(p) : convert(float(Int), p),
    first_element_product_dictionary,
    first_element_residual_dictionary
end

function unique_group_mz_ab(mztable, transitions, colmz, colab, gf = nothing)
    if isnothing(gf)
        gfv = [zeros(length(x)) for x in transitions]
    else
        gfv = map(gf, transitions)
    end
    gids = [groupfind(x -> x[begin:i], gfv) for i in eachindex(first(gfv))]
    ma = map(enumerate(gids)) do (i, gid)
        map(gid) do v 
            c = [x[begin:i] for x in transitions[v]]
            u = v[[findfirst(x -> x == t, c) for t in unique(c)]]
            ab = getproperty(mztable, colab[i])[u]
            mean(getproperty(mztable, colmz[i])[u], weights(ab)), sum(ab)
        end
    end
    gid = last(gids)
    gt = map(keys(gid)) do k
        ks = [k[begin:i] for i in eachindex(gids)]
        ms = [ma[i][ks[i]] for i in eachindex(gids)]
        (; [c => first(m) for (m, c) in zip(ms, colmz)]..., [c => last(m) for (m, c) in zip(ms, colab)]...)
    end
    gid, gt
end

function gf_parent_isotope(isotope = "[13C]")
    isotope_unit = elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]]
    x -> [chemicalparent(m) => _isotopomerstate(isotopomersisotopes(m), isotope_unit) for m in x]
end

function _isotopomerstate(isotopes, isotope_unit)
    ds = 0
    for (e, n) in isotopes
        ds += ustrip((elements_mass()[e] - elements_mass()[elements_parents()[e]]) * n)
    end
    round(Int, ds / ustrip(isotope_unit))
end

function swap_elements(product, residual, swap)
    isempty(swap) && return product
    product2 = copy(product)
    for x in union(keys(product), keys(residual))
        if parent_element(x) in swap
            set!(product2, x, get(residual, x, 0))
        end
    end
    product2
end

get_isotope_vec(input_element::Vector) = filter(x -> !iselement(first(x)), input_element)
get_isotope_vec(input_element::Dictionary) = [k => v for (k, v) in pairs(input_element) if !iselement(k)]
get_element_dictinonary(input_element::Dictionary) = get_element_dictinonary(pairs(input_element))
get_element_dictinonary_fixmass(input_element::Dictionary) = get_element_dictinonary_fixmass(pairs(input_element))
function get_element_dictinonary(input_element)
    element_dictionary = Dictionary{String, Int}()
    for (e, n) in input_element
        if iselement(e)
            get!(element_dictionary, e, 0)
            element_dictionary[e] += n
        end
    end
    element_dictionary
end

function get_element_dictinonary_fixmass(input_element)
    element_dictionary = Dictionary{String, Int}()
    msfix = 0.0u"g"
    for (e, n) in input_element
        if iselement(e)
            get!(element_dictionary, e, 0)
            element_dictionary[e] += n
        else
            msfix += elements_mass()[e] * n
        end
    end
    element_dictionary, ustrip(msfix)
end

const MIN_ABUNDANCE = 1e-4
min_abundance() = MIN_ABUNDANCE
# ==========================================================================================================================
# Mid level MS1
isotopologues_elements(x::AbstractString, abundance, abtype, threshold; table = true, normalize = true) = 
    isotopologues_elements(chemicalelements(x), abundance, abtype, threshold; table, normalize)
isotopologues_elements(input_element::Vector, abundance, abtype, threshold; table = true, normalize = true) = 
    isotopologues_elements(get_element_dictinonary_fixmass(input_element)..., abundance, abtype, threshold; table, normalize)
function isotopologues_elements(element_dictionary::Dictionary, msfix, abundance, abtype, threshold; table = true, normalize = true)
    # element => isoptope pairs
    # remove major
    isempty(element_dictionary) && return (table ? Table(; Element = [get_isotope_vec(element_dictionary)], Mass = [mmi(element_dictionary)], Abundance = [abundance]) : [get_isotope_vec(element_dictionary)])
    element_isotope_pair = mapreduce(vcat, collect(keys(element_dictionary))) do e
        v = map(get(elements_isotopes(), e, e)) do x
            (e, x)
        end
        deleteat!(v, 1)
    end
    sort!(element_isotope_pair; by = x -> elements_abundunce()[last(x)], rev = true)
    first_proportion = isotopicabundance(element_dictionary)
    base_abundance_cutoff = minimum(makecrit_value(crit(threshold), abundance)) / abundance
    abundance_cutoff = base_abundance_cutoff * first_proportion
    if first_proportion < min_abundance()
        max_dictionary = maximal_elements(element_dictionary)
        max_proportion = isotopicabundance(max_dictionary)
        if abtype == :input && abundance_cutoff / max_proportion < min_abundance() 
            throw(ArgumentError("Isotopic abundance of input chemical is too small; try use `abtype` other than `:input` or larger threshold`"))
        elseif abtype != :input 
            abundance_cutoff = base_abundance_cutoff * max_proportion
            first_proportion = isotopicabundance(max_dictionary)
            element_chemical, abundance_chemical = rec_addminusisotopes!([max_dictionary], [abundance * first_proportion], max_dictionary, element_isotope_pair, 1, abundance * first_proportion, abundance_cutoff, true, true)
        else
            element_chemical, abundance_chemical = rec_addisotopes!([element_dictionary], [abundance * first_proportion], element_dictionary, element_isotope_pair, 1, abundance * first_proportion, abundance_cutoff)
        end
    else
        element_chemical, abundance_chemical = rec_addisotopes!([element_dictionary], [abundance * first_proportion], element_dictionary, element_isotope_pair, 1, abundance * first_proportion, abundance_cutoff)
    end
    abundance_chemical_normalize = normalize_abundance(abundance_chemical, abundance, abtype, [:max, :input, :list, :total])
    mass_chemical = map(mmi, element_chemical) .+ msfix
    id = sortperm(mass_chemical)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(abundance_chemical_normalize)))
    filter!(x -> >=(abundance_chemical_normalize[x], abundance_cutoff), id)
    element_chemical = map(get_isotope_vec, @view element_chemical[id]) 
    abundance_chemical = normalize ? abundance_chemical_normalize[id] : abundance_chemical[id]
    mass_chemical = mass_chemical[id]
    table ? Table(; Element = element_chemical, Mass = mass_chemical, Abundance = abundance_chemical) : element_chemical
end

# ==========================================================================================================================
# Mid level MS2
function isotopologues_elements_ms2(it1, element_product, element_residual, abundance, abtype, threshold; gain = false, loss = false, normalize = true, table = true)
    if isempty(element_product)
        data = map(it1) do r 
            (; Element = [element_product], 
            Mass = [0], 
            Abundance = [r.Abundance])
        end
    elseif gain 
        element_residual_dictionary = get_element_dictinonary(element_residual)
        msfix = mmi(element_product) - mmi(element_residual_dictionary)
        proportion_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * isotopicabundance(element_residual_dictionary) / abundance
        it2 = isotopologues_elements(element_residual_dictionary, 0, 1, :total, proportion_cutoff)
        data = map(it1) do r
            (; Element = it2.Element, 
            Mass = it2.Mass, 
            Abundance = it2.Abundance .* r.Abundance)
        end
    else
        element_product_dictionary, msfix = get_element_dictinonary_fixmass(element_product)
        element_residual_dictionary = get_element_dictinonary(element_residual)
        i = findfirst(isempty, it1.Element)
        p = 1
        if isnothing(i) 
            _, i = findmax(it1.Abundance)
            p = maximal_proportion(it1.Element[i], element_product_dictionary, element_residual_dictionary)
        end
        base_proportion_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * it1.Abundance[i] * p / abundance
        data = map(it1) do r
            proportion_cutoff = base_proportion_cutoff / r.Abundance
            element_product_is, proportion_product_is = isotopes_proportion(r.Element, element_product_dictionary, element_residual_dictionary, proportion_cutoff)
            mass_product_is = map(mmi, element_product_is) 
            proportion_product_is .*= r.Abundance
            id = sortperm(mass_product_is)
            el = map(get_isotope_vec, @view element_product_is[id])
            (; Element = loss ? [collect(pairs(loss_elements(r.Element, x))) for x in el] : el, 
            Mass = mass_product_is[id], 
            Abundance = proportion_product_is[id])
        end
    end
    element_product = vcat((getproperty(x, :Element) for x in data)...)
    mass_product = vcat((getproperty(x, :Mass) for x in data)...)
    abundance_pair = vcat((getproperty(x, :Abundance) for x in data)...)
    id_pair = vcat(([i for _ in eachindex(getproperty(x, :Abundance))] for (i, x) in enumerate(data))...)
    abundance_pair_normalize = normalize_abundance(abundance_pair, abundance, abtype, [:max, :list, :total])
    # spectrum specific threshold ?
    abundance_cutoff = isempty(abundance_pair_normalize) ? 0 : minimum(makecrit_value(crit(threshold), maximum(abundance_pair_normalize)))
    id = findall(>=(abundance_cutoff), abundance_pair_normalize)
    id_pair = id_pair[id]
    element_product = element_product[id]
    mass_product = mass_product[id]
    mass_product .+= msfix
    abundance_pair = normalize ? abundance_pair_normalize[id] : abundance_pair[id]
    table ? Table(; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = abundance_pair) : element_product
end

function isotopes_proportion(element_precursor_dictionary::Dictionary, element_product_dictionary::Dictionary, element_residual_dictionary::Dictionary, threshold)
    first_element_product_dictionary = copy(element_product_dictionary)
    first_element_residual_dictionary = copy(element_residual_dictionary)
    swap = String[]
    for (e, n) in pairs(first_element_product_dictionary)
        get!(first_element_residual_dictionary, e, 0)
        if first_element_residual_dictionary[e] < n
            push!(swap, e)
            first_element_product_dictionary[e] = first_element_residual_dictionary[e]
            first_element_residual_dictionary[e] = n
        end
    end
    swap_element_product_dictionary = copy(first_element_product_dictionary)
    swap_element_residual_dictionary = copy(first_element_product_dictionary)
    for (e, n) in pairs(element_precursor_dictionary)
        p = parent_element(e)
        if p != e 
            get!(first_element_residual_dictionary, e, 0)
            first_element_residual_dictionary[e] += n
            get!(first_element_residual_dictionary, p, 0)
            first_element_residual_dictionary[p] -= n
        end 
    end
    # element => isoptope pairs
    # remove first
    element_isotope_pair = mapreduce(vcat, collect(keys(element_product_dictionary))) do e
        v = map(get(elements_isotopes(), e, e)) do x
            (e, x)
        end
        deleteat!(v, 1)
    end
    sort!(element_isotope_pair; by = x -> elements_abundunce()[last(x)], rev = true)
    first_proportion = initial_proportion(first_element_product_dictionary, first_element_residual_dictionary)
    bidir = false
    if first_proportion < min_abundance()
        bidir = true
        first_proportion, first_element_product_dictionary, first_element_residual_dictionary = maximal_proportion_elements(element_precursor_dictionary, swap_element_product_dictionary, swap_element_residual_dictionary)
    end
    first_element_product_dictionary2 = swap_elements(first_element_product_dictionary, first_element_residual_dictionary, swap)
    bidir ? rec_exchangeisotopes!([first_element_product_dictionary2], [first_proportion], first_element_product_dictionary, first_element_residual_dictionary, element_isotope_pair, 1, first_proportion, threshold, swap, true, true) : 
        rec_moveisotopes!([first_element_product_dictionary2], [first_proportion], first_element_product_dictionary, first_element_residual_dictionary, element_isotope_pair, 1, first_proportion, threshold, swap)
end

# ==========================================================================================================================
# Low level MS1
function rec_addisotopes!(
        element_vec::Vector, 
        abundance_vec::Vector, 
        element_dictionary::Dictionary, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_abundance, 
        threshold
    )
    new_isotope_position = -1
    for (e, i) in @inbounds @views element_isotope_pair[isotope_position:end]
        ne = get(element_dictionary, e, 0)
        new_isotope_position += 1
        if ne > 0
            new_element_dictionary = copy(element_dictionary)
            new_element_dictionary[e] -= 1
            get!(new_element_dictionary, i, 0)
            abundance = update_abundance(prev_abundance, e, i, element_dictionary[e], new_element_dictionary[i], 1)
            new_element_dictionary[i] += 1
            if abundance >= threshold  
                push!(element_vec, new_element_dictionary)
                push!(abundance_vec, abundance)
            end
            if abundance >= threshold || abundance >= prev_abundance
                rec_addisotopes!(element_vec, abundance_vec, new_element_dictionary, element_isotope_pair, isotope_position + new_isotope_position, abundance, threshold)
            end
        end
    end
    element_vec, abundance_vec
end

function rec_addminusisotopes!(
        element_vec::Vector, 
        abundance_vec::Vector, 
        element_dictionary::Dictionary, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_abundance, 
        threshold, 
        backward,
        forward
    )
    new_isotope_position = -1
    noe = @inbounds last(element_isotope_pair[isotope_position])
    for (e, i) in @inbounds @views element_isotope_pair[isotope_position:end]
        ne = get(element_dictionary, e, 0)
        ni = get(element_dictionary, i, 0)
        new_isotope_position += 1
        if (i != noe || backward) && ni > 0 
            new_element_dictionary = copy(element_dictionary)
            new_element_dictionary[i] -= 1
            get!(new_element_dictionary, e, 0)
            abundance = update_abundance(prev_abundance, i, e, element_dictionary[i], new_element_dictionary[e], 1)
            new_element_dictionary[e] += 1
            if abundance >= threshold  
                push!(element_vec, new_element_dictionary)
                push!(abundance_vec, abundance)
            end
            if abundance >= threshold || abundance >= prev_abundance
                rec_addminusisotopes!(element_vec, abundance_vec, new_element_dictionary, element_isotope_pair, isotope_position + new_isotope_position, abundance, threshold, true, false)
            end
        end
        if (i != noe || forward) && ne > 0 
            new_element_dictionary = copy(element_dictionary)
            new_element_dictionary[e] -= 1
            get!(new_element_dictionary, i, 0)
            abundance = update_abundance(prev_abundance, e, i, element_dictionary[e], new_element_dictionary[i], 1)
            new_element_dictionary[i] += 1
            if abundance >= threshold  
                push!(element_vec, new_element_dictionary)
                push!(abundance_vec, abundance)
            end
            if abundance >= threshold || abundance >= prev_abundance
                rec_addminusisotopes!(element_vec, abundance_vec, new_element_dictionary, element_isotope_pair, isotope_position + new_isotope_position, abundance, threshold, false, true)
            end
        end
    end
    element_vec, abundance_vec
end

function update_abundance(prev_abundance, old_element, new_element, nold, nnew, delta)
    x = get(elements_abundunce(), old_element, 1)
    y = get(elements_abundunce(), new_element, 1)
    (x == 1 || y == 1) && return prev_abundance
    # take care of large numbers
    f = try 
        factorial(nold, nold - delta) / factorial(nnew + delta, nnew)
    catch 
        factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)
    end
    prev_abundance * f * (y / x) ^ delta
end

# ==========================================================================================================================
# Low level MS2
function rec_moveisotopes!(
        element_vec::Vector, 
        abundance_vec::Vector, 
        element_product_dictionary::Dictionary, 
        element_residual_dictionary::Dictionary, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_proportion,
        threshold,
        swap
    )
    new_isotope_position = -1
    for (e, i) in @inbounds @views element_isotope_pair[isotope_position:end]
        np = get(element_product_dictionary, e, 0)
        nr = get(element_residual_dictionary, i, 0)
        new_isotope_position += 1
        if np > 0 && nr > 0
            new_element_product_dictionary = copy(element_product_dictionary)
            new_element_residual_dictionary = copy(element_residual_dictionary)
            new_element_residual_dictionary[i] -= 1
            new_element_product_dictionary[e] -= 1
            get!(new_element_product_dictionary, i, 0)
            get!(new_element_residual_dictionary, e, 0)
            if prev_proportion == 0
                new_element_product_dictionary[i] += 1
                new_element_residual_dictionary[e] += 1
                proportion = initial_proportion(new_element_product_dictionary, new_element_residual_dictionary)
            else
                proportion = update_proportion(prev_proportion, element_product_dictionary[e], new_element_product_dictionary[i], 1)
                proportion = update_proportion(proportion, element_residual_dictionary[i], new_element_residual_dictionary[e], 1)
                new_element_product_dictionary[i] += 1
                new_element_residual_dictionary[e] += 1
            end
            if proportion >= threshold  
                push!(element_vec, swap_elements(new_element_product_dictionary, new_element_residual_dictionary, swap))
                push!(abundance_vec, proportion)
            end
            if proportion >= threshold || proportion >= prev_proportion
                rec_moveisotopes!(element_vec, abundance_vec, new_element_product_dictionary, new_element_residual_dictionary, element_isotope_pair, isotope_position + new_isotope_position, proportion, threshold, swap)
            end
        end
    end
    element_vec, abundance_vec
end

function rec_exchangeisotopes!(
        element_vec::Vector, 
        abundance_vec::Vector, 
        element_product_dictionary::Dictionary, 
        element_residual_dictionary::Dictionary, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_proportion,
        threshold,
        swap,
        backward, 
        forward
    )
    new_isotope_position = -1
    noe = @inbounds last(element_isotope_pair[isotope_position])
    for (e, i) in @inbounds @views element_isotope_pair[isotope_position:end]
        new_isotope_position += 1
        if (noe != i || backward) && get(element_residual_dictionary, e, 0) > 0 && get(element_product_dictionary, i, 0) > 0
            new_element_product_dictionary = copy(element_product_dictionary)
            new_element_residual_dictionary = copy(element_residual_dictionary)
            new_element_residual_dictionary[e] -= 1
            new_element_product_dictionary[i] -= 1
            get!(new_element_residual_dictionary, i, 0)
            get!(new_element_product_dictionary, e, 0)
            if prev_proportion == 0
                new_element_residual_dictionary[i] += 1
                new_element_product_dictionary[e] += 1
                proportion = initial_proportion(new_element_product_dictionary, new_element_residual_dictionary)
            else
                proportion = update_proportion(prev_proportion, element_residual_dictionary[e], new_element_residual_dictionary[i], 1)
                proportion = update_proportion(proportion, element_product_dictionary[i], new_element_product_dictionary[e], 1)
                new_element_residual_dictionary[i] += 1
                new_element_product_dictionary[e] += 1
            end
            if proportion >= threshold 
                push!(element_vec, swap_elements(new_element_product_dictionary, new_element_residual_dictionary, swap))
                push!(abundance_vec, proportion)
            end
            if proportion >= threshold || proportion >= prev_proportion
                rec_exchangeisotopes!(element_vec, abundance_vec, new_element_product_dictionary, new_element_residual_dictionary, element_isotope_pair, isotope_position + new_isotope_position, proportion, threshold, swap, true, false)
            end
        end
        if (noe != i || forward) && get(element_product_dictionary, e, 0) > 0 && get(element_residual_dictionary, i, 0) > 0
            new_element_product_dictionary = copy(element_product_dictionary)
            new_element_residual_dictionary = copy(element_residual_dictionary)
            new_element_product_dictionary[e] -= 1
            new_element_residual_dictionary[i] -= 1
            get!(new_element_product_dictionary, i, 0)
            get!(new_element_residual_dictionary, e, 0)
            if prev_proportion == 0
                new_element_product_dictionary[i] += 1
                new_element_residual_dictionary[e] += 1
                proportion = initial_proportion(new_element_product_dictionary, new_element_residual_dictionary)
            else
                proportion = update_proportion(prev_proportion, element_product_dictionary[e], new_element_product_dictionary[i], 1)
                proportion = update_proportion(proportion, element_residual_dictionary[i], new_element_residual_dictionary[e], 1)
                new_element_product_dictionary[i] += 1
                new_element_residual_dictionary[e] += 1
            end
            if proportion >= threshold 
                push!(element_vec, swap_elements(new_element_product_dictionary, new_element_residual_dictionary, swap))
                push!(abundance_vec, proportion)
            end
            if proportion >= threshold || proportion >= prev_proportion
                rec_exchangeisotopes!(element_vec, abundance_vec, new_element_product_dictionary, new_element_residual_dictionary, element_isotope_pair, isotope_position + new_isotope_position, proportion, threshold, swap, false, true)
            end
        end
    end
    element_vec, abundance_vec
end

function update_proportion(prev_proportion, nold, nnew, delta; precise = false) 
    f = precise ? factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew) : try 
        factorial(nold, nold - delta) / factorial(nnew + delta, nnew)
    catch 
        factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)
    end
    precise ? prev_proportion * float(f) : prev_proportion * convert(float(Int), f)
end

# function distribute_isotopes!(element_precursor_dictionary::Dictionary, element_product_dictionary::Dictionary, isotope_product_dictionary::Dictionary, element_residual_dictionary::Dictionary, isotope_residual_dictionary::Dictionary)
#     elements = unique(map(parent_element, collect(keys(element_precursor_dictionary))))
#     element_product_vec = [element_product_dictionary]
#     isotope_product_vec = [isotope_product_dictionary]
#     element_residual_vec = [element_residual_dictionary]
#     isotope_residual_vec = [isotope_residual_dictionary]
#     for e in elements
#         isotopes = reverse(filter(i -> i != e && get(element_precursor_dictionary, i, 0) > 0, elements_isotopes()[e]))
#         nisotopes_vec = map(x -> element_precursor_dictionary[x], isotopes)
#         nisotopes_total = sum(nisotopes_vec)
#         if get(last(element_residual_vec), e, 0) >= nisotopes_total 
#             for element_residual in element_residual_vec
#                 get!(element_residual, e, 0)
#                 element_residual[e] -= nisotopes_total 
#             end
#             for isotope in isotopes
#                 for isotope_residual in isotope_residual_vec
#                     get!(isotope_residual, isotope, 0)
#                     isotope_residual[isotope] += element_precursor_dictionary[isotope]
#                 end
#             end
#             continue
#         end
#         new_element_product_vec = empty(element_product_vec)
#         new_isotope_product_vec = empty(isotope_product_vec)
#         new_element_residual_vec = empty(element_residual_vec)
#         new_isotope_residual_vec = empty(isotope_residual_vec)
#         for nisotopes_residual_vec in multiexponents(length(isotopes), get(last(element_residual_vec), e, 0)) 
#             any(x -> <(x...), zip(nisotopes_vec, nisotopes_residual_vec)) && continue 
#             delta = nisotopes_vec .- nisotopes_residual_vec
#             for (isotope_residual, element_residual) in zip(isotope_residual_vec, element_residual_vec)
#                 push!(new_isotope_residual_vec, copy(isotope_residual))
#                 for (x, isotope) in zip(nisotopes_residual_vec, isotopes)
#                     get!(last(new_isotope_residual_vec), isotope, 0)
#                     last(new_isotope_residual_vec)[isotope] += x
#                 end
#                 push!(new_element_residual_vec, copy(element_residual))
#                 get!(last(new_element_residual_vec), e, 0)
#                 last(new_element_residual_vec)[e] = 0
#             end
#             for (isotope_product, element_product) in zip(isotope_product_vec, element_product_vec)
#                 push!(new_isotope_product_vec, copy(isotope_product))
#                 for (d, isotope) in zip(delta, isotopes)
#                     get!(last(new_isotope_product_vec), isotope, 0)
#                     last(new_isotope_product_vec)[isotope] += d
#                 end
#                 push!(new_element_product_vec, copy(element_product))
#                 get!(last(new_element_product_vec), e, 0)
#                 last(new_element_product_vec)[e] -= nisotopes_total - get(last(element_residual_vec), e, 0) - get(last(isotope_residual_vec), e, 0)
#             end
#         end
#         element_product_vec = new_element_product_vec
#         isotope_product_vec = new_isotope_product_vec 
#         element_residual_vec = new_element_residual_vec
#         isotope_residual_vec = new_isotope_residual_vec
#     end
#     element_product_vec, isotope_product_vec, element_residual_vec, isotope_residual_vec
# end