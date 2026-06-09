"""
    detectedchemicaldata(raw_product; precursor)

Detected chemical and dictionary of elements.
"""
function detectedchemicaldata(raw_product; precursor)
    prod = detectedchemical(raw_product; precursor)
    prod, dictionary_elements(filter(x -> last(x) != 0, chemicalelements(prod)))
end

"""
    serieschemicaldata(raw_product; precursor)

Series detected chemical and dictionaries of elements.
"""
function serieschemicaldata(input_chemical)
    precursor = seriesanalyzedchemical(input_chemical)
    v = map(precursor) do p 
        elements_precursor = chemicalelements(p)
        p, dictionary_elements(filter(x -> last(x) != 0, elements_precursor))
    end
    for (p, d) in v 
        all(>=(0), values(d)) || throw(ArgumentError("Product can only contain elements restricted by precursor."))
    end
    v
end

"""
    seriesisotopomerize(transitions::Vector{<: AbstractChemicalsSchema}, els::Vector{<: Vector{<: Dict}})

Serial isotomoperize `transitions` with detected isotpic replacements `els`.
"""
function seriesisotopomerize(transitions::Vector{<: AbstractChemicalsSchema}, els::Vector{<: Vector{<: Dict}})
    @inbounds map(els) do el
        ChemicalTransition(map(enumerate(transitions)) do (i, trans)
            (islossscheme(trans) || isgainscheme(trans)) && i == firstindex(transitions) && throw(ArgumentError("$(typeof(trans)) cannot be input chemical."))
            islossscheme(trans) ? isotopomerize(trans, collect(loss_elements(el[i - 1], el[i]))) : isotopomerize(trans, el[i])
        end
        )
    end
end

"""
    maximal_elements(elements::Dict) -> Dict

Dictionary of elements of maximal abundance.
"""
function maximal_elements(elements)
    max_dictionary = copy(elements)
    for (e, m) in elements
        x = minor_isotope(e)
        n = floor(Int, m * (1 - elements_abundance()[e]) / (elements_abundance()[e] + elements_abundance()[x]))
        if (max_dictionary[e] - n) / (n + 1) * elements_abundance()[x] / elements_abundance()[e] > 1
            n += 1
        end
        max_dictionary[e] -= n 
        max_dictionary[minor_isotope(e)] = n
    end
    max_dictionary
end

"""
    initial_proportion(element_precursor::Dict, element_product::Dict; precise = false) -> AbstractFloat

Proportion of isotopologue `element_product` relative to all possible isotopologues fragmented from `element_precursor`. 
"""
function initial_proportion(element_precursor::Dict, element_product::Dict; precise = false)
    any(<(0), values(element_product)) && return 0
    p = 1
    for (e, n) in element_precursor
        n - get(element_product, e, 0) < 0 && return 0
    end
    for (e, n) in element_precursor
        iselement(e) || continue
        np = get(element_product, e, 0)
        nr = n - np 
        vp = [get(element_product, i, 0) for i in elements_isotopes()[e]]
        f = precise ? multinomial(big(np), vp...) : try 
            multinomial(np, vp...)
        catch
            multinomial(big(np), vp...)
        end
        p *= f
        v = [get(element_precursor, i, 0) for i in elements_isotopes()[e]]
        f = precise ? multinomial(big(n), v...) : try 
            multinomial(n, v...)
        catch
            multinomial(big(n), v...)
        end
        p /= f
        vr = v .- vp
        f = precise ? multinomial(big(nr), vr...) : try 
            multinomial(nr, vr...)
        catch
            multinomial(big(nr), vr...)
        end
        p *= f
    end
    precise ? float(p) : convert(float(Int), p)
end

"""
    maximal_proportion(element_precursor::Dict, element_product::Dict; precise = false) -> AbstractFloat

Proportion of maximal isotopologue of `element_product` relative to all possible isotopologues fragmented from `element_precursor`. 
"""
function maximal_proportion(element_precursor, element_product; precise = false)
    el = group(parent_element, keys(element_precursor))
    p = 1
    for (e, v) in pairs(el) 
        pn = get(element_product, e, 0)
        v = filter(!iselement, v)
        pre = [element_precursor[x] for x in v]
        en = get(element_precursor, e, 0) + sum(pre)
        rn = en - pn
        pro = [round(Int, element_precursor[x] * pn / en) for x in v]
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
            f = precise ? multinomial(big(pn) - sum(pro), pro...) * multinomial(big(rn) - sum(res), res...) / multinomial(big(en) - sum(pre), pre...) : try 
                multinomial(pn - sum(pro), pro...) * multinomial(rn - sum(res), res...) / multinomial(en - sum(pre), pre...)
            catch
                multinomial(big(pn) - sum(pro), pro...) * multinomial(big(rn) - sum(res), res...) / multinomial(big(en) - sum(pre), pre...)
            end
            p *= f
        end
    end
    precise ? float(p) : convert(float(Int), p)
end

"""
    maximal_proportion_elements(element_precursor::Dict, element_product::Dict; precise = false) 

Dictionary of elements and proportion of maximal isotopologue of `element_product` relative to all possible isotopologues fragmented from `element_precursor`. 
"""
function maximal_proportion_elements(element_precursor, element_product; precise = false)
    el = group(parent_element, keys(element_precursor))
    p = 1
    first_element_product_dictionary = copy(element_product)
    for (e, v) in pairs(el) 
        pn = get(element_product, e, 0)
        v = filter(!iselement, v)
        pre = [element_precursor[x] for x in v]
        en = get(element_precursor, e, 0) + sum(pre)
        rn = en - pn
        pro = [round(Int, element_precursor[x] * pn / en) for x in v]
        for (i, x) in enumerate(v) 
            get!(first_element_product_dictionary, x, 0)
            first_element_product_dictionary[x] += pro[i]
            get!(first_element_product_dictionary, e, 0)
            first_element_product_dictionary[e] -= pro[i]
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
            f = precise ? multinomial(big(pn) - sum(pro), pro...) * multinomial(big(rn) - sum(res), res...) / multinomial(big(en) - sum(pre), pre...) : try 
                multinomial(pn - sum(pro), pro...) * multinomial(rn - sum(res), res...) / multinomial(en - sum(pre), pre...)
            catch
                multinomial(big(pn) - sum(pro), pro...) * multinomial(big(rn) - sum(res), res...) / multinomial(big(en) - sum(pre), pre...)
            end
            p *= f
        end
    end
    precise ? float(p) : convert(float(Int), p),
    first_element_product_dictionary
end

"""
    unique_group_mz_ab(mztable, transitions, colmz, colab, gf = nothing)

Group `transitions` by `gf`, create maps from group to id, and extract m/z values and abundance from `matable`.
"""
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

"""
    gf_parent_isotope(isotope = "[13C]")

Create function for grouping isotopologues by parent chemical and isotopomer state defining by `isotope`.
"""
function gf_parent_isotope(isotope = "[13C]")
    isotope_unit = elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]]
    x -> [chemicalparent(m) => _isotopomerstate(isotopomersisotopes(m), isotope_unit) for m in x]
end

function _isotopomerstate(isotopes, isotope_unit)
    ds = 0
    for (e, n) in isotopes
        ds += (elements_mass()[e] - elements_mass()[elements_parents()[e]]) * n
    end
    round(Int, ds / isotope_unit)
end

"""
    swap_elements(product, precursor, swap)

Swap elements of `products` and elements of residuals of `precursor` and `product` for eleemnts in `swap`.
"""
function swap_elements(product, precursor, swap)
    isempty(swap) && return copy(product)
    product2 = copy(product)
    for (k, v) in precursor
        if parent_element(k) in swap
            product2[k] = v - get(product, k, 0)
        end
    end
    product2
end

"""
    get_isotope_vec(input_element::Vector{Pair{String, Int}})
    get_isotope_vec(input_element::Dict{String, Int})

Create isotopic replacement vectors.
"""
get_isotope_vec(input_element::Vector{Pair{String, Int}}) = filter(x -> !iselement(first(x)), input_element)
get_isotope_vec(input_element::Dict{String, Int}) = [k => v for (k, v) in input_element if !iselement(k)]

"""
    get_element_dictinonary(input_element)

Element dictionary without isotopes.
"""
function get_element_dictinonary(input_element)
    element_dictionary = Dict{String, Int}()
    for (e, n) in input_element
        if iselement(e)
            get!(element_dictionary, e, 0)
            element_dictionary[e] += n
        end
    end
    element_dictionary
end

"""
    get_element_dictinonary(input_element)

Element dictionary without isotopes, and mass of isotopes
"""
function get_element_dictinonary_fixmass(input_element)
    element_dictionary = Dict{String, Int}()
    msfix = 0.0
    for (e, n) in input_element
        if iselement(e)
            get!(element_dictionary, e, 0)
            element_dictionary[e] += n
        else
            msfix += elements_mass()[e] * n
        end
    end
    element_dictionary, msfix
end

"""
    element_isotope_pairs(element_dictionary; sort = true)

Elements and isotopes pairs.
"""
function element_isotope_pairs(element_dictionary; sort = true)
    element_isotope_pair = vcat(map(collect(keys(element_dictionary))) do e
        v = get(elements_isotopes(), e, nothing)
        isnothing(v) && return Pair{String, String}[]
        length(v) < 1 && return Pair{String, String}[]
        map(v[begin + 1:end]) do x
            (e, x)
        end
    end...)
    sort ? sort!(element_isotope_pair; by = pair_last_abundance, rev = true) : element_isotope_pair
end

pair_last_abundance(x) = elements_abundance()[last(x)]

const MIN_PROPORTION = 1e-6
min_propotion() = MIN_PROPORTION
# ==========================================================================================================================
# Mid level MS1
isotopologues_elements(x::AbstractString, abundance, abtype, threshold; normalize = true) = 
    isotopologues_elements(chemicalelements(x), abundance, abtype, threshold; normalize)
isotopologues_elements(input_element::Vector, abundance, abtype, threshold; normalize = true) = 
    isotopologues_elements(get_element_dictinonary_fixmass(input_element)..., abundance, abtype, threshold; normalize)
function isotopologues_elements(element_dictionary::Dict, msfix, abundance, abtype, threshold; normalize = true)
    abtype = abtyped(abtype)
    isempty(element_dictionary) && return (; Element = [get_isotope_vec(element_dictionary)], Mass = [mmi(element_dictionary)], Abundance = [abundance]) 
    element_isotope_pair = element_isotope_pairs(element_dictionary)
    first_proportion = isotopicabundance(element_dictionary)
    max_dictionary = maximal_elements(element_dictionary)
    max_proportion = isotopicabundance(max_dictionary)
    total, abundance_cutoff = abundance_threshold(abtype, abundance, threshold, first_proportion, max_proportion)
    if abtype == Input() && first_proportion / max_proportion < min_propotion() 
        throw(ArgumentError("Isotopic abundance of input chemical is too small; try use `abtype` other than `Input()` or larger threshold`"))
    end
    if max_proportion * abundance_cutoff / total / first_proportion ^ 2 > 1 
        # abundance_cutoff / total << first_proportion
        element_chemical = [max_dictionary]
        abundance_chemical = [total * max_proportion]
        rec_addminusisotopes!(element_chemical, abundance_chemical, max_dictionary, element_isotope_pair, 1, first(abundance_chemical), abundance_cutoff, true, true)
    else
        element_chemical = [element_dictionary]
        abundance_chemical = [total * first_proportion]
        rec_addisotopes!(element_chemical, abundance_chemical, element_dictionary, element_isotope_pair, 1, first(abundance_chemical), abundance_cutoff)
    end
    abundance_chemical_normalize = normalize_abundance(abundance_chemical, abundance, preabtype(abtype), [Max(), Input(), Total()])
    mass_chemical = map(mmi, element_chemical) .+ msfix
    id = sortperm(mass_chemical)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(abundance_chemical_normalize)))
    filter!(x -> >=(abundance_chemical_normalize[x], abundance_cutoff), id)
    element_chemical = map(get_isotope_vec, element_chemical[id]) 
    if normalize && dopostnormalize(abtype)
        abundance_chemical = normalize_abundance(abundance_chemical_normalize[id], abundance, abtype, [Max(), Input(), List(), Total()])
    elseif normalize 
        abundance_chemical = abundance_chemical_normalize[id]
    else
        abundance_chemical = abundance_chemical[id]
    end
    mass_chemical = mass_chemical[id]
    (; Element = element_chemical, Mass = mass_chemical, Abundance = abundance_chemical) 
end

# ==========================================================================================================================
# Mid level MS2
function isotopologues_elements_ms2(it1, element_precursor_dictionary, element_product, abundance, abtype, proportion, threshold; gain = false, loss = false, normalize = true)
    if isempty(element_product)
        throw(ArgumentError("Product chemical must contain at least one element."))
    elseif gain 
        msfix = mmi(element_product) - mmi(element_precursor_dictionary)
        proportion_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * isotopicabundance(element_precursor_dictionary) / abundance
        it2 = isotopologues_elements(element_precursor_dictionary, 0, 1, Total(), proportion_cutoff)
        data = map(it1.Abundance) do r
            (; Element = it2.Element, 
            Mass = it2.Mass, 
            Abundance = it2.Abundance .* (r * proportion))
        end
    else
        element_product_dictionary, msfix = get_element_dictinonary_fixmass(element_product)
        i = findfirst(x -> all(iselement, keys(x)), it1.Element)
        p = 1
        if isnothing(i) 
            _, i = findmax(it1.Abundance)
            p = maximal_proportion(it1.Element[i], element_product_dictionary)
        end
        proportion_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * it1.Abundance[i] * p / abundance * proportion
        data = map(zip(it1.Element, it1.Abundance)) do (element, abundance)
            element_product_is, proportion_product_is = isotopes_proportion(element, element_product_dictionary, element_precursor_dictionary, abundance * proportion, proportion_cutoff)
            mass_product_is = map(mmi, element_product_is) 
            id = sortperm(mass_product_is)
            (; Element = loss ? [filter(!iselement ∘ first, loss_elements(element, x)) for x in element_product_is[id]] : element_product_is[id], 
            Mass = mass_product_is[id], 
            Abundance = proportion_product_is[id])
        end
    end
    element_product = vcat((getproperty(x, :Element) for x in data)...)
    mass_product = vcat((getproperty(x, :Mass) for x in data)...)
    abundance_pair = vcat((getproperty(x, :Abundance) for x in data)...)
    id_pair = vcat(([i for _ in eachindex(getproperty(x, :Abundance))] for (i, x) in enumerate(data))...)
    abtype = abtyped(abtype)
    abundance_pair_normalize = normalize_abundance(abundance_pair, abundance, abtype, [Total()])
    # spectrum specific threshold ?
    abundance_cutoff = isempty(abundance_pair_normalize) ? 0 : minimum(makecrit_value(crit(threshold), maximum(abundance_pair_normalize)))
    id = findall(>=(abundance_cutoff), abundance_pair_normalize)
    id_pair = id_pair[id]
    element_product = map(get_isotope_vec, element_product[id])
    mass_product = mass_product[id]
    mass_product .+= msfix
    abundance_pair = normalize ? abundance_pair_normalize[id] : abundance_pair[id]
    (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = abundance_pair) 
end

function isotopes_proportion(precursor_dictionary::Dict, element_product_dictionary::Dict, element_precursor_dictionary::Dict, abundance_factor, threshold)
    first_element_product_dictionary = copy(element_product_dictionary)
    swap = String[]
    for (e, n) in first_element_product_dictionary
        m = element_precursor_dictionary[e] - n
        if m < n
            push!(swap, e)
            first_element_product_dictionary[e] = m
        end
    end
    element_isotope_pair = element_isotope_pairs(element_product_dictionary; sort = false)
    first_proportion = initial_proportion(precursor_dictionary, first_element_product_dictionary)
    bidir = false
    if first_proportion < min_propotion()
        bidir = true
        first_proportion, first_element_product_dictionary = maximal_proportion_elements(precursor_dictionary, first_element_product_dictionary)
    end
    first_element_product_dictionary2 = swap_elements(first_element_product_dictionary, precursor_dictionary, swap)
    element_vec = [first_element_product_dictionary2]
    abundance_vec = [first_proportion * abundance_factor]
    bidir ? rec_exchangeisotopes!(element_vec, abundance_vec, first_element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(abundance_vec), threshold, swap, true, true) : 
        rec_moveisotopes!(element_vec, abundance_vec, first_element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(abundance_vec), threshold, swap)
    element_vec, abundance_vec
end

# ==========================================================================================================================
# Mid level MSn
function isotopologues_elements_msn(element_dictionary_vec, msfix_vec, abundance, abtype, threshold) 
    max_dictionary_vec = map(maximal_elements, element_dictionary_vec)
    max_proportion_vec = map(isotopicabundance, max_dictionary_vec)
    total, abundance_cutoff = abundance_threshold_msn(abtype, abundance, threshold, max_proportion_vec, element_dictionary_vec) 
    tbls = map(isotopologues_elements_msx, element_dictionary_vec, msfix_vec, max_dictionary_vec, max_proportion_vec, [total for _ in eachindex(element_dictionary_vec)], [abundance_cutoff for _ in eachindex(element_dictionary_vec)])
    first(tbls).Abundance .*= total
    els = Vector{Dict}[]
    abv = float(Int)[]
    mass = Vector{float(Int)}[]
    rec_vec_ab!(els, abv, mass, tbls, abundance_cutoff, prod(first(tbl.Abundance) for tbl in tbls), [1 for _ in eachindex(tbls)], 1)
    idm = sortperm(mass)
    mass = mass[idm]
    abv = abv[idm]
    if dopostnormalize(abtype)
        abv = normalize_abundance(abv, abundance, abtype)
    end
    (; Element = els, Mass = mass, Abundance = abv)
end

function isotopologues_elements_msx(element_dictionary, msfix, max_dictionary, max_proportion, abundance_factor, proportioon_cutoff)
    isempty(element_dictionary) && return (; Element = [dictionary_elements(get_isotope_vec(element_dictionary))], Mass = [mmi(element_dictionary) + msfix], Abundance = [float(1)])
    element_isotope_pair = element_isotope_pairs(element_dictionary)
    element_chemical = [max_dictionary]
    abundance_chemical = [max_proportion * abundance_factor]
    rec_addminusisotopes!(element_chemical, abundance_chemical, max_dictionary, element_isotope_pair, 1, first(abundance_chemical), proportioon_cutoff, true, true)
    mass_chemical = map(mmi, element_chemical) .+ msfix
    aid = sortperm(abundance_chemical; rev = true)
    element_chemical = map(x -> filter!(!iselement ∘ first, x), element_chemical[aid]) 
    mass_chemical = mass_chemical[aid]
    abundance_chemical = [abundance_chemical[i] / abundance_factor for i in aid]
    (; Element = element_chemical, Mass = mass_chemical, Abundance = abundance_chemical) 
end

# ==========================================================================================================================
# Low level MS1
function rec_addisotopes!(
        element_vec::Vector, 
        abundance_vec::Vector, 
        element_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_abundance, 
        threshold
    )
    current_state = false
    rec_state = false
    next_state = true
    max_abundance = prev_abundance
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair)
        next_state, next_abundance = rec_addisotopes!(element_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position + 1, prev_abundance, threshold)
        max_abundance = max(max_abundance, next_abundance)
    end
    ne = get(element_dictionary, e, 0)
    if ne > 0
        abundance = update_abundance(prev_abundance, e, i, ne, get(element_dictionary, i, 0), 1)
        max_abundance = max(max_abundance, abundance)
        current_state = (abundance >= threshold)
        rec_state = current_state || rec_state || (abundance >= prev_abundance)
        if !rec_state && (prev_abundance >= abundance > 0)
            rec_state = max_abundance / prev_abundance * abundance >= threshold 
        end
        if rec_state
            new_element_dictionary = copy(element_dictionary)
            new_element_dictionary[e] -= 1
            get!(new_element_dictionary, i, 0)
            new_element_dictionary[i] += 1
            if current_state
                push!(element_vec, new_element_dictionary)
                push!(abundance_vec, abundance)
            end
            rec_state, rec_abundance = rec_addisotopes!(element_vec, abundance_vec, new_element_dictionary, element_isotope_pair, isotope_position, abundance, threshold)
            current_state = rec_state || current_state || next_state
            max_abundance = max(max_abundance, rec_abundance)
        end
    end
    current_state, max_abundance
end

function rec_addminusisotopes!(
        element_vec::Vector, 
        abundance_vec::Vector, 
        element_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_abundance, 
        threshold, 
        backward,
        forward
    )
    next_state = true
    max_abundance = prev_abundance
    backward_state = false
    rec_state = false
    backward_abundance = prev_abundance
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair)
        next_state, next_abundance = rec_addminusisotopes!(element_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position + 1, prev_abundance, threshold, true, true)
        max_abundance = max(max_abundance, next_abundance)
    end
    ne = get(element_dictionary, e, 0)
    ni = get(element_dictionary, i, 0)
    if backward && ni > 0 
        abundance = update_abundance(prev_abundance, i, e, ni, ne, 1)
        backward_abundance = max(max_abundance, abundance)
        backward_state = (abundance >= threshold)
        rec_state = backward_state || rec_state || (abundance >= prev_abundance)
        if !rec_state && (prev_abundance >= abundance > 0)
            rec_state = backward_abundance / prev_abundance * abundance >= threshold 
        end
        if rec_state
            new_element_dictionary = copy(element_dictionary)
            new_element_dictionary[i] -= 1
            get!(new_element_dictionary, e, 0)
            new_element_dictionary[e] += 1
            if backward_state
                push!(element_vec, new_element_dictionary)
                push!(abundance_vec, abundance)
            end
            rec_state, rec_abundance = rec_addminusisotopes!(element_vec, abundance_vec, new_element_dictionary, element_isotope_pair, isotope_position, abundance, threshold, true, false)
            backward_state = rec_state || backward_state 
            backward_abundance = max(backward_abundance, rec_abundance)
        end
    end
    forward_state = false
    rec_state = false
    forward_abundance = prev_abundance
    if forward && ne > 0 
        abundance = update_abundance(prev_abundance, e, i, ne, ni, 1)
        forward_abundance = max(max_abundance, abundance)
        forward_state = (abundance >= threshold)
        rec_state = forward_state || rec_state || (abundance >= prev_abundance)
        if !rec_state && (prev_abundance >= abundance > 0)
            rec_state = forward_abundance / prev_abundance * abundance >= threshold 
        end
        if rec_state
            new_element_dictionary = copy(element_dictionary)
            new_element_dictionary[e] -= 1
            get!(new_element_dictionary, i, 0)
            new_element_dictionary[i] += 1
            if forward_state
                push!(element_vec, new_element_dictionary)
                push!(abundance_vec, abundance)
            end
            rec_state, rec_abundance = rec_addminusisotopes!(element_vec, abundance_vec, new_element_dictionary, element_isotope_pair, isotope_position, abundance, threshold, false, true)
            forward_state = rec_state || forward_state 
            forward_abundance = max(forward_abundance, rec_abundance)
        end
    end
    forward_state || backward_state || next_state, max(forward_abundance, backward_abundance)
end

function update_abundance(prev_abundance, old_element, new_element, nold, nnew, delta)
    x = get(elements_abundance(), old_element, 1)
    y = get(elements_abundance(), new_element, 1)
    (x == 1 || y == 1) && return prev_abundance
    if delta == 1
        return prev_abundance * nold / (nnew + 1) * (y / x)
    end
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
        element_product_dictionary::Dict, 
        element_precursor_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_proportion,
        threshold,
        swap
    )
    current_state = false
    rec_state = false
    next_state = true
    max_proportion = prev_proportion
    (e, i) = element_isotope_pair[isotope_position]
    np = get(element_product_dictionary, e, 0)
    ni = get!(element_product_dictionary, i, 0)
    nr = get(element_precursor_dictionary, i, 0) - ni
    if isotope_position < lastindex(element_isotope_pair) 
        next_state, next_proportion = rec_moveisotopes!(element_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_proportion, threshold, swap)
        max_proportion = max(max_proportion, next_proportion)
    end
    if np > 0 && nr >= 0
        if prev_proportion == 0
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
            proportion = initial_proportion(element_precursor_dictionary, element_product_dictionary)
        else
            proportion = update_proportion(prev_proportion, np, ni, 1)
            proportion = update_proportion(proportion, nr, element_precursor_dictionary[e] - np, 1)
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
        end
        if proportion >= threshold  
            push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
            push!(abundance_vec, proportion)
            current_state = true
        end
        max_proportion = max(max_proportion, proportion)
        rec_state = current_state || rec_state || (proportion >= prev_proportion) 
        if !rec_state && (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
            ne, ni = element_isotope_pair[isotope_position + 1]
            # same parent replacement
            if ne == e && get(element_product_dictionary, ni, 0) > 0
                rec_state = max_proportion >= threshold 
            else
                # Other parent peaks
                rec_state = max_proportion / prev_proportion * proportion >= threshold 
            end
        end
        if rec_state
            rec_state, rec_proportion = rec_moveisotopes!(element_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, proportion, threshold, swap)
            current_state = rec_state || current_state || next_state
            max_proportion = max(max_proportion, rec_proportion)
        end
        element_product_dictionary[e] += 1
        element_product_dictionary[i] -= 1
    end
    current_state, max_proportion
end

function rec_exchangeisotopes!(
        element_vec::Vector, 
        abundance_vec::Vector, 
        element_product_dictionary::Dict, 
        element_precursor_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_proportion,
        threshold,
        swap,
        backward, 
        forward
    )
    next_state = true
    max_proportion = prev_proportion
    backward_state = false
    rec_state = false
    backward_proportion = prev_proportion
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair) 
        next_state, next_proportion = rec_exchangeisotopes!(element_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_proportion, threshold, swap, true, true)
        max_proportion = max(max_proportion, next_proportion)
    end
    pis = get!(element_product_dictionary, i, 0)
    pel = get!(element_product_dictionary, e, 0)
    rel = get(element_precursor_dictionary, e, 0) - pel
    ris = get(element_precursor_dictionary, i, 0) - pis
    if backward && pis > 0 && rel >= 0
        if prev_proportion == 0
            element_product_dictionary[i] -= 1
            element_product_dictionary[e] += 1
            proportion = initial_proportion(element_precursor_dictionary, element_product_dictionary)
        else
            proportion = update_proportion(prev_proportion, rel, ris, 1)
            proportion = update_proportion(proportion, pis, pel, 1)
            element_product_dictionary[i] -= 1
            element_product_dictionary[e] += 1
        end
        if proportion >= threshold 
            push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
            push!(abundance_vec, proportion)
            backward_state = true 
        end
        backward_proportion = max(max_proportion, proportion)
        rec_state = backward_state || rec_state || (proportion >= prev_proportion) 
        if !rec_state && (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
            ne, ni = element_isotope_pair[isotope_position + 1]
            # same parent replacement
            if ne == e && get(element_precursor_dictionary, ni, 0) - get(element_product_dictionary, ni, 0) > 0
                rec_state = max_proportion >= threshold 
            else
                # Other parent peaks
                rec_state = max_proportion / prev_proportion * proportion >= threshold 
            end
        end
        if rec_state
            rec_state, rec_proportion = rec_exchangeisotopes!(element_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, proportion, threshold, swap, true, false) 
            backward_state = rec_state || backward_state 
            backward_proportion = max(backward_proportion, rec_proportion)
        end
        element_product_dictionary[i] += 1
        element_product_dictionary[e] -= 1
    end
    forward_state = false
    rec_state = false
    forward_proportion = prev_proportion
    if forward && pel > 0 && ris >= 0
        if prev_proportion == 0
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
            proportion = initial_proportion(element_precursor_dictionary, element_product_dictionary)
        else
            proportion = update_proportion(prev_proportion, pel, pis, 1)
            proportion = update_proportion(proportion, ris, rel, 1)
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
        end
        if proportion >= threshold  
            push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
            push!(abundance_vec, proportion)
            forward_state = true
        end
        backward_proportion = max(max_proportion, proportion)
        rec_state = forward_state || rec_state || (proportion >= prev_proportion) 
        if !rec_state && (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
            ne, ni = element_isotope_pair[isotope_position + 1]
            # same parent replacement
            if ne == e && get(element_product_dictionary, ni, 0) > 0
                rec_state = max_proportion >= threshold 
            else
                # Other parent peaks
                rec_state = max_proportion / prev_proportion * proportion >= threshold 
            end
        end
        if rec_state
            rec_state, rec_proportion = rec_exchangeisotopes!(element_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, proportion, threshold, swap, false, true) 
            forward_state = rec_state || forward_state 
            backward_proportion = max(backward_proportion, rec_proportion)
        end
        element_product_dictionary[e] += 1
        element_product_dictionary[i] -= 1
    end
    forward_state || backward_state || next_state, max(forward_proportion, backward_proportion)
end

function update_proportion(prev_proportion, nold, nnew, delta; precise = false) 
    if delta == 1
        return precise ? prev_proportion * float(nold) / float(nnew + 1) : prev_proportion * (nold / (nnew + 1))
    end
    f = precise ? factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew) : try 
        factorial(nold, nold - delta) / factorial(nnew + delta, nnew)
    catch 
        factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)
    end
    precise ? prev_proportion * float(f) : prev_proportion * convert(float(Int), f)
end

# ==========================================================================================================================
# Low level MSn
function rec_vec_ab!(els, abv, mass, tbls, abundance_cutoff, maxab, id, msn)
    msn == lastindex(tbls) && return rec_vec_ab_end!(els, abv, mass, tbls, abundance_cutoff, maxab, id)
    maxab /= first(tbls[msn].Abundance)
    pass = false
    @inbounds for (i, a) in enumerate(tbls[msn].Abundance)
        if a <= 0 
            break
        end
        id[msn] = i
        next_pass = rec_vec_ab!(els, abv, mass, tbls, abundance_cutoff, maxab * a, id, msn + 1)
        pass = pass || next_pass
        next_pass ? continue : break
    end
    id[msn] = 1
    pass
end

function rec_vec_ab_end!(els, abv, mass, tbls, abundance_cutoff, maxab, id)
    maxab /= first(tbls[end].Abundance)
    pass = false
    @inbounds for (i, a) in enumerate(tbls[end].Abundance)
        ab = maxab * a
        ab < abundance_cutoff && break
        id[end] = i
        el = [copy(tbl.Element[i]) for (tbl, i) in zip(tbls, id)]
        prev = empty(el[end])
        for i in Iterators.reverse(eachindex(el))
            prev = gain_elements!(el[i], prev)
        end
        push!(els, el)
        push!(abv, ab)
        push!(mass, reverse(cumsum(tbl.Mass[i] for (tbl, i) in Iterators.reverse(zip(tbls, id)))))
        pass = true
    end
    id[end] = 1
    pass
end
