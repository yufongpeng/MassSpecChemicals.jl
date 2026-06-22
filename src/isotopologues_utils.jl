const MIN_PROPORTION = 1e-3
min_propotion() = MIN_PROPORTION

"""
    detectedchemicaldata(raw_product; precursor)

Detected chemical and dictionary of elements.
"""
function detectedchemicaldata(raw_product; precursor)
    prod = detectedchemical(raw_product; precursor)
    prod, unique_elements(Dict, chemicalelements(prod))
end

"""
    serieschemicaldata(raw_product; precursor)

Series detected chemical and dictionaries of elements.
"""
function serieschemicaldata(input_chemical)
    precursor = seriesanalyzedchemical(input_chemical)
    v = map(precursor) do p 
        elements_precursor = chemicalelements(p)
        p, unique_elements(Dict, elements_precursor)
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
function seriesisotopomerize(transitions::Vector{<: AbstractChemicalsSchema}, els::Vector{<: Vector{Vector{Pair{String, Int}}}})
    @inbounds map(els) do el
        ChemicalTransition(map(enumerate(transitions)) do (i, trans)
            (islossscheme(trans) || isgainscheme(trans)) && i == firstindex(transitions) && throw(ArgumentError("$(typeof(trans)) cannot be input chemical."))
            islossscheme(trans) ? isotopomerize(trans, loss_elements(el[i], el[i - 1])) : isotopomerize(trans, el[i])
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
        abe = elements_abundance()[e]
        abx = elements_abundance()[x]
        n = floor(Int, m * (1 - abe) / (abe + abx))
        if (max_dictionary[e] - n) / (n + 1) * abx / abe > 1
            n += 1
        end
        max_dictionary[e] -= n 
        max_dictionary[x] = n
    end
    max_dictionary
end

"""
    isotopologue_proportion(element_precursor::Dict, element_product::Dict; precise = false) -> AbstractFloat

Proportion of isotopologue `element_product` relative to all possible isotopologues fragmented from `element_precursor`. 
"""
function isotopologue_proportion(element_precursor::Dict, element_product::Dict; precise = false)
    any(<(0), values(element_product)) && return 0
    el = group(parent_element, keys(element_precursor))
    p = 1
    for (e, v) in pairs(el)
        np = get(element_product, e, 0)
        np < 0 && return 0
        n = get(element_precursor, e, 0)
        n - np < 0 && return 0
        filter!(!iselement, v)
        va = [get(element_precursor, i, 0) for i in v]
        vp = [get(element_product, i, 0) for i in v]
        if precise || check_overflow_multinomial(n, va...)
            p *= multinomial(big(np), vp...) / multinomial(big(n), va...) * multinomial(big(n - np), (va .- vp)...)
        else
            p *= multinomial(np, vp...) / multinomial(n, va...) * multinomial(n - np, (va .- vp)...)
        end
    end
    precise ? p : convert(float(Int), p)
end

function isotopologue_combination(element_precursor::Dict, element_product::Dict, p = 1; precise = false)
    any(<(0), values(element_product)) && return 0
    el = group(parent_element, keys(element_precursor))
    for (e, v) in pairs(el)
        np = get(element_product, e, 0)
        np < 0 && return 0
        n = get(element_precursor, e, 0)
        n - np < 0 && return 0
        filter!(!iselement, v)
        va = [get(element_precursor, i, 0) for i in v]
        vp = [get(element_product, i, 0) for i in v]
        if precise || (np * 2 > n ? check_overflow_multinomial(np, vp...) : check_overflow_multinomial(n - np, (va .- vp)...))
            p *= multinomial(big(np), vp...) 
            p *= multinomial(big(n - np), (va .- vp)...)
        else
            p *= safe_multinomial(np, vp...) 
            p *= safe_multinomial(n - np, (va .- vp)...)
        end
    end
    precise ? p : convert(float(Int), p)
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
        filter!(!iselement, v)
        pre = [get(element_precursor, x, 0) for x in v]
        spre = sum(pre)
        en = get(element_precursor, e, 0) + spre
        rn = en - pn
        pro = isempty(pre) ? Int[] : [round(Int, x * pn / en) for x in pre]
        spro = sum(pro)
        if abs(pn + spro * 2 - rn - spre) > en - spre
            p *= safe_multinomial(pn, rn) 
            for (i, j) in zip(pre, pro)
                p *= safe_multinomial(i - j, j) 
            end
        else
            res = pre .- pro
            if precise || check_overflow_multinomial(en - spre, pre...)
                p *= multinomial(big(pn - spro), pro...) / multinomial(big(en - spre), pre...) * multinomial(big(rn - spre + spro), res...)
            else
                p *= multinomial(pn - spro, pro...) / multinomial(en - spre, pre...) * multinomial(rn - spre + spro, res...)
            end
        end
    end
    precise ? p : convert(float(Int), p)
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
        filter!(!iselement, v)
        pre = [get(element_precursor, x, 0) for x in v]
        spre = sum(pre)
        en = get(element_precursor, e, 0) + spre
        rn = en - pn
        pro = isempty(pre) ? Int[] : [round(Int, x * pn / en) for x in pre]
        spro = sum(pro)
        for (x, n) in zip(v, pro) 
            get!(first_element_product_dictionary, x, 0)
            first_element_product_dictionary[x] += n
            get!(first_element_product_dictionary, e, 0)
            first_element_product_dictionary[e] -= n
        end
        if abs(pn + spro - rn - spre + spro) > en - spre
            p *= safe_multinomial(pn, rn) 
            for (i, j) in zip(pre, pro)
                p *= safe_multinomial(i - j, j) 
            end
        else
            res = pre .- pro
            if precise || check_overflow_multinomial(en - spre, pre...)
                p *= multinomial(big(pn - spro), pro...) / multinomial(big(en - spre), pre...) * multinomial(big(rn - spre + spro), res...)
            else
                p *= multinomial(pn - spro, pro...) / multinomial(en - spre, pre...) * multinomial(rn - spre + spro, res...)
            end
        end
    end
    precise ? p : convert(float(Int), p),
    first_element_product_dictionary
end

function maximal_combination(element_precursor, element_product, p = 1; precise = false)
    el = group(parent_element, keys(element_precursor))
    for (e, v) in pairs(el) 
        pn = get(element_product, e, 0)
        filter!(!iselement, v)
        pre = [get(element_precursor, x, 0) for x in v]
        spre = sum(pre)
        en = get(element_precursor, e, 0) + spre
        rn = en - pn
        pro = isempty(pre) ? Int[] : [round(Int, x * pn / en) for x in pre]
        spro = sum(pro)
        res = pre .- pro
        if precise || (pn > rn ? check_overflow_multinomial(pn - spro, pro...) : check_overflow_multinomial(rn - spre + spro, res...))
            p *= multinomial(big(pn - spro), pro...) 
            p *= multinomial(big(rn - spre + spro), res...)
        else
            p *= safe_multinomial(pn - spro, pro...) 
            p *= safe_multinomial(rn - spre + spro, res...)
        end
    end
    precise ? p : convert(float(Int), p)
end

function maximal_combination_elements(element_precursor, element_product, p = 1; precise = false)
    el = group(parent_element, keys(element_precursor))
    first_element_product_dictionary = copy(element_product)
    for (e, v) in pairs(el) 
        pn = get(element_product, e, 0)
        filter!(!iselement, v)
        pre = [get(element_precursor, x, 0) for x in v]
        spre = sum(pre)
        en = get(element_precursor, e, 0) + spre
        rn = en - pn
        pro = isempty(pre) ? Int[] : [round(Int, x * pn / en) for x in pre]
        spro = sum(pro)
        res = pre .- pro
        for (x, n) in zip(v, pro) 
            get!(first_element_product_dictionary, x, 0)
            first_element_product_dictionary[x] += n
            get!(first_element_product_dictionary, e, 0)
            first_element_product_dictionary[e] -= n
        end
        if precise || (pn > rn ? check_overflow_multinomial(pn - spro, pro...) : check_overflow_multinomial(rn - spre + spro, res...))
            p *= multinomial(big(pn - spro), pro...) 
            p *= multinomial(big(rn - spre + spro), res...)
        else
            p *= safe_multinomial(pn - spro, pro...) 
            p *= safe_multinomial(rn - spre + spro, res...)
        end
    end
    precise ? p : convert(float(Int), p), 
    first_element_product_dictionary
end

function isotopologue_inverse_combination(elements::Dict; ignore_isotopes = false, precise = false)
    ignore_isotopes && return 1
    element_dict = group(parent_element ∘ first, elements)
    abundance_sum = 1
    for v in element_dict
        abundance_sum /= safe_multinomial(precise, last.(v))
    end
    precise ? abundance_sum : convert(float(Int), abundance_sum)
end

function isotopologue_inverse_combination(element_precursor, element_product, swap; precise = false)
    element_dict = group(parent_element ∘ first, element_precursor)
    abundance_sum = 1
    for (e, v) in pairs(element_dict)
        if e in swap 
            v = [n - get(element_product, i, 0) for (i, n) in v]
            abundance_sum /= safe_multinomial(precise, v)
        else
            abundance_sum /= safe_multinomial(precise, [get(element_product, i, 0) for (i, n) in v])
        end
    end
    precise ? abundance_sum : convert(float(Int), abundance_sum)
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
            c = [y[begin:i] for y in transitions[v]]
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
    x -> [chemicalparent(m) => isotopomerstate(m; isotope_unit) for m in x]
end

"""
    swap_elements(product, precursor, swap)

Swap elements of `products` and elements of residuals of `precursor` and `product` for eleemnts in `swap`.
"""
swap_elements(product, precursor, swap) = isempty(swap) ? get_isotope_vec(product) : [k => nswap(k, v, swap, product) for (k, v) in precursor if !iselement(k)]

nswap(k, v, swap, product) = parent_element(k) in swap ? v - get(product, k, 0) : get(product, k, 0) 

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
get_element_dictinonary(input_element::Dict) = filter(iselement ∘ first, input_element)

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

element_mass_delta(old_element, new_element) = elements_mass()[new_element] - elements_mass()[old_element]
parent_mass_delta(e) = elements_mass()[e] - elements_mass()[parent_element(e)]
deltammi(x) = isempty(x) ? float(0) : sum(k -> last(k) * parent_mass_delta(first(k)), x)

function update_abundance(prev_abundance, old_element, new_element, nold, nnew, delta; precise = false)
    x = get(elements_abundance(), old_element, 1)
    y = get(elements_abundance(), new_element, 1)
    (x == 1 || y == 1) && return prev_abundance
    delta == 1 && return precise ? big(prev_abundance) * nold / (nnew + 1) * (big(y) / x) : prev_abundance * nold / (nnew + 1) * (y / x)
    # take care of large numbers
    if precise || (nold > nnew + delta ? check_overflow_factorial(nold, nold - delta) : check_overflow_factorial(nnew + delta, nnew))
        prev_abundance *= factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)
    else
        prev_abundance *= safe_factorial(nold, nold - delta) / safe_factorial(nnew + delta, nnew) 
    end
    precise ? prev_abundance * (big(y) / x) ^ delta : convert(float(Int), prev_abundance * (y / x) ^ delta)
end

function update_proportion(prev_proportion, nold, nnew, delta; precise = false) 
    delta == 1 && return precise ?  prev_proportion * (big(nold) / (nnew + 1)) : prev_proportion * (nold / (nnew + 1))
    f = precise ? factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew) : try 
        factorial(nold, nold - delta) / factorial(nnew + delta, nnew)
    catch 
        factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)
    end
    precise ? prev_proportion * float(f) : prev_proportion * convert(float(Int), f)
end


function update_inverse_proportion(prev_proportion, nold, nnew, delta; precise = false)
    delta == 1 && return precise ? prev_proportion * (big(nnew + 1) / nold) : prev_proportion * ((nnew + 1) / nold)
    # take care of large numbers
    if precise || (nold > nnew + delta ? check_overflow_factorial(nold, nold - delta) : check_overflow_factorial(nnew + delta, nnew))
        prev_proportion *= factorial(big(nnew + delta), nnew) / factorial(big(nold), nold - delta) 
    else
        prev_proportion *= safe_factorial(nnew + delta, nnew) / safe_factorial(nold, nold - delta) 
    end
    precise ? prev_proportion : convert(float(Int), prev_proportion)
end