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

return_abundance(::Val{true}, x) = x 
return_abundance(::Val{false}, x) = convert(float(Int), x) 

"""
    isotopologue_proportion(precise::Val, element_precursor::Dict, element_product::Dict, p = 1.0) -> AbstractFloat

Proportion of isotopologue `element_product` relative to all possible isotopologues fragmented from `element_precursor`. 
"""
function isotopologue_proportion(precise::Val, element_precursor::Dict, element_product::Dict, p::T = 1.0) where T
    any(<(0), values(element_product)) && return zero(T)
    for (e, v) in pairs(group(parent_element, keys(element_precursor)))
        np = get(element_product, e, 0)
        np < 0 && return zero(T)
        n = get(element_precursor, e, 0)
        n - np < 0 && return zero(T)
        filter!(!iselement, v)
        va = [get(element_precursor, i, 0) for i in v]
        vp = [get(element_product, i, 0) for i in v]
        p = update_isotopologue_proportion(precise, p, n, va, np, vp)
    end
    return_abundance(precise, p)
end

update_isotopologue_proportion(::Val{true}, p, n, va, np, vp) = 
    p * multinomial(big(np), vp...) / multinomial(big(n), va...) * multinomial(big(n - np), (va .- vp)...)

function update_isotopologue_proportion(::Val{false}, p, n, va, np, vp) 
    if check_overflow_multinomial(n, va...) 
        p * multinomial(big(np), vp...) / multinomial(big(n), va...) * multinomial(big(n - np), (va .- vp)...)
    else
        p * multinomial(np, vp...) / multinomial(n, va...) * multinomial(n - np, (va .- vp)...)
    end
end

"""
    isotopologue_combination(precise::Val, element_precursor::Dict, element_product::Dict, p = 1.0) -> AbstractFloat

The number of combination of product isotopologue `element_product` and precursor `element_precursor`. 
"""
function isotopologue_combination(precise::Val, element_precursor::Dict, element_product::Dict, p::T = 1.0) where T
    any(<(0), values(element_product)) && return zero(T)
    for (e, v) in pairs(group(parent_element, keys(element_precursor)))
        np = get(element_product, e, 0)
        np < 0 && return zero(T)
        n = get(element_precursor, e, 0)
        n - np < 0 && return zero(T)
        filter!(!iselement, v)
        va = [get(element_precursor, i, 0) for i in v]
        vp = [get(element_product, i, 0) for i in v]
        p = update_isotopologue_combination(precise, p, n, va, np, vp)
    end
    return_abundance(precise, p)
end

update_isotopologue_combination(::Val{true}, p, n, va, np, vp) = 
    p * multinomial(big(np), vp...) * multinomial(big(n - np), (va .- vp)...)

function update_isotopologue_combination(::Val{false}, p, n, va, np, vp) 
    if (np * 2 > n ? check_overflow_multinomial(np, vp...) : check_overflow_multinomial(n - np, (va .- vp)...))
        p * multinomial(big(np), vp...) * multinomial(big(n - np), (va .- vp)...)
    else
        p * safe_multinomial(np, vp...) * safe_multinomial(n - np, (va .- vp)...)
    end
end

"""
    maximal_proportion(precise::Val, element_precursor::Dict, element_product::Dict, p = 1.0) -> AbstractFloat

Estimate maximal product isotopologue of `element_product`, and compute the proportion relative to all possible isotopologues fragmented from `element_precursor`. 
"""
function maximal_proportion(precise::Val, element_precursor, element_product, p = 1.0)
    for (e, v) in pairs(group(parent_element, keys(element_precursor))) 
        pn = get(element_product, e, 0)
        filter!(!iselement, v)
        pre = [get(element_precursor, x, 0) for x in v]
        spre = sum(pre)
        en = get(element_precursor, e, 0) + spre
        pro = isempty(pre) ? Int[] : [round(Int, x * pn / en) for x in pre]
        spro = sum(pro)
        p = update_maximal_proportion(precise, p, en, spre, pre, pn, spro, pro)
    end
    return_abundance(precise, p)
end

update_maximal_proportion(::Val{true}, p, en, spre, pre, pn, spro, pro) = 
    p * multinomial(big(pn - spro), pro...) / multinomial(big(en - spre), pre...) * multinomial(big(en - pn - spre + spro), (pre .- pro)...)

function update_maximal_proportion(::Val{false}, p, en, spre, pre, pn, spro, pro) 
    if check_overflow_multinomial(en - spre, pre...)
        p * multinomial(big(pn - spro), pro...) / multinomial(big(en - spre), pre...) * multinomial(big(en - pn - spre + spro), (pre .- pro)...)
    else
        p * multinomial(pn - spro, pro...) / multinomial(en - spre, pre...) * multinomial(en - pn - spre + spro, (pre .- pro)...)
    end
end

"""
    maximal_proportion_elements(precise::Val, element_precursor::Dict, element_product::Dict, p = 1.0)

Estimate maximal product isotopologue of `element_product`, and compute `Dictionary` of elements and proportion relative to all possible isotopologues fragmented from `element_precursor`. 
"""
function maximal_proportion_elements(precise::Val, element_precursor::Dict, element_product::Dict, p = 1.0)
    first_element_product_dictionary = copy(element_product)
    for (e, v) in pairs(group(parent_element, keys(element_precursor))) 
        pn = get(element_product, e, 0)
        filter!(!iselement, v)
        pre = [get(element_precursor, x, 0) for x in v]
        spre = sum(pre)
        en = get(element_precursor, e, 0) + spre
        pro = isempty(pre) ? Int[] : [round(Int, x * pn / en) for x in pre]
        spro = sum(pro)
        for (x, n) in zip(v, pro) 
            get!(first_element_product_dictionary, x, 0)
            first_element_product_dictionary[x] += n
            get!(first_element_product_dictionary, e, 0)
            first_element_product_dictionary[e] -= n
        end
        p = update_maximal_proportion(precise, p, en, spre, pre, pn, spro, pro)
    end
    return_abundance(precise, p), first_element_product_dictionary
end

"""
    maximal_combination(precise::Val, element_precursor::Dict, element_product::Dict, p = 1.0) -> AbstractFloat

Estimate maximal product isotopologue of `element_product`, and compute the number of combinations. 
"""
function maximal_combination(precise::Val, element_precursor::Dict, element_product::Dict, p = 1.0)
    for (e, v) in pairs(group(parent_element, keys(element_precursor))) 
        pn = get(element_product, e, 0)
        filter!(!iselement, v)
        pre = [get(element_precursor, x, 0) for x in v]
        spre = sum(pre)
        en = get(element_precursor, e, 0) + spre
        rn = en - pn
        pro = isempty(pre) ? Int[] : [round(Int, x * pn / en) for x in pre]
        spro = sum(pro)
        res = pre .- pro
        p = update_maximal_combination(precise, p, en, spre, pre, pn, spro, pro)
    end
    return_abundance(precise, p)
end

"""
    maximal_combination_elements(precise::Val, element_precursor::Dict, element_product::Dict, p = 1.0)

Estimate maximal product isotopologue of `element_product`, and compute `Dictionary` of elements and the number of combinations. 
"""
function maximal_combination_elements(precise::Val, element_precursor, element_product, p = 1.0)
    first_element_product_dictionary = copy(element_product)
    for (e, v) in pairs(group(parent_element, keys(element_precursor))) 
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
        p = update_maximal_combination(precise, p, en, spre, pre, pn, spro, pro)
    end
    return_abundance(precise, p), first_element_product_dictionary
end

update_maximal_combination(::Val{true}, p, en, spre, pre, pn, spro, pro) = 
    p * multinomial(big(pn - spro), pro...) * multinomial(big(en - pn - spre + spro), (pre .- pro)...)

function update_maximal_combination(::Val{false}, p, en, spre, pre, pn, spro, pro) 
    res = pre .- pro
    if (2 * pn > en ? check_overflow_multinomial(pn - spro, pro...) : check_overflow_multinomial(en - pn - spre + spro, res...))
        p * multinomial(big(pn - spro), pro...) * multinomial(big(en - pn - spre + spro), res...)
    else
        p * safe_multinomial(pn - spro, pro...) * safe_multinomial(en - pn - spre + spro, res...)
    end
end

"""
    isotopologue_inverse_combination(precise::Val, elements::Dict)
    isotopologue_inverse_combination(precise::Val, element_precursor::Dict, element_product::Dict, swap::Vector{String})

Compute inverse of combination.
"""
function isotopologue_inverse_combination(precise::Val, elements::Dict)
    # ignore_isotopes && return 1.0
    return_abundance(precise, mapfoldl(/, group(parent_element ∘ first, elements); init = 1.0) do v 
        safe_multinomial(precise, last.(v))
    end)
end

function isotopologue_inverse_combination(precise::Val, element_precursor::Dict, element_product::Dict, swap::Vector{String})
    return_abundance(precise, mapfoldl(/, pairs(group(parent_element ∘ first, element_precursor)); init = 1.0) do (e, v) 
        if e in swap 
            v = [n - get(element_product, i, 0) for (i, n) in v]
            safe_multinomial(precise, v)
        else
            safe_multinomial(precise, [get(element_product, i, 0) for (i, n) in v])
        end
    end)
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
    swap_elements(product::Dict, precursor::Dict, swap::Vector{String})

Swap elements of `products` and elements of residuals of `precursor` and `product` for eleemnts in `swap`.
"""
swap_elements(product::Dict, precursor::Dict, swap::Vector{String}) = isempty(swap) ? get_isotope_vec(product) : [k => nswap(k, v, swap, product) for (k, v) in precursor if !iselement(k)]

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
    msfix = elements_mass()[""]
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
    element_isotope_pair = mapreduce(vcat, collect(keys(element_dictionary))) do e
        v = get(elements_isotopes(), e, nothing)
        isnothing(v) && return Pair{String, String}[]
        length(v) < 1 && return Pair{String, String}[]
        map(v[begin + 1:end]) do x
            (e, x)
        end
    end
    sort ? sort!(element_isotope_pair; by = pair_last_abundance, rev = true) : element_isotope_pair
end

pair_last_abundance(x) = elements_abundance()[last(x)]

"""
    element_mass_delta(old_element, new_element)

Mass of `new_element` minus mass of `old_element`.
"""
element_mass_delta(old_element, new_element) = elements_mass()[new_element] - elements_mass()[old_element]

"""
    parent_mass_delta(element)

Mass of `element` minus mass of parent element of `element`.
"""
parent_mass_delta(e) = elements_mass()[e] - elements_mass()[parent_element(e)]

"""
    parent_mass_delta(elements::Vector{Pair{String, Int}})

Sum of mass of each element in `elements` minus mass of their parent elements.
"""
deltammi(x) = isempty(x) ? float(0) : sum(k -> last(k) * parent_mass_delta(first(k)), x)

"""
    update_abundance(precise::Val, prev_abundance, old_element, new_element, nold, nnew, n)

Update of `prev_abundance` after `n` `old_element` replaced by `new_element`. `nold` and `nnew` are number of elements before replacing.
"""
function update_abundance(precise::Val, prev_abundance::T, old_element, new_element, nold, nnew, delta) where T
    x = get(elements_abundance(), old_element, one(T))
    y = get(elements_abundance(), new_element, one(T))
    if (x == one(T) || y == one(T))
        prev_abundance
    elseif delta == 1 
        update_abundance1(precise, prev_abundance, x, y, nold, nnew)
    else
        update_abundancen(precise, prev_abundance, x, y, nold, nnew, delta)
    end
end

update_abundance1(::Val{true}, prev_abundance, x, y, nold, nnew) = prev_abundance * (big(nold) / (nnew + 1)) * (big(y) / x)
update_abundance1(::Val{false}, prev_abundance, x, y, nold, nnew) = prev_abundance * (nold / (nnew + 1)) * (y / x)
update_abundancen(::Val{true}, prev_abundance, x, y, nold, nnew, delta) = 
    prev_abundance * (factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)) * (big(y) / x) ^ delta
function update_abundancen(::Val{false}, prev_abundance, x, y, nold, nnew, delta) 
    if (nold > nnew + delta ? check_overflow_factorial(nold, nold - delta) : check_overflow_factorial(nnew + delta, nnew))
        prev_abundance *= factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)
    else
        prev_abundance *= safe_factorial(nold, nold - delta) / safe_factorial(nnew + delta, nnew) 
    end
    convert(float(Int), prev_abundance * (y / x) ^ delta)
end

"""
    update_proportion(precise::Val, prev_proportion, nold, nnew, n) 
    
Update of `prev_proportion` after `n` elements replacing. `nold` and `nnew` are number of elements before replacing.
"""
function update_proportion(precise::Val, prev_proportion, nold, nnew, delta) 
    if delta == 1
        update_proportion1(precise, prev_proportion, nold, nnew)
    else
        update_proportionn(precise, prev_proportion, nold, nnew, delta)
    end
end

update_proportion1(::Val{true}, prev_proportion, nold, nnew) = prev_proportion * (big(nold) / (nnew + 1)) 
update_proportion1(::Val{false}, prev_proportion, nold, nnew) = prev_proportion * (nold / (nnew + 1)) 
update_proportionn(::Val{true}, prev_proportion, nold, nnew, delta) = 
    prev_proportion * (factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)) 
function update_proportionn(::Val{false}, prev_proportion, nold, nnew, delta) 
    if (nold > nnew + delta ? check_overflow_factorial(nold, nold - delta) : check_overflow_factorial(nnew + delta, nnew))
        convert(float(Int), prev_proportion * (factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)))
    else
        convert(float(Int), prev_proportion * (safe_factorial(nold, nold - delta) / safe_factorial(nnew + delta, nnew)))
    end
end

"""
    update_inverse_proportion(precise::Val, prev_inverse_proportion, nold, nnew, delta)

Update of `prev_inverse_proportion` after `n` elements replacing. `nold` and `nnew` are number of elements before replacing.
"""
function update_inverse_proportion(precise::Val, prev_proportion, nold, nnew, delta)
    if delta == 1 
        update_inverse_proportion1(precise, prev_proportion, nold, nnew)
    else
        update_inverse_proportionn(precise, prev_proportion, nold, nnew, delta)
    end
end

update_inverse_proportion1(::Val{true}, prev_proportion, nold, nnew) = prev_proportion * (big(nnew + 1) / nold) 
update_inverse_proportion1(::Val{false}, prev_proportion, nold, nnew) = prev_proportion * ((nnew + 1) / nold) 
update_inverse_proportionn(::Val{true}, prev_proportion, nold, nnew, delta) = 
    prev_proportion * (factorial(big(nnew + delta), nnew) / factorial(big(nold), nold - delta) ) 
function update_inverse_proportionn(::Val{false}, prev_proportion, nold, nnew, delta) 
    if (nold > nnew + delta ? check_overflow_factorial(nold, nold - delta) : check_overflow_factorial(nnew + delta, nnew))
        convert(float(Int), prev_proportion * (factorial(big(nnew + delta), nnew) / factorial(big(nold), nold - delta)))
    else
        convert(float(Int), prev_proportion * (safe_factorial(nnew + delta, nnew) / safe_factorial(nold, nold - delta)))
    end
end