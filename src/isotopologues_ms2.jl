# ==========================================================================================================================
# Mid level MS2
function isotopologues_elements_ms2(it1, element_precursor_dictionary, element_product, abundance, abtype, proportion, threshold; iter = true, gain = false, loss = false, normalize = true, precise = false)
    if isempty(element_product)
        throw(ArgumentError("Product chemical must contain at least one element."))
    elseif gain 
        msfix = mmi(element_product) - mmi(element_precursor_dictionary)
        proportion_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * isotopicabundance(element_precursor_dictionary) / abundance
        it2 = isotopologues_elements(dictionary_elements(element_precursor_dictionary), msfix, 1, Total(), proportion_cutoff)
        data = map(it1.Abundance) do r
            (; Element = it2.Element, 
            Mass = it2.Mass, 
            Abundance = it2.Abundance .* (r * proportion))
        end
    elseif hasproperty(it1, :Preab)
        element_product_dictionary = get_element_dictinonary(element_product)
        i = findfirst(x -> all(iselement, keys(x)), it1.Element)
        use_max = false
        if isnothing(i) 
            _, i = findmax(it1.Abundance)
            p = maximal_combination(it1.Element[i], element_product_dictionary, it1.Abundance[i] * it1.Preab[i] * proportion; precise)
        else
            p = it1.Abundance[i] * it1.Preab[i] * proportion
        end
        use_max = isotopologue_combination(first(it1.Element), element_product_dictionary, first(it1.Preab); precise) < min_propotion()
        msfix = mmi(element_product)
        element_isotope_pair = element_isotope_pairs(element_product_dictionary; sort = false)
        swap = String[]
        for (e, n) in element_product_dictionary
            m = element_precursor_dictionary[e] - n
            if m < n
                push!(swap, e)
                element_product_dictionary[e] = m
            end
        end
        proportion_cutoff = p * minimum(makecrit_value(crit(threshold), abundance)) / abundance
        @inbounds data = map(eachindex(it1.Element)) do idp
            element_product_is, mass_product_is, proportion_product_is, preab_product_is, use_max = isotopologues_proportion(it1.Element[idp], element_product_dictionary, swap, element_isotope_pair, msfix, (it1.Abundance[idp] * proportion, it1.Preab[idp]), proportion_cutoff, use_max, true, iter, precise)
            id = sortperm(mass_product_is)
            if loss 
                isotopes_precursor = it1.Isotope[idp]
                element_product_is = [loss_elements(element_product_is[i], isotopes_precursor) for i in id] 
            else
                element_product_is = element_product_is[id]
            end
            if iter 
                (; Element = element_product_is, 
                Mass = mass_product_is[id], 
                Abundance = proportion_product_is[id],
                Preab = preab_product_is[id])
            else
                (; Element = element_product_is, 
                Mass = mass_product_is[id], 
                Abundance = proportion_product_is[id]
                )
            end
        end
    else
        element_product_dictionary = get_element_dictinonary(element_product)
        i = findfirst(x -> all(iselement, keys(x)), it1.Element)
        p = 1
        use_max = false
        if isnothing(i) 
            _, i = findmax(it1.Abundance)
            p = maximal_proportion(it1.Element[i], element_product_dictionary; precise)
        end
        use_max = isotopologue_proportion(first(it1.Element), element_product_dictionary; precise) < min_propotion()
        msfix = mmi(element_product)
        element_isotope_pair = element_isotope_pairs(element_product_dictionary; sort = false)
        swap = String[]
        for (e, n) in element_product_dictionary
            m = element_precursor_dictionary[e] - n
            if m < n
                push!(swap, e)
                element_product_dictionary[e] = m
            end
        end
        proportion_cutoff = p * it1.Abundance[i] * minimum(makecrit_value(crit(threshold), abundance)) / abundance * proportion
        @inbounds data = map(eachindex(it1.Element)) do idp
            element_product_is, mass_product_is, proportion_product_is, preab_product_is, use_max = isotopologues_proportion(it1.Element[idp], element_product_dictionary, swap, element_isotope_pair, msfix, it1.Abundance[idp] * proportion, proportion_cutoff, use_max, false, iter, precise)
            id = sortperm(mass_product_is)
            if loss 
                isotopes_precursor = it1.Isotope[idp]
                element_product_is = [loss_elements(element_product_is[i], isotopes_precursor) for i in id] 
            else
                element_product_is = element_product_is[id]
            end
            if iter 
                (; Element = element_product_is, 
                Mass = mass_product_is[id], 
                Abundance = proportion_product_is[id],
                Preab = preab_product_is[id])
            else
                (; Element = element_product_is, 
                Mass = mass_product_is[id], 
                Abundance = proportion_product_is[id]
                )
            end
        end
    end
    v = getproperty.(data, :Element)
    element_product = ChainedVector(v)
    mass_product = ChainedVector(getproperty.(data, :Mass))
    abundance_pair = ChainedVector(getproperty.(data, :Abundance))
    id_pair = ChainedVector([[i for _ in eachindex(x)] for (i, x) in enumerate(v)])
    abtype = abtyped(abtype)
    abundance_pair_normalize = normalize_abundance(abundance_pair, abundance, abtype, [Total()])
    # spectrum specific threshold ?
    if isempty(abundance_pair_normalize) 
        if iter
            (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = normalize ? abundance_pair_normalize : abundance_pair, Preab = ChainedVector(getproperty.(data, :Preab)))
        else
            (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = normalize ? abundance_pair_normalize : abundance_pair)
        end
    else
        abundance_cutoff = isempty(abundance_pair_normalize) ? 0 : minimum(makecrit_value(crit(threshold), maximum(abundance_pair_normalize)))
        id = abundance_pair_normalize .>= abundance_cutoff
        if all(id) 
            if iter 
                (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = normalize ? abundance_pair_normalize : abundance_pair, Preab = ChainedVector(getproperty.(data, :Preab)))
            else
                (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = normalize ? abundance_pair_normalize : abundance_pair) 
            end
        elseif iter
            (; ID = id_pair[id], Element = element_product[id], Mass = mass_product[id], Abundance = normalize ? abundance_pair_normalize[id] : abundance_pair[id], Preab = ChainedVector(getproperty.(data, :Preab))[id])
        else
            (; ID = id_pair[id], Element = element_product[id], Mass = mass_product[id], Abundance = normalize ? abundance_pair_normalize[id] : abundance_pair[id]) 
        end
    end
end

function isotopologues_proportion(precursor_dictionary::Dict, element_product_dictionary::Dict, swap, element_isotope_pair, msfix, abundance_factor, threshold, use_max, preiter, iter, precise)
    if preiter
        if !use_max
            first_proportion = isotopologue_combination(precursor_dictionary, element_product_dictionary, *(abundance_factor...); precise)
        end
        if use_max || first_proportion / first(abundance_factor) < min_propotion()
            use_max = true
            first_proportion, element_product_dictionary = maximal_combination_elements(precursor_dictionary, element_product_dictionary, *(abundance_factor...); precise)
        end
    else
        if !use_max
            first_proportion = isotopologue_proportion(precursor_dictionary, element_product_dictionary; precise)
        end
        if use_max || first_proportion < min_propotion()
            use_max = true
            first_proportion, element_product_dictionary = maximal_proportion_elements(precursor_dictionary, element_product_dictionary; precise)
        end
        first_proportion *= abundance_factor
    end
    element_vec = [swap_elements(element_product_dictionary, precursor_dictionary, swap)]
    mass_vec = [deltammi(first(element_vec)) + msfix]
    abundance_vec = [first_proportion]
    if iter 
        preab_vec = [isotopologue_inverse_combination(precursor_dictionary, element_product_dictionary, swap; precise)]
        use_max ? rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, swap, true, true, precise) : 
            rec_moveisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, swap, precise)
        element_vec, mass_vec, abundance_vec, preab_vec, use_max
    else
        use_max ? rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), threshold, swap, true, true, precise) : 
            rec_moveisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), threshold, swap, precise)
        element_vec, mass_vec, abundance_vec, nothing, use_max
    end
end

# ==========================================================================================================================
# Low level MS2
function rec_moveisotopes!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        element_product_dictionary::Dict, 
        element_precursor_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_proportion,
        threshold,
        swap, 
        precise
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
        next_state, next_proportion = rec_moveisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, threshold, swap, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    if np > 0 && nr >= 0
        if prev_proportion == 0
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
            proportion = isotopologue_proportion(element_precursor_dictionary, element_product_dictionary; precise)
        else
            proportion = update_proportion(prev_proportion, np, ni, 1; precise)
            proportion = update_proportion(proportion, nr, element_precursor_dictionary[e] - np, 1; precise)
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
        end
        new_mass = prev_mass + (e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i))
        if proportion >= threshold 
            push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
            push!(abundance_vec, proportion)
            push!(mass_vec, new_mass)
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
            rec_state, rec_proportion = rec_moveisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, threshold, swap, precise)
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
        mass_vec::Vector,
        abundance_vec::Vector, 
        element_product_dictionary::Dict, 
        element_precursor_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_proportion,
        threshold,
        swap,
        backward, 
        forward, 
        precise
    )
    next_state = true
    max_proportion = prev_proportion
    backward_state = false
    rec_state = false
    backward_proportion = prev_proportion
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair) 
        next_state, next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, threshold, swap, true, true, precise)
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
            proportion = isotopologue_proportion(element_precursor_dictionary, element_product_dictionary; precise)
        else
            proportion = update_proportion(prev_proportion, rel, ris, 1; precise)
            proportion = update_proportion(proportion, pis, pel, 1; precise)
            element_product_dictionary[i] -= 1
            element_product_dictionary[e] += 1
        end
        new_mass = prev_mass + (e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e))
        if proportion >= threshold 
            push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
            push!(abundance_vec, proportion)
            push!(mass_vec, new_mass)
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
            rec_state, rec_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, threshold, swap, true, false, precise) 
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
            proportion = isotopologue_proportion(element_precursor_dictionary, element_product_dictionary; precise)
        else
            proportion = update_proportion(prev_proportion, pel, pis, 1; precise)
            proportion = update_proportion(proportion, ris, rel, 1; precise)
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
        end
        new_mass = prev_mass + (e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i))
        if proportion >= threshold  
            push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
            push!(abundance_vec, proportion)
            push!(mass_vec, new_mass)
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
            rec_state, rec_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, threshold, swap, false, true, precise) 
            forward_state = rec_state || forward_state 
            backward_proportion = max(backward_proportion, rec_proportion)
        end
        element_product_dictionary[e] += 1
        element_product_dictionary[i] -= 1
    end
    forward_state || backward_state || next_state, max(forward_proportion, backward_proportion)
end

function rec_moveisotopes_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_product_dictionary::Dict, 
        element_precursor_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_proportion,
        prev_preab,
        threshold,
        swap, 
        precise
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
        next_state, next_proportion = rec_moveisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, prev_preab, threshold, swap, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    if np > 0 && nr >= 0
        if prev_proportion == 0
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
            proportion = isotopologue_proportion(element_precursor_dictionary, element_product_dictionary; precise)
            preab = isotopologue_inverse_combination(element_precursor_dictionary, element_product_dictionary, swap; precise)
        else
            proportion = update_proportion(prev_proportion, np, ni, 1; precise)
            proportion = update_proportion(proportion, nr, element_precursor_dictionary[e] - np, 1; precise)
            preab = e in swap ? update_inverse_proportion(prev_preab, nr, element_precursor_dictionary[e] - np, 1; precise) : update_inverse_proportion(prev_preab, np, ni, 1; precise)
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
        end
        new_mass = prev_mass + (e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i))
        if proportion >= threshold 
            push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
            push!(abundance_vec, proportion)
            push!(preab_vec, preab)
            push!(mass_vec, new_mass)
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
            rec_state, rec_proportion = rec_moveisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, precise)
            current_state = rec_state || current_state || next_state
            max_proportion = max(max_proportion, rec_proportion)
        end
        element_product_dictionary[e] += 1
        element_product_dictionary[i] -= 1
    end
    current_state, max_proportion
end

function rec_exchangeisotopes_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_product_dictionary::Dict, 
        element_precursor_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_proportion,
        prev_preab,
        threshold,
        swap,
        backward, 
        forward, 
        precise
    )
    next_state = true
    max_proportion = prev_proportion
    backward_state = false
    rec_state = false
    backward_proportion = prev_proportion
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair) 
        next_state, next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, prev_preab, threshold, swap, true, true, precise)
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
            proportion = isotopologue_proportion(element_precursor_dictionary, element_product_dictionary; precise)
            preab = isotopologue_inverse_combination(element_precursor_dictionary, element_product_dictionary, swap; precise)
        else
            proportion = update_proportion(prev_proportion, rel, ris, 1; precise)
            proportion = update_proportion(proportion, pis, pel, 1; precise)
            preab = e in swap ? update_inverse_proportion(prev_preab, rel, ris, 1; precise) : update_inverse_proportion(prev_preab, pis, pel, 1; precise)
            element_product_dictionary[i] -= 1
            element_product_dictionary[e] += 1
        end
        new_mass = prev_mass + (e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e))
        if proportion >= threshold 
            push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
            push!(abundance_vec, proportion)
            push!(preab_vec, preab)
            push!(mass_vec, new_mass)
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
            rec_state, rec_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, true, false, precise) 
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
            proportion = isotopologue_proportion(element_precursor_dictionary, element_product_dictionary; precise)
            preab = isotopologue_inverse_combination(element_precursor_dictionary, element_product_dictionary, swap; precise)
        else
            proportion = update_proportion(prev_proportion, pel, pis, 1; precise)
            proportion = update_proportion(proportion, ris, rel, 1; precise)
            preab = e in swap ? update_inverse_proportion(prev_preab, ris, rel, 1) : update_inverse_proportion(prev_preab, pel, pis, 1; precise)
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
        end
        new_mass = prev_mass + (e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i))
        if proportion >= threshold  
            push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
            push!(abundance_vec, proportion)
            push!(preab_vec, preab)
            push!(mass_vec, new_mass)
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
            rec_state, rec_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, false, true, precise) 
            forward_state = rec_state || forward_state 
            backward_proportion = max(backward_proportion, rec_proportion)
        end
        element_product_dictionary[e] += 1
        element_product_dictionary[i] -= 1
    end
    forward_state || backward_state || next_state, max(forward_proportion, backward_proportion)
end