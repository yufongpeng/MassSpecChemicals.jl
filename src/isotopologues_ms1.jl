# ==========================================================================================================================
# Mid level MS1
isotopologues_elements(x::AbstractString, abundance, abtype, threshold; normalize = true, precise = false) = 
    isotopologues_elements(chemicalelements(x), abundance, abtype, threshold; normalize, precise)
isotopologues_elements(input_element::Vector, abundance, abtype, threshold; normalize = true, precise = false) = 
    isotopologues_elements(get_element_dictinonary_fixmass(input_element)..., abundance, abtype, threshold; normalize, precise)
function isotopologues_elements(element_dictionary::Dict, msfix, abundance, abtype, threshold; normalize = true, precise = false)
    abtype = abtyped(abtype)
    isempty(element_dictionary) && return (; Element = [element_dictionary], Mass = [mmi(element_dictionary)], Abundance = [abundance]) 
    element_isotope_pair = element_isotope_pairs(element_dictionary)
    first_proportion = isotopicabundance(element_dictionary; precise)
    max_dictionary = maximal_elements(element_dictionary)
    max_proportion = isotopicabundance(max_dictionary; precise)
    total, abundance_cutoff = abundance_threshold(abtype, abundance, threshold, first_proportion, max_proportion)
    if abtype == Input() && first_proportion / max_proportion < min_propotion() 
        throw(ArgumentError("Isotopic abundance of input chemical is too small; try use `abtype` other than `Input()` or larger threshold`"))
    end
    if max_proportion * abundance_cutoff / total / first_proportion ^ 2 > 1 
        # abundance_cutoff / total << first_proportion
        element_chemical = [get_isotope_vec(max_dictionary)]
        mass_chemical = [mmi(max_dictionary) + msfix]
        abundance_chemical = [total * max_proportion]
        rec_addminusisotopes!(element_chemical, mass_chemical, abundance_chemical, max_dictionary, element_isotope_pair, 1, first(mass_chemical), first(abundance_chemical), abundance_cutoff, true, true, precise)
    else
        element_chemical = [get_isotope_vec(element_dictionary)]
        mass_chemical = [mmi(element_dictionary) + msfix]
        abundance_chemical = [total * first_proportion]
        rec_addisotopes!(element_chemical, mass_chemical, abundance_chemical, element_dictionary, element_isotope_pair, 1, first(mass_chemical), first(abundance_chemical), abundance_cutoff, precise)
    end
    abundance_chemical_normalize = normalize_abundance(abundance_chemical, abundance, preabtype(abtype), [Max(), Input(), Total()])
    id = sortperm(mass_chemical)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(abundance_chemical_normalize)))
    filter!(x -> >=(abundance_chemical_normalize[x], abundance_cutoff), id)
    if normalize && dopostnormalize(abtype)
        abundance_chemical = normalize_abundance(abundance_chemical_normalize[id], abundance, abtype, [Max(), Input(), List(), Total()])
    elseif normalize 
        abundance_chemical = abundance_chemical_normalize[id]
    else
        abundance_chemical = abundance_chemical[id]
    end
    (; Element = element_chemical[id], Mass = mass_chemical[id], Abundance = abundance_chemical) 
end

isotopologues_elements_iter(x::AbstractString, abundance, abtype, threshold; normalize = true, precise = false) = 
    isotopologues_elements_iter(chemicalelements(x), abundance, abtype, threshold; normalize, precise)
isotopologues_elements_iter(input_element::Vector, abundance, abtype, threshold; normalize = true, precise = false) = 
    isotopologues_elements_iter(get_element_dictinonary_fixmass(input_element)..., abundance, abtype, threshold; normalize, precise)
function isotopologues_elements_iter(element_dictionary::Dict, msfix, abundance, abtype, threshold; normalize = true, precise = false)
    abtype = abtyped(abtype)
    isempty(element_dictionary) && return (; Element = [element_dictionary], Mass = [mmi(element_dictionary)], Abundance = [abundance], Preab = [1.0]) 
    element_isotope_pair = element_isotope_pairs(element_dictionary)
    first_proportion = isotopicabundance(element_dictionary; precise)
    max_dictionary = maximal_elements(element_dictionary)
    max_proportion = isotopicabundance(max_dictionary; precise)
    total, abundance_cutoff = abundance_threshold(abtype, abundance, threshold, first_proportion, max_proportion)
    if abtype == Input() && first_proportion / max_proportion < min_propotion() 
        throw(ArgumentError("Isotopic abundance of input chemical is too small; try use `abtype` other than `Input()` or larger threshold`"))
    end
    if max_proportion * abundance_cutoff / total / first_proportion ^ 2 > 1 
        # abundance_cutoff / total << first_proportion
        element_chemical = [get_isotope_vec(max_dictionary)]
        mass_chemical = [mmi(max_dictionary) + msfix]
        abundance_chemical = [total * max_proportion]
        preab_chemical = [isotopologue_inverse_combination(max_dictionary; precise)]
        rec_addminusisotopes_iter!(element_chemical, mass_chemical, abundance_chemical, preab_chemical, max_dictionary, element_isotope_pair, 1, first(mass_chemical), first(abundance_chemical), first(preab_chemical), abundance_cutoff, true, true, precise)
    else
        element_chemical = [get_isotope_vec(element_dictionary)]
        mass_chemical = [mmi(element_dictionary) + msfix]
        abundance_chemical = [total * first_proportion]
        preab_chemical = [1.0]
        rec_addisotopes_iter!(element_chemical, mass_chemical, abundance_chemical, preab_chemical, element_dictionary, element_isotope_pair, 1, first(mass_chemical), first(abundance_chemical), first(preab_chemical), abundance_cutoff, precise)
    end
    abundance_chemical_normalize = normalize_abundance(abundance_chemical, abundance, preabtype(abtype), [Max(), Input(), Total()])
    id = sortperm(mass_chemical)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(abundance_chemical_normalize)))
    filter!(x -> >=(abundance_chemical_normalize[x], abundance_cutoff), id)
    if normalize && dopostnormalize(abtype)
        abundance_chemical = normalize_abundance(abundance_chemical_normalize[id], abundance, abtype, [Max(), Input(), List(), Total()])
    elseif normalize 
        abundance_chemical = abundance_chemical_normalize[id]
    else
        abundance_chemical = abundance_chemical[id]
    end
    (; Element = element_chemical[id], Mass = mass_chemical[id], Abundance = abundance_chemical, Preab = preab_chemical[id]) 
end

# ==========================================================================================================================
# Low level MS1
function rec_addisotopes!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        element_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_abundance, 
        threshold,
        precise
    )
    current_state = false
    rec_state = false
    next_state = true
    max_abundance = prev_abundance
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair)
        next_state, next_abundance = rec_addisotopes!(element_vec, mass_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_abundance, threshold, precise)
        max_abundance = max(max_abundance, next_abundance)
    end
    ne = get(element_dictionary, e, 0)
    if ne > 0
        abundance = update_abundance(prev_abundance, e, i, ne, get(element_dictionary, i, 0), 1; precise)
        max_abundance = max(max_abundance, abundance)
        current_state = (abundance >= threshold)
        rec_state = current_state || rec_state || (abundance >= prev_abundance)
        if !rec_state && (prev_abundance >= abundance > 0)
            rec_state = max_abundance / prev_abundance * abundance >= threshold 
        end
        if rec_state
            element_dictionary[e] -= 1
            get!(element_dictionary, i, 0)
            element_dictionary[i] += 1
            new_mass = prev_mass + element_mass_delta(e, i)
            if current_state
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(abundance_vec, abundance)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_abundance = rec_addisotopes!(element_vec, mass_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position, new_mass, abundance, threshold, precise)
            element_dictionary[e] += 1
            element_dictionary[i] -= 1
            current_state = rec_state || current_state || next_state
            max_abundance = max(max_abundance, rec_abundance)
        end
    end
    current_state, max_abundance
end

function rec_addminusisotopes!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        element_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_abundance, 
        threshold, 
        backward,
        forward, 
        precise
    )
    next_state = true
    max_abundance = prev_abundance
    backward_state = false
    rec_state = false
    backward_abundance = prev_abundance
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair)
        next_state, next_abundance = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_abundance, threshold, true, true, precise)
        max_abundance = max(max_abundance, next_abundance)
    end
    ne = get(element_dictionary, e, 0)
    ni = get(element_dictionary, i, 0)
    if backward && ni > 0 
        abundance = update_abundance(prev_abundance, i, e, ni, ne, 1; precise)
        backward_abundance = max(max_abundance, abundance)
        backward_state = (abundance >= threshold)
        rec_state = backward_state || rec_state || (abundance >= prev_abundance)
        if !rec_state && (prev_abundance >= abundance > 0)
            rec_state = backward_abundance / prev_abundance * abundance >= threshold 
        end
        if rec_state
            element_dictionary[i] -= 1
            get!(element_dictionary, e, 0)
            element_dictionary[e] += 1
            new_mass = prev_mass + element_mass_delta(i, e)
            if backward_state
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(abundance_vec, abundance)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_abundance = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position, new_mass, abundance, threshold, true, false, precise)
            element_dictionary[i] += 1
            element_dictionary[e] -= 1
            backward_state = rec_state || backward_state 
            backward_abundance = max(backward_abundance, rec_abundance)
        end
    end
    forward_state = false
    rec_state = false
    forward_abundance = prev_abundance
    if forward && ne > 0 
        abundance = update_abundance(prev_abundance, e, i, ne, ni, 1; precise)
        forward_abundance = max(max_abundance, abundance)
        forward_state = (abundance >= threshold)
        rec_state = forward_state || rec_state || (abundance >= prev_abundance)
        if !rec_state && (prev_abundance >= abundance > 0)
            rec_state = forward_abundance / prev_abundance * abundance >= threshold 
        end
        if rec_state
            element_dictionary[e] -= 1
            get!(element_dictionary, i, 0)
            element_dictionary[i] += 1
            new_mass = prev_mass + element_mass_delta(e, i)
            if forward_state
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(abundance_vec, abundance)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_abundance = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position, new_mass, abundance, threshold, false, true, precise)
            element_dictionary[e] += 1
            element_dictionary[i] -= 1
            forward_state = rec_state || forward_state 
            forward_abundance = max(forward_abundance, rec_abundance)
        end
    end
    forward_state || backward_state || next_state, max(forward_abundance, backward_abundance)
end

function rec_addisotopes_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_abundance, 
        prev_preab,
        threshold,
        precise
    )
    current_state = false
    rec_state = false
    next_state = true
    max_abundance = prev_abundance
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair)
        next_state, next_abundance = rec_addisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_abundance, prev_preab, threshold, precise)
        max_abundance = max(max_abundance, next_abundance)
    end
    ne = get(element_dictionary, e, 0)
    if ne > 0
        abundance = update_abundance(prev_abundance, e, i, ne, get(element_dictionary, i, 0), 1; precise)
        preab = update_inverse_proportion(prev_preab, ne, get(element_dictionary, i, 0), 1; precise)
        max_abundance = max(max_abundance, abundance)
        current_state = (abundance >= threshold)
        rec_state = current_state || rec_state || (abundance >= prev_abundance)
        if !rec_state && (prev_abundance >= abundance > 0)
            rec_state = max_abundance / prev_abundance * abundance >= threshold 
        end
        if rec_state
            element_dictionary[e] -= 1
            get!(element_dictionary, i, 0)
            element_dictionary[i] += 1
            new_mass = prev_mass + element_mass_delta(e, i)
            if current_state
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(abundance_vec, abundance)
                push!(preab_vec, preab)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_abundance = rec_addisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_dictionary, element_isotope_pair, isotope_position, new_mass, abundance, preab, threshold, precise)
            element_dictionary[e] += 1
            element_dictionary[i] -= 1
            current_state = rec_state || current_state || next_state
            max_abundance = max(max_abundance, rec_abundance)
        end
    end
    current_state, max_abundance
end

function rec_addminusisotopes_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_abundance, 
        prev_preab, 
        threshold, 
        backward,
        forward, 
        precise
    )
    next_state = true
    max_abundance = prev_abundance
    backward_state = false
    rec_state = false
    backward_abundance = prev_abundance
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair)
        next_state, next_abundance = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_abundance, prev_preab, threshold, true, true, precise)
        max_abundance = max(max_abundance, next_abundance)
    end
    ne = get(element_dictionary, e, 0)
    ni = get(element_dictionary, i, 0)
    if backward && ni > 0 
        abundance = update_abundance(prev_abundance, i, e, ni, ne, 1; precise)
        preab = update_inverse_proportion(prev_preab, ni, ne, 1; precise)
        backward_abundance = max(max_abundance, abundance)
        backward_state = (abundance >= threshold)
        rec_state = backward_state || rec_state || (abundance >= prev_abundance)
        if !rec_state && (prev_abundance >= abundance > 0)
            rec_state = backward_abundance / prev_abundance * abundance >= threshold 
        end
        if rec_state
            element_dictionary[i] -= 1
            get!(element_dictionary, e, 0)
            element_dictionary[e] += 1
            new_mass = prev_mass + element_mass_delta(i, e)
            if backward_state
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(abundance_vec, abundance)
                push!(preab_vec, preab)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_abundance = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_dictionary, element_isotope_pair, isotope_position, new_mass, abundance, preab, threshold, true, false, precise)
            element_dictionary[i] += 1
            element_dictionary[e] -= 1
            backward_state = rec_state || backward_state 
            backward_abundance = max(backward_abundance, rec_abundance)
        end
    end
    forward_state = false
    rec_state = false
    forward_abundance = prev_abundance
    if forward && ne > 0 
        abundance = update_abundance(prev_abundance, e, i, ne, ni, 1; precise)
        preab = update_inverse_proportion(prev_preab, ne, ni, 1; precise)
        forward_abundance = max(max_abundance, abundance)
        forward_state = (abundance >= threshold)
        rec_state = forward_state || rec_state || (abundance >= prev_abundance)
        if !rec_state && (prev_abundance >= abundance > 0)
            rec_state = forward_abundance / prev_abundance * abundance >= threshold 
        end
        if rec_state
            element_dictionary[e] -= 1
            get!(element_dictionary, i, 0)
            element_dictionary[i] += 1
            new_mass = prev_mass + element_mass_delta(e, i)
            if forward_state
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(abundance_vec, abundance)
                push!(preab_vec, preab)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_abundance = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_dictionary, element_isotope_pair, isotope_position, new_mass, abundance, preab, threshold, false, true, precise)
            element_dictionary[e] += 1
            element_dictionary[i] -= 1
            forward_state = rec_state || forward_state 
            forward_abundance = max(forward_abundance, rec_abundance)
        end
    end
    forward_state || backward_state || next_state, max(forward_abundance, backward_abundance)
end