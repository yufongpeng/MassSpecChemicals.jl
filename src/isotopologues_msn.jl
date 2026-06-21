# ==========================================================================================================================
# Mid level MSn
function isotopologues_elements_msn(element_dictionary_vec, msfix_vec, abundance, abtype, threshold; precise = false) 
    max_dictionary_vec = map(maximal_elements, element_dictionary_vec)
    max_proportion_vec = precise ? map(x -> isotopicabundance(x; precise), max_dictionary_vec) : map(isotopicabundance, max_dictionary_vec)
    total, abundance_cutoff = abundance_threshold_msn(abtype, abundance, threshold, max_proportion_vec, element_dictionary_vec) 
    tbls = map(isotopologues_elements_msx, element_dictionary_vec, msfix_vec, max_dictionary_vec, max_proportion_vec, [total for _ in eachindex(element_dictionary_vec)], [abundance_cutoff for _ in eachindex(element_dictionary_vec)], [precise for _ in eachindex(element_dictionary_vec)])
    first(tbls).Abundance .*= total
    els = Vector{Vector{Pair{String, Int}}}[]
    abv = float(Int)[]
    mass = Vector{float(Int)}[]
    rec_vec_ab!(els, abv, mass, tbls, abundance_cutoff, prod(first(tbl.Abundance) for tbl in tbls), [1 for _ in eachindex(tbls)], 1)
    idm = sortperm(mass)
    abv = abv[idm]
    if dopostnormalize(abtype)
        abv = normalize_abundance(abv, abundance, abtype)
    end
    (; Element = els[idm], Mass = mass[idm], Abundance = abv)
end

function isotopologues_elements_msx(element_dictionary, msfix, max_dictionary, max_proportion, abundance_factor, proportioon_cutoff, precise)
    isempty(element_dictionary) && return (; Element = [get_isotope_vec(element_dictionary)], Mass = [mmi(element_dictionary) + msfix], Abundance = [float(1)])
    element_isotope_pair = element_isotope_pairs(element_dictionary)
    element_chemical = [get_isotope_vec(max_dictionary)]
    abundance_chemical = [max_proportion * abundance_factor]
    mass_chemical = [mmi(max_dictionary) + msfix]
    rec_addminusisotopes!(element_chemical, mass_chemical, abundance_chemical, max_dictionary, element_isotope_pair, 1, first(mass_chemical), first(abundance_chemical), proportioon_cutoff, true, true, precise)
    id = sortperm(abundance_chemical; rev = true)
    (; Element = element_chemical[id], Mass = mass_chemical[id], Abundance = [abundance_chemical[i] / abundance_factor for i in id]) 
end

# ==========================================================================================================================
# Low level MSn
function rec_vec_ab!(els, abv, mass, tbls, abundance_cutoff, maxab, id, msn)
    msn == lastindex(tbls) && return rec_vec_ab_end!(els, abv, mass, tbls, abundance_cutoff, maxab, id)
    maxab /= first(tbls[msn].Abundance)
    pass = false
    @inbounds for (i, a) in enumerate(tbls[msn].Abundance)
        a <= 0 && break
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