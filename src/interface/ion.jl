in(x::AbstractIon, y::IonCluster) = any(i -> isionequal(i, x), y)
length(ions::IonCluster) = length(ions.ions)

abundantion(ion::AbstractIon) = ion
abundantion(ions::IonCluster) = first(ions.ions)
isionequal(x::IonCluster, y::IonCluster) = all(isionequal(a, b) for (a, b) in zip(x.ions, y.ions)) && all(isapprox(a, b) for (a, b) in zip(x.abundance, y.abundance))
isionequal(x::IonCluster, y::AbstractIon) = length(x) == 1 && isionequal(abundantion(x), y)
isionequal(x::AbstractIon, y::IonCluster) = length(y) == 1 && isionequal(x, abundantion(y))
isionequal(x::AbstractIon, y::AbstractIon) = isadductequal(ionadduct(x), ionadduct(y)) && ischemicalequal(ioncore(x), ioncore(y))

chemicalname(ion::AbstractIon) = string(ionadduct(ion), " of ", chemicalname(ioncore(ion)))
function chemicalformula(ion::AbstractIon)
    el = parse_compound(interpret_isotope(chemicalformula(ioncore(ion))))
    nm = kmer(ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    ael = adductelement(ion)
    for e in ael
        id = findfirst(x -> first(x) == first(e), el)
        if isnothing(id)
            push!(el, e)
        else
            el[id] = first(el[id]) => (last(el[id]) + last(e))
        end
    end
    mapreduce(*, el) do e
        last(e) == 1 ? first(e) : string(first(e), last(e))
    end
end

ioncore(ion::Ion) = ion.core
ionadduct(ion::Ion) = ion.adduct

kmer(ion::AbstractIon) = kmer(ionadduct(ion))
charge(ion::AbstractIon) = charge(ionadduct(ion))
adductelement(ion::AbstractIon) = adductelement(ionadduct(ion))

rt(ion::AbstractIon) = rt(ioncore(ion))
rt(ion::IonCluster) = sum(map(i -> rt(i), ion.ions) .* ion.abundance) / sum(ion.abundance)

function transform_ions(ion::AbstractIon, formulas::Vector{String}, abundance; names = :a)
    if isnothing(names)
        IonCluster([Ion(Chemical(chemicalname(ioncore(ion)), f, collect(infopairs(ioncore(ion)))), ionadduct(ion)) for f in formulas], abundance)
    elseif names == :n
        el = unique_elements(parse_compound(interpret_isotope(chemicalformula(ioncore(ion)))))
        i = findfirst(x -> ==(unique_elements(parse_compound(interpret_isotope(x))), el), formulas)
        if isnothing(i)
            newname = chemicalname(ioncore(ion)) .* string.(1:length(formulas))
            IonCluster([Ion(Chemical(n, f, collect(infopairs(ioncore(ion)))), ionadduct(ion)) for (n, f) in zip(newname, formulas)], abundance)
        else
            newname = chemicalname(ioncore(ion)) .* string.(1:length(formulas) - 1)
            insert!(newname, i, chemicalname(ioncore(ion)))
            IonCluster([Ion(Chemical(n, f, collect(infopairs(ioncore(ion)))), ionadduct(ion)) for (n, f) in zip(newname, formulas)], abundance)
        end
    elseif names == :a
        n = 1
        oldname = chemicalname(ioncore(ion))
        oldformula = chemicalformula(ioncore(ion))
        el = unique_elements(parse_compound(interpret_isotope(oldformula)))
        re = [Regex(string(first(e), "(?!\\])\\d*")) => "" for e in el]
        newname = map(formulas) do f
            replace(f, re...)
        end
        for (i, x) in enumerate(newname)
            if isempty(x)
                if unique_elements(parse_compound(interpret_isotope(formulas[i]))) == el
                    newname[i] = oldname
                else
                    newname[i] = string(oldname, "_", n)
                    n += 1
                end
            else
                newname[i] = string(oldname, " ", x)
            end
        end
        IonCluster([Ion(Chemical(n, f, collect(infopairs(ioncore(ion)))), ionadduct(ion)) for (n, f) in zip(newname, formulas)], abundance)
    else
        IonCluster([Ion(Chemical(n, f, collect(infopairs(ioncore(ion)))), ionadduct(ion)) for (n, f) in zip(names, formulas)], abundance)
    end
end

function transform_ions(ion::AbstractIon, elements::Vector{<: Union{<: Vector{<: Pair}, <: Dict}}, abundance; names = :a)
    if isnothing(names)
        IonCluster([Ion(Chemical(chemicalname(ioncore(ion)), chemicalformula(e), collect(infopairs(ioncore(ion)))), ionadduct(ion)) for e in elements], abundance)
    elseif names == :n
        el = unique_elements(parse_compound(interpret_isotope(chemicalformula(ioncore(ion)))))
        i = findfirst(x -> ==(unique_elements(x), el), elements)
        if isnothing(i)
            newname = chemicalname(ioncore(ion)) .* string.(1:length(elements))
            IonCluster([Ion(Chemical(n, chemicalformula(e), collect(infopairs(ioncore(ion)))), ionadduct(ion)) for (n, e) in zip(newname, elements)], abundance)
        else
            newname = chemicalname(ioncore(ion)) .* string.(1:length(elements) - 1)
            insert!(newname, i, chemicalname(ioncore(ion)))
            IonCluster([Ion(Chemical(n, chemicalformula(e), collect(infopairs(ioncore(ion)))), ionadduct(ion)) for (n, e) in zip(newname, elements)], abundance)
        end
    elseif names == :a
        n = 1
        oldname = chemicalname(ioncore(ion))
        oldformula = chemicalformula(ioncore(ion))
        el = unique_elements(parse_compound(interpret_isotope(oldformula)))
        newname = map(elements) do e
            chemicalformula(Dict(k => v for (k, v) in e if !haskey(el, k)))
        end
        for (i, x) in enumerate(newname)
            if isempty(x)
                if unique_elements(elements[i]) == el
                    newname[i] = oldname
                else
                    newname[i] = string(oldname, "_", n)
                    n += 1
                end
            else
                newname[i] = string(oldname, " ", x)
            end
        end
        IonCluster([Ion(Chemical(n, chemicalformula(e), collect(infopairs(ioncore(ion)))), ionadduct(ion)) for (n, e) in zip(newname, elements)], abundance)
    else
        IonCluster([Ion(Chemical(n, chemicalformula(e), collect(infopairs(ioncore(ion)))), ionadduct(ion)) for (n, e) in zip(names, elements)], abundance)
    end
end