ischemicalequal(x::AbstractChemical, y::AbstractChemical) = isequal(x, y)
ischemicalequal(x::Chemical, y::Chemical) = isequal(chemicalname(x), chemicalname(y)) || isequal(unique_elements(parse_compound(interpret_isotope(chemicalformula(x)))), unique_elements(parse_compound(interpret_isotope(chemicalformula(y)))))

chemicalname(m::Chemical) = m.name
chemicalformula(m::Chemical) = m.formula
chemicalformula(m::AbstractChemical, adduct::AbstractAdduct) = chemicalformula(Ion(m, adduct))
chemicalformula(m::AbstractChemical, adduct::AbstractString) = chemicalformula(Ion(m, parse_adduct(adduct)))

infonames(m::Chemical) = first.(m.info)
function getinfo(m::Chemical, i::Symbol)
    for (p, v) in m.info
        p == i && return v
    end
    return nothing
end
infopairs(m::Chemical) = m.info

function rt(m::Chemical)
    for (p, v) in m.info
        p == :rt && return v
    end
    return nothing
end

# nothing, nochange
# :a, append
# :n, number
# vector
function transform_chemicals(chemical::AbstractChemical, formulas::Vector{String}; names = :a)
    if isnothing(names)
        [Chemical(chemicalname(chemical), f, collect(infopairs(chemical))) for f in formulas]
    elseif names == :n
        el = unique_elements(parse_compound(interpret_isotope(chemicalformula(chemical))))
        i = findfirst(x -> ==(unique_elements(parse_compound(interpret_isotope(x))), el), formulas)
        if isnothing(i)
            newname = chemicalname(chemical) .* string.(1:length(formulas))
            [Chemical(n, f, collect(infopairs(chemical))) for (n, f) in zip(newname, formulas)]
        else
            newname = chemicalname(chemical) .* string.(1:length(formulas) - 1)
            insert!(newname, i, chemicalname(ioncore(ion)))
            [Chemical(n, f, collect(infopairs(chemical))) for (n, f) in zip(newname, formulas)]
        end
    elseif names == :a
        n = 1
        oldname = chemicalname(chemical)
        oldformula = chemicalformula(chemical)
        el = unique_elements(parse_compound(interpret_isotope(oldformula)))
        re = [Regex(string(k, "(?!\\])\\d*")) => "" for (k, _) in el]
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
        [Chemical(n, f, collect(infopairs(chemical))) for (n, f) in zip(newname, formulas)]
    else
        [Chemical(n, f, collect(infopairs(chemical))) for (n, f) in zip(names, formulas)]
    end
end

function transform_chemicals(chemical::AbstractChemical, elements::Vector{<: Union{<: Vector{<: Pair}, <: Dict}}; names = :a)
    if isnothing(names)
        [Chemical(chemicalname(chemical), chemicalformula(e), collect(infopairs(chemical))) for e in elements]
    elseif names == :n
        el = unique_elements(parse_compound(interpret_isotope(chemicalformula(chemical))))
        i = findfirst(e -> ==(unique_elements(e), el), elements)
        if isnothing(i)
            newname = chemicalname(chemical) .* string.(1:length(elements))
            [Chemical(n, chemicalformula(e), collect(infopairs(chemical))) for (n, e) in zip(newname, elements)]
        else
            newname = chemicalname(chemical) .* string.(1:length(elements) - 1)
            insert!(newname, i, chemicalname(ioncore(ion)))
            [Chemical(n, chemicalformula(e), collect(infopairs(chemical))) for (n, e) in zip(newname, elements)]
        end
    elseif names == :a
        n = 1
        oldname = chemicalname(chemical)
        oldformula = chemicalformula(chemical)
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
        [Chemical(n, chemicalformula(e), collect(infopairs(chemical))) for (n, e) in zip(newname, elements)]
    else
        [Chemical(n, chemicalformula(e), collect(infopairs(chemical))) for (n, e) in zip(names, elements)]
    end
end
