# Coelution(separation, exp, lib, table)
# table: index of target in mztable => mztable
# (modify abundance)
# Add mz
# ms analyze of coeluting chemicals 
# ...
# last ms scan 
# peak_table
# error_table
function CoelutingIsobars(separation, ms, exp::Table, lib::Table; error_sep = [[value_error, percentage_error] for x in separation], error_mz = [[value_error, ppm_error] for x in ms])
    exptransitions = chemicaltransitions.(exp.Chemical)
    exp_icol = findall(x -> startswith(x, "Abundance"), string.(propertynames(exp)))
    exp_col = sort(propertynames(exp)[exp_icol]; by = x -> parse(Int, last(split(string(x)))))
    exp_ab = [getproperty(exp, x) for x in exp_col]
    libtransitions = chemicaltransitions.(lib.Chemical)
    lib_icol = findall(x -> startswith(x, "Abundance"), string.(propertynames(lib)))
    lib_col = sort(propertynames(lib)[lib_icol]; by = x -> parse(Int, last(split(string(x)))))
    lib_ab = [getproperty(lib, x) for x in lib_col]
    exp_sep_value = [[fn(first(c)) for c in exptransitions] for (fn, _) in separation]
    lib_sep_value = [[fn(first(c)) for c in libtransitions] for (fn, _) in separation]
    exp_mz_value = [[mz(c[i]) for c in exptransitions] for i in eachindex(first(exptransitions))]
    lib_mz_value = [[mz(c[i]) for c in libtransitions] for i in eachindex(first(libtransitions))]
    quals = sep_qual!(separation, exp_sep_value, lib_sep_value)
    mz_qual!(ms, exp_mz_value, exp_ab, lib_mz_value, lib_ab, quals)
    cols = Pair{Symbol, Vector}[]
    id = Pair[]
    for (i, qual) in enumerate(quals)
        all(!, qual) && continue 
        push!(id, i => findall(qual))
    end
    push!(cols, :Chemical => vcat(([exp.Chemical[first(i)] for _ in last(i)] for i in id)...))
    push!(cols, :Isobar => vcat((lib.Chemical[last(i)] for i in id)...))
    for (i, fns) in enumerate(error_sep)
        sep_value = vcat(([(exp_sep_value[i][first(j)], lib_sep_value[i][k]) for k in last(j)] for j in id)...)
        for fn in fns 
            push!(cols, Symbol(colname(first(separation[i]), fn)) => [fn(x...) for x in sep_value])
        end
    end
    for (i, fns) in enumerate(error_mz)
        mz_value = vcat(([(exp_mz_value[i][first(j)], lib_mz_value[i][k]) for k in last(j)] for j in id)...)
        for fn in fns 
            push!(cols, Symbol(colname(string("MZ", i), fn)) => [fn(x...) for x in mz_value]) 
        end 
    end 
    Table(; cols...)
end

function mz_qual!(ms::AbstractVector, exp_mz_value, exp_ab, lib_mz_value, lib_ab, quals = [trues(length(first(lib_mz_value)) for i in eachindex(first(exp_mz_value)))])
    all(x -> all(!, x), quals) && return quals
    if isempty(exp_ab)
        exp_ab = [repeat([1], length(exp_mz_value)) for _ in eachindex(first(exp_mz_value))]
    end
    if isempty(lib_ab)
        lib_ab = [repeat([1], length(lib_mz_value)) for _ in eachindex(first(lib_mz_value))]
    end
    i = 0
    while i < lastindex(ms)
        i += 1 
        i1 = min(lastindex(exp_mz_value), i)
        i2 = min(lastindex(exp_ab), i)
        i3 = min(lastindex(lib_mz_value), i)
        i4 = min(lastindex(lib_ab), i)
        for j in eachindex(quals)
            all(!, quals[j]) && continue 
            imz_qual!(
                ms[i], 
                exp_mz_value[i1][j],
                exp_ab[i2][j],
                lib_mz_value[i3],
                lib_ab[i4],
                quals[j]
            )
        end
        any(any, quals) && continue 
    end
    quals
end

function imz_qual!(ms::Criteria, target_mz, target_ab, lib_mz, lib_ab, quals)
    tol = union(makecrit_delta(ms, target_mz)...)
    quals .*= [x in tol for x in lib_mz]
end

function sep_qual!(sep, exp_sep_value, lib_sep_value, quals = [trues(length(first(lib_sep_value))) for i in eachindex(first(exp_sep_value))])
    all(x -> all(!, x), quals) && return quals
    i = 0
    while i < lastindex(sep)
        i += 1 
        i1 = min(lastindex(exp_sep_value), i)
        i2 = min(lastindex(lib_sep_value), i)
        for j in eachindex(quals)
            all(!, quals[j]) && continue 
            isep_qual!(
                last(sep[i]), 
                exp_sep_value[i1][j],
                lib_sep_value[i2],
                quals[j]
            )
        end
        any(any, quals) && continue 
    end
    quals
end

function isep_qual!(sep::Criteria, target_sep, lib_sep, quals)
    tol = union(makecrit_delta(sep, target_sep)...)
    quals .*= [x in tol for x in lib_sep]
end
# lex(v1, v2, i = 1) = i > lastindex(v1) || (i <= lastindex(v2) && (v1[i] == v2[i] ? lex(v1, v2, i + 1) : v1[i] < v2[i])) 
# separation dim: indivisual fwhm / crit
# [rt => Criteria]
# [rt => Resolution(fwhm) => cutoff]
# ms: analyzer dependent / crit
# [msanalyzer => Criteria (abundance)], [Criteria]
# error, percentage error, ppm error, resolution