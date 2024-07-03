const ADDUCT_ALIAS = Dict(
    "FA-H" => "HCOO",
    "OAc" => "CH3COO",
    "AcO" => "CH3COO",
    "-Serine+H2O" => "-C3H5NO2",
    "Serine" => "C3H7NO3",
    "Me" => "CH3"
)

function chemicalformula(x::Union{Vector{<: Pair}, Dict})
    join(begin
        e = get(ELEMENTS, k, k)
        k == e ? string(k, v > 1 ? v : "") : string("[", round(Int, ustrip(MW[k])), e, "]", v > 1 ? v : "")
    end for (k, v) in x)
end

function unique_elements(x::Vector{<: Pair{X, Y}}) where {X, Y}
    d = Dict{String, Y}()
    for (k, v) in x
        get!(d, k, 0)
        d[k] += v
    end
    d
end
unique_elements(x::Dict) = x
function add_elements!(d::Dict, y::Union{Vector{<: Pair}, Dict})
    for (k, v) in y
        get!(d, k, 0)
        d[k] += v
    end
    d
end
function loss_elements!(d::Dict, y::Union{Vector{<: Pair}, Dict})
    for (k, v) in y
        get!(d, k, 0)
        d[k] -= v
    end
    d
end

function interpret_isotope(formula::AbstractString)
    f = string(formula)
    f2 = f
    for i in eachmatch(r"\[(\d*)([^\]]*)\]", f)
        m, e = i
        delta = isempty(m) ? 0 : (parse(Int, m) - round(Int, ustrip(MW[e])))
        e = delta > 0 ? string(e, "it") * "n" ^ delta :
            delta < 0 ? string(e, "it") * "p" ^ delta : e
        f2 = replace(f2, i.match => e)
    end
    f2
end

function parse_adduct(adduct::AbstractString)
    adduct = replace(adduct, ADDUCT_ALIAS...)
    ion, charge = split(adduct, "]")
    pos = endswith(charge, "+")
    charge = split(charge, "+", keepempty = false)
    charge = isempty(charge) ? 1 : begin
        charge = split(only(charge), "-", keepempty = false)
        isempty(charge) ? 1 : parse(Int, only(charge))
    end
    ion = replace(ion, "[" => "")
    nm, ion = split(ion, "M")
    nm = isempty(nm) ? 1 : parse(Int, nm)
    pos ? PosAdduct(nm, ion, charge) : NegAdduct(nm, ion, charge)
end
