const ADDUCT_ABBR = Dict{Regex, String}()

function push_adduct_abbr!(abbr::AbstractString, fm::AbstractString)
    rabbr = Regex(string("(?<=[+-])", abbr, "(?=[+-\\]])"))
    push!(ADDUCT_ABBR, rabbr => fm)
end
push_adduct_abbr!(abbr::Regex, fm::AbstractString) = push!(ADDUCT_ABBR, abbr => fm)
for (abbr, fm) in [("FA-H", "HCOO"), ("FA", "HCOOH"), ("OFo", "HCOO"), ("FoO", "HCOO"), ("AcA-H", "CH3COO"), ("AcA", "CH3COOH"), ("OAc", "CH3COO"), ("AcO", "CH3COO"), ("Me", "CH3")]
    push_adduct_abbr!(abbr, fm)
end

const ADDUCT_NAME = Dict{String, AbstractAdduct}(
    "[M]+"          => LossElectron(),
    "[M+H]+"        => Protonation(),
    "[M+H-H2O]+"    => ProtonationNLH2O(),
    "[M-H2O+H]+"    => ProtonationNLH2O(),
    "[M+H-2H2O]+"   => ProtonationNL2H2O(),
    "[M-2H2O+H]+"   => ProtonationNL2H2O(),
    "[M+H-3H2O]+"   => ProtonationNL3H2O(),
    "[M-3H2O+H]+"   => ProtonationNL3H2O(),
    "[M+2H]2+"      => DiProtonation(),
    "[M+3H]3+"      => TriProtonation(),
    "[M+NH4]+"      => AddNH4(),
    "[M+H+NH4]2+"   => AddHNH4(),
    "[M+NH4+H]2+"   => AddHNH4(),
    "[M+2NH4]2+"    => Add2NH4(),
    "[M+Na]+"       => Sodization(),
    "[M+Na+H]2+"    => SodizationProtonation(),
    "[M+H+Na]2+"    => SodizationProtonation(),
    "[M+2Na]2+"     => DiSodization(),
    "[M+Na+NH4]2+"  => SodizationAddNH4(),
    "[M+NH4+Na]2+"  => SodizationAddNH4(),
    "[M+K]+"        => Potassiation(),
    "[M+K+H]2+"     => PotassiationProtonation(),
    "[M+H+K]2+"     => PotassiationProtonation(),
    "[M+2K]2+"      => DiPotassiation(),
    "[M+Li]+"       => Lithiation(),
    "[M+Li+H]2+"    => LithiationProtonation(),
    "[M+H+Li]2+"    => LithiationProtonation(),
    "[M+2Li]2+"     => DiLithiation(),
    "[M]-"          => AddElectron(),
    "[M-H]-"        => Deprotonation(),
    "[M-H-H2O]-"    => DeprotonationNLH2O(),
    "[M-H2O-H]-"    => DeprotonationNLH2O(),
    "[M-2H]2-"      => DiDeprotonation(),
    "[M-3H]3-"      => TriDeprotonation(),
    "[M+CH3COO]-"   => AddOAc(),
    "[M+OAc]-"      => AddOAc(),
    "[M+AcO]-"      => AddOAc(),
    "[M+HCOO]-"     => AddOFo(),
    "[M+OFo]-"      => AddOFo(),
    "[M+FoO]-"      => AddOFo(),
    # "[M-CH2O]-"     => LossCH2O(),
    # "[M+O]-"        => AddO(),
    # "[M+C2H2O]-"    => AddC2H2O(),
    # "[M-CH8NO]-"    => LossCH8NO(),
    # "[M-C2H8NO]-"   => LossC2H8NO(),
    # "[M+C3H5NO]-"   => AddC3H5NO(),
    # "[M+C2H5NO]-"   => AddC2H5NO(),
    # "[M-CH3]-"      => LossCH3(),
    # "[M-Me]-"       => LossCH3(),
    # "[M-H-Serine+H2O]-" => DeprotonationLossSerineAddH2O(),
    # "[M-Serine+H2O-H]-" => DeprotonationLossSerineAddH2O(),
    # "[M+H2O-H-Serine]-" => DeprotonationLossSerineAddH2O(),
    # "[M-Serine-H+H2O]-" => DeprotonationLossSerineAddH2O(),
    # "[M-H+H2O-Serine]-" => DeprotonationLossSerineAddH2O(),
    # "[M+H2O-Serine-H]-" => DeprotonationLossSerineAddH2O(),
    "[M+F]-"        => Fluoridation(),
    "[M+Cl]-"       => Chloridation(),
    "[M-Me]-"       => Demethylation()

)

"""
    push_adduct_name!(nm::AbstractString, adduct::AbstractAdduct)

Push the pair `nm` and `adduct` into `ADDUCT_NAME` for parsing specific adduct string.

For custumized adduct type, this function is required to make `nm` parsed into `adduct` by `parse_adduct`. 
"""
push_adduct_name!(nm::AbstractString, adduct::AbstractAdduct) = push!(ADDUCT_NAME, nm => adduct)

function chemicalformula(x::Union{Vector{<: Pair}, Dict}; delim = "")
    join((begin
        e = get(ELEMENTS, k, k)
        k == e ? string(k, v > 1 ? v : "") : string("[", round(Int, ustrip(MW[k])), e, "]", v > 1 ? v : "")
    end for (k, v) in x), delim)
end

function unique_elements(x::Vector{<: Pair{X, Y}}) where {X, Y}
    d = Dict{String, Y}()
    for (k, v) in x
        get!(d, k, 0)
        d[k] += v
    end
    d
end
unique_elements(x::Dict) = deepcopy(x)
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

function transform_isotope_repr(formula::AbstractString)
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

"""
    parse_adduct(::AbstractString) -> AbstractAdduct 

Parse string into `AbstractAdduct`. This function searches string in `ADDUCT_NAME` first, then destructs the input string and constructs a `PosAdduct` OR `NegAdduct`.
"""
parse_adduct(adduct::AbstractString) = get(ADDUCT_NAME, string(adduct), _parse_adduct(adduct))
function _parse_adduct(adduct::AbstractString)
    adduct = replace(adduct, ADDUCT_ABBR...)
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
parse_adduct(adduct::AbstractAdduct) = adduct
