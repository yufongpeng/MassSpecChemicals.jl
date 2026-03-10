const ADDUCT_ABBR = Dict{Regex, String}()

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
    "[M+H+NH4]2+"   => AddNH4Protonation(),
    "[M+NH4+H]2+"   => AddNH4Protonation(),
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

adducts_doc = """
    adducts_name()
    adducts_abbr()

Access constants related to adducts. 
* `adducts_name`: adduct type for an adduct expression.
* `adducts_abbr`: conversion of abbreviation to valid chemical formula.
"""
@doc adducts_doc
adducts_name() = ADDUCT_NAME

@doc adducts_doc
adducts_abbr() = ADDUCT_ABBR

"""
    set_addabbr!(abbr::AbstractString, fm::AbstractString)

Set `abbr` to be an abbreviation of adduct formula `fm`. 
"""
function set_addabbr!(abbr::AbstractString, fm::AbstractString)
    rabbr = Regex(string("(?<=[+-])", abbr, "(?=[+-\\]])"))
    adducts_abbr()[rabbr] = fm
end
set_addabbr!(abbr::Regex, fm::AbstractString) = push!(adducts_abbr(), abbr => fm)

for (abbr, fm) in [("FA-H", "HCOO"), ("FA", "HCOOH"), ("OFo", "HCOO"), ("FoO", "HCOO"), ("AcA-H", "CH3COO"), ("AcA", "CH3COOH"), ("OAc", "CH3COO"), ("AcO", "CH3COO"), ("Me", "CH3")]
    set_addabbr!(abbr, fm)
end

"""
    set_adduct!(nm::AbstractString, adduct::AbstractAdduct)

Set `nm` to be `adduct` for parsing adduct.

For custumized adduct type, this function is required to make `nm` parsed into `adduct` by `parse_adduct`. 
"""
set_adduct!(nm::AbstractString, adduct::AbstractAdduct) = push!(adducts_name(), nm => adduct)

"""
    chemicalformula(elements::Dictionary; delim = "", unique = true) -> String
    chemicalformula(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}; delim = "", unique = true) -> String

Create chemical formula using given element-number pairs. 

* `delim` assigns the delimiter between each element.
* `unique` determines whether combines the elements to become unique or not.
"""
function chemicalformula(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}; delim = "", unique = true)
    unique ? chemicalformula(unique_elements(elements); delim) : join((v == 1 ? k : string(k, v) for (k, v) in elements if v != 0), delim)
end
function chemicalformula(elements::Dictionary; delim = "", unique = true)
    join((v == 1 ? k : string(k, v) for (k, v) in pairs(elements) if v != 0), delim)
end

"""
    unique_elements(elements::Vector{<: Pair}) -> Dictionary
    unique_elements(elements::Dictionary) -> Dictionary

Create a `Dictionary` from element-number pairs. As elements can be duplicated in the original vector, the new dictionary is convenient for updating elements number.
"""
function unique_elements(elements::Vector{<: Pair{X, Y}}) where {X, Y}
    d = Dictionary{String, Y}()
    for (k, v) in elements
        v == 0 && continue
        get!(d, k, 0)
        d[k] += v
    end
    filter!(!=(0), d)
end
unique_elements(elements::Dictionary) = deepcopy(elements)

"""
    add_elements(elements, y) -> Dictionary

Add elements in `y` to deepcopied `elements`.
"""
add_elements(elements, y) = add_elements!(deepcopy(elements), y)
add_elements(elements::Vector, y) = add_elements!(unique_elements(elements), y)

"""
    add_elements!(elements, y) -> Dictionary

Add elements in `y` to `elements`.
"""
add_elements!(elements::Dictionary, y::Dictionary) = add_elements!(elements, pairs(y))
function add_elements!(elements::Dictionary, y::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary})
    for (k, v) in y
        get!(elements, k, 0)
        elements[k] += v
    end
    filter!(!=(0), elements)
end

"""
    loss_elements(elements, y) -> Dictionary

Substract elements in `y` from deepcopied `elements`.
"""
loss_elements(elements, y) = loss_elements!(deepcopy(elements), y)
loss_elements(elements::Vector, y) = loss_elements!(unique_elements(elements), y)

"""
    loss_elements!(elements, y) -> Dictionary

Substract elements in `y` from `elements`.
"""
loss_elements!(elements::Dictionary, y::Dictionary) = loss_elements!(elements, pairs(y))
function loss_elements!(elements::Dictionary, y::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary})
    for (k, v) in y
        get!(elements, k, 0)
        elements[k] -= v
    end
    filter!(!=(0), elements)
end

function encode_isotopes(formula::AbstractString)
    f = string(formula)
    f2 = f
    for i in eachmatch(r"\[(\d*)([^\]]*)\]", f)
        m, e = i
        delta = isempty(m) ? 0 : (parse(Int, m) - round(Int, ustrip(elements_mass()[e])))
        e = delta > 0 ? string(e, "it") * "n" ^ delta :
            delta < 0 ? string(e, "it") * "p" ^ abs(delta) : string(e, "itz")
        f2 = replace(f2, i.match => e)
    end
    f2
end

"""
    chemicalelements(formula::AbstractString) -> Vector{Pair{String, Int}}

Create element-number pairs from chemical formula. 
"""
function chemicalelements(formula::AbstractString; kwargs...)
    formula = startswith(formula, "-") ? formula[begin + 1:end] : formula
    [elements_decodes()[k] => v for (k, v) in parse_compound(encode_isotopes(formula))]
end

"""
    parse_adduct(::AbstractString) -> AbstractAdduct 

Parse string into `AbstractAdduct`. This function searches string in `ELEMENTS[:NAME]` first, then destructs the input string and constructs a `PosAdduct` or `NegAdduct`.
"""
parse_adduct(adduct::AbstractString) = get(adducts_name(), string(adduct), _parse_adduct(adduct))
function _parse_adduct(adduct::AbstractString)
    adduct = replace(adduct, adducts_abbr()...)
    ion, charge = split(adduct, "]")
    pos = endswith(charge, "+") ? 1 : -1
    charge = split(charge, "+", keepempty = false)
    charge = isempty(charge) ? 1 : begin
        charge = split(only(charge), "-", keepempty = false)
        isempty(charge) ? 1 : parse(Int, only(charge))
    end
    charge *= pos
    ion = replace(ion, "[" => "")
    nm, ion = split(ion, "M"; limit = 2)
    nm = isempty(nm) ? 1 : parse(Int, nm)
    Adduct(nm, ion, charge) 
end
parse_adduct(adduct::AbstractAdduct) = adduct
