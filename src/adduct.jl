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
