const SCHEME_ABBR = Dict{String, AbstractChemical}(
    "H"         => Proton(), 
    "H2O"       => Water(), 
    "NH3"       => Ammonia(), 
    "NH4"       => Ammonium(), 
    "Na"        => Sodium(), 
    "K"         => Potassium(), 
    "Li"        => Lithium(), 
    "Ag"        => Silver(), 
    "CH3COO"    => Acetate(), 
    "OAc"       => Acetate(), 
    "AcO"       => Acetate(), 
    "HCOO"      => Formate(), 
    "OFo"       => Formate(), 
    "FoO"       => Formate(), 
    "F"         => Fluoride(),
    "Cl"        => Chloride(),
    "Me"        => Methenium(),
    "HOAc"      => AceticAcid(),
    "AcOH"      => AceticAcid(),
    "AcA"       => AceticAcid(),
    "CH3COOH"   => AceticAcid(),
    "HOFo"      => FormicAcid(),
    "FoOH"      => FormicAcid(),
    "FoA"       => FormicAcid(),
    "FA"        => FormicAcid(),
    "HCOOH"     => FormicAcid(),
    "MeOAc"     => MethylAcetate(),
    "CH3OAc"    => MethylAcetate(),
    "AcOMe"     => MethylAcetate(),
    "AcOCH3"    => MethylAcetate(),
    "CH3COOCH3" => MethylAcetate(),
    "MeOFo"     => MethylFormate(),
    "CH3OFo"    => MethylFormate(),
    "FoOMe"     => MethylFormate(),
    "FoOCH3"    => MethylFormate(),
    "HCOOCH3"   => MethylFormate()
)

const SCHEME_NAME = Dict{String, AbstractScheme}()

scheme_doc = """
    scheme_name()
    scheme_abbr()

Access constants related to scheme. 
* `scheme_name`: scheme type for an scheme expression.
* `scheme_abbr`: conversion of abbreviation to chemical for scheme.
"""
@doc scheme_doc
scheme_name() = SCHEME_NAME

@doc scheme_doc
scheme_abbr() = SCHEME_ABBR

"""
    set_schabbr!(abbr::AbstractString, chemical::AbstractChemical)

Set `abbr` to be an abbreviation of scheme chemical `fm` for default `AdductParser`. 
"""
set_schabbr!(abbr::AbstractString, chemical::AbstractChemical) = push!(scheme_abbr(), abbr => chemical)

for (abbr, fm) in [
        ("FA", FormicAcid()), 
        ("FoOH", FormicAcid()), 
        ("FoA", FormicAcid()), 
        ("HOFo", FormicAcid()), 
        ("OFo", Formate()), 
        ("FoO", Formate()), 
        ("AcA", AceticAcid()), 
        ("AcOH", AceticAcid()), 
        ("HOAc", AceticAcid()), 
        ("OAc", Acetate()), 
        ("AcO", Acetate()), 
        ("Me", Methenium())]
    set_schabbr!(abbr, fm)
end

"""
    set_scheme!(nm::AbstractString, scheme::AbstractScheme)

Set `nm` to be `scheme` for default `SchemeParser`.

For custumized scheme type, this function is required to make `nm` parsed into `scheme` by `parse_scheme`. 
"""
set_scheme!(nm::AbstractString, scheme::AbstractScheme) = push!(scheme_name(), nm => scheme)

"""
    parse_adduct([adductparser::AbstractAdductParser,] adduct::AbstractString; args = false, kwargs...) 

Parse string into `NamedTuple` or `Tuple` using `adductparser`. 
The default `AdductParser` determines the number of core first (`"[2M+...]"` for 2 cores), and then searches string in `scheme_name()`, or parses the string using `adductparser.chemicalparser`.

# Keyword Arguments
* `args::Bool`: determines whether returning a tuple (`false`) or a named tuple (`true`). Named tuple is useful for user input from a data table and ionization using `ionize`; tuple can be used for fixed adduct ion constructor interface.
"""
parse_adduct(adduct::AbstractString; args = false, kwargs...) = parse_adduct(AdductParser(), adduct; args, kwargs...)
parse_adduct(adduct::AbstractScheme; args = false, kwargs...) = args ? (adduct, 1) : (; adduct, ncore = 1)
function parse_adduct(adductparser::AdductParser, adduct::AbstractString; args = false, kwargs...) 
    nm, sch = split(adduct, "M"; limit = 2)
    nm = nm[begin + 1:end]
    ncore = isempty(nm) ? 1 : parse(Int, nm)
    sch = string("[", sch)
    adduct = haskey(scheme_name(), sch) ? scheme_name()[sch] : parse_chemical(adductparser.chemicalparser, sch)
    args ? (adduct, ncore) : (; adduct, ncore)
end

parse_adduct(x::Tuple; kwargs...) = x
parse_adduct(x::NamedTuple; kwargs...) = x

# parse_adduct(adduct::AbstractChemicalScheme) = (adduct, 1)
# parse_adduct(::AbstractAdductParser, adduct::AbstractChemicalScheme) = (adduct, 1)
