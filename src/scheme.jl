const SCHEME_ABBR = Dict{String, AbstractChemical}(
    "H"         => StandardChemical[Proton], 
    "H2O"       => StandardChemical[Water], 
    "NH3"       => StandardChemical[Ammonia],
    "NH4"       => StandardChemical[Ammonium],
    "Na"        => StandardChemical[Sodium],
    "K"         => StandardChemical[Potassium],
    "Li"        => StandardChemical[Lithium],
    "Ag"        => StandardChemical[Silver],
    "CH3COO"    => StandardChemical[Acetate],
    "OAc"       => StandardChemical[Acetate],
    "AcO"       => StandardChemical[Acetate],
    "HCOO"      => StandardChemical[Formate],
    "OFo"       => StandardChemical[Formate],
    "FoO"       => StandardChemical[Formate],
    "F"         => StandardChemical[Fluoride],
    "Cl"        => StandardChemical[Chloride],
    "Me"        => StandardChemical[Methenium],
    "HOAc"      => StandardChemical[AceticAcid],
    "AcOH"      => StandardChemical[AceticAcid],
    "AcA"       => StandardChemical[AceticAcid],
    "CH3COOH"   => StandardChemical[AceticAcid],
    "HOFo"      => StandardChemical[FormicAcid],
    "FoOH"      => StandardChemical[FormicAcid],
    "FoA"       => StandardChemical[FormicAcid],
    "FA"        => StandardChemical[FormicAcid],
    "HCOOH"     => StandardChemical[FormicAcid],
    "MeOAc"     => StandardChemical[MethylAcetate],
    "CH3OAc"    => StandardChemical[MethylAcetate],
    "AcOMe"     => StandardChemical[MethylAcetate],
    "AcOCH3"    => StandardChemical[MethylAcetate],
    "CH3COOCH3" => StandardChemical[MethylAcetate],
    "MeOFo"     => StandardChemical[MethylFormate],
    "CH3OFo"    => StandardChemical[MethylFormate],
    "FoOMe"     => StandardChemical[MethylFormate],
    "FoOCH3"    => StandardChemical[MethylFormate],
    "HCOOCH3"   => StandardChemical[MethylFormate]
)

const SCHEME_NAME = Dict{String, AbstractScheme}()

scheme_doc = """
    scheme_name()
    scheme_abbr()

Access constants related to a scheme.
* `scheme_name`: the scheme type for a scheme expression.
* `scheme_abbr`: conversion from an abbreviation to a chemical in a scheme.
"""
@doc scheme_doc
scheme_name() = SCHEME_NAME

@doc scheme_doc
scheme_abbr() = SCHEME_ABBR

"""
    set_schabbr!(abbr::AbstractString, chemical::AbstractChemical)

Set `abbr` as the abbreviation of scheme `chemical` for the default [`AdductParser`](@ref).

# Examples
```julia
julia> serine = Chemical("Serine", "C3H5NO2"; abbreviation = "Ser") 
Serine

julia> set_schabbr!("Ser", serine)
Dict{String, AbstractChemical} with 37 entries:
  "HCOOH"     => Formic Acid
  "Me"        => Methenium
  ⋮           => ⋮

julia> parse_adduct("[M-H-Ser]-")
(adduct = Loss_Proton|Loss_Serine, ncore = 1)
```
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
    set_scheme!(name::AbstractString, scheme::AbstractScheme)

Set `name` as the scheme name for the default [`AdductParser`](@ref).

For custom scheme types, this function is required so that `name` is parsed into `scheme` by [`parse_adduct`](@ref).

# Examples
```julia
julia> struct NewScheme <: AbstractCompleteScheme end 

julia> set_scheme!("[+NSC]+", NewScheme())
Dict{String, AbstractScheme} with 1 entry:
  "[+NSC]+" => NewScheme

julia> parse_adduct("[2M+NSC]+")
(adduct = NewScheme, ncore = 2)
```
"""
set_scheme!(name::AbstractString, scheme::AbstractScheme) = push!(scheme_name(), name => scheme)