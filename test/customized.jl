# Interface AbstractChemical
abstract type Hexose <: AbstractChemical end
struct Glucose <: Hexose
    chirality::String
    nD::Int
    n13C::Int
    retentiontime::Float64
end
struct Galactose <: Hexose
    chirality::String
    nD::Int
    n13C::Int
    retentiontime::Float64
end
struct Mannose <: Hexose 
    chirality::String
    nD::Int
    n13C::Int
    retentiontime::Float64
end

repr_isotope(m::Hexose) = repr_isotope(m.nD, m.n13C)
function repr_isotope(nD, n13C)
    if nD == 0 && n13C == 0
        ""
    elseif nD == 0
        string("[13C", n13C, "]")
    elseif n13C == 0
        string("[D", nD, "]")
    else
        string("[13C", n13C, ",D", nD, "]")
    end
end

chemicalname(m::T; kwargs...) where {T <: Hexose} = string(m.chirality, "-", string(T), repr_isotope(m))
function chemicalelements(m::Hexose; kwargs...) 
    [
        "C"     => 6 - m.n13C,
        "[13C]" => m.n13C,
        "H"     => 12 - m.nD, 
        "D"     => m.nD,
        "O"     => 6
    ]
end
# function getchemicalattr(m::Hexose, ::Val{:formula}) 
#     nd = m.nD 
#     n13c = m.n13C
#     nc = 6 - n13c
#     nh = 12 - nd
#     no = 6
#     rc = string("C", nc > 1 ? nc : "", n13c > 0 ? string("[13C]", n13c > 1 ? n13c : "") : "") 
#     rh = string("H", nh > 1 ? nh : "", nd > 0 ? string("[2H]", nd > 1 ? nd : "") : "") 
#     string(rc, rh, "O", no)
# end

chemicalabbr(m::Hexose; kwargs...) = string(first(chemicalname(m), 5), repr_isotope(m))
chemicalabbr(m::Glucose; kwargs...) = string(first(chemicalname(m), 4), "c", repr_isotope(m))
chemicalsmiles(m::Hexose; kwargs...) = ""

# Fatty acid 
struct FattyAcid <: AbstractChemical
    ncb::Int 
    ndb::Int
    nD::Int 
    n13C::Int 
end

# PS     
struct DiacylPS <: AbstractChemical
    headgroup::Tuple{Int, Int}
    fa1::FattyAcid
    fa2::FattyAcid
    retentiontime::Float64
end
repr_headgroup(m::DiacylPS) = string("PS", repr_isotope(m.headgroup...))
repr_fa(fa::FattyAcid) = string(fa.ncb, ":", fa.ndb, repr_isotope(fa.nD, fa.n13C))
chemicalname(m::FattyAcid; kwargs...) = string("FA", " ", repr_fa(m))
chemicalname(m::DiacylPS; kwargs...) = string(repr_headgroup(m), " ", repr_fa(m.fa1), "/", repr_fa(m.fa2))
function chemicalelements(m::FattyAcid; kwargs...)
    filter(x -> last(x) > 0, [
        "C"     => chemicalncb(m) - chemicaln13C(m), 
        "[13C]" => chemicaln13C(m),
        "H"     => chemicalncb(m) * 2 - chemicalndb(m) * 2 - chemicalnD(m),
        "D"     => chemicalnD(m),
        "O"     => 2
    ])
end
chemicaln13C(m::FattyAcid; kwargs...) = getchemicalproperty(m, :n13C)
chemicalncb(m::FattyAcid; kwargs...) = getchemicalproperty(m, :ncb)
chemicalnD(m::FattyAcid; kwargs...) = getchemicalproperty(m, :nD)
chemicalndb(m::FattyAcid; kwargs...) = getchemicalproperty(m, :ndb)

# function getchemicalattr(m::DiacylPS, ::Val{:elements}; kwargs...)
#     filter(x -> last(x) > 0, [
#         "C"     => 6 + getchemicalattr(m, :ncb) - getchemicalattr(m, :n13C), 
#         "[13C]" => getchemicalattr(m, :n13C),
#         "H"     => 12 + getchemicalattr(m, :ncb) * 2 - getchemicalattr(m, :ndb) * 2 - getchemicalattr(m, :nD),
#         "D"     => getchemicalattr(m, :nD),
#         "N"     => 1, 
#         "O"     => 10, 
#         "P"     => 1
#     ])
# end
function chemicalelements(m::DiacylPS; kwargs...)
    fa1 = chemicalelements(m.fa1)
    fa2 = chemicalelements(m.fa2)
    i = findfirst(x -> ==(first(x), "H"), fa1)
    fa1[i] = "H" => (last(fa1[i]) - 1)
    i = findfirst(x -> ==(first(x), "H"), fa2)
    fa2[i] = "H" => (last(fa2[i]) - 1)
    i = findfirst(x -> ==(first(x), "O"), fa1)
    fa1[i] = "O" => (last(fa1[i]) - 1)
    i = findfirst(x -> ==(first(x), "O"), fa2)
    fa2[i] = "O" => (last(fa2[i]) - 1)
    filter(x -> last(x) > 0, vcat(["C" => (6 - m.headgroup[2]), "[13C]" => m.headgroup[2], "H" => (12 - m.headgroup[1]), "D" => m.headgroup[1], "N" => 1, "O" => 8, "P" => 1], fa1, fa2))
end
chemicalformula(m::DiacylPS; unique = false, kwargs...) = chemicalformula(chemicalelements(m); unique, kwargs...)
# getchemicalattr(m::DiacylPS, ::Val{:formula}) = 
#     string("C", 6 + getchemicalattr(m, :ncb) - getchemicalattr(m, :n13C), 
#         getchemicalattr(m, :n13C) > 0 ? string("[13C]", getchemicalattr(m, :n13C)) : "",
#         "H", 12 + getchemicalattr(m, :ncb) * 2 - getchemicalattr(m, :ndb) * 2 - getchemicalattr(m, :nD),
#         getchemicalattr(m, :nD) > 0 ? string("D", getchemicalattr(m, :nD)) : "",
#         "N", "O", 10, "P"
#     )
chemicalsmiles(m::FattyAcid; kwargs...) = ""
chemicalsmiles(m::DiacylPS; kwargs...) = ""
chemicalnD(m::DiacylPS; kwargs...) = +(m.headgroup[1], chemicalnD(m.fa1; kwargs...), chemicalnD(m.fa2; kwargs...))
chemicaln13C(m::DiacylPS; kwargs...) = +(m.headgroup[2], chemicaln13C(m.fa1; kwargs...), chemicaln13C(m.fa2; kwargs...))
chemicalncb(m::DiacylPS; kwargs...) = chemicalncb(m.fa1; kwargs...) + chemicalncb(m.fa2; kwargs...)
chemicalndb(m::DiacylPS; kwargs...) = chemicalndb(m.fa1; kwargs...) + chemicalndb(m.fa2; kwargs...)

# Interface AbstractAdduct
struct DeSerine <: AbstractAdduct end
set_adduct!("[M-Ser]-", DeSerine())
adductelements(::DeSerine) = ["C" => -3, "H" => -6, "N" => -1, "O" => -2]
adductformula(::DeSerine) = "-Ser"
kmer(::DeSerine) = 1
charge(::DeSerine) = -1

struct Halfprotonation <: AbstractAdduct end
set_adduct!("[2M+H]+", Halfprotonation())
adductelements(::Halfprotonation) = ["H" => 1]
adductformula(::Halfprotonation) = "+H"
kmer(::Halfprotonation) = 2
charge(::Halfprotonation) = 1

# Interface AbstractIon
adductisotopes(ion::AdductIon{DiacylPS, DeSerine}) = ["H" => ioncore(ion).headgroup[1], "D" => -ioncore(ion).headgroup[1], "C" => ioncore(ion).headgroup[2], "[13C]" => -ioncore(ion).headgroup[2]]

glc = Glucose("D", 0, 0, 1.5)
gld = Glucose("D", 6, 0, 1.5)
ps = DiacylPS((0, 0), FattyAcid(18, 0, 0, 0), FattyAcid(20, 4, 0, 0), 7.8)
psi1 = DiacylPS((3, 3), FattyAcid(18, 0, 0, 0), FattyAcid(20, 4, 0, 0), 7.8)
psi2 = DiacylPS((0, 0), FattyAcid(18, 0, 5, 0), FattyAcid(20, 4, 0, 0), 7.8)
iglc = [AdductIon(glc, Protonation())]
igld = [AdductIon(gld, "[M+H]+")]
ips = [AdductIon(ps, "[M-Ser]-"), AdductIon(ps, "[2M+H]+")]
ipsi1 = [AdductIon(psi1, "[M-Ser]-"), AdductIon(psi1, "[2M+H]+")]
ipsi2 = [AdductIon(psi2, "[M-Ser]-"), AdductIon(psi2, "[2M+H]+")]
it3 = Isotopologues(ipsi2[1]; abtype = :total, threshold = crit(1e-3, 1e-3))
git3 = group_isotopologues(it3)