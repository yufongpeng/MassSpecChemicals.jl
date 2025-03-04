using MassSpecChemicals
const MSC = MassSpecChemicals
import MassSpecChemicals: getchemicalattr, adductelement, adductformula, kmer, charge, attrpairs
using Test

# Interface AbstractChemical
abstract type Hexose <: AbstractChemical end
struct Glucose <: Hexose
    chirality::String
    nD::Int
    n13C::Int
    rt::Float64
end
struct Galactose <: Hexose     
    chirality::String
    nD::Int
    n13C::Int
    rt::Float64
end
struct Mannose <: Hexose 
    chirality::String
    nD::Int
    n13C::Int
    rt::Float64
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

getchemicalattr(m::T, ::Val{:name}) where {T <: Hexose} = string(m.chirality, "-", string(T), repr_isotope(m))
function getchemicalattr(m::Hexose, ::Val{:formula}) 
    nd = m.nD 
    n13c = m.n13C
    nc = 6 - n13c
    nh = 12 - nd
    no = 6
    rc = string("C", nc > 1 ? nc : "", n13c > 0 ? string("[13C]", n13c > 1 ? n13c : "") : "") 
    rh = string("H", nh > 1 ? nh : "", nd > 0 ? string("[2H]", nd > 1 ? nd : "") : "") 
    string(rc, rh, "O", no)
end

getchemicalattr(m::Hexose, ::Val{:abbreviation}) = string(first(chemicalname(m), 5), repr_isotope(m))
getchemicalattr(m::Glucose, ::Val{:abbreviation}) = string(first(chemicalname(m), 4), "c", repr_isotope(m))
getchemicalattr(m::Hexose, ::Val{:SMILES}) = ""
attrpairs(m::Hexose) = [s => getchemicalattr(m, s) for s in [:abbreviation, :SMILES, :rt, :chirality]]

struct DiacylPS <: AbstractChemical
    headgroup::Tuple{Int, Int}
    fa1::Tuple{Int, Int, Int, Int}
    fa2::Tuple{Int, Int, Int, Int}
    rt::Float64
end
repr_headgroup(m::DiacylPS) = string("PS", repr_isotope(m.headgroup...))
repr_fa(fa::Tuple) = string(fa[1], ":", fa[2], repr_isotope(fa[3], fa[4]))
getchemicalattr(m::DiacylPS, ::Val{:name}) = string(repr_headgroup(m), " ", repr_fa(m.fa1), repr_fa(m.fa2))
getchemicalattr(m::DiacylPS, ::Val{:formula}) = 
    string("C", 6 + getchemicalattr(m, :ncb) - getchemicalattr(m, :n13C), 
        getchemicalattr(m, :n13C) > 0 ? string("[13C]", getchemicalattr(m, :n13C)) : "",
        "H", 12 + getchemicalattr(m, :ncb) * 2 - getchemicalattr(m, :ndb) * 2 - getchemicalattr(m, :nD),
        getchemicalattr(m, :nD) > 0 ? string("D", getchemicalattr(m, :nD)) : "",
        "O", 10, "P", 1
    )
getchemicalattr(m::DiacylPS, ::Val{:SMILES}) = ""
getchemicalattr(m::DiacylPS, ::Val{:nD}) = +(m.headgroup[1], m.fa1[3], m.fa2[3])
getchemicalattr(m::DiacylPS, ::Val{:n13C}) = +(m.headgroup[2], m.fa1[4], m.fa2[4])
getchemicalattr(m::DiacylPS, ::Val{:ncb}) = +(m.fa1[1], m.fa2[1])
getchemicalattr(m::DiacylPS, ::Val{:ndb}) = +(m.fa1[2], m.fa2[2])
attrpairs(m::DiacylPS) = [s => getchemicalattr(m, s) for s in [:SMILES, :rt, :nD, :n13C, :ncb, :ndb]]

# Interface AbstractAdduct
struct DeSerine <: AbstractNegAdduct end
push_adduct_name!("[M-Ser]-", DeSerine())
adductelement(::DeSerine) = ["C" => 3, "H" => 6, "N" => 1, "O" => 2]
adductformula(::DeSerine) = "-Ser"
kmer(::DeSerine) = 1
charge(::DeSerine) = -1

struct Halfprotonation <: AbstractPosAdduct end
push_adduct_name!("[2M+H]+", Halfprotonation())
adductelement(::Halfprotonation) = ["H" => 1]
adductformula(::Halfprotonation) = "+H"
kmer(::Halfprotonation) = 2
charge(::Halfprotonation) = 1

# Interface AbstractIon
adductisotope(ion::Ion{DiacylPS, DeSerine}) = ["H" => ioncore(ion).headgroup[1], "D" => -ioncore(ion).headgroup[1], "C" => ioncore(ion).headgroup[2], "[13C]" => -ioncore(ion).headgroup[2]]

@testset "MassSpecChemicals.jl" begin
    # Default chemical, adduct
    cglc = Chemical("Glucose", "C6H12O6"; rt = 1.5, abbreviation = "Glc", SMILES = "")
    cgld = parse_chemical("Glucose-d6", "C6H6D6O6"; rt = 1.5, abbreviation = "Glc[D6]", SMILES = "")
    cps = Chemical("PS 18:0/20:4(5Z,8Z,11Z,14Z)", "C44H80NO10P"; rt = 7.8)
    cpsi1 = Chemical("PS[D3,13C3] 18:0/20:4(5Z,8Z,11Z,14Z)", "C41[13C]3H77D3NO10P"; rt = 7.8)
    cpsi2 = Chemical("PS 18:0[D5]/20:4(5Z,8Z,11Z,14Z)", "C44H75D5NO10P"; rt = 7.8)
    lossserine = NegAdduct(1, "-C3H6NO2", 1)
    lossserinei = NegAdduct(1, "-[13C]3H3D3NO2", 1)
    dimh = PosAdduct(2, "+H", 1)
    # test all default adduct 
    # mw, mz, formula
    icglcall = [Ion(cglc, k) for k in keys(MSC.ADDUCT_NAME)]
    icglc = [Ion(cglc, Protonation())]
    icgld = [Ion("Glucose-d6", "C6H6D6O6", "[M+H]+"; rt = 1.5), Ion("Glucose-d6", "C6H6D6O6", Protonation(); rt = 1.5), Ion(cgld, Protonation())]
    icps = [Ion(cps, lossserine), Ion(cps, dimh)]
    icpsi1 = [Ion(cpsi1, lossserinei), Ion(cpsi1, dimh)]
    icpsi2 = [Ion(cpsi2, lossserine), Ion(cpsi2, dimh)]
    @test chemicalname(cglc) == "Glucose"
    @test chemicalname(icgld[1]) == "Glucose-d6[M+H]+"
    @test chemicalformula(cglc) == "C6H12O6"
    @test chemicalformula(icgld[1]) == "C6H7D6O6"
    @test isapprox(mz(icpsi1[1]), mz(icps[1]))
    @test isapprox((mw(icps[2]) * 2 + mw("H")) / 2, mz(icps[2]))
    @test isapprox(rt(icpsi2[1]), 7.8)
    @test chemicalabbr(cglc) == "Glc"
    @test chemicalsmiles(cgld) == ""
    @test ischemicalequal(icgld[1], icgld[2])
    # ionvariants, chemicalvariants
    iv1 = ionvariants(icglc[1], ["C6H12O6", "C6H11DO6", "C6H10D2O6", "C6H9D3O6", "C6H8D4O6", "C6H7D5O6", "C6H6D6O6"])
    @test chemicalname(first(iv1)) == "Glucose[1][M+H]+"
    @test chemicalname(last(iv1)) == "Glucose[[2H]6][M+H]+"
    @test chemicalformula(iv1[2]) == "C6H12DO6"
    # isotopes, isobars
    it = isotopetable(icglc[1], 1e5; threshold = crit(1e3, 1e-2))
    # User defined chemical, adduct
    glc = Glucose("D", 0, 0, 1.5)
    glcd = Glucose("D", 6, 0, 1.5)
    ps = DiacylPS((0, 0), (18, 0, 0, 0), (20, 4, 0, 0), 7.8)
    psi1 = DiacylPS((3, 3), (18, 0, 0, 0), (20, 4, 0, 0), 7.8)
    psi2 = DiacylPS((0, 0), (18, 0, 5, 0), (20, 4, 0, 0), 7.8)
    # mw, mz, formula
    # ionvariants, chemicalvariants
end
