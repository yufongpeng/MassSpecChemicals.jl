using MassSpecChemicals
using Statistics, StatsBase, TypedTables
const MSC = MassSpecChemicals
import MassSpecChemicals: AbstractChemical, getchemicalattr, adductelements, adductisotopes, adductformula, kmer, charge
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

getchemicalattr(m::T, ::Val{:name}; kwargs...) where {T <: Hexose} = string(m.chirality, "-", string(T), repr_isotope(m))
function getchemicalattr(m::Hexose, ::Val{:elements}; kwargs...) 
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

getchemicalattr(m::Hexose, ::Val{:abbreviation}; kwargs...) = string(first(chemicalname(m), 5), repr_isotope(m))
getchemicalattr(m::Glucose, ::Val{:abbreviation}; kwargs...) = string(first(chemicalname(m), 4), "c", repr_isotope(m))
getchemicalattr(m::Hexose, ::Val{:SMILES}; kwargs...) = ""

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
    rt::Float64
end
repr_headgroup(m::DiacylPS) = string("PS", repr_isotope(m.headgroup...))
repr_fa(fa::FattyAcid) = string(fa.ncb, ":", fa.ndb, repr_isotope(fa.nD, fa.n13C))
getchemicalattr(m::FattyAcid, ::Val{:name}; kwargs...) = string("FA", " ", repr_fa(m))
getchemicalattr(m::DiacylPS, ::Val{:name}; kwargs...) = string(repr_headgroup(m), " ", repr_fa(m.fa1), "/", repr_fa(m.fa2))
function getchemicalattr(m::FattyAcid, ::Val{:elements}; kwargs...)
    filter(x -> last(x) > 0, [
        "C"     => getchemicalattr(m, :ncb) - getchemicalattr(m, :n13C), 
        "[13C]" => getchemicalattr(m, :n13C),
        "H"     => getchemicalattr(m, :ncb) * 2 - getchemicalattr(m, :ndb) * 2 - getchemicalattr(m, :nD),
        "D"     => getchemicalattr(m, :nD),
        "O"     => 2
    ])
end
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
function getchemicalattr(m::DiacylPS, ::Val{:elements}; kwargs...)
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
getchemicalattr(m::DiacylPS, ::Val{:formula}; unique = false, kwargs...) = chemicalformula(chemicalelements(m); unique, kwargs...)
# getchemicalattr(m::DiacylPS, ::Val{:formula}) = 
#     string("C", 6 + getchemicalattr(m, :ncb) - getchemicalattr(m, :n13C), 
#         getchemicalattr(m, :n13C) > 0 ? string("[13C]", getchemicalattr(m, :n13C)) : "",
#         "H", 12 + getchemicalattr(m, :ncb) * 2 - getchemicalattr(m, :ndb) * 2 - getchemicalattr(m, :nD),
#         getchemicalattr(m, :nD) > 0 ? string("D", getchemicalattr(m, :nD)) : "",
#         "N", "O", 10, "P"
#     )
getchemicalattr(m::FattyAcid, ::Val{:SMILES}; kwargs...) = ""
getchemicalattr(m::DiacylPS, ::Val{:SMILES}; kwargs...) = ""
getchemicalattr(m::FattyAcid, ::Val{:nD}; kwargs...) = m.nD
getchemicalattr(m::DiacylPS, ::Val{:nD}; kwargs...) = +(m.headgroup[1], getchemicalattr(m.fa1, :nD), getchemicalattr(m.fa2, :nD))
getchemicalattr(m::FattyAcid, ::Val{:n13C}; kwargs...) = m.n13C
getchemicalattr(m::DiacylPS, ::Val{:n13C}; kwargs...) = +(m.headgroup[2], getchemicalattr(m.fa1, :n13C), getchemicalattr(m.fa2, :n13C))
getchemicalattr(m::FattyAcid, ::Val{:ncb}; kwargs...) = m.ncb
getchemicalattr(m::DiacylPS, ::Val{:ncb}; kwargs...) = getchemicalattr(m.fa1, :ncb) + getchemicalattr(m.fa2, :ncb)
getchemicalattr(m::FattyAcid, ::Val{:ndb}; kwargs...) = m.ndb
getchemicalattr(m::DiacylPS, ::Val{:ndb}; kwargs...) = getchemicalattr(m.fa1, :ndb) + getchemicalattr(m.fa2, :ndb)

# Interface AbstractAdduct
struct DeSerine <: AbstractNegAdduct end
set_adduct_name!("[M-Ser]-", DeSerine())
adductelements(::DeSerine) = ["C" => -3, "H" => -6, "N" => -1, "O" => -2]
adductformula(::DeSerine) = "-Ser"
kmer(::DeSerine) = 1
charge(::DeSerine) = -1

struct Halfprotonation <: AbstractPosAdduct end
set_adduct_name!("[2M+H]+", Halfprotonation())
adductelements(::Halfprotonation) = ["H" => 1]
adductformula(::Halfprotonation) = "+H"
kmer(::Halfprotonation) = 2
charge(::Halfprotonation) = 1

# Interface AbstractIon
adductisotopes(ion::AdductIon{DiacylPS, DeSerine}) = ["H" => ioncore(ion).headgroup[1], "D" => -ioncore(ion).headgroup[1], "C" => ioncore(ion).headgroup[2], "[13C]" => -ioncore(ion).headgroup[2]]

test_show(x) = show(IOBuffer(), x)
macro test_noerror(x)
    return quote
        try 
            $x
            true
        catch e
            false
        end
    end
end

@testset "MassSpecChemicals.jl" begin
    @testset "Generic interface" begin
        # Default chemical, adduction
        global cglc = Chemical("Glucose", ["C" => 6, "H" => 12, "O" => 6]; rt = 1.5, abbreviation = "Glc", SMILES = "")
        global cgld = parse_chemical("Glucose-d6"; formula = "C6H6D6O6", rt = 1.5, abbreviation = "Glc[D6]", SMILES = "")
        global cps = Chemical("PS 18:0/20:4(5Z,8Z,11Z,14Z)", "C44H78NO10P"; rt = 7.8)
        global cpsi1 = Chemical("PS[D3,13C3] 18:0/20:4(5Z,8Z,11Z,14Z)", "C41[13C]3H75D3NO10P"; rt = 7.8)
        global cpsi2 = Chemical("PS 18:0[D5]/20:4(5Z,8Z,11Z,14Z)", "C44H73D5NO10P"; rt = 7.8)
        global lossserine = NegAdduct(1, "-C3H6NO2", 1)
        global lossserinei = NegAdduct(1, "-[13C]3H3D3NO2", 1)
        global dimh = PosAdduct(2, "+H", 1)
        # test all default adduct 
        global icglcall = [AdductIon(cglc, k) for k in keys(MSC.ADDUCT_NAME)]
        @test @test_noerror test_show(cglc)
        @test @test_noerror [test_show(k) for k in icglcall]
        @test @test_noerror [chemicalformula(k) for k in icglcall]
        global icglc = [AdductIon(cglc, Protonation())]
        global icgld = [AdductIon("Glucose-d6", "[M+H]+"; formula = "C6H6D6O6", rt = 1.5), AdductIon("Glucose-d6", Protonation(); formula = "C6H6D6O6", rt = 1.5), AdductIon(cgld, Protonation())]
        global icps = [AdductIon(cps, lossserine), AdductIon(cps, dimh)]
        global icpsi1 = [AdductIon(cpsi1, lossserinei), AdductIon(cpsi1, dimh)]
        global icpsi2 = [AdductIon(cpsi2, lossserine), AdductIon(cpsi2, dimh)]
        global cp1 = ChemicalPair(icps[1], AdductIon(Chemical("FA 20:4", "C20H32O2"; rt = 7.78), Deprotonation()))
        # name, formula, elements
        @test chemicalname(cglc) == "Glucose"
        @test chemicalname(icgld[1]) == "[Glucose-d6+H]+"
        @test chemicalformula(cglc) == "C6H12O6"
        @test chemicalformula(icgld[1]) == "C6H7D6O6"
        @test chemicalelements(icpsi2[2]) == vcat([a => b * 2 for (a, b) in chemicalelements(cpsi2)], adductelements(icpsi2[2]))
        # mmi, mz 
        @test isnan(mz(cglc))
        @test isapprox(mz(icpsi1[1]), mz(icps[1]))
        @test isapprox(mmi(ioncore(icps[2])) * 2 + mmi("H") - MSC.ustrip(MSC.ME), mz(icps[2]))
        # other attr
        @test chemicalabbr(cglc) == "Glc"
        @test chemicalabbr(icgld[3]) == "[Glc[D6]+H]+"
        @test chemicalsmiles(cgld) == ""
        @test ischemicalequal(ioncore(icpsi1[1]), cpsi1)
        @test ionadduct(icpsi2[1]) == lossserine
        @test isadductequal(PosAdduct(2, "+Na+H-H2O", 1), PosAdduct(2, "-H2O+H+Na", 1))
        @test isadductequal(NegAdduct(2, "-H-H2O", 1), NegAdduct(2, "-H2O-H", 1))
        @test !isadductequal(PosAdduct(2, "+Na+H-H2O", 1), NegAdduct(2, "-H2O-H", 1))
        @test kmer(icgld[1]) == 1
        @test charge(icgld[1]) == 1
        @test ncharge(cglc) == 0
        @test ischemicalequal(abundantchemical(icgld[1]), abundantchemical(icgld[2]))
        @test ischemicalequal(cp1, cp1)
        @test isapprox(rt(icpsi2[1]), 7.8) 
    end
    @testset "Customized interface" begin 
        # User defined chemical, adduct
        global glc = Glucose("D", 0, 0, 1.5)
        global gld = Glucose("D", 6, 0, 1.5)
        global ps = DiacylPS((0, 0), FattyAcid(18, 0, 0, 0), FattyAcid(20, 4, 0, 0), 7.8)
        global psi1 = DiacylPS((3, 3), FattyAcid(18, 0, 0, 0), FattyAcid(20, 4, 0, 0), 7.8)
        global psi2 = DiacylPS((0, 0), FattyAcid(18, 0, 5, 0), FattyAcid(20, 4, 0, 0), 7.8)
        global iglc = [AdductIon(glc, Protonation())]
        global igld = [AdductIon(gld, "[M+H]+")]
        global ips = [AdductIon(ps, "[M-Ser]-"), AdductIon(ps, "[2M+H]+")]
        global ipsi1 = [AdductIon(psi1, "[M-Ser]-"), AdductIon(psi1, "[2M+H]+")]
        global ipsi2 = [AdductIon(psi2, "[M-Ser]-"), AdductIon(psi2, "[2M+H]+")]
        global clossserine = ChemicalLoss(Chemical("Serine", "C3H5NO2"))
        # name, formula, elements 
        @test chemicalname(ps) == "PS 18:0/20:4"
        @test chemicalname(igld[1]) == "[D-Glucose[D6]+H]+"
        @test chemicalname(ips[1]) == "[(PS 18:0/20:4)-Ser]-"
        @test chemicalformula(ps) != chemicalformula(cps) 
        @test chemicalformula(ipsi1[2]; unique = true) == chemicalformula(icpsi1[2]) 
        @test MSC.unique_elements(chemicalelements(ipsi1[1])) ==  MSC.unique_elements(vcat(chemicalelements(ioncore(ipsi1[1])), adductelements(ipsi1[1]))) 
        # mmi, mz
        @test isapprox(mz(glc, Protonation()), mz(iglc[1]))
        @test isapprox(mz(ips[2]), mz(ips[2], dimh))
        @test isapprox(mz(ipsi1[1]), mmi(icpsi1[1]))
        # other attr
        @test chemicalabbr(igld[1]) == "[D-Glc[D6]+H]+"
        @test kmer(ipsi1[2]) == 2
        @test charge(ipsi1[2]) == 1
        @test ischemicalequal(abundantchemical(ioncore(ipsi1[1])), abundantchemical(ioncore(ipsi1[2])))
        @test isapprox(rt(ips[1]), 7.8)
    end
    @testset "Isotopologues" begin 
        # isotopologues
        global it1 = isotopologues_table(icglc[1], 1e5; threshold = crit(1e1, 1e-2))
        global it2 = isotopologues_table(ioncore(icglc[1]), 1e5; abtype = :total, threshold = crit(1e-2, 1e-2))
        it4 = isotopologues_table("C6H12O6", 1e5; abtype = :total, threshold = crit(1e-2, 1e-2))
        @test all(>(1e2), mapreduce(x -> x.abundance, vcat, it1.Isotopologues))
        @test all(ischemicalequal.(isotopologues(icglc[1], 1e5; threshold = crit(1e1, 1e-2)), it1.Isotopologues))
        # isotopologues MS/MS 
        itit1 = isotopologues_table(cp1, 1e5; threshold = crit(1e1, 1e-2))
        itit2 = isotopologues_table(chemicalformula(icps[1]) => chemicalformula(cp1.product), 1e5; threshold = crit(1e1, 1e-2))
        itit4 = isotopologues_table(chemicalformula(icps[1]) => chemicalformula(cp1.product), 1e5; threshold = crit(1e1, 1e-2), isobaric = false)
        itit5 = isotopologues_table(chemicalformula(AdductIon(cps, "[M-H]-")) => "-C3H5NO2", 1e5; threshold = crit(1e1, 1e-2), isobaric = false, net_charge = (-1, 0))
        itit7 = isotopologues_table("S8" => "S7"; abtype = :total, isobaric = false, net_charge = (2, 1))
        itit8 = isotopologues_table("C2H5[13C]OO" => "C2H5"; abtype = :total, isobaric = false, net_charge = (-1, -1))
        itit9 = isotopologues_table("C2H5[13C]OO" => "-[13C]O2"; abtype = :total, isobaric = false, net_charge = (-1, 0))
        itit10 = isotopologues_table("C2H5[12C]OO" => "-[12C]O2"; abtype = :total, isobaric = false, net_charge = (-1, 0))
        @test all(>(1e1), mapreduce(x -> x.abundance, vcat, itit1.Isotopologues))
        @test all(ischemicalequal.(isotopologues(cp1, 1e5; threshold = crit(1e1, 1e-2)), itit1.Isotopologues))
        @test isapprox(itit1.Abundance[2], itit2.Abundance[2])
        @test isapprox(itit7.Abundance[10],
            ELEMENTS[:ABUNDANCE]["S"] ^ 6 * ELEMENTS[:ABUNDANCE]["[33S]"] * ELEMENTS[:ABUNDANCE]["[34S]"] * 
            factorial(7, 5) 
        )
        @test isapprox(itit8.Abundance[2],
            ELEMENTS[:ABUNDANCE]["[13C]"] ^ 1 * ELEMENTS[:ABUNDANCE]["C"] ^ 1 * ELEMENTS[:ABUNDANCE]["H"] ^ 5 * ELEMENTS[:ABUNDANCE]["O"] ^ 2 * 2
        )
        @test isapprox(itit8.Abundance[3], itit9.Abundance[3])
        @test isapprox(itit10.Abundance[4], itit9.Abundance[4])
        # isotopologues
        global it3 = isotopologues_table(ipsi2[1], 1; abtype = :total, threshold = crit(1e-3, 1e-3), isobaric = false)
        d = MSC.unique_elements(chemicalelements(ipsi2[1]))
        @test isapprox(isotopicabundance(glc), isotopicabundance(MSC.unique_elements(chemicalelements(glc))))
        @test isapprox(it3.Abundance[6], 
            factorial(d["C"], d["C"] - 2) / factorial(2) * ELEMENTS[:ABUNDANCE]["C"] ^ (d["C"] - 2)* ELEMENTS[:ABUNDANCE]["[13C]"] ^ 2 * 
            ELEMENTS[:ABUNDANCE]["H"] ^ d["H"] * 
            ELEMENTS[:ABUNDANCE]["N"] ^ get(d, "N", 0) * 
            ELEMENTS[:ABUNDANCE]["O"] ^ d["O"] * 
            ELEMENTS[:ABUNDANCE]["P"] ^ d["P"] 
            )
        @test ischemicalequal(ipsi2[1], it3.Isotopologues[1])
        # isotopologues MS/MS
        itit3 = isotopologues_table(ChemicalPair(ips[1], AdductIon(ioncore(ips[1]).fa1, "[M-H]-")), 1; abtype = :total, threshold = crit(1e-8, 1e-8), isobaric = false)
        itit6 = isotopologues_table(ChemicalPair(AdductIon(ps, "[M-H]-"), clossserine), 1; abtype = :total, threshold = crit(1e-8, 1e-8))
        d1 = MSC.unique_elements(chemicalelements(itit3.Isotopologues[14].precursor))
        d2 = MSC.unique_elements(chemicalelements(itit3.Isotopologues[14].product))
        @test isapprox(itit3.Abundance[1], isotopicabundance(ips[1]))
        @test isapprox(itit3.Abundance[14], 
            ELEMENTS[:ABUNDANCE]["C"] ^ (d1["C"]) * ELEMENTS[:ABUNDANCE]["[13C]"] ^ d1["[13C]"] * 
            ELEMENTS[:ABUNDANCE]["H"] ^ d1["H"] * 
            ELEMENTS[:ABUNDANCE]["N"] ^ get(d1, "N", 0) * 
            ELEMENTS[:ABUNDANCE]["O"] ^ (d1["O"]) * ELEMENTS[:ABUNDANCE]["[17O]"] ^ d1["[17O]"] * 
            ELEMENTS[:ABUNDANCE]["P"] ^ d1["P"] * 
            factorial(d2["C"] + d2["[13C]"], d2["C"]) / factorial(d2["[13C]"]) * 
            factorial(d1["O"] + d1["[17O]"] - d2["O"] - get(d2, "[17O]", 0), d1["O"] - d2["O"]) / factorial(d1["[17O]"] - get(d2, "[17O]", 0))
        )
    end
    @testset "Isobars" begin
        # isobars
        # name, formula, elements
        @test @test_noerror test_show(it1.Isotopologues[2])
        @test all(chemicalformula(it1.Isotopologues[2]) .== ["C5H13O6[13C]", "C6H13O5[17O]", "C6H12O6D"])
        @test chemicalname(it1.Isotopologues[1]) == "Isobars[[Glucose+H]+]"
        @test chemicalname(it2.Isotopologues[3]) == "Isobars[Glucose[18O], Glucose[13C2]]"
        @test MSC.unique_elements(chemicalelements(it1.Isotopologues[1])[1]) == MSC.unique_elements(chemicalelements(icglc[1]))
        # mmi, mz
        @test all(isapprox.(mmi.(isotopologues(ioncore(icglc[1]), 1; abtype = :total, threshold = crit(1e-2, 1e-2))), it2.Mass))
        @test all(isapprox.(mz.(it1.Isotopologues), it1.MZ))
        @test isapprox(mean(mz.(it1.Isotopologues[2].chemicals), weights(it1.Isotopologues[2].abundance)), mz(it1.Isotopologues[2], "[M+H]+"))
        @test isapprox(molarmass(cglc), 6 * molarmass("C") + 12 * molarmass("H") + 6 * molarmass("O"))
        # other attr
        @test chemicalabbr(it2.Isotopologues[1]) == "Isobars[Glc]"
        @test chemicalsmiles(it2.Isotopologues[1]) == "Isobars[]"
        @test isapprox(charge(it1.Isotopologues[2]), 1)
        @test getchemicalattr(it1.Isotopologues[2], :chemicals) == it1.Isotopologues[2].chemicals
        @test getchemicalattr(it1.Isotopologues[2], :abundance) == it1.Isotopologues[2].abundance
        @test ischemicalequal(abundantchemical(it1.Isotopologues[2]), first(it1.Isotopologues[2].chemicals))
        @test isapprox(rt(it2.Isotopologues[1]), 1.5)
    end
    @testset "Isotopomers" begin 
        # isotopomers
        # name, formula, elements 
        @test @test_noerror test_show(it3.Isotopologues[2])
        @test chemicalname(it3.Isotopologues[2]) == "[(PS 18:0[D5]/20:4)-Ser]-[13C]"
        # mmi, mz
        # other attrs
        @test ischemicalequal(abundantchemical(igld[1]), igld[1])
        @test ischemicalequal(abundantchemical(glc), glc)
        @test ischemicalequal(getchemicalattr(abundantchemical(it3.Isotopologues[1]), :parent), ipsi2[1])
        @test MSC.unique_elements(getchemicalattr(it3.Isotopologues[2], :isotopes)) ==  MSC.unique_elements(["[13C]" => 1])
    end
    @testset "ChemicalPair and ChemicalLoss" begin 
        # chemicalpair
        @test @test_noerror test_show(cp1)
        @test chemicalname(cp1) == string(chemicalname(getchemicalattr(cp1, :precursor)), " -> ", chemicalname(getchemicalattr(cp1, :product)))
        @test chemicalabbr(cp1) == string(chemicalabbr(getchemicalattr(cp1, :precursor)), " -> ", chemicalabbr(getchemicalattr(cp1, :product)))
        @test chemicalformula(cp1) == Pair(chemicalformula.([chemicalelements(cp1)...])...)
        @test isapprox(rt(cp1), 7.78)
        @test charge(cp1) == (-1 => -1)
        # chemicalloss
        @test @test_noerror test_show(clossserine)
        @test chemicalname(clossserine) == "Loss_Serine"
        @test chemicalabbr(clossserine) == "Loss_Serine"
        @test chemicalformula(clossserine) == chemicalformula(chemicalelements(clossserine))
    end
    @testset "RT" begin 
        exp1 = Table(;
            AdductIon = AdductIon.([Chemical("C1", "C6H12O6"; rt = 1.5), Chemical("C2", "C6H12O6"; rt = 2), Chemical("C3", "C6H13O5N"; rt = 2.5)], "[M+H]+"),
            FWHM = [0.2, 0.21, 0.25]
        )
        exp2 = isotopologues_table(AdductIon(Chemical("B1", "C6H12O7"; rt = 1.5), "[M-H2O+H]+"))
        exp2 = Table(exp2; FWHM = [(1 + log(1e6, 1)) / 5 for x in exp2.Abundance])
        lib = Table(;
            AdductIon = AdductIon.([
                Chemical("Fructose", "C6H12O6"; rt = 1.25),
                Chemical("Glucose", "C6H12O6"; rt = 1.51), 
                Chemical("Galactose", "C6H12O6"; rt = 2.005), 
                Chemical("Glucosamine", "C6H13O5N"; rt = 2.49),
                Chemical("Galactosamine", "C6H13O5N"; rt = 2.91)
            ], "[M+H]+")
        )
        ir1 = isobars_rt(exp1.AdductIon[1], lib.AdductIon)
        ir2 = isobars_rt(exp1.AdductIon[1], lib; libmz = nothing, librt = nothing, libfwhm = nothing)
        irt1 = isobars_rt_table(exp1.AdductIon, lib.AdductIon)
        irt2 = isobars_rt_table(exp1.AdductIon, lib; libmz = nothing, librt = nothing, libfwhm = nothing)
        irt3 = isobars_rt_table(exp2, lib.AdductIon; expchemical = :Isotopologues, expmz = nothing, exprt = nothing)
        irt4 = isobars_rt_table(exp2, lib; libmz = nothing, librt = nothing, libfwhm = nothing, expchemical = :Isotopologues, expmz = nothing, exprt = nothing, mz_tol = acrit(0.1))
        @test ir1 == lib.AdductIon[2:2]
        @test ir2 == ir1
        @test all(splat(isapprox), zip(irt1.ΔRT, sort(rt.(exp1.AdductIon) .- rt.(lib.AdductIon[2:4]))))
        @test irt2.LibID == [2, 3, 4]
        @test isempty(irt3)
        @test irt4.ExpID == [3]
    end
    @testset "New elements" begin 
        m = [23.98504168, 24.985836966, 25.982592972]
        a = [0.78965, 0.1001, 0.11025]
        set_elements!("Mg", m , a)
        @test isapprox(mmi("Mg"), m[1])
        @test isapprox(molarmass("Mg"), m'a)
    end
    @testset "Utils" begin 
        global ct1 = crit(10)
        global ct2 = rcrit(0.2)
        global ct3 = crit(10, 0.2)
        @test @test_noerror test_show(ct1)
        @test @test_noerror test_show(ct2)
        @test @test_noerror test_show(ct3)
        qualified_peak1(x, x̂, ct) = all(c -> in(x, c), makecrit_delta(ct, x̂))
        qualified_peak2(x, x̂, ct) = any(c -> x >= c, makecrit_value(ct, x̂))
        @test !qualified_peak1(85, 100, ct1)
        @test qualified_peak1(85, 100, ct2)
        @test !qualified_peak2(8, 40, ct1)
        @test qualified_peak2(8, 40, ct2)
        @test union(ri"[1,3)", ri"[2,4)") == ri"[1,4)"
        @test 2 in ri"(-Inf, 5]"
        @test ri"(-5, 5]" * (-10) / 2 + 10 - 5 == ri"[-20, 30)"
        @test ri"(-∞, ∞)" * (-10) / 2 + 10 - 5 == ri"(-∞, ∞)"
        @test ri"(-5, 5]" * Inf64 == ri"(-5, 5]" / 0
        @test ri"(-5, 5]" + Inf64 == ri"(-5, 5]" - Inf64
    end
end
