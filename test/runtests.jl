using MassSpecChemicals
using Statistics, StatsBase
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

struct DiacylPS <: AbstractChemical
    headgroup::Tuple{Int, Int}
    fa1::Tuple{Int, Int, Int, Int}
    fa2::Tuple{Int, Int, Int, Int}
    rt::Float64
end
repr_headgroup(m::DiacylPS) = string("PS", repr_isotope(m.headgroup...))
repr_fa(fa::Tuple) = string(fa[1], ":", fa[2], repr_isotope(fa[3], fa[4]))
getchemicalattr(m::DiacylPS, ::Val{:name}; kwargs...) = string(repr_headgroup(m), " ", repr_fa(m.fa1), "/", repr_fa(m.fa2))
function getchemicalattr(m::DiacylPS, ::Val{:elements}; kwargs...)
    [
        "C"     => 6 + getchemicalattr(m, :ncb) - getchemicalattr(m, :n13C), 
        "[13C]" => getchemicalattr(m, :n13C),
        "H"     => 12 + getchemicalattr(m, :ncb) * 2 - getchemicalattr(m, :ndb) * 2 - getchemicalattr(m, :nD),
        "D"     => getchemicalattr(m, :nD),
        "N"     => 1, 
        "O"     => 10, 
        "P"     => 1
    ]
end
# getchemicalattr(m::DiacylPS, ::Val{:formula}) = 
#     string("C", 6 + getchemicalattr(m, :ncb) - getchemicalattr(m, :n13C), 
#         getchemicalattr(m, :n13C) > 0 ? string("[13C]", getchemicalattr(m, :n13C)) : "",
#         "H", 12 + getchemicalattr(m, :ncb) * 2 - getchemicalattr(m, :ndb) * 2 - getchemicalattr(m, :nD),
#         getchemicalattr(m, :nD) > 0 ? string("D", getchemicalattr(m, :nD)) : "",
#         "N", "O", 10, "P"
#     )
getchemicalattr(m::DiacylPS, ::Val{:SMILES}; kwargs...) = ""
getchemicalattr(m::DiacylPS, ::Val{:nD}; kwargs...) = +(m.headgroup[1], m.fa1[3], m.fa2[3])
getchemicalattr(m::DiacylPS, ::Val{:n13C}; kwargs...) = +(m.headgroup[2], m.fa1[4], m.fa2[4])
getchemicalattr(m::DiacylPS, ::Val{:ncb}; kwargs...) = +(m.fa1[1], m.fa2[1])
getchemicalattr(m::DiacylPS, ::Val{:ndb}; kwargs...) = +(m.fa1[2], m.fa2[2])

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
    # Default chemical, adduction
    cglc = Chemical("Glucose", "C6H12O6"; rt = 1.5, abbreviation = "Glc", SMILES = "")
    cgld = parse_chemical("Glucose-d6"; formula = "C6H6D6O6", rt = 1.5, abbreviation = "Glc[D6]", SMILES = "")
    cps = Chemical("PS 18:0/20:4(5Z,8Z,11Z,14Z)", "C44H80NO10P"; rt = 7.8)
    cpsi1 = Chemical("PS[D3,13C3] 18:0/20:4(5Z,8Z,11Z,14Z)", "C41[13C]3H77D3NO10P"; rt = 7.8)
    cpsi2 = Chemical("PS 18:0[D5]/20:4(5Z,8Z,11Z,14Z)", "C44H75D5NO10P"; rt = 7.8)
    lossserine = NegAdduct(1, "-C3H6NO2", 1)
    lossserinei = NegAdduct(1, "-[13C]3H3D3NO2", 1)
    dimh = PosAdduct(2, "+H", 1)
    # test all default adduct 
    global icglcall = [AdductIon(cglc, k) for k in keys(MSC.ADDUCT_NAME)]
    @test @test_noerror [test_show(k) for k in icglcall]
    @test @test_noerror [chemicalformula(k) for k in icglcall]
    icglc = [AdductIon(cglc, Protonation())]
    icgld = [AdductIon("Glucose-d6", "[M+H]+"; formula = "C6H6D6O6", rt = 1.5), AdductIon("Glucose-d6", Protonation(); formula = "C6H6D6O6", rt = 1.5), AdductIon(cgld, Protonation())]
    icps = [AdductIon(cps, lossserine), AdductIon(cps, dimh)]
    icpsi1 = [AdductIon(cpsi1, lossserinei), AdductIon(cpsi1, dimh)]
    icpsi2 = [AdductIon(cpsi2, lossserine), AdductIon(cpsi2, dimh)]
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
    @test kmer(icgld[1]) == 1
    @test charge(icgld[1]) == 1
    @test ncharge(cglc) == 0
    @test ischemicalequal(abundantchemical(icgld[1]), abundantchemical(icgld[2]))
    @test isapprox(rt(icpsi2[1]), 7.8) 
    # isotopologues
    it1 = isotopologues_table(icglc[1], 1e5; threshold = crit(1e1, 1e-2))
    it2 = isotopologues_table(ioncore(icglc[1]), 1; abtype = :all, threshold = crit(1e-2, 1e-2))
    @test all(>(1e2), mapreduce(x -> x.abundance, vcat, it1.Isotopologues))
    @test all(ischemicalequal.(isotopologues(icglc[1], 1e5; threshold = crit(1e1, 1e-2)), it1.Isotopologues))
    # isobars
    # name, formula, elements
    @test all(chemicalformula(it1.Isotopologues[2]) .== ["C5H13O6[13C]", "C6H13O5[17O]", "C6H12O6D"])
    @test chemicalname(it1.Isotopologues[1]) == "Isobars[[Glucose+H]+]"
    @test chemicalname(it2.Isotopologues[3]) == "Isobars[Glucose[18O], Glucose[13C2]]"
    @test MSC.unique_elements(chemicalelements(it1.Isotopologues[1])[1]) == MSC.unique_elements(chemicalelements(icglc[1]))
    # mmi, mz
    @test all(isapprox.(mmi.(isotopologues(ioncore(icglc[1]), 1; abtype = :all, threshold = crit(1e-2, 1e-2))), it2.Mass))
    @test all(isapprox.(mz.(it1.Isotopologues), it1.MZ))
    @test isapprox(mean(mz.(it1.Isotopologues[2].chemicals), weights(it1.Isotopologues[2].abundance)), mz(it1.Isotopologues[2], "[M+H]+"))
    # other attr
    @test chemicalabbr(it2.Isotopologues[1]) == "Isobars[Glc]"
    @test chemicalsmiles(it2.Isotopologues[1]) == "Isobars[]"
    @test isapprox(charge(it1.Isotopologues[2]), 1)
    @test getchemicalattr(it1.Isotopologues[2], :chemicals) == it1.Isotopologues[2].chemicals
    @test getchemicalattr(it1.Isotopologues[2], :abundance) == it1.Isotopologues[2].abundance
    @test ischemicalequal(abundantchemical(it1.Isotopologues[2]), first(it1.Isotopologues[2].chemicals))
    @test isapprox(rt(it2.Isotopologues[1]), 1.5)
    # User defined chemical, adduct
    glc = Glucose("D", 0, 0, 1.5)
    gld = Glucose("D", 6, 0, 1.5)
    ps = DiacylPS((0, 0), (18, 0, 0, 0), (20, 4, 0, 0), 7.8)
    psi1 = DiacylPS((3, 3), (18, 0, 0, 0), (20, 4, 0, 0), 7.8)
    psi2 = DiacylPS((0, 0), (18, 0, 5, 0), (20, 4, 0, 0), 7.8)
    iglc = [AdductIon(glc, Protonation())]
    igld = [AdductIon(gld, "[M+H]+")]
    ips = [AdductIon(ps, "[M-Ser]-"), AdductIon(ps, "[2M+H]+")]
    ipsi1 = [AdductIon(psi1, "[M-Ser]-"), AdductIon(psi1, "[2M+H]+")]
    ipsi2 = [AdductIon(psi2, "[M-Ser]-"), AdductIon(psi2, "[2M+H]+")]
    # name, formula, elements 
    @test chemicalname(ps) == "PS 18:0/20:4"
    @test chemicalname(igld[1]) == "[D-Glucose[D6]+H]+"
    @test chemicalname(ips[1]) == "[(PS 18:0/20:4)-Ser]-"
    @test chemicalformula(ps) == chemicalformula(cps) 
    @test chemicalformula(ipsi1[2]) == chemicalformula(icpsi1[2]) 
    @test  MSC.unique_elements(chemicalelements(ipsi1[1])) ==  MSC.unique_elements(vcat(chemicalelements(ioncore(ipsi1[1])), adductelements(ipsi1[1]))) 
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
    # isotopologues
    it3 = isotopologues_table(ipsi2[1], 1; abtype = :all, threshold = crit(1e-3, 1e-3), isobaric = false)
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
    # isotopomers
    # name, formula, elements 
    @test chemicalname(it3.Isotopologues[2]) == "[(PS 18:0[D5]/20:4)-Ser]-[13C]"
    # mmi, mz
    # other attrs
    @test ischemicalequal(abundantchemical(igld[1]), igld[1])
    @test ischemicalequal(abundantchemical(glc), glc)
    @test ischemicalequal(getchemicalattr(abundantchemical(it3.Isotopologues[1]), :parent), ipsi2[1])
    @test MSC.unique_elements(getchemicalattr(it3.Isotopologues[2], :isotopes)) ==  MSC.unique_elements(["[13C]" => 1])
    @testset "Utils" begin 
        ct1 = crit(10)
        ct2 = rcrit(0.2)
        qualified_peak1(x, x̂, ct) = all(c -> in(x, c), makecrit_delta(ct, x̂))
        qualified_peak2(x, x̂, ct) = any(c -> x >= c, makecrit_value(ct, x̂))
        @test !qualified_peak1(85, 100, ct1)
        @test qualified_peak1(85, 100, ct2)
        @test !qualified_peak2(8, 40, ct1)
        @test qualified_peak2(8, 40, ct2)
        @test union(ri"[1,3)", ri"[2,4)") == ri"[1,4)"
        @test 2 in ri"(-Inf, 5]"
    end
end
