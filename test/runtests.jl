using MassSpecChemicals
using Statistics, StatsBase, TypedTables, DataPipes
const MSC = MassSpecChemicals
import MassSpecChemicals: AbstractChemical, chemicalname, chemicalelements, chemicalformula, chemicalabbr, chemicalsmiles, adductelements, adductisotopes, adductformula, kmer, charge
using Test

include("generic.jl")
include("customized.jl")
include("gq1.jl")
include("coeluting.jl")
include("utils.jl")

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
        @test @test_noerror test_show(cglc)
        @test @test_noerror [test_show(k) for k in icglcall]
        @test @test_noerror [chemicalformula(k) for k in icglcall]
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
        @test isadductequal(Adduct(2, "+Na+H-H2O", 1), Adduct(2, "-H2O+H+Na", 1))
        @test isadductequal(Adduct(2, "-H-H2O", -1), Adduct(2, "-H2O-H", -1))
        @test !isadductequal(Adduct(2, "+Na+H-H2O", 1), Adduct(2, "-H2O-H", -1))
        @test kmer(ComposedAdduct([Protonation(), AddNH4()])) == 1
        @test charge(ComposedAdduct([Protonation(), AddNH4()])) == 2
        @test adductformula(ComposedAdduct([Protonation(), AddNH4()])) == "+H+NH4"
        @test adductelements(ComposedAdduct([Protonation(), AddNH4()])) == ["H" => 1, "N" => 1, "H" => 4]
        @test kmer(icgld[1]) == 1
        @test charge(icgld[1]) == 1
        @test ncharge(cglc) == 0
        @test ischemicalequal(cp1, cp1)
        @test isapprox(retentiontime(icpsi2[1]), 7.8) 
        # chemicalstructure  
        @test ischemicalequal(chemicalentity(icgld[1]), chemicalentity(icgld[2]))
        @test ischemicalequal(first(chemicalspecies(icgld[1])), icgld[2])
        @test ischemicalequal(first(chemicaltransition(icgld[1])), icgld[2])
        @test ischemicalequal(chemicalparent(icgld[1]), icgld[2])
        @test isempty(isotopomersisotopes(icgld[1]))
        @test ischemicalequal(analyzedchemical(icgld[1]), icgld[2])
        @test ischemicalequal(detectedchemical(icgld[1]), icgld[2])
        @test ischemicalequal(inputchemical(icgld[1]), icgld[2])
        @test ischemicalequal(outputchemical(icgld[1]), icgld[2])
        @test ischemicalequal(first(seriesanalyzedchemical(icgld[1])), first(chemicaltransition(icgld[2])))
        @test isempty(detectedisotopes(icgld[1]))
        @test isempty(last(seriesanalyzedisotopes(icgld[1])))
        @test detectedcharge(icgld[1]) == last(seriesanalyzedcharge(icgld[1]))
        @test detectedelements(icgld[1]) == last(seriesanalyzedelements(icgld[1]))
        # ChemicalSeries
        @test ischemicalequal(ChemicalSeries(cglc), cglc)
        @test ischemicalequal(ChemicalSeries("C6H12O6"), FormulaChemical("C6H12O6"))
        @test ischemicalequal(ChemicalSeries("[C6H12O6]+"), AdductIon(FormulaChemical("C6H12O6"), LossElectron()))
        @test ischemicalequal(ChemicalSeries("[C6H12O6+H]+"), AdductIon(FormulaChemical("C6H12O6"), Protonation()))
        @test ischemicalequal(ChemicalSeries("[C6H12O6+H]+" => "-H2O"), ChemicalSeries(AdductIon(FormulaChemical("C6H12O6"), Protonation()), ChemicalLoss(FormulaChemical("H2O"))))
        @test ischemicalequal(ChemicalSeries("[C6H12O6+2H]2+" => "C6H10O5"; charge = 1), ChemicalSeries(AdductIon(FormulaChemical("C6H12O6"), DiProtonation()) => AdductIon(FormulaChemical("C6H10O5"), LossElectron())))
        @test ischemicalequal(ChemicalSeries("[C6H12O6+2H]2+" => "-H3O" => "-[CO2]"; charge = 1, loss = 1), ChemicalSeries(ChemicalSeries(AdductIon(FormulaChemical("C6H12O6"), DiProtonation()), ChemicalSeries(ChemicalLoss(AdductIon(FormulaChemical("H3O"), LossElectron())), ChemicalLoss(FormulaChemical("CO2"))))))
        @test ischemicalequal(ChemicalSeries("[C6H12O6+2H]2+" => "-H3O" => "-[CO2]"; charge = 1, loss = 1), ChemicalSeries(AdductIon(FormulaChemical("C6H12O6"), DiProtonation()) => ChemicalSeries(ChemicalLoss(AdductIon(FormulaChemical("H3O"), LossElectron())), ChemicalLoss(FormulaChemical("CO2")))))
    end
    @testset "Customized interface" begin 
        # User defined chemical, adduct
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
        @test ischemicalequal(chemicalentity(ioncore(ipsi1[1])), chemicalentity(ioncore(ipsi1[2])))
        @test isapprox(retentiontime(ips[1]), 7.8)
    end
    @testset "Isotopologues" begin 
        # isotopologues
        it1 = Isotopologues(icglc[1]; abundance = 1e5, threshold = crit(1e1, 1e-2))
        it2 = Isotopologues(ioncore(icglc[1]); abundance = 1e5, abtype = :total, threshold = crit(1e1, 1e-2))
        it3 = Isotopologues(ipsi2[1]; abtype = :total, threshold = crit(1e-3, 1e-3))
        it4 = Isotopologues("C6H12O6"; abundance = 1e5, abtype = :total, threshold = crit(1e1, 1e-2))
        @test all(>(1e1), it1.Abundance1)
        @test isapprox(first(it2.Abundance1), first(first(it4.Abundance1)))
        # @test all(ischemicalequal.(isotopologues(icglc[1], 1e5; threshold = crit(1e1, 1e-2)), it1.Chemical))
        # isotopologues MS/MS 
        itit1 = Isotopologues(cp1; abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.3)
        itit2 = Isotopologues(chemicalformula(icps[1]) => chemicalformula(last(chemicaltransition(cp1))); abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.3, charge = -1)
        itit3 = Isotopologues(ChemicalSeries(AdductIon(cps, "[M-H]-"), clossserine); abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.5)
        itit4 = Isotopologues(string("[", chemicalformula(cps), "-H]-") => "-C3H5NO2"; abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.5)
        itit5 = Isotopologues(ChemicalSeries(ips[1], AdductIon(ioncore(ips[1]).fa1, "[M-H]-")); abtype = :total, threshold = crit(1e-8, 1e-8), transmission = 0.3)
        itit6 = Isotopologues(ChemicalSeries(AdductIon(ps, "[M-H]-"), clossserine); abtype = :total, threshold = crit(1e-8, 1e-8), transmission = 0.5)
        itit7 = Isotopologues("[S8]2-" => "[S7]-"; abtype = :total)
        itit8 = Isotopologues("[C2H5[13C]OO]-" => "[C2H5]-"; abtype = :total)
        itit9 = Isotopologues("[C2H5[13C]OO]-" => "-[13C]O2"; abtype = :total)
        itit10 = Isotopologues("[C2H5[12C]OO]-" => "-[12C]O2"; abtype = :total)
        @test all(>(1e1), itit1.Abundance2)
        # @test all(ischemicalequal.(isotopologues(cp1, 1e5; threshold = crit(1e1, 1e-2)), itit1.Chemical))
        @test isapprox(itit1.Abundance2[2], itit2.Abundance2[2])
        @test isapprox(itit7.Abundance2[10],
            MSC.elements_abundunce()["S"] ^ 6 * MSC.elements_abundunce()["[33S]"] * MSC.elements_abundunce()["[34S]"] * 
            factorial(7, 5) 
        )
        @test isapprox(itit8.Abundance2[2],
            MSC.elements_abundunce()["[13C]"] ^ 1 * MSC.elements_abundunce()["C"] ^ 1 * MSC.elements_abundunce()["H"] ^ 5 * MSC.elements_abundunce()["O"] ^ 2 * 2
        )
        @test isapprox(itit8.Abundance2[3], itit9.Abundance2[3])
        @test isapprox(itit10.Abundance2[4], itit9.Abundance2[4])
        # isotopologues
        d = MSC.dictionary_elements(chemicalelements(ipsi2[1]))
        @test isapprox(isotopicabundance(glc), isotopicabundance(MSC.unique_elements(chemicalelements(glc))))
        @test isapprox(it3.Abundance1[6], 
            factorial(d["C"], d["C"] - 2) / factorial(2) * MSC.elements_abundunce()["C"] ^ (d["C"] - 2)* MSC.elements_abundunce()["[13C]"] ^ 2 * 
            MSC.elements_abundunce()["H"] ^ d["H"] * 
            MSC.elements_abundunce()["N"] ^ get(d, "N", 0) * 
            MSC.elements_abundunce()["O"] ^ d["O"] * 
            MSC.elements_abundunce()["P"] ^ d["P"] 
            )
        @test ischemicalequal(ipsi2[1], it3.Chemical[1])
        # isotopologues MS/MS
        d1 = MSC.dictionary_elements(chemicalelements(inputchemical(itit5.Chemical[14])))
        d2 = MSC.dictionary_elements(chemicalelements(outputchemical(itit5.Chemical[14])))
        @test isapprox(itit5.Abundance2[1], isotopicabundance(ips[1]))
        @test isapprox(itit5.Abundance2[14], 
            MSC.elements_abundunce()["C"] ^ (d1["C"]) * MSC.elements_abundunce()["[13C]"] ^ d1["[13C]"] * 
            MSC.elements_abundunce()["H"] ^ d1["H"] * 
            MSC.elements_abundunce()["N"] ^ get(d1, "N", 0) * 
            MSC.elements_abundunce()["O"] ^ (d1["O"]) * MSC.elements_abundunce()["[17O]"] ^ d1["[17O]"] * 
            MSC.elements_abundunce()["P"] ^ d1["P"] * 
            factorial(d2["C"] + d2["[13C]"], d2["C"]) / factorial(d2["[13C]"]) * 
            factorial(d1["O"] + d1["[17O]"] - d2["O"] - get(d2, "[17O]", 0), d1["O"] - d2["O"]) / factorial(d1["[17O]"] - get(d2, "[17O]", 0))
        )
        # TaandemIsotopologues
    end
    @testset "Spectrum" begin 
        @test @test_noerror test_show(spec1)
        @test @test_noerror plot_spectrum((948.1, 951.5), spec2)
        @test @test_noerror plot_spectrum!((948.1, 951.5), spec2; deconvolution = true)
        @test @test_noerror plot_window(MSC.GaussianWindow())
        @test @test_noerror plot_window!(MSC.TukeyWindow(0.2))
        @test @test_noerror plot_window!(MSC.CosineWindow())
        @test @test_noerror plot_window!(MSC.SuperGaussianWindow(4))
        @test @test_noerror plot_resolving_power((0, 1000), TOF())
        @test @test_noerror plot_resolving_power!((0, 1000), Orbitrap())
    end
    @testset "Isobars" begin
        # isobars
        # name, formula, elements
        @test @test_noerror test_show(pt1.Chemical[2])
        @test chemicalformula(pt1.Chemical[2]) == "C5H13O6[13C]"
        @test chemicalname(pt1.Chemical[1]) == "Isobars[[Glucose+H]+]"
        @test chemicalname(pt1.Chemical[3]) == "Isobars[[Glucose+H]+[18O], [Glucose+H]+[13C2], [Glucose+H]+[13C,17O]]"
        @test MSC.unique_elements(chemicalelements(pt1.Chemical[1])) == MSC.unique_elements(chemicalelements(icglc[1]))
        # mmi, mz
        @test all(isapprox.(mmi.(pt1.Chemical), pt1.MZ; rtol = 20e-6))
        @test all(isapprox.(mz.(pt1.Chemical), pt1.MZ; rtol = 20e-6))
        @test isapprox(molarmass(pt1.Chemical[1]), 6 * molarmass("C") + 13 * molarmass("H") + 6 * molarmass("O") - MSC.ustrip(MSC.ME))
        # other attr
        @test chemicalabbr(pt1.Chemical[1]) == "Isobars[[Glc+H]+]"
        @test chemicalsmiles(pt1.Chemical[1]) == ""
        @test isapprox(charge(pt1.Chemical[2]), 1)
        @test getchemicalproperty(pt1.Chemical[2], :chemicals) == pt1.Chemical[2].chemicals
        @test getchemicalproperty(pt1.Chemical[2], :abundance) == pt1.Chemical[2].abundance
        @test ischemicalequal(chemicalentity(pt1.Chemical[2]), first(pt1.Chemical[2].chemicals))
        @test isapprox(retentiontime(pt1.Chemical[1]), 1.5)
    end
    @testset "Isotopomers" begin 
        # isotopomers
        # name, formula, elements 
        @test @test_noerror test_show(it3.Chemical[2])
        @test chemicalname(it3.Chemical[2]) == "[(PS 18:0[D5]/20:4)-Ser]-[13C]"
        # mmi, mz
        @test all(isapprox.(mmi.(it3.Chemical), it3.MZ1; rtol = 20e-6))
        @test all(isapprox.(mz.(it3.Chemical), it3.MZ1; rtol = 20e-6))
        # other attrs
        @test ischemicalequal(chemicalentity(igld[1]), igld[1])
        @test ischemicalequal(chemicalentity(glc), glc)
        @test ischemicalequal(chemicalparent(chemicalentity(it3.Chemical[1])), ipsi2[1])
        @test MSC.unique_elements(isotopomersisotopes(it3.Chemical[2])) ==  MSC.unique_elements(["[13C]" => 1])
    end
    @testset "Groupedisotopomers" begin 
        @test @test_noerror test_show(git3.Chemical[2])
        @test chemicalname(git3.Chemical[2]) == "[(PS 18:0[D5]/20:4)-Ser]-(+1)"
        @test all(isapprox.(mmi.(git3.Chemical), git3.MZ1; rtol = 20e-6))
        @test all(isapprox.(mz.(git3.Chemical), git3.MZ1; rtol = 20e-6))
    end
    @testset "ChemicalTransition, ChemicalLoss, and ChemicalGain" begin 
        # chemicaltransition
        @test @test_noerror test_show(cp1)
        @test chemicalname(cp1) == string(chemicalname(inputchemical(cp1)), " -> ", chemicalname(outputchemical(cp1)))
        @test chemicalabbr(cp1) == string(chemicalabbr(inputchemical(cp1)), " -> ", chemicalabbr(outputchemical(cp1)))
        @test chemicalformula(cp1) == chemicalformula(chemicalelements(cp1))
        @test isapprox(retentiontime(cp1), 7.8)
        @test charge(cp1) == -1
        # chemicalloss
        @test @test_noerror test_show(clossserine)
        @test chemicalname(clossserine) == "Loss_Serine"
        @test chemicalabbr(clossserine) == "Loss_Serine"
        @test chemicalformula(clossserine) == chemicalformula(chemicalelements(clossserine))
        # chemicalgain
        @test @test_noerror test_show(gainno)
        @test chemicalname(gainno) == "+NO"
        @test chemicalabbr(gainno) == "+NO"
        @test chemicalformula(gainno) == chemicalformula(chemicalelements(gainno))
        # other attr 
        # chemicalstructure  
        @test ischemicalequal(chemicalentity(cp2), chemicalentity(cp3))
        @test ischemicalequal(analyzedchemical(cp2), analyzedchemical(cp3))
        @test ischemicalequal(detectedchemical(cp2), detectedchemical(cp3))
        @test ischemicalequal(inputchemical(cp2), inputchemical(cp3))
        @test ischemicalequal(outputchemical(cp1), outputchemical(cp3))
        @test ischemicalequal(seriesanalyzedchemical(cp2)[2], Chemical(string(chemicalname(AdductIon(cps, "[M-H]-")), chemicalname(clossserine)), MSC.loss_elements(chemicalelements(AdductIon(cps, "[M-H]-")), chemicalelements(clossserine))))
        @test detectedisotopes(cp5) == last(seriesanalyzedisotopes(cp5))
        @test detectedcharge(cp5) == last(seriesanalyzedcharge(cp5))
        @test MSC.dictionary_elements(seriesanalyzedelements(cp5)[2]) == MSC.dictionary_elements(MSC.loss_elements(chemicalelements(cp5), chemicalelements(chemicaltransition(cp5)[2])))
        @test MSC.dictionary_elements(seriesanalyzedelements(cp6)[2]) == MSC.dictionary_elements(MSC.gain_elements(chemicalelements(cp6), chemicalelements(chemicaltransition(cp6)[2])))
    end
    @testset "CoelutingIsobars" begin 
        @test @test_noerror test_show(ci1)
        @test retentiontime(ci1.target.Chemical[4]) - retentiontime(ci1.tables[4].Chemical[1]) <= 0.15
        @test isapprox(sum(ci1.tables[4].Abundance1[2:3]), ib1.var"Abundance1(%)"[4] * ci1.target.Abundance1[4] / 100)
        @test length(ib2) == 6
        @test all(x -> abs(x) < 0.1, ib2.var"ΔMZ1")
        @test all(>(100 * 1e-4), ib2.var"Abundance2(%)")
    end
    @testset "New elements" begin 
        m = [23.98504168, 24.985836966, 25.982592972]
        a = [0.78965, 0.1001, 0.11025]
        set_element!("Mg", m , a)
        @test isapprox(mmi("Mg"), m[1])
        @test isapprox(molarmass("Mg"), m'a)
    end
    @testset "Utils" begin 
        @test @test_noerror test_show(ct1)
        @test @test_noerror test_show(ct2)
        @test @test_noerror test_show(ct3)
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
