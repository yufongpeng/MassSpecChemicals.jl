using MassSpecChemicals
using Statistics, StatsBase, TypedTables, DataPipes
const MSC = MassSpecChemicals
import MassSpecChemicals: chemicalname, chemicalelements, chemicalformula, chemicalabbr, chemicalsmiles, completescheme, adductionscheme, charge
using Test

include("generic.jl")
include("customized.jl")
include("large.jl")
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
        # name, formula, elements
        @test chemicalname(cglc) == "Glucose"
        @test chemicalname(icgld[1]) == "[Glucose_d6+H]+"
        @test chemicalformula(cglc) == "C6H12O6"
        @test chemicalformula(icgld[1]) == "C6H7D6O6"
        @test chemicalelements(icpsi2[2]) == vcat([a => b * 2 for (a, b) in chemicalelements(cpsi2)], chemicalelements(ionadduct(icpsi2[2]); loss = false))
        # chemicalloss
        @test @test_noerror test_show(lossserine)
        @test chemicalname(lossserine) == "Loss_Serine"
        @test chemicalabbr(lossserine) == "-Ser"
        @test chemicalformula(lossserine) == string("-", chemicalformula(chemicalelements(lossserine; loss = true)))
        # chemicalgain
        @test @test_noerror test_show(ChemicalGain(Ammonia()))
        @test chemicalname(ChemicalGain(Ammonia())) == "Gain_Ammonia"
        @test chemicalabbr(ChemicalGain(Ammonia())) == "+NH3"
        @test chemicalformula(ChemicalGain(Ammonia())) == string("+", chemicalformula(chemicalelements(ChemicalGain(Ammonia()); loss = false)))
        # mix 
        @test @test_noerror test_show(losshserine)
        @test chemicalname(losshserine) == "Loss_[Serine+H]+"
        @test chemicalabbr(losshserine) == "-[Ser+H]+"
        @test chemicalformula(losshserine) == string("-", chemicalformula(chemicalelements(losshserine; loss = true)))

        # mmi, mz 
        @test isnan(mz(cglc))
        @test isapprox(mz(icpsi1[1]), mz(icps[1]))
        @test isapprox(mmi(ioncore(icps[2])) * 2 + mmi("H") - MSC.ustrip(MSC.ME), mz(icps[2]))
        # other attr
        @test chemicalabbr(cglc) == "Glc"
        @test chemicalabbr(icgld[3]) == "[Glc[D6]+H]+"
        @test chemicalsmiles(cgld) == ""
        @test ischemicalequal(ioncore(icpsi1[1]), cpsi1)
        @test ischemicalequal(icpsi1[1], AdductIon(cpsi1, losshserine))
        @test ischemicalequal(parse_chemical("[+Na+H-H2O]2+"), parse_chemical("[-H2O+Na+H]2+"))
        @test ischemicalequal(parse_chemical("-[2H]2+"), parse_chemical("[-2H]2-"))
        @test ischemicalequal(parse_chemical("+[NH4]+"), parse_chemical("[+NH4]+"))
        @test ncore(icpsi1[2]) == 2
        @test charge(icpsi1[1]) == -1
        @test ncore(icgld[1]) == 1
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
        @test ischemicalequal(parse_chemical(cglc), cglc)
        @test ischemicalequal(parse_chemical("C6H12O6"), FormulaChemical("C6H12O6"))
        @test ischemicalequal(parse_chemical("[C6H12O6]+"), AdductIon(FormulaChemical("C6H12O6"), "[M]+"))
        @test ischemicalequal(parse_chemical("[C6H12O6+H]+"), AdductIon(FormulaChemical("C6H12O6"), ChemicalGain(Proton())))
        @test ischemicalequal(parse_chemical("[C6H12O6+H]+" => "-H2O"), parse_chemical(AdductIon(FormulaChemical("C6H12O6"), ChemicalGain(Proton())), ChemicalLoss(Water())))
        @test ischemicalequal(detectedchemical(parse_chemical("[C6H12O6+H]+" => "-H2O")), parse_chemical("[C6H12O6+H-H2O]+"))
        @test ischemicalequal(parse_chemical(ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1)), "[C6H12O6+2H]2+" => "C6H10O5"), parse_chemical(AdductIon(FormulaChemical("C6H12O6"), "[M+2H]2+") => AdductIon(FormulaChemical("C6H10O5"), "[M]+")))
        @test ischemicalequal(parse_chemical(ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 1)), "[C6H12O6+2H]2+ -> -H3O ->  -[CO2]"), parse_chemical(ChemicalSeries(AdductIon(FormulaChemical("C6H12O6"), "[M+2H]2+"), ChemicalSeries(ChemicalLoss(AdductIon(FormulaChemical("H3O"), "[M]+")), ChemicalLoss(FormulaChemical("CO2"))))))
        @test ischemicalequal(parse_chemical(ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 1)), "[C6H12O6+2H]2+" => "-H3O" => "-[CO2]"), parse_chemical(AdductIon(FormulaChemical("C6H12O6"), "[M+2H]2+") => ChemicalSeries(ChemicalLoss(AdductIon(FormulaChemical("H3O"), "[M]+")), ChemicalLoss(FormulaChemical("CO2")))))
        @test Chemical("Glucose", ["C" => 6, "H" => 12, "O" => 6]; retentiontime = 1.5, abbreviation = "Glc", SMILES = "") == cglc 
        @test hash(Chemical("Glucose", ["C" => 6, "H" => 12, "O" => 6]; retentiontime = 1.5, abbreviation = "Glc", SMILES = "")) == hash(cglc) 
        @test lossserine == ChemicalLoss(cserine)
        @test hash(lossserine) == hash(ChemicalLoss(cserine))
        @test structuralscheme(ionadduct(icpsi1[2])) == ChemicalGain(Proton())
        @test hash(structuralscheme(ionadduct(icpsi1[2]))) == hash(ChemicalGain(Proton()))
        @test last(icglc) == AdductIon(cglc, "[M+H-H2O]+")
        @test hash(last(icglc)) == hash(AdductIon(cglc, "[M+H-H2O]+"))
        @test cp3 == ChemicalSeries(AdductIon(cps, "[M-H]-"), cp1)
        @test hash(cp3) == hash(ChemicalSeries(AdductIon(cps, "[M-H]-"), cp1))
        @test pt1.Chemical[3] == peak_table(MSScan(Isotopologues(icglc[1]; abundance = 1e5, threshold = crit(1e1, 1e-2)))).Chemical[3]
        @test hash(pt1.Chemical[3]) == hash(peak_table(MSScan(Isotopologues(icglc[1]; abundance = 1e5, threshold = crit(1e1, 1e-2)))).Chemical[3])
        @test git3.Chemical[2] == group_isotopologues(Isotopologues(ipsi2[1]; abtype = :total, threshold = crit(1e-3, 1e-3))).Chemical[2]
        @test hash(git3.Chemical[2]) == hash(group_isotopologues(Isotopologues(ipsi2[1]; abtype = :total, threshold = crit(1e-3, 1e-3))).Chemical[2])
    end
    @testset "Customized interface" begin 
        # User defined chemical, adduct
        # name, formula, elements 
        @test chemicalname(ps) == "PS 18:0/20:4"
        @test chemicalname(igld[1]) == "[D-Glucose[D6]+H]+"
        @test chemicalname(ips[1]) == "[(PS 18:0/20:4)-Ser-H]-"
        @test chemicalname(ionadduct(first(ipsi1))) == "Loss_[Serine[D3,13C3]+H]+"
        @test chemicalabbr(ionadduct(first(ipsi1))) == "-[Ser[D3,13C3]+H]+"
        @test chemicalformula(ionadduct(first(ipsi1))) == chemicalformula(chemicalelements(ionadduct(first(ipsi1)); loss = true))
        @test chemicalformula(ps) != chemicalformula(cps) 
        @test chemicalformula(ipsi1[2]; unique = true) == chemicalformula(icpsi1[2]) 
        @test MSC.unique_elements(chemicalelements(ipsi1[1])) ==  MSC.unique_elements(vcat(chemicalelements(ioncore(ipsi1[1])), chemicalelements(ionadduct(ipsi1[1]); loss = false))) 
        @test parse_chemical("[-Ser-H]-") == LossProtonSerine()
        @test parse_chemical("[-Ser-H-H2O]-") == ChemicalSchema(ChemicalLoss(Serine()), ChemicalLoss(Proton()), ChemicalLoss(Water()))
        # mmi, mz
        @test isapprox(mz(glc, "[M+H]+"), mz(iglc[1]))
        @test isapprox(mz(ips[2]), mz(ioncore(ips[2]), "[2M+H]+"))
        @test isapprox(mz(ipsi1[1]), mmi(icpsi1[1]))
        # other attr
        @test chemicalabbr(igld[1]) == "[D-Glc[D6]+H]+"
        @test ncore(ipsi1[2]) == 2
        @test charge(ipsi1[2]) == 1
        @test ischemicalequal(chemicalentity(ioncore(ipsi1[1])), chemicalentity(ioncore(ipsi1[2])))
        @test isapprox(retentiontime(ips[1]), 7.8)
    end
    @testset "ChemicalTransition" begin 
        # chemicaltransition
        @test @test_noerror test_show(cp1)
        @test chemicalname(cp1) == string(chemicalname(inputchemical(cp1)), " -> ", chemicalname(outputchemical(cp1)))
        @test chemicalabbr(cp1) == string(chemicalabbr(inputchemical(cp1)), " -> ", chemicalabbr(outputchemical(cp1)))
        @test chemicalformula(cp1) == chemicalformula(chemicalelements(cp1))
        @test isapprox(retentiontime(cp1), 7.8)
        @test charge(cp1) == -1
        # other attr 
        # chemicalstructure  
        @test !ischemicalequal(ChemicalGain(Proton()), lossserine)
        @test !ischemicalequal(cp1, lossserine)
        @test ischemicalequal(chemicalentity(cp2), chemicalentity(cp3))
        @test ischemicalequal(analyzedchemical(cp2), analyzedchemical(cp3))
        @test ischemicalequal(detectedchemical(cp2), detectedchemical(cp3))
        @test ischemicalequal(inputchemical(cp2), inputchemical(cp3))
        @test ischemicalequal(outputchemical(cp1), outputchemical(cp3))
        @test ischemicalequal(seriesanalyzedchemical(cp2)[2], icps[1])
        @test ischemicalequal(seriesanalyzedchemical(sp2)[2], ipsi1[1])
        @test ischemicalequal(detectedchemical(sp4), AdductIon(psi2.fa1, ChemicalLoss(Proton())))
        @test detectedisotopes(cp5) == last(seriesanalyzedisotopes(cp5))
        @test detectedcharge(cp5) == last(seriesanalyzedcharge(cp5))
        @test MSC.dictionary_elements(seriesanalyzedelements(cp5)[2]) == MSC.dictionary_elements(MSC.loss_elements(chemicalelements(cp5), chemicalelements(chemicaltransition(cp5)[2])))
        @test MSC.dictionary_elements(seriesanalyzedelements(sp5)[2]) == MSC.dictionary_elements(chemicalelements(isotopomerize(completescheme(ipsi1[1], outputchemical(sp1)), ["[13C]" => 5])))
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
        # Tandemisotopologues / Isotopologues MS/MS
        itit1 = TandemIsotopologues(cp1; abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.3)
        itit2 = TandemIsotopologues(chemicalformula(icps[1]) => chemicalformula(last(chemicaltransition(cp1))); abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.3, chemicalparser = ChemicalTransitionParser())
        itit3 = TandemIsotopologues(ChemicalSeries(AdductIon(cps, "[M-H]-"), lossserine); abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.5)
        itit4 = TandemIsotopologues(string("[", chemicalformula(cps), "-H]-") => "-C3H5NO2"; abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.5)
        itit5 = TandemIsotopologues(ChemicalSeries(ips[1], AdductIon(ioncore(ips[1]).fa1, "[M-H]-")); abtype = :total, threshold = crit(1e-8, 1e-8), transmission = 0.3)
        itit6 = TandemIsotopologues(ChemicalSeries(AdductIon(ps, "[M-H]-"), lossserine); abtype = :total, threshold = crit(1e-8, 1e-8), transmission = 0.5)
        itit7 = Isotopologues("[S8]2-" => "[S7]-"; abtype = :total)
        itit8 = Isotopologues("[C2H5[13C]OO]-" => "[C2H5]-"; abtype = :total)
        itit9 = Isotopologues("[C2H5[13C]OO]-" => "-[13C]O2"; abtype = :total)
        itit10 = Isotopologues("[C2H5[12C]OO]-" => "-[12C]O2"; abtype = :total)
        itit11 = Isotopologues("[C2H5[12C]OO]-" => "+H2O"; abtype = :total)
        @test all(>(1e1), itit1.Abundance2)
        # @test all(ischemicalequal.(isotopologues(cp1, 1e5; threshold = crit(1e1, 1e-2)), itit1.Chemical))
        @test isapprox(itit1.Abundance2[2], itit2.Abundance2[2])
        @test isapprox(itit3.Abundance2[2], itit4.Abundance2[2])
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
        @test isapprox(itl.Abundance1[begin], isotopicabundance(itl.Chemical[begin]); rtol = 1e-6)
        # TaandemIsotopologues
    end
    @testset "Spectrum" begin 
        @test @test_noerror test_show(spec1)
        @test isapprox(pt2.Abundance2[1], 20000)
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
        @test !ischemicalequal(pt1.Chemical[1], pt1.Chemical[2])
        @test ischemicalequal(pt1.Chemical[2], pt1.Chemical[2])
        @test ischemicalequal(chemicalentity(pt1.Chemical[1]), pt1.Chemical[1])
        @test ischemicalequal(chemicalentity(pt1.Chemical[2]), first(pt1.Chemical[2].chemicals))
        @test isapprox(retentiontime(pt1.Chemical[1]), 1.5)
    end
    @testset "Isotopomers" begin 
        # isotopomers
        # name, formula, elements 
        @test @test_noerror test_show(it3.Chemical[2])
        @test chemicalname(it3.Chemical[2]) == "[(PS 18:0[D5]/20:4)-Ser-H]-[13C]"
        # mmi, mz
        @test all(isapprox.(mmi.(it3.Chemical), it3.MZ1; rtol = 20e-6))
        @test all(isapprox.(mz.(it3.Chemical), it3.MZ1; rtol = 20e-6))
        # other attrs
        @test ischemicalequal(chemicalentity(igld[1]), igld[1])
        @test ischemicalequal(chemicalentity(glc), glc)
        @test ischemicalequal(it3.Chemical[1], ipsi2[1])
        @test MSC.unique_elements(isotopomersisotopes(it3.Chemical[2])) ==  MSC.unique_elements(["[13C]" => 1])
    end
    @testset "Groupedisotopomers" begin 
        @test @test_noerror test_show(git3.Chemical[2])
        @test chemicalname(git3.Chemical[2]) == "[(PS 18:0[D5]/20:4)-Ser-H]-(+1)"
        @test all(isapprox.(mmi.(git3.Chemical), git3.MZ1; rtol = 20e-6))
        @test all(isapprox.(mz.(git3.Chemical), git3.MZ1; rtol = 20e-6))
        @test ischemicalequal(git3.Chemical[2], git3.Chemical[2])
        @test ischemicalequal(git3.Chemical[1], git3.Chemical[1].parent)
        @test ischemicalequal(MSC.Groupedisotopomers(git3.Chemical[1].parent, 1, "[13C]", [["D" => 1]], [1.0]), Isotopomers(git3.Chemical[1].parent, ["D" => 1]))
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
        @test MSC.lastcolnum(["MZ1", "Abundance2", "MZ3", "MZ2", "Abundance1", "Abundance3"], "MZ") == :MZ3
        @test MSC.ithcolnum(["MZ1", "Abundance2", "MZ3", "MZ2", "Abundance1", "Abundance3"], "MZ", 1) == :MZ1
        @test MSC.allcolnum(["MZ1", "Abundance2", "MZ3", "MZ2", "Abundance1", "Abundance3"], "Abundance") == [:Abundance1, :Abundance2, :Abundance3]
        @test MSC.normalize_abundance([1, 2, 3, 4], 1, :raw) == [1, 2, 3, 4]
        @test MSC.normalize_abundance([1, 2, 3, 4], 1, :total) == [1, 2, 3, 4]
        @test MSC.normalize_abundance([1, 2, 3, 4], 1, :input) == [1, 2, 3, 4] ./ 1
        @test all(isapprox.(MSC.normalize_abundance([1, 2, 3, 4], 1, :max), [1, 2, 3, 4] ./ 4))
        @test all(isapprox.(MSC.normalize_abundance([1, 2, 3, 4], 1, :list), [1, 2, 3, 4] ./ 10))
        @test match_chemical(AdductIon.([Chemical("Fructose", "C6H12O6"), Chemical("Glucose", "C6H12O6"), Chemical("Galactose", "C6H12O6")], "[M+H]+"), exp1).LibID == [1, 2, 3]
    end
end
