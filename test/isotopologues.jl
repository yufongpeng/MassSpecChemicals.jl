@info "Running tests for Isotopologues and TandemIsotopologues"

@testset "Isotopologues and TandemIsotopologues" begin 
    # isotopologues

    @test all(>(1e1), it1.Abundance1)
    @test isapprox(first(it2.Abundance1), first(first(it4.Abundance1)))
    # @test all(ischemicalequal.(isotopologues(icglc[1], 1e5; threshold = crit(1e1, 1e-2)), it1.Chemical))
    @test all(>(1e1), itit1.Abundance2)
    # @test all(ischemicalequal.(isotopologues(cp1, 1e5; threshold = crit(1e1, 1e-2)), itit1.Chemical))
    @test isapprox(itit1.Abundance2[2], itit2.Abundance2[2])
    @test isapprox(itit3.Abundance2[2], itit4.Abundance2[2])
    @test isapprox(itit7.Abundance2[10],
        MSC.elements_abundance()["S"] ^ 6 * MSC.elements_abundance()["[33S]"] * MSC.elements_abundance()["[34S]"] * 
        factorial(7, 5) 
    )
    @test isapprox(itit8.Abundance2[2],
        MSC.elements_abundance()["[13C]"] ^ 1 * MSC.elements_abundance()["C"] ^ 1 * MSC.elements_abundance()["H"] ^ 5 * MSC.elements_abundance()["O"] ^ 2 * 2
    )
    @test isapprox(itit8.Abundance2[3], itit9.Abundance2[3])
    @test isapprox(itit10.Abundance2[4], itit9.Abundance2[4])
    @test isapprox(itit12.Abundance3[10], itit11.Abundance3[10]*itit12.Abundance3[1])
    # isotopologues
    @test isapprox(isotopicabundance(glc), isotopicabundance(MSC.unique_elements(chemicalelements(glc))))
    @test isapprox(it3.Abundance1[6], 
        factorial(d["C"], d["C"] - 2) / factorial(2) * MSC.elements_abundance()["C"] ^ (d["C"] - 2)* MSC.elements_abundance()["[13C]"] ^ 2 * 
        MSC.elements_abundance()["H"] ^ d["H"] * 
        MSC.elements_abundance()["N"] ^ get(d, "N", 0) * 
        MSC.elements_abundance()["O"] ^ d["O"] * 
        MSC.elements_abundance()["P"] ^ d["P"] 
        )
    @test ischemicalequal(ipsi2[1], it3.Chemical[1])
    # isotopologues MS/MS
    @test isapprox(itit5.Abundance2[1], isotopicabundance(ips[1]))
    @test isapprox(itit5.Abundance2[14], 
        MSC.elements_abundance()["C"] ^ (d1["C"]) * MSC.elements_abundance()["[13C]"] ^ d1["[13C]"] * 
        MSC.elements_abundance()["H"] ^ d1["H"] * 
        MSC.elements_abundance()["N"] ^ get(d1, "N", 0) * 
        MSC.elements_abundance()["O"] ^ (d1["O"]) * MSC.elements_abundance()["[17O]"] ^ d1["[17O]"] * 
        MSC.elements_abundance()["P"] ^ d1["P"] * 
        factorial(d2["C"] + d2["[13C]"], d2["C"]) / factorial(d2["[13C]"]) * 
        factorial(d1["O"] + d1["[17O]"] - d2["O"] - get(d2, "[17O]", 0), d1["O"] - d2["O"]) / factorial(d1["[17O]"] - get(d2, "[17O]", 0))
    )
    @test isapprox(itit17.Abundance2[14], 
        big(MSC.elements_abundance()["C"]) ^ (d1["C"]) * big(MSC.elements_abundance()["[13C]"]) ^ d1["[13C]"] * 
        big(MSC.elements_abundance()["H"]) ^ d1["H"] * 
        big(MSC.elements_abundance()["N"]) ^ get(d1, "N", 0) * 
        big(MSC.elements_abundance()["O"]) ^ (d1["O"]) * big(MSC.elements_abundance()["[17O]"]) ^ d1["[17O]"] * 
        big(MSC.elements_abundance()["P"]) ^ d1["P"] * 
        factorial(big(d2["C"] + d2["[13C]"]), d2["C"]) / factorial(big(d2["[13C]"])) * 
        factorial(big(d1["O"] + d1["[17O]"] - d2["O"] - get(d2, "[17O]", 0)), d1["O"] - d2["O"]) / factorial(big(d1["[17O]"] - get(d2, "[17O]", 0)));
        rtol = 1e-10
    )
    @test isapprox(itl.Abundance1[begin], isotopicabundance(itl.Chemical[begin]); rtol = 1e-6)
end