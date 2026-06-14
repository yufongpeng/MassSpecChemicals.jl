@info "Running tests for Spectrum and CoelutingIsobars"

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
    @test @test_noerror plot_resolving_power!((0, 1000), FTICR())
    @test MSC.fwhm_mz(QIT(), 100.0) == MSC.fwhm_mz(LIT(), 100.0)
end

@testset "CoelutingIsobars" begin 
    @test @test_noerror test_show(ci1)
    @test retentiontime(ci1.target.Chemical[4]) - retentiontime(ci1.tables[4].Chemical[1]) <= 0.15
    @test isapprox(sum(ci1.tables[4].Abundance1[2:3]), ib1.var"Abundance1(%)"[4] * ci1.target.Abundance1[4] / 100)
    @test all(x -> abs(x) < 0.1, ib2.var"ΔMZ1")
    @test all(>(100 * 1e-4), ib2.var"Abundance2(%)")
end
