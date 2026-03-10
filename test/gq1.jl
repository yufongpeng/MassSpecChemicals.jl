# GQ1 
GQ1 = Chemical("GQ1b 18:1;O2/18:0", "C106H182N6O55")
GD3 = Chemical("GD3 18:1;O2/18:0", "C70H125N3O29")
nana = Chemical("N-acetylneuraminic acid", "C11H19NO9")
galnac = Chemical("GalNAc", "C8H15NO6")
f366 = Chemical("GalGalNAc", "C6H12O6C8H13NO5")
f657 = Chemical("NeuAcGalGalNAc", "C11H19NO9C6H10O5C8H13NO5")
f948 = Chemical("NeuAc2GalGalNAc", "C22H36N2O17C6H10O5C8H13NO5")
pt = Table(; 
    Chemical = [
        AdductIon(GQ1, "[M+2H]2+"), 
        AdductIon(GQ1, "[M+H+Na]2+"), 
        AdductIon(GQ1, "[M+2H-H2O]2+"), 
        AdductIon(GD3, "[M+H]+"), 
        AdductIon(f948, "[M+H-H2O]+"), 
        AdductIon(f657, "[M+H-H2O]+"), 
        AdductIon(f366, "[M+H-H2O]+"), 
        AdductIon(galnac, "[M+H-H2O]+"), 
        AdductIon(nana, "[M+H-H2O]+"), 
        AdductIon(nana, "[M+H-2H2O]+")], 
    Product = [
        [
            AdductIon(f948, "[M+H-H2O]+"), 
            AdductIon(f657, "[M+H-H2O]+"), 
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-2H2O]+"),
            AdductIon(GQ1, "[M+2H]2+"), 
            AdductIon(GQ1, "[M+2H-H2O]2+"), 
            AdductIon(GD3, "[M+H]+")
        ],
        [
            AdductIon(f948, "[M+H-H2O]+"), 
            AdductIon(f657, "[M+H-H2O]+"), 
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+"),
            AdductIon(GQ1, "[M+H+Na]2+"), 
            AdductIon(GD3, "[M+H]+")
        ],
        [
            AdductIon(f948, "[M+H-H2O]+"), 
            AdductIon(f657, "[M+H-H2O]+"), 
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+"),
            AdductIon(GQ1, "[M+2H-H2O]2+"), 
            AdductIon(GD3, "[M+H]+")
        ],
        [
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+"),
            AdductIon(GD3, "[M+H]+")
        ],
        [
            AdductIon(f948, "[M+H-H2O]+"), 
            AdductIon(f657, "[M+H-H2O]+"), 
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+")
        ],
        [
            AdductIon(f657, "[M+H-H2O]+"), 
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+")
        ],
        [
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
        ],
        [
            AdductIon(galnac, "[M+H-H2O]+")
        ],
        [
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+")
        ],
        [
            AdductIon(nana, "[M+H-2H2O]+")
        ],
    ], 
    Proportion = [
        [0.3, 0.1, 0.1, 0.05, 0.1, 0.2, 0.3, 0.05, 0.1], 
        [0.35, 0.12, 0.12, 0.07, 0.1, 0.2, 0.2], 
        [0.35, 0.12, 0.2, 0.07, 0.12, 0.3, 0.2], 
        [0.1, 0.4, 0.5], 
        [0.05, 0.15, 0.15, 0.3, 0.1, 0.3], 
        [0.05, 0.1, 0.4, 0.15, 0.4], 
        [0.1, 0.7], 
        [0.2, 0.8], [0.8]
    ])
IGQ1b = AdductIon(GQ1, "[M+2H]2+")
spec1 = @p Isotopologues(IGQ1b; abtype = :list, abundance = 100000) |> MSScan
spec2 = @p spec1 |> TargetIon(Quadrupole(1210.588728; fwhm = 1.3, offset = 0.3)) |> Fragmentation(pt) |> MSScan
spec3 = @p spec2 |> TargetIon(Quadrupole(948.3303; fwhm = 1.3, offset = 0.3)) |> Fragmentation(pt) |> MSScan
