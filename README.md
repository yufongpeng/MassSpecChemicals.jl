# MassSpecChemicals

[![Build Status](https://github.com/yufongpeng/MassSpecChemicals.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/yufongpeng/MassSpecChemicals.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/yufongpeng/MassSpecChemicals.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yufongpeng/MassSpecChemicals.jl)

A package for representing molecules or ions formed in mass spectrometers. 

# Basic type
1. `Chemical`: unstructured chemicals, storing name, formula, and other information; `Chemical(name::String, formula::String, info::Vector{Pair{Symbol, Any}})`.
2. `Ion`: charged `Chemical` with specific adduct or molecule loss; `Ion(core::AbstractChemical, adduct::AbstractAdduct)`.
3. `IonCluster`: multiple `Ion`s with similar m/z; `IonCluster(ions::Vector{<: AbstractIon}, abundance::Vector{Float64})`.

# Formula
Currently supported elements:

C, H, O, N, P, S, Li, Na, K, F, Cl, Ag, [13C], [2H] (D), [3H] (T), [17O], [18O], [15N], [33S], [34S], [35S], [6Li], [40K], [41K], [37Cl], [109Ag]

# Adduct
Prefefined adducts:

* `LossElectron`: [M]+
* `Protonation`: [M+H]+
* `ProtonationNLH2O`: [M+H-H2O]+
* `ProtonationNL2H2O`: [M+H-2H2O]+
* `ProtonationNL3H2O`: [M+H-3H2O]+
* `DiProtonation`: [M+2H]2+
* `TriProtonation`: [M+3H]3+
* `AddNH4`: [M+NH4]+
* `AddHNH4`: [M+H+NH4]2+
* `Add2NH4`: [M+2NH4]2+
* `Sodization`: [M+Na]+
* `SodizationProtonation`: [M+Na+H]2+
* `DiSodization`: [M+2Na]2+
* `AddElectron`: [M]-
* `Deprotonation`: [M-H]-
* `DeprotonationNLH2O`: [M-H-H2O]-
* `DiDeprotonation`: [M-2H]2-
* `TriDeprotonation`: [M-3H]3-
* `AddOAc`: [M+CH3COO]-
* `AddHCOO`: [M+HCOO]-
* `LossCH2O`: [M-CH2O]-
* `AddO`: [M+O]-
* `AddC2H2O`: [M+C2H2O]-
* `LossCH8NO`: [M-CH8NO]-
* `LossC2H8NO`: [M-C2H8NO]-
* `AddC3H5NO`: [M+C3H5NO]-
* `AddC2H5NO`: [M+C2H5NO]-
* `LossCH3`: [M-CH3]-
* `DeprotonationLossSerineAddH2O`: [M-H-Serine+H2O]-

User can custumize adduct by `PosAdduct` and `NegAdduct`.

For example, [2M+H]+, and [M-2H-2H2O]2- can be created by `PosAdduct(2, "+H", 1)`, and `NegAdduct(1, "-2H-2H2O", 2)` respectively.

# Interface
* `chemicalname`: name of a chemical (`AbstractChemical`).
* `chemicalformula`: formula of a chemical.
* `ischemicalequal`: whether two chemicals equal.
* `ioncore`: core of an ion.
* `ionadduct`: adduct of an ion.
* `kmer`: number of "M" of an adduct or ion.
* `charge`: number of charges of an adduct or ion.
* `adductelement`: changes of elements of an adduct or ion.
* `adductformula`: formula of an adduct.
* `abundantion`: most abundant ion of an `IonCluster`.
* `isionequal`: whether two ions equal.
* `isionequal`: whether two adducts equal.
* `rt`: retention time.

# Other function
* `mw`: molecular weight of a chemical.
* `mz`: m/z of an ion.
* `abundance`: abundance of a chemical among natural isotopes.
* `isotopic_abundance`: formula and abundance of possible natural isotopes of a chemical.
* `getisotopes`: formula of possible natural isotopes of a chemical.
* `isobar_rt_mz1`: coeluting isobars as a table.
* `getisobars`: coeluting isobars as a vector.
* `acrit`: create absolute criterion.
* `rcrit`: create relative criterion.
* `crit`: create both absolute and relative criterion.
* `@ri_str`: real number interval.

# BioChemicals
This submodule defines multiple biology related molecules, including aminoacids, lipids, glycans, etc. (still in active development)