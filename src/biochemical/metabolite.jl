module Metabolites
using Reexport
@reexport using ..BioChemicals
using ..MassSpecChemicals: AbstractChemical
import ..BioChemicals: originalmolecule, leavinggroup, conjugation
export Metabolite, 
       Ethanolamine, 
       Nmethylethanolamine, 
       NNdimethylethanolamine, 
       GABA, 
       Dopamine, 
       Taurine,
       Tauryl, # T/Tau
       Carnitine, 
       Choline, 
       CoA, 
       GlycolicAcid, 
       Glycoyl, # Gc
       LacticAcid, 
       Lactyl, # Lt
       PyruvicAcid, 
       Pyruvyl # Py

abstract type Metabolite <: AbstractChemical end
struct Ethanolamine <: Metabolite end
struct Nmethylethanolamine <: Metabolite end
struct NNdimethylethanolamine <: Metabolite end
struct GABA <: Metabolite end
struct Dopamine <: Metabolite end
struct Taurine <: Metabolite end # T
struct Tauryl <: FunctionalGroup{Taurine, Dehydrogen} end # T/Tau
struct Carnitine <: Metabolite end
struct Choline <: Metabolite end
struct CoA <: Metabolite end
struct GlycolicAcid <: Metabolite end
struct Glycoyl <: FunctionalGroup{GlycolicAcid, Dehydroxy} end # Gc
struct LacticAcid <: Metabolite end
struct Lactyl <: FunctionalGroup{LacticAcid, Dehydroxy} end # Lt
struct PyruvicAcid <: Metabolite end
struct Pyruvyl <: FunctionalGroup{PyruvicAcid, Dehydroxy} end # Py
conjugation(t::Taurine) = Tauryl()

end