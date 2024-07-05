module AminoAcids
using Reexport
@reexport using ..BioChemicals
using ..MassSpecChemicals: AbstractChemical
import ..BioChemicals: originalmolecule, leavinggroup, conjugation
export AminoAcid, Glycine, Serine, Ornithine
abstract type AminoAcid <: AbstractChemical end

struct Alanine <: AminoAcid end # G
struct Alanyl <: FunctionalGroup{Alanine, Dehydroxy} end # Ala
struct Glycine <: AminoAcid end # G
struct Glycyl <: FunctionalGroup{Glycine, Dehydroxy} end # G/Gly
struct Serine <: AminoAcid end # S
struct Ornithine <: AminoAcid end # Orn
conjugation(aa::AminoAcid) = Substituent{Ndehydrogen}(aa)

end