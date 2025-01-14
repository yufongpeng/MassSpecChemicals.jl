using MassSpecChemicals
using MassSpecChemicals.BioChemicals
using MassSpecChemicals.BioChemicals.Lipids
const MBL = MassSpecChemicals.BioChemicals.Lipids
using Test, JSON3

function test_annotationlevel(x)
    tal = MBL.annotationlevel(x.object; partial = true, additional = true, pass = true)
    length(tal) == length(x.annotationlevel) && all(t -> in(t, x.annotationlevel), tal)
end

@testset "MassSpecChemicals.jl" begin
    # Write your tests here.
    @testset "MassSpecChemicals.BioChemicals.Lipids" begin
        # use json
        test_lipid_js = JSON3.read(joinpath("data", "test_lipid.json"))
        @testset "io" begin
            global test_lipid = Dict{UnionAll, Dict}()
            for (c, s) in test_lipid_js
                @testset string(c) begin
                    dict = Dict{String, Any}()
                    for (l, a) in s
                        l = string(l)
                        lipid = MBL.parse_lipid(l)
                        show(lipid)
                        print(", ")
                        push!(dict, l => (object = lipid, annotationlevel = [eval(Meta.parse(aa)) for aa in a]))
                    end
                    push!(test_lipid, eval(c) => dict)
                end
            end
        end
        @testset "annotationlevel" begin
            
        end
    end
end
