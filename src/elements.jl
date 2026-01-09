# AME 2020
const MASS = Dict(
    ""          => 0u"g",
    "C"         => 12u"g",
    "[12C]"     => 12u"g",
    "[13C]"     => 13.003354835336u"g",
    # "[14C]"   => 14.003241989u"g",
    "H"         => 1.007825031898u"g",
    "[1H]"      => 1.007825031898u"g",
    "D"         => 2.014101777844u"g",
    "[2H]"      => 2.014101777844u"g",
    "O"         => 15.994914619257u"g",
    "[16O]"     => 15.994914619257u"g",
    "[17O]"     => 16.999131755953u"g",
    "[18O]"     => 17.999159612136u"g",
    "N"         => 14.003074004251u"g",
    "[14N]"     => 14.003074004251u"g",
    "[15N]"     => 15.000108898266u"g",
    "P"         => 30.973761997677u"g",
    "[31P]"     => 30.973761997677u"g",
    "S"         => 31.9720711735u"g",
    "[32S]"     => 31.9720711735u"g",
    "[33S]"     => 32.9714589086u"g",
    "[34S]"     => 33.96786701u"g",
    "[36S]"     => 35.96708069u"g",
    "Li"        => 7.016003434u"g",
    "[7Li]"     => 7.016003434u"g",
    "[6Li]"     => 6.0151228874u"g",
    "Na"        => 22.9897692820u"g",
    "[23Na]"    => 22.9897692820u"g",
    "K"         => 38.963706485u"g",
    "[39K]"     => 38.963706485u"g",
    "[40K]"     => 39.96399817u"g",
    "[41K]"     => 40.961825256u"g",
    "F"         => 18.9984031621u"g",
    "[19F]"     => 18.9984031621u"g",
    "Cl"        => 34.96885269u"g",
    "[35Cl]"    => 34.96885269u"g",
    "[37Cl]"    => 36.96590257u"g",
    "Ag"        => 106.9050915u"g",
    "[108Ag]"   => 106.9050915u"g",
    "[109Ag]"   => 108.9047558u"g",
    "Se"        => 79.916521761u"g",
    "[74Se]"    => 73.922475933u"g",
    # "[75Se]"    => 74.922522870u"g",
    "[76Se]"    => 75.919213702u"g",
    "[77Se]"    => 76.919914150u"g",
    "[78Se]"    => 77.917309244u"g",
    "[80Se]"    => 79.916521761u"g",
    "[82Se]"    => 81.916699531u"g"
    # "Mg"        => 23.98504168u"g",
    # "[24Mg]"    => 23.98504168u"g",
    # "[25Mg]"    => 24.985836966u"g",
    # "[26Mg]"    => 25.982592972u"g"
)

# CIAAW
const ABUNDANCE = Dict(
    ""          => 1,
    "C"         => 0.989165,
    "[12C]"     => 0.989165,
    "[13C]"     => 0.010835,    # C3: [0.010 674, 0.010 827], C4: [0.010 881, 0.010 958], use mean
    # "[14C]"     => 1e-12,
    "H"         => 0.9998576,
    "[1H]"      => 0.9998576,
    "D"         => 0.0001424,   # Non-marine organisms: [0.000 1188, 0.000 1660]
    "[2H]"      => 0.0001424,   # Non-marine organisms: [0.000 1188, 0.000 1660]
    "O"         => 0.9975835,
    "[16O]"     => 0.9975835,
    "[17O]"     => 0.0003835,   # [0.000 367, 0.000 400]
    "[18O]"     => 0.0020330,   # Cellulose, lipids, and tissue: [0.001 9918, 0.002 0742]
    "N"         => 0.99637,
    "[14N]"     => 0.99637,
    "[15N]"     => 0.00363,     # Plants and animals: [0.003 484, 0.003 776]
    "P"         => 1,
    "[31P]"     => 1,
    "S"         => 0.9500205, 
    "[32S]"     => 0.9500205, 
    "[33S]"     => 0.00763,     # [0.007 29, 0.007 97]
    "[34S]"     => 0.0421915,   # Animals: [0.041571, 0.042812]
    "[36S]"     => 0.000158,    # [0.000 129, 0.000 187]
    "Li"        => 0.92409,     # LSVEC
    "[7Li]"     => 0.92409,     # LSVEC
    "[6Li]"     => 0.07591, 
    "Na"        => 1,
    "[23Na]"    => 1,
    "K"         => 0.932581,
    "[39K]"     => 0.932581,
    "[40K]"     => 0.000117,
    "[41K]"     => 0.067302,
    "F"         => 1,
    "[19F]"     => 1,
    "Cl"        => 0.75773,
    "[35Cl]"    => 0.75773,
    "[37Cl]"    => 0.24227,     # SMOC
    "Ag"        => 0.51839,
    "[108Ag]"   => 0.51839,
    "[109Ag]"   => 0.48161,
    "Se"        => 0.49803,
    "[74Se]"    => 0.008393,
    # "[75Se]"    => 74.922522870u"g",
    "[76Se]"    => 0.09237,
    "[77Se]"    => 0.07607,
    "[78Se]"    => 0.236922,
    "[80Se]"    => 0.49803,
    "[82Se]"    => 0.088215
    # "Mg"        => 0.78965, # [0.7888, 0.7905]
    # "[24Mg]"    => 0.78965, # [0.7888, 0.7905]
    # "[25Mg]"    => 0.1001, # [0.099 88, 0.100 34]
    # "[26Mg]"    => 0.11025 # [0.1096, 0.1109]
)

const ISOTOPES = Dict(
    ""   => [""],
    "C"  => ["[12C]", "[13C]"],
    "H"  => ["[1H]", "D"],
    "O"  => ["[16O]", "[18O]", "[17O]"],
    "N"  => ["[14N]", "[15N]"],
    "P"  => ["[31P]"],
    "S"  => ["[32S]", "[34S]", "[33S]", "[36S]"],
    "Li" => ["[7Li]", "[6Li]"],
    "Na" => ["[23Na]"],
    "K"  => ["[39K]", "[41K]", "[40K]"],
    "F"  => ["[19F]"],
    "Cl" => ["[35Cl]", "[37Cl]"],
    "Ag" => ["[108Ag]", "[109Ag]"],
    "Se" => ["[80Se]", "[78Se]", "[76Se]", "[82Se]", "[77Se]", "[74Se]"]
    # "Mg" => ["[24Mg]", "[26Mg]", "[25Mg]"]
) 

const PARENTS = Dict(
    ""          => "",
    "C"         => "C",
    "[12C]"     => "C",
    "[13C]"     => "C",
    # "[14C]"   => "C",
    "H"         => "H",
    "[1H]"       => "H",
    "D"         => "H",
    "[2H]"      => "H",
    "O"         => "O",
    "[16O]"     => "O",
    "[17O]"     => "O",
    "[18O]"     => "O",
    "N"         => "N",
    "[14N]"     => "N",
    "[15N]"     => "N",
    "P"         => "P",
    "[31P]"     => "P",
    "S"         => "S", 
    "[32S]"     => "S", 
    "[33S]"     => "S",
    "[34S]"     => "S",
    "[36S]"     => "S",
    "Li"        => "Li",
    "[7Li]"     => "Li",
    "[6Li]"     => "Li", 
    "Na"        => "Na",
    "[23Na]"    => "Na",
    "K"         => "K",
    "[39K]"     => "K",
    "[40K]"     => "K",
    "[41K]"     => "K",
    "F"         => "F",
    "[19F]"     => "F",
    "Cl"        => "Cl",
    "[335Cl]"   => "Cl",
    "[37Cl]"    => "Cl",
    "Ag"        => "Ag",
    "[108Ag]"   => "Ag",
    "[109Ag]"   => "Ag",
    "Se"        => "Se",
    "[74Se]"    => "Se",
    # "[75Se]"    => "Se",
    "[76Se]"    => "Se",
    "[77Se]"    => "Se",
    "[78Se]"    => "Se",
    "[80Se]"    => "Se",
    "[82Se]"    => "Se"
    # "Mg"        => "Mg",
    # "[24Mg]"    => "Mg",
    # "[25Mg]"    => "Mg",
    # "[26Mg]"    => "Mg"
)

const DECODES = Dict(
    ""          => "",
    "C"         => "C",
    "Citz"      => "[12C]",
    "Citn"      => "[13C]",
    # "Citnn"     => "C",
    "H"         => "H",
    "Hitz"      => "[1H]",
    "D"         => "D",
    "Hitn"      => "D",
    "O"         => "O",
    "Oitz"      => "[16O]",
    "Oitn"      => "[17O]",
    "Oitnn"     => "[18O]",
    "N"         => "N",
    "Nitz"      => "[14N]",
    "Nitn"      => "[15N]",
    "P"         => "P",
    "Pitz"      => "[31P]",
    "S"         => "S", 
    "Sitz"      => "[32S]",
    "Sitn"      => "[33S]",
    "Sitnn"     => "[34S]",
    "Sitnnnn"   => "[36S]",
    "Li"        => "Li",
    "Liitz"     => "[7Li]", 
    "Liitp"     => "[6Li]", 
    "Na"        => "Na",
    "Naitz"     => "[23Na]",
    "K"         => "K",
    "Kitz"      => "[39K]",
    "Kitn"      => "[40K]",
    "Kitnn"     => "[41K]",
    "F"         => "F",
    "Fitz"      => "[19F]",
    "Cl"        => "Cl",
    "Clitz"     => "[35Cl]",
    "Clitnn"    => "[37Cl]",
    "Ag"        => "Ag",
    "Agitz"     => "[108Ag]",
    "Agitn"     => "[109Ag]",
    "Se"        => "Se",
    "Seitpppppp"    => "[74Se]",
    # "[75Se]"    => "Se",
    "Seitpppp"  => "[76Se]",
    "Seitppp"   => "[77Se]",
    "Seitpp"    => "[78Se]",
    "Seitz"     => "[80Se]",
    "Seitnn"    => "[82Se]"
    # "Mg"        => "Mg",
    # "Mgitz"     => "[24Mg]",
    # "Mgitn"     => "[25Mg]",
    # "Mgitnn"    => "[26Mg]"
)

"""
    const ELEMENTS

Constants related to elements. It is a dictionary with five keys
* `:MASS`: atomic mass.
* `:ABUNDANCE`: natrural abundance.
* `:ISOTOPES`: possible isotopes.
* `:PARENTS`: parent element. 
* `:DECODES`: decode encoded string for `parse_compound`.

# Parent Elements and Major isotopes
|Symbol|Major isotopes|Atomic number|Mass number|
|--------|-------------|-----------|------------------|
|C|[12C]|6|12|
|H|[1H]|1|1|
|O|[16H]|8|16|
|N|[14N]|7|14|
|P|[31P]|15|31|
|S|[32S]|16|32|
|Li|[7Li]|3|7|
|Na|[23Na]|11|23|
|K|[39K]|19|39|
|F|[19F]|9|19|
|Cl|[35Cl]|17|35|
|Ag|[108Ag]|47|108|
|Se|[80Se]|34|80|

# Minor Isotopes
|Minor isotopes|Atomic number|Mass number|Alternative symbol|
|--------|-------------|-----------|------------------|
|[13C]|6|13||
|D|1|2|[2H]|
|[17O]|8|17||
|[18O]|8|18||
|[15N]|7|15||
|[33S]|16|33||
|[34S]|16|34||
|[36S]|16|36||
|[6Li]|3|6||
|[40K]|19|40||
|[41K]|19|41||
|[37Cl]|17|37||
|[109Ag]|47|109||
|[74Se]|34|74||
|[76Se]|34|76||
|[77Se]|34|77||
|[78Se]|34|78||
|[82Se]|34|82||

By default, parent elements are considered as major isotopes possibly replaced by minor isotopes. For instance,
* CO2 has a carbon-12 and two oxygen-16, but any minor isotopes replacements are possible.
* [13C][16O]O has a carbon-13, an oxygen-16, and an oxygen-16 possibly replaced by other minor isotopes.

One exception is that in `parent` chemical of `Isotopomers`, parent elements are major isotopes, and the number of replacement is restricted by field `isotopes`. 
"""
const ELEMENTS = Dict{Symbol, Dict}(
    :MASS       => deepcopy(MASS),
    :ABUNDANCE  => deepcopy(ABUNDANCE),
    :ISOTOPES   => deepcopy(ISOTOPES),
    :PARENTS    => deepcopy(PARENTS),
    :DECODES    => deepcopy(DECODES)
)

"""
    set_elements!(element::AbstractString, mass, abundance; minor_name = nothing)

Update or insert `elements` in `ELEMENTS`. 

* `mass`: atomic mass of all isotopes.
* `abundance`: natural abundance of all isotopes.
* `minor_name`: customized minor element names.
"""
function set_elements!(element::AbstractString, mass, abundance; minor_name = nothing)
    isotopes = isnothing(minor_name) ? map(mass) do m 
        n = round(Int, m)
        string("[", n, element, "]")
    end : minor_name
    _, i = findmax(abundance)
    ELEMENTS[:PARENTS][element] = element
    ELEMENTS[:DECODES][element] = element
    ELEMENTS[:MASS][element] = mass[i]
    ELEMENTS[:ABUNDANCE][element] = abundance[i]
    ee = encode_isotopes.(isotopes)
    for (e, i, m, a) in zip(ee, isotopes, mass, abundance) 
        ELEMENTS[:DECODES][e] = i 
        ELEMENTS[:MASS][i] = m 
        ELEMENTS[:ABUNDANCE][i] = a 
        ELEMENTS[:PARENTS][i] = element
    end
    id = sortperm(abundance; rev = true)
    ELEMENTS[:ISOTOPES][element] = isotopes[id]
    ELEMENTS
end

const ME = 0.00054857990924u"g"