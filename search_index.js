var documenterSearchIndex = {"docs":
[{"location":"#PLAYA","page":"Home","title":"PLAYA","text":"","category":"section"},{"location":"#Documentation","page":"Home","title":"Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To add Package use:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> Pkg.add(\"https://github.com/tungli/PLAYA.jl.git\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"A try to replicate a study of air kinetics from [1] is in demo/air-kinetics.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = PLAYA","category":"page"},{"location":"","page":"Home","title":"Home","text":"parse_reactions\nsafe_replace\nparse_bolos_output\nkeep_fixed!","category":"page"},{"location":"#PLAYA.parse_reactions","page":"Home","title":"PLAYA.parse_reactions","text":"parse_reactions(filename, env; to_replace=[])\n\nParse reactions from a file and output a Catalyst ReactionSystem and the rate values.\n\nInputs\n\nFile name.\nEnvironment of variables that get substituted into rates. \nOptional replacement list (Array{Pair{String, String}}) that translates the\n\nspecies names (mainly to replace symbols that Catalyst complains about).\n\nReaction file\n\nOn each line there are two  seperators: ':' and '!' and three data fields. At the end of the line a comment starting with '#' may be inserted. Each line is one (forward) reaction:\n\nid: reaction ! rate expression (# optional comment)\n\nThe reaction is specified same as a reaction in Catalyst but without the \"rate part\". This imposes limitations on the use of parenthesis and symbols that represent functions in Julia.\n\nRate expression can be one either a function or a number. Evaluations are done using the env::Dict which acts as an environment for variables.\n\nAdditionally, to avoid long functions, a \"subdict\" with functions can be used with the following syntax:\n\nsubdict { key }\n\nwhich evaluates to:\n\nenv[Symbol(subdict)](key::String, env)\n\nFor example:\n\nBOLOS { A -> B } * 1e6 \n\nwill search for the BOLOS dict under the keyword \"A -> B\", evaluate the found function and substitute the corresponding value and evaluate the multiplication.\n\n\n\n\n\n","category":"function"},{"location":"#PLAYA.safe_replace","page":"Home","title":"PLAYA.safe_replace","text":"safe_replace(s, l)\n\nReplaces multiple strings \"safely\" – longer strings are replace first, and a string cannot be replaced twice.\n\n\n\n\n\n","category":"function"},{"location":"#PLAYA.parse_bolos_output","page":"Home","title":"PLAYA.parse_bolos_output","text":"Parses one block of BOLOS output file. Return a Dict{String, Interpolation} where the String is the reaction.\n\n\n\n\n\n","category":"function"},{"location":"#PLAYA.keep_fixed!","page":"Home","title":"PLAYA.keep_fixed!","text":"keep_fixed!(odesys::ModelingToolkit.AbstractODESystem, species_name)\n\nModify the system of ODEs, such that the concentration of species will be fixed, i.e. its time derivative will be zero.\n\nInputs\n\nodesys - system of ODEs (from DifferentialEquations.jl)\nspecies_name - array of species identifiers\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"[1]: A. Obrusník, P. Bílek, T. Hoder, M. Šimek, Z. Bonaventura, Plasma Sources Sci. Technol., 2018, DOI:10.1088/1361-6595/aad663. ","category":"page"}]
}