"""
Apply multiple replacements upon a `String`.
The replacements can be specified as an `Array{Pair{String, String}}`.
"""
function replace_seq(s::AbstractString, list::Array{Pair{String, String}})
    reduce((acc, x) -> replace(acc, x), list, init=s)
end


"""
    safe_replace(s, l)

Replaces multiple strings "safely" -- longer strings are replace first, and a
string cannot be replaced twice.

# Examples

```jldoctest
julia> safe_replace("NO(s) -> N2 + N", Dict(["N" => "M", "NO" => "N2", "N2" => "A2"]))
"N2(s) => A2 + M"
```
"""
function safe_replace(s, to_replace)
    ps = collect(to_replace)
    sorted = sort(ps, by= x-> length(x.first), rev=true)
    sorted_strings = map(x -> x.first, sorted)
    subs1 = map(x -> x => string(uuid1()), sorted_strings)

    dummies = map(x -> x.second, subs1)
    @assert allunique(dummies) "Not expected! (uuid not unique)"
    new_strings = map(x -> x.second, sorted)

    s2 = replace_seq(s, subs1)
    subs2 = [ i => j for (i, j) in zip(dummies, new_strings) ]
    replace_seq(s2, subs2)
end

"""
Parses one line from the reaction file.
"""
function parse_line(line)
    num, rest = split(line, ":", limit=2)
    num = strip(num)
    rest = strip(rest)
    reaction, rate = split(rest, "!", limit=2)
    reaction = strip(reaction)
    rate = strip(rate)
    [num, reaction, rate]
end

"""
Transform a line from the reaction file into a `Catalyst` specified reaction.
"""
function interpret_line(line, env, rates_list)
    n, r, k = parse_line(line)
    k = replace_num_d(k)
    k = replace_exp(k)
    rate = interpret_rate(k, env)
    rate_param_name = "p$n"
    push!(rates_list, (rate_param_name, rate))
    "$rate_param_name, $r    # $k"
end

"""
Replaces FORTRAN doubles for julia to parse, i.e. "7.d3" => "7.e3"
"""
function replace_num_d(s)
    for m in eachmatch(r"([0-9]|\.)d([0-9]|-)", s)
        j = m.offset + 1
        @assert(s[j] == 'd')
        s = join([ i == j ? 'e' : c for (i, c) in enumerate(s) ])
    end
    s
end

function replace_exp(s)
    replace(s, "**" => "^")
end

"""
Interprets a rate expression resulting in a `Real`.
"""
function interpret_rate(rate_expr::String,
                        env
                       )::Real
    # load env
    for (name, val) in env
        eval(:($name = $val))
    end

    m = match(r"(\w+)\s*\{(.*)\}", rate_expr)
    if m !== nothing
        key = replace(m[2], " " => "")
        fun = env[Symbol(m[1])]
        rate = fun(key, env)
        new_expr = replace(rate_expr, m.match => "$rate")
        println(new_expr)
        interpret_rate(new_expr, env)
    else
        eval(Meta.parse(rate_expr))
    end
end

"""
Parse reactions from a file and output a
[Catalyst](https://catalyst.sciml.ai/stable/) `ReactionSystem` and the rate
values.

# Inputs

* File name.
* Environment of variables that get substituted into rates. 
* Optional replacement list (`Array{Pair{String, String}}`) that translates the
species names (mainly to replace symbols that Catalyst complains about).

# Reaction file

On each line there are two  seperators: ':' and '!' and three data fields.
At the end of the line a comment starting with '#' may be inserted.
Each line is one (forward) reaction:

> `id`: `reaction` ! `rate expression` (# `optional comment`)

The `reaction` is specified same as a reaction in
[Catalyst](https://catalyst.sciml.ai/stable/tutorials/basics/) but without the
"rate part". This imposes limitations on the use of parenthesis and symbols
that represent functions in Julia.

Rate expression can be one either a function or a number. Evaluations are done
using the `env::Dict` which acts as an environment for variables.

Additionally, to avoid long functions, a "subdict" with functions can be used
with the following syntax:

> `subdict` { `key` }

which evaluates to:
```
env[Symbol(subdict)](key::String, env)
```

For example:
> BOLOS { A -> B } * 1e6 
will search for the `BOLOS` dict under the keyword `"A -> B"`, evaluate the
found function and substitute the corresponding value and evaluate the
multiplication.
"""
function parse_reactions(filename, env; to_replace=[])
    rates_list = Vector{Tuple{String, Float64}}(undef, 0)
    reactions = open(filename) do file
        [ interpret_line(line, env, rates_list) for line in eachline(file) ]
    end
    reactions = map(x -> safe_replace(x, to_replace), reactions)
    s = join(reactions, "\n")
    println(s)

    parameter_string = join(map(x -> x[1], rates_list), " ")
    catalyst_reactions = join(["@reaction_network begin\n", s, "\nend ", parameter_string ])
    reaction_system = eval(Meta.parse(catalyst_reactions))
    rate_values = map(x -> x[2], rates_list)
    (reaction_system, rate_values)
end


"""
Modify the system of ODEs, such that the concentration of `species` will be
fixed, i.e. its time derivative will be zero.

# Inputs
* `odesys` - system of ODEs (from DifferentialEquations.jl)
* `species_name` - array of species identifiers
"""
function keep_fixed!(odesys::ModelingToolkit.AbstractODESystem, species_name)
    species_num = findfirst(x -> string(x) == species_name, odesys.states)
    eq = odesys.eqs[species_num]
    odesys.eqs[species_num] = Equation(eq.lhs, Operation(0))
end




