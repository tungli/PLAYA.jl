function extract_reaction(s)
    reaction = match(r"\{(.*)\}", s)[1]
    @assert occursin("->", reaction)
    replace(reaction, " " => "") # remove whitespace
end

"""
Parses one block of BOLOS output file.
Returns a reaction string and an interpolation.
"""
function parse_bolos_block(block)
    table = split(strip(block), "\n")
    (reaction, table) = if table[1][1:8] == "E/N (Td)"
        y_quantity_name = strip(split(replace(table[1], "E/N (Td)" => ""), "(")[1])
        (y_quantity_name, table[2:end])
    elseif table[2][1:8] == "E/N (Td)"
        reaction = extract_reaction(table[1])
        (reaction, table[3:end])
    else
        error("`E/N (Td)` not found!")
    end
    table = map(split, table)
    table = [parse(Float64, x[i]) for x in table, i in 1:2]
    x = table[:, 1]
    y = table[:, 2]
    itp = interpolate((x,), y, Gridded(Linear()))
    if !isnothing(reaction)
        println("Loaded $(reaction) from BOLOS file")
    end

    (reaction, itp)
end

"""
Parses one block of BOLOS output file.
Return a `Dict{String, Interpolation}` where the `String` is the reaction.
"""
function parse_bolos_output(filename)
    content = open(filename) do file
        read(file, String) 
    end
    content = replace(content, "\r\n" => "\n")
    Dict([ parse_bolos_block(i) for i in split(strip(content), "\n\n") ])
end
