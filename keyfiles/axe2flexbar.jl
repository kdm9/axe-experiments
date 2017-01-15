#!/usr/bin/env julia


open(ARGS[1]) do infile
    for line in eachline(infile)
        bcd, name = split(rstrip(line))
        println(">$name")
        println(bcd)
    end
end
