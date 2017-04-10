#!/usr/bin/env julia

Pkg.update()
PKGS = ["Bio", "JSON", "DataStructures", "Combinatorics", "Iterators", "YAML"]
for p in PKGS
    Pkg.add(p)
    eval(Expr(:import, Symbol(p))) # and import it so that it gets precompiled
end

