using Iterators

function i2axe(fn::String)
    bcd = readdlm(fn)[2:end, 1]
    width = Int(ceil(log(length(bcd))/log(26)))
    lbl = sort(collect(map(x -> join([x...]),
                           product(repeated('A':'Z', width)...))))
    open("($fn).axe", "w") do ofh
        for (i, b) in enumerate(bcd)
            println(ofh, "$b\t$(lbl[i])")
        end
    end
end
