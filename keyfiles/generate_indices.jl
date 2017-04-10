#!/usr/bin/env julia
module BarcodeFactory

import Combinatorics
import Iterators
using Bio.Seq
import DataStructures
import YAML

immutable Barcode
	name::String
	r1barcode::DNASequence
	r2barcode::Nullable{DNASequence}
	Barcode(n, s1::DNASequence, s2::DNASequence) = new(n, s1, Nullable{DNASequence}(s2))
	Barcode(n, s) = new(n, s, Nullable{DNASequence}())
end

abstract Demuxer
immutable Axe <: Demuxer end
immutable AdapterRemoval <: Demuxer end
immutable Fastx <: Demuxer end
immutable Flexbar <: Demuxer end

function Base.write(io::IO, ::Type{Axe}, bcds::Vector{Barcode})
	for bcd in bcds
		r2bcd = isnull(bcd.r2barcode) ? "" : "\t$(get(bcd.r2barcode))"
		println(io, "$(bcd.r1barcode)$r2bcd\t$(bcd.name)")
	end
end

function Base.write{T<:Union{AdapterRemoval, Fastx}}(io::IO, ::Type{T},
                                                     bcds::Vector{Barcode})
	for bcd in bcds
		r2bcd = isnull(bcd.r2barcode) ? "" : "\t$(get(bcd.r2barcode))"
		println(io, "$(bcd.name)\t$(bcd.r1barcode)$r2bcd")
	end
end

function Base.write(io::IO, ::Type{Flexbar}, bcds::Vector{Barcode})
	for bcd in bcds
		println(io, ">$(bcd.name)\n$(bcd.r1barcode)")
	end
end

Base.write{T<:Demuxer}(::Type{T}, bcds) = write(STDOUT, T, bcds)


"""
Computes the cloud of words, w \in A^k, such that d(w_i, w_j) >= dist for all i
and j
"""
function hamming_space(k, dist)
    accepted = Vector{DNASequence}()
    for candidate_nt in Iterators.product(repeated(ACGT, k)...)
        candidate = DNASequence(DNANucleotide[candidate_nt...])
        acceptable = true
        for j in accepted
            if mismatches(candidate, j) < dist
                acceptable = false
                break
            end
        end
        if acceptable
            push!(accepted, candidate)
        end
    end
    return accepted
end

"Generate a list of alphabetic codes of length >= n"
function alphacode(n)
    A = collect('A':'Z')
    len = Int(ceil(log(n)/log(26)))
    codes = String[]
    for code in Iterators.product(repeated(A, len)...)
        push!(codes, join(code))
    end
    return codes
end


function generate_index_set(num_indices::Vector{Int}, idx_length::Int,
                            distance::Int=3)
    barcodes = Barcode[]
    all_index_seqs = hamming_space(idx_length, distance)
    shuffle!(all_index_seqs)

    totaln = prod(num_indices)
    maxn = maximum(num_indices)
    if maxn > length(all_index_seqs)
        error("More indicies requested than available with"*
              "length $idx_length and distance $distance")
    end

    idx_names = alphacode(totaln)

    for (i, index) in enumerate(Iterators.product(map(x -> (1:x), num_indices)...))
        idx_seqs = collect(map(x -> all_index_seqs[x], index))
        push!(barcodes, Barcode(idx_names[i], idx_seqs...))
    end
    return barcodes
end

const DEMUXERS = Dict{String, Type}(
        "axe" => Axe,
        "ar" => AdapterRemoval,
        "flexbar" => Flexbar,
        "fastx" => Fastx,
)

function generate_from_yaml(yamlfile::String, root::String)
    sets = YAML.load_file(yamlfile)
    for set in sets
        setname = set["name"]
        indices = generate_index_set(set["num_indices"], set["length"],
                                     set["dist"])
        for demuxer in set["demuxers"]
            ext = ".$demuxer"
            if demuxer == "flexbar"
                ext = "_flexbar.fasta"
            end
            open("$root/$setname$ext", "w") do outf
                write(outf, DEMUXERS[demuxer], indices)
            end
        end
    end
end

end # module BarcodeFactory

using ArgParse

function main()
    ap = ArgParseSettings()
    @add_arg_table ap begin
        "settings"
            help="YAML-encoding settings file"
            required=true
            arg_type=String
        "rootdir"
            help="Output directory"
            default="."
            required=false
            arg_type=String
    end
    args = parse_args(ap)

    BarcodeFactory.generate_from_yaml(args["settings"], args["rootdir"])
end

main()
