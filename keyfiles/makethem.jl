#!/usr/bin/env julia

import Combinatorics
import Iterators
using Bio.Seq
import Bio.Seq: ACGTN
import DataStructures

# From Bio.jl
# Calculate the Hamming distance between `seq1` and `seq2`.
function hamming_distance(seq1, seq2)
    @assert length(seq1) == length(seq2)
    n = 0
    for (x, y) in zip(seq1, seq2)
        n += x != y
    end
    return n
end


immutable Barcode
	name::String
	r1barcode::DNASequence
	r2barcode::Nullable{DNASequence}
	Barcode(n, s1::DNASequence, s2::DNASequence) = new(n, s1, Nullable{DNASequence}(s2))
	Barcode(n, s) = new(n, s, Nullable{DNASequence}())
end


function write_axe(io::IO, bcds::Vector{Barcode})
	for bcd in bcds
		r2bcd = isnull(bcd.r2barcode) ? "" : "\t$(get(bcd.r2barcode))"
		println(io, "$(bcd.name)\t$(bcd.r1barcode)$r2bcd")
	end
end

function all_kmers(k)


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
            if hamming_distance(candidate, j) < dist
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


function generate_index_set(name::String, idx_length::Int, num_indices::Vector{Int};
                            combinatorial::Bool=false, distance::Int=3)
    barcodes = Barcode[]
    all_indices = hamming_space(idx_length, distance)
    shuffle!(all_indices)

    n = prod(num_indices)
    if n > length(all_indices)
        error("More indicies requested than available with length $idx_length and distance $distance")
    end

    idx_names = alphacode(n)

    for (i, index) in enumerate(Iterators.product(map(x -> (1:x), num_indices)...))
        idx_seqs = collect(map(x -> all_indices[x], index))
        push!(barcodes, Barcode(idx_names[i], idx_seqs...))
    end
    return barcodes
end

g = generate_index_set("test", 9, Int[96, 12], combinatorial=true, distance=3)
