#!/usr/bin/env julia
using Bio.Seq
import JSON
using DataStructures


function readfile_sampleids(filename::String)
    ctr = counter(String)
    open(FASTQReader, filename, quality_encoding=:sanger) do reader
        rec = FASTQSeqRecord{DNASequence}()
        while !eof(reader)
            read!(reader, rec)
            rdata = JSON.parse(rec.metadata.description)
            sample = rdata["id"]
            push!(ctr, sample)
        end
    end
    return ctr
end

immutable Assessment
    total::Int
    correct::Int
    incorrect::Int
    unassigned::Int
    Assessment() = new(0,0,0,0)
end

function assess_samples(truthfile, unknown_file, filenames)
    truecounts = JSON.parsefile(truthfile)
    unknown_samples = readfile_sampleids(unknown_file)
    sample_stats = Dict{String, Assessment}()
    for fname in filenames
        sample, _ = splitext(basename(fname))
        got = readfile_sampleids(fname)
        nread_expect = truecounts[sample]
        nread_got_total = sum(values(got))
        nread_got_sample = got[sample]
        nread_unassigned = unknown_samples[sample]
        # Wrong is those not unassigned, and not in correct file. Must be
        # elsewhere
        nread_got_wrong = nread_expect - nread_unassigned 
        a = Assessment(nread_expect, nread_got_sample, nread_got_wrong,
                       nread_unassigned)
        sample_stats[sample] = a
    end
    return sample_stats
end

function append!(a::Assessment, b::Assessment)
    a.total += b.total
    a.correct += b.correct
    a.incorrect += b.incorrect
    a.unassigned += b.unassigned
end

function summarise_assessment(assessment)
    total = Assessment()
    for (s, a) in assessment
        append!(total, a)
    end
    return total
end
    
function dump_assessment(file, assess)
    open(file, "w") do f
        prinln(f, "Sample\tTotal\tCorrect\tIncorrect\tUnassigned")
        for (s, a) in assess
            println(f, "$s\t$(a.total)\t$(a.correct)\t$(a.incorrect)\t$(a.unassigned)")
        end
    end
end


##################
#  END OF FUNCS  #
##################


const afile = ARGS[1]
const seed = ARGS[2]
const barcode_set = ARGS[3]
const demuxer = ARGS[4]
const truename = ARGS[5]
const unknown = ARGS[6]
const filenames = ARGS[7:end]

const assessment = assess_samples(truename, unknown, filenames)
dump_assessment(afile, assessment)

const tot = summarise_assessment(assessment)

println("$seed\t$demuxer\t$barcode_set\t$(tot.correct)\t$(tot.incorrect)\t$(tot.unassigned)")
