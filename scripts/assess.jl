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
            if endswith(rec.name, "/2")
                # Skip R2, otherwise counts are inflated
                continue
            end
            rdata = JSON.parse(rec.metadata.description)
            sample = rdata["id"]
            push!(ctr, sample)
        end
    end
    return ctr
end

type Assessment
    total::Int
    correct::Int
    incorrect::Int
    unassigned::Int
    Assessment() = new(0,0,0,0)
    Assessment(t, c, i, u) = new(t, c, i, u)
end

function assess_samples(truthfile, sampledir)
    # Read true counts
    samplecounts = JSON.parsefile(truthfile)

    # Find unknown file and read it
    unknowncounts = readfile_sampleids(joinpath(sampledir, "unknown.fastq"))

    sample_stats = Dict{String, Assessment}()
    for sample in keys(samplecounts)
        fname = joinpath(sampledir, "$(sample).fastq")
        got = readfile_sampleids(fname)
        nread_expect = samplecounts[sample]
        nread_got_sample = got[sample]
        nread_unassigned = unknowncounts[sample]
        # Wrong is those not unassigned, and not in correct file. Must be
        # elsewhere
        nread_got_wrong = nread_expect - nread_unassigned -  nread_got_sample
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
        println(f, "Sample\tTotal\tCorrect\tIncorrect\tUnassigned")
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
const truthfile = ARGS[5]
const sampledir = ARGS[6]

const assessment = assess_samples(truthfile, sampledir)
dump_assessment(afile, assessment)

const tot = summarise_assessment(assessment)

println("$seed\t$demuxer\t$barcode_set\t$(tot.correct)\t$(tot.incorrect)\t$(tot.unassigned)")
