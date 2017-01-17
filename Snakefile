import simhelpers
import random

random.seed(3301)
NREPS = 1
SEEDS = random.sample(range(10000), NREPS)
# FIXME: Increase this to something more reasonable for any publication version
NREADS = 100000
KEEP_READS = False
READLEN = 101 # for benchmarking purposes
BC_DM = [
    ('8bp-se', 'axe'),
    ('8bp-se', 'flexbar'),
    ('8bp-se', 'fastx'),
    ('nextera-all-se', 'axe'),
    ('nextera-all-se', 'flexbar'),
    ('nextera-all-se', 'fastx'),
    ('nextera-all-pe', 'axe'),
    ('bvzlab-pst1-gbs-pe', 'axe'),
    ('gbs-se', 'axe'),
    ('gbs-se', 'flexbar'),
    ('gbs-se', 'fastx'),
    ('nested-se', 'axe'),
    ('nested-se', 'fastx'),
#    ('nested-pe', 'axe'),
]
BARCODE_SETS = list(set([bd[0] for bd in BC_DM]))
BARCODE_NAMES = {bcd: simhelpers.keyfile_names("keyfiles/{}.axe".format(bcd))
                 for bcd in BARCODE_SETS}

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

rule all:
    input:
        expand("data/reads/{seed}_{barcode}.fastq.gz", seed=SEEDS,
               barcode=BARCODE_SETS) if KEEP_READS else [],
        ["data/assessments/{seed}_{bc}_{dm}.tsv".format(dm=dm, bc=bc, seed=seed)
            for bc, dm in BC_DM
            for seed in SEEDS],
        "data/stats/all.tsv"


rule clean:
    shell:
        'rm -rf data'

rule global_assess:
    input:
        ["data/stats/{seed}_{bc}_{dm}.tsv".format(dm=dm, bc=bc, seed=seed)
            for bc, dm in BC_DM
            for seed in SEEDS],
    output:
        "data/stats/all.tsv"
    run:
        with open(output[0], "w") as ofh:
            print("Seed\tDemuxer\tBarcodeSet\tCorrect\tIncorrect\tUnassigned", file=ofh)
            for infn in input:
                with open(infn) as ifh:
                    ofh.write(ifh.read())

rule assess:
    input:
        stats="data/bcd_stats/{seed}_{barcode}.json",
        reads=dynamic('data/demuxed.{demuxer}/{seed}_{barcode}/{sample}.fastq'),
    output:
        afile='data/assessments/{seed}_{barcode}_{demuxer}.tsv',
        summary=temp('data/stats/{seed}_{barcode}_{demuxer}.tsv'),
    shell:
        "scripts/assess.jl"
        "   {output.afile}"
        "   {wildcards.seed}"
        "   {wildcards.barcode}"
        "   {wildcards.demuxer}"
        "   {input.stats}"
        "   data/demuxed.{wildcards.demuxer}/{wildcards.seed}_{wildcards.barcode}"
        "> {output.summary}"


rule axe:
    input:
        kf='keyfiles/{barcode}.axe',
        reads="data/reads/{seed}_{barcode}.fastq.gz",
    output:
        dynamic('data/demuxed.axe/{seed}_{barcode}/{sample}.fastq')
    params:
        outdir=lambda w: 'data/demuxed.axe/{}_{}/'.format(w.seed, w.barcode),
        combo=lambda w: '-c' if 'pe' in w.barcode else '',
        inmode=lambda w: '-f' if 'se' in w.barcode else '-i',
        outmode=lambda w: '-F' if 'se' in w.barcode else '-I',
    log:
        'data/log/axe/{seed}_{barcode}.log'
    benchmark:
        'data/benchmarks/axe/{seed}_{barcode}.txt'
    shell:
        '(axe-demux'
        '   {params.combo}'
        '   -b {input.kf}'
        '   {params.inmode} {input.reads}'
        '   {params.outmode} {params.outdir}'
        ' && pushd {params.outdir}'
        " && rename 's/_(il|R1|R2)//' *.fastq"
        ' && popd'
        ") >{log} 2>&1"

rule fastx_dm:
    input:
        kf='keyfiles/{barcode}.fastx',
        reads="data/reads/{seed}_{barcode}.fastq.gz",
    output:
        dynamic('data/demuxed.fastx/{seed}_{barcode}/{sample}.fastq')
    params:
        outdir=lambda w: 'data/demuxed.fastx/{}_{}/'.format(w.seed, w.barcode),
        combo=lambda w: '-c' if 'se' in w.barcode else '',
    log:
        'data/log/fastx/{seed}_{barcode}.log'
    benchmark:
        'data/benchmarks/fastx/{seed}_{barcode}.txt'
    shell:
        '(zcat {input.reads} | '
        '   fastx_barcode_splitter.pl'
        '   --bol'
        '   --bcfile {input.kf}'
        '   --suffix .fastq'
        '   --prefix {params.outdir} &&'
        ' pushd {params.outdir} &&'
        ' mv unmatched.fastq unknown.fastq &&'
        ' popd'
        ') >{log} 2>&1'


rule flexbar_dm:
    input:
        kf='keyfiles/{barcode}_flexbar.fasta',
        reads="data/reads/{seed}_{barcode}.fastq.gz",
    output:
        dynamic('data/demuxed.flexbar/{seed}_{barcode}/{sample}.fastq')
    params:
        outdir=lambda w: 'data/demuxed.flexbar/{}_{}/'.format(w.seed, w.barcode),
    log:
        'data/log/flexbar/{seed}_{barcode}.log'
    benchmark:
        'data/benchmarks/flexbar/{seed}_{barcode}.txt'
    shell:
        '(flexbar'
        '   -be LEFT'  # 5' barcodes
        '   -bu'
        '   --min-read-length 1'
        '   -b {input.kf}'
        '   -r {input.reads}'
        '   -t {params.outdir} &&'
        " pushd {params.outdir} &&"
        " ls -lahF . &&"
        " rename 's/^.*_barcode_([\w]+.fastq)$/$1/' *.fastq &&"
        " mv unassigned.fastq unknown.fastq &&"
        " popd"
        ") >{log} 2>&1"


rule bcdreads:
    input:
        barcode="keyfiles/{barcode}.axe",
        r1="data/tmp/{seed}-R1.prebcd.fastq",
        r2="data/tmp/{seed}-R2.prebcd.fastq",
    params:
        re_site=lambda w: '-r TGCAG' if 'gbs' in w.barcode else '',
        singleend=lambda w: '-S' if 'se' in w.barcode else '',
    output:
        reads="data/reads/{seed}_{barcode}.fastq.gz",
        stats="data/bcd_stats/{seed}_{barcode}.json",
    log:
        "data/log/bcdreads/{seed}_{barcode}.log"
    shell:
        "(scripts/add-barcodes.py -s {wildcards.seed}"
        "   {params.re_site}"
        "   {params.singleend}"
        "   -y {output.stats}"
        "   {input.barcode}"
        "   {input.r1}"
        "   {input.r2} |"
        "  gzip -4 >{output.reads}) 2>{log}"


rule reads:
    input:
        gen="data/genomes/{seed}.fa",
    output:
        r1=temp("data/tmp/{seed}-R1.prebcd.fastq"),
        r2=temp("data/tmp/{seed}-R2.prebcd.fastq"),
    log:
        "data/log/reads/{seed}.log"
    shell:
        "mason_simulator --seed {wildcards.seed}"
        "   --illumina-read-length {READLEN}"
        "   -n {NREADS}"
        "   -ir {input.gen}"
        "   -o {output.r1}"
        "   -or {output.r2} >{log} 2>&1"


rule randgenome:
    output:
        "data/genomes/{seed}.fa"
    params:
        genome_size=100000
    log:
        "data/log/genome/{seed}.log"
    shell:
        "scripts/randfa.py"
        "   {params.genome_size}"
        "   {wildcards.seed}"
        " > {output}"
        " 2> {log}"
