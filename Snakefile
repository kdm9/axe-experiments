import simhelpers
import random

random.seed(3301)
NREPS = 2
SEEDS = random.sample(range(10000), NREPS)
# FIXME: Increase this to something more reasonable for any publication version
NREADS = 100000
KEEP_READS = False
READLEN = 101 # for benchmarking purposes
BC_DM = [
    ('4bp-se', 'axe'),
    ('nested-se', 'axe'),
    ('gbs-se', 'axe'),
    ('4bp-pe', 'axe'),
    ('nested-pe', 'axe'),
    ('4bp-se', 'flexbar'),
    ('nested-se', 'flexbar'),
    ('gbs-se', 'flexbar'),
    ('bvzlab-pst1-gbs-pe', 'axe'),
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
        ["data/assessments/{seed}_{bc}_{dm}.tab".format(dm=dm, bc=bc, seed=seed)
            for bc, dm in BC_DM
            for seed in SEEDS],


rule clean:
    shell:
        'rm -rf data'


rule assess:
    input:
        stats="data/bcd_stats/{seed}_{barcode}.json",
        reads=dynamic('data/demuxed.{demuxer}/{seed}_{barcode}/{sample}.fastq'),
    output:
        'data/assessments/{seed}_{barcode}_{demuxer}.tab'
    shell:
        "scripts/assess-barcodes.py {input} > {output}"


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

rule flexbar_dm:
    input:
        kf='keyfiles/{barcode}_flexbar.fasta',
        reads="data/reads/{seed}_{barcode}.fastq.gz",
    output:
        dynamic('data/demuxed.flexbar/{seed}_{barcode}/{sample}.fastq')
        #lambda wc: expand('data/demuxed.flexbar/{{seed}}_{{barcode}}/{s}.fastq',
        #                  s=BARCODE_NAMES[wc.barcode])
    params:
        outdir=lambda w: 'data/demuxed.flexbar/{}_{}/'.format(w.seed, w.barcode),
        combo=lambda w: '-c' if 'se' in w.barcode else '',
    log:
        'data/log/flexbar/{seed}_{barcode}.log'
    benchmark:
        'data/benchmarks/flexbar/{seed}_{barcode}.txt'
    shell:
        '(flexbar'
        '   -be LEFT'  # 5' barcodes
        '   -bu'
        '   -b {input.kf}'
        '   -r {input.reads}'
        '   -t {params.outdir} &&'
        " pushd {params.outdir} &&"
        " ls -lahF . &&"
        " rename 's/^.*_barcode_([\w]+.fastq)$/$1/' *.fastq &&"
        " mv unassigned.fastq unkown.fastq &&"
        " popd"
        ") >{log} 2>&1"


rule bcdreads:
    input:
        barcode="keyfiles/{barcode}.axe",
        r1="data/tmp/{seed}_{barcode}-R1.prebcd.fastq",
        r2="data/tmp/{seed}_{barcode}-R2.prebcd.fastq",
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
        barcode="keyfiles/{barcode}.axe",
    output:
        r1=temp("data/tmp/{seed}_{barcode}-R1.prebcd.fastq"),
        r2=temp("data/tmp/{seed}_{barcode}-R2.prebcd.fastq"),
    log:
        "data/log/reads/{seed}_{barcode}.log"
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
