from collections import defaultdict
import random
import yaml

random.seed(int(config.get('seed', 3301)))
SEEDS = random.sample(range(10000), config.get('num_reps', 1))

with open("keyfiles/index_sets.yaml") as fh:
    INDEX_SETS = yaml.load(fh)


DEMUXER_SETS = defaultdict(list)
for s in INDEX_SETS:
    whitelist = config.get("set_whitelist", None)
    blacklist = config.get("set_blacklist", None)
    if whitelist is not None and s["name"] not in whitelist:
        continue
    if blacklist is not None and s["name"] in blacklist:
        continue
    for demuxer in s["demuxers"]:
        DEMUXER_SETS[demuxer].append(s["name"])

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

def all_assessments(w):
    return ["data/stats/accuracy/{seed}_{bc}_{dm}.tsv".format(dm=dm, bc=bc, seed=seed)
            for dm, bcds in DEMUXER_SETS.items()
            for bc in bcds
            for seed in SEEDS]

rule all:
    input:
        expand("data/reads/{seed}_{barcode}.fastq.gz", seed=SEEDS,
               barcode=[s["name"] for s in INDEX_SETS]) if config.get('keep_reads', False) else [],
        "data/stats/accuracy_summary.tsv",
        "data/stats/timing_summary.tsv",
       #all_assessments,


rule clean:
    shell:
        'rm -rf data'

rule global_assess:
    input:
        all_assessments,
    output:
        acc="data/stats/accuracy_summary.tsv",
        time="data/stats/timing_summary.tsv"
    run:
        import glob
        from os.path import dirname, basename
        with open(output.acc, "w") as ofh:
            print("Seed\tDemuxer\tBarcodeSet\tCorrect\tIncorrect\tUnassigned", file=ofh)
            for infn in input:
                with open(infn) as ifh:
                    ofh.write(ifh.read())

        with open(output.time, "w") as ofh:
            print("Seed\tDemuxer\tBarcodeSet\tSec", file=ofh)
            benchmarks = glob.glob('data/benchmarks/*/*.txt')
            for infn in benchmarks:
                dm = basename(dirname(infn))
                if dm not in set(DEMUXER_SETS.keys()):
                    continue
                sb = basename(infn)
                sb = sb[:-len(".txt")]
                seed, bcd = sb.split("_")
                with open(infn) as ifh:
                    next(ifh)  # Skip header
                    for line in ifh:
                        s, tm = line.strip().split()
                        print(seed, dm, bcd, s, sep="\t", file=ofh)


rule assess:
    input:
        stats="data/bcd_stats/{seed}_{barcode}.json",
        reads='data/demuxed/{demuxer}/{seed}_{barcode}/reads.tar',
    output:
        afile='data/assessments/{seed}_{barcode}_{demuxer}.tsv',
        summary=temp('data/stats/accuracy/{seed}_{barcode}_{demuxer}.tsv'),
    log:
        "data/log/assess/{seed}_{barcode}_{demuxer}.log"
    shell:
        "tmp=$(mktemp -d -p data/);"
        "(tar xf {input.reads} -C $tmp --strip-components=4 &&"
        " scripts/assess.jl"
        "   {output.afile}"
        "   {wildcards.seed}"
        "   {wildcards.barcode}"
        "   {wildcards.demuxer}"
        "   {input.stats}"
        "   $tmp > {output.summary}"
        ") >{log} 2>&1;"
        "rm -rf $tmp"


rule axe:
    input:
        kf='keyfiles/{barcode}.axe',
        reads="data/reads/{seed}_{barcode}.fastq.gz",
    output:
        'data/demuxed/axe/{seed}_{barcode}/reads.tar'
    params:
        outdir=lambda w: 'data/demuxed/axe/{}_{}/'.format(w.seed, w.barcode),
        combo=lambda w: '-c' if 'pe' in w.barcode else '',
        inmode=lambda w: '-f' if 'se' in w.barcode else '-i',
        outmode=lambda w: '-F' if 'se' in w.barcode else '-I',
    log:
        'data/log/axe/{seed}_{barcode}.log'
    benchmark:
        'data/benchmarks/axe/{seed}_{barcode}.txt'
    threads: 1
    shell:
        '(axe-demux'
        '   {params.combo}'
        '   -b {input.kf}'
        '   {params.inmode} {input.reads}'
        '   {params.outmode} {params.outdir}'
        ' && pushd {params.outdir}'
        " && rename 's/_(il|R1|R2)//' *.fastq"
        ' && popd'
        " && tar cvf {output} {params.outdir}"
        " && rm -vf {params.outdir}/*.fastq"
        ") >{log} 2>&1"

rule fastx_dm:
    input:
        kf='keyfiles/{barcode}.fastx',
        reads="data/reads/{seed}_{barcode}.fastq.gz",
    output:
        'data/demuxed/fastx/{seed}_{barcode}/reads.tar'
    params:
        outdir=lambda w: 'data/demuxed/fastx/{}_{}/'.format(w.seed, w.barcode),
        combo=lambda w: '-c' if 'se' in w.barcode else '',
    log:
        'data/log/fastx/{seed}_{barcode}.log'
    benchmark:
        'data/benchmarks/fastx/{seed}_{barcode}.txt'
    threads: 1
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
        " && tar cvf {output} {params.outdir}"
        " && rm -vf {params.outdir}/*.fastq"
        ') >{log} 2>&1'


rule ar_dm:
    input:
        kf='keyfiles/{barcode}.ar',
        reads="data/reads/{seed}_{barcode}.fastq.gz",
    output:
        "data/demuxed/ar/{seed}_{barcode}/reads.tar"
    params:
        outdir=lambda w: 'data/demuxed/ar/{}_{}/'.format(w.seed, w.barcode),
        il=lambda w: "--interleaved" if 'se' not in w.barcode else ''
    log:
        'data/log/ar/{seed}_{barcode}.log'
    benchmark:
        'data/benchmarks/ar/{seed}_{barcode}.txt'
    threads: 1
    shell:
        '( AdapterRemoval'
        '   {params.il}'
        '   --file1 {input.reads}'
        '   --demultiplex-only'
        '   --barcode-list {input.kf}'
        '   --basename {params.outdir}ar &&'
        " pushd {params.outdir} &&"
        " rename 's/^ar.([\w]+)(.paired)?.fastq$/$1.fastq/' ar*.fastq &&"
        " mv unidentified.fastq unknown.fastq &&"
        " popd"
        " && tar cvf {output} {params.outdir}"
        " && rm -vf {params.outdir}/*.fastq"
        ") >{log} 2>&1"

rule flexbar_dm:
    input:
        kf='keyfiles/{barcode}_flexbar.fasta',
        reads="data/reads/{seed}_{barcode}.fastq.gz",
    output:
        "data/demuxed/flexbar/{seed}_{barcode}/reads.tar"
    params:
        outdir=lambda w: 'data/demuxed/flexbar/{}_{}/'.format(w.seed, w.barcode),
    log:
        'data/log/flexbar/{seed}_{barcode}.log'
    benchmark:
        'data/benchmarks/flexbar/{seed}_{barcode}.txt'
    threads: 1
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
        " rename 's/^.*_barcode_([\w]+.fastq)$/\\1/' *.fastq &&"
        " mv unassigned.fastq unknown.fastq &&"
        " popd"
        " && tar cvf {output} {params.outdir}"
        " && rm -vf {params.outdir}/*.fastq"
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
    params:
        nread=config.get('num_reads', 10000),
        readlen=config.get('read_length', 10),
    shell:
        "mason_simulator --seed {wildcards.seed}"
        "   --illumina-read-length {params.readlen}"
        "   -n {params.nread}"
        "   -ir {input.gen}"
        "   -o {output.r1}"
        "   -or {output.r2} >{log} 2>&1"


rule randgenome:
    output:
        "data/genomes/{seed}.fa"
    params:
        genome_size=config.get('genome_size', 100000)
    log:
        "data/log/genome/{seed}.log"
    shell:
        "scripts/randfa.py"
        "   {params.genome_size}"
        "   {wildcards.seed}"
        " > {output}"
        " 2>{log}"
