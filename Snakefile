GENOMES = [
    'random',
]
# FIXME: Increase this to something more reasonable for any publication version
NREADS = 10000
SEED = 198712
READLEN = 20
BC_DM = [
    ('4bp-se', 'axe'),
    ('nested-se', 'axe'),
    ('4bp-pe', 'axe'),
    ('nested-pe', 'axe'),
]
BARCODE_SETS = list(set([bd[0] for bd in BC_DM]))

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

rule all:
    input:
        expand("data/reads/{gen}_{barcode}.fastq.gz", gen=GENOMES,
               barcode=BARCODE_SETS),
        ["data/demuxed.{dm}/{gen}_{bc}/unknown_il.fastq".format(dm=dm, bc=bc, gen=gen)
            for bc, dm in BC_DM
            for gen in GENOMES],
        ["data/assessments/{gen}_{bc}_{dm}.tab".format(dm=dm, bc=bc, gen=gen)
            for bc, dm in BC_DM
            for gen in GENOMES],


rule clean:
    shell:
        'rm -rf data/reads'


rule assess:
    input:
        expand('data/demuxed.{{demuxer}}/{{genome}}_{{barcode}}/{s}_il.fastq',
               s=['A', 'B', 'C', 'D', 'unknown'])
    output:
        'data/assessments/{genome}_{barcode}_{demuxer}.tab'
    shell:
        "scripts/assess-barcodes.py {input} > {output}"


rule axe:
    input:
        kf='keyfiles/{barcode}.axe',
        reads="data/reads/{genome}_{barcode}.fastq.gz",
    output:
        ['data/demuxed.axe/{{genome}}_{{barcode}}/{s}_il.fastq'.format(s=s)
            for s in ['A', 'B', 'C', 'D', 'unknown']]
    params:
        outdir=lambda w: 'data/demuxed.axe/{}_{}/'.format(w.genome, w.barcode),
        combo=lambda w: '-c' if 'pe' in w.barcode else ''
    log:
        'data/log/axe_{genome}_{barcode}.log'
    shell:
        'axe-demux'
        '   {params.combo}'
        '   -b {input.kf}'
        '   -i {input.reads}'
        '   -I {params.outdir}'
        '    >{log} 2>&1'

rule sabre_dm:
    input:
        kf='keyfiles/{barcode}.sabre',
        r1="data/reads/{genome}_{barcode}_R1.fastq.gz",
        r2="data/reads/{genome}_{barcode}_R2.fastq.gz",
    output:
        ['data/demuxed.sabre/{{genome}}_{{barcode}}/{s}_il.fastq'.format(s=s)
            for s in ['A', 'B', 'C', 'D', 'unknown']]
    params:
        outdir='data/demuxed.flexbar',
        combo=lambda w: 'se' if 'se' in w.barcode else '',
    log:
        'data/log/flexbar_{genome}_{barcode}.log'
    shell:
        '(p=$(pwd); '
        't=$(mktemp -d); '
        'pushd $t; '
        '~/prog/bio/forks/sabre/sabre pe '
        '   -f $p/{input.r1}'
        '   -r $p/{input.r2}'
        '   -b $p/{input.kf}'
        '   -u unknown_R1.fastq'
        '   -w unknown_R2.fastq'
        '&& rename "s/^/{wildcards.genome}_{wildcards.barcode}/" *.fastq'
        '&& gzip *.fastq'
        '&& mv *.fastq.gz $p/{params.outdir} )'
        '    >{log} 2>&1'

rule flexbar_dm:
    input:
        kf='keyfiles/{barcode}.flexbar',
        reads="data/reads/{genome}_{barcode}.fastq.gz",
    output:
        ['data/demuxed.flexbar/{{genome}}_{{barcode}}/{s}_il.fastq'.format(s=s)
            for s in ['A', 'B', 'C', 'D', 'unknown']]
    params:
        outdir=lambda w: 'data/demuxed.flexbar/{}_{}'.format(w.genome, w.barcode),
        combo=lambda w: '-c' if 'se' in w.barcode else '',
    log:
        'data/log/flexbar_{genome}_{barcode}.log'
    shell:
        '(flexbar'
        '   -be LEFT'
        '   -b {input.kf}'
        '   -r {input.reads}'
        '   -t {params.outdir} &&'
        'rename "s/_barcode//" {params.outdir}*.fastq )'
        '    >{log} 2>&1'


rule bcdreads:
    input:
        barcode="keyfiles/{barcode}.axe",
        r1="tmp/{genome}_{barcode}-R1.prebcd.fastq",
        r2="tmp/{genome}_{barcode}-R2.prebcd.fastq",
    params:
        re_site=lambda w: '-r TGCAG' if 'gbs' in w.barcode else ''
    output:
        "data/reads/{genome}_{barcode}.fastq.gz"
    log:
        "data/log/bcdreads.{genome}_{barcode}.log"
    shell:
        "(scripts/add-barcodes.py -s {SEED}"
        "   {params.re_site}"
        "   {input.barcode}"
        "   {input.r1}"
        "   {input.r2} |"
        "  gzip -4 >{output}) 2>{log}"


rule reads:
    input:
        gen="data/genomes/{genome}.fa",
        barcode="keyfiles/{barcode}.axe",
    output:
        r1=temp("tmp/{genome}_{barcode}-R1.prebcd.fastq"),
        r2=temp("tmp/{genome}_{barcode}-R2.prebcd.fastq"),
    log:
        "data/log/reads.{genome}_{barcode}.log"
    shell:
        "mason_simulator --seed {SEED}"
        "   --illumina-read-length {READLEN}"
        "   -n {NREADS}"
        "   -ir {input.gen}"
        "   -o {output.r1}"
        "   -or {output.r2} >{log} 2>&1"


rule ath:
    output:
        "data/genomes/athaliana.fa"
    shell:
        "wget -c ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr1.fas "
        "        ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr2.fas "
        "        ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr3.fas "
        "        ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr4.fas "
        "        ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr5.fas "
        "        ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chrC.fas "
        "        ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chrM.fas && "
        "cat TAIR10_chr*.fas >{output} && "
        "rm -f TAIR10_chr*.fas*"

rule ecoli:
    output:
        "data/genomes/ecoli.fa"
    shell:
        "wget -c ftp://ftp.ensemblgenomes.org/pub/bacteria/release-30/fasta/bacteria_9_collection/escherichia_coli_str_k_12_substr_dh10b/dna/Escherichia_coli_str_k_12_substr_dh10b.GCA_000019425.1.30.dna.genome.fa.gz &&"
        "zcat Escherichia_coli_*.fa.gz >{output} &&"
        "rm -f Escherichia_coli_*.fa.gz*"

rule yeast:
    output:
        "data/genomes/scerevesiae.fa"
    shell:
        "wget -c http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz &&"
        "tar xf S288C_reference_genome_Current_Release.tgz &&"
        "mv S288C_reference_genome*/S288C_reference_sequence*.fsa {output} &&"
        "rm -rf S288C_referenc*.tgz* S288C_reference_genome*"

rule randgenome:
    output:
        "data/genomes/random.fa"
    shell:
        "scripts/randfa.py 1000000 {SEED} > {output}"
