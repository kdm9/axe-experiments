GENOMES = [
    #'athaliana',
    'ecoli',
    #'scerevesiae',
]
# FIXME: Increase this to something more reasonable for any publication version
NREADS = 1000
SEED = 198712
READLEN = 20
BARCODE_SETS = [
#    'gbs',
#    'gbs-se',
    '8bp-se',
]
DEMUXERS = [
    'axe',
    'flexbar',
]

rule all:
    input:
        expand("data/reads/{gen}_{barcode}.fastq.gz", gen=GENOMES,
               barcode=BARCODE_SETS),
        #expand("data/demuxed.{demuxer}/{gen}_{barcode}_unknown_il.fastq", gen=GENOMES,
        #       barcode=BARCODE_SETS, demuxer=DEMUXERS)
        #expand("data/assessments/{gen}_{barcode}_{demuxer}.tab", gen=GENOMES,
        #       barcode=BARCODE_SETS, demuxer=DEMUXERS)


rule clean:
    shell:
        'rm -rf data/reads'

rule axe_combo:
    input:
        kf='keyfiles/{barcode}.axe',
        reads="data/reads/{genome}_{barcode}.fastq.gz",
    output:
        'data/demuxed.axe/{genome}_{barcode}/{genome}_{barcode}_unknown_il.fastq'
    params:
        outdir=lambda w: 'data/demuxed.axe/{}_{}'.format(w.genome, w.barcode),
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
        kf='keyfiles/{barcode}.flexbar',
        r1="data/reads/{genome}_{barcode}_R1.fastq.gz",
        r2="data/reads/{genome}_{barcode}_R2.fastq.gz",
    output:
        'data/demuxed.flexbar/{genome}_{barcode}_unknown_il.fastq'
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
        'data/demuxed.flexbar/{genome}_{barcode}_unknown_il.fastq'
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
        "(./add-barcodes.py -s {SEED}"
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
