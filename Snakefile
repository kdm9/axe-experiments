GENOMES = [
    'athaliana',
    'ecoli',
    'scerevesiae',
]
NREADS = 1000000
SEED = 198712
BARCODE_SETS = [
    'gbs',
    'gbs-se',
    'nested',
]

rule all:
    input:
        expand("data/genomes/{gen}.fa", gen=GENOMES),
        expand("data/reads/{gen}_{barcode}.fastq.gz", gen=GENOMES, barcode=BARCODE_SETS)

rule clean:
    shell:
        'rm -rf data/'

rule reads:
    input:
        gen="data/genomes/{genome}.fa",
        barcode="{barcode}.axe",
    output:
        "data/reads/{genome}_{barcode}.fastq.gz"
    params:
        tmpbase=lambda w: "/dev/shm/tmp_{}_{}".format(w.genome, w.barcode)
    shell:
        "mason_simulator --seed {SEED}"
        "   -n {NREADS}"
        "   -ir {input.gen}"
        "   -o {params.tmpbase}_R1.fq"
        "   -or {params.tmpbase}_R2.fq &&"
        "( ./add-barcodes.py -s {SEED}"
        "   {input.barcode}"
        "   {params.tmpbase}_R1.fq"
        "   {params.tmpbase}_R2.fq |"
        "  gzip -4 >{output} ) && "
        "rm -f {params.tmpbase}_R?.fq"


rule ath:
    output:
        "data/genomes/athaliana.fa"
    shell:
        "rm -f TAIR10_chr*.fas* && "
        "wget -nv ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr1.fas "
        "         ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr2.fas "
        "         ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr3.fas "
        "         ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr4.fas "
        "         ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr5.fas "
        "         ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chrC.fas "
        "         ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chrM.fas && "
        "cat TAIR10_chr*.fas >{output} && "
        "rm -f TAIR10_chr*.fas*"

rule ecoli:
    output:
        "data/genomes/ecoli.fa"
    shell:
        "rm -f Escherichia_coli_*.fa* &&"
        "wget -nv ftp://ftp.ensemblgenomes.org/pub/bacteria/release-30/fasta/bacteria_9_collection/escherichia_coli_str_k_12_substr_dh10b/dna/Escherichia_coli_str_k_12_substr_dh10b.GCA_000019425.1.30.dna.genome.fa.gz &&"
        "zcat Escherichia_coli_*.fa.gz >{output} &&"
        "rm -f Escherichia_coli_*.fa.gz*"

rule yeast:
    output:
        "data/genomes/scerevesiae.fa"
    shell:
        "rm -f S288C_referenc*.tgz* &&"
        "wget -nv http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz &&"
        "tar xvf S288C_reference_genome_Current_Release.tgz &&"
        "mv S288C_reference_genome*/S288C_reference_sequence*.fsa {output} &&" 
        "rm -rf S288C_referenc*.tgz* S288C_reference_genome*"
