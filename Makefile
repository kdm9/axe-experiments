GENOMES = athaliana ecoli scerevesiae
NREADS = 1000000
SEED = 198712
BARCODE_SETS = gbs gbs_se nested

.PHONY: barcode_reads
prebarcode_reads: $(foreach gen,$(GENOMES),data/reads/$(gen)_barcode.fastq.gz)

.PHONY: prebarcode_reads
prebarcode_reads: $(foreach gen,$(GENOMES),data/reads/$(gen)_prebarcode.fastq.gz)

data/reads/%.fastq.gz: data/genomes/athaliana.fa
	mason_simulator --seed $(SEED) \
		-n $(NREADS) \
	       	-ir $^ \
		-o data/tmp_$*_R1.fq \
		-or data/tmp_$*_R2.fq
	./add-barcodes.py -s $(SEED) \

		data/tmp_$*_R1.fq \
		data/tmp_$*_R2.fq | \
		gzip -4 >$@
	rm -f data/tmp_$*_R?.fq


.PHONY: genomes
genomes: $(foreach gen,$(GENOMES),data/genomes/$(gen).fa)

data/genomes/athaliana.fa: | dirs
	rm -f TAIR10_chr*.fas*
	wget -nv $(foreach i,1 2 3 4 5 C M,ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr$(i).fas)
	cat TAIR10_chr*.fas >$@
	rm -f TAIR10_chr*.fas*

data/genomes/ecoli.fa: | dirs
	rm -f Escherichia_coli_str_k_12_substr_dh10b.GCA_000019425.1.30.dna.genome.fa.gz*
	wget -nv ftp://ftp.ensemblgenomes.org/pub/bacteria/release-30/fasta/bacteria_9_collection/escherichia_coli_str_k_12_substr_dh10b/dna/Escherichia_coli_str_k_12_substr_dh10b.GCA_000019425.1.30.dna.genome.fa.gz
	gunzip Escherichia_coli_str_k_12_substr_dh10b.GCA_000019425.1.30.dna.genome.fa.gz
	mv Escherichia_coli_str_k_12_substr_dh10b.GCA_000019425.1.30.dna.genome.fa $@
	rm -f Escherichia_coli_str_k_12_substr_dh10b.GCA_000019425.1.30.dna.genome.fa.gz*

data/genomes/scerevesiae.fa: | dirs
	rm -f S288C_referenc*.tgz*
	wget -nv http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz
	tar xvf S288C_reference_genome_Current_Release.tgz
	mv S288C_reference_genome*/S288C_reference_sequence*.fsa $@
	rm -rf S288C_referenc*.tgz* S288C_reference_genome*

dirs:
	mkdir -p data/genomes
	mkdir -p data/reads
