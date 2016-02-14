GENOMES = athaliana ecoli scerevesiae

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
