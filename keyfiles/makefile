SETS = 8bp 4bp
EXTS = _flexbar.fasta .sabre .adapterremoval
EXTS_PE = .adapterremoval

SE = $(foreach ext,$(EXTS),$(foreach set,$(SETS),$(set)-se$(ext)))
PE = $(foreach ext,$(EXTS_PE),$(foreach set,$(SETS),$(set)-pe$(ext)))
all: $(SE) $(PE)

%-se_flexbar.fasta: %-se.axe
	grep -v '^[#;]' $< | grep -iv '^Barcode' | \
	   	awk '{printf(">%s\n%s\n", $$2, $$1);}' > $@

%-se.sabre: %-se.axe
	grep -v '^[#;]' $< | grep -iv '^Barcode' | \
	   	awk '{printf("%s\t%s_R1.fastq\t%s_R2.fastq\n", $$1, $$2, $$2);}' > $@

%-se.adapterremoval: %-se.axe
	grep -v '^[#;]' $< | grep -iv '^Barcode' | \
		awk '{printf("%s %s\n", $$2, $$1);}' > $@

%-pe.adapterremoval: %-pe.axe
	grep -v '^[#;]' $< | grep -iv '^Barcode' | \
	   	awk '{printf("%s %s %s\n", $$3, $$1, $$2);}' > $@
