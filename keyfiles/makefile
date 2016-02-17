SETS = 8bp-se gbs-se gbs nested
EXTS = flexbar sabre

all: $(foreach ext,$(EXTS),$(foreach set,$(SETS),$(set).$(ext)))

%.flexbar: %.axe
	grep -v '^[#;]' $< | grep -iv '^Barcode' | awk '{printf(">%s\n%s\n", $$2, $$1);}' > $@

%.sabre: %.axe
	grep -v '^[#;]' $< | grep -iv '^Barcode' | awk '{printf("%s\t%s_R1.fastq\t%s_R2.fastq\n", $$1, $$2, $$2);}' > $@
