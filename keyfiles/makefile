SETS = 8bp-se gbs-se gbs nested
EXTS = flexbar

all: $(foreach ext,$(EXTS),$(foreach set,$(SETS),$(set).$(ext)))

%.flexbar: %.axe
	grep -v '^[#;]' $< | grep -iv '^Barcode' | awk '{printf(">%s\n%s\n", $$2, $$1);}' > $@

