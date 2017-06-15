
.PHONY: all
all: plot_results.html

%.html: %.Rmd ../data/stats/accuracy_summary.tsv ../data/stats/timing_summary.tsv
	Rscript -e 'rmarkdown::render("$<", "html_document")'
