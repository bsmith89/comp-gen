PROJ_DIRS = $(shell find . \( -name ".git" \) -prune -o -type d -print)

all:   docs figs
figs:

# Don't delete any intermediate files
.SECONDARY:

.PHONY: all docs figs



# --------------
#  Data Recipes
# --------------

# This is very not-editable, since it's just line ranges for the mapping
# and line ranges for the tree.
meta/bacteria.names.tsv: raw/RAxML_3473_1_flux_trim4_10.19.14.nex
	cat $> \
		| sed -n '18,1164p' \
		| awk '{print $$1"\t"$$2}' \
		| sed 's:[,;]::' \
		> $@

tre/bacteria.nwk: raw/RAxML_3473_1_flux_trim4_10.19.14.nex
	cat $> \
		| sed -n '1165p' \
		| cut -d ' ' -f 5 \
		> $@

tre/bacteria.names.nwk: scripts/rename_tree.py tre/bacteria.nwk \
                        meta/bacteria.names.tsv
	$^ > $@

meta/K0.tsv: ./scripts/clean_tables.py ./raw/2.6.15_KEGG_K0.csv
	$^ > $@

meta/M0.tsv: ./scripts/clean_tables.py ./raw/2.6.15_KEGG_M0.csv
	$^ > $@

# -----------------------
#  Analysis Recipes
# -----------------------

res/%.pairwise_log16S_corr.tsv: scripts/pairwise_log16S_corr.py \
                                tre/bacteria.names.nwk meta/%.tsv
	$^ > $@

# -----------------------
#  Documentation Recipes
# -----------------------

ALL_DOCS_MD = $(foreach d,${PROJ_DIRS}, $(wildcard ${d}/*.md))
ALL_DOCS_HTML = $(patsubst %.md,%.html, ${ALL_DOCS_MD})

docs: ${ALL_DOCS_HTML}

pandoc_recipe_md2html = \
cat $< \
| pandoc -f markdown -t html5 -s \
			--highlight-style pygments --mathjax \
			--toc --toc-depth=4 \
			--css static/main.css \
> $@
%.html: %.md
	${pandoc_recipe_md2html}

# ------------------------
#  Initialization Recipes
# ------------------------
SEMAPHORE = .initialized

ifeq (,$(wildcard ${SEMAPHORE}))  # The semaphor file does NOT exist
touch_semaphore:
	# Create the initialization semaphore.
	touch ${SEMAPHORE}
base_init: touch_semaphore
	git submodule update --init --recursive
	# Configure IPYNB output filtering
	git config --local filter.dropoutput_ipynb.clean scripts/utils/ipynb_output_filter.py
	git config --local filter.dropoutput_ipynb.smudge cat
init_from_project: base_init
init_from_template: base_init
	# Remove the template remote
	git remote remove origin
	# Link README to project notes, instead of template notes.
	unlink README.md
	ln -s NOTE.md README.md
	# Remove all of the commits after the first, leaving files intact,
	# add files created/changed during initialization,
	# and amend the first commit with everything else.
	# TL;DR Squash everything to a single first commit.
	git reset --soft $(git rev-list --max-parents=0 HEAD)
	git add -A
	git commit --amend -em "Clean project.  Let's get started!"
else
base_init:
	$(error "This directory contains a file '${SEMAPHORE}' indicating \
			 that it has already been initialized. \
			 Delete this file if you would like to run the initialization \
			 scripts again (which may cause problems).")
init_from_project: base_init
init_from_template: base_init
endif

# -----------------
#  Cleanup Recipes
# -----------------

clean:
	rm -f ${ALL_DOCS_HTML}
