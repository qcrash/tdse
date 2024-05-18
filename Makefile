lectures := $(patsubst %.tex,%.pdf,$(wildcard lecture*.tex))
assignments := $(patsubst %.tex,%.pdf,$(wildcard assignment*.tex))
# You want latexmk to *always* run, because make does not have all the info.
# Also, include non-file targets in .PHONY so they are run regardless of any
# file of the given name existing.
.PHONY: $(lectures) $(assignments) $(handout) syllabus.pdf all clean

# The first rule in a Makefile is the one executed by default ("make"). It
# should always be the "all" rule, so that "make" and "make all" are identical.
all: $(lectures) $(assignments) $(handout) syllabus.pdf

# CUSTOM BUILD RULES

# In case you didn't know, '$@' is a variable holding the name of the target,
# and '$<' is a variable holding the (first) dependency of a rule.
# "raw2tex" and "dat2tex" are just placeholders for whatever custom steps
# you might have.

%.tex: %.raw
	./raw2tex $< > $@

%.tex: %.dat
	./dat2tex $< > $@

# MAIN LATEXMK RULE

# -pdf tells latexmk to generate PDF directly (instead of DVI).
# -pdflatex="" tells latexmk to call a specific backend with specific options.
# -use-make tells latexmk to call make for generating missing files.

# -interaction=nonstopmode keeps the pdflatex backend from stopping at a
# missing file reference and interactively asking you for an alternative.

$(lectures): %.pdf: %.tex
	latexmk -lualatex -shell-escape $<
	latexmk -bibtex- -lualatex -shell-escape \
	-usepretex='\PassOptionsToClass{handout}{beamer}' \
	-outdir=handout $<  

$(assignments): %.pdf: %.tex
	latexmk -lualatex -shell-escape $<

syllabus.pdf: syllabus.tex
	latexmk -lualatex $<

clean:
	latexmk -C
	latexmk -C -outdir=handout
	rm -f lecture*.nav lecture*.snm handout/*.nav handout/*.snm
