objects := dump_psi.o propagate.o scalar.o ham_psi.o propagate_fft.o\
	segal_bargmann.o helloworld.o propagate_green.o propagate_ab.o\
	propagate_trap.o variance.o propagate_central.o pval.o wigner.o\
	propagate_convert.o pvar.o xval.o
sources := dump_psi.f90 propagate.f90 scalar.f90 ham_psi.f90 propagate_fft.f90\
	segal_bargmann.f90 helloworld.f90 propagate_green.f90 propagate_ab.f90\
	propagate_trap.f90 variance.f90 propagate_central.f90 pval.f90 wigner.f90\
	propagate_convert.f90 pvar.f90 xval.f90
progs := tdse.f90
progname := tdse
MKLROOT := /modfac/apps/Intel/compilers_and_libraries_2020.1.217/linux/mkl

# You want latexmk to *always* run, because make does not have all the info.
# Also, include non-file targets in .PHONY so they are run regardless of any
# file of the given name existing.
.PHONY: $(objects) $(progs) $(progname) all clean

# The first rule in a Makefile is the one executed by default ("make"). It
# should always be the "all" rule, so that "make" and "make all" are identical.
all: $(progname) $(objects) 

# CUSTOM BUILD RULES

# In case you didn't know, '$@' is a variable holding the name of the target,
# and '$<' is a variable holding the (first) dependency of a rule.
# "raw2tex" and "dat2tex" are just placeholders for whatever custom steps
# you might have.

# %.tex: %.dat
# 	./dat2tex $< > $@

# MAIN LATEXMK RULE

# -pdf tells latexmk to generate PDF directly (instead of DVI).
# -pdflatex="" tells latexmk to call a specific backend with specific options.
# -use-make tells latexmk to call make for generating missing files.

# -interaction=nonstopmode keeps the pdflatex backend from stopping at a
# missing file reference and interactively asking you for an alternative.

$(progname): $(progs) $(objects) 
	gfortran $< -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -o $(progname)

$(objects): $(sources)
	gfortran -c $< -o $@
	@echo made

clean:
	rm -f *.o $(progname)
