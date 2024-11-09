objects := dump_psi.o propagate.o scalar.o ham_psi.o \
	segal_bargmann.o propagate_green.o propagate_ab.o\
	propagate_trap.o variance.o propagate_central.o pval.o wigner.o\
	propagate_convert.o pvar.o xval.o
        #propagate_fft.o\
sources := $(objects:.o=.f90)
progname := tdse
progs := $(progname).f90
modules := #mkl_dfti.mod
modsources := $(modules:.mod=.f90)
FFLAGS := -fdefault-integer-8  -m64  -I"${MKLROOT}/include"
LDFLAGS := -m64  -Wl,--start-group ${MKLROOT}/lib/libmkl_gf_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

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

$(progname): $(objects) 
	gfortran $(progs) $^ $(LDFLAGS) -o $(progname)

$(objects): $(modules) $(sources)
	gfortran $(FFLAGS) -c $(subst .o,.f90,$@) -o $@

$(modules): $(modsources)
	gfortran $(FFLAGS) -c $(subst .mod,.f90,$@)

clean:
	rm -f $(objects) $(progname) $(modules)
