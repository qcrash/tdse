bin_PROGRAMS = tdse
tdse_SOURCES = tdse.f90 dump_psi.f90 propagate.f90 scalar.f90 ham_psi.f90 \
	segal_bargmann.f90 propagate_green.f90 propagate_ab.f90\
	propagate_trap.f90 variance.f90 propagate_central.f90 pval.f90 wigner.f90\
	propagate_convert.f90 pvar.f90 xval.f90 getcli.f90
sources := $(objects:.o=.f90)
progname := tdse
progs := $(progname).f90
modules := #mkl_dfti.mod
modsources := $(modules:.mod=.f90)
FFLAGS := -fdefault-integer-8  -m64  -I"${MKLROOT}/include"
LDFLAGS := -fdefault-integer-8 -m64  -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

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
