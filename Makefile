SQLITE=no
OMP=no
DEBUG=no
PROG=euclid
ODIR=src
OBJS=main.o pari_min_c.o graph.o basef.o
LDFLAGS=-lpari -lm -ldl
CFLAGS=-Wall -fcommon
ifeq "$(DEBUG)" "yes"
CFLAGS+=-g -pedantic -O0
else
CFLAGS+=-O3
endif
ifeq "$(SQLITE)" "yes"
OBJS+= datah.o sqlite3.o
CFLAGS+=-D SQLITE
LDFLAGS+= -lpthread
ifeq "$(OMP)" "yes"
CFLAGS+=-fopenmp -D OPENMP
endif
else
ifeq "$(OMP)" "yes"
CFLAGS+=-fopenmp -D OPENMP
LFLAGS+=-lpthread
endif
endif

OBJ = $(patsubst %,$(ODIR)/%,$(OBJS))


euclid: $(OBJ)
	gcc $(CFLAGS) -o $(PROG) $(OBJ) $(LDFLAGS)

$(ODIR)/%.o: %.c header.h
	gcc $(CFLAGS) -o $@ $< -c

%.c: %.gp
	gp2c -S -g > $@ $<

clean: clean-chapel
	rm -f $(ODIR)/*.o $(PROG) result

tags:
	ctags *.c *.h

# ======================================================================
# Chapel implementation (see CHAPEL.md): build and test automation.
# The C tool above (`make euclid`, the default target) is unaffected and
# remains this Makefile's ground-truth oracle for `make test`.
# ======================================================================
CHPL=chpl
CHPL_BIN=euclid_chpl
CHPL_DEBUG_BIN=euclid_chpl_debug
# --fast disables Chapel's runtime bounds/overflow checking in exchange for
# a large speedup; use `make chapel-debug` (no --fast) when tracking down a
# correctness bug, since that keeps the checks enabled and gives better
# error messages/backtraces.
CHPL_FLAGS=--fast
CHPL_MODULES=NumberField.chpl SmallElements.chpl Sieve.chpl Certify.chpl
TEST_BUILD_DIR=tests/build
SMOKE_TESTS=test_numberfield test_smallelements test_sieve test_sieve_debug test_certify

# name/polynomial pairs used by `make fixtures` (kept in sync with the
# FIELDS list in tests/validate.sh). Quoted and fed through `read` below
# rather than iterated as shell words, so the literal '*' and '^' in some
# polynomials are never glob-expanded.
define FIELD_LIST
x2-2      x^2-2
x2-61     x^2-61
x3+x2-1   x^3+x^2-1
x3-3x-1   x^3-3*x-1
x4-4x2+2  x^4-4*x^2+2
x5-x-1    x^5-x-1
endef
export FIELD_LIST

.PHONY: all chapel chapel-debug smoke test check fixtures clean-chapel

all: euclid chapel

# Optimized build (default): what you want for actually running the sieve.
chapel: $(CHPL_BIN)

# Force the standard CPU ("flat") locale model for this project's builds.
# If CHPL_LOCALE_MODEL=gpu is set in the ambient environment (e.g. a Chapel
# toolchain configured for GPU offload work), chpl 2.9.0 hits an internal
# compiler error (COD-CG--XPR-03263) analyzing Sieve.chpl's forall loops for
# GPU eligibility; this code doesn't use any GPU-specific features, so
# there's nothing to gain from that mode here, only risk.
CHPL_ENV=CHPL_LOCALE_MODEL=flat

$(CHPL_BIN): main.chpl $(CHPL_MODULES)
	$(CHPL_ENV) $(CHPL) $(CHPL_FLAGS) -M . main.chpl -o $(CHPL_BIN)

# Unoptimized build: bounds/overflow checks stay enabled; use this while
# debugging correctness issues. Produces ./euclid_chpl_debug so it can
# coexist with the optimized ./euclid_chpl.
chapel-debug: $(CHPL_DEBUG_BIN)

$(CHPL_DEBUG_BIN): main.chpl $(CHPL_MODULES)
	$(CHPL_ENV) $(CHPL) -M . main.chpl -o $(CHPL_DEBUG_BIN)

$(TEST_BUILD_DIR):
	mkdir -p $(TEST_BUILD_DIR)

$(TEST_BUILD_DIR)/%: tests/%.chpl $(CHPL_MODULES) | $(TEST_BUILD_DIR)
	$(CHPL_ENV) $(CHPL) -M . $< -o $@

# Compiles and runs the standalone per-module smoke tests under tests/.
smoke: $(addprefix $(TEST_BUILD_DIR)/,$(SMOKE_TESTS))
	@for t in $(SMOKE_TESTS); do \
	  echo "=== $$t ==="; \
	  ./$(TEST_BUILD_DIR)/$$t || exit 1; \
	done

# Generates any tests/fixtures/*.txt that don't already exist yet, via the
# one-time Pari/gp field setup step (requires `gp` on PATH).
fixtures:
	@echo "$$FIELD_LIST" | while read -r name poly; do \
	  [ -z "$$name" ] && continue; \
	  if [ ! -f tests/fixtures/$$name.txt ]; then \
	    echo "generating tests/fixtures/$$name.txt from $$poly"; \
	    POLY="$$poly" gp -q generate_field.gp >/dev/null && cp field_data.txt tests/fixtures/$$name.txt; \
	  fi; \
	done

# Runs the Chapel pipeline against the reference C tool over tests/fixtures/.
test: euclid chapel fixtures
	tests/validate.sh

check: smoke test

clean-chapel:
	rm -rf $(CHPL_BIN) $(CHPL_BIN)_real $(CHPL_DEBUG_BIN) $(CHPL_DEBUG_BIN)_real $(TEST_BUILD_DIR) *.dSYM
