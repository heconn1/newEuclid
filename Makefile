SQLITE=no
OMP=no
DEBUG=no
PROG=euclid
ODIR=src
OBJS=main.o pari_min_c.o graph.o basef.o
LDFLAGS=-lpari -lm -ldl
CFLAGS=-Wall
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
	gcc $(CFLAGS) -fcommon -o $@ $< -c

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
CHPL_MODULES=NumberField.chpl SmallElements.chpl Sieve.chpl Certify.chpl
TEST_BUILD_DIR=tests/build
SMOKE_TESTS=test_numberfield test_smallelements test_sieve test_sieve_debug test_certify

# name/polynomial pairs used by `make fixtures` (kept in sync with the
# FIELDS list in tests/validate.sh). Quoted and fed through `read` below
# rather than iterated as shell words, so the literal '*' and '^' in some
# polynomials are never glob-expanded.
define FIELD_LIST
x2-2     x^2-2
x2-61    x^2-61
x3+x2-1  x^3+x^2-1
x3-3x-1  x^3-3*x-1
x5-x-1   x^5-x-1
endef
export FIELD_LIST

.PHONY: all chapel smoke test check fixtures clean-chapel

all: euclid chapel

chapel: $(CHPL_BIN)

$(CHPL_BIN): main.chpl $(CHPL_MODULES)
	$(CHPL) -M . main.chpl -o $(CHPL_BIN)

$(TEST_BUILD_DIR):
	mkdir -p $(TEST_BUILD_DIR)

$(TEST_BUILD_DIR)/%: tests/%.chpl $(CHPL_MODULES) | $(TEST_BUILD_DIR)
	$(CHPL) -M . $< -o $@

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
	rm -rf $(CHPL_BIN) $(CHPL_BIN)_real $(TEST_BUILD_DIR) *.dSYM
