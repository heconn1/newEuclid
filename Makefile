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
	gcc $(CFLAGS) -o $@ $< -c

%.c: %.gp
	gp2c -S -g > $@ $<

clean:
	rm -f $(ODIR)/*.o $(PROG) result

tags:
	ctags *.c *.h
