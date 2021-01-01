MF=	Makefile

CC=mpicc -cc=icc
CFLAGS=	-O3

LFLAGS= $(CFLAGS)

EXE=	percolate

INC= \
	src/src_h/percolate.h\
	src/src_h/Changevalue.h\
	src/src_h/parallel.h

SRC= \
	src/percolate.c \
	src/percwritedynamic.c \
	src/uni.c\
	src/ChangeGridValue.c\
	src/parallelpart.c
SRC0= \
	percolate.o \
	percwritedynamic.o \
	uni.o\
	ChangeGridValue.o\
	parallelpart.o

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(OBJ):	$(INC)

$(EXE):	$(OBJ)
	$(CC) $(LFLAGS) -o $@ $(SRC0)

$(OBJ):	$(MF)

clean:
	rm -f $(EXE) $(SRC0) core
	rm -f map.pgm
