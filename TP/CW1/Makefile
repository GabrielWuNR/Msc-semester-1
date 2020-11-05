
CC=     icc -O3 -qopenmp -std=c99
LIB=    -lm

OBJ= loop.o

#
# Compile
#
loop:   $(OBJ)
	$(CC) -o $@ $(OBJ) $(LIB)

.c.o:
	$(CC) -c $<

#
# Clean out object files and the executable.
#
clean:
	rm *.o loop

