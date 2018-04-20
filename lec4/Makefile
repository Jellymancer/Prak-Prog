CFLAGS = -Wall -std=gnu99 -O
LDLIBS = -lm

.PHONEY: all clean

all: out.txt
	cat $<

out.txt: main
	./$< > $@


main: main.o komplex.o


main.o komplex.o: komplex.h



clean: 
	$(RM) main main.o komplex.o out.txt
