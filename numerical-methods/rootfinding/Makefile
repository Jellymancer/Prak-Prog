CFLAGS += -Wall -std=c99 `gsl-config --cflags`
LDLIBS += `gsl-config --libs` -lm

all:out.txt
	cat $<

out.txt: main
	./$< > $@

main:main.o newton.o print.o func.o GSL.o

clean:
	rm -f main *.o out*
