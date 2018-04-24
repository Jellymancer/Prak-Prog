CFLAGS += -Wall -std=c11 `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

all:out.txt
	cat $<

out.txt: main
	./$< > $@

main:main.o newton.o print.o func.o

clean:
	rm -f main *.o out*
