CFLAGS += -Wall -std=c11 `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
QRDIR = ../qr.gramschmidt.gsl
CFLAGS += -I$(QRDIR)

all:out.txt
	cat $<

out.txt: main
	./$< > $@

main:main.o newton.o print.o

clean:
	rm -f main *.o out*
