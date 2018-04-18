CFLAGS += -Wall -std=gnu1x
LDLIBS += $$(gsl-config --libs) 


all: main
	xdg-open airy.svg
	cat out.txt
main: airy.gpi airy.data linalg
	gnuplot $<
	./lin > out.txt

airy.data: airy
	./airy > airy.data



airy: airy.o airy.c
	gcc airy.c -o airy $(LDLIBS) 


linalg: linalg.c linalg.o
	gcc linalg.c -o lin $(LDLIBS)

.PHONEY: clean

clean:
	rm *.data *.o airy lin *.svg *.txt
