CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)



print: main
	./main > out.txt
	cat out.txt

main: main.o qrdecomp.o myprint.o qrgssolve.o


.PHONEY: clean

clean:
	rm *.o *.txt
