CFLAGS += -Wall -std=c11 `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

all: out.txt as1.svg plotairy.svg
	cat out.txt
	xdg-open as1.svg
	xdg-open plotairy.svg

as1.svg: out.txt Makefile
	echo '\
        set term svg background rgb "white";\
        set out "$@";\
        set xlabel "x";\
        set ylabel "y";\
        set tics out;\
        set key out;\
        set title "Solution to dydx=3x+4" ;\
        plot "as1.txt" index 0 with lines title "y(x)";\
        '|gnuplot

plotairy.svg: out.txt Makefile
	echo '\
        set term svg size 1200,800 background rgb "white";\
        set out "$@";\
        set xlabel "x";\
        set ylabel "y";\
        set tics out;\
        set key out;\
        set title "Airy equation solution" ;\
        plot "airy.txt" using 1:2 with lines title "Numeric-airy1", "airy.txt" using 1:3 with lines title "Numeric-airy2",\
	"GSLairy.txt" using 1:2 with points title "GSL-airy1","GSLairy.txt" using 1:3 with points title "GSL-airy";\
        '|gnuplot

out.txt: main
	./main > out.txt

main: main.o print.o rkstepX.o driver.o


clean:
	rm -f main *.o out*
