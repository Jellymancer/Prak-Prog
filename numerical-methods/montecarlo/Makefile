CFLAGS += -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS += `gsl-config --libs` -lm

all: out.txt plot.svg
	cat out.txt
	xdg-open plot.svg

plot.svg: out.txt Makefile
	echo ' \
	set term svg background rgb "white" size 1200,800; \
	set out "$@"; \
	set grid; \
	set xlabel "N"; \
	set ylabel "error"; \
	set title "Error of the cube integration as function of the number of random points N"; \
	f(x) = a/sqrt(x); \
	FIT_LIMIT = 1e-6; \
	fit f(x) "errcheck.txt" using 1:2 via a; \
	plot "errcheck.txt" index 0 with points title "Error" \
	,a/sqrt(x) with lines title "best fit g(x)=C/sqrt(x)" ; \
	' | gnuplot


out.txt: main
	./main > out.txt

main: main.o mc.o

clean:
	rm -f main *.o out* *.svg *.txt



