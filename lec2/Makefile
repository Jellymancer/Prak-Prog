CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
RM = -rm

all: ima.txt real.txt opg2.txt
	cat real.txt
	cat ima.txt
	cat opg2.txt
ima.txt: opg1ima
	./opg1ima > ima.txt

real.txt: opg1
	./opg1 > real.txt

opg2.txt: opg2
	./opg2 > opg2.txt

opg1: opg1.o

opg1ima: opg1ima.o

opg2: opg2.o

.PHONY:clean
clean:
	$(RM) -f opg1.o opg1ima.o opg1 opg1ima real.txt ima.txt
