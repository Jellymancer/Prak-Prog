#CC = tcc
CC = gcc -std=c99
CFLAGS += -Werror -O2
LDLIBS += -lm

all: out.A.data.txt
	cat out.A.data.txt

out.A.data.txt : main ; ./$< > $@

main: main.o ode_driver.o rkstep12.o rkstep23.o

clean : ; rm -f main *.o *.png out*
