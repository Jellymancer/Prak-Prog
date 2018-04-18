CFLAGS = -Wall -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS = `gsl-config --libs`
ECHO = /bin/echo


show: ros.txt
	cat ros.txt

ros.txt: ros
	./ros > ros.txt

ros: ros.c
	gcc ros.c -o ros $(CFLAGS) $(LDLIBS)


.PHONEY:clean

clean:
	rm  ros ros.txt
