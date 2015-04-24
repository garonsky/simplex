dlvl = ./.
include $(dlvl)/../Makefile.in


all: simplex

simplex: 
	$(CC) -c $(CFLAGS) -I../src simplex2.c
	$(LOADER) -Wall -g -o simplex simplex2.o $(CBLIB) $(BLLIB)

cleanall:
	rm -f *.o simplex 
