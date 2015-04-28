dlvl = ./.
include $(dlvl)/../Makefile.in


all: simplex

simplex: 
	$(CC) -c $(CFLAGS) -I../src simplex3.c
	$(LOADER) -Wall -g -o simplex simplex3.o $(CBLIB) $(BLLIB)

cleanall:
	rm -f *.o simplex 
