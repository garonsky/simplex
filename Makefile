dlvl = ./.
include $(dlvl)/../Makefile.in


all: simplex

simplex: 
	$(CC) -c $(CFLAGS) -I../src simplex2.c
	$(LOADER) -g -o simplex simplex.o $(CBLIB) $(BLLIB)

cleanall:
	rm -f *.o simplex 
