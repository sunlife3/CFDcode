CC   = gfortran
EULERSRC = subprog.f90 plot.f90 InitCond.f90 BdyCond.f90  evalFlux.f90  AUSMDV.f90 HLL.f90 main.f90
NSSRC    = subprog.f90 plot.f90 InitCond.f90 BdyCond.f90 evalFlux.f90  viscFlux.f90 AUSMDV.f90 HLL.f90 main.f90 
EULEROBJ = $(EULERSRC:%.f90=%.o)
NSOBJ = $(NSSRC:%.f90=%.o)
TARGET = a.out

.PHONY: buildeuler buildNS 
buildeuler:
	$(CC) -c $(EULERSRC)
	$(CC) -o $(TARGET) $(EULEROBJ)

clean:
	rm *.o
	mv *.json result/
	mv *.dat result/


buildNS:
	$(CC) -c $(NSSRC)
	$(CC) -o $(TARGET) $(NSOBJ)