buildeuler: InitCond.f90 BdyCond.f90 subprog.f90 evalFlux.f90 SD-SLAU.f90 AUSMDV.f90 HLL.f90 main.f90 \
		InitCond.o BdyCond.o subprog.o SD-SLAU.o AUSMDV.f90 HLL.o evalFlux.o plot.o main.o
	gfortran -c subprog.f90
	gfortran -c plot.f90
	gfortran -c main.f90
	gfortran -c InitCond.f90
	gfortran -c BdyCond.f90
	gfortran -c SD-SLAU.f90
	gfortran -c AUSMDV.f90
	gfortran -c HLL.f90
	gfortran -c evalFlux.f90
	gfortran -o a.exe InitCond.o BdyCond.o subprog.o SD-SLAU.o AUSMDV.o HLL.o evalFlux.o plot.o main.o

buildNS: InitCond.f90 BdyCond.f90 subprog.f90 evalFlux.f90 viscFlux.f90 SD-SLAU.f90 AUSMDV.f90 HLL.f90 main.f90 \
	 InitCond.o BdyCond.o subprog.o evalFlux.o viscFlux.o SD-SLAU.o AUSMDV.o HLL.o main.o
	ifort -c main.f90
	ifort -c InitCond.f90
	ifort -c BdyCond.f90
	ifort -c subprog.f90
	ifort -c SD-SLAU.f90
	ifort -c AUSMDV.f90
	ifort -c HLL.f90
	ifort -c evalFlux.f90
	ifort -c viscFlux.f90
	ifort -c plot.f90
	ifort -o solver InitCond.o BdyCond.o subprog.o SD-SLAU.o AUSMDV.o HLL.o evalFlux.o viscFlux.o plot.o main.o	

ramprun: 
	b.exe ramp

Bowrun: 
	a.exe Bow

SNSrun: 
	a.exe SNS 

SOrun: 
	a.exe SO 

SWBLIrun:
	./solver SWBLI > log.dat &  

platerun: 
	a.exe visplate
