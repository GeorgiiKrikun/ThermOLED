.PHONY: all execute clean cleanbackup

all: main

constants.o: constants.f90
	mpif90 -ffpe-trap=invalid -fcheck=all -c constants.f90 -o constants.o

SGI.o: SGI.f90 constants.o
	mpif90 -ffpe-trap=invalid -fcheck=all -c SGI.f90 -o SGI.o

PoisonHeatIteration.o: PoisonHeatIteration.f90 constants.o
	mpif90 -c -ffpe-trap=invalid -fcheck=all PoisonHeatIteration.f90 -o PoisonHeatIteration.o

main.o: main.f90 constants.o PoisonHeatIteration.o SGI.o
	mpif90 -ffpe-trap=invalid -fcheck=all -c main.f90 -o main.o

main: main.o constants.o PoisonHeatIteration.o SGI.o
	mpif90 -ffpe-trap=invalid -fcheck=all main.o constants.o PoisonHeatIteration.o SGI.o -o thermoled

clean: 
	rm -fv *.o *.mod thermoled

cleanbackup:
	sh scripts/clearbackup.sh

