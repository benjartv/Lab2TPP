INC = -I./include
CFLAGS = -fopenmp -Wall -lm
CC = gcc

main: ./build/det.o ./build/detfunct.o
	@ $(CC) -o ./bin/det ./build/det.o ./build/detfunct.o $(CFLAGS)
	@ echo "The compilation has been completed successfully"

build/det.o: ./src/det.c
	@ echo "src/det.c"
	@ $(CC) -c ./src/det.c $(INC) -o ./build/det.o $(CFLAGS)

build/detfunct.o: ./src/detfunct.c
	@ echo "src/detfunct.c"
	@ $(CC) -c ./src/detfunct.c $(INC) -o ./build/detfunct.o $(CFLAGS)

clean:
	@ clear
	@ echo "Cleaning folders.."
	@ rm -rf build/*
	@ rm -rf bin/*
	@ rm -rf output/*

run:
	@ clear
	@ echo "Det OPENMPI"
	@ ./bin/det -i ./input/matriz1212.raw -N 12 -H 2
