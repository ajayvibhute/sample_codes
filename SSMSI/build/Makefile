NVCC = nvcc
CC=mpicc
INCLUDEDIR = ../include
LIBDIR = ../lib
SRCDIR = ../src/
all : svd_fit Main

sim:
	nvcc -c ../src/svd_fit.cu  -o ../src/svd_fit.o -I ../include/ -lm   -D __SIM__
	mpicc  ../src/main.c ../src/svd_fit.o -lcudart -L /usr/local/cuda/lib64/ -I /usr/local/cuda/incl:../include -lm -o ../bin/ssmsi
svd_fit:
	nvcc -c ../src/svd_fit.cu  -o ../src/svd_fit.o -I ../include/   -lm

Main :
	mpicc ../src/main.c ../src/svd_fit.o -lcudart -L /usr/local/cuda/lib64/ -I /usr/local/cuda/incl:../include  -lm -o ../bin/ssmsi 
clean:
	rm svd
	rm gensha

# mpicc  main.c New_Main.o  -lcudart -L /usr/local/cuda/lib64/ -I /usr/local/cuda/incl
