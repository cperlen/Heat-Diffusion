CXX = g++
CXXFLAGS = -g -Wall 
SHELL:=/bin/bash


all: heat_mpi heat_omp heat_serial  
	module load openmpi

heat_mpi : heat_mpi.cc 
	mpic++ heat_mpi.cc $(CXXFLAGS) -o heat_mpi

heat_omp : heat_omp.cc 
	$(CXX) heat_omp.cc $(CXXFLAGS) -fopenmp -o heat_omp

heat_serial : heat_serial.cc 
	$(CXX) heat_serial.cc $(CXXFLAGS) -o heat_serial

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
