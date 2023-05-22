INCLUDE_DIR = ./include
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -fopenmp #-O3 -march=native #-Wall
VPATH = ./src

CC = gcc
CCI = icc

TARGET = SPIN

DEPENDENCIES =	dynamical_system.c \
				aux_vmo.c \
				periodic_orbit.c \
				spin_orbit.c

.PHONY: run run_intel compile compile_intel clean
.SILENT: run clean

run: compile
	./$(TARGET)

run_intel: compile_intel
	./$(TARGET)

clean:
	-rm -f $(TARGET)

python_requirements:
	pip install -r python_tools/requirements.txt

#prerequisites
compile: main.c $(DEPENDENCIES) -lgsl -lgslcblas -lm  
		   @$(CC) $(CFLAGS) -o $(TARGET) $^

compile_intel: main.c $(DEPENDENCIES) -lgsl -lgslcblas -lm  
		   @$(CCI) $(CFLAGS) -o $(TARGET) $^