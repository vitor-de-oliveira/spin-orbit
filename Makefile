INCLUDE_DIR = ./include
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -fopenmp #-O3 -march=native #-Wall
VPATH = ./src

CC = gcc

TARGET = SPIN

DEPENDENCIES =	dynamical_system.c \
				aux_vmo.c \
				periodic_orbit.c \
				spin_orbit.c

.PHONY: run compile clean
.SILENT: run clean

run: compile
	./$(TARGET)

clean:
	-rm -f $(TARGET)

python_requirements:
	pip install -r python_tools/requirements.txt

#prerequisites
compile: main.c $(DEPENDENCIES) -lgsl -lgslcblas -lm  
		   @$(CC) $(CFLAGS) -o $(TARGET) $^