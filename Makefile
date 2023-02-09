INCLUDE_DIR = ./include
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE \
		 -fopenmp #-Wall
VPATH = ./src

CC = gcc

TARGET = SPIN

DEPENDENCIES =	dynamical_system.c \
				aux_vmo.c \
				periodic_orbit.c \
				spin_orbit.c

TARGET_TEST = TEST

# on this moment, june 17 2022, there are
# changes that I did in the dynamical system
# lib that I did not update on pendulum.c
# DEPENDENCIES_TEST =	dynamical_system.c \
					aux_vmo.c \
					spin_orbit.c \
					pendulum.c \
					test_spin_orbit.c

DEPENDENCIES_TEST =	dynamical_system.c \
					aux_vmo.c \
					periodic_orbit.c \
					spin_orbit.c \
					test_spin_orbit.c

.PHONY: run compile test compile_test clean
.SILENT: run test clean

run: compile
	./$(TARGET)

test: compile_test
	./$(TARGET_TEST)

clean:
	-rm -f $(TARGET)
	-rm -f $(TARGET_TEST)

#prerequisites
compile: main.c $(DEPENDENCIES) -lgsl -lgslcblas -lm  
		   @$(CC) $(CFLAGS) -o $(TARGET) $^

compile_test: ../dbg/test.c $(DEPENDENCIES_TEST) -lgsl -lgslcblas -lm  
				@$(CC) $(CFLAGS) -o $(TARGET_TEST) $^
