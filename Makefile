INCLUDE_DIR = ./include
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE \
		 -fopenmp #-Wall
VPATH = ./src

CC = gcc

TARGET = SPIN

TARGET_TEST = TEST

DEPENDENCIES = main_functions.c auxiliar_functions.c \
			   auxiliar_functions_gsl.c auxiliar_functions_vmo.c \
			   dynamical_system.c kepler.c

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

compile_test: main_test.c $(DEPENDENCIES) -lgsl -lgslcblas -lm  
				@$(CC) $(CFLAGS) -o $(TARGET_TEST) $^
