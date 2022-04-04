INCLUDE_DIR = ./include
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE #-fopenmp #-Wall
VPATH = ./src

CC = gcc

TARGET = SPIN

TARGET_TEST = TEST

DEPENDENCIES = main_functions.c auxiliar_functions.c \
			   auxiliar_functions_gsl.c auxiliar_functions_vmo.c \
			   dynamical_system.c kepler.c

build: $(TARGET)

$(TARGET): main.c $(DEPENDENCIES) -lgsl -lgslcblas -lm  
		   $(CC) $(CFLAGS) -o $(TARGET) $^

test: $(TARGET_TEST)

$(TARGET_TEST): main_test.c $(DEPENDENCIES) -lgsl -lgslcblas -lm  
				$(CC) $(CFLAGS) -o $(TARGET_TEST) $^


.PHONY: clean
clean:
	-rm -f $(TARGET)
	-rm -f $(TARGET_TEST)
