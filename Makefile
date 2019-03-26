# Use Prefix to define the home folder of the source code.
# It can be different from the folder in which you want to compile and run the mesh generator. 
# In the current directory ./ you only need to have the main.cpp and this Makefile

TRILIBDEFS = -DSINGLE -DTRILIBRARY -DANSI_DECLARATORS

INC_GL = /opt/X11/include/
LIB_GL = /opt/X11/lib/

INC_OMP = /usr/local/opt/libomp/include
LIB_OMP = /usr/local/opt/libomp/lib
OMPLIB = -lfftw3f_threads -lfftw3f #-lfftw3

# Define compiler and optimizer's flags
CC = clang++ -Xpreprocessor -fopenmp -O3 -Wno-c++11-extensions #-Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize 

CFLAGS = -I$(INC_GL) -L$(LIB_GL) -I$(INC_OMP) -L$(LIB_OMP)

LIBS = -lGLU -lGL -lglut $(OMPLIB) -lm -lomp

TARGET = incomp

# List of objects
OBJ_SRC = main.o\
	point.o\
	triangle.o

OBJ = $(OBJ_SRC)

all: $(TARGET)

$(TARGET): $(OBJ) 
		$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) -lm $(LIBS) 

main.o: main.cpp ib.h display.h defines.h globals.h
	$(CC) $(CFLAGS) -c $< -o $@

point.o: point.cpp point.h 
	$(CC) $(CFLAGS) -c $< -o $@

triangle.o: triangle.c triangle.h
	$(CC) $(CSWITCHES) $(TRILIBDEFS) -c $< -o $@

clean:
		rm -rf $(TARGET) $(OBJ)
