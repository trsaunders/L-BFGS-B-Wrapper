# Linux settings.
MATLAB = /usr/local/MATLAB/R2011b
MEX         = $(MATLAB)/bin/mex
CXX         = g++4.3
FC          = gfortran4.3
CFLAGS      = -O3 -fPIC -pthread 
FFLAGS      = -O3 -fPIC -Wall -fbounds-check -g -Wno-uninitialized

TARGET = lbfgsb_
OBJS   = lbfgsb_mex.o lbfgsb.o linpack.o blas.o timer.o
CFLAGS += -Wall -ansi -DMATLAB_MEXFILE

all: $(TARGET)

%.o: %.cpp
	$(CXX) $(CFLAGS) -I$(MATLAB)/extern/include -o $@ -c $^

%.o: %.f
	$(FC) $(FFLAGS) -o $@ -c $^

$(TARGET): $(OBJS)
	$(MEX) -cxx CXX=$(CXX) CC=$(CXX) FC=$(FCC) LD=$(CXX) -lgfortran -lm \
        -O -output $@ $^
	rm -f *.o

clean:
	rm -f *.o $(TARGET).mex*