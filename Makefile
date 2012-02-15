# Linux settings.
MEX         = /usr/local/MATLAB/R2011a/bin/mex
MEXSUFFIX   = 
MATLAB_HOME = /usr/local/MATLAB/R2011a
CXX         = g++4.3
FC          = gfortran4.3
CFLAGS      = -O3 -fPIC -pthread 
FFLAGS      = -O3 -fPIC -Wall -fbounds-check -g -Wno-uninitialized

TARGET = lbfgsb_
OBJS   = lbfgsb_mex.o lbfgsb.o linpack.o blas.o timer.o

CFLAGS += -Wall -ansi -DMATLAB_MEXFILE

all: $(TARGET)

%.o: %.cpp
	$(CXX) $(CFLAGS) -I$(MATLAB_HOME)/extern/include -o $@ -c $^

%.o: %.f
	$(FC) $(FFLAGS) -o $@ -c $^

$(TARGET): $(OBJS)
	$(MEX) -cxx CXX=$(CXX) CC=$(CXX) FC=$(FCC) LD=$(CXX) -lgfortran -lm \
        -O -output $@ $^

clean:
	rm -f *.o $(TARGET)
