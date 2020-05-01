EXEC   = ReadSnapshot

OPTIMIZE += $(OPTS)

OBJS   = main.o read_snapshot.o constants.o io.o

CC     = g++

INCL   = read_snapshot.h constants.h io.h

LIBS   = -lm -lgsl -lgslcblas 

# HDF5 Location in Lux
HDF5_INCL = -I/cm/shared/apps/hdf5/1.10.6/include
HDF5_LIBS = -L/cm/shared/apps/hdf5/1.10.6/lib -lhdf5


CFLAGS  = $(OPTIMIZE) $(HDF5_INCL) -I/home/brvillas/code/gsl/include/
LDFLAGS = $(HDF5_LIBS) -L/home/brvillas/code/gsl/lib

$(EXEC): $(OBJS) 
	 $(CC) $(OBJS) $(CFLAGS) $(LDFLAGS) $(LIBS) -o $(EXEC)   

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

