FC      = gfortran
FFLAGS  = -O3 -fbounds-check -Wall -Wextra -Wunused-parameter -Wunused-variable
LDFLAGS = 

SRC  = $(wildcard *.f90)
OBJ  = $(patsubst %.f90,%.o,$(SRC))

TARGET = afd_1d

$(TARGET): $(OBJ) 
	$(FC) -o $(TARGET) $(OBJ) $(LDFLAGS)

variables.mod: variables.f90
	$(FC) -c $(FFLAGS) variables.f90

%.o: %.f90 variables.mod
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.mod *.o $(OBJ)
