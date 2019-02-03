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

amplification_factor.mod: amplification_factor.f90
	$(FC) -c $(FFLAGS) amplification_factor.f90

finite_difference.mod: finite_difference.f90
	$(FC) -c $(FFLAGS) finite_difference.f90

build_matrix.mod: build_matrix.f90
	$(FC) -c $(FFLAGS) build_matrix.f90

tools.mod: tools.f90
	$(FC) -c $(FFLAGS) tools.f90

spectral_analysis.mod: spectral_analysis.f90
	$(FC) -c $(FFLAGS) spectral_analysis.f90

%.o: %.f90 variables.mod amplification_factor.mod finite_difference.mod tools.mod build_matrix.mod spectral_analysis.mod 
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.mod *.o $(OBJ)
