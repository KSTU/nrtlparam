all: nrtl

# Define Fortran compiler
#FC= /opt/intel/composerxe-2011.5.220/bin/ia32/ifort
#FC=/opt/pgi/linux86/10.6/bin/pgfortran  # gfortran
FC=gfortran

nrtl: main.f90
	$(FC) -o nrtl main.f90

#cuda_add.o: cuda_add.cu
#	/usr/local/cuda/bin/nvcc -c cuda_add.cu     #-deviceemu

clean: 
	rm nrtl main.o
