all: build run

coor_flag=-I/usr/include -lfftw3 -lfftw3f -lm
f2py_flag= ${coor_flag} -lgomp  -liomp5
gfor_flag= ${coor_flag} -fopenmp -O3 -lgomp  -liomp5
file=./src/corr_fft.f90

build:
	# ifort correlate.f90 -qopenmp -lgomp -O3
	# f2py -c correlate.f90 -m correlate --f90flags="-qopenmp  -lgomp -O3" -lgomp  --fcompiler=intelem -liomp5
	# gfortran correlate.f90 -fopenmp -lgomp -O3
	# f2py -c correlate.f90 -m correlate --f90flags="-fopenmp  -lgomp -O3" -lgomp  -liomp5
	
	gfortran ${file} ${gfor_flag} -o benchmark.out
	f2py -c ${file} -m corrFFT --f90flags="${gfor_flag}"  ${f2py_flag} 
	mv *.so ./lib
	mv *.out ./bin

clean:
	rm -rf *.so *.out