all: sing_asc_lmm.so

sing_asc_lmm.so: sing_asc_lmm.cpp
	Rscript --vanilla -e "TMB::compile('sing_asc_lmm.cpp')"
	rm -f *.o
clean:
	rm -f *.o *.so
