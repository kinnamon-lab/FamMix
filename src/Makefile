all: FamModel.so

FamModel.so: FamModel.cpp sing_asc_lmm.hpp
	"$(R_HOME)/bin/Rscript" --vanilla \
		-e "TMB::compile('FamModel.cpp', " \
		-e "safebounds = FALSE, safeunload = FALSE)"
	rm -f *.o
clean:
	rm -f *.o *.so
