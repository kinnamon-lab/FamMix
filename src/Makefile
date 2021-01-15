all: sing_asc_lmm.so

sing_asc_lmm.so: sing_asc_lmm.cpp
	"$(R_HOME)/bin/Rscript" --vanilla \
		-e "TMB::compile('sing_asc_lmm.cpp', " \
		-e "safebounds = FALSE, safeunload = FALSE)"
	rm -f *.o
clean:
	rm -f *.o *.so
