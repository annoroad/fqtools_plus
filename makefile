tmpdir=$(dir  $(abspath $(firstword $(MAKEFILE_LIST))))
BIN=$(tmpdir)
gcc=/usr/bin/gcc
ar=/usr/bin/ar
zip=$(BIN)/speedy_gzip
install:
	cd $(BIN)
	$(gcc) -Wall -c yarn.c -O3
	$(ar) cr libyarn.a yarn.o
	$(gcc) -Wall -c hash.c  -O3
	$(ar) cr libhash.a hash.o
	$(gcc) -Wall -c cutadapt.c -lz -O3
	$(ar) cr libcutadapt.a cutadapt.o
	$(gcc) -Wall -O3 -o fqtools_plus fqtools.c -lz -lm -lpthread -L. -lyarn -lhash -lcutadapt -Dspeedy_gzip_software="\"$(zip)\""
	@echo "done"
