all: libdeflate xz bzip2

libdeflate:
	cd libdeflate && make -j4

#xz:
#	cd xz && ./configure --disable-shared --enable-static && make -j4

bzip2:
	cd bzip2 && make -f Makefile-libbz2_so

clean:
	cd libdeflate && make clean
	cd xz && make clean
	cd bzip2 && make clean