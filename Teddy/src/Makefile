.PHONY: all deps tools clean

all: deps tools

deps:
	cd ./deps/libdeflate && rm -rf build && mkdir build && cd build && cmake .. && make
	cd ./deps/xz && ./configure --prefix=$(PWD)/deps/xz/install && make && make install
	cd ./deps/bzip2 && make libbz2.a
	cd ./deps/htslib && make libhts.a

tools: tools-old tools-new

tools-old:
	cd ./stringtie && make release \
	LIBS="../deps/htslib/libhts.a ../deps/libdeflate/build/libdeflate.a ../deps/xz/install/lib/liblzma.a ../deps/bzip2/libbz2.a -lz -llzma -lbz2 -pthread" \
	CFLAGS="-I../deps/htslib -I../deps/libdeflate -I../deps/xz/install/include -I../deps/bzip2"

tools-new:
	cd ./stringtie3 && make clean && make release \
	LIBS="../deps/htslib/libhts.a ../deps/libdeflate/build/libdeflate.a ../deps/xz/install/lib/liblzma.a ../deps/bzip2/libbz2.a -lz -llzma -lbz2 -pthread" \
	CFLAGS="-I../deps/htslib -I../deps/libdeflate -I../deps/xz/install/include -I../deps/bzip2"

clean:
	cd deps/libdeflate/build && make clean || true
	cd deps/xz && make clean || true
	cd deps/bzip2 && make clean || true
	cd stringtie && make clean || true