.PHONY: all tools-old clean distclean

all: tools-old

tools-old:
	cd ./deps/libdeflate && make clean && make
	cd ./deps/bzip2 && make clean && make
	cd ./deps/xz && make clean && ./configure --prefix=$(PWD)/install && make -j1 && make install
	cd ./deps/htslib && make clean && make libhts.a
	$(MAKE) -C ./stringtie clean
	$(MAKE) -C ./stringtie \
		LIBS="../deps/htslib/libhts.a ../deps/libdeflate/build/libdeflate.a ../deps/xz/install/lib/liblzma.a ../deps/bzip2/libbz2.a -lz -llzma -lbz2 -pthread" \
		CFLAGS="-I../deps/htslib -I../deps/libdeflate -I../deps/xz/install/include -I../deps/bzip2"

clean:
	$(MAKE) -C ./stringtie clean || true
	cd ./deps/htslib && make clean || true
	cd ./deps/libdeflate && make clean || true
	cd ./deps/bzip2 && make clean || true
	cd ./deps/xz && make clean || true
	rm -rf ./deps/xz/install

distclean: clean
	find ./deps -name '*.o' -delete
	find ./deps -name '*.a' -delete
	find ./deps -name '.DS_Store' -delete
	find ./deps -name '*.dylib' -delete
	find ./deps -name '*.so' -delete
	find ./deps -name '*.lo' -delete
	find ./deps -name '*.la' -delete
	find ./deps -name '.libs' -type d -prune -exec rm -rf {} +
