.PHONY: all deps tools tools-old tools-new tools-gffcompare install-tools clean distclean

all: deps tools install-tools

deps:
	rm -rf ./deps/libdeflate/build
	cmake -S ./deps/libdeflate -B ./deps/libdeflate/build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF
	cmake --build ./deps/libdeflate/build -j1
	cd ./deps/xz && make clean || true && ./configure --prefix=$(PWD)/install && make -j1 && make install
	cd ./deps/bzip2 && make clean || true && make libbz2.a
	cd ./deps/htslib && make clean || true && make libhts.a

tools: tools-old tools-new tools-gffcompare

tools-old:
	$(MAKE) -C ./stringtie clean || true
	$(MAKE) -C ./stringtie release

tools-new:
	$(MAKE) -C ./stringtie3 clean || true
	$(MAKE) -C ./stringtie3 release

tools-gffcompare:
	$(MAKE) -C ./gffcompare clean || true
	$(MAKE) -C ./gffcompare release

install-tools:
	mkdir -p ../inst
	cp ./stringtie/stringtie ../inst/stringtie
	mkdir -p ../inst/stringtie3
	cp ./stringtie3/stringtie ../inst/stringtie3/stringtie
	cp ./gffcompare/gffcompare ../inst/gffcompare

clean:
	$(MAKE) -C ./stringtie clean || true
	$(MAKE) -C ./stringtie3 clean || true
	$(MAKE) -C ./gffcompare clean || true
	cd ./deps/htslib && make clean || true
	cd ./deps/bzip2 && make clean || true
	cd ./deps/xz && make clean || true
	rm -rf ./deps/xz/install
	rm -rf ./deps/libdeflate/build
	rm -f ../inst/stringtie
	rm -rf ../inst/stringtie3
	rm -f ../inst/gffcompare

distclean: clean
	find ./deps -name '*.o' -delete
	find ./deps -name '*.a' -delete
	find ./deps -name '.DS_Store' -delete
	find ./deps -name '*.dylib' -delete
	find ./deps -name '*.so' -delete
	find ./deps -name '*.lo' -delete
	find ./deps -name '*.la' -delete
	find ./deps -name '.libs' -type d -prune -exec rm -rf {} + || true
