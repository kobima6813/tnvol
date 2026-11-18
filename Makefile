# Optional standalone build of the shared library (pip install . is the usual path).
CFLAGS ?= -O3 -mavx2 -mfma -ffp-contract=off -DNPAIR=16 -fPIC
CC ?= gcc

src/tnvol/libwlc.so: src/tnvol/wlc.c
	$(CC) $(CFLAGS) -shared -o $@ $< -lm

clean:
	rm -f src/tnvol/libwlc.so src/tnvol/_libwlc*.so
	rm -rf build dist *.egg-info src/*.egg-info
.PHONY: clean
