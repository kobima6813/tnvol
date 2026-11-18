# Optional standalone build of the shared library (pip install . is the usual path).
# `make` builds the portable library everywhere, and also the AVX2 library on x86-64.
NPAIR ?= 16
FPFLAGS ?= -ffp-contract=off
CFLAGS_PORTABLE ?= -O3 $(FPFLAGS) -DNPAIR=$(NPAIR) -fPIC
CFLAGS_AVX2 ?= -O3 -mavx2 -mfma $(FPFLAGS) -DNPAIR=$(NPAIR) -fPIC
CC ?= cc

UNAME_S := $(shell uname -s 2>/dev/null || echo unknown)
UNAME_M := $(shell uname -m 2>/dev/null || echo unknown)

ifeq ($(OS),Windows_NT)
  LIB = src/tnvol/libwlc.dll
  LIBAVX = src/tnvol/libwlc_avx2.dll
  LIBS =
else ifeq ($(UNAME_S),Darwin)
  LIB = src/tnvol/libwlc.dylib
  LIBAVX = src/tnvol/libwlc_avx2.dylib
  LIBS = -lm
else
  LIB = src/tnvol/libwlc.so
  LIBAVX = src/tnvol/libwlc_avx2.so
  LIBS = -lm
endif

X86_64 := $(filter x86_64 amd64 x64,$(UNAME_M))

.PHONY: all avx2 clean
all: $(LIB)
ifneq ($(X86_64),)
all: $(LIBAVX)
endif

$(LIB): src/tnvol/wlc.c
	$(CC) $(CFLAGS_PORTABLE) -shared -o $@ $< $(LIBS)

$(LIBAVX): src/tnvol/wlc.c
	$(CC) $(CFLAGS_AVX2) -shared -o $@ $< $(LIBS)

avx2: $(LIBAVX)

clean:
	rm -f src/tnvol/libwlc.so src/tnvol/libwlc.dylib src/tnvol/libwlc.dll \
	      src/tnvol/libwlc_avx2.so src/tnvol/libwlc_avx2.dylib src/tnvol/libwlc_avx2.dll \
	      src/tnvol/_libwlc*.so src/tnvol/*.exp src/tnvol/*.lib src/tnvol/*.obj
	rm -rf build dist *.egg-info src/*.egg-info
