# platform: linux, macports, mingw_dynamic, mingw_static
# build: debug, release
PLATFORM = linux
BUILD    = release
INSTALL_PATH=/usr/local


TARGET_SOURCES = \
	hdr_create.cc \
	hdr_fix_ca.cc \
	hdr_denoise.cc \
	hdr_squish.cc \
	hdr_adjust.cc \
	hdr_convert.cc

DEPENDENCY_SOURCES = \
	image.cc \
	fix_chromatic_aberration.cc \
	fix_noise.cc \
	squish_luminance.cc \
	solve_poisson.cc

DEPENDENCY_HEADERS = \
	image.hh \
	fix_chromatic_aberration.hh \
	fix_noise.hh \
	squish_luminance.hh \
	solve_poisson.hh \
	shared_array.hh \
	helpers.hh \
	headers.hh


ifeq (${PLATFORM},linux)
CFLAGS = -pthread -fopenmp -Wall -msse2 -mfpmath=sse -I/usr/include/ImageMagick -I/usr/include/OpenEXR
LDFLAGS = -lIlmImf -lMagick++ -ltiff -lpopt -lstdc++ -lm
EXECUTABLES = ${TARGET_SOURCES:%.cc=%}
endif

ifeq ($(SYSTEM), osx)
CFLAGS = -msse2 -mfpmath=sse -I/usr/local/include/ImageMagick-6 -I/usr/local/include/OpenEXR
LDFLAGS = -lIlmImf -lHalf -ltiff -L/usr/local/Cellar/imagemagick/6.9.2-7/lib -lMagick++-6.Q16 -lpopt -lstdc++ -lm
EXECUTABLES = ${TARGET_SOURCES:%.cc=%}
endif

ifeq  ($(PLATFORM),mingw_dynamic)
CFLAGS = -msse2 -mfpmath=sse -I/mingw/include/ImageMagick -I/mingw/include/OpenEXR
LDFLAGS = -enable-auto-import -lIlmImf -lHalf -lMagick++ -ltiff -lpopt -lstdc++ -lm
EXECUTABLES = ${TARGET_SOURCES:%.cc=%.exe}
endif

ifeq  ($(PLATFORM),mingw_static)
CFLAGS = -msse2 -mfpmath=sse -I/mingw/include/ImageMagick -I/mingw/include/OpenEXR
LDFLAGS = -static -static-libgcc -lIlmImf -lIlmThread -lIex -lHalf -lMagick++ -lMagickWand -lMagickCore -ljpeg -ltiff -lgdi32 -lpopt -lintl -liconv -lz -lltdl -lstdc++ -lm
EXECUTABLES = ${TARGET_SOURCES:%.cc=%.exe}
endif


ifeq ($(BUILD),debug)
CFLAGS := -g3 $(CFLAGS)
endif

ifeq ($(BUILD),release)
CFLAGS := -O3 -funroll-loops -fomit-frame-pointer -DNDEBUG $(CFLAGS)
endif


.PHONY : all
all : ${TARGET_SOURCES:%.cc=%}

% : %.cc $(DEPENDENCY_HEADERS) $(DEPENDENCY_SOURCES) Makefile
	@echo "----  building \"$@\" from \"$<\" ($(PLATFORM)/$(BUILD))  ----"
	g++ $(CFLAGS) $< $(DEPENDENCY_SOURCES) -o $@ $(LDFLAGS)
	@echo


.PHONY : strip
strip : all
	strip $(EXECUTABLES)


.PHONY : clean
clean :
	rm -f $(EXECUTABLES)

.PHONY : install
install:
	cp $(EXECUTABLES) $(INSTALL_PATH)/bin/