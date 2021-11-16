PYTHON=  python3
PACKAGE= _core.abi3.so
LIB=     libalouette.so

BUILD_DIR= build

CC=     gcc
LD=     $(CC) -shared
CFLAGS= -O3 -g -Wall
SHARED= -fPIC


.PHONY: lib
lib: lib/$(LIB)

.PHONY: package
package: alouette/$(PACKAGE)


# Bundle TAUOLA
TAUOLAPP_VERSION= 1.1.8
TAUOLAPP_TARBALL= TAUOLA.$(TAUOLAPP_VERSION)-LHC.tar.gz

TAUOLA_DIR=      src/TAUOLA
TAUOLA_F_DIR=    $(TAUOLA_DIR)/tauola-fortran/tauola-modified
TAUOLA_I_DIR=    $(TAUOLA_DIR)/src/tauolaFortranInterfaces
TAUOLA_INCLUDES= $(TAUOLA_F_DIR)/new-currents/RChL-currents/funct_declar.inc \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/parameter.inc
TAUOLA_SOURCES=  $(TAUOLA_F_DIR)/curr_cleo.f \
                 $(TAUOLA_F_DIR)/formf.f \
                 $(TAUOLA_F_DIR)/tauola.f \
                 $(TAUOLA_F_DIR)/f3pi.f \
                 $(TAUOLA_F_DIR)/pkorb.f \
                 $(TAUOLA_F_DIR)/new-currents/other-currents/frho_pi_belle.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_3pi/f3pi_rcht.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_3pi/funct_3pi.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_common/FA1RCHL.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_common/ffwid3pi.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_common/funct_rpt.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_common/gaus_integr.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_common/gfact.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_common/initA1Tab.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_common/initA1TabKKpi.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_common/value_parameter.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_common/wid_a1_fit.f \
                 $(TAUOLA_F_DIR)/new-currents/RChL-currents/rcht_common/wid_a1_fitKKpi.f \
                 $(TAUOLA_I_DIR)/tauola_extras.f
TAUOLA_OBJS=     $(BUILD_DIR)/tauola.o
TAUOLA_FFLAGS=   $(CFLAGS) -w \
                 -fno-automatic -fno-backslash -ffixed-line-length-132 \
                 -J$(BUILD_DIR)

$(BUILD_DIR)/tauola.o: src/tauola.f | build
	$(CC) -o $@ $(TAUOLA_FFLAGS) $(SHARED) -c $<

src/tauola.f: src/wrap-tauola.py
	$(PYTHON) $< -s $(TAUOLA_SOURCES) -i $(TAUOLA_INCLUDES) -w $@

$(TAUOLA_DIR):
	cd src && \
	wget https://tauolapp.web.cern.ch/tauolapp/resources/TAUOLA.$(TAUOLAPP_VERSION)/$(TAUOLAPP_TARBALL) && \
	tar -xf $(TAUOLAPP_TARBALL) && \
	rm -f $(TAUOLAPP_TARBALL)


# Build the library
ALOUETTE_INCLUDES= include/alouette.h
ALOUETTE_OBJS=     $(BUILD_DIR)/alouette.o

$(BUILD_DIR)/alouette.o: src/alouette.c $(ALOUETTE_INCLUDES) | build
	$(CC) -o $@ $(CFLAGS) $(SHARED) -Iinclude -c $<


lib/$(LIB): $(ALOUETTE_OBJS) $(TAUOLA_OBJS) | libdir
	$(LD) -o $@ $^ -lm


# Build the package
alouette/$(PACKAGE): lib/$(LIB) $(ALOUETTE_INCLUDES) src/build-alouette.py setup.py
	@cd alouette && ln -fs ../lib/$(LIB)
	$(PYTHON) setup.py build --build-lib .


# Build examples
EXAMPLES_CFLAGS= $(CFLAGS) -Iinclude -Llib -Wl,-rpath $(PWD)/lib
EXAMPLES_LDFLAGS= -lalouette -lm

examples: bin/example-forward bin/example-backward

bin/example-%: examples/%.c lib | bin
	$(CC) -o $@ $(EXAMPLES_CFLAGS) $< $(EXAMPLES_LDFLAGS)


.PHONY: bin
bin:
	@mkdir -p bin

.PHONY: build
build:
	@mkdir -p build

.PHONY: libdir
libdir:
	@mkdir -p lib

.PHONY: clean
clean:
	rm -rf bin build lib alouette/libalouette.*

.PHONY: dist-clean
dist-clean:
	rm -rf .eggs alouette.egg-info dist/*.whl
