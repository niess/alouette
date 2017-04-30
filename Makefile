CFLAGS := -O2 -std=c99 -pedantic -Wall
FFLAGS := -O2 -fno-second-underscore -fno-backslash -fno-automatic -ffixed-line-length-132

RCHLCURRENTS  = new-currents/RChL-currents
OTHERCURRENTS = new-currents/other-currents
TAUOLA_OBJS =  $(addprefix src/tauola/, \
	formf.o tauola.o curr_cleo.o pkorb.o f3pi.o tauola_extras.o \
	$(RCHLCURRENTS)/rcht_3pi/f3pi_rcht.o \
	$(RCHLCURRENTS)/rcht_3pi/funct_3pi.o \
	$(RCHLCURRENTS)/rcht_common/funct_rpt.o \
	$(RCHLCURRENTS)/rcht_common/value_parameter.o \
	$(RCHLCURRENTS)/rcht_common/FA1RCHL.o \
	$(RCHLCURRENTS)/rcht_common/ffwid3pi.o \
	$(RCHLCURRENTS)/rcht_common/initA1TabKKpi.o \
	$(RCHLCURRENTS)/rcht_common/wid_a1_fit.o \
	$(RCHLCURRENTS)/rcht_common/initA1Tab.o \
	$(RCHLCURRENTS)/rcht_common/wid_a1_fitKKpi.o  \
	$(RCHLCURRENTS)/rcht_common/gaus_integr.o \
	$(RCHLCURRENTS)/rcht_common/gfact.o \
	$(OTHERCURRENTS)/frho_pi_belle.o)

.PHONY: lib examples clean

lib: lib/libalouette.so
	@rm -f *.o $(TAUOLA_OBJS)

examples: bin/example-forward bin/example-backward

clean:
	@rm -rf bin lib *.o $(TAUOLA_OBJS)

lib/lib%.so: src/%.c include/%.h $(TAUOLA_OBJS)
	@mkdir -p lib
	@gcc -o $@ $(CFLAGS) -fPIC -Iinclude -shared $< $(TAUOLA_OBJS) -lm -lgfortran

src/tauola/%.o: src/tauola/%.f
	@gfortran -o $@ $(FFLAGS) -fPIC -c $<

bin/example-%: examples/example-%.c lib
	@mkdir -p bin
	@gcc -o $@ $(CFLAGS) -Iinclude $< -Llib -lalouette -lm
