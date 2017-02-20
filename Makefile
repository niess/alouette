CFLAGS := -O2 -std=c99 -pedantic -Wall
CPPFLAGS := -O2 -pedantic -Wall
INCLUDE := -Iinclude
LIBS := -lTauolaCxxInterface -lTauolaFortran -lm -lstdc++

.PHONY: lib examples clean

lib: lib/libtauola-c.so
	@rm -f *.o

examples: bin/example-basic

clean:
	@rm -rf bin lib *.o

lib/lib%.so: src/%.cpp include/%.h
	@mkdir -p lib
	@g++ -o $@ $(CPPFLAGS) -fPIC $(INCLUDE) -shared $< $(LIBS)

bin/example-%: examples/example-%.c lib
	@mkdir -p bin
	@gcc -o $@ $(CFLAGS) $(INCLUDE) $< -Llib -ltauola-c $(LIBS)
