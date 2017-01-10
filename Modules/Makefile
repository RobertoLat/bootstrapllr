MKDIR=../Make
include $(MKDIR)/MkFlags  

SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o)
LIBOBJECTS = ap.o alglibinternal.o alglibmisc.o interpolation.o optimization.o linalg.o solvers.o specialfunctions.o integration.o
INCLUDE=-I ../Include -I LibAlg/src

DEPS:=$(SOURCES:%.cc=%.d)


default: $(LIBOBJECTS) $(OBJECTS)
$(OBJECTS): $(SOURCES)

$(DEPS): $(SOURCES)

-include $(DEPS) 


%.o: %.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $<

%.d : %.cc
	$(CC) $(MM) $(INCLUDE) -o $@ $<

ap.o: LibAlg/src/ap.cpp LibAlg/src/ap.h
	$(CC) $(CFLAGS) -w LibAlg/src/ap.cpp

alglibinternal.o: LibAlg/src/alglibinternal.cpp LibAlg/src/alglibinternal.h
	$(CC) $(CFLAGS) -w LibAlg/src/alglibinternal.cpp

alglibmisc.o: LibAlg/src/alglibmisc.cpp LibAlg/src/alglibmisc.h
	$(CC) $(CFLAGS) -w LibAlg/src/alglibmisc.cpp

interpolation.o: LibAlg/src/interpolation.cpp LibAlg/src/interpolation.h
	$(CC) $(CFLAGS) -w LibAlg/src/interpolation.cpp

optimization.o: LibAlg/src/optimization.cpp LibAlg/src/optimization.h
	$(CC) $(CFLAGS) -w LibAlg/src/optimization.cpp

linalg.o: LibAlg/src/linalg.cpp LibAlg/src/linalg.h
	$(CC) $(CFLAGS) -w LibAlg/src/linalg.cpp

solvers.o: LibAlg/src/solvers.cpp LibAlg/src/solvers.h
	$(CC) $(CFLAGS) -w LibAlg/src/solvers.cpp

specialfunctions.o: LibAlg/src/specialfunctions.cpp LibAlg/src/specialfunctions.h
	$(CC) $(CFLAGS) -w LibAlg/src/specialfunctions.cpp

integration.o: LibAlg/src/integration.cpp LibAlg/src/integration.h
	$(CC) $(CFLAGS) -w LibAlg/src/integration.cpp




#datasample.o: datasample.cc
#	$(CC) -c $(CFLAGS) $(INCLUDE) -o datasample.o datasample.cc

#datasample.d : datasample.cc
#	$(CC) $(MM) $(INCLUDE) -o $@ $<


.PHONY clean:
	rm -f *.o *.d *~
