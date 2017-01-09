SHELL := /bin/bash
MKDIR=../Make
include $(MKDIR)/MkFlags
PROF= #-pg
MOD=../Modules/
EXE = a bs_corr logZ prob
SRC= a.cc bs_corr.cc logZ.cc prob.cc

OBJ = $(SRC:.cc=.o)
DEPS=$(SRC:.cc=.d)
INC=-I ../Include -I ../Modules/LibAlg/src
MODOBJ=$(wildcard $(MOD)/*.o)

all: $(MOD) $(EXE)
.PHONY: mod $(MOD)

mod:
	$(MAKE) -C $(MOD)

$(OBJ): $(SRC)

$(DEPS): $(SRC)


$(EXE): $(OBJ)

.PHONY cleandir:
	$(MAKE) clean -C $(MOD)

clean:cleandir
	rm -f $(EXE) *.o *~ *.d

%.o : %.cc
	$(CC) $(CFLAGS) $(INC) -o $@ $<

%.d : %.cc
	$(CC) $(MM) $(INC) -o $@ $<

% : %.o mod
	$(CC) $(MODOBJ) $< $(INC) -lm -o $@ $(PROF) 
