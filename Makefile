CC      = gcc -fcommon
#CC      = x86_64-linux-gnu-gcc -g
CC2     = g++ -std=c++11 -fcommon
#CC2     = x86_64-linux-gnu-g++ -g

LDIR    = lib
DEPS    = $(LDIR)/mcce.h
LIB     = $(LDIR)/mcce.a
AR      = ar
ARFLAGS = rvs
STEP6   = $(LDIR)/analysis_adv.o

SRC     = $(wildcard $(LDIR)/*.c)
OBJ     = $(SRC:.c=.o)

BIN     = bin
DELPHI  = $(BIN)/delphi
MCCE    = $(BIN)/mcce

# default: Build mcce
default: $(MCCE)

# mcce: Build mcce only
mcce: $(MCCE)

# delphi: Build delphi only
delphi: $(DELPHI)

# ngpb = ngpb container
ngpb:
	@echo "Building Apptainer Image for NGPB container (.sif) ..."
	sudo apptainer build $(BIN)/NextGenPB_MCCE4.sif $(BIN)/recipe_MCCE.def

# all = Build  mcce, delphi, ngpb
all: $(MCCE) $(DELPHI) ngpb

# bin/mcce = core solver
$(MCCE): mcce.c $(LIB) $(DEPS) $(STEP6)
	$(CC2) -o $(MCCE) mcce.c $(STEP6) $(LIB) -lm

$(STEP6): lib/analysis_adv.cpp $(DEPS)
	cd $(LDIR) && $(CC2) -c -o analysis_adv.o ../lib/analysis_adv.cpp $(CFLAGS)

$(LIB): $(OBJ)
	cd $(LDIR) && $(AR) $(ARFLAGS) mcce.a *.o

$(LDIR)/%.o: $(LDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# bin/delphi = copied after external build
$(DELPHI): $(LDIR)/delphi/delphi
	cp $(LDIR)/delphi/delphi $(DELPHI)

$(LDIR)/delphi/delphi:
	$(MAKE) -C $(LDIR)/delphi

# Clean
clean:
	-rm -f $(MCCE) $(DELPHI) $(LIB) $(LDIR)/*.o
	$(MAKE) clean -C $(LDIR)/delphi

cleanbin/mcce:
	-rm -f $(MCCE) $(DELPHI) $(LIB) $(LDIR)/*.o
	$(MAKE) clean -C $(LDIR)/delphi
