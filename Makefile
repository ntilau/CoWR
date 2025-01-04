.SUFFIXES: .cpp .o

ARCH = $(shell uname -m)
ARCH = x86_64
PLAT = $(shell uname -s | tr '[:upper:]' '[:lower:]')-gnu
PLAT = w64-mingw32
#EXTRA = -bundle -undefined dynamic_lookup
EXTRA = -lpsapi -liphlpapi

BIN = core

CC = $(ARCH)-$(PLAT)-g++
INCDIR = -I./dep/include
LIBDIR = -L./bin/$(ARCH)-$(PLAT)/ 
#-L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/current/
BINDIR  = ./bin/$(ARCH)-$(PLAT)
OBJDIR  = ./obj/$(ARCH)-$(PLAT)
SRCDIR  = ./src

CFLAGS = $(INCDIR) -std=c++17 -O2 -DTETLIBRARY -DTRILIBRARY
LFLAGS = $(LIBDIR) -std=c++17 \
	-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord \
	-ltet -ltriangle \
	-larpack -lopenblas -lgfortran -lquadmath -lpthread \
	$(EXTRA)
LFLAGS = $(LIBDIR) -std=c++17 \
	-lmumps \
	-larpack -lopenblas -lgfortran -lquadmath -lpthread \
	$(EXTRA)

SRCS=$(wildcard  $(SRCDIR)/*.cpp)
OBJS=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))

all: $(BINDIR) $(OBJDIR) $(BIN)

$(BIN): $(OBJS)
	$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS)

$(BINDIR):
	if [ ! -d "$(BINDIR)" ]; then mkdir -p $(BINDIR); fi

$(OBJDIR):
	if [ ! -d "$(OBJDIR)" ]; then mkdir -p $(OBJDIR); fi

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c  $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/$(BIN)

.PHONY: test
test:
	cd $(BINDIR)
	$(BINDIR)/$(BIN) $(BINDIR)/../data/Strip.core
