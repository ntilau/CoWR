.SUFFIXES: .cpp .o

ARCH = $(shell uname -m)
#ARCH = aarch64
PLAT = $(shell uname -s | tr '[:upper:]' '[:lower:]')-gnu
# PLAT = w64-mingw32
EXTRA = -bundle -undefined dynamic_lookup
# EXTRA = -lpsapi -liphlpapi

BIN = core

CC = g++
INCDIR = -I./dep/include
LIBDIR = -L./dep/lib/$(ARCH)-$(PLAT)/ -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/current/

BINDIR  = ./bin/$(ARCH)-$(PLAT)
OBJDIR  = ./obj/$(ARCH)-$(PLAT)
SRCDIR  = ./src

CFLAGS = $(INCDIR) -std=c++17 -O2 -DTETLIBRARY -DTRILIBRARY
LFLAGS = $(LIBDIR) -std=c++17 \
	-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord \
	-ltet -ltriangle \
	-larpack -lopenblas -lgfortran -lquadmath -lpthread

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
	$(BINDIR)/$(BIN) $(BINDIR)/../data/Strip.core
