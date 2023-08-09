.SUFFIXES: .cpp .o

ARCH = $(shell uname -m)
#ARCH = aarch64
PLAT = $(shell uname -s | tr '[:upper:]' '[:lower:]')-gnu
# PLAT = w64-mingw32
EXTRA = 
# EXTRA = -lpsapi -liphlpapi

BIN = fes$(VER)

CC = g++ -w

INCDIR = -I./dep/include
LIBDIR = -L./dep/lib/$(ARCH)-$(PLAT)/

BINDIR  = ./bin/$(ARCH)-$(PLAT)
OBJDIR  = ./obj/$(ARCH)-$(PLAT)
SRCDIR  = ./src$(VER)

CFLAGS = $(INCDIR) -std=c++17 -fopenmp -O2 -DTETLIBRARY -DTRILIBRARY
LFLAGS = $(LIBDIR) -std=c++17 -fopenmp -static \
	-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord \
	-ltet -ltriangle \
	-larpack -lopenblas -lgfortran -lquadmath \
	-Wl,--whole-archive -lpthread -Wl,--no-whole-archive
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
	$(BINDIR)/$(BIN) ../models/WR10.hfss