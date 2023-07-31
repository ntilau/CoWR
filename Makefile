.SUFFIXES: .cpp .o

ARCH = $(shell uname -m)
#ARCH = aarch64
PLAT = $(shell uname -s | tr '[:upper:]' '[:lower:]')-gnu
# PLAT = w64-mingw32
EXTRA = 
# EXTRA = -lpsapi -liphlpapi
LEGACY = 
# LEGACY = _legacy

BIN = fes$(LEGACY)

CC = $(ARCH)-$(PLAT)-g++ -w

INCDIR = -I./dep/include -I./dep/include/arma/include -DTETLIBRARY
LIBDIR = -L./dep/lib/$(ARCH)-$(PLAT)/

BINDIR  = ./bin/$(ARCH)-$(PLAT)
OBJDIR  = ./obj/$(ARCH)-$(PLAT)
SRCDIR  = ./src$(LEGACY)
TEST    = ./test

ARGS = $(TEST)/RectangularWG.poly

CFLAGS = $(INCDIR) -std=c++17 -O3 -fopenmp -static
LFLAGS = $(LIBDIR) -std=c++17 -fopenmp -s \
	-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq_seq -lpord \
	-ltet -ltriangle \
	-larpack -llapack -lblas -lgfortran -lquadmath \
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

.PHONY: test
test:
	$(BINDIR)/$(BIN) $(ARGS)

.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/$(BIN)


