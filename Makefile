.SUFFIXES: .cpp .o

ARCH = $(shell uname -m)
#ARCH = aarch64
PLAT = $(shell uname -s | tr '[:upper:]' '[:lower:]')-gnu
# PLAT = w64-mingw32
EXTRA = 
# EXTRA = -lpsapi -liphlpapi

BIN = fes

CC = -$(ARCH)-$(PLAT)-g++
INCDIR = -I./ext/include
LIBDIR = -L./ext/lib/$(ARCH)-$(PLAT)/

BINDIR  = ./bin/$(ARCH)-$(PLAT)
OBJDIR  = ./obj/$(ARCH)-$(PLAT)
SRCDIR  = ./src

ARGS = $(BIN)/RectangularWG.poly

CFLAGS = $(INCDIR) -std=gnu++17 -O3 -fopenmp -static
LFLAGS = $(LIBDIR) -fopenmp -static -s \
	-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord \
	-ltet -ltriangle \
	-larpack -llapack -lblas -lgfortran -lquadmath \
	$(EXTRA)
SRCS=$(wildcard  $(SRCDIR)/*.cpp)
OBJS=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))
#OBJS = $(addprefix $(OBJDIR)/, main.o model.o project.o solver.o)

all: $(BINDIR) $(OBJDIR) $(BIN)

$(BIN): $(OBJS)
	$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS)

$(BINDIR):
	if [ ! -d "$(BINDIR)" ]; then mkdir $(BINDIR); fi

$(OBJDIR):
	if [ ! -d "./obj" ]; then mkdir ./obj; fi
	if [ ! -d "$(OBJDIR)" ]; then mkdir $(OBJDIR); fi

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c  $< -o $@

.PHONY: test
test:
	$(BINDIR)/$(BIN) $(ARGS)

.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/$(BIN)


