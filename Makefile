.SUFFIXES: .cpp .o

ARCH = $(shell uname -m)
#ARCH = aarch64
PLAT = $(shell uname -s | tr '[:upper:]' '[:lower:]')-gnu
# PLAT = w64-mingw32
EXTRA = #./bin/$(ARCH)-$(PLAT)/libtet.so  #\
	./bin/$(ARCH)-$(PLAT)/libmumps_common.so \
	./bin/$(ARCH)-$(PLAT)/libpord.so \
	./bin/$(ARCH)-$(PLAT)/libmpiseq.so \
	./bin/$(ARCH)-$(PLAT)/libsmumps.so \
	./bin/$(ARCH)-$(PLAT)/libdmumps.so \
	./bin/$(ARCH)-$(PLAT)/libcmumps.so \
	./bin/$(ARCH)-$(PLAT)/libzmumps.so \
	./bin/$(ARCH)-$(PLAT)/libmpiseq.so

# EXTRA = -lpsapi -liphlpapi

BIN = fes$(VER)

CC = $(ARCH)-$(PLAT)-g++ -w

INCDIR = -I./dep/include #-I./dep/include/mumps_seq
LIBDIR = #-L./dep/lib/$(ARCH)-$(PLAT)/

BINDIR  = ./bin/$(ARCH)-$(PLAT)
OBJDIR  = ./obj/$(ARCH)-$(PLAT)
SRCDIR  = ./src$(VER)

CFLAGS = $(INCDIR) -std=c++17 -O2 -fopenmp -DTETLIBRARY -DTRILIBRARY
LFLAGS = $(LIBDIR) -L$(BINDIR) -std=c++17 -fopenmp -Wl,-Bdynamic \
	-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord \
	-ltet -ltriangle \
	-larpack -lopenblas -lgfortran # -lquadmath -lm -lc -lgomp -lstdc++ -lgcc_s
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


