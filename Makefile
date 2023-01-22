.SUFFIXES: .cpp .o

BIN = fes

CC = g++ #
INCDIR = 
LIBDIR = 

BINDIR = ./bin
SRCDIR = ./src
OBJDIR = ./obj

ARGS = $(BINDIR)/RectangularWG +freq 1e10

CFLAGS = $(INCDIR) -std=gnu++11 -m64 -O2 -fopenmp -static
LFLAGS = $(LIBDIR) -m64 -fopenmp -static -s \
	-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq_seq -lpord \
	-lopenblas -larpack -lgfortran -lquadmath
OBJS = $(addprefix $(OBJDIR)/, main.o model.o project.o solver.o)

all: $(OBJDIR) $(BIN)

$(BIN): $(OBJS)
	$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS)
	$(BINDIR)/$(BIN) $(ARGS)

$(OBJDIR):
	if [ ! -d "$(OBJDIR)" ]; then mkdir $(OBJDIR); fi

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c  $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/$(BIN)


