.SUFFIXES: .cpp .o

CC=g++ #

BINDIR = ./bin
SRCDIR = ./src
OBJDIR =  /tmp
MODELS = ./models

INCDIR = 
LIBDIR = 

CFLAGS = $(INCDIR) -std=gnu++11 -m64 -O2 -fopenmp -static
LFLAGS = $(LIBDIR) -m64 -fopenmp -static -s -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq_seq -lpord -lopenblas -larpack -lgfortran -lquadmath
 
OBJS = $(addprefix $(OBJDIR)/, main.o model.o project.o solver.o)

fes: $(OBJS)
	$(CC) -o $(BINDIR)/fes $(OBJS) $(LFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c  $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o


