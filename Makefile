.SUFFIXES: .cpp .o

CC=g++

BINDIR= ./bin
SRCDIR= ./src
OBJDIR= ./src

INCDIR = -I./src -I./dep/include
LIBDIR-WIN = -L./dep/lib/x86_64-w64-mingw32/
LIBDIR-LNX = -L./dep/lib/x86_64-linux-gnu/

ifdef OS
   RM = del /F /S /Q
   FixPath = $(subst /,\,$1)
   LIBDIR = $(LIBDIR-WIN)
else
   ifeq ($(shell uname), Linux)
      RM = rm -f
      FixPath = $1
      LIBDIR = $(LIBDIR-LNX)
   endif
endif

CFLAGS = $(INCDIR) -std=gnu++11 -m64 -O2 -fopenmp -static
LFLAGS = $(LIBDIR) -m64 -fopenmp -static -s -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord -lopenblas -larpack -lgfortran -lquadmath
 
OBJS = $(addprefix $(OBJDIR)/, main.o model.o project.o solver.o)

all: $(OBJS)
	$(CC) -o $(BINDIR)/fes $(OBJS) $(LFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c  $< -o $@

.PHONY: clean
clean:
	$(RM) $(call FixPath,$(OBJDIR)/*.o)


