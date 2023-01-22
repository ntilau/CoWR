.SUFFIXES: .cpp .o

CC=x86_64-w64-mingw32-g++-posix

BINDIR= ./bin
SRCDIR= ./src
OBJDIR= ./obj

#INCDIR = -I./src -I./dep/include
INCDIR = -I./src -I./dep/include -I./dep/include/vtk -I./dep/lib/x86_64-w64-mingw32/mswu
#INCDIR = -I./dep/include/vtk -I/usr/lib/x86_64-linux-gnu/wx/include/gtk3-unicode-3.0 -I/usr/include/wx-3.0 -D_FILE_OFFSET_BITS=64 -DWXUSINGDLL -D__WXGTK__ -pthread
LIBDIR-WIN = -L./dep/lib/x86_64-w64-mingw32
LIBDIR-LNX = -L./dep/lib/x86_64-linux-gnu

LIBWIN = -lkernel32 -luser32 -lgdi32 -lwinspool -lcomdlg32 -ladvapi32 -lshell32 -lole32 -loleaut32 -luuid -lcomctl32 -lwsock32 -lodbc32 -lshlwapi -lpsapi -liphlpapi -lversion -luxtheme -loleacc -lopengl32 -mwindows
LIBALG = -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord -lopenblas -larpack -lgfortran -lquadmath
LIBWX = -lpthread -lwxmsw31u -liconv -lwxjpeg -lwxpng -lwxexpat -lwxregexu -lwxscintilla -lwxzlib -lwxtiff
LIBVTKALL = -lvtkChartsCore -lvtkCommonColor -lvtkCommonComputationalGeometry -lvtkCommonCore -lvtkCommonDataModel -lvtkCommonExecutionModel -lvtkCommonMath -lvtkCommonMisc -lvtkCommonSystem -lvtkCommonTransforms -lvtkDICOMParser -lvtkDomainsChemistry -lvtkdoubleconversion -lvtkexodusII -lvtkexpat -lvtkFiltersAMR -lvtkFiltersCore -lvtkFiltersExtraction -lvtkFiltersFlowPaths -lvtkFiltersGeneral -lvtkFiltersGeneric -lvtkFiltersGeometry -lvtkFiltersHybrid -lvtkFiltersHyperTree -lvtkFiltersImaging -lvtkFiltersModeling -lvtkFiltersParallel -lvtkFiltersParallelImaging -lvtkFiltersPoints -lvtkFiltersProgrammable -lvtkFiltersSelection -lvtkFiltersSMP -lvtkFiltersSources -lvtkFiltersStatistics -lvtkFiltersTexture -lvtkFiltersTopology -lvtkFiltersVerdict -lvtkfreetype -lvtkGeovisCore -lvtkgl2ps -lvtkglew -lvtkhdf5 -lvtkhdf5_hl -lvtkImagingColor -lvtkImagingCore -lvtkImagingFourier -lvtkImagingGeneral -lvtkImagingHybrid -lvtkImagingMath -lvtkImagingMorphological -lvtkImagingSources -lvtkImagingStatistics -lvtkImagingStencil -lvtkInfovisCore -lvtkInfovisLayout -lvtkInteractionImage -lvtkInteractionStyle -lvtkInteractionWidgets -lvtkIOAMR -lvtkIOAsynchronous -lvtkIOCityGML -lvtkIOCore -lvtkIOEnSight -lvtkIOExodus -lvtkIOExport -lvtkIOExportGL2PS -lvtkIOExportPDF -lvtkIOGeometry -lvtkIOImage -lvtkIOImport -lvtkIOInfovis -lvtkIOLegacy -lvtkIOLSDyna -lvtkIOMINC -lvtkIOMotionFX -lvtkIOMovie -lvtkIONetCDF -lvtkIOOggTheora -lvtkIOParallel -lvtkIOParallelXML -lvtkIOPLY -lvtkIOSegY -lvtkIOSQL -lvtkIOTecplotTable -lvtkIOVeraOut -lvtkIOVideo -lvtkIOXML -lvtkIOXMLParser -lvtkjpeg -lvtkjsoncpp -lvtklibharu -lvtklibproj -lvtklibxml2 -lvtkloguru -lvtklz4 -lvtklzma -lvtkmetaio -lvtknetcdf -lvtkogg -lvtkParallelCore -lvtkParallelDIY -lvtkpng -lvtkpugixml -lvtkRenderingAnnotation -lvtkRenderingContext2D -lvtkRenderingCore -lvtkRenderingFreeType -lvtkRenderingGL2PSOpenGL2 -lvtkRenderingImage -lvtkRenderingLabel -lvtkRenderingLOD -lvtkRenderingOpenGL2 -lvtkRenderingSceneGraph -lvtkRenderingUI -lvtkRenderingVolume -lvtkRenderingVolumeOpenGL2 -lvtkRenderingVtkJS -lvtksqlite -lvtksys -lvtkTestingRendering -lvtktheora -lvtktiff -lvtkverdict -lvtkViewsContext2D -lvtkViewsCore -lvtkViewsInfovis -lvtkWrappingTools -lvtkzlib
LIBVTK = -lvtkRenderingOpenGL2 -lvtkRenderingVolumeOpenGL2 -lvtkRenderingVolume -lvtkRenderingUI -lvtkRenderingFreeType -lvtkRenderingAnnotation -lvtkInteractionStyle -lvtkRenderingCore -lvtkCommonColor -lvtkFiltersCore -lvtkFiltersSources -lvtkFiltersGeneral -lvtkImagingCore -lvtkCommonExecutionModel -lvtkCommonDataModel -lvtkCommonMath -lvtkCommonTransforms -lvtkCommonSystem -lvtkCommonCore -lvtkCommonMisc -lvtksys -lvtkloguru -lvtkglew -lvtkfreetype -lvtkzlib

WXCFLAGS = $(INCDIR) -std=gnu++11 -m64 -O2 -pthread -static -fopenmp -DNDEBUG -D_FILE_OFFSET_BITS=64 -DWX_PRECOMP -fno-common -fpermissive -Wall -Wundef -Wunused-parameter -Wno-ctor-dtor-privacy -Woverloaded-virtual -Wno-deprecated-declarations
WXLFLAGS = $(LIBDIR) -m64 -fopenmp -static -s $(LIBALG) $(LIBWX) $(LIBVTK) $(LIBWIN)
#WXLFLAGS = $(LIBDIR) -m64 -fopenmp -static -s -L/usr/lib/x86_64-linux-gnu -pthread   -lwx_gtk3u_xrc-3.0 -lwx_gtk3u_html-3.0 -lwx_gtk3u_qa-3.0 -lwx_gtk3u_adv-3.0 -lwx_gtk3u_core-3.0 -lwx_baseu_xml-3.0 -lwx_baseu_net-3.0 -lwx_baseu-3.0


WXOBJS = $(addprefix $(OBJDIR)/, wxfes.o model.o project.o solver.o wxVTKRenderWindowInteractor.o)

LIBDIR = $(LIBDIR-WIN)

ifdef OS
   RM = del /F /S /Q
   FixPath = $(subst /,\,$1)
else
   ifeq ($(shell uname), Linux)
      RM = rm -f
      FixPath = $1
   endif
endif

CFLAGS = $(INCDIR) -std=gnu++11 -m64 -O2 -fopenmp -static
LFLAGS = $(LIBDIR) -m64 -fopenmp -static -s -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord -lopenblas -larpack -lgfortran -lquadmath
 
OBJS = $(addprefix $(OBJDIR)/, main.o model.o project.o solver.o)

fes: $(OBJS)
	$(CC) -o $(BINDIR)/fes $(OBJS) $(WXLFLAGS)

wxfes: $(WXOBJS)
	$(CC) -o $(BINDIR)/wxfes $(WXOBJS) $(WXLFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(WXCFLAGS) -c  $< -o $@

.PHONY: clean
clean:
	$(RM) $(call FixPath,$(OBJDIR)/*.o)


