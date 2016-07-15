include Makefile.inc
#FC = ifort
#F90 = ifort
#MPIF90 = ifort
#FFLAGS = -check uninit 

ifeq "$(wildcard $(OBJ) )" ""
-include $(shell mkdir $(OBJ)) $(wildcard $(OBJ)/*)
endif
ifeq "$(wildcard $(MODDIR) )" ""
-include $(shell mkdir $(MODDIR)) $(wildcard $(MODDIR)/*)
endif

ELWAVE = elwave

INTERFACE_OBJ = $(OBJ)/modules.o $(OBJ)/interface.o
FD_OBJ = $(OBJ)/asmble.o $(OBJ)/calcwts.o \
	 $(OBJ)/fdfd3d.o $(OBJ)/fixed3d.o \
	 $(OBJ)/pfdfd3d.o $(OBJ)/pml3d.o \
	 $(OBJ)/setcoefs3d.o $(OBJ)/zeroas.o
COMMON_OBJ = $(OBJ)/calc1d.o $(OBJ)/dispersion.o \
	     $(OBJ)/h5_cinter.o $(OBJ)/h5_finter.o \
	     $(OBJ)/hetfldsrc.o $(OBJ)/info.o $(OBJ)/iniparser.o \
	     $(OBJ)/memory.o $(OBJ)/makeuh.o $(OBJ)/modelio.o \
	     $(OBJ)/nulluest.o $(OBJ)/planewvexp.o \
	     $(OBJ)/readini.o $(OBJ)/recsub.o $(OBJ)/srcreg.o $(OBJ)/srcsetup.o \
	     $(OBJ)/srcsub.o $(OBJ)/sort.o $(OBJ)/test1d.o \
	     $(OBJ)/tpfree.o $(OBJ)/zero.o
SPARSE_OBJ = $(OBJ)/sparse.o
HASKELL_OBJ = $(OBJ)/haskgrn.o
SOLVER_OBJ = $(OBJ)/strumpack_finter.o
GRAPH_OBJ = $(OBJ)/srgraphs.o $(OBJ)/metis_finter.o
ELWAVE_OBJ = $(INTERFACE_OBJ) $(COMMON_OBJ) \
	     $(FD_OBJ) $(HASKELL_OBJ) $(SPARSE_OBJ) \
	     $(GRAPH_OBJ) $(SOLVER_OBJ) $(OBJ)/elwave.o

all: $(ELWAVE)

$(ELWAVE): $(ELWAVE_OBJ)
	$(MPIF90) $(FFLAGS) $(INCF) -o $(ELWAVE) $(ELWAVE_OBJ) $(LIBALL)

$(OBJ)/asmble.o: asmble.f90
	$(F90) $(FFLAGS) $(INCF) -c asmble.f90 -o $(OBJ)/asmble.o

$(OBJ)/calc1d.o: calc1d.f
	$(FC) $(FFLAGS) $(INCF) -c calc1d.f -o $(OBJ)/calc1d.o

$(OBJ)/calcwts.o: calcwts.f
	$(FC) $(FFLAGS) $(INCF) -c calcwts.f -o $(OBJ)/calcwts.o

$(OBJ)/dispersion.o: dispersion.f
	$(FC) $(FFLAGS) $(INCF) -c dispersion.f -o $(OBJ)/dispersion.o

$(OBJ)/fdfd3d.o: fdfd3d.f
	$(FC) $(FFLAGS) $(INCF) -c fdfd3d.f -o $(OBJ)/fdfd3d.o

$(OBJ)/fixed3d.o: fixed3d.f
	$(FC) $(FFLAGS) $(INCF) -c fixed3d.f -o $(OBJ)/fixed3d.o

$(OBJ)/elwave.o: elwave.f90
	$(MPIF90) $(FFLAGS) $(INCF) -c elwave.f90 -o $(OBJ)/elwave.o

$(OBJ)/h5_cinter.o: h5_cinter.c
	$(CC) $(CFLAGS) $(INCC) $(INI_HDF5) -c h5_cinter.c -o $(OBJ)/h5_cinter.o

$(OBJ)/h5_finter.o: h5_finter.f90
	$(F90) $(FFLAGS) $(INCF) $(MODULE) -c h5_finter.f90 -o $(OBJ)/h5_finter.o

$(OBJ)/haskgrn.o: haskgrn.f
	$(FC) $(FFLAGS) $(INCF) -c haskgrn.f -o $(OBJ)/haskgrn.o

$(OBJ)/hetfldsrc.o: hetfldsrc.f90
	$(F90) $(FFLAGS) $(INCF) -c hetfldsrc.f90 -o $(OBJ)/hetfldsrc.o

$(OBJ)/info.o: info.f
	$(FC) $(FFLAGS) $(INCF) -c info.f -o $(OBJ)/info.o

$(OBJ)/iniparser.o: iniparser.c
	$(CC) $(CFLAGS) $(INC_INI) -c iniparser.c -o $(OBJ)/iniparser.o

$(OBJ)/interface.o: interface.f90
	$(FC) $(FFLAGS) $(INCF) $(MODULE) -c interface.f90 -o $(OBJ)/interface.o

$(OBJ)/makeuh.o: makeuh.f
	$(FC) $(FFLAGS) $(INCF) -c makeuh.f -o $(OBJ)/makeuh.o

$(OBJ)/memory.o: memory.f90
	$(F90) $(FFLAGS) $(INCF) -c memory.f90 -o $(OBJ)/memory.o

$(OBJ)/metis_finter.o: metis_finter.c
	$(CC) $(CFLAGS) $(INCC) $(INC_METIS) -c metis_finter.c -o $(OBJ)/metis_finter.o

$(OBJ)/modelio.o: modelio.f90
	$(F90) $(FFLAGS) $(INCF) -c modelio.f90 -o $(OBJ)/modelio.o

$(OBJ)/modules.o: modules.F90
	$(MPIF90) $(FFLAGS) $(INCF) $(MODULE) -c modules.F90 -o $(OBJ)/modules.o

$(OBJ)/nulluest.o: nulluest.f
	$(FC) $(FFLAGS) $(INCF) -c nulluest.f -o $(OBJ)/nulluest.o

$(OBJ)/pfdfd3d.o: pfdfd3d.f
	$(FC) $(FFLAGS) $(INCF) -c pfdfd3d.f -o $(OBJ)/pfdfd3d.o

$(OBJ)/planewvexp.o: planewvexp.f90
	$(F90) $(FFLAGS) $(INCF) -c planewvexp.f90 -o $(OBJ)/planewvexp.o

$(OBJ)/pml3d.o: pml3d.f
	$(FC) $(FFLAGS) $(INCF) -c pml3d.f -o $(OBJ)/pml3d.o

$(OBJ)/setcoefs3d.o: setcoefs3d.f
	$(FC) $(FFLAGS) $(INCF) -c setcoefs3d.f -o $(OBJ)/setcoefs3d.o

$(OBJ)/sort.o: sort.c
	$(CC) $(CFLAGS) $(INCC) -c sort.c -o $(OBJ)/sort.o

$(OBJ)/sparse.o: sparse.f
	$(FC) $(FFLAGS) $(INCF) -c sparse.f -o $(OBJ)/sparse.o

$(OBJ)/srcreg.o: srcreg.f
	$(FC) $(FFLAGS) $(INCF) -c srcreg.f -o $(OBJ)/srcreg.o

$(OBJ)/srcsetup.o: srcsetup.f
	$(FC) $(FFLAGS) $(INCF) -c srcsetup.f -o $(OBJ)/srcsetup.o

$(OBJ)/srcsub.o: srcsub.f
	$(FC) $(FFLAGS) $(INCF) -c srcsub.f -o $(OBJ)/srcsub.o

$(OBJ)/srgraphs.o: srgraphs.f
	$(FC) $(FFLAGS) $(INCF) -c srgraphs.f -o $(OBJ)/srgraphs.o

$(OBJ)/strumpack_finter.o: strumpack_finter.c
	$(MPICC) $(CFLAGS) $(INCC) $(INC_STRUMPACK) -c strumpack_finter.c -o $(OBJ)/strumpack_finter.o

$(OBJ)/readini.o: readini.f90
	$(F90) $(FFLAGS) $(INCF) -c readini.f90 -o $(OBJ)/readini.o

$(OBJ)/recsub.o: recsub.f
	$(FC) $(FFLAGS) $(INCF) -c recsub.f -o $(OBJ)/recsub.o

$(OBJ)/test1d.o: test1d.f
	$(FC) $(FFLAGS) $(INCF) -c test1d.f -o $(OBJ)/test1d.o

$(OBJ)/tpfree.o: tpfree.f
	$(FC) $(FFLAGS) $(INCF) -c tpfree.f -o $(OBJ)/tpfree.o

$(OBJ)/zero.o: zero.f
	$(FC) $(FFLAGS) $(INCF) -c zero.f -o $(OBJ)/zero.o

$(OBJ)/zeroas.o: zeroas.f
	$(FC) $(FFLAGS) $(INCF) -c zeroas.f -o $(OBJ)/zeroas.o

clean:
	$(RM) -rf $(OBJ)/*.o $(MODDIR)/*.mod

