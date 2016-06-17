include Makefile.inc

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
COMMON_OBJ = $(OBJ)/dispersion.o $(OBJ)/h5_cinter.o $(OBJ)/h5_finter.o \
	     $(OBJ)/hetfldsrc.o $(OBJ)/iniparser.o \
	     $(OBJ)/memory.o $(OBJ)/makeuh.o $(OBJ)/modelio.o \
	     $(OBJ)/readini.o $(OBJ)/srcreg.o $(OBJ)/srcsetup.o \
	     $(OBJ)/test1d.o $(OBJ)/tpfree.o $(OBJ)/zero.o
HASKELL_OBJ = $(OBJ)/haskgrn.o
GRAPH_OBJ = $(OBJ)/srgraphs.o
ELWAVE_OBJ = $(INTERFACE_OBJ) $(COMMON_OBJ) \
	     $(FD_OBJ) $(HASKELL_OBJ) $(GRAPH_OBJ) $(OBJ)/elwave.o

all: $(ELWAVE)

$(ELWAVE): $(ELWAVE_OBJ)
	$(MPIF90) $(FFLAGS) $(INCF) -o $(ELWAVE) $(ELWAVE_OBJ) $(LIBALL)

$(OBJ)/asmble.o: asmble.f
	$(FC) $(FFLAGS) $(INCF) -c asmble.f -o $(OBJ)/asmble.o

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

$(OBJ)/hetfldsrc.o: hetfldsrc.f
	$(FC) $(FFLAGS) $(INCF) -c hetfldsrc.f -o $(OBJ)/hetfldsrc.o

$(OBJ)/iniparser.o: iniparser.c
	$(CC) $(CFLAGS) $(INC_INI) -c iniparser.c -o $(OBJ)/iniparser.o

$(OBJ)/interface.o: interface.f90
	$(FC) $(FFLAGS) $(INCF) $(MODULE) -c interface.f90 -o $(OBJ)/interface.o

$(OBJ)/makeuh.o: makeuh.f
	$(FC) $(FFLAGS) $(INCF) -c makeuh.f -o $(OBJ)/makeuh.o

$(OBJ)/memory.o: memory.f90
	$(F90) $(FFLAGS) $(INCF) -c memory.f90 -o $(OBJ)/memory.o

$(OBJ)/modelio.o: modelio.f90
	$(F90) $(FFLAGS) $(INCF) -c modelio.f90 -o $(OBJ)/modelio.o

$(OBJ)/modules.o: modules.F90
	$(MPIF90) $(FFLAGS) $(INCF) $(MODULE) -c modules.F90 -o $(OBJ)/modules.o

$(OBJ)/pfdfd3d.o: pfdfd3d.f
	$(FC) $(FFLAGS) $(INCF) -c pfdfd3d.f -o $(OBJ)/pfdfd3d.o

$(OBJ)/pml3d.o: pml3d.f
	$(FC) $(FFLAGS) $(INCF) -c pml3d.f -o $(OBJ)/pml3d.o

$(OBJ)/srcsetup.o: srcsetup.f
	$(FC) $(FFLAGS) $(INCF) -c srcsetup.f -o $(OBJ)/srcsetup.o

$(OBJ)/setcoefs3d.o: setcoefs3d.f
	$(FC) $(FFLAGS) $(INCF) -c setcoefs3d.f -o $(OBJ)/setcoefs3d.o

$(OBJ)/srcreg.o: srcreg.f
	$(FC) $(FFLAGS) $(INCF) -c srcreg.f -o $(OBJ)/srcreg.o

$(OBJ)/srgraphs.o: srgraphs.f
	$(FC) $(FFLAGS) $(INCF) -c srgraphs.f -o $(OBJ)/srgraphs.o

$(OBJ)/readini.o: readini.f90
	$(F90) $(FFLAGS) $(INCF) -c readini.f90 -o $(OBJ)/readini.o

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

