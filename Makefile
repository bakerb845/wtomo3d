include Makefile.inc

ifeq "$(wildcard $(OBJ) )" ""
-include $(shell mkdir $(OBJ)) $(wildcard $(OBJ)/*)
endif
ifeq "$(wildcard $(MODDIR) )" ""
-include $(shell mkdir $(MODDIR)) $(wildcard $(MODDIR)/*)
endif

ELWAVE = elwave

INTERFACE_OBJ = $(OBJ)/modules.o $(OBJ)/interface.o
COMMON_OBJ = $(OBJ)/asmble.o $(OBJ)/calcwts.o $(OBJ)/dispersion.o \
	     $(OBJ)/fdfd3d.o $(OBJ)/fixed3d.o \
	     $(OBJ)/hetfldsrc.o $(OBJ)/pfdfd3d.o \
	     $(OBJ)/pml3d.o $(OBJ)/setcoefs3d.o \
	     $(OBJ)/tpfree.o $(OBJ)/zero.o $(OBJ)/zeroas.o
HASKELL_OBJ = $(OBJ)/haskgrn.o
GRAPH_OBJ = $(OBJ)/test1d.o
ELWAVE_OBJ = $(OBJ)/elwave.o $(INTERFACE_OBJ) $(COMMON_OBJ) $(HASKELL_OBJ) $(GRAPH_OBJ)

all: $(ELWAVE)

$(ELWAVE): $(ELWAVE_OBJ)
	$(MPIF90) $(FFLAGS) $(INCF) -o $(ELWAVE) $(ELWAVE_OBJ)

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

$(OBJ)/haskgrn.o: haskgrn.f
	$(FC) $(FFLAGS) $(INCF) -c haskgrn.f -o $(OBJ)/haskgrn.o

$(OBJ)/hetfldsrc.o: hetfldsrc.f
	$(FC) $(FFLAGS) $(INCF) -c hetfldsrc.f -o $(OBJ)/hetfldsrc.o

$(OBJ)/interface.o: interface.f90
	$(FC) $(FFLAGS) $(INCF) $(MODULE) -c interface.f90 -o $(OBJ)/interface.o

$(OBJ)/modules.o: modules.F90
	$(MPIF90) $(FFLAGS) $(INCF) $(MODULE) -c modules.F90 -o $(OBJ)/modules.o

$(OBJ)/pfdfd3d.o: pfdfd3d.f
	$(FC) $(FFLAGS) $(INCF) -c pfdfd3d.f -o $(OBJ)/pfdfd3d.o

$(OBJ)/pml3d.o: pml3d.f
	$(FC) $(FFLAGS) $(INCF) -c pml3d.f -o $(OBJ)/pml3d.o

$(OBJ)/setcoefs3d.o: setcoefs3d.f
	$(FC) $(FFLAGS) $(INCF) -c setcoefs3d.f -o $(OBJ)/setcoefs3d.o

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

