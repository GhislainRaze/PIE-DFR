# --- CHOICE OF COMPILER --- #
F90 = gfortran
F90_CPPFLAG = -cpp $(P_OPT)

# --- COMPILER OPTIONS --- $
# classical options:
#OPT = -O3 -Wall
# debug options:
OPT = -O0 -Wall -ffpe-trap=overflow,underflow -fbounds-check

# --- SUFFIXES --- #
.SUFFIXES: .out .o .f .h .f90

#-> GEOM
GEOM_O=
LIST_GEOM= $(GEOM_O)

#-> INPUTS
INPUTS_O=GeneralInput.o InitComp.o AsciiTreatment.o CreateTemplate.o GeneralInputBnd.o
LIST_INPUTS = $(INPUTS_O)

#-> IOLIB
IOLIB_O=IoGMSH.o
LIST_IOLIB = $(IOLIB_O)


#-> STRUCTURE
STRUCTURE_O=General.o Structure.o 
LIST_STRUCTURE=$(STRUCTURE_O)


FLAGS = $(OPT)



# --- directory GEOM    ---
###TransformCGNS.o:$(GEOM)TransformCGNS.F90 Structure.o
###	@echo "Compiling TransformCGNS.F90"
###	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c $(GEOM)TransformCGNS.F90)
###
###IsoParamTransform.o:$(GEOM)IsoParamTransform.f90 TransformCGNS.o Structure.o
###	@echo "Compiling IsoParamTransform.f90"
###	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c $(GEOM)IsoParamTransform.f90)
###
###FaceConnectivity.o:$(GEOM)FaceConnectivity.F90 Structure.o
###	@echo "Compiling FaceConnectivity.F90"
###	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c $(GEOM)FaceConnectivity.F90)
###
#### --- directory INPUTS  ---
###AsciiTreatment.o:$(INPUTS)AsciiTreatment.F90 Structure.o
###	@echo "Compiling AsciiTreatment.F90"
###	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c $(INPUTS)AsciiTreatment.F90)
###
###CreateTemplate.o:$(INPUTS)CreateTemplate.F90 Structure.o
###	@echo "Compiling CreateTemplate.F90"
###	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c $(INPUTS)CreateTemplate.F90)
###
###GeneralInput.o:$(INPUTS)GeneralInput.F90 IoCGNS.o IoGMSH.o AsciiTreatment.o \
###	Structure.o
###	@echo "Compiling GeneralInput.F90"
###	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c $(INPUTS)GeneralInput.F90)
###
###GeneralInputBnd.o:$(INPUTS)GeneralInputBnd.F90 Structure.o
###	@echo "Compiling GeneralInputBnd.F90"
###	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c $(INPUTS)GeneralInputBnd.F90)
###
###InitComp.o:$(INPUTS)InitComp.f90 TransformCGNS.o Structure.o
###	@echo "Compiling InitComp.f90"
###	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c $(INPUTS)InitComp.f90)

ReadMesh.o:ReadMesh.f90 structure.o
	@echo "Compiling ReadMesh.f90"
	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c ReadMesh.f90)


# --- directory MAIN    ---
Main.o:Main.f90 $(LIST_STRUCTURE) $(LIST_IOLIB) $(LIST_INPUTS) \
	$(LIST_EOS) $(LIST_GEOM) $(LIST_BOCO) $(LIST_FCONV) $(LIST_TIME) \
	Rhs.o $(LIST_OUTPUTS) $(LIST_JOIN)
	@echo "Compiling Main.F90"
	@($(F90) $(F90_CPPFLAG) $(FLAGS) $(P_OPT) -c Main.90)

constants.o: constants.F90
	@echo "Compiling constants.F90"
	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c constants.f90)

structure.o: structure.f90 constants.o
	@echo "Compiling structure.f90"
	@($(F90) $(F90_CPPFLAG) $(FLAGS) -c structure.f90)

MAIN.exe: Main.o
	@echo "Building MAIN.exe"
	@($(F90) $(OPT) -o MAIN.exe Main.o \
	$(LIST_STRUCTURE) \
	$(LIST_IOLIB) \
	$(LIST_INPUTS) \
	$(LIST_JOIN) \
	$(LIST_GEOM) \
	$(LIST_BOCO) \
	$(LIST_FCONV) \
	$(LIST_OUTPUTS) \
	$(LIST_RHS) \
	$(LIST_SOURCE) \
	$(LIST_TIME) \
	$(LIST_EOS) \
	$(CGNSLIB)  )

clean: 
	@(rm -f *.o *.mod *.kmo *.d work.* ifc* MAIN.exe *~ )

