
include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${EMPIRICALGMDIR} -I${EMPIRICALIEDIR} -I${DATAREADINDICESDIR}

include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

OBJECTS = \
        ModPlanet.o\
	ModCrcm.o\
	ModCrcmPlot.o\
	ModFieldTrace.o\
	ModGmCRCM.o\
	ModGrid.o\
	ModIeCrcm.o\
	ModImSat.o\
	ModImTime.o\
	ModInitialize.o\
	ModRestart.o\
	ModCrcmBoundary.o\
	ModTsyInput.o\
	ModPrerunField.o\
	crcm.o\
	set_parameters.o\
	plot_fieldline.o\
	trace_dipole.o\
	geopack.o

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS}


MY_LIB = libCRCM.a

LIB: DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} 

# Object files not included into the library
OBJECTS_EXE = crcm_main.o

# Other requireed libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a
LIBINDICES  = ${LIBDIR}/libINDICES.a
LIBEMPIRICALIE   = ${LIBDIR}/libEMPIRICALIE.a
LIBEMPIRICALGM = ${LIBDIR}/libEMPIRICALGM.a

# Libraries should be compiled first, because modules are used in main
${OBJECTS_EXE}: ${LIBSHARE} ${LIBTIMING} ${MY_LIB}   

EXE = ${BINDIR}/crcm.exe

CRCM:
	@echo ' '
	@make ${EXE}
	@echo "Program ${EXE} has been brought up to date."
	@echo ' '

${EXE}: ${OBJECTS_EXE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBSHARE}; \
		ar -x ${LIBINDICES}; \
		ar -x ${LIBEMPIRICALGM};\
		ar -x ${LIBEMPIRICALIE};
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o ${LBLAS} ${Lflag1}
	rm -rf Tmp_

CRCM_orig:
	${COMPILE.f90} ${Lflag1} -o ${EXE} crcm_main.f90 ${OBJECTS} \
	-L${LIBDIR} -lTIMING -lSHARE -lNOMPI -lINDICES

distclean: clean
	rm -f ModPlanet.f90 ModGrid.f90
