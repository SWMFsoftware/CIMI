default : CRCM

include Makefile.def

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcSAMI3/Makefile.DEPEND \
		srcSAMI3/Makefile.RULES \
		srcInterface/Makefile.DEPEND


install: 
	touch ${INSTALLFILES}
	./Config.pl -EarthHO -GridDefault

#
#       General Housekeeping
#

CRCM:
	@cd ${SHAREDIR};  	make LIB
	@cd ${NOMPIDIR};	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd ${DATAREADINDICESDIR};make LIB
	@cd src;	make LIB
	@cd src;	make CRCM

CIMI_SAMI:
	@cd ${SHAREDIR};  	make LIB
	@cd ${NOMPIDIR};	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd ${DATAREADINDICESDIR};make LIB
	@cd srcSAMI3;	make LIB
	@cd src;	make LIB
	@cd src;	make CIMI_SAMI


SAMI3:
	@cd ${SHAREDIR};  	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd srcSAMI3;	make LIB
	@cd srcSAMI3;	make SAMI3

NOMPI:
	cd util/NOMPI/src; $(MAKE) LIB

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

TESTDIR = run_test

test:
	@echo "test_compile..." > test.diff
	make   test_compile
	@echo "test_rundir..." >> test.diff
	make   test_rundir
	@echo "test_run..."    >> test.diff
	make   test_run
	@echo "test_check..."  >> test.diff
	make   test_check
	ls -l test_*.diff

test_all:
	@echo "test_compile..." > test.diff
	make   test_compile
	@echo "test_rundir..." >> test.diff
	make   test_rundir
	@echo "test_run..."    >> test.diff
	make   test_run
	@echo "test_check..."  >> test.diff
	make   test_check
	@echo "test_rundir_dipole..." >> test.diff
	make   test_rundir_dipole
	@echo "test_run..."    >> test.diff
	make   test_run
	@echo "test_check_dipole..."  >> test.diff
	make   test_check_dipole
	@echo "test_rundir_Prerun..." >> test.diff
	make   test_rundir_Prerun
	@echo "test_run..."    >> test.diff
	make   test_run
	@echo "test_check_Prerun..."  >> test.diff
	make   test_check_Prerun
	ls -l test_*.diff

test_dipole:
	@echo "test_compile..." > test.diff
	make   test_compile
	@echo "test_rundir_dipole..." >> test.diff
	make   test_rundir_dipole
	@echo "test_run..."    >> test.diff
	make   test_run
	@echo "test_check_dipole..."  >> test.diff
	make   test_check_dipole
	ls -l test_*.diff

test_flux:
	@echo "test_compile..." > test.diff
	make   test_compile
	@echo "test_rundir_flux..." >> test.diff
	make   test_rundir_flux
	@echo "test_run..."    >> test.diff
	make   test_run
	@echo "test_check_flux..."  >> test.diff
	make   test_check_flux
	@echo "test_check_eq..."  >> test.diff
	make   test_check_eq
	ls -l test_*.diff

test_drift:
	@echo "test_compile..." > test.diff
	make   test_compile
	@echo "test_rundir_drift..." >> test.diff
	make   test_rundir_drift
	@echo "test_run..."    >> test.diff
	make   test_run
	@echo "test_check_eq..."  >> test.diff
	make   test_check_eq
	@echo "test_check_drift..."  >> test.diff
	make   test_check_drift
	ls -l test_*.diff

test_Prerun:
	@echo "test_compile..." > test.diff
	make   test_compile
	@echo "test_rundir_Prerun..." >> test.diff
	make   test_rundir_Prerun
	@echo "test_run..."    >> test.diff
	make   test_run
	@echo "test_check_Prerun..."  >> test.diff
	make   test_check_Prerun
	ls -l test_*.diff

test_compile:
	./Config.pl -EarthHO -GridDefault -show
	make CRCM

test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/*.dat ${TESTDIR}/
	cp input/testfiles/PARAM.in ${TESTDIR}/PARAM.in

test_rundir_flux:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/*.dat ${TESTDIR}/
	cp input/testfiles/PARAM.in.test.flux ${TESTDIR}/PARAM.in

test_rundir_drift:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/*.dat ${TESTDIR}/
	cp input/testfiles/PARAM.in.test.drift ${TESTDIR}/PARAM.in

test_rundir_dipole:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/*.dat ${TESTDIR}/
	cp input/testfiles/PARAM.in.test.dipole ${TESTDIR}/PARAM.in

test_rundir_Prerun:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/imf.dat.Prerun ${TESTDIR}/imf.dat
	cp input/testfiles/Indices.dat.Prerun ${TESTDIR}/Indices.dat
	cp input/testfiles/Prerun/*.dat ${TESTDIR}/IM/.
	cp input/testfiles/PARAM.in.test.Prerun ${TESTDIR}/PARAM.in

test_run:
	cd ${TESTDIR}; ./crcm.exe > runlog 

test_check:
	make test_check_flux
	make test_check_eq
	make test_check_drift

test_check_flux:
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmFlux_h.fls \
		output/CrcmFlux_h.fls \
		> test_h_fls.diff
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmFlux_o.fls \
		output/CrcmFlux_o.fls \
		> test_o_fls.diff
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmFlux_e.fls \
		output/CrcmFlux_e.fls \
		> test_e_fls.diff

test_check_drift:
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmDrift_h.vp \
		output/CrcmDrift_h.vp \
		> test_h_vp.diff
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmDrift_o.vp \
		output/CrcmDrift_o.vp \
		> test_o_vp.diff
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmDrift_e.vp \
		output/CrcmDrift_e.vp \
		> test_e_vp.diff

test_check_dipole:
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmDrift_h.vp \
		output/CrcmDrift_h.vp.dipole \
		> test_h_vp_dipole.diff
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmDrift_o.vp \
		output/CrcmDrift_o.vp.dipole \
		> test_o_vp_dipole.diff
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmDrift_e.vp \
		output/CrcmDrift_e.vp.dipole \
		> test_e_vp_dipole.diff

test_check_eq:
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CRCMeq.outs \
		output/CRCMeq.outs \
		> test_eq.diff

test_check_Prerun:
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmFlux_h.fls \
		output/CrcmFlux_h.fls.Prerun \
		> test_h_fls_Prerun.diff
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmFlux_o.fls \
		output/CrcmFlux_o.fls.Prerun \
		> test_o_fls_Prerun.diff
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CrcmFlux_e.fls \
		output/CrcmFlux_e.fls.Prerun \
		> test_e_fls_Prerun.diff
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CRCMeq.outs \
		output/CRCMeq.outs.Prerun \
		> test_eq_Prerun.diff


clean:
	@touch ${INSTALLFILES}
	@cd src; make clean
	@cd srcSAMI3; make clean
	@cd srcInterface; make clean
	@(if [ -d util ];  then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean:
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	cd src; make distclean
	cd srcInterface; make distclean
	rm -f config.log *~

#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/IM
	@(cd ${RUNDIR}; \
		if [ ! -e "EIE/README" ]; then \
			ln -s ${EMPIRICALIEDIR}/data EIE;\
		fi;)
	cd ${RUNDIR}/IM; \
		cp ${IMDIR}/input/quiet*fin . ;\
		cp ${IMDIR}/input/WaveData/*dat . ;\
		mkdir plots restartIN restartOUT
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/crcm.exe .   ; \
		touch core ; chmod 444 core;\
	fi);

rundir_cimi_sami:
	mkdir -p ${RUNDIR}/IM
	@(cd ${RUNDIR}; \
		if [ ! -e "EIE/README" ]; then \
			ln -s ${EMPIRICALIEDIR}/data EIE;\
		fi;)
	cd ${RUNDIR}/IM; \
		cp ${IMDIR}/input/quiet*fin . ;\
		cp ${IMDIR}/input/WaveData/*dat . ;\
		mkdir plots restartIN restartOUT plotsSAMI3
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp input/testfiles/*.dat ${RUNDIR}/ ;\
		cp input/testfiles/PARAM.in ${RUNDIR}/;\
		cp srcSAMI3/sami3_mpi-1.98.namelist ${RUNDIR}/ ;\
		cp srcSAMI3/*.inp ${RUNDIR}/ ;\
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/cimi_sami.exe .   ; \
		touch core ; chmod 444 core;\
	fi);


rundir_sami:
	mkdir -p ${RUNDIR}/IM
	cd ${RUNDIR}/IM; \
		mkdir restartIN restartOUT plotsSAMI3
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp srcSAMI3/sami3_mpi-1.98.namelist ${RUNDIR}/ ;\
		cp srcSAMI3/*.inp ${RUNDIR}/ ;\
		ln -s srcSAMI3/sami3.exe ${RUNDIR}/SAMI3.exe ;\
		cd ${RUNDIR} ; \
		touch core ; chmod 444 core;\
	fi);
