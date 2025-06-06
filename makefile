#
# This is a simple makefile for building spin-squared calculation code.
#
MQCDir       = $(mqcinstall)
MQCMODS      = $(MQCDir)/NVidia/mod
MQCLIB       = $(MQCDir)/NVidia/lib
LIBS         = -llapack -lblas -L$(MQCLIB)
F03Flags     = 
RunF         = pgfortran -i8 -r8 -Mallocatable=03 -mp
#RunF         = pgfortran -i8 -r8
#
#
# The 'all' rule.
#
#hph all: basisCounting.exe gbs.exe
all: gbs.exe pad.exe basisCounting.exe

#
# Generic rules for building module (*.mod) and object (*.o) files.
#
mqc_integrals1.mod: mqc_integrals1.F03 $(MQCLIB)/libmqc.a
	$(RunF) -I$(MQCMODS) -c $*.F03 $(MQCLIB)/libmqc.a

%.o: %.f90
	$(RunF) -I$(MQCMODS) -c $*.f90

%.o: %.f03
	$(RunF) $(F03Flags) -I$(MQCMODS) -c $*.f03

#
# Generic rule for building general executable program (*.exe) from a standard
# f03 source (*.f03) file.
#

%.exe: %.f03 mqc_integrals1.mod gbs_mod.f03 $(MQCLIB)/libmqc.a
	$(RunF) $(LIBS) $(Prof) -I$(MQCMODS) -o $*.exe $*.f03 mqc_integrals1.o $(MQCLIB)/libmqc.a

#
# Some clean rules.
#
clean:
	(rm -f *.o ; rm -f *.mod ; rm -f *.exe)	
cleanGTests:
	(cd GTests ; rm -f *.chk ; rm -f *.mat ; rm -f *.log)	
