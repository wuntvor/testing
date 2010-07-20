#-------------------------------------------------------------------------
#									
#               M A K E F I L E   F O R   L B C		        
#									
#-------------------------------------------------------------------------

# compile on  Asama/A1/ONTAKE for SX 
# --------------------------------
#PLATFORM?=$(shell uname -m)
#HOSTNAME?=$(shell uname -n|cut --delimiter=. --fields=1)
#

#---------------------
# Define Machine Name here:

#HPC_TYPE=sx8
#HPC_TYPE=sx9
#HPC_TYPE=bwgrid
#HPC_TYPE=cray
HPC_TYPE=nehalem
#HPC_TYPE=bluegene

#---------------------
# Compiler:

# Add _debug for debugging
#COMPILER=gnu
#COMPILER=sx
#COMPILER=vampir
#COMPILER=intel_debug
COMPILER=intel
#COMPILER=ibmxl
#COMPILER=pgi_debug
#COMPILER=cray
#COMPILER=g95




# IBM Compiler: Pass Arguments with additional switch -wf,-D...
ifeq ($(findstring ibmxl,$(COMPILER)),ibmxl) 
 FDEF = -WF,
endif

#-------------------------------------------------------------------------
# Possible Arguments for Parallel

# (arguments are set using -DARG1 -DARG2

#-----------------------------------------------------
# arg                             # description

# For parallel either use USE_MPI with suboptions or use USE_CAF
PARALLEL += $(FDEF)-DUSE_MPI    # Activate MPI:              
#PARALLEL += $(FDEF)-DDEBUG_MPI 
#PARALLEL += $(FDEF)-DVERBOSE   # Output detailed information about parallel execution 
# MPI Method
#PARALLEL += $(DEF)-DPERSISTENT   # Persistent MPI Communication. Must be used together with BUFFER_EACH 
#PARALLEL += $(DEF)-DPERSISTENT2   # Persistent MPI Communication. Must be used together with BUFFER_EACH 
PARALLEL += $(FDEF)-DISEND_IRECV       # use MPI_IRECV and MPI_ISEND
#PARALLEL += $(FDEF)-DSENDRECV          # use MPI_SENDRECV

# MPI options
#PARALLEL += $(DEF)-DMPI_SUBARRAY      # use a mpi derived type subarray for the exchange
PARALLEL += $(DEF)-DBUFFER_EACH       # Use separate buffers for each direction 
PARALLEL += $(FDEF)-DCOMM_REDUCED      # Only Exchange necessary densities
#PARALLEL += $(FDEF)-DSYNCHRONIZE_COMM  # mpl_barrier before communication
#                              
# Use ADCL. USE_MPI has to be activated 
#PARALLEL += $(FDEF)-DUSE_ADCL          

# CoArray Fortran Usage (use instead of mpi)
#PARALLEL += $(FDEF)-DUSE_CAF           
# use  
#PARALLEL += $(FDEF)-DCAF_DER_TYPE
#PARALLEL += $(FDEF)-DCAF_BENCH=2       # instead of communication routine, run a little caf benchmark 


#-----------------------------------------------------
# Model arguments 
# Lattice Boltzmann Model: 
# Choose ONE of the following
# 2D:  use two-dimensional model with nine velocities
LBMODEL += $(FDEF)-DD2Q9
# 3D:  use three-dimensional model with nineteen velocities
#LBMODEL += $(FDEF)-DD3Q19

# Relaxation modes (optional)
# Standard is BGK

# MRT: multi-relaxation time model (not working yet)
#LBMODEL += $(FDEF)-DMRT
#LBMODEL += $(FDEF)-DTEST_MRT
# Check mass and momentum conservation in MRT routine
#LBMODEL += $(FDEF)-DMRT_MASS_MOMENTUM
#LBMODEL += $(FDEF)-DMRT_BGK_PARAMS
# Debug option with sample data for the MRT routine
#LBMODEL += $(FDEF)-DCHECK_MRT

# Use the incompressible model for higher stability
#LBMODEL += $(FDEF)-DINCOMPRESSIBLE
#LBMODEL += $(FDEF)-DSPONGE
# Use the turbulence model LES Smagorinsky for higher Reynolds numbers
#LBMODEL += $(FDEF)-DLES_SMAGORINSKY
#LBMODEL += $(FDEF)-DDEBUG_LES
# TRT: two-relaxation time mode
#LBMODEL += $(FDEF)-DTRT

#Special options for co-rotating vortex pair (averaging)
#LBMODEL += $(FDEF)-DCRVP

# Performance issues (optional)
# Combine stream and collide routine, split relaxation 
#LBMODEL += $(FDEF)-DCOMBINE_STREAM_COLLIDE
#LBMODEL += $(FDEF)-DSPLITLOOPS

# Memory layout:
# Choose ONE of the following
# save array in memory as x-,y-,z- component and l as density (1...9 / 1...19)
LBMODEL += $(FDEF)-DLAYOUT_IJKL
#  save array in memory as l as density (1...9 / 1...19) and x-,y-,z- component 
#LBMODEL += $(FDEF)-DLAYOUT_LIJK

# Initialization of geometry and boundaries (for parallel execution)
#  First process assigns geometry and distributes to other processes
#  Gathers data to output from other processes
LBMODEL += $(FDEF)-DINIT_WITH_ROOT  
#
# Timestep information
#LBMODEL += $(FDEF)-DWRITE_TSTATUS=2 #    Write info about the timestep duration every ith step
# Use Fortran 2003 statements
#LBMODEL += $(FDEF)-DF2003_ALLOC_DERIVED_TYPE # use allocatables in derived type
#LBMODEL += $(FDEF)-DCHECK_CAF
LBMODEL += $(FDEF)-DTEC_OUT # Write Tecplot Binary output 
#LBMODEL += $(FDEF)-DVERBOSE # write out extra information 


CPPDEFINES=$(LBMODEL) $(PARALLEL)

#-----------------------------------------------------
# Use Paramesh AMR library?
USE_PARAMESH = 0
ifeq ($(USE_PARAMESH),1)
   LBMODEL += -DAMR
endif


#----------------------------------------------------------------
# Machine specific 
#----------------------------------------------------------------

#--------------------------------
# Vector Machine

ifeq ($(findstring sx,$(HPC_TYPE)),sx) #vector machine
   COMP=sxmpif90
   CCOMP=sxmpicc

   F90_CFLAGS_DEBUG=-Chopt  -g -EP  -Wf,'-M zdiv,flovf,flunf,fxovf,inv' -Wf,'-init stack=nan,heap=nan'
   F90_LFLAGS_DEBUG=-g -Cdebug

   #F90CFLAGS_OPT=-float0 -f4 -Pmulti -ew -EP -p -f2003 #-Wf"-ptr byte-pvctl noassume -pvctl chgpwr -pvctl loopcnt=200000" #-Wl"-Z 60G"
   F90_CFLAGS_OPT=-Chopt -EP -ftrace #-pvctl fullmsg
   F90_LFLAGS_OPT=-Chopt -ftrace
   ifeq ($(findstring F2003,$(HPC_TYPE)),F2003) #vector machine
      F90_CFLAGS_OPT+= -f2003
      F90_LFLAGS_OPT+= -f2003
      F90_CFLAGS_DEBUG+= -f2003
      F90_LFLAGS_DEBUG+= -f2003
   endif
#   ifeq ($(HPC_TYPE),sx9)
#   #-ftrace -mpitrace
#   else ($(HPC_TYPE),sx8)
#   #
#   endif

else 

#--------------------------------
# no vector machine

#--------------------------------
# IBM Compiler BlueGene/P
   ifeq ($(findstring ibmxl,$(COMPILER)),ibmxl)
      COMP=mpixlf90
      CCOMP=mpixlc

#      CPPDEFINES=-WF,"$(LBMODEL) $(PARALLEL)"

      F90_CFLAGS_OPT=-O3 -qstrict -qarch=450 -qtune=450
      F90_LFLAGS_OPT=-O3 -qstrict -qarch=450 -qtune=450
      ifeq ($(findstring vamp,$(COMPILER)),vamp)
      COMP=bgvtf90 -vt:f90 mpixlf90
      CCOMP=bgvtcc -vt:cc mpixlc
      endif
      ifeq ($(findstring scalasca,$(COMPILER)),scalasca)
      COMP=scalasca -instrument mpixlf90
      CCOMP=scalasca -instrument mpixlc
      endif
        #--------------------------------
      # PROFILING 
#      ifeq ($(findstring profile,$(COMPILER)),profile)
#         F90_CFLAGS_DEBUG+=-p -g
#         F90_LFLAGS_DEBUG+=-p
#         F90_CFLAGS_OPT+=-p -g
#         F90_LFLAGS_OPT+=-p
#      endif
#      ifeq ($(findstring vamp,$(COMPILER)),vamp)
#      COMP=vtf90 -vt:f90 mpif90
#      CCOMP=vtcc -vt:cc mpicc
#      endif
   else

#--------------------------------
# INTEL 
   ifeq ($(findstring intel,$(COMPILER)),intel)
cflags_debug= -i4 -r8 -Mpreprocess -C -g -Mbounds,chkptr,dclchk -O0  
      COMP=mpif90
      CCOMP=mpicc

      F90_CFLAGS_DEBUG=-fpp -O0 -g -check all -r8 -i4 -traceback -implicitnone -warn all -keep
      F90_LFLAGS_DEBUG=-g -O0 -r8 -i4  
      #F90_CFLAGS_DEBUG+=-warn alignments,declarations,errors,general,ignore_loc,interfaces,truncated_source,uncalled,nounused,usage

      F90_CFLAGS_OPT=-fpp -r8 -i4 -O3 -i-dynamic 
      F90_LFLAGS_OPT=-r8  -i4 -O3 -i-dynamic 
      #--------------------------------
      # PROFILING 
      ifeq ($(findstring profile,$(COMPILER)),profile)
         F90_CFLAGS_DEBUG+=-p -g
         F90_LFLAGS_DEBUG+=-p
         F90_CFLAGS_OPT+=-p -g
         F90_LFLAGS_OPT+=-p
      endif
      ifeq ($(findstring vamp,$(COMPILER)),vamp)
      COMP=vtf90 -vt:f90 mpif90
      CCOMP=vtcc -vt:cc mpicc
      endif
   endif

#--------------------------------
# G95
   ifeq ($(findstring g95,$(COMPILER)),g95)
      COMP=g95   
      CCOMP=mpicc

      F90_CFLAGS_DEBUG=-xf95-cpp-input -O0 -g -fbounds-check -Wall
      F90_LFLAGS_DEBUG=-O0 -g -fbounds-check -Wall -pedantic

      F90_CFLAGS_OPT=-std=f2003 -O3
      F90_LFLAGS_OPT=-std=f2003 -O3
      #--------------------------------
      # PROFILING 
      ifeq ($(findstring profile,$(COMPILER)),profile)
         F90_CFLAGS_DEBUG+=-pg
         F90_LFLAGS_DEBUG+=-pg
         F90_CFLAGS_OPT+=-pg
         F90_LFLAGS_OPT+=-pg
      endif
   endif

#--------------------------------
# Marmot

   ifeq ($(findstring marmot,$(COMPILER)),marmot)
      COMP=marmotf90
      CCOMP=marmotcc
   endif

#--------------------------------
# Vampir

   ifeq ($(findstring vampir,$(COMPILER)),vampir)
      COMP=vtf90 -vt:f90 mpif90
      CCOMP=vtcc -vt:cc mpicc

      F90_CFLAGS_DEBUG=-xf95-cpp-input -O0 -g -fbounds-check -Wall -pedantic -fcray-pointer 
      F90_LFLAGS_DEBUG=-O0 -g -fbounds-check -Wall -pedantic 

#      F90_CFLAGS_OPT=-xf95-cpp-input -std=f2003 -O3 -fcray-pointer 
#      F90_LFLAGS_OPT=-std=f2003 -O3

      F90_CFLAGS_OPT=-warn all -fpp -r8 -i4 -O3 -i-dynamic 
      F90_LFLAGS_OPT=-warn all -r8 -i4 -O3 -i-dynamic 
 
      ifeq ($(findstring profile,$(COMPILER)),profile)
         F90_CFLAGS_DEBUG+=-pg
         F90_LFLAGS_DEBUG+=-pg
         F90_CFLAGS_OPT+=-pg
         F90_LFLAGS_OPT+=-pg
      endif
   endif

#--------------------------------
# GNU

   ifeq ($(findstring gnu,$(COMPILER)),gnu)
      COMP=mpif90
      CCOMP=mpicc

      F90_CFLAGS_DEBUG=-xf95-cpp-input -O0 -g -fbounds-check -Wall -pedantic -fcray-pointer 
      F90_LFLAGS_DEBUG=-O0 -g -fbounds-check -Wall -pedantic -fcray-pointer 

      F90_CFLAGS_OPT=-xf95-cpp-input -O3 -fcray-pointer 
#-std=f2003  
      F90_LFLAGS_OPT= -O3 -fcray-pointer 
      #--------------------------------
      # PROFILING 
      ifeq ($(findstring profile,$(COMPILER)),profile)
         F90_CFLAGS_DEBUG+=-pg
         F90_LFLAGS_DEBUG+=-pg
         F90_CFLAGS_OPT+=-pg
         F90_LFLAGS_OPT+=-pg
      endif
   endif
   
#--------------------------------
# CRAY 

   ifeq ($(findstring cray,$(HPC_TYPE)),cray) 
      ifeq ($(findstring cray,$(COMPILER)),cray)
         COMP=ftn 
         CCOMP=cc
   
         # -eZ preprocessing
   
         # -ec Interface checking
         # -eC call site checking
         # -eD all debugging information (same as -O0, -g, -m2, -R bcdspi, and -rl)
         # -eI all variables must be typed
         # -en report nonstandard Fortran usage
         # -g  debugging support (same as -G0)
         # -m 0 messages of type Error, warning, caution, note, and comment
         # -M 300 supress messages of type 300
         # -O0 no optimization
         # -R abcdps runtime checks
   
         F90_CFLAGS_DEBUG= -eZ -ec -eC -eD -eI -en -g  -m 0 -O0 -R abcdps -F 
# -s real64 integer32 
         F90_LFLAGS_DEBUG= -ec -eC -eD -eI -en -g  -m 0 -O0 -R abcdps -F 
   
         # -eo display optimizations options to stderr
         # -G2 full optimization
         # -h msgs write optimization messages to stderr.
         # -h negmsgs write messages explaining why optimizations did not occur in a given instance to stderr
         # -O3 optimization
         # -O ipa5 Aggressive interprocedural analysis (IPA)
         # #-O aggress compiler treats a program unit as a single optimization region => increases compile time
   
         F90_CFLAGS_OPT=-eZ -eo -G2 -h msgs -h negmsgs -O3 -O ipa5 -F
         F90_LFLAGS_OPT=-eo -G2 -h msgs -h negmsgs -O3 -O ipa5 -F
         ifeq ($(findstring CAF,$(PARALLEL)),CAF)
           F90_CFLAGS_OPT   += -hcaf
           F90_LFLAGS_OPT   += -hcaf
           F90_CFLAGS_DEBUG += -hcaf
           F90_LFLAGS_DEBUG += -hcaf
         endif
      endif
   endif

#--------------------------------
# PGI Compiler

   ifeq ($(findstring pgi,$(COMPILER)),pgi)
      COMP=mpif90 
      CCOMP=cc

      # -Mpreprocess preprocessing

      # #-help=debug display information
      # -C Exposes Ansi warnings only.
      # -g Includes debugging information in the object module.
      # -Ktrap=fp,denorm,unf,inexact
      # -M bounds, array bound checking
      #    chkptr,  NULL pointers 
      #       dclchk,  all variables declared? 
      # -O0 no optimization

      F90_CFLAGS_DEBUG= -i4 -r8 -Mpreprocess -C -g -Mbounds,chkptr,dclchk -O0  #crashes with: -Ktrap=fp,denorm,unf,inexact  
      F90_LFLAGS_DEBUG=  -i4 -r8 -C -g -Mbounds,chkptr,dclchk -O0  #crashes with: -Ktrap=fp,denorm,unf,inexact 

      # #-help=opt display information
      # -fast    Generally optimal set of flags for targets that support SSE capability.
     # -fastsse Generally optimal set of flags for targets that include SSE/SSE2 capability.
      # -M info  prints info messages
	 #    ipa  inter-procedural optimization:w
      # -nfast Generally optimal set of flags for the target. Doesn't use SSE.
	 # -O4  Specifies code optimization level where <level> is 0, 1, 2, 3, or 4.

      F90_CFLAGS_OPT= -i4 -r8 -Mpreprocess -fast -fastsse -Minfo,ipa -O4
      F90_LFLAGS_OPT= -i4 -r8 -fast -fastsse -Minfo,ipa -O4
   endif

#--------------------------------
# Scalasca: simply add scalasca -instrument before command

   ifeq ($(findstring scalasca,$(COMPILER)),scalasca)
      COMP=scalasca -instrument $(COMP) 
      CCOMP=scalasca -instrument $(CCOMP) 
   endif


   ifeq ($(findstring cray,$(HPC_TYPE)),cray) 
      COMP=ftn
   endif

#--------------------------------
# LIBRARIES
   ifeq ($(findstring vampir,$(COMPILER)),vampir) 
      PROFEXT+=_vampir
   endif

   endif
endif

   LIBEXT=_$(HPC_TYPE)$(PROFEXT)
   ifeq ($(HPC_TYPE),nehalem)
      LIBEXT=_$(HPC_TYPE)_x86_64$(PROFEXT)
   else 
      ifeq ($(HPC_TYPE),bwgrid)
         LIBEXT=_$(HPC_TYPE)_x86_64$(PROFEXT)
      endif
   endif



ifeq ($(findstring debug,$(COMPILER)),debug)
F90CFLAGS= $(F90_CFLAGS_DEBUG) 
F90LFLAGS= $(F90_LFLAGS_DEBUG) 
else
F90CFLAGS= $(F90_CFLAGS_OPT) 
F90LFLAGS= $(F90_LFLAGS_OPT) 
endif

#--------------------------------
# share compiler and flag variables with sub-Makefiles

export COMP      
export F90CFLAGS 

#----------------------------------------------------------------
#


#--------------------------------
# links and includes 

BASELIB_DIR = ./lib
INCLUDE_DIR = ./include
PARAMESH_INC = amr_headers

LINK_BASELIB=-L$(BASELIB_DIR) 


ifeq ($(findstring ADCL,$(PARALLEL)),ADCL)
ADCL_INCLUDES=$(INCLUDE_DIR)
LINK_BASELIB+=-ladcl$(LIBEXT)
INCLUDES=-I$(BASELIB_DIR) -I$(ADCL_INCLUDES)
endif

LIBS= $(LINK_BASELIB)
ifeq ($(findstring USE_CAF,$(PARALLEL)),USE_CAF)
   LIB_CAF= 
#-L/opt/xt-pe/2.1.50HD/lib_TV/snos64
#   LIBS+=$(LIB_CAF)
else
   LIB_CAF=
endif
ifeq ($(USE_PARAMESH),1)
   INCLUDES += -I./$(PARAMESH_INC)
   LIBS+= -lmpi_paramesh
   LIBS+= -lparamesh
   LIBS+= -lmodules
   export PARAMESH_INC
endif
ifeq ($(findstring TEC_OUT,$(LBMODEL)),TEC_OUT)
   LIBS+= -lstdc++
   TEC_LIB=./lib/tecio.a
else
   TEC_LIB=
endif



# paramesh
# Shell:
# ------
SHELL=/bin/bash

sources := nrt_lib.F90 timing.F90 lbmodel.F90 tools.F90 mpl_set.F90
sources += function_set.F90 lb_geo.F90 lb_init.F90 lb_bc.F90 lbm_functions.F90 
ifeq ($(USE_PARAMESH),1)
   sources += amr_1blk_bcset.F90 amr_initial_soln.F90 amr_test_refinement.F90 advance_soln.F90 amr_timestep.F90  
endif

# Object-Files:
# -------------
OBJS := $(sources:.F90=.o)

main := lbc.F90
mainobject := $(main:.F90=.o) 

CMD := lbc

# Suffixes:
# ---------
#.SUFFIXES:
#.SUFFIXES:.o .c .f .f90 .a

# Defult rules:
# -------------
#.f90.o:
#	$(COMP) $(F90CFLAGS) $(CPPDEFINES) $(INCLUDES) -c $<	
#.c.o:
#	$(CCOMP) $(CFLAGS) $(CPPDEFINES) $(INCLUDES) -c $<	
# sparate preprocessor:
# ---------------------
CPP=cpp 
#-traditional-cpp
%.o: %.F90
	$(COMP) $(CPPDEFINES) $(F90CFLAGS) -c $< $(INCLUDES) $(LIB_CAF)  

# Targets:
# --------
#echo:
#	@echo $(OBJS)
ifeq ($(USE_PARAMESH),1)
all: amr_headers amr_source amr_mpi_source $(CMD) 
else
all: $(CMD) 
endif


$(CMD): $(OBJS) $(mainobject)
	$(COMP) $(F90LFLAGS) -o  $@  $^ $(TEC_LIB) $(INCLUDES) $(LIBS)
	@echo -e '\033[32m                                                       '
	@echo '      '
	@echo '         oooo   .o8                       '
	@echo '          888   888                       '
	@echo '          888   888oooo.   .ooooo.        '
	@echo '          888   d88   88b d88    Y8       '
	@echo '          888   888   888 888             '
	@echo '          888   888   888 888   .o8       '
	@echo '         o888o  `Y8bod8P  `Y8bod8P        '
	@echo '      '
	@echo -e '\033[0m                                                        '
	@echo '==============================================================='
	@echo '=  L B C   C O M P I L E D                                    ='
	@echo '==============================================================='

.PHONY: $(LIBS)

.PHONY: amr_headers
amr_headers:
	$(MAKE) -C ./$@ -f Makefile.gnu
	cp -f ./$@/libmodules.a $(BASELIB_DIR)
	@echo '==============================================================='
	@echo '=  A M R   headers compiled                                   ='
	@echo '==============================================================='

amr_mpi_source: amr_headers
	$(MAKE) -C ./$@ -f Makefile.gnu
	cp -f ./$@/libmpi_paramesh.a $(BASELIB_DIR)
	@echo '==============================================================='
	@echo '=  A M R   mpi_source compiled                                ='
	@echo '==============================================================='

amr_source: amr_headers
	$(MAKE) -C ./$@ -f Makefile.gnu
	cp -f ./$@/libparamesh.a $(BASELIB_DIR)
	@echo '==============================================================='
	@echo '=  A M R   source compiled                                    ='
	@echo '==============================================================='

clean:
	- rm -f *.o *.mod *.i lbc
	@echo '==============================================================='
	@echo '=  c l e a n e d                                              ='
	@echo '==============================================================='

cleanall:
	- rm -f *.o *.mod *.i *.i90 lbc
	for dir in amr_headers amr_mpi_source amr_source; do \
	$(MAKE) -C ./$$dir -f Makefile.gnu clean; \
	done
	@echo '==============================================================='
	@echo '=  c l e a n e d     u p    a l l                             ='
	@echo '==============================================================='

cleandata:
	- rm -f *.plt *.res *.inf fort.* 
	@echo '==============================================================='
	@echo '=  c l e a n e d     d a t a                                  ='
	@echo '==============================================================='

touch:
	- touch *.F90
#display:
#	@echo
#"==============================================================================="
#	@echo "Plattform:        $(PLATFORM)"
#	@echo "HPC_TYPE:         $(HPC_TYPE)"
#	@echo "Fortran-Compiler: $(COMP)"
#	@echo "Fortran-Optionen: $(F90CFLAGS)"
#	@echo "C-Compiler:       $(CCOMP)"
#	@echo "Lib-Extension:    $(LIBEXT)"
#	@echo "BASELIB_DIR:      $(BASELIB_DIR)"
#	@echo
#"================================================================================"
	


# Abhaengigkeiten
(lbmodel.o): (nrt_lib.o)
(tools.o): (nrt_lib.o) (lbmodel.o) 
(mpl_set.o): (nrt_lib.o) (lbmodel.o) (tools.o) 
(function_set.o): (nrt_lib.o) (timing.o) (lbmodel.o) (mpl_set.o) 
(init.o): (nrt_lib.o) (timing.o) (lbmodel.o) (mpl_set.o) (lb_geo.o)
(lb_bc.o): (nrt_lib.o) (timing.o) (lbmodel.o) (mpl_set.o) (lb_geo.o)
(lbm_functions.o): (nrt_lib.o) (timing.o) (lbmodel.o) (mpl_set.o) (lb_geo.o) (lb_bc.o)
(lbc.o): (nrt_lib.o) (timing.o) (lbmodel.o) (mpl_set.o) (lbm_functions.o) (function_set.o) 
#(amr_initial_soln.o) (amr_1blk_bcset.o)
