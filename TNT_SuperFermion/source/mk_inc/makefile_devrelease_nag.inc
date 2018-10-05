# Include this file to use the NAG MKL libraries and compile with the Intel compilers on the ARCUS cluster.

# Define the architecture:
ARCH := glnxa64

# Root Directory for application
APP_RT := /system/software/linux-x86_64

# Directory for each of the libraries required
DIR_CDF := $(APP_RT)/lib/netcdf/4.4__intel-2013/lib
DIR_MAT := $(APP_RT)/matlab/R2015a
DIR_TNT := $(TNT_ROOT)/dev-release
DIR_NAG := $(APP_RT)/lib/nag/fortran/fll6i25dcl

# Define the linked MATLAB libraries and MATLAB flags
MATLIBS := -L$(DIR_MAT)/bin/$(ARCH) -L$(TMW_ROOT)/sys/os/$(ARCH) -lmat -lmx
RPATH := -Wl,-rpath-link,$(DIR_MAT)/bin/$(ARCH),-rpath-link,$(DIR_MAT)/sys/os/$(ARCH)
MATFLAGS := -I$(DIR_MAT)/extern/include 

# Define the linked netcdf library
CDFLIBS := -L$(DIR_CDF) -lnetcdf

# Define the linked TNT libraries and TNT flags
TNTLIBS := $(DIR_TNT)/lib/arcusnag/libtnt.a $(DIR_TNT)/lib/arcusnag/libtntMps.a
TNTFLAGS := -I$(DIR_TNT)/include  

# Define the NAG flags
NAGFLAGS := -DLONG_CHARLEN -DRETURN_COMPLEX_PARAM 

# Define all linked libraries
ALLLIBS := $(MATLIBS) $(CDFLIBS) -Wl,--start-group $(DIR_NAG)/lib/libnag_mkl.a $(TNTLIBS) -Wl,--end-group
ALLLIBS := $(ALLLIBS) -liomp5 -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -lpthread -nofor-main
ALLLIBS := $(ALLLIBS) $(RPATH) -DTNT_SAFESVD -DTNT_RDYNAMIC -DTNTNAGAVAIL

# Define all library flags
LIBFLAGS := $(MATFLAGS) $(TNTFLAGS) $(NAGFLAGS)

# Specify the compilers to use and any additional compiler flags
CC := icc
FF := ifort

CFLAGS := $(LIBFLAGS)
FFLAGS := $(ALLLIBS)

# Extra flags for compiling OpenMP parallelism into TNT 
CFLAGS += -DTNT_OMP -openmp
FFLAGS += -DTNT_OMP -openmp

# Flag to say arpack is available and to use SVD as backup to SDD
CFLAGS += -DTNT_SAFESVD
FFLAGS += -DTNT_SAFESVD
