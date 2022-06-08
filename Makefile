SRCDIR = src

TARGET1 = libcmils.so
SRCS1 = $(SRCDIR)/gfunc/cmils.f90 $(wildcard $(SRCDIR)/lib/**/*.f)
OBJS1 = $(addsuffix .o, $(basename $(SRCS1)))
TARGET2 = libils.so
SRCS2 = $(SRCDIR)/gfunc/ils.f90 $(wildcard $(SRCDIR)/lib/**/*.f)
OBJS2 = $(addsuffix .o, $(basename $(SRCS2)))
TARGET3 = libics.so
SRCS3 = $(SRCDIR)/gfunc/ics.f90 $(wildcard $(SRCDIR)/lib/**/*.f)
OBJS3 = $(addsuffix .o, $(basename $(SRCS3)))

FC = gfortran
#FC = ifort

FFLAGS =
LDFLAGS =
LIBS =

SOFLAGS = -shared

# for gfortran
ifeq ($(FC),gfortran)
	FFLAGS += -fPIC -fopenmp
	FFLAGS += -O3 -fall-intrinsics
	LIBS += -llapack -lblas
endif

# for ifort
ifeq (${FC},ifort)
	FFLAGS += -fpic -fast
	LIBS += -qmkl -qopenmp
endif

all: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1) : $(OBJS1)
	$(FC) $(SOFLAGS) $(FFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^
	@echo $@ created!

$(TARGET2) : $(OBJS2)
	$(FC) $(SOFLAGS) $(FFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^
	@echo $@ created!

$(TARGET3) : $(OBJS3)
	$(FC) $(SOFLAGS) $(FFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^
	@echo $@ created!

%.o:%.f90
	$(FC) $(FFLAGS) $(LDFLAGS) $(LIBS) -o $@ -c $<

%.o:%.f
	$(FC) $(FFLAGS) $(LDFLAGS) $(LIBS) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) $(TARGET1) $(OBJS1) $(TARGET2) $(OBJS2) $(TARGET3) $(OBJS3) *.mod
test:
	@echo $(SRCS1)
