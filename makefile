-include ../petscdir.mk
CFLAGS     =
CPPFLAGS   =
LIBFILES   = 
TARGET     = main
OBJFILES   = main.o 
CLEANFILES = $(TARGET)
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

all: $(TARGET)
$(TARGET) : $(OBJFILES)
	${CLINKER} -o $(TARGET) $(OBJFILES) ${PETSC_KSP_LIB}
