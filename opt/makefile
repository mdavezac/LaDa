CC     := gcc
CXX    := g++
F95    := g95
LD     := g++
AR     := ar rc
RANLIB := ranlib
DEBUG = NO

F95FLAGS   := 
CFLAGS     := 
CXXFLAGS   := -malign-double -ffriend-injection 
DEFS       := 


DEBUG_CFLAGS      := -Wall -Wno-format -g -O0 -Wno-unknown-pragmas -fbounds-check
DEBUG_CXXFLAGS    := ${DEBUG_CFLAGS} -D_DEBUG_LADA_ 
DEBUG_LDFLAGS     := -g 
DEBUG_F95FLAGS    := -g 

RELEASE_CXXFLAGS := ${RELEASE_CFLAGS}
RELEASE_LDFLAGS   := 
RELEASE_CFLAGS    := -Wall -Wno-unknown-pragmas -Wno-format -O3

LIBS:=  
INCS := -I atat -I /home/davezac/usr/include 


ifeq (YES, ${DEBUG})
   CFLAGS      := ${CFLAGS} ${DEBUG_CFLAGS} 
   CXXFLAGS    := ${CXXFLAGS} ${DEBUG_CXXFLAGS}
   LDFLAGS     := ${LDFLAGS} ${DEBUG_LDFLAGS}
   F95FLAGS    := ${F95FLAGS} ${DEBUG_F95FLAGS}
else
   F95FLAGS    := ${F95FLAGS} -O3
   CFLAGS      := ${CFLAGS} ${RELEASE_CFLAGS} -O3
   CXXFLAGS    := ${CXXFLAGS} ${RELEASE_CXXFLAGS} -O3 
   LDFLAGS     := ${LD_FLAGS} ${RELEASE_LDFLAGS} -O3
endif

CFLAGS   := ${CFLAGS}   ${DEFS}
CXXFLAGS := ${CXXFLAGS} ${DEFS}

OUTPUT := libopt.a

all: $(OUTPUT) 

SRCS := dbrent.F90  df1dim.F90  frprmn.F90  linmin.F90

LIBOBJS := $(addsuffix .o,$(basename ${SRCS}))

.PHONY: clean cleanall

${OUTPUT}: ${SRCS}
	${F95} -c ${SRCS}
	ar rvu ${OUTPUT} ${LIBOBJS}
	ranlib ${OUTPUT}


clean:
	- rm -f ${LIBOBJS}
	- rm -f ${OUTPUT}


cleanall: clean




