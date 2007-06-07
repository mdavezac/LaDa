CC     := gcc
CXX    := gcc
F77    := g77
LD     := gcc
AR     := ar rc
RANLIB := ranlib
DEBUG = NO

FFLAGS     := 
CFLAGS     := 
CXXFLAGS   := -malign-double -ffriend-injection 
DEFS       := 


DEBUG_CFLAGS     := -Wall -Wno-format -g -O0 -Wno-unknown-pragmas -fbounds-check
DEBUG_CXXFLAGS    := ${DEBUG_CFLAGS} -D_DEBUG_LADA_ 
DEBUG_LDFLAGS     := -g 

RELEASE_CXXFLAGS := ${RELEASE_CFLAGS}
RELEASE_LDFLAGS   := 
RELEASE_CFLAGS    := -Wall -Wno-unknown-pragmas -Wno-format -O3

LIBS:=   -lm  -lstdc++  -L /home/davezac/usr/lib/  -llamarck  -latat  -ltinyxml 
INCS := -I /home/davezac/usr/include 


ifeq (YES, ${DEBUG})
   CFLAGS      := ${CFLAGS} ${DEBUG_CFLAGS} 
   CXXFLAGS    := ${CXXFLAGS} ${DEBUG_CXXFLAGS}
   LDFLAGS     := ${LDFLAGS} ${DEBUG_LDFLAGS}
else
   CFLAGS      := ${CFLAGS} ${RELEASE_CFLAGS} -O3
   CXXFLAGS    := ${CXXFLAGS} ${RELEASE_CXXFLAGS} -O3 
   LDFLAGS     := ${LD_FLAGS} ${RELEASE_LDFLAGS} -O3
endif

CFLAGS   := ${CFLAGS}   ${DEFS}
CXXFLAGS := ${CXXFLAGS} ${DEFS}

OUTPUT := liblamarck.a

all: atat $(OUTPUT) test 

SRCS := main.cc

LIBSRCS := constituent_strain.cc polynome.cc structure.cc harmonic.cc cluster.cc \
	lattice.cc functional_builder.cc convex_hull.cc 

ATATSRCS := atat/xtalutil.cc atat/integer.cc atat/findsym.cc 

OBJS := $(addsuffix .o,$(basename ${SRCS}))

LIBOBJS := $(addsuffix .o,$(basename ${LIBSRCS}))

ATATOBJS := $(addsuffix .o,$(basename ${ATATSRCS}))

.PHONY: clean cleanall

${OUTPUT}: lib atat 


${OBJS} : 
	${CXX} -c ${CXXFLAGS} ${INCS} $< -o $@


${ATATOBJS} : 
	${CXX} -c ${CXXFLAGS}  ${INCS} $< -o $@


${LIBOBJS} : 
	${CXX} -c ${CXXFLAGS}  ${INCS} $< -o $@

test: ${OBJS} lib atat
	${LD} ${LDFLAGS} -o lamarck ${OBJS}  ${LIBS} ${EXTRALIBS}


lib: ${LIBOBJS}
	ar rvu ${OUTPUT} ${LIBOBJS} 
	ranlib ${OUTPUT}


atat: ${ATATOBJS}
	ar rvu atat/lib$@.a ${ATATOBJS} 
	ranlib atat/lib$@.a


clean:
	- rm -f ${OBJS} ${LIBOBJS}
	- rm -f ${OUTPUT}



cleanall:
	- rm -f ${OBJS} ${LIBOBJS} ${ATATOBJS} ${TIXMLOBJS}
	- rm -f ${OUTPUT}
	- rm -f atat/libatat.a ${OUTPUT} 



include .dependencies

