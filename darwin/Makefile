CXX    := gcc
F77    := g77
LD     := gcc
AR     := ar rc
RANLIB := ranlib
DEBUG = YES

CXXFLAGS   := 
DEFS       := -D _G_HAVE_BOOL -D ANSI_HEADERS -D HAVE_SSTREAM -D _NOLAUNCH


ifneq ( ,$(findstring gcc, $(CXX)) )
  CXXFLAGS   += -malign-double -ffriend-injection 
endif
ifeq ( $(CXX), g++ )
  CXXFLAGS   += -malign-double -ffriend-injection 
endif

DEBUG_CFLAGS     := -Wall -Wno-format -g -O0 -Wno-unknown-pragmas -fbounds-check
DEBUG_CXXFLAGS    := ${DEBUG_CFLAGS} -D_DEBUG_LADA_ 
DEBUG_LDFLAGS     := -g 

RELEASE_CXXFLAGS := ${RELEASE_CFLAGS}
RELEASE_LDFLAGS   := 
RELEASE_CFLAGS    := -Wall -Wno-unknown-pragmas -Wno-format 

LIBS:=   -lm  -lstdc++  -L /home/davezac/usr/lib/  -llamarck  -latat  -ltinyxml    -lga  -leoutils  -leo 
INCLUDES := -I ../
INCS := -I /opt/mpich/include -I /home/davezac/usr/include -I /home/davezac/usr/include/eo

ifeq ( ,$(findstring _MPI, $(CXXFLAGS)) )
  LIBS   += -L/opt/mpich/ch-p4/lib/   -lpmpich++  -lpmpich  -lmpiobject
endif

ifeq (YES, ${DEBUG})
   DEFS        += -D _DEBUG
   CXXFLAGS    := ${CXXFLAGS} ${DEBUG_CXXFLAGS} 
   LDFLAGS     := ${LDFLAGS} ${DEBUG_LDFLAGS}
else
   CXXFLAGS    += ${RELEASE_CXXFLAGS} -O3 
   LDFLAGS     += ${LD_FLAGS} ${RELEASE_LDFLAGS} 
endif

CXXFLAGS := ${CXXFLAGS} ${DEFS}

depends_dir := .dependencies/

sources := $(filter-out main.cc, $(wildcard *.cc))
objects := $(sources:.cc=.o)
refs := $(addprefix ${depends_dir}, $(sources:.cc=.d) )
refs += ${depends_dir}cemain.d
refs += ${depends_dir}pescanmain.d

vpath %.d
vpath %.d ${depends_dir}
vpath %.cc ../lamarck ../vff
vpath %.h ../lamarck

all: ce pescan

${depends_dir}%.d: CXXFLAGS += -MM
${depends_dir}%.d: %.cc
	@set -e; rm -f $@; \
	$(CXX) ${CXXFLAGS} ${INCLUDES} $< > $@.$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$ > $@; \
	rm -f $@.$$

ce: cemain.o $(filter-out pescan.o, ${objects})
	${LD} ${LDFLAGS} -o $@ $^ ${LIBS} ${EXTRALIBS}
cemain.o: CXXFLAGS += -D _CE
cemain.o: main.cc 
	${CXX} -c ${CXXFLAGS} ${INCLUDES} ${INCS} $< -o $@

${depends_dir}cemain.d: CXXFLAGS += -D_CE -MM
${depends_dir}cemain.d: main.cc
	@set -e; rm -f $@; \
	$(CXX) ${CXXFLAGS} ${INCLUDES} $< > $@.$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$ > $@; \
	sed 's,main.o,cemain.o,g' < $@ > $@.$$; \
	mv $@.$$ $@

pescan: LIBS += -lpescan -lvff -lphysics -lgslcblas -lgsl 
pescan: main.o $(filter-out ce.o, ${objects})
	${LD} ${LDFLAGS} -o $@ $^ ${LIBS} ${EXTRALIBS}
main.o: CXXFLAGS += -D _PESCAN
main.o: main.cc 
	${CXX} -c ${CXXFLAGS} ${INCLUDES} ${INCS} $< -o $@

${depends_dir}pescanmain.d: CXXFLAGS += -D_PESCAN -MM
${depends_dir}pescanmain.d: main.cc
	@set -e; rm -f $@; \
	$(CXX) ${CXXFLAGS} ${INCLUDES} $< > $@.$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$ > $@; \
	sed 's,main.o,pescanmain.o,g' < $@ > $@.$$; \
	mv $@.$$ $@ 

${objects} : ${sources}
	${CXX} -c ${CXXFLAGS} ${INCLUDES} ${INCS} $(basename $@).cc -o $@


.PHONY: clean cleanall 


clean:
	- rm -f ${objects} main.o $(refs)
	- rm -f pescan ce


commit: 
	-sed -i 's/\(const t_unsigned svn_revision = \).*\;/\1$(shell svn info | grep Revision | awk '{print $$2+1}' );/' darwin.impl.h 
	-svn ci
	-svn update .


include ${refs}
# include ${depends_dir}mainpescan.d
# include ${depends_dir}maince.d


