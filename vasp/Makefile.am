#
# Version: $Id$
#
if DOPYTHON
include $(top_srcdir)/Files

subpkgpyexec_PYTHON = __init__.py _vasp.py relaxation.py mesh.py
subpkgpyexecdir = $(pkgpyexecdir)/vasp
__init__.py: 
	@echo "from _vasp import *"  > __init__.py; 
clean-local:
	-rm -f __init__.py
endif
