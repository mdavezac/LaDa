#! /bin/bash

aclocal-1.7 
libtoolize --force
autoconf
autoheader
automake-1.7 -a
