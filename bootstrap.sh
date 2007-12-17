#! /bin/bash

aclocal 
libtoolize --force
autoconf
autoheader
automake-1.7 -a
