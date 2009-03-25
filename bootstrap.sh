#! /bin/bash

aclocal
libtoolize --force
autoconf
autoheader
automake -a
