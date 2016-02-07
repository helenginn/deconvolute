# what flags you want to pass to the Fortran compiler & linker
AM_FFLAGS = -ffixed-form -ffixed-line-length-132 -fd-lines-as-comments -fimplicit-none -m64 -O3 -pedantic -std=f2003 -Wall

bin_PROGRAMS = fsmelt
fsmelt_SOURCES = fsmelt.f
fsmelt_LDADD = ../gapsrc/libgw14.a ../libccp4/fortran/.libs/libccp4f.a ../libccp4/ccp4/.libs/libccp4c.a
