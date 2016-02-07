# what flags you want to pass to the Fortran compiler & linker
AM_FFLAGS = -ffixed-line-length-132 -fd-lines-as-comments -fimplicit-none -m64 -O3 -Wall

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
lib_LIBRARIES = libgw14.a
libgw14_a_SOURCES = getword14.f gw14.inc g2.f g2.inc p2-linux32.inc p2-linux64.inc m2-p2.inc m2.inc p2_mem.inc p2.f p2_sys.f p2_io.f p2_str.f p2_mem-@target@.f m2_io-@target@.f m2_mem-@target@.f m2_sys-@target@.f c2-@target@.c

bin_PROGRAMS = gap
gap_SOURCES = gap6.5.0-rme.for average.fcm
gap_LDADD = libgw14.a ../libccp4/fortran/.libs/libccp4f.a ../libccp4/ccp4/.libs/libccp4c.a
