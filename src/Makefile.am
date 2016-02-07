# what flags you want to pass to the C compiler & linker
AM_CPPFLAGS = -std=c++0x -O3 -I$(BOOST_DIST)
# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = deconvolute
deconvolute_SOURCES = FileParser.cpp FileReader.cpp InputFileParser.cpp StatisticsManager.cpp Holder.cpp Miller.cpp main.cpp MtzManager.cpp Matrix.cpp Vector.cpp FileParser.h FileReader.h InputFileParser.h StatisticsManager.h Holder.h Miller.h MtzManager.h Matrix.h Vector.h definitions.h misc.cpp misc.h parameters.h Deconvoluter.h Deconvoluter.cpp
deconvolute_LDADD = -lccp4io -lcctbx