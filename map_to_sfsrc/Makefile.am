# what flags you want to pass to the Fortran compiler & linker
AM_FFLAGS = -ffixed-form -m64 -O3 -Wall

bin_PROGRAMS = map_to_sf
map_to_sf_SOURCES = change_sg.f90 check_full_symm.f90 check_symmetry.f90 extract.f90 ffl_ortho_2.f90 ffl_p1.f90 ffl_p1p.f90 ffl_p2p_te.f90 fftrr.for fftutils.for map_to_sf.f90 mapdefs.f90 mapread.for read_options.for space_group.f90 space_group_setup.f90 write_mtz.for 
space_group.o space_group.mod: space_group.f90
mapdefs.o mapdefs.mod: space_group.o space_group.mod mapdefs.f90
ffl_p1.o: mapdefs.o mapdefs.mod ffl_p1.f90
ffl_p1p.o: ffl_p1p.f90
ffl_p2p_dxml.o: ffl_p2p_dxml.f90
ffl_p2p_te.o: ffl_p2p_te.f90
fftrr.o: fftrr.for
fftutils.o: fftutils.for
mapread.o:  mapdefs.o mapdefs.mod mapread.for
read_options.o: read_options.for
space_group_setup.o: space_group.o space_group.mod mapdefs.o mapdefs.mod space_group_setup.f90
write_mtz.o:  mapdefs.o mapdefs.mod write_mtz.for
change_sg.o: space_group.o space_group.mod mapdefs.o mapdefs.mod change_sg.f90
check_full_symm.o: space_group.o space_group.mod mapdefs.o mapdefs.mod \
check_full_symm.f90
check_symmetry.o: space_group.o space_group.mod mapdefs.o mapdefs.mod check_symmetry.f90
extract.o: space_group.o space_group.mod mapdefs.o mapdefs.mod extract.f90
ffl_ortho_2.o: mapdefs.o mapdefs.mod ffl_ortho_2.f90
map_to_sf_LDADD = ../libccp4/fortran/.libs/libccp4f.a ../libccp4/ccp4/.libs/libccp4c.a
