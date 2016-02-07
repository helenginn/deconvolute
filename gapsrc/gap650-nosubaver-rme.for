      PROGRAM GAP
      IMPLICIT NONE
c
c******************************************************************************
c*                                                                            *
c*                       *=======*=========*=======                           *
c*                       General Averaging Program.                           *
c*                       *=======*=========*=======                           *
c*                                                                            *
c******************************************************************************
c*                                                                            *
c*      Initial code: Dave Stuart, Lab of Mol Biophysics, Oxford. (Nov-1991)  *
c*                                                                            *
c*      Version 1.0:  D.S. 24-Febuary-1992                                    *
c*                                                                            *
c*      Version 2.0:  D.S./Jonathan Grimes  6-March-1992                      *
c*                    Incorporate crystallographic symm./matrix refinement    *
c*                    Further changes 16-March-1992, D.S./J.G.                *
c*                                                                            *
c*      Version 3.0:  More changes to matrix refinement   8-May-1992,         *
c*                     D.S./J.G.                                              *
c*                                                                            *
c*      Version 3.1:  Crystallographic symmetry implemented.                  *
c*                       25-May-1992, D.S./J.G.                               *
c*                                                                            *
c*      Version 3.2:  Fold option added for cyclic refinement of phases.      *
c*                       12-July-1992, D.S./J.G.                              *
c*                    Fold option implemented. 16-July-1992, D.S./J.G.        *
c*                                                                            *
c*      Version 3.3:  Routine to calculate ncs operators from                 *
c*                    heavy atom positions added.  17-July-1992, D.S./J.G.    * 
c*                                                                            *
c*      Version 4.0:  Automatic envelope determination added, various options *
c*                    for map flattening and truncation added. 19-July-1992,  *
c*                    D.S./J.G.                                               *
c*                    31-July_1992   Options work!                            *
c*                                                                            *
c*      Version 4.1:  1-August-1992  Some new options added - print mask etc  *
c*                    Subroutines reordered. Some input passed up to top level*
c*                    Documentation updated. Internal data streams added.     *
c*                    8 -August-1992 Pipeing of streams added + other options *
c*                    24-August-1992 DMAP option implemented.                 *
c*                    29-August-1992 Various bits added.                      *
c*                    29-August-1992 AUTO added to MENV option - not working  *
c*                    4-December-1992 OUTX option implemented.                *
c*                                                                            *
c*      Version 4.2:  30-Febuary-1993 ATOM option in MENV implemented.        *
c*                    8-March-1993 TIDY option added                          *
c*                                                                            *
c*      Version 4.3:  End of March OENV, STAT and RESE options implemented    *
c*                                                                            *
c*      Version 4.4:  November-1993, incorporate RE speed-up of de-bricking   *
c*                                                                            *
c*      Version 4.5:  Nov-1993, Changes to parser for machine compatibility   *
c*                    Dec-1993  Numerous small changes. Listing and writing   *
c*                              inverse ops.  Option to scale up 1 map or     *
c*                              scaling 2 maps together. Reset shift and scale*
c*                              factors....etc.                               *
c*                                                                            *
c*      Version 4.6:  Dec-1993  Simplified averaging option and labelled      *
c*                              envelope option intiated. Verbose, terse      *
c*                              options added. Error messages standardized.   *
c*                              Output reworked, limited to 80 columns.       *
c*                              Various statistics now work in doubleprecision*
c*                                                                            *
c*      Version 4.7:  Jan-1994  Average-many-maps option. Map labels expanded *
c*                              from 4 to eight characters.                   *
c*                                                                            *
c*      Version 4.8:  Aug-1994  FLIP option incorporated JG                   *
c*                                                                            *
c*      Version 4.9:  Sep-1995  New options in TOSS. More double precision    *
c*                    Jan-1996  DMAP option corrected, Correlation mapping    *
c*                              implemented in AVER routine                   *
c*                                                                            *
c*      Version 5.0   Feb-1996  Dynamic memory allocation added               *
c*                              Rapid compressed envelope read/write added    *
C*                              Reconstruction added (for double interp)      *
C*                              Map & mask printing options extended          *
C*                              Printed output compressed                     *
C*                        CPU timing option added                             *
C*                                                                            *
C*      Version 6.0 through to 6.3.                                           *
C*                    1997-1999 Various additions introduced by Jon Diprose   *
C*                                                                            *
C*      Version 6.4.0 Jan-2000  Additions to convert a map to a set of        *
C*                              coordinates (M2AT).  Updates to subroutines   *
C*                              go_cut and spew_map.                          *
C*                                                                            *
C*      Version 6.4.1 Feb-2000  Addition of R.E subroutine to idealise a      *
C*                              set of proper NCS symmetry operators.         *
C*                        Addition to make matrices pure rotations            *
C*                                                                            *
C*      Version 6.4.2 Apr-2002  JD changes, and from DIS to analyse peaks     *
C*                                                                            *
C*      Version 6.5.0 Apr-2010  JD changes for linux gfortran 4.1.2           *
C*                              compatibility                                 *
C*                                                                            *
c******************************************************************************
c*                                                                            *
c*      OVERVIEW OF PROGRAM FUNCTIONS                             AUGUST-1992 *
c*      =============================                             =========== *
c*                                                                            *
c*      The philosophy is simple, have a stack of maps (electron density      *
c*      or envelopes or whatever, we are egalitarian - all get the same       *
c*      treatment) stored in a scratch array that is dynamically allocated,   *
c*      this forms a database of maps. We also have a data base of non-       *
c*      crystallographic symmetry operators. Then we have a number of smaller *
c*      controller arrays that contain all the necessary definitions for each *
c*      map.                                                                  *
c*                                                                            *
c*      All the subroutines needed are included (we use ASCII version of      *
c*      the CCP4 map format to avoid libraries, machine specificity etc),     *
c*      with the exception of the parser which ....                           *
c*      As of 24-feb-92 is Robert Esnouf's parser:       GETWORD              *
c*                                         calls :       GETOPEN              *
c*                                                       KEYWORD              *
c*                                                       INTVAL               *
c*                                                       REALVAL              *
c*                                                       CHARVAL              *
c*                                                       G2CHRU               *
c*                                                       G2EQ                 *
c*                                                       P2DATE               *
c*                                                       P2TIME               *
c*                                                       P2AL                 *
c*                                                       P2DAL                *
c*                                                                            *
c*                                                                            *
c* The program is driven by keywords, these may be nested but the major       *
c* commands are all invoked from the base level keywords. The presently       *
c* available base level keywords are (in alphabetical order):                 *
c*                                                                            *
c* ADCG - (for ncs refinement) refine the rotations about a local centre of   *
c*        gravity.                                                            *
c* AENV - define and envelope as the logical .and. of two envelopes.          *
c* ASSI - assign a map in the data base to an internal stream.                *
c* ASYM - add a non-crystallographic symmetry operator (NCS) to the list of   *
c*        such operators.                                                     *
c* AVER - average a map in the data base.                                     *
c* AVEM - averages many maps from the data base.                              *
c* CELL - redefine the cell dimensions for a certain map in the data base.    *
c* CMAP - copy one map directly into another (defined or undefined).          *
c* CONV - convolute electron density with a Gaussian smearing function.       *
c* CPU  - print cpu time used, both overall and since last call.              *
c* CUTG - truncate electron density to lie between bounds CUTL and CUTU.      *
c* CUTL - define lower cutoff for electron density.                           *
c* CUTP - print cutoffs on electron density.                                  *
c* CUTU - define upper cutoff for electron density.                           *
c* DEFA - *******
c* DMAP - delete a map from the data base and free up the space it occupied.  *
c* EXIT - end the program.                                                    *
c* FLAT - set electron density to its mean value.                             *
c* FLIP - flip electron density value.                                        *
c* FMAP - list a map presently in the data base (print a selected portion)    *
c*        on a fine grid.                                                     *
c* FMAS - list a map presently in the data base (print as a mask) on a        *
c*          fine grid.                                                        *
c* FOLD - take a volume of density and try to fill the o/p map using the      *
c*        crystallographic symmetry operators.                                *
c* FSOL - calculate fractional solvent content for each pixel                 *
c* HIST - form histogram of electron density and use it to define envelope.   *
c* IDEA - idealise a set of NCS operators back to proper symmetry             *
c* IMPR - flag to use improper NCS                                            *
c* INIT - initialize the main arrays, this will delete the entire data base   *
c*        of maps.                                                            *
c* LENV - generate up a labelled envelope using NCS operators.  This sets     *
c*        the envelope value to be the value for the NCS operator used to     *
c*        to generate it.  To be used in conjunction with UNPA in double      *
c*        interpolation averaging.                                            *
c* LISY - ******
c* LMAP - list a map presently in the data base (print a selected portion).   *
c* LMAS - list a map presently in the data base (print as a mask).            *
c* LSYM - list currently defined NCS operators.                               *
c* MENV - make an envelope (create a map in the data base using certain       *
c*        geometrically constrained logical operators).                       *
c* MCSY - switch on or off crystallographic symmetry operators for a          *
c*        particular map.                                                     *
c* M2AT - convert a map to a pdb file of pseudo-atoms                         *
c* NCSY - no, don't use the crystallographic symmetry operators.              *
c* NOPI - disable piping (qv)                                                 *
c* NORE - switch off recording of printed output in a file                    *
c* NOXP - write out NCS operators for entry into GAP                          *
c* OENV - ********
c* OMAT - orthogonalise a set of NCS matrices.                                *
c* OUTX - write out all non-crystallographically related molecules to         *
c*        files with flags, for input to XPLOR                                *
c* PEAK - analyse map recursively to find peaks in decreasing order of height *
c* PIPE - pipe output of one operation to input of next                       *
c* PROP - flag to use proper NCS                                              *
c* QASS - show current assignments to internal streams                        *
c* QMAP - list the local names of all the maps in the data base.              *
c* QMEM - show how much memory remains free in the data base                  *
c* QSYM - list all currently active ncs operators                             *
c* RECO - record all printed output in a file                                 *
c* REFI - refine some of the non-crystallographic symmetry operators.         *
c* RESE - update the rholims in map header                                    *
c* RMAP - read a map into the data base from a file.                          *
c* RMAF - read a map into the data base from a file with a defined scale      *
C* RZIP - read a compressed unformatted envelope.                             *
c* SCUT - cut an electron density map using limits defined as multiples of    *
c*        the sigma level.  The density is set to its absolute value before   *
c*        the cuts are applied.                                               *
c* SCMA - *******
c* SETM - set electron density to some defined value.                         *
c* SGRP - change space group number.                                          *
c* SOLV - set solvent density to some defined value                           *
c* SHMA - *******
c* SMAP - *******
c* STAT - *******
c* STOP - *******
c* TERS - limit printed output.                                               *
C* TIDY - tidy up an envelope. Resets certain pixels solvent/protein          *
c* TIME - print time and date                                                 *
c* TITL - define a local title for this program run.                          *
c* TOSS - read in coordinates for ncs related atoms and from them derive      *
c*        least squares best values for some ncs operators                    *
c* UNPA - unpacks the averaged density for a protomer into the other NCS      *
c*        related protomers using the labelled envelope made by LENV          *
c* VERB - produce more than usual amount of printed output.                   *
c* WANG - make wang envelope.                                                 *
c* WISY - ******
c* WMAP - write a map from the data base into a file.                         *
c* WSYM - write currently defined NCS operators to a file.                    *
C* WZIP - write a compressed unformatted envelope.                            *
c* XPLO - write out NCS operators for direct entry into XPLOR                 *
c* YCSY - yes, use the crystallographic symmetry operators.                   *
c* ZSYM - cancel all the currently defines NCS operators.                     *
c*                                                                            *
c*                                                                            *
c* THE PARSER                                                                 *
c* ==========                                                                 *
c*    (1) Keywords can be upper or lower case. They can be abbreviated as     *
c*        long as they stay unique.                                    *
c*    (2) A shreik :-  !  causes rest of the current input line to be ignored *
c*    (3) Stuff enlosed in curly brackets .. { ignore me } .. is ignored.     *
c*    (4) Any equals signs are ignored.                                       *
c*    (5) The @ convention for directing input from another file works as     *
c*        expected.                                                           *
c*    (6) A command preceded by a hash character is passed through to the     *
c*        operating system.                                                   *
c*    (7) Spaces and blank lines are ignored.                                 *
c*    (8) For more information (use of parameters etc) see GETWORD.FOR source *
c*        code.                                                               *
c*                                                                            *
c******************************************************************************
c*                                                                            *
c*      PROGRAMMING NOTES                                         AUGUST-1992 *
c*      =================                                         =========== *
c*                                                                            *
c*      This program has been put together with pieces of code from           *
c*      a number of other programs written by probably more scrupulous        *
c*      workers, in particular A.Brunger, W.Kabsch, M.Rossmann. Thank you,    *
c*      - this does not of course mean the resultant programme is up to       *
c*      their standards !                                                     *
c*                                                                            *
c*      Subroutines used in the program (in alphabetical order):              *
c*      NAME            FUNCTION                                              *
c*      add_cofg   --   I have run out of steam, finish this another time!    *
c*      add_env
c*      addoxt
c*      add_sym
c*      alloc
c*      av_ass
c*      average
c*      average_m
c*      av_setup
c*      av_setup_m
c*      brickel_1
c*      brickel_3
c*      brickit
c*      cell
c*      convert
c*      con_setup
c*      convolv
c*      copy
c*      copyc
c*      copy_head
c*      cop_map
c*      cop_set
c*      del_map
c*      det
c*      dir_map
c*      domain
c*      exit_it
c*      flat_setup
c*      flatten
c*      flip_setup
c*      flip_map
c*      fold_av
c*      fold_set
c*      frac_prot
c*      frac_prot_setup
c*      frac_sol
c*      frac_sol_setup
c*      f_stat
c*      f_to_j1
c*      f_to_j3
c*      f_to_o
c*      g_c_n_i
c*      gcut_setup
c*      go_cut
c*      g_p_i_i
c*      g_p_i_i_x
c*      g_p_n_i
c*      gp3
c*      gp3_a
c*      gulp_map
c*      header
c*      histo
c*      hist_set
c*      i_box
c*      iinout
c*      index
c*      init
c*      init_map
c*      interp
c*      invers
c*      invrt
c*      irscale
c*      i_stat
c*      j1_to_j3
c*      j1_to_o
c*      j3_to_f
c*      j3_to_j1
c*      lass_add
c*      lass_get
c*      lass_init
c*      lass_q
c*      list_map
c*      list_sym
c*      list_map
c*      list_mask
c*      make_env
c*      map_atom
c*      map_atom_setup
c*      matmul
c*      matmul1
c*      matm33
c*      matrot
c*      minv
c*      mmtply
c*      mod_csy
c*      mvmul
c*      normal
c*      o_to_f
c*      o_to_j1
c*      o_to_j3
c*      outxpl
c*      pip
c*      pipe_set
c*      polar
c*      p_p_n_i
c*      pp3
c*      pp3_a
c*      pri_cut
c*      prop
c*      prtsym
c*      qkfit
c*      qmem
c*      qsym
c*      read_headers
c*      read_map
c*      read_pdb
c*      reco
c*      rem
c*      reset_av
c*      reset
c*      r_stat
c*      rtrue
c*      rigset
c*      rotmat
c*      same_grid
c*      same_map
c*      set_cutl
c*      set_cutu
c*      set_map
c*      set_sollev
c*      sharp_setup
c*      sharp
c*      sig_cut
c*      sj33ev
c*      spew_map
c*      spin
c*      spin_in
c*      spinn
c*      spin_out
c*      spin_p
c*      spin_r
c*      sterefi
c*      subaver
c*      syminp
c*      tidy_setup
c*      tidy
c*      time_get
c*      title
c*      to_from
c*      tosser
c*      vadd
c*      vdot
c*      vlen
c*      vmtply
c*      vsub
c*      vunit
c*      vvmtply
c*      write_headers
c*      write_map
c*      write_sym
c*      wrtsym
c*      wtfsec
c*      x_copy
c*      xrfrac
c*      xr_isym
c*      xr_symsy
c*      x_trans
c*      zero_sym
c*                                                                            *
c******************************************************************************
c*                                                                            *
c*       Input/output streams:                                                 *
c*      stream                                                                *
c*      1     - map i/o , coordinate i/o                                      *
c*      6     - standard output                                               *
c*      60    - copy of output (record file)                                  *
c*                                                                            *
c******************************************************************************
c*                                                                            *
c*      map geometry conventions                                              *
c*      ========================                                              *
c*                                                                            *
c*      there are some naming conventions used in the code below:             *
c*                                                                            *
c*      o - orthogonal coords in Ang.                                         *
c*      f - fractional coords                                                 * 
c*      3 - coord in 3-d map grid points                                      *
c*      1 - coord in 1-d map address space.                                   *
c*                                                                            * 
c*      XO,YO,ZO - orthogonal axes                                            *
c*      X,Y,Z    - fractional axes                                            *
c*      J3       - 3-d array pointers                                         *
c*      AJ3      - 3-d array pointers (non-integral)                          *
c*      JSINGLE  - pointer to 1-d array.                                      *
c*      MAP      - map number.                                                *
c*      Since the array is stored bricked in the scratch array there are      *
c*      further functions etc to extract brick number, element within the     *
c*      brick, etc.  We also need to take into account the offset of the map  *
c*      within the scratch array.                                             *
c*                                                                            *
c*      Some of the arrays manipulate coords while some do actually store     *
c*      or retrieve the contents of the map.                                  *
c*                                                                            *
c*                                                                            *
c******************************************************************************
c
      INCLUDE 'average.fcm'
c     =====================
c
c---- that defines all the common blocks, it should be put with every
c---- subroutine and function for safety.
c
c----LOCAL VARIABLES
c
      integer j, nopt, keyword
c
c---- open parser - second argument non-zero for initialisation
c
      call getopen(' ',1)
c
c---- print headers
c
      call header(.true.)
c
c---- initialize scratch storage etc
c
      call init
c
c---- initialize internal streams
c
      call lass_init
c
c ccp4 init
c
      call ccpfyp
c
c---set up some values for some common block variables
c
      verbose =.false.
      terse   =.false.
      nrec    =60
      record  =.false.
      pipe    =.true.
      proper  =.false.
      xplor   =.false.
      shi_map =0
      sca_map =1
c
c---- new for vms, unknown for unix 
c
      cstatus ='unknown'
c
c---- intialise centre of gravity
c
      cg_check=.true.
      do j=1,3
        cg(j)=0.0
      enddo
c
c---- treat a '=' sign in input as if it were a space
c
      call G2EQ
c
c---------------------------------------------------------------------------
c
100   nopt=keyword('AV_ROOT',
     & 'INIT!RMAP!WMAP!MENV!ZSYM!ASYM!LSYM!LMAP!TITL!AVER'//
     &'!CELL!EXIT!QMAP!YCSY!NCSY!MCSY!REFI!ADCG!FOLD!TOSS'//
     &'!CUTL!CUTU!CUTP!CUTG!FLAT!SETM!CONV!HIST!QMEM!QSYM'//
     &'!RECO!NORE!LMAS!ASSI!QASS!WANG!FMAP!FMAS!AENV!CMAP'//
     &'!TIME!PIPE!NOPI!PROP!IMPR!SCUT!DMAP!WSYM!OUTX!TIDY'//
     &'!STAT!OENV!RESE!SHMA!SCMA!SMAP!CSEA!LISY!WISY!LENV'//
     &'!TERS!VERB!DEFA!AVEM!SGRP!FLIP!XPLO!NOXP!RZIP!WZIP'//
     &'!UNPA!CPU !WBIN!PEXP!SOLV!FSOL!FPRO!RMAF!WMRH!SHMT'//
     &'!GEN !RTNC!GROW!M2AT!IDEA!OMAT!FLOO!ADDM!SHMS!CUTS'//
     &'!WBIS!PEAK',0)
c
      if(nopt.eq.1)then
        call init
        write(6,*)'INITIALIZE SCRATCH AREAS'
        if(record)write(nrec,*)'INITIALIZE SCRATCH AREAS'
      endif
      if(nopt.eq.2 )call read_map_check(.false.)
      if(nopt.eq.3 )call write_map('slow')
      if(nopt.eq.4 )call make_env
      if(nopt.eq.5 )call zero_sym
      if(nopt.eq.6 )call add_sym
      if(nopt.eq.7 )call list_sym
      if(nopt.eq.8 )call list_map(.false.)
      if(nopt.eq.9 )call title
      if(nopt.eq.10)call average
      if(nopt.eq.11)call cell
      if(nopt.eq.12)call exit_it
      if(nopt.eq.13)call dir_map
      if(nopt.eq.14)then
        write(6,200)
        if(record)write(nrec,200)
200     format(' CRYSTALLOGRAPHIC SYMMETRY SWITCHED ON')
        do j=1,mmaps
          xr_con(j)=.true.
        enddo
      endif
      if(nopt.eq.15)then
        write(6,300)
        if(record)write(nrec,300)
300     format(' CRYSTALLOGRAPHIC SYMMETRY SWITCHED OFF')
        do j=1,mmaps
          xr_con(j)=.false.
        enddo
      endif
      if(nopt.eq.16)call mod_csy
      if(nopt.eq.17)call sterefi
      if(nopt.eq.18)call add_cofg
      if(nopt.eq.19)call fold_set
      if(nopt.eq.20)call tosser
      if(nopt.eq.21)call set_cutl
      if(nopt.eq.22)call set_cutu
      if(nopt.eq.23)call pri_cut
      if(nopt.eq.24)call go_cut
      if(nopt.eq.25)call flat_setup
      if(nopt.eq.26)call set_map
      if(nopt.eq.27)call convolv
      if(nopt.eq.28)call histo
      if(nopt.eq.29)call qmem
      if(nopt.eq.30)call qsym
      if(nopt.eq.31)call reco(.true.)
      if(nopt.eq.32)call reco(.false.)
      if(nopt.eq.33)call list_mask(.false.)
      if(nopt.eq.34)call lass_add
      if(nopt.eq.35)call lass_q
      if(nopt.eq.36)then
        call sig_cut
        call convolv
        call histo
        call pip('ENVO','MAPI',.true.)
        call pip('    ','ENVO',.false.)
        call list_mask(.true.)
      endif
      if(nopt.eq.37)call list_map(.true.)
      if(nopt.eq.38)call list_mask(.true.)
      if(nopt.eq.39)call add_env(.true.)
      if(nopt.eq.40)call cop_map
      if(nopt.eq.41)call time_get
      if(nopt.eq.42)call pipe_set(.true.)
      if(nopt.eq.43)call pipe_set(.false.)
      if(nopt.eq.44)call prop(.true.)
      if(nopt.eq.45)call prop(.false.)
      if(nopt.eq.46)call sig_cut
      if(nopt.eq.47)call del_map
      if(nopt.eq.48)call write_sym
      if(nopt.eq.49)call outxpl
      if(nopt.eq.50)call tidy
      if(nopt.eq.51)call av_stat
      if(nopt.eq.52)call add_env(.false.)
      if(nopt.eq.53)call reset
      if(nopt.eq.54)call map_shift
      if(nopt.eq.55)call map_scale
      if(nopt.eq.56)call scale_map
      if(nopt.eq.57)call cell_search
      if(nopt.eq.58)call ilist_sym
      if(nopt.eq.59)call write_sym_inv
      if(nopt.eq.60)call label_env
      if(nopt.eq.61)then
        write(6,*)'OUTPUT WILL BE TERSE'
        if(record)write(nrec,*)'OUTPUT WILL BE TERSE'
        terse=.true.
        verbose=.false.
      endif
      if(nopt.eq.62)then
        write(6,*)'OUTPUT WILL BE VERBOSE'
        if(record)write(nrec,*)'OUTPUT WILL BE VERBOSE'
        terse=.false.
        verbose=.true.
      endif
      if(nopt.eq.63)then
        write(6,*)'OUTPUT WILL BE AT THE DEFAULT LEVEL'
        if(record)write(nrec,*)'OUTPUT WILL BE AT THE DEFAULT LEVEL'
        terse=.false.
        verbose=.false.
      endif
      if(nopt.eq.64) call average_m
      if(nopt.eq.65) call sgrp
      if(nopt.eq.66) call flip_setup
      if(nopt.eq.67) then
        write(6,*)'OUTPUT MATRICES WILL BE IN XPLOR FORMAT'
        if(record)write(nrec,*)'OUTPUT MATRICES WILL BE IN XPLOR FORMAT'
        xplor=.true.
      endif
      if(nopt.eq.68) then
        write(6,*)'OUTPUT MATRICES WILL BE IN GAP FORMAT'
        if(record)write(nrec,*)'OUTPUT MATRICES WILL BE IN GAP FORMAT'
        xplor=.false.
      endif
      if(nopt.eq.69) call read_map1('fast')
      if(nopt.eq.70) call write_map('fast')
      if(nopt.eq.71) call reconst
      if(nopt.eq.72) call p2pcpu
      if(nopt.eq.73) call write_mapb(.false.)
      if(nopt.eq.74) call point_exp
      if(nopt.eq.75) call set_sollev
      if(nopt.eq.76) call frac_sol
      if(nopt.eq.77) call frac_prot
      if(nopt.eq.78) call read_map_check(.true.)
      if(nopt.eq.79) call write_mean
      if(nopt.eq.80) call map_shift_true
      if(nopt.eq.81) call gen
      if(nopt.eq.82) call rtnc
      if(nopt.eq.83) call grow
      if(nopt.eq.84) call map_atom
      if(nopt.eq.85) call idealise_mats
      if(nopt.eq.86) call ortho_mats
      if(nopt.eq.87) call flood_setup
      if(nopt.eq.88) call add_map
      if(nopt.eq.89) call map_shift_sol
      if(nopt.eq.90) call set_cut_sol
      if(nopt.eq.91) call write_mapb(.true.)
      if(nopt.eq.92) call peak_setup
c
c
c
      goto 100
      end
c
c---- end main routine
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE add_cofg
c     ===================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      integer i
      real realval
c
      write(6,*)'CENTRE OF GRAVITY REMOVAL INVOKED'
      if(record)write(nrec,*)' CENTRE OF GRAVITY REMOVAL INVOKED'
      if(verbose)then
        write(6,*)
        write(6,*) 'Please look at your map and decide on a vector'
        write(6,*) 'that will position approximately the origin of'
        write(6,*) 'of coordinate space on c.of.g of the object...'
        write(6,*) 'this will hopefully reduce errors when applying'
        write(6,*) 'rotation matrices. Input now as 3 real numbers'
        write(6,*) '- Angs wrt orthogonal system.'
      endif
c

c
      do i=1,3
        cg(i)=0.0
        cg(i)=realval('AV_C_OF_G_VEC:',-100000.0,100000.0,0.0)
      enddo
c
      write(6,*)' C of G :',cg
      if(record)write(nrec,*)' C of G :',cg
c
      cg_check=.true.
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE add_env(flag)
c     ========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c
      integer     map1,map2
      character*8 label1,label2
c
      integer     iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
      character*8 G2CHRU,ch_temp
      integer     n_reset
      integer*2   irho,gp3_a
      logical     lass_get,flag
      logical     same,same_grid,same_map
      integer     iindex(3)
      real        aindex(3),x(3),x1(3)
      external    gp3, gp3_a, same_grid, same_map
c
c      equivalence  (i_temp,ch_temp)
c
      same=.true.
c
      if(.not.lass_get(ch_temp,'ENVI'))then
        ch_temp=G2CHRU('ADD_ENV:ENTER_ENV_LABEL',1,8,0)
        if(verbose)then
          write(6,5) ch_temp
          if(record)write(nrec,5) ch_temp
5         format(' First envelope,  label:', a)
        endif
        call lass_set(ch_temp,'ENVI')
      endif
      label1=ch_temp
c
      if(.not.lass_get(ch_temp,'ENVO'))then
        ch_temp=G2CHRU('ADD_ENV:ENTER_ENV_LABEL',1,84,0)
        if(verbose)then
          write(6,10) ch_temp
          if(record)write(nrec,10) ch_temp
10        format(' Second envelope, label:', a)
        endif
        call lass_set(ch_temp,'ENV2')
      endif
      label2=ch_temp
c
c
      write(6,20)label1,label2
      if(record)write(nrec,20)label1,label2
20    format(' ADD TWO ENVELOPES, logical names ',a,' and ',a)
c
      call av_ass(map1, label1)
c
      call av_ass(map2, label2)
c
c
      if(.not.defined(map1))then
        write(6,40)
        if(record)write(nrec,40)
40      format(' %ADD_ENV-ERR: First envelope empty')
        goto 100
      endif
c
      if(.not.defined(map2))then
        write(6,42)
        if(record)write(nrec,42)
42      format(' %ADD_ENV-ERR: Second envelope empty')
        goto 100
      endif
c
      n_reset=0
c
      same=same_grid(map1,map2)
c
      if(.not.same)then
        write(6,44)
        if(record)write(nrec,44)
44      format(' %ADD_ENV-ERR: Envelopes not on same grid')
        goto 100
      endif
c
      same=same_map(map1,map2)
c
      if(.not.same)then
        write(6,46)
        if(record)write(nrec,46)
46      format(' **WARNING** envelopes do not cover the same volume')
      endif
c
      if(flag)then
        write(6,51)
        if(record)write(nrec,51)
51      format(' Forming the logical .AND. of envelopes')
      else
        write(6,52)
        if(record)write(nrec,52)
52      format(' Forming the logical .OR. of envelopes')
      endif
c        
c
c---- now skip thro' resetting
c---- loop over bricks
c
      do iz=1,n_brick(3,map2)
c
        jz=(iz-1)*brick(3,map2)
c
        do iy=1,n_brick(2,map2)
c
          jy=(iy-1)*brick(2,map2)
c
          do ix=1,n_brick(1,map2)
c
            jx=(ix-1)*brick(1,map2)
c
c---- now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map2)
c
              if(jjz.le.nx(3,map2))then
c
                do jjy=jy+1,jy+brick(2,map2)
c
                  if(jjy.le.nx(2,map2))then
c
                    do jjx=jx+1,jx+brick(1,map2)
c               
                      if(jjx.le.nx(1,map2))then
c
      iindex(1)=jjx
      iindex(2)=jjy
      iindex(3)=jjz
c
c----if maps not with same vol go and get mapped index for map1
c----gp3 and gp3_a cover for pixels out of volume and print 1 warning
c
      if(.not.same)then
        call j3_to_f(iindex,x1,map2)
        call f_to_o(x1,x,map2)
        call o_to_j3(x,aindex,iindex,map1)
      endif
c
      irho=gp3_a(map1,iindex)
c
c
c-- this is logical .and.
      if(flag)then
        if(irho.eq.0)then
          call pp3(map2,jjx,jjy,jjz,irho)
          n_reset=n_reset+1
        endif
      endif
c
c this is logical .or.
      if(.not.flag)then
        if(irho.ne.0)then
          call pp3(map2,jjx,jjy,jjz,irho)
          n_reset=n_reset+1
        endif
      endif
c
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c
      write(6,*)'Envelope calculated,',n_reset,' pixels were reset'
      if(record)then
        write(nrec,*)'Envelope calculated,',n_reset,
     &      ' pixels were reset'
      endif
c
100   return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE add_map
c     ========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c
      integer     map1,map2,env1,env2,nsummed,isum,nover,nunder,
     &            iz,jz,iy,jy,ix,jx,jjz,jjy,jjx,iindex(3)
      integer*2   irho1,irho2,isum2,gp3_a
      real        aindex(3),x(3),x1(3)
      character*8 g2chru,lmap1,lmap2,lenv1,lenv2
      logical     lass_get,icheck,ocheck,same,same_grid,same_map
c
c
c
      if(.not.lass_get(lmap1,'MAPI'))then
        lmap1=G2CHRU('ADD_MAP:ENTER_MAP_LABEL',1,8,0)
        if(verbose)then
          write(6,10)lmap1
          if(record)write(nrec,10)lmap1
        endif
        call lass_set(lmap1,'MAPI')
      endif
c
c
c
      if(.not.lass_get(lmap2,'MAPO'))then
        lmap2=G2CHRU('ADD_MAP:ENTER_MAP_LABEL',1,8,0)
        if(verbose)then
          write(6,20)lmap2
          if(record)write(nrec,20)lmap2
        endif
        call lass_set(lmap2,'MAPO')
      endif
c
c
c
      if(.not.lass_get(lenv1,'ENVI'))then
        lenv1=G2CHRU('ADD_MAP:ENTER_ENV_LABEL',1,8,0)
        if(verbose)then
          write(6,30)lenv1
          if(record)write(nrec,30)lenv1
        endif
        call lass_set(lenv1,'ENVI')
      endif
c
c
c
      if(.not.lass_get(lenv2,'ENVO'))then
        lenv2=G2CHRU('ADD_MAP:ENTER_ENV_LABEL',1,8,0)
        if(verbose)then
          write(6,40)lenv2
          if(record)write(nrec,40)lenv2
        endif
        call lass_set(lenv2,'ENVO')
      endif
c
10    format(' Input map,       label:', a)
20    format(' Output map,      label:', a)
30    format(' Input envelope,  label:', a)
40    format(' Output envelope, label:', a)
c
c
c
      write(6,50)lmap1,lmap2,lenv1,lenv2
      if(record)write(nrec,50)lmap1,lmap2,lenv1,lenv2
50    format(' ADD TWO MAPS',/,' MAPI: ',a,' MAPO: ',a,
     &            ' ENVI: ',a,' ENVO: ',a)
c
c
c
      call av_ass(map1,lmap1)
      if(.not.defined(map1))then
        write(6,60)
        if(record)write(nrec,60)
        return
      endif
c
c
c
      call av_ass(map2,lmap2)
      if(.not.defined(map2))then
        write(6,70)
        if(record)write(nrec,70)
        return
      endif
c
c
c
      icheck=(lenv1.ne.'OFF')
      if(icheck)then
        call av_ass(env1,lenv1)
        if(.not.defined(env1))then
          write(6,80)
          if(record)write(nrec,80)
          icheck=.false.
        endif
      endif
c
c
c
      ocheck=(lenv2.ne.'OFF')
      if(ocheck)then
        call av_ass(env2,lenv2)
        if(.not.defined(env2))then
          write(6,90)
          if(record)write(nrec,90)
          ocheck=.false.
        endif
      endif
c
60    format(' %ADD_MAP-ERR: Input map empty - Exiting')
70    format(' %ADD_MAP-ERR: Output map empty - Exiting')
80    format(' %ADD_MAP-ERR: Input env empty - Switching envi off')
90    format(' %ADD_MAP-ERR: Output env empty- Switching envo off')
c
c
c
      same=same_grid(map1,map2)
      if(.not.same)then
        write(6,100)
        if(record)write(nrec,100)
        return
      endif
100   format(' %ADD_MAP-ERR: Maps not on same grid - Exiting')
c
c
c
      same=same_map(map1,map2)
      if(.not.same)then
        write(6,110)
        if(record)write(nrec,110)
      endif
110   format(' **WARNING** envelopes do not cover the same volume')
c
c
c
      nsummed=0
c
c--- Loop over bricks in map2
c
      do iz=1,n_brick(3,map2)
        jz=(iz-1)*brick(3,map2)
        do iy=1,n_brick(2,map2)
          jy=(iy-1)*brick(2,map2)
          do ix=1,n_brick(1,map2)
            jx=(ix-1)*brick(1,map2)
c
c--- Now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map2)
              if(jjz.le.nx(3,map2))then
                do jjy=jy+1,jy+brick(2,map2)
                  if(jjy.le.nx(2,map2))then
                    do jjx=jx+1,jx+brick(1,map2)
                      if(jjx.le.nx(1,map2))then
c
c
c
      iindex(1)=jjx
      iindex(2)=jjy
      iindex(3)=jjz
c
c---- If not within envo, next pixel
c
      if(ocheck)then
        if(gp3_a(env2,iindex).le.0)goto200
      endif
c
c
c
      irho2=gp3_a(map2,iindex)
c
c---- If maps not with same vol go and get mapped index for map1
c
      if(.not.same)then
        call j3_to_f(iindex,x1,map2)
        call f_to_o(x1,x,map2)
        call o_to_j3(x,aindex,iindex,map1)
      endif
c
c---- If not within envi, next pixel
c
      if(icheck)then
        if(gp3_a(env1,iindex).le.0)goto200
      endif
c
c
c
      irho1=gp3_a(map1,iindex)
c
c---- Check for possible overflow
c
      isum=irho1+irho2
      if(isum.gt.32767)then
        nover=nover+1
        isum2=32767
      elseif(isum.lt.-32767)then
        nunder=nunder+1
        isum2=-32767
      else
        isum2=isum
      endif
c
c---- Put the sum - use jjx,jjy,jjz in case index has changed
c
      call pp3(map2,jjx,jjy,jjz,isum)
      nsummed=nsummed+1
c
c
c
200                   endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c
c
      write(6,300)nsummed
      if(record)write(nrec,300)nsummed
300   format(i1,' pixels were summed')
c
c
c
      if(nover.gt.0)then
        write(6,310)nover
        if(record)write(6,310)nover
      endif
      if(nunder.gt.0)then
        write(6,320)nunder
        if(record)write(6,320)nunder
      endif
310   format(' **** WARNING **** ',i1,' overflows occured')
320   format(' **** WARNING **** ',i1,' underflows occured')
c
c
c
      return
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE add_sym
c     ==================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c     =====================
c
      real    R_local(3,3),t(3),rinv(3,3),tinv(3),det,realval
      integer i,j,k
c
      do i=1,3
        do j=1,3
          r_local(i,j)=realval('AV_ADD_SYMM:ENTER_ROT_MAT',-100.0,
     &        100.0,0)
        enddo
      enddo
c
      do i=1,3
        t(i)=realval('AV_ADD_SYMM:ENTER_T_VECTOR',-100000.0,100000.0,
     &      0)
      enddo
c
      nsym=nsym+1
c
      call invers(r_local,rinv,det)
c
      if(abs(abs(det)-1.0).gt.0.99)then
        write(6,*)'**WARNING** det .ne. 1'
        if(record)write(nrec,*)'**WARNING** det .ne. 1'
      endif
c
      do j=1,3
        tinv(j)=-t(j)
      enddo
      call copy(9,r_local,ops(1,1,nsym))
      call copy(9,rinv,ops_inv(1,1,nsym))
      call copy(3,t,vecs(1,nsym))
      if(.not.terse)then
        write(6,100)nsym,((r_local(j,k),k=1,3),j=1,3),t
c     &                  ,((rinv(j,k),k=1,3),j=1,3),tinv
        if(record)write(nrec,100)nsym,((r_local(j,k),k=1,3),j=1,3),t
c     &                  ,((rinv(j,k),k=1,3),j=1,3),tinv
100     format(10x,i6,' Symm ops now in use',/,2(3(10x,3f11.5,/),
     &           10x,3f11.2))
      endif
c
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE ADDOXT
c     =================
      IMPLICIT NONE
c*******************************************************************************
c
C        - A little utility program written by R. ESNOUF to add the 
c           terminal oxygen to a peptide chain.....we have altered it
c           for our needs.
C
C
C       written by Robert Esnouf, August 1989
c
C******************************************************************************
c
      INCLUDE 'average.fcm'
c     ==================== 
c
c
C The contents of these variables and their types follows:
c
c
      CHARACTER*5     C5
      REAL*4          V1(3),V2(3),V3(3)
      COMMON /SILENT/ SILENT
      LOGICAL         SILENT
      integer nca, nc, no, n
      real*4  fac,vdot
      external vdot
c
c
      SILENT=.TRUE.
c
c
      IF (P_NA.EQ.0) THEN
        WRITE(6,1000) 'No PDB data available'
        RETURN
      ENDIF
c
C Find the numbers of the CA, C and O atoms
c
      NCA=0
      NC=0
      NO=0
      N=P_NA
c
10    IF (N.EQ.0.OR.P_RNUM(N).NE.P_RNUM(P_NA)) GOTO 20
      IF (P_ANAM(N).EQ.' CA  ') NCA=N
      IF (P_ANAM(N).EQ.' C   ') NC=N
      IF (P_ANAM(N).EQ.' O   ') NO=N
      IF ((P_ANAM(N).EQ.' OXT ').OR.(P_ANAM(N).EQ.' OT  ')) THEN
        if(.not.terse)then
          WRITE(6,1000)'Terminal oxygen already present'
          if(record)WRITE(nrec,1000)'Terminal oxygen already present'
        endif
        RETURN
      ENDIF
      N=N-1
      GOTO 10
c
20    IF (NCA.EQ.0) THEN
        WRITE(6,1000) 'Can''t find final alpha-C'
        if(record)WRITE(nrec,1000) 'Can''t find final alpha-C'
        RETURN
      ELSE IF (NC.EQ.0) THEN
        WRITE(6,1000) 'Can''t find final carbonyl-C'
        if(record)WRITE(nrec,1000) 'Can''t find final carbonyl-C'
        RETURN
      ELSE IF (NO.EQ.0) THEN
        WRITE(6,1000) 'Can''t find final carbonyl-O'
        if(record)WRITE(nrec,1000) 'Can''t find final carbonyl-O'
        RETURN
      ELSE IF (P_RECTYP(P_NA).NE.'ATOM  ') THEN
        WRITE(6,1000) 'Last residue is HETATM'
        if(record)WRITE(nrec,1000) 'Last residue is HETATM'
        RETURN
      ENDIF
c
      P_NA=P_NA+1
      P_RATM(P_NR)=P_RATM(P_NR)+1
      C5=P_ANUM(P_NA-1)
      P_RECTYP(P_NA)=P_RECTYP(P_NA-1)
      P_ANAM(P_NA)=' OXT '
      P_ANAM(P_NA-1)=' O   '
      P_RNAM(P_NA)=P_RNAM(P_NA-1)
      P_RNUM(P_NA)=P_RNUM(P_NA-1)
      P_OCC(P_NA)=1.0
      P_BFAC(P_NA)=20.0
c
C Give the atom a new 'number'
c
      N=5
30    IF (N.GT.1.AND.C5(N:N).EQ.'Z') THEN
        C5(N:N)='A'
        N=N-1
        GOTO 30
      ELSE IF (N.GT.1.AND.C5(N:N).EQ.'9') THEN
        C5(N:N)='0'
        N=N-1
        GOTO 30
      ELSE IF (C5(N:N).EQ.' ') THEN
        C5(N:N)='1'
      ELSE
        C5(N:N)=CHAR(ICHAR(C5(N:N))+1)
      ENDIF
      P_ANUM(P_NA)=C5
c
C Now calculate the new coordinates
c
      CALL VSUB(V1,P_XYZ(1,NO ),P_XYZ(1,NC))
      CALL VSUB(V2,P_XYZ(1,NCA),P_XYZ(1,NC))
      CALL VUNIT(V2,V2)
      FAC=2.0*VDOT(V1,V2)
      V2(1)=FAC*V2(1)
      V2(2)=FAC*V2(2)
      V2(3)=FAC*V2(3)
      CALL VSUB(V3,V2,V1)
      CALL VADD(P_XYZ(1,P_NA),P_XYZ(1,NC),V3)
c
c
      RETURN
c
1000  FORMAT(' %ADDOXT-ERR: **',A,'**')
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE alloc(map)
c     =====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      integer map
C RME: 23/4/2010: Made index to scratch space explicitly integer*8
      integer*8 j8
c
      if(defined(map)) goto 999
c
c
c  allocate memory for this map. Use P2AL from parser library.
      print *,'Allocating memory'
      ns(map) = P2AL(scratch(1),scratch(2),nx_tot(map))
      print *,'Returned from allocating memory'
      if(ns(map).eq.0) goto 888
c
c       nxtpix=nxtpix + nx_tot(map)
      ne(map)=ns(map)+nx_tot(map)
      defined(map)=.true.
c
      do j8=ns(map),ne(map)
        scratch(j8)=0
      enddo
c
      if(.not.terse)then
        write(6,1000)npix(map),map
        if(record)write(nrec,1000)npix(map),map
1000    format(' Dynamically allocated',i10,' pixels for map',i4)
      endif
      return
c
888   write(6,1010)npix(map)
      if(record)write(nrec,1010)npix(map)
1010  format(' %ALLOC-ERR: Unable to allocate memory',i10,
     &  ' elements required')
      return
c
999   write(6,1020)
      if(record)write(nrec,1020)
1020  format(' %ALLOC-ERR: Space already defined')
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE av_ass(mapn,mapl)
c     ============================
      IMPLICIT NONE
c
      integer mapn,j
      character*(*) mapl
      logical ok
c
      INCLUDE 'average.fcm'
c     =====================
c
      ok = .false.
      if(nmaps.le.0) goto 5
      do j=1,nmaps
        if(mapl.eq.labels(j))then
          ok=.true.
          mapn=j
        endif
      enddo
5     if(.not.ok)then
        if(.not.terse)then
          write(6,10)mapl
          if(record)write(nrec,10)mapl
10        format(' No existing map labelled ',a)
        endif
c---- put down a marker for a new map
        nmaps=nmaps+1
        defined(nmaps)=.false.
        labels(nmaps)=mapl
        mapn=nmaps
        write(6,15)mapn
        if(record)write(nrec,15)mapn
15      format(' New map made, number:',i5)
      else
c---- assign an old map
        write(6,20)mapn,mapl
        if(record)write(nrec,20)mapn,mapl
20      format(' Map #',i4,' assigned to ',a)
      endif
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE av_setup(map1,map2,env1,env2,nsym1,nsym2,i_check,
     &   o_check,use_out,no_trans,trans,itype,quit)
c     ============================================================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- Local variables
c
      integer      map1,map2,env1,env2,nsym1,nsym2,itype
      logical      i_check,o_check,use_out,no_trans,trans,quit
c
      real         Xrvol
      integer      nopt,keyword,j,intval,ics_local
      logical      lass_get
      character*8  G2CHRU,imap,omap,ienv,oenv
      character*13 mess(0:1)
      data mess/'on all pairs ','from the mean'/
c
c--- Initialise variables
c
      nsym1=1
      nsym2=nsym
      i_check=.false.
      o_check=.false.
      use_out=.false.
      no_trans=.false.
      trans=.true.
      itype=11
      quit=.false.
c
      correl=.false.
      istat=1
c
c--- Get map labels
c
      if(.not.lass_get(imap,'MAPI'))then
        imap=G2CHRU('AV_AVER:INPUT_MAP',1,8,0)
        call lass_set(imap,'MAPI')
      endif
c
      if(.not.lass_get(omap,'MAPO'))then
        omap=G2CHRU('AV_AVER:OUTPUT_MAP',1,8,0)
        call lass_set(omap,'MAPO')
      endif
c
      if(.not.lass_get(ienv,'ENVI'))then
        ienv=G2CHRU('AV_AVER:INPUT_ENV',1,8,0)
        call lass_set(ienv,'ENVI')
      endif
c
      if(.not.lass_get(oenv,'ENVO'))then
        oenv=G2CHRU('AV_AVER:OUTPUT_ENV',1,8,0)
        call lass_set(oenv,'ENVO')
      endif
c
c--- Report map and envelope labels
c
      if(verbose)then
        write(6,100) imap,omap,ienv,oenv
        if(record)write(nrec,100) imap,omap,ienv,oenv
      else
        write(6,110) imap,omap,ienv,oenv
        if(record)write(nrec,110) imap,omap,ienv,oenv
      endif
c
100   format(/' AVERAGE ELECTRON DENSITY',/,
     &        ' ========================'
     &       /' Input map:       ',a,
     &       /' Output map:      ',a,
     &       /' Input envelope:  ',a,
     &       /' Output envelope: ',a)
110   format(' AVERAGE ELECTRON DENSITY',/,
     &       ' I/P map:',a,'  O/P map:',a,'  I/P env:',a,
     &       '  O/P env:',a)
c
c--- Setup envelope logicals
c
      i_check=ienv.ne.'OFF'
      o_check=oenv.ne.'OFF'
c
c--- Report on envelope status
c
      if(verbose)then
        if(i_check) write(6,120)
        if(.not.i_check)write(6,130)
        if(o_check) write(6,140)
        if(.not.o_check)write(6,150)
c
        if(record) then
          if(i_check) write(nrec,120)
          if(.not.i_check)write(nrec,130)
          if(o_check) write(nrec,140)
          if(.not.o_check)write(nrec,150)
        endif
      endif
c
120   format(' Input map envelope filtered')
130   format(' Input map not filtered')
140   format(' Output map envelope filtered')
150   format(' Output map not filtered')
c
c--- Get slot numbers for all maps
c--- map1 should exist by now!
c
      call av_ass(map1,imap)
      if(.not.defined(map1))then
        write(6,160)
        if(record)write(nrec,160)
        return
      endif
c
c--- Check status of output file
c
      call av_ass(map2,omap)
      if(.not.terse)then
        if(defined(map2))write(6,170)
        if(.not.defined(map2))write(6,180)
c
        if(record)then
          if(defined(map2))write(nrec,170)
          if(.not.defined(map2))write(nrec,180)
        endif
      endif
c
c---- Set up allocations for envelopes
c
      if(i_check)then
        call av_ass(env1,ienv)
        if(.not.defined(env1))then
          i_check=.false.
          write(6,190)
          if(record)write(6,190)
        endif
      endif
c
      if(o_check)then
        call av_ass(env2,oenv)
        if(.not.defined(env2))then
          o_check=.false.
          write(6,200)
          if(record)write(6,200)
        endif
      endif
c
160   format(' %AVERAGE-ERR: Input map empty')
170   format(' Output map exists: will overwrite')
180   format(' Make output from i/p template')
190   format(' %AVERAGE-ERR: Input envelope empty - switched off')
200   format(' %AVERAGE-ERR: Output envelope empty - switched off')
c
c--- Set up output map header information
c
      if(.not.defined(map2))then
        call copy_head(map1,map2)
c
c---- Check if we want to redefine the output grid (eg coarser than input)
c
300     nopt=keyword
     &   ('AV_AVERAGE_UPDATE','UPDA!RESE!NOTR!USEO!BACK!QLIM!GO  ',0)
c
c---- UPDA
c
        if(nopt.eq.1) then
          do j=1,3
            nstart(j,map2)=intval('AV_AVERAGE:NSTART',-10000,10000,0)
            nend(j,map2)  =intval('AV_AVERAGE:NEND',-10000,10000,0)
            nunit(j,map2) =intval('AV_AVERAGE:NUNIT',-10000,10000,0)
c
            nx(j,map2)=nend(j,map2)-nstart(j,map2)+1
            norg(j,map2)=1-nstart(j,map2)
c
            if(.not.terse)then
              write(6,310)
     &          nstart(j,map2),nend(j,map2),nx(j,map2),nunit(j,map2)
              if(record)write(nrec,310)
     &          nstart(j,map2),nend(j,map2),nx(j,map2),nunit(j,map2)
            endif
          enddo
c
          npix(map2)=nx(1,map2)*nx(2,map2)*nx(3,map2)
c
          write(6,320)npix(map2)
          if(record)write(nrec,320)npix(map2)
        endif
c
310     format(' ',4i5)
320     format(' There will be',i10,' pixels in the output map')
c
c---- RESE
c
        if(nopt.eq.2)then
          do j=1,3
            iuvw(j,map2)=intval('AV_AVERAGE:ENTER_AXIS',1,3,0)
          enddo
        endif
c
c---- NOTR
c
        if(nopt.eq.3) then
          no_trans=.true.
          trans=.false.
        endif
c
c---- USEO
c
        if(nopt.eq.4) then
          use_out=.true.
        endif
c
c---- BACK
c
        if(NOPT.eq.5)then
          quit=.true.
          return
        endif
c
c---- QLIM
c
        if(NOPT.eq.6)then
          write(6,330)(nstart(j,map2),nend(j,map2),nunit(j,map2),
     &      nx(j,map2),norg(j,map2),j=1,3)
        endif
c
330     format( 3(' START:',i5,' END:',i5,' GRID:',i5,
     &            ' POINTS:',i5,' ORIGIN:',i5,/) )
c
c---- GO
c
        if(nopt.ne.7)goto 300
c
c---- End of user interaction about output map
c---- Set up cell data
c
        call XRfrac(XRcell(1,map2),XRtr(1,1,map2),XRintr(1,1,map2),
     &              XRvol,.true.)
c
c---- Set up symmetry data
c
        call XR_symsy(XRcell(1,map2),lgrp(map2),ics_local,
     &                n_xsym(map2),XR_sym(1,map2) )
c
c---- Set up integerised CSO's:
c
        call XR_isym(map2)
c
c---- Set up pointers indicating which symm op to try first
c
        do j=1,maxsym
          last_sym(j,map2)=1
        enddo
c
c---- Set up bricking data
c
        call brickit(map2)
c
c---- Allocate space for output map
c
        call alloc(map2)
c
      endif
c
c--- Set up symmetry operators to be used
c
400   nopt=keyword('AV_AVERAGE_SYMM',
     &  'ALL !STAR!END !GO  !INTP!CORQ!CORA!USEO!CORM!NCOR',0)
c
      if(nopt.eq.1)then
        nsym1=1
        nsym2=nsym
      endif
      if(nopt.eq.2)nsym1=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.eq.3)nsym2=intval('AV_AVERAGE_SYMM',1,nsym,0)
c
      if(nopt.eq.5)itype=intval('AV_AVERAGE_INTP',1,64,0)
c
      if(nopt.eq.6)istat=1
      if(nopt.eq.7)istat=0
c
      if(NOPT.eq.8) then
        use_out=.true.
      endif
c
      if(nopt.eq.9) then
        correl=.true.
        write(6,410)
        if(record)write(nrec,410)
      endif
410   format('CORRELATION MAPPING SWITCHED ON')
c
      if(nopt.eq.10)then
        correl=.false.
        write(6,420)
        if(record)write(nrec,420)
      endif
420   format('CORRELATION MAPPING SWITCHED OFF')
c
      if(nopt.ne.4)goto 400
c
c--- Status reports and sanity checks
c--- Status report on sym ops and whether we're using the output map
c
      if(use_out)then
        write(6,500)nsym1,nsym2
        if(record)write(nrec,500)nsym1,nsym2
      else
        write(6,510)nsym1,nsym2
        if(record)write(nrec,510)nsym1,nsym2
      endif
c
500   format(' Get densities from i/p map using ncs ops',
     &             i4,' to',i4,' and average with o/p map')
510   format(' Form average density from i/p map using ncs ops',
     &             i4,' to',i4)
c
c--- Sanity check on sym ops
c
      if( (nsym1.gt.nsym2).or.(nsym2.gt.nsym) )then
        write(6,520)
        if(record)write(nrec,520)
      endif
520   format(' %AVERAGE-ERR: In choice of symmetry operators ')
c
c--- Show sym ops
c
      if(verbose)then
        do j=nsym1,nsym2
          call prtsym(6,j)
        enddo
      endif
c
c--- Sanity check on interpolation value
c
      if(itype.ne.11.and.itype.ne.8.and.itype.ne.1.and.itype.ne.64)
     &   itype=11
c
c--- Report on interpolation and cc types
c
      write(6,530)itype,mess(istat)
      if(record)write(nrec,530)itype,mess(istat)
530   format(i3,' point interpolation. Correlation coeff calculated ',
     &   a)
c
c--- And its done...
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE av_setup_m(maps,map2,envs,env2,nsym1,nsym2,i_checks,
     &   o_check,use_out,no_trans,trans,ITYPE,quit)
c     ===============================================================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c     =====================
c
c----local variables:
      integer      maps(maxsym),map2,envs(maxsym),env2,nsym1,nsym2,itype
      logical      i_checks(maxsym),o_check,use_out,no_trans,trans
c
      character*8  G2CHRU,ch_temp,omap,oenv
      real         xrvol
      integer      nopt,keyword,j,intval,ics_local,n_ops(mmaps)
      integer      n_read,n_map_1,keylist,jj
      logical      lass_get,quit
      character*13 mess(0:1)
      data mess/'on all pairs ','from the mean'/
c
c
      no_trans=.false.
      trans=.true.
      use_out=.false.
      correl=.false.
      sigma=.false.
      n_map_1=0
c
c
      n_read=intval('AV_AVER:NUMBER_OF_OPS',1,MMAPS,0)
      do j=1,n_read
        n_ops(j)=intval('AV_AVER:SYMM_OP_NUMBER',1,MAXSYM,0)
        if(j.eq.1)n_map_1=n_ops(j)
        maps(n_ops(j))=keylist('AV_AVER:MAP_LABEL',LABELS,NMAPS,0)
        envs(n_ops(j))=0
c
c----check if envelope to be used
c
        if('ENV'.eq.G2CHRU('AV_AVER:ENV/NOE',0,3,0))
     &    envs(n_ops(j))=keylist('AV_AVER:ENV_LABEL',LABELS,NMAPS,0)
        i_checks(n_ops(j))=envs(n_ops(j)).gt.0
      enddo
c
c
      write(6,1003)'Number of operators to be used:',n_read
      if(record)write(nrec,1003)'Number of operators to be used:',n_read
c
c
      do j=1,n_read
        jj=n_ops(j)
        if(i_checks(jj))write(6,1001)jj,maps(jj),envs(jj)
        if(.not.i_checks(jj))
     &    write(6,1002)jj,maps(jj)
        if(record)then
          if(i_checks(jj))write(6,1001)jj,maps(jj),envs(jj)
          if(.not.i_checks(jj))
     &      write(6,1002)jj,maps(jj)
        endif
      enddo
c
c
1001  format(' Operator',i4,' will use map',i3,' and Input env',i3)
1002  format(' Operator',i4,' will use map',i3,' and no Input env')
1003  format(' ',a,i4)
c
c
      if(.not.lass_get(ch_temp,'MAPO'))then
c
        ch_temp=G2CHRU('AV_AVER:OUTPUT_MAP',1,8,0)
c
        call lass_set(ch_temp,'MAPO')
      endif
      omap=ch_temp
c
c
      if(.not.lass_get(ch_temp,'ENVO'))then
c
        ch_temp=G2CHRU('AV_AVER:OUTPUT_ENV',1,8,0)
c
        call lass_set(ch_temp,'ENVO')
      endif
      oenv=ch_temp
c
c
      if(verbose)then
        write(6,110) omap,oenv
        if(record)write(nrec,110) omap,oenv
110     format(/' AVERAGE ELECTRON DENSITY',/,
     &          ' ========================'
     &         /' Output map:      ',a,
     &         /' Output envelope: ',a)
      else
        write(6,111) omap,oenv
        if(record)write(nrec,111) omap,oenv
111     format(' AVERAGE ELECTRON DENSITY',
     &    ' O/P map:',a,'  O/P env:',a)
      endif
C
      o_check=oenv.ne.'OFF'
c
      if(verbose)then
        if(o_check) write(6,140)
        if(.not.o_check)write(6,150)
c
        if(record) then
          if(o_check) write(nrec,140)
          if(.not.o_check)write(nrec,150)
        endif
      endif
c
140   format(' Output map envelope filtered')
150   format(' Output map not filtered')
c
c
c       get slot numbers for all maps
c
c---- check status of output file
      call av_ass(map2,omap)
      if(.not.terse)then
        if(defined(map2))write(6,170)
        if(.not.defined(map2))write(6,180)
c
        if(record)then
          if(defined(map2))write(nrec,170)
          if(.not.defined(map2))write(nrec,180)
        endif
      endif
c
170   format(' Output map exists: will overwrite')
180   format(' Make output from i/p template')
c
c---- set up output map header information
c
      if(.not.defined(map2))then
c----take the first of the input maps as the one to use to define the o/p
c----...there is a limit to the flexibility......
        call copy_head(maps(n_map_1),map2)
c
c---- check if we want to redefine the output grid (eg coarser than input)
c
        quit=.false.
c
20      nopt=keyword
     &  ('AV_AVERAGE_UPDATE','UPDA!RESE!NOTR!USEO!BACK!QLIM!GO  ',0)
c
        if(nopt.eq.1) then
          do j=1,3
            nstart(j,map2)=intval('AV_AVERAGE:NSTART',-3000,3000,0)
            nend(j,map2)  =intval('AV_AVERAGE:NEND',-3000,3000,0)
            nunit(j,map2) =intval('AV_AVERAGE:NUNIT',-3000,3000,0)
            nx(j,map2)=nend(j,map2)-nstart(j,map2)+1
            norg(j,map2)=1-nstart(j,map2)
            if(.not.terse)then
              write(6,200)nstart(j,map2),nend(j,map2),
     &                  nx(j,map2),nunit(j,map2)
              if(record)
     &          write(nrec,200)nstart(j,map2),nend(j,map2),
     &                  nx(j,map2),nunit(j,map2)
200           format(' ',4i5)
            endif
          enddo
          npix(map2)=nx(1,map2)*nx(2,map2)*nx(3,map2)
          write(6,210)npix(map2)
          if(record)write(nrec,210)npix(map2)
210       format(' There will be',i10,' pixels in the output map')
        endif
c
        if(nopt.eq.2)then
          do j=1,3
            iuvw(j,map2)=intval('AV_AVERAGE:ENTER_AXIS',1,3,0)
          enddo
        endif
c
        if(NOPT.EQ.3) then
          no_trans=.true.
          trans=.false.
        endif
c
        if(NOPT.eq.4) then
          use_out=.true.
        endif
c
        if(NOPT.eq.5)then
          quit=.true.
          return
        endif
c
        if(NOPT.eq.6)then
          write(6,214)(nstart(j,map2),nend(j,map2),nunit(j,map2),
     &      nx(j,map2),norg(j,map2),j=1,3)
214       format( 3(' START:',i4,' END:',i4,' GRID:',i4,
     &              ' POINTS:',i4,' ORIGIN:',i4,/) )
        endif
c
        if(nopt.ne.7)goto 20
c
c---- set up cell data
c
        call XRfrac(XRcell(1,map2),XRtr(1,1,map2),XRintr(1,1,map2),
     &    XRvol,.true.)
c
c---- set up symmetry data
c
        call XR_symsy(XRcell(1,map2),lgrp(map2),ics_local,
     &    n_xsym(map2),XR_sym(1,map2) )
c
c---- and integerized CSO's:
c
        call XR_isym(map2)
c
c---- also set up pointers indicating which symm op to try first
c
        do j=1,maxsym
          last_sym(j,map2)=1
        enddo
c
c
c---- set up bricking data
c
        call brickit(map2)
c
c---- allocate space for output map
c
        call alloc(map2)
c
      endif
c
c---- set up allocations for envelopes
c
      if(o_check) call av_ass(env2,oenv)
c
c
c
c---- list symmetry operators to be used
c
      itype=11
      istat=1
      nsym1=1
      nsym2=nsym
c
10    nopt=keyword
     &  ('AV_AVERAGE_SYMM','ALL !STAR!END !GO  !INTP!CORQ!CORA'//
     &          '!USEO!CORM!SIGM!DENS',0)
c
      if(nopt.eq.2)nsym1=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.eq.3)nsym2=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.eq.5)itype=intval('AV_AVERAGE_INTP',1,64,0)
      if(nopt.eq.6)istat=1
      if(nopt.eq.7)istat=0
      if(NOPT.eq.8) then
        use_out=.true.
      endif
      if(nopt.eq.9) then
        correl=.true.
        sigma=.false.
        write(6,219)
        if(record)write(nrec,219)
      endif
219   format('CORRELATION MAPPING SWITCHED ON')
      if(nopt.eq.10) then
        correl=.false.
        sigma=.true.
        write(6,218)
        if(record)write(nrec,218)
      endif
218   format('SIGMA MAPPING SWITCHED ON')
      if(nopt.eq.11) then
        correl=.false.
        sigma=.false.
        write(6,217)
        if(record)write(nrec,217)
      endif
217   format('DENSITY MAPPING SWITCHED ON')
c
      if(nopt.ne.4)goto 10
c
      if(use_out)then
        write(6,220)nsym1,nsym2
        if(record)write(nrec,220)nsym1,nsym2
220     format(' Get densities from i/p map using ncs ops',
     &             i4,' to',i4,' and average with o/p map')
      else
        write(6,221)nsym1,nsym2
        if(record)write(nrec,221)nsym1,nsym2
221     format(' Form average density from i/p map using ncs ops',
     &             i4,' to',i4)
      endif
      if( (nsym1.gt.nsym2).or.(nsym2.gt.nsym) )then
        write(6,230)
        if(record)write(nrec,230)
230     format(' %AVERAGE-ERR: In choice of symmetry operators ')
      endif
c
      if(verbose)then
        do j=nsym1,nsym2
          call prtsym(6,j)
        enddo
      endif
C----
c
      if(itype.ne.11.and.itype.ne.8.and.itype.ne.1.and.itype.ne.64)
     &    itype=11
      write(6,240)itype,mess(istat)
      if(record)write(nrec,240)itype,mess(istat)
240   format(i3,' point interpolation. Correlation coeff calculated ',
     &    a)
c
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE average
c      ==================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c      =====================
c
c--- Local variables
c
       integer   map1, map2, env1, env2, nout
       integer   iindex(3),new_index(3)
       real      aindex(3),x(3),x1(3)
       integer   np,nc
       logical   i_check,o_check,use_in,no_trans,trans,check,use_out
       logical   iiinout,iinout,quit
       integer*2 irho
       external  gp3,gp3_a
       integer*2 gp3,gp3_a
c
       integer   nsym1,nsym2,itype,ngadd
       integer   iz,jz,iy,jy,ix,jx,jjz,jjy,jjx,j
c
c--- Defaults
c
       use_in=.true.
       check=.true.
       iiinout=.false.
c
c--- Get user input
c
        call av_setup(map1,map2,env1,env2,nsym1,nsym2,i_check,o_check,
     &                use_out,no_trans,trans,itype,quit)
c
c--- Check for BACK card
c
       if(quit)then
         quit=.false.
         write(6,*)'Quitting...'
         return
       endif
c
c---- the averaging process is very simple, we are driven by the
c---- output map, for each pixel in that map (accessed in a sequential
c---- fashion, but since map is bricked we move slowly through space)
c---- we pick up all the symm related pixels in the input map and take
c---- the average.
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
       call init_map(map2)
c
c---- initialize r-factors, corr coefs, etc
c
       ngadd=0
       np=0
       nc=0
       nout=0
c
       call i_stat
c
       if(o_check)call init_map(env2)
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
       do iz=1,n_brick(3,map2)
c
         jz=(iz-1)*brick(3,map2)
c
         do iy=1,n_brick(2,map2)
c
           jy=(iy-1)*brick(2,map2)
c
           do ix=1,n_brick(1,map2)
c
             jx=(ix-1)*brick(1,map2)
c
c---- now loop within each brick
c
             do jjz=jz+1,jz+brick(3,map2)
c
             if(jjz.le.nx(3,map2))then
c
               do jjy=jy+1,jy+brick(2,map2)
c
               if(jjy.le.nx(2,map2))then
c
                 do jjx=jx+1,jx+brick(1,map2)
c                    
                 if(jjx.le.nx(1,map2))then
c
c---- Check that this pixel is in required region of the output map
c
                   if(o_check)then
                     if(gp3(env2,jjx,jjy,jjz).ne.1)goto 200
                   endif
c
c---- Get the orthogonal co-ords of this pixel
c
                   iindex(1)=jjx
                   iindex(2)=jjy
                   iindex(3)=jjz
c
                   call j3_to_f(iindex,x1,map2)
                   call f_to_o(x1,x,map2)
c
c---- Zero the counters
c----  in_temp = #pixels contributing to this output pixel
c----  sum_temp = runningdensity total, and later the average density
c
                   in_temp=0
                   sum_temp=zero
c
c---- Loop over all the required NCS operations:
c
                   do j=nsym1,nsym2
c
c---- Get the appropriate postion for this copy:
c
                     if(trans)    call spin (x,x1,ops(1,1,j),vecs(1,j))
c
                     if(no_trans) call vmtply(ops(1,1,j),x,x1)
c
                     call o_to_j3(x1,aindex,iindex,map1)
c
c---- If required check that this is in a permitted region of the Input map:
c
                     if(i_check)then
c
c---- This code allows for the envelope to be limited to the cryst. aysm.
c---- portion of the map (irrespective of what area the actual maps cover)
c---- Not entirely sure its being used, though... JMD 7/6/2001
c
                       iiinout=iinout(env1,iindex(1),iindex(2),iindex(3)
     &                            ,new_index,check,j)
c
c---- Use the envelope pixel that is closest to non-integral grid point
c---- Pixels outside envelope for i/p counted
c---- Surely this should use new_index???
c
                       if((.not.iiinout).or.
     &                          (gp3_a(env1,iindex).eq.0))then
                         nout=nout+1
                         goto 300
                       endif
c
                     endif
c
c---- Now get the interpolated electron density:
c---- If not using xtal symm....much quicker
                     if(xr_con(map1))then
c
c----- itype specifies interp: 1, 8 or 11 point.
c
                       call g_p_i_i_x(map1,irho,aindex,iindex,itype,
     &                                j,use_in,np,nc)
c
                     else
                       call g_p_i_i(map1,irho,aindex,iindex,itype,
     &                              use_in,np,nc)
                     endif
c
c---- if the pixel was found ok add it into the sums:
c
                     if(.not.use_in) then
                       nout=nout+1
                     else
                       sum_temp=sum_temp+irho
                       in_temp=in_temp+1
                       ir_temp(in_temp)=irho
                       ngadd=ngadd+1
                     endif
c
300                  continue
                   enddo
c
c---- That the end of the loop over the ncs relationships in map1.
c---- We may want to add in a component from map2, the output map.
c
                   if(use_out) then
c
c---- Get this pixel from map2
c
                     irho=gp3(map2,jjx,jjy,jjz)
                     sum_temp=sum_temp+irho
                     in_temp=in_temp+1
                     ir_temp(in_temp)=irho
                     ngadd=ngadd+1
                   endif
c
c---- Now we have gathered all the contributions.
c---- a_stat evaluates mean rho and accumulates stats for this output pixel:
c
                   call a_stat
c
                   irho=nint(sum_temp)
c
c---- Check to see if we are doing correlation mapping
c
                   if(correl) then
                     irho=0
                     if(cc_bot.ne.0.0)
     &                 irho=nint( (cc_top/cc_bot)*10000.0 )
                   endif
c
c---- And now finally store the pixel in the output map:
c
                   call pp3(map2,jjx,jjy,jjz,irho)
c
200                continue
                 endif
                 enddo
               endif
               enddo
             endif
             enddo
           enddo
         enddo
       enddo
c
c--- Phew, averaged, lets see how we did ...
c
       if(verbose)then
         write(6,1000)ngadd
         if(record)write(nrec,1000)ngadd
       else
         write(6,1010)ngadd
         if(record)write(nrec,1010)ngadd
       endif
c
       write(6,1020)np,nout,nc
       if(record)write(nrec,1020)np,nout,nc
c
1000   format(' AVERAGING FINISHED',/,
     &        ' ==================',/,
     &        i12,' pixels were inspected')
1010   format(' AVERAGING FINISHED, ',i12,' pixels were inspected')
1020   format(' # interpolated:',i9,
     &         '  outside i/p map:',i9,
     &         '  not interpolated:',i9)
c
c--- Print the statistics
c
       call p_stat
c
c--- Pipe, if we're piping
c
       if(pipe)then
         call pip('MAPO','MAPI',.true.)
         call pip('    ','MAPO',.false.)
       endif
c
c--- And its all over...
c
       return
       end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE average_m
c     ====================
      IMPLICIT NONE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c This version allows any input map to be associated with any of the NCS
c operators. This permits averaging between many different maps/crystal forms.
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- local variables
c
      integer   map2,env2,nout
      integer   iindex(3),new_index(3)
      real      aindex(3),x(3),x1(3)
      integer   np,nc
      logical   o_check,i_checks(maxsym),use_in,no_trans,trans,check,
     &          use_out
      logical   iiinout,iinout,quit
      integer*2 irho
      external  gp3,gp3_a
      integer*2 gp3,gp3_a
c
      integer   nsym1,nsym2,itype,ngadd
      integer   iz,jz,iy,jy,ix,jx,jjz,jjy,jjx,j
      integer   maps(maxsym),envs(maxsym)
c
c--- defaults
c
      use_in=.true.
      check=.true.
      iiinout=.false.
c
c--- get user input
c
      call av_setup_m(maps,map2,envs,env2,nsym1,nsym2,i_checks,o_check,
     &                use_out,no_trans,trans,ITYPE,quit)
c
c--- check for a BACK card
c
      if(quit)then
        quit=.false.
        return
      endif
c
c---- the averaging process is very simple, we are driven by the
c---- output map, for each pixel in that map (accessed in a sequential
c---- fashion, but since map is bricked we move slowly through space)
c---- we pick up all the symm related pixels in the input map and take
c---- the average.
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
      call init_map(map2)
c
c---- initialize r-factors, corr coefs, etc
c
      ngadd=0
      np=0
      nc=0
      nout=0
c
      call i_stat
c
      if(o_check)call init_map(env2)
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
      do iz=1,n_brick(3,map2)
c
        jz=(iz-1)*brick(3,map2)
c
        do iy=1,n_brick(2,map2)
c
          jy=(iy-1)*brick(2,map2)
c
          do ix=1,n_brick(1,map2)
c
            jx=(ix-1)*brick(1,map2)
c
c---- now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map2)
c
              if(jjz.le.nx(3,map2))then
c
                do jjy=jy+1,jy+brick(2,map2)
c
                  if(jjy.le.nx(2,map2))then
c
                    do jjx=jx+1,jx+brick(1,map2)
c               
                      if(jjx.le.nx(1,map2))then
c
c---- check that this pixel is in required region of the output map
c
      if(o_check)then
        if(gp3(env2,jjx,jjy,jjz).ne.1)goto 200
      endif
c
      iindex(1)=jjx
      iindex(2)=jjy
      iindex(3)=jjz
c
      call j3_to_f(iindex,x1,map2)
      call f_to_o(x1,x,map2)
c
c
      in_temp=0
      sum_temp=zero
c
c
c---- loop over all the required NCS operations:
c
      do j=nsym1,nsym2
c
c---- get the appropriate postion for this copy:
c
        if(trans) call spin (x,x1,ops(1,1,j),vecs(1,j))
c
        if(no_trans) call vmtply(ops(1,1,j),x,x1)
c
        call o_to_j3(x1,aindex,iindex,maps(j))
c
c---- if required check that this is in a permitted region of the Input map:
c
        if(i_checks(j))then
c
c---- this code allows for the envelope to be limited to the cryst. aysm.
c---- portion of the map (irrespective of what area the actual maps cover)
c
          iiinout=iinout(envs(j),iindex(1),iindex(2),
     &                   iindex(3),new_index,check,j)
c
c          (use the envelope pixel that is closest to non-integral grid point)
c
          if(.not.iiinout.or.gp3_a(envs(j),iindex).eq.0)
     &                      then
            nout=nout+1
c---- pixels outside envelope for i/p counted
            goto 300
          endif
c
        endif
c
c---- now get the interpolated electron density:
c
c---- if not using xtal symm....much quicker
        if(xr_con(maps(j)))then
          call g_p_i_i_x(maps(j),irho,aindex,iindex,itype,j,use_in
     &                                    ,np,nc)
c---- ( itype specifies interp: 1, 8 or 11 point. )
        else
          call g_p_i_i(maps(j),irho,aindex,iindex,itype,use_in,np,nc)
        endif
c
c---- if the pixel was found ok add it into the sums:
c
        if(.not.use_in) then
          nout=nout+1
        else
          sum_temp=sum_temp+irho
          in_temp=in_temp+1
          ir_temp(in_temp)=irho
          ngadd=ngadd+1
        endif
300     continue
      enddo
c
c---- thats ended the loop over the ncs relationships in input maps, we may
c---- want to add in a component from map2, the output map, for present
c---- assume scale factor of 1 between maps (use map scale option if needed).
c
      if(use_out) then
c---- check the sym ops for map 2
c
        irho=gp3(map2,jjx,jjy,jjz)
        sum_temp=sum_temp+irho
        in_temp=in_temp+1
        ir_temp(in_temp)=irho
        ngadd=ngadd+1
      endif
c
c---- now we have gathered all the contributions.
c---- a_stat evaluates mean rho and accumualtes stats for this output pixel:
c
      call a_stat
c
      irho=nint(sum_temp)
c
c---- check to see we are not doing correlation mapping 
c
      if(correl) then
        irho=0
        if(cc_bot.ne.0.0)
     &    irho=nint( (cc_top/cc_bot)*100. )
      endif
c
c---- check to see if we are doing sigma mapping 
c
      if(sigma)irho=nint(rms_delta*10.)
c
c---- and now finally store the pixel in the Output map:
c
      call pp3(map2,jjx,jjy,jjz,irho)
c
200   continue
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c----phew, averaged, lets see how we did ...
c
c999   continue
      if(verbose)then
        write(6,998)ngadd
        if(record)write(nrec,998)ngadd
998     format(' AVERAGING FINISHED',/,
     &         ' ==================',/,
     &      i9,' pixels were inspected')
      else
        write(6,997)ngadd
        if(record)write(nrec,997)ngadd
997     format(' AVERAGING FINISHED, ',i10,' pixels were inspected')
      endif
c
      write(6,1000)np,nout,nc
      if(record)write(nrec,1000)np,nout,nc
1000  format(' # interpolated:',i9,
     &       '  outside i/p map:',i9,
     &       '  not interpolated:',i9)
c
c---- print some stats
c
      call p_stat
c
      if(pipe)then
        call pip('MAPO','MAPI',.true.)
        call pip('    ','MAPO',.false.)
      endif
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


c
      FUNCTION brickel_1(j1,map)
c     ==========================
      IMPLICIT NONE
      INCLUDE 'average.fcm'
c     =====================
      integer map,j1,j3(3),brickel_1,brickel_3
c JMD!PORT: Not really external?
c      external brickel_3
c
      call j1_to_j3(j1,j3,nx(1,map))
      brickel_1=brickel_3(j3,map)
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      FUNCTION brickel_3(j3,map)
c     ==========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      integer map,j3(3),brickel_3,i3(3),i1,nb,rem3(3)
      integer n_off,b_off,j
c
c---- (i1 and i3 are brick indices not pixel indices)
c
c---- new speedier version, 23-feb-92
c
      do j=1,3
        nb     = brick(j,map)
        i3(j)  = (j3(j)-1)/nb
        rem3(j)= j3(j)-i3(j)*nb
      enddo
c
      i1= ( i3(3)*n_brick(2,map)+i3(2) )*n_brick(1,map) + i3(1)
      n_off = i1*npix_b(map)
      b_off = ( ( rem3(3)-1 )*brick(2,map) + rem3(2)-1
     &              )*brick(1,map) + rem3(1)
      brickel_3 = n_off + b_off
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE brickit(map)
c     =======================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      integer j,map
c
      npix_b(map)=ibsize**3
c
      do j=1,3
c
c---- (
c---- cheat and do counters for cryst. symm. checking, 10*size of map to
c---- ensure it is always possible to simply take remainder to place in positive
c---- portion of cell.
        n10(j,map)=10*nunit(j,map)
c---- )
        brick(j,map)=ibsize
        n_brick(j,map)=nx(j,map)/brick(j,map)
        if(n_brick(j,map)*brick(j,map).ne.nx(j,map))n_brick(j,map)=
     &     n_brick(j,map)+1
c---- bricking involves padding in memeory ... record padded map size:
        nx_pad(j,map)=n_brick(j,map)*brick(j,map)
      enddo
      nx_tot(map)=nx_pad(1,map)*nx_pad(2,map)*nx_pad(3,map)
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE cell
c     ===============
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      character*8 G2CHRU,ch_temp,mapl
      real        realval,Xrvol
c
c----update cell information
c
      integer mapn,j
c
      ch_temp=G2CHRU('AV_CELL:ENTER_MAP_NAME',1,8,0)
c
      mapl=ch_temp
c
      call av_ass(mapn,mapl)
c
      write(6,10)mapn
      if(record)write(nrec,10)mapn
10    format(' Update cell information for map:',i4)
c
      do j=1,6
        xrcell(j,mapn)=realval('AV_CELL:ENTER_CELL',0.0,100000.0,0)
      enddo
c
      call XRfrac(XRcell(1,mapn),XRtr(1,1,mapn),XRintr(1,1,mapn),
     &            XRvol,.true.)
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE cell_search
c     ======================
      IMPLICIT NONE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- local variables
c
      integer      map1,map2,env1,env2
      character*8  imap,omap,ienv,oenv
      logical      i_check,o_check,no_trans,use_out
c
      character*8  G2CHRU,ch_temp,c_string
      real         realval,Xrvol
      integer      nopt,keyword,j,intval,k
      integer      itype,nsym1,nsym2,ic_param(6),n_step
      logical      lass_get,r_check
      real         rot(3,3),c_step,c_start,c_end,ce_cs,c_store(3)
      real         c_now,r_chec
      character*13 mess(0:1)
      data mess/'on all pairs ','from the mean'/
c
c
      r_check  = .false.
      no_trans = .false.
      use_out  = .true.
      r_chec    = 1.0
c
c
      if(.not.lass_get(ch_temp,'MAPI'))then
c
        ch_temp=G2CHRU('AV_CSEARCH:INPUT_MAP',1,8,0)
c
        call lass_set(ch_temp,'MAPI')
      endif
      imap=ch_temp
c
c
      if(.not.lass_get(ch_temp,'MAPO'))then
c
        ch_temp=G2CHRU('AV_CSEARCH:OUTPUT_MAP',1,8,0)
c
        call lass_set(ch_temp,'MAPO')
      endif
      omap=ch_temp
c
c
      if(.not.lass_get(ch_temp,'ENVI'))then
c
        ch_temp=G2CHRU('AV_CSEARCH:INPUT_ENV',1,8,0)
c
        call lass_set(ch_temp,'ENVI')
      endif
      ienv=ch_temp
c
c
      if(.not.lass_get(ch_temp,'ENVO'))then
c
        ch_temp=G2CHRU('AV_CSEARCH:OUTPUT_ENV',1,8,0)
c
        call lass_set(ch_temp,'ENVO')
      endif
      oenv=ch_temp
c
c
      if(verbose)then
        write(6,110) imap,omap,ienv,oenv
        if(record)write(nrec,110) imap,omap,ienv,oenv
110     format(/' REFINE CELL PARAMETERS',/,
     &          ' ======================'
     &      /' Input map:       ',a,
     &      /' Output map:      ',a,
     &      /' Input envelope:  ',a,
     &      /' Output envelope: ',a)
      else
        write(6,111) imap,omap,ienv,oenv
        if(record)write(nrec,111) imap,omap,ienv,oenv
111     format(' REFINE CELL PARAMETERS',/,
     &    ' I/P map:',a,'  O/P map:',a,'  I/P env:',a,
     &    '  O/P env:',a)
      endif
C
      i_check=ienv.ne.'OFF'
      o_check=oenv.ne.'OFF'
c
      if(verbose)then
        if(i_check) write(6,120)
        if(.not.i_check)write(6,130)
        if(o_check) write(6,140)
        if(.not.o_check)write(6,150)
c
        if(record) then
          if(i_check) write(nrec,120)
          if(.not.i_check)write(nrec,130)
          if(o_check) write(nrec,140)
          if(.not.o_check)write(nrec,150)
        endif
      endif
c
120   format(' Input map envelope filtered')
130   format(' Input map not filtered')
140   format(' Output map envelope filtered')
150   format(' Output map not filtered')
c
c
c       get slot numbers for all maps
c
      call av_ass(map1,imap)
c----map 1 should exist by now !
      if(.not.defined(map1)) then
        write(6,160)
        if(record)write(nrec,160)
160     format(' %CELL_SEARCH-ERR: Input map empty')
        return
      endif
c
c---- check status of output file
c
      call av_ass(map2,omap)
c----map 2 should exist by now !
      if(.not.defined(map2)) then
        write(6,170)
        if(record) write(nrec,170)
170     format(' %CELL_SEARCH-ERR: Input map2 empty')
        return
      endif
c
c
c---- set up allocations for envelopes
c
      if(i_check) call av_ass(env1,ienv)
      if(o_check) call av_ass(env2,oenv)
c
c
c
c---- list symmetry operators to be used
c
      itype=11
      istat=1
      nsym1=1
      nsym2=nsym
c
10    nopt=keyword
     &  ('AV_CSEARCH_SYMM','ALL !STAR!END !GO  !INTP!CORQ!CORA',0)
c
      if(nopt.eq.2)nsym1=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.eq.3)nsym2=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.eq.5)itype=intval('AV_AVERAGE_INTP',1,64,0)
      if(nopt.eq.6)istat=1
      if(nopt.eq.7)istat=0
      if(nopt.ne.4)goto 10
c
      write(6,220)nsym1,nsym2
      if(record)write(nrec,220)nsym1,nsym2
220   format(' Will use symmetry operators from ',
     &             i4,' to',i4)
      if( (nsym1.gt.nsym2).or.(nsym2.gt.nsym) )then
        write(6,230)
        if(record)write(nrec,230)
230     format(
     &    ' %CELL_SEARCH-ERR: In choice of symmetry operators ')
        return
      endif
c
      if(verbose)then
        do j=nsym1,nsym2
          call prtsym(6,j)
        enddo
      endif
C----
c
      if(itype.ne.11.and.itype.ne.8.and.itype.ne.1.and.itype.ne.64)
     &    itype=11
      write(6,240)itype,mess(istat)
      if(record)write(nrec,240)itype,mess(istat)
240   format(i3,' point interpolation. Correlation coeff calculated ',a)
c
c
c---- set up map1 and map2. map2 will have the grid search on cell
c---- parameters done.  i would envisage map1 being the em cryo
c---- reconstructions.  map1s cell params will be altered and cc calced.
c---- at the mo refine a,b,c,ab or abc. no angle refinement.
c---- first ask for c_string to define what part parameter you want to refine
c
C---- intialise ic_param(j)
      do j=1,6
        ic_param(j)=0
      enddo
c       
250   c_string=G2CHRU('REFINE_CELL_PARAMETER:',1,4,0)
c
      if(c_string.eq.'A') ic_param(1)=1
      if(c_string.eq.'B') ic_param(2)=1
      if(c_string.eq.'C') ic_param(3)=1
      if(c_string.eq.'AB')then
        ic_param(1)=1
        ic_param(2)=1
      endif
      if(c_string.eq.'ABC')then
        ic_param(1)=1
        ic_param(2)=1
        ic_param(3)=1
      endif
      if((c_string.ne.'A').and.(c_string.ne.'B').and.(c_string.ne.'C')
     &    .and.(c_string.ne.'AB').and.(c_string.ne.'ABC')) goto 250
c
c---- save cell params for map2...for mo
c
      c_start = realval('AV_CSEARCH:C_MIN',0.0,100000.0,0)
      c_end   = realval('AV_CSEARCH:C_MAX',c_start,100000.0,0)
      ce_cs   = c_end-c_start
      c_step  = realval('AV_CSEARCH:C_STEP',0.0,ce_cs,0)
c
c
      c_store(1)=xrcell(1,map2)
      c_store(2)=xrcell(2,map2)
      c_store(3)=xrcell(3,map2)
c
c--- set up number of steps in grid search...
c
      n_step= nint( (ce_cs)/c_step)
c
      do j= 0,n_step
c
        c_now = c_start + float(j)*c_step
c
c---set up current cell parameter
c
        do k = 1,6
          if(ic_param(k).ne.0)xrcell(k,map2)=c_now
        enddo
c
c--- recalc frac to ortho matrices and inverse
c
        call i_stat
c
        call XRfrac(XRcell(1,map2),XRtr(1,1,map2),XRintr(1,1,map2),
     &              XRvol,.false.)
c
c
        call subaver(map1,map2,env1,env2,nsym1,nsym2,
     &    i_check,o_check,use_out,no_trans,ITYPE,rot(1,1),r_check)
c
        call g_stat
c
        if(r.lt.r_chec)then
          r_chec=r
c
          do k = 1,6
            if(ic_param(k).ne.0)c_store(k)=c_now
          enddo
c
          write(6,100),'For cell parameter size',c_now,' corr coeff=',
     &      cc,' r_fac=',r,' r_top=',r1,' **'
          if(record)write(nrec,100)
     &       ,'For cell parameter size ',c_now,' corr coeff=',cc,
     &        ' r_fac=',r,' r_top=',r1,' **'
100       format(a,f7.2,a,f6.3,a,f5.3,a,f5.3,a)
c
        else
          if(.not.terse)then
            write(6,100),'For cell parameter size',c_now,' corr coeff=',
     &                cc,' r_fac=',r,' r_top=',r1,' '
            if(record)write(nrec,100)
     &        ,'For cell parameter size ',c_now,' corr coeff=',cc,
     &         ' r_fac=',r,' r_top=',r1,' '
          endif
        endif
c
c--- finish n_steps
      enddo
c
c----  should we just use if(ic_param(1).ne. etc for 2 and 3 and
c
c--- also should we do gsearch on alpha beta gamma...if so trivial
c--- to add in
c
      do k = 1,6
        if(ic_param(k).ne.0)xrcell(k,map2)=c_store(k)
      enddo
c
      return
      end

C=============================================================================

      INTEGER FUNCTION COMP(SEG1, SEG2)
      IMPLICIT NONE

C Compare the rotation angles of two operators

      COMMON /ANGS/ ANG
C Modified subscript to match later definition, RME 22/4/2010
C     DOUBLE PRECISION ANG(1)
      INTEGER MAXNCS
      PARAMETER (MAXNCS=24)
      DOUBLE PRECISION ANG(MAXNCS)
      SAVE /ANGS/

      INTEGER SEG1, SEG2

      IF (ANG(SEG1) .LT. ANG(SEG2)) THEN
        COMP = -1
      ELSE IF (ANG(SEG1) .GT. ANG(SEG2)) THEN
        COMP = 1
      ELSE
        COMP = 0
      END IF

      RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE con_setup(map1,env1,i_check)
c     =======================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c----local variables:
      integer     map1, env1
      character*8 imap,ienv
c
      character*8 G2CHRU,ch_temp
      logical     lass_get,i_check
c
c      equivalence  (i_temp,ch_temp)
c
c
      if(.not.lass_get(ch_temp,'MAPI'))then
c
        ch_temp=G2CHRU('AV_CONVOLV:INPUT_MAP',1,8,0)
c
        call lass_set(ch_temp,'MAPI')
      endif
      imap=ch_temp
c
c
      if(.not.lass_get(ch_temp,'ENVI'))then
c
        ch_temp=G2CHRU('AV_CONVOLV:INPUT_ENV',1,8,0)
c
        call lass_set(ch_temp,'ENVI')
      endif
      ienv=ch_temp
c
c
      if(verbose)then
        write(6,110) imap,ienv
        if(record)write(nrec,110) imap,ienv
110     format(/' CONVOLUTE ELECTRON DENSITY',/,
     &          ' =========================='
     &      /' Input map:       ',a
     &      /' Input envelope:  ',a)
      else
        write(6,111) imap,ienv
        if(record)write(nrec,111) imap,ienv
111     format(' CONVOLUTE ELECTRON DENSITY',
     &    '  I/P map:',a'  I/P env:',a)
      endif
C
      i_check=ienv.ne.'OFF'
c
      if(verbose)then
        if(i_check) write(6,120)
        if(.not.i_check)write(6,130)
c
        if(record)then
          if(i_check) write(nrec,120)
          if(.not.i_check)write(nrec,130)
        endif
      endif
c
120   format(' Input map envelope filtered')
130   format(' Input map not filtered')
c
c
c       get slot numbers for all maps
c
      call av_ass(map1,imap)
c----map 1 should exist by now !
      if(.not.defined(map1)) then
        write(6,140)
        if(record)write(nrec,140)
140     format(' %CON_SETUP-ERR: Input map empty')
        return
      endif
c---- check status of output file
c
c---- set up allocations for envelopes
c
      if(i_check) call av_ass(env1,ienv)
c
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE convolv
c     ==================
      IMPLICIT NONE
c
c--- Perform convolution of map with a smearing function.
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- Local variables:
c
C RME: 23/4/2010: Made index to scratch space explicitly integer*8
      integer*8 j8
C
      integer     map1,env1,map2
      integer     iindex1(3),new_i(3)
      real        x11(3),x1(3),x22(3),x2(3)
      logical     i_check
      integer*2   irho
      integer*4   irsum
      integer     intval,nstep_now,n2,nstep_now_sq
      integer     ndiagadd,nxadd,nyadd,nzadd
      integer     iz,jz,iy,jy,ix,jx,jjz,jjy,jjx,jjz_z,jjy_y,jjx_x
      external    gp3,gp3_a
      integer*2   gp3,gp3_a
c
      character*8 ch_temp
      real        realval,B,pi,u2,u,fudge
      real        diag_weight,this_diag_weight,x_weight,this_x_weight
      real        z_weight,y_weight,this_y_weight,this_z_weight
      integer     nstep,nopt,keyword,np_err,edge,prev_edge
      integer*8   i8temp,offset
      logical     iiinout,iinout,check,no_symm
c
c--- Initialise
c
      no_symm=.false.
      check=.true.
      iiinout=.false.
      pi=4.0*atan(1.0)
c
c--- Defaults
c

      nstep=2
      B=50.0
      np_err=0
c
c--- Setup
c
      call con_setup(map1,env1,i_check)
c
10    nopt=keyword('AV_CONV','STEP!BFAC!GO  !BACK!NSYM',0)
c
      if(nopt.eq.1)nstep=intval('AV_CONV:NSTEP',1,30000,0)
      if(nopt.eq.2)B  =realval('AV_CONV:BFAC',-300000.0,300000.0,0)
      if(nopt.eq.4)return
      if(nopt.eq.5)NO_SYMM=.true.
      if(nopt.ne.3) goto 10
c
      if(no_symm)check=.false.
c
      if(.not.no_symm)then
        write(6,*)
     &    'Crystallographic symmetry used to "wrap-around" edges'
        if(record)write(nrec,*)
     &    'Crystallographic symmetry used to "wrap-around" edges'
      endif
c
      if(verbose)then
        if(no_symm)write(6,*)
     &    'Beware! No crystallographic symmetry used'
        write(6,100)nstep,b
        if(record.and.no_symm)
     &    write(nrec,*)
     &    'Beware! No crystallographic symmetry used'
        if(record)write(nrec,100)nstep,b
      else
        if(no_symm)write(6,*)
     &    'Beware! No crystallographic symmetry used'
        write(6,101)nstep,b
        if(record.and.no_symm)
     &    write(nrec,*)'Beware! No crystallographic symmetry used'
        if(record)write(nrec,101)nstep,b
      endif
100   format(' Largest step size of which to smear is:',i4,/,
     &       ' B factor to weight by is              :',f7.1)
101   format(' Largest step size:',i4,'   B factor:',f7.1)
c
c--- Report on Guassian
c--- This is based on Biso = 8(pi**2) * (<u>**2)/3
c
      u2=3.0*b/(8.0*pi**2)
      if(u2.ge.0)then
        u=sqrt(u2)
      else
        u=-sqrt(-u2)
      endif
      if(verbose)then
        write(6,110)u
        if(record)write(nrec,110)u
      endif
110   format(' Width of Gaussian:',f7.1)
c
c--- Report on function range
c
      if(.not.terse)then
        write(6,120)nstep
        if(record)write(nrec,120)nstep
120     format(' For each of the ',i4,' passes through the map',
     &    ' statistics will be printed')
      endif
c
c--- precalculate weights for our weighting scheme
c--- NB: For each nstep_now, this weight should be multiplied
c--- by nstep_now**2 to make the math correct 
c
c--- get map1 111 as xyz
c
      iindex1(1)=1
      iindex1(2)=1
      iindex1(3)=1
      call j3_to_f(iindex1,x11,map1)
      call f_to_o(x11,x1,map1)
c
c--- calculate weight on diagonals for nstep_now=1
c
      iindex1(1)=2
      iindex1(2)=2
      iindex1(3)=2
      call j3_to_f(iindex1,x22,map1)
      call f_to_o(x22,x2,map1)
      diag_weight=( (x2(1)-x1(1))**2 + (x2(2)-x1(2))**2
     &            + (x2(3)-x1(3))**2 ) / u2
c
c--- calculate weight for x neighbours when nstep_now=1
c
      iindex1(1)=2
      iindex1(2)=1
      iindex1(3)=1
      call j3_to_f(iindex1,x22,map1)
      call f_to_o(x22,x2,map1)
      x_weight=( (x2(1)-x1(1))**2 + (x2(2)-x1(2))**2
     &            + (x2(3)-x1(3))**2 ) / u2
c
c--- calculate weight for y neighbours when nstep_now=1
c
      iindex1(1)=1
      iindex1(2)=2
      iindex1(3)=1
      call j3_to_f(iindex1,x22,map1)
      call f_to_o(x22,x2,map1)
      y_weight=( (x2(1)-x1(1))**2 + (x2(2)-x1(2))**2
     &            + (x2(3)-x1(3))**2 ) / u2
c
c--- calculate weight for z neighbours when nstep_now=1
c
      iindex1(1)=1
      iindex1(2)=1
      iindex1(3)=2
      call j3_to_f(iindex1,x22,map1)
      call f_to_o(x22,x2,map1)
      z_weight=( (x2(1)-x1(1))**2 + (x2(2)-x1(2))**2
     &            + (x2(3)-x1(3))**2 ) / u2
c
c--- get our dummy map
c
      ch_temp="        "
      call av_ass(map2,ch_temp)
      call copy_head(map1,map2)
      call brickit(map2)
      call alloc(map2)
      offset=ns(map2)-ns(map1)
      do j8=ns(map2),ne(map2)
        scratch(j8)=scratch(j8-offset)
      enddo
c
c--- outer loop
c
      do nstep_now=1,nstep
        write(6,*)'  Step: ',nstep_now
        if(record)write(nrec,*)'   step: ',nstep_now
c
c---- we need double this step size for the loops later
c
        n2 = nstep_now*2
c
c---- Fudge factor to take account of the pixels we aren't bothering
c---- to check - we look at 14 out of however big this shell is.
c---- Assume evenly spread between x, y, z and diag pixels
c
        edge=n2 + 1
        prev_edge=n2 - 1
        fudge=float(edge**3 - prev_edge**3) / 14.0
c
c---- Precalculate our weighting scheme
c
        nstep_now_sq = nstep_now**2
        this_diag_weight = fudge*exp( -diag_weight * nstep_now_sq )
        this_x_weight    = fudge*exp( -x_weight * nstep_now_sq )
        this_y_weight    = fudge*exp( -y_weight * nstep_now_sq )
        this_z_weight    = fudge*exp( -z_weight * nstep_now_sq )
        write(6,*)'  Weights (diag,x,y,z):',
     &    this_diag_weight,this_x_weight,this_y_weight,this_z_weight
        if(record)write(nrec,*)'  Weights (diag,x,y,z):',
     &    this_diag_weight,this_x_weight,this_y_weight,this_z_weight
c
c---- Firstly set up pointers for the maps that drive the process
c---- (backwards!)
c
        call init_map(map1)
        if(i_check)call init_map(env1)
c
c---- Initialize r-factors, corr coefs, etc
c
        call f_stat
c
c---- Off to work...
c---- =============
c
c---- Loop over bricks
c
        do iz=1,n_brick(3,map1)
c
          jz=(iz-1)*brick(3,map1)
c
          do iy=1,n_brick(2,map1)
c
            jy=(iy-1)*brick(2,map1)
c
            do ix=1,n_brick(1,map1)
c
              jx=(ix-1)*brick(1,map1)
c
c---- Now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map1)
c
                if(jjz.le.nx(3,map1))then
c
                  do jjy=jy+1,jy+brick(2,map1)
c
                    if(jjy.le.nx(2,map1))then
c
                      do jjx=jx+1,jx+brick(1,map1)
c                      
                        if(jjx.le.nx(1,map1))then
c
      if(i_check)then
        if(gp3(env1,jjx,jjy,jjz).eq.0)goto 200
      endif
c
c----- Zero counters
c
      ndiagadd=0
      nxadd=0
      nyadd=0
      nzadd=0
      irsum=0
c
c----- Diag related pixels
c
      do jjz_z= jjz-nstep_now,jjz+nstep_now,n2
c
        do jjy_y= jjy-nstep_now,jjy+nstep_now,n2
c
          do jjx_x= jjx-nstep_now,jjx+nstep_now,n2
c
            iiinout=iinout(map1,jjx_x,jjy_y,jjz_z,new_i,check,1)
            if(.not.iiinout)then
              np_err=np_err+1
              if(np_err.lt.2)write(6,*)
     &          '%CONVOLV-ERR: Pixel can not be found'
            else
              irsum=irsum+(this_diag_weight*gp3_a(map1,new_i))
              ndiagadd=ndiagadd+1
            endif
c
          enddo
        enddo
      enddo
c
c----- X related pixels
c
      do jjx_x= jjx-nstep_now,jjx+nstep_now,n2
c
        iiinout=iinout(map1,jjx_x,jjy,jjz,new_i,check,1)
        if(.not.iiinout)then
          np_err=np_err+1
          if(np_err.lt.2)write(6,*)
     &      '%CONVOLV-ERR: Pixel can not be found'
        else
          irsum=irsum+(this_x_weight*gp3_a(map1,new_i))
          nxadd=nxadd+1
        endif
c
      enddo
c
c----- Y related pixels
c
      do jjy_y= jjy-nstep_now,jjy+nstep_now,n2
c
        iiinout=iinout(map1,jjx,jjy_y,jjz,new_i,check,1)
        if(.not.iiinout)then
          np_err=np_err+1
          if(np_err.lt.2)write(6,*)
     &      '%CONVOLV-ERR: Pixel can not be found'
        else
          irsum=irsum+(this_y_weight*gp3_a(map1,new_i))
          nyadd=nyadd+1
        endif
c
      enddo
c
c----- Z related pixels
c
      do jjz_z= jjz-nstep_now,jjz+nstep_now,n2
c
        iiinout=iinout(map1,jjx,jjy,jjz_z,new_i,check,1)
        if(.not.iiinout)then
          np_err=np_err+1
          if(np_err.lt.2)write(6,*)
     &      '%CONVOLV-ERR: Pixel can not be found'
        else
          irsum=irsum+(this_z_weight*gp3_a(map1,new_i))
          nzadd=nzadd+1
        endif
c
      enddo
c
      irho=gp3(map1,jjx,jjy,jjz)
      irho=nint( (irho+irsum) / ( 1.0 +
     &           (ndiagadd*this_diag_weight) +
     &           (nxadd*this_x_weight) +
     &           (nyadd*this_y_weight) +
     &           (nzadd*this_z_weight) ) )
c
c----- a_stat accumulates stats for this pass
c
      call f_a_stat(irho)
c
c----- Finally store the new pixel value in the output map
c
      call pp3(map2,jjx,jjy,jjz,irho)
c        
200   continue
                        endif
                      enddo
                    endif
                  enddo
                endif
              enddo
            enddo
          enddo
        enddo
c
c
c--- Print stats gathered during this cycle. NB This shows what went in, not
c--- what came out.
c
        call f_p_stat(.false.)
c
c--- Copy new map over original map
c
        i8temp=ns(map1)
        ns(map1)=ns(map2)
        ns(map2)=i8temp
c
        i8temp=ne(map1)
        ne(map1)=ne(map2)
        ne(map2)=i8temp
c
c--- And next step
c
      enddo
c
      if(np_err.gt.0)write(6,*)'%CONVOLV-ERR: ',
     &  np_err,'pixels outside map ignored'
c
c--- Delete temporary map - use parser deallocation routine
c
      if(.not.P2DAL(scratch(ns(map2)))) then
        write(6,*)' %DEL_MAP-ERROR: memory deallocation fails'
        if(record)write(nrec,*)
     &    ' %DEL_MAP-ERROR: memory deallocation fails'
      endif
c
c--- Reset first free pointer
c
      nmaps=nmaps-1
c
c--- And its all over
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE cop_map
c       ==================
      IMPLICIT NONE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      INCLUDE 'average.fcm'
c       =====================
c
c---- local variables:
c
        integer         map1,map2
        integer         iindex(3),new_index(3)
        real            aindex(3),x(3),x1(3)
        integer         ngadd,nout
        logical         iiinout,iinout,check
        integer*2       irho
        external        gp3,gp3_a
        integer*2       gp3,gp3_a
        integer         iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
c
c---- initialize logicals
c
        check=.true.
        iiinout=.false.
c
c---- get input and ouput map numbers
c
        call cop_set(map1,map2)
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
        call init_map(map2)
c
c---- initialize counters
c
        nout=0
        ngadd=0
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
        do iz=1,n_brick(3,map2)
c
          jz=(iz-1)*brick(3,map2)
c
          do iy=1,n_brick(2,map2)
c
            jy=(iy-1)*brick(2,map2)
c
            do ix=1,n_brick(1,map2)
c
              jx=(ix-1)*brick(1,map2)
c
c---- now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map2)
c
              if(jjz.le.nx(3,map2))then
c
                do jjy=jy+1,jy+brick(2,map2)
c
                if(jjy.le.nx(2,map2))then
c
                  do jjx=jx+1,jx+brick(1,map2)
c               
                  if(jjx.le.nx(1,map2))then
c
                    iindex(1)=jjx
                    iindex(2)=jjy
                    iindex(3)=jjz
c
c---- find equivalent position in map1 space
c
                    call j3_to_f(iindex,x1,map2)
                    call f_to_o(x1,x,map2)
                    call o_to_j3(x,aindex,iindex,map1)
c
c---- can we get to the volume map1 actually covers by using the CSOs?
c
                    iiinout=iinout(map1,iindex(1),iindex(2),iindex(3)
     &                      ,new_index,check,1)
                    if(iiinout)then
c
c---- yes, lets get on with it
c
                      irho=gp3_a(map1,new_index)
                      call pp3(map2,jjx,jjy,jjz,irho)
                      ngadd=ngadd+1
                    else
c
c---- not this time
c
                      nout=nout+1
                    endif
                  endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
        enddo
c
c---- phew, copied, lets see how we did ...
c
        if(verbose)then
          write(6,998)ngadd,nout
          if(record)write(nrec,998)ngadd,nout
        else
          write(6,997)ngadd,nout
          if(record)write(nrec,997)ngadd,nout
        endif
c
c---- take care of the piping
c
        if(pipe)then
          call pip('MAPO','MAPI',.true.)
          call pip('    ','MAPO',.false.)
        endif
c
c---- and its all over
c
        return
c
c---- format statements
c
997     format(' Copying finished,',i10,' pixels were copied,',
     &             i10,' rejected')
998     format(' COPYING FINISHED',/,
     &            ' ================',/,
     &             i10,' pixels were copied,',i10,' rejected')
c
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE cop_set(map1,map2)
c       =============================
      IMPLICIT NONE
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      INCLUDE 'average.fcm'
c       =====================
c
c--- local variables
c
        integer     map1,map2
        real        xrvol
        character*8 G2CHRU,ch_temp,imap,omap
        integer     nopt,keyword,j,intval,ics_local
        logical     lass_get,quit
c
c
c
        if(.not.lass_get(ch_temp,'MAPI'))then
          ch_temp=G2CHRU('AV_COPY_MAP:INPUT_MAP',1,8,0)
          call lass_set(ch_temp,'MAPI')
        endif
        imap=ch_temp
c
c
        if(.not.lass_get(ch_temp,'MAPO'))then
          ch_temp=G2CHRU('AV_COPY_MAP:OUTPUT_MAP',1,8,0)
          call lass_set(ch_temp,'MAPO')
        endif
        omap=ch_temp
c
c
        if(verbose)then
          write(6,110) imap,omap
          if(record)write(nrec,110) imap,omap
110       format(/' COPY ELECTRON DENSITY',/,
     &           ' ====================='
     &      /' Input map:       ',a,
     &      /' Output map:      ',a)
          else
          write(6,111) imap,omap
          if(record)write(nrec,111) imap,omap
111       format(' COPY ELECTRON DENSITY',
     &    '  I/P map:',a,'  O/P map:',a)
          endif
C
c
        call av_ass(map1,imap)
c----map 1 should exist by now !
        if(.not.defined(map1)) then
          write(6,160)
160       format('  %COP_SET-ERR: Input map empty')
          return
          endif
c---- check status of output file
        call av_ass(map2,omap)
        if(defined(map2)) then
          write(6,170)
          if(record)write(nrec,170)
170       format(' **WARNING** ',
     &' output map already defined: cannot extend')
        return
        endif
c
c---- set up output map header information
c
          call copy_head(map1,map2)
c
c---- check if we want to redefine the output grid (eg coarser than input)
c
20      nopt=keyword('AV_COPY_UPDATE','UPDA!RESE!BACK!QLIM!GO  ',0)
c
        if(nopt.eq.1) then
          do j=1,3
            nstart(j,map2)=intval('AV_AVERAGE:NSTART',-3000,3000,0)
            nend(j,map2)  =intval('AV_AVERAGE:NEND',-3000,3000,0)
            nunit(j,map2) =intval('AV_AVERAGE:NUNIT',-3000,3000,0)
            nx(j,map2)=nend(j,map2)-nstart(j,map2)+1
            norg(j,map2)=1-nstart(j,map2)
            if(verbose)then
              write(6,200)nstart(j,map2),nend(j,map2),
     &                  nx(j,map2),nunit(j,map2)
              if(record)
     &          write(nrec,200)nstart(j,map2),nend(j,map2),
     &                  nx(j,map2),nunit(j,map2)
200           format(' ',4i5)
              endif
          enddo
          npix(map2)=nx(1,map2)*nx(2,map2)*nx(3,map2)
          if(.not.terse)then
            write(6,210)npix(map2)
            if(record)write(nrec,210)npix(map2)
210         format(' There will be',i9,' pixels in the output map')
            endif
        endif
c
        if(nopt.eq.2)then
          do j=1,3
            iuvw(j,map2)=intval('AV_AVERAGE:ENTER_AXIS',1,3,0)
          enddo
        endif
c
        if(NOPT.eq.3)then
           quit=.true.
           return
           endif
c
        if(NOPT.eq.4)then
          write(6,214)(nstart(j,map2),nend(j,map2),nunit(j,map2),
     &      nx(j,map2),norg(j,map2),j=1,3)
214       format( 3(' START:',i4,' END:',i4,' GRID:',i4,
     &              ' POINTS:',i4,' ORIGIN:',i4,/) )
          endif
c
        if(nopt.ne.5)goto 20
c
c---- set up cell data
c
          call XRfrac(XRcell(1,map2),XRtr(1,1,map2),XRintr(1,1,map2),
     &    XRvol,.true.)
c
c---- set up symmetry data
c
          call XR_symsy(XRcell(1,map2),lgrp(map2),ics_local,
     &                  n_xsym(map2),XR_sym(1,map2) )
c
c---- and integerized CSO's:
c
          call XR_isym(map2)
c
c---- also set up pointers indicating which symm op to try first
c
           do j=1,maxsym
             last_sym(j,map2)=1
           enddo
c
c
c---- set up bricking data
c
          call brickit(map2)
c
c---- allocate space for output map
c
          call alloc(map2)
c
        return
        end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE copy(n,a,b)
c       ======================
      IMPLICIT NONE
        integer n,j
        real a(n),b(n)
        do j=1,n
        b(j)=a(j)
        enddo
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE copy_head(map1,map2)
c       ===============================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
        integer map1,map2
c
        call copyi(1,npix(map1),npix(map2))
        call copyi(1,npix_b(map1),npix_b(map2))
        call copyi(3,n_brick(1,map1),n_brick(1,map2))
        call copyi(3,brick(1,map1),brick(1,map2))
        call copyi(3,nend(1,map1),nend(1,map2))
        call copyi(3,nstart(1,map1),nstart(1,map2))
        call copyi(3,iuvw(1,map1),iuvw(1,map2))
        call copyi(3,nunit(1,map1),nunit(1,map2))
        call copyi(3,norg(1,map1),norg(1,map2))
        call copyi(3,nx(1,map1),nx(1,map2))
        call copy(9,XRtr(1,1,map1),XRtr(1,1,map2))
        call copy(9,XRintr(1,1,map1),XRintr(1,1,map2))
        call copy(9,XRcell(1,map1),XRcell(1,map2))
        call copyc(20,titles(map1),titles(map2))
        call copyi(1,mtype(map1),mtype(map2))
        call copyi(1,lgrp(map1),lgrp(map2))
        call copy(4,rholim(1,map1),rholim(1,map2))
        call copy(1,mscale(map1),mscale(map2))
        call copy(1,moffset(map1),moffset(map2))
        sollev(map2)=sollev(map1)
        if(.not.terse)then
          write(6,10)
          if(record)write(nrec,10)
10        format(' Headers copied from input to output map')
          endif
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE copyc(n,a,b)
c       ======================
      IMPLICIT NONE
        integer n,j
        character*1 a(n),b(n)
        do j=1,n
        b(j)=a(j)
        enddo
        return
        end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE copyi(n,a,b)
c       ======================
      IMPLICIT NONE
        integer n,a(n),b(n),j
        do j=1,n
        b(j)=a(j)
        enddo
        return
        end
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE d_copy(reel,duble,num)
      IMPLICIT NONE
        integer num,j
        double precision duble(num)
        real reel(num)
        do j=1,num
        duble(j)=dble(reel(j))
        enddo
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE del_map
c       ==================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c
        character*8 G2CHRU,ch_temp,label
        integer     j,n_del,k,k3,l3
c
c----   read label of map to be deleted, its too dangerous to use stream MAPI !
        ch_temp=G2CHRU('AV_DEL_MAP:ENTER_LABEL',1,8,0)
        label=ch_temp
c
c----   search through
        do j=1,nmaps
        n_del=j
        if (label.eq.labels(j)) goto 100
        enddo
        write(6,10)label
        if(record)write(nrec,10)label
10      format(' %DEL_MAP-ERROR: 'a,
     &    ' is not in the current list of map labels')
        return
c
c----  n_del is now the number of the slot to be deleted
100     continue
c
        write(6,110)label,n_del
        if(record)write(nrec,110)label,n_del
110     format(' MAP ',a,', NUMBER,'i4,' WILL BE DELETED FROM MEMORY')
c
c----   use parser deallocation routine
        if(.not.P2DAL(scratch(ns(n_del)))) then
         write(6,*)' %DEL_MAP-ERROR: memory deallocation fails'
         if(record)write(nrec,*)
     &     ' %DEL_MAP-ERROR: memory deallocation fails'
        endif
c
c----   pack things
c----   copy maps themselves
c----   increase free space counter
c----   set first free pointer
        nmaps=nmaps-1
        write(6,120)nmaps
        if(record)write(nrec,120)nmaps
120     format(i5,' maps remain')
c
c----   if no maps left, quit.
        if(nmaps.lt.1)return
c
c----   if map to be deleted is at end of list we are done !
        if(n_del.gt.nmaps)return
c
c----   now we must shuffle all the pointers, maps and control data along
c----   to fill the gap.
c
        do k=n_del,nmaps
c
c---- now compress all pointer arrays
          ns(k)     = ns(k+1)
          ne(k)     = ne(k+1)
          defined(k)= defined(k+1)
          npix(k)   = npix(k+1)
          npix_b(k) = npix_b(k+1)
          nx_tot(k) = nx_tot(k+1)
c
          labels(k) = labels(k+1)
          titles(k) = titles(k+1)
          mtype(k)  = mtype(k+1)
          mscale(k) = mscale(k+1)
          moffset(k)= moffset(k+1)
c
          pix_now(k) = pix_now(k+1)
          word_now(k)= word_now(k+1)
c
          lgrp(k)    = lgrp(k+1)
          n_xsym(k)  = n_xsym(k+1)
c
          do k3=1,3
            nx(k3,k)     = nx(k3,k+1)
            norg(k3,k)   = norg(k3,k+1)
            iuvw(k3,k)   = iuvw(k3,k+1)
            nunit(k3,k)  = nunit(k3,k+1)
            nstart(k3,k) = nstart(k3,k+1)
            nend(k3,k)   = nend(k3,k+1)
            n10(k3,k)    = n10(k3,k+1)
            nx_pad(k3,k) = nx_pad(k3,k+1)
            brick(k3,k)  = brick(k3,k+1)
            n_brick(k3,k)= n_brick(k3,k+1)
c
            do l3=1,3
              xrtr(k3,l3,k)  =xrtr(k3,l3,k+1)
              xrintr(k3,l3,k)=xrintr(k3,l3,k+1)
            enddo
c
          enddo
c
          do k3=1,4
            rholim(k3,k)=rholim(k3,k+1)
          enddo
c
          do k3=1,9
            xrcell(k3,k)=xrcell(k3,k+1)
          enddo
c
c
          do k3=1,maxxrc
            xr_sym(k3,k)  =xr_sym(k3,k+1)
            xr_insym(k3,k)=xr_insym(k3,k+1)
          enddo
c
          do k3=1,maxsym
            last_sym(k3,k)=last_sym(k3,k+1)
          enddo
c
          sollev(k)=sollev(k+1)
c
        enddo
c
        defined(nmaps+1)=.false.
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      FUNCTION DET(RWORK)
c     ===================
      IMPLICIT NONE
c
c*******************************************************************************
c
      real det,rwork
      DIMENSION RWORK(3,3)
      DET= RWORK(1,1)*( RWORK(2,2)*RWORK(3,3) - RWORK(2,3)*RWORK(3,2))
     &  -RWORK(1,2)*( RWORK(2,1)*RWORK(3,3) - RWORK(2,3)*RWORK(3,1) )
     &  +RWORK(1,3)*( RWORK(2,1)*RWORK(3,2) - RWORK(2,2)*RWORK(3,1) )
      RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE dir_map
c       ==================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        integer j
c
        if(nmaps.lt.1)then
          write(6,10)
          if(record)write(nrec,10)
10        format(' No maps in memory')
        else
          do j=1,nmaps
            write(6,30)j,labels(j)
            if(record)write(nrec,30)j,labels(j)
30          format(' Map number',i4,'  label: ',a)
          enddo
        endif
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE domain (n,r1,r2)
c       ===========================
      IMPLICIT NONE
c*******************************************************************************
c
c---- routine to manipulate matrices to get the difference in polar angles 
c
c
c*****************************************************************************
        common/av_rec /record,nrec
c*****************************************************************************
        logical record
        integer nrec,n
c
c
c       character*80 title
        real r1,r2,r3,r4,phi,psi,chi
        dimension r1(3,3),r2(3,3),r3(3,3),r4(3,3)
c
c
        call invrt(r1,r4,3)
c
        call matm33(r2,r4,r3)
c
C----
      CALL  POLAR  ( r3,PHI,PSI,CHI )
C----
      PHI= PHI*57.296
      PSI= PSI*57.296
      CHI= CHI*57.296
C----
      WRITE(6,400) n,CHI,PHI,PSI
      if(record)WRITE(nrec,400) n,CHI,PHI,PSI
400   FORMAT(' Refinement of NCS op',i3,' led to a ',f7.1,
     & ' deg rotation about an axis at phi=',F7.1,' psi=',F7.1)
c
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE exit_it
c       ==================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c
c
        character*20 date
        character*18 time
        integer      ilen
c
c
        call P2TIME(time,ilen)
        call P2DATE(date,ilen)
        write(6,10),time,date
        if(record)write(nrec,10),time,date
10      format(' PROGRAM FINISHES AT ',a,' ON ',a)
        stop
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE f_stat
c      =================
      IMPLICIT NONE
c
c--- These subroutines accumulate statistics on the map before flattening.
c--- Max, min, mean and sigma are calculated.
c
      INCLUDE 'average.fcm'
c      =====================
c
c--- Local variables
c
       integer   num_pixels
       integer*2 irho,irho_max,irho_min
       real*8    firho,sum_r,sum_sq,sav_tempf,ssign,ssig
       logical   full
c
       save num_pixels,sum_r,sum_sq,irho_max,irho_min
c
c--- Initialise counters
c
       av_temp=zero
       sigm=zero
       sum_r=zero
       sum_sq=zero
       num_pixels=0
       irho_max=-32000
       irho_min=32000
c
       return
c
c--- Accumulate statistics
c       
       entry f_a_stat(irho)
c      ====================
c
c--- Add an observation
c
       firho=float(irho)
c
       num_pixels=num_pixels+1
       sum_r=sum_r+firho
       sum_sq=sum_sq+(firho**2)
c
c--- Rho min/maxes
c
       if(irho.gt.irho_max)irho_max=irho
       if(irho.lt.irho_min)irho_min=irho
c
       return
c
c--- Print statistics
c
       entry f_p_stat(full)
c      ====================
c
c--- Sanity check - we should have some pixels
c
       if(num_pixels.le.0) then
         write(6,10)
         if(record)write(nrec,10)
         return
       endif
c
10     format(' %F_STAT-ERR: No pixels to analyse')
c
c--- Calculate statistics - av_temp, sigm
c--- NB av_temp and sigm are in one of the average.fcm common blocks
c
       av_temp=sum_r/float(num_pixels)
c
       sav_tempf=float(num_pixels)*(av_temp**2)
       ssign=sav_tempf + sum_sq - (2.0d0 * av_temp * sum_r)
       ssig = ssign/float(num_pixels)
       if(ssig.gt.0.0)sigm=sqrt(ssig)
c
c--- Print a summary of the statistics
c
       if(full.and.(.not.terse))write(6,20)
     &   irho_max,irho_min,av_temp,sigm
       if((.not.full).or.terse)  write(6,30)
     &   irho_max,irho_min,av_temp,sigm
c
       if(record) then
         if(full.and.(.not.terse))write(nrec,20)
     &     irho_max,irho_min,av_temp,sigm
         if((.not.full).or.terse)write(nrec,30)
     &     irho_max,irho_min,av_temp,sigm
       endif
c
20     format(' Statistics on map pixels used for current operation'/,
     &     ' Max rho',i6,
     &  '    min rho',i6,
     &  '    mean rho',f8.1,
     &  '    sigma',f8.1)
c
30     format(
     &     ' Max rho',i6,
     &  '    min rho',i6,
     &  '    mean rho',f8.1,
     &  '    sigma',f8.1)
c
       return
       end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE f_to_j3(x,aj3,j3,map)
c       ================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        real x(3),aj3(3)
        integer map,j3(3),j
c       
        do j=1,3
            aj3(j) = x(j)*float(nunit(j,map)) + float(norg(j,map))
            j3(j)  = nint(aj3(j))
        enddo
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE f_to_o(x,xo,map)
c       ===========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        real xo(3),x(3)
        integer map
c
        call vmtply(xrintr(1,1,map),x,xo)
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE flat_setup
c       =====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c----local variables:
        integer map1, env1
        character*8 imap,ienv
        integer     nopt,keyword,intval
        integer*2 i_solv
        external  gp3,gp3_a,same_grid
        integer*2 gp3,gp3_a
c
        character*8 G2CHRU,ch_temp
        real        realval,n_solv
        logical lass_get,i_check,r_solv
        logical     same,same_grid
c
c       equivalence  (i_temp,ch_temp)
c
c--- set flag so that map and mean solv level stays as it is calced...
c--- solvent level not set to i_solv
        r_solv=.false.
c
        if(.not.lass_get(ch_temp,'MAPI'))then
c
          ch_temp=G2CHRU('AV_FLATTEN:INPUT_MAP',1,8,0)
c
          call lass_set(ch_temp,'MAPI')
        endif
          imap=ch_temp
c
c
        if(.not.lass_get(ch_temp,'ENVI'))then
c
        ch_temp=G2CHRU('AV_FLATTEN:INPUT_ENV',1,84,0)
c
          call lass_set(ch_temp,'ENVI')
        endif
        ienv=ch_temp
c
c
        if(verbose)then
          write(6,110) imap,ienv
          if(record)write(nrec,110) imap,ienv
110       format(/' FLATTEN ELECTRON DENSITY',/,
     &           ' ========================'
     &      /' Input map:       ',a
     &      /' Input envelope:  ',a)
          else
          write(6,111) imap,ienv
          if(record)write(nrec,111) imap,ienv
111       format(' FLATTEN ELECTRON DENSITY',
     &    '  I/P map:',a,'  I/P env:',a)
          endif
C
        i_check=ienv.ne.'OFF'
c
c
        if(verbose)then
          if(i_check) write(6,120)
          if(.not.i_check)write(6,130)
c
          if(record)then
            if(i_check) write(nrec,120)
            if(.not.i_check)write(nrec,130)
            endif
          endif
c
120     format(' Input map envelope filtered')
130     format(' Input map not filtered')
c
c
c       get slot numbers for all maps
c
        call av_ass(map1,imap)
c----map 1 should exist by now !
        if(.not.defined(map1)) then
          write(6,140)
          if(record)write(nrec,140)
140       format(' %FLAT_SETUP-ERR: Input map empty')
          return
          endif
c---- check status of output file
c
c---- set up allocations for envelopes
c
        if(i_check) call av_ass(env1,ienv)
c
        if(i_check)then
        same=same_grid(map1,env1)
c
        if(.not.same)then
                write(6,49)
                 if(record)write(nrec,49)
49      format(' %FLAT_SETUP-ERR: Mapi and envi not on same grid')
                goto 700
        endif
        endif
c
c---- option to set solvent level to 0.  with float one sets i_solv
c---- to be -30000 to 30000....
c
600     nopt=keyword('AV_FLAT_FLOAT:','FLOT!FLOS!GO  !BACK',0)
c
        if(nopt.eq.1)then
          r_solv=.true.
          i_solv=intval('AV_FLAT_FLOAT:ENTER_SOLV_LEVEL',-32000,32000,0)
        endif
c
        if(nopt.eq.2)then
          r_solv=.true.
          n_solv=realval('AV_FLAT_FLOAT:ENTER_SOLV_LEVEL',-100000.0,
     &      100000.0,0)
          i_solv=nint(n_solv*mscale(map1))
        endif
c
        if(nopt.eq.4) goto 700
        if(nopt.ne.3) goto 600


        call flatten(map1,env1,i_check,r_solv,i_solv)
c
700     return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE flatten(map1,env1,i_check,r_solv,i_solv)
c       ===================================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c----local variables:
        integer map1,  env1
        logical i_check,r_solv
        integer*2 irho,av_dens,i_solv
        external  gp3,gp3_a
        integer*2 gp3,gp3_a
c
        integer   ngadd
        integer iz,jz,iy,jy,ix,jx,jjz,jjy,jjx

c
c       equivalence  (i_temp,ch_temp)
c
        if(.not.terse)then
          if(.not.i_check)then
            write(6,*)'**WARNING** input map in peril'
            if(record)write(nrec,*)'**WARNING** input map in peril'
            endif
          endif
c
c
c
c
c---- the averaging process is very simple, we are driven by the
c---- output map, for each pixel in that map (accessed in a sequential
c---- fashion, but since map is bricked we move slowly through space)
c---- flatten sets all pixels above and below certain values to zero.
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
        call init_map(map1)
c
c---- initialize r-factors, corr coefs, etc
c
        ngadd=0
c
        call f_stat
c
        if(i_check)call init_map(env1)
c       if(i_check)write(6,*) ' i_check true '
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
        do iz=1,n_brick(3,map1)
c
          jz=(iz-1)*brick(3,map1)
c
          do iy=1,n_brick(2,map1)
c
            jy=(iy-1)*brick(2,map1)
c
            do ix=1,n_brick(1,map1)
c
              jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                if(jjy.le.nx(2,map1))then
c
                  do jjx=jx+1,jx+brick(1,map1)
c               
                  if(jjx.le.nx(1,map1))then
c
c---- check that this pixel is in required region of the output map
c
                    if(i_check)then
                      if(gp3(env1,jjx,jjy,jjz).eq.0)then
c                       write(6,*) ' get solv irho '
                        irho=gp3(map1,jjx,jjy,jjz)
                        ngadd=ngadd+1
                        call f_a_stat(irho)
                      endif
                    else
                      irho=gp3(map1,jjx,jjy,jjz)
                      ngadd=ngadd+1
c                       write(6,*) ' get irho '
                      call f_a_stat(irho)
                    endif
c
                  endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
        enddo
c
c---- print some stats
c
        call f_p_stat(.true.)
c
c-----we put in if r_solv=false do below blah blah
c
        av_dens=nint(av_temp)
c
        if(.not.r_solv)then
c
c---- make header
c
        sollev(map1)=av_dens
c
        write(6,100)av_dens
        if(record)write(nrec,100)av_dens
100     format(' Solvent region will be set to:',i5)
c
c----  now work through the map setting all pixel points below av_dens
c----    to av_dens....flatten map to average density.
c
        call init_map(map1)
c
        if(i_check)call init_map(env1)
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
        do iz=1,n_brick(3,map1)
c
          jz=(iz-1)*brick(3,map1)
c
          do iy=1,n_brick(2,map1)
c
            jy=(iy-1)*brick(2,map1)
c
            do ix=1,n_brick(1,map1)
c
              jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                if(jjy.le.nx(2,map1))then
c
                  do jjx=jx+1,jx+brick(1,map1)
c               
                  if(jjx.le.nx(1,map1))then
c
c---- check that this pixel is in required region of the output map
c
c
                    if(i_check)then
                      if(gp3(env1,jjx,jjy,jjz).eq.0)call
     &                       pp3(map1,jjx,jjy,jjz,av_dens)
                    else
                      call pp3(map1,jjx,jjy,jjz,av_dens)
                    endif
c
                  endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
        enddo
c
        else
c
c----- make header
c
        sollev(map1)=i_solv
c
        write(6,101)i_solv
        if(record)write(nrec,101)i_solv
101     format(' Solvent region to be set to',i5
     &  ,' and map floated to this level.')
c
c----  now work through the map setting all pixel points below av_dens
c----    to av_dens....flatten map to average density.
c
        call init_map(map1)
c
        if(i_check)call init_map(env1)
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
        do iz=1,n_brick(3,map1)
c
          jz=(iz-1)*brick(3,map1)
c
          do iy=1,n_brick(2,map1)
c
            jy=(iy-1)*brick(2,map1)
c
            do ix=1,n_brick(1,map1)
c
              jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                if(jjy.le.nx(2,map1))then
c
                  do jjx=jx+1,jx+brick(1,map1)
c               
                  if(jjx.le.nx(1,map1))then
c
c---- check that this pixel is in required region of the output map
c
c
                    if(i_check)then
                      if(gp3(env1,jjx,jjy,jjz).eq.0)call
     &                       pp3(map1,jjx,jjy,jjz,av_dens)
                    else
                      call pp3(map1,jjx,jjy,jjz,av_dens)
                    endif
c
c--- now set level to irho - av_dens + float
c
                irho=gp3(map1,jjx,jjy,jjz)
c
                irho=irho - av_dens + i_solv
c
                call pp3(map1,jjx,jjy,jjz,irho)
c               call f_a_stat(irho)
c
                  endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
        enddo
c
c
        endif
c
c999     continue
        if(verbose)then
          write(6,998)ngadd
          if(record)write(nrec,998)ngadd
998       format(/,' FLATTENING FINISHED,',/,
     &            ' ===================',/,
     &             i9,' pixels were inspected')
c
          if(.not.i_check)then
            write(6,*)'(that was the whole map)'
            if(record)write(nrec,*)'(that was the whole map)'
            endif
          else
          write(6,997)ngadd
          if(record)write(nrec,997)ngadd
997       format(' FLATTENING FINISHED',
     &    i10,' pixels were inspected')
          endif
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE flip_map(map1,env1,i_check)
c       ======================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c--- local variables
c
        integer   map1,env1
        logical   i_check
c
        integer   ngadd,iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
        integer*2 irho,av_dens,i_solv
        external  gp3,gp3_a
        integer*2 gp3,gp3_a
c
c
        if(.not.terse)then
          if(.not.i_check)then
            write(6,*)'**WARNING** input map in peril'
            if(record)write(nrec,*)'**WARNING** input map in peril'
            endif
          endif
c
c---- the averaging process is very simple, we are driven by the
c---- output map, for each pixel in that map (accessed in a sequential
c---- fashion, but since map is bricked we move slowly through space)
c---- flatten sets all pixels above and below certain values to zero.
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
        call init_map(map1)
c
c---- initialize r-factors, corr coefs, etc
c
        ngadd=0
c
        call f_stat
c
        if(i_check)call init_map(env1)
c       if(i_check)write(6,*) ' i_check true '
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
        do iz=1,n_brick(3,map1)
c
          jz=(iz-1)*brick(3,map1)
c
          do iy=1,n_brick(2,map1)
c
            jy=(iy-1)*brick(2,map1)
c
            do ix=1,n_brick(1,map1)
c
              jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                if(jjy.le.nx(2,map1))then
c
                  do jjx=jx+1,jx+brick(1,map1)
c               
                  if(jjx.le.nx(1,map1))then
c
c---- check that this pixel is in required region of the output map
c
                    if(i_check)then
                      if(gp3(env1,jjx,jjy,jjz).eq.0)then
c
                         irho=gp3(map1,jjx,jjy,jjz)
c
                         ngadd=ngadd+1
                         call f_a_stat(irho)
                      endif
                    else
                        irho=gp3(map1,jjx,jjy,jjz)
                        ngadd=ngadd+1
c                       write(6,*) ' get irho '
                        call f_a_stat(irho)
                    endif
c
                  endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
        enddo
c
c---- print some stats
c
        call f_p_stat(.true.)
c
c-----we put in if r_solv=false do below blah blah
c
        av_dens= 2 * (nint(av_temp))
c
c       write(6,*) '2 * avdens ', av_dens
c
        write(6,100)
        if(record)write(nrec,100)
100     format(' Solvent region will be flipped around mean value ')
c
c
        call init_map(map1)
c
c---- back through map to flip solvent density values around mean
c
        if(i_check)call init_map(env1)
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
        do iz=1,n_brick(3,map1)
c
          jz=(iz-1)*brick(3,map1)
c
          do iy=1,n_brick(2,map1)
c
            jy=(iy-1)*brick(2,map1)
c
            do ix=1,n_brick(1,map1)
c
              jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                if(jjy.le.nx(2,map1))then
c
                  do jjx=jx+1,jx+brick(1,map1)
c               
                  if(jjx.le.nx(1,map1))then
c
c---- check that this pixel is in required region of the output map
c
                    if(i_check)then
                      if(gp3(env1,jjx,jjy,jjz).eq.0)then
c
                         irho=gp3(map1,jjx,jjy,jjz)
c
c--- flip irho about mean....ie (2*av_dens - irho)
c
                         i_solv = av_dens - irho
                        call pp3(map1,jjx,jjy,jjz,i_solv)
c
                      endif
                    else
                         irho=gp3(map1,jjx,jjy,jjz)
                         i_solv = av_dens - irho
                         call pp3(map1,jjx,jjy,jjz,i_solv)
                    endif
c
                  endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
        enddo
c
c
c
        if(verbose)then
          write(6,998)ngadd
          if(record)write(nrec,998)ngadd
998       format(/,' FLIPPING FINISHED,',/,
     &            ' =================',/,
     &             i9,' pixels were inspected')
c
          if(.not.i_check)then
            write(6,*)'(that was the whole map)'
            if(record)write(nrec,*)'(that was the whole map)'
            endif
          else
          write(6,997)ngadd
          if(record)write(nrec,997)ngadd
997       format(' FLIPPING FINISHED',
     &    i10,' pixels were inspected')
          endif
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE flip_setup
c       =====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c--- local variables
c
        integer     map1,env1
        character*8 G2CHRU,ch_temp,imap,ienv
        logical     lass_get,i_check,same,same_grid
c
c--- set flag so that map and mean solv level stays as it is calced...
c--- solvent level not set to i_solv
c
        if(.not.lass_get(ch_temp,'MAPI'))then
c
          ch_temp=G2CHRU('AV_FLIP:INPUT_MAP',1,8,0)
c
          call lass_set(ch_temp,'MAPI')
        endif
          imap=ch_temp
c
c
        if(.not.lass_get(ch_temp,'ENVI'))then
c
        ch_temp=G2CHRU('AV_FLIP:INPUT_ENV',1,84,0)
c
          call lass_set(ch_temp,'ENVI')
        endif
        ienv=ch_temp
c
c
        if(verbose)then
          write(6,110) imap,ienv
          if(record)write(nrec,110) imap,ienv
110       format(/' FLIP ELECTRON DENSITY',/,
     &           ' ====================='
     &      /' Input map:       ',a
     &      /' Input envelope:  ',a)
          else
          write(6,111) imap,ienv
          if(record)write(nrec,111) imap,ienv
111       format(' FLIP ELECTRON DENSITY',
     &    '  I/P map:',a,'  I/P env:',a)
          endif
C
        i_check=ienv.ne.'OFF'
c
c
        if(verbose)then
          if(i_check) write(6,120)
          if(.not.i_check)write(6,130)
c
          if(record)then
            if(i_check) write(nrec,120)
            if(.not.i_check)write(nrec,130)
            endif
          endif
c
120     format(' Input map envelope filtered')
130     format(' Input map not filtered')
c
c
c       get slot numbers for all maps
c
        call av_ass(map1,imap)
c----map 1 should exist by now !
        if(.not.defined(map1)) then
          write(6,140)
          if(record)write(nrec,140)
140       format(' %FLIP_SETUP-ERR: Input map empty')
          return
          endif
c---- check status of output file
c
c---- set up allocations for envelopes
c
        if(i_check) call av_ass(env1,ienv)
c
        if(i_check)then
        same=same_grid(map1,env1)
c
        if(.not.same)then
                write(6,49)
                 if(record)write(nrec,49)
49      format(' %FLIP_SETUP-ERR: Mapi and envi not on same grid')
                goto 700
        endif
        endif
c

        call flip_map(map1,env1,i_check)
c
700     return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE flood(map, inew, ostart, ispace)
c     ===========================================
      IMPLICIT NONE
c
c Always called from flood_setup
c
c This is a simple flood-fill algorithm:
c
c  Allocate ispace*3 space for a stack of indices
c  Use the passed orthogonal coordinates to get the first index
c  Get the test value from the first point
c  Set the first point to the passed new value
c  Push its index to the stack
c
c  Pull an index off the stack
c  Get values from its 6 actual neighbours (not the diagonals)
c  If value matches test, set it to the new value and push its
c    coordinates to the stack
c  Loop until there is nothing left on the stack
c
c Should make some effort to follow xtal symmetry when desired
c
c Modified (read blatantly stolen) from Robert Esnouf's Volumes Program
c JMD Jul2000
c
      INCLUDE 'average.fcm'
c
      integer    map,ispace
      integer*2  inew
      real       ostart(3)
c
      integer    ipixels_set,stack(2),iindex(3),jx,jy,jz,joffset,isym
      integer*8  bottom_of_stack,top_of_stack
      integer*2  itest,gp3_a
      real       aindex(3)
      logical    i_box,spin_in
c
c--- Allocate stack space
c
      ispace=ispace*3
      bottom_of_stack=p2al(stack(1),stack(2),ispace)
      if(bottom_of_stack.eq.0)then
        write(6,*)'AV_FLOOD-ERR% Allocate failed for stack'
        if(record)write(6,*)'AV_FLOOD-ERR% Allocate failed for stack'
      endif
c
c--- Find the starting gridpoint
c
      call o_to_j3(ostart,aindex,iindex,map)
c
c--- Get the test value from the starting gridpoint
c
      itest=gp3_a(map,iindex)
      write(6,*)'Test value = ',itest
      if(record)write(nrec,*)'Test value = ',itest
c
c--- If the test value equals the new value we're in the poo, so run away
c
      if(itest.eq.inew)then
        write(6,*)'AV_FLOOD-ERR% Test value equals fill value',
     &   ' - stopping'
        if(record)write(nrec,*)'AV_FLOOD-ERR% Test value equals',
     &   ' fill value - stopping'
        return
      endif
c
c--- Set it to the new value
c
      call pp3_a(map,iindex,inew)
      ipixels_set=1
c
c--- Push it to the stack
c
      top_of_stack=bottom_of_stack+3
      stack(top_of_stack-2)=iindex(1)
      stack(top_of_stack-1)=iindex(2)
      stack(top_of_stack)=iindex(3)
c
c--- Now the main loop....
c
c--- Pop the next pixel off the stack and check its neighbours
c
100   jx=stack(top_of_stack-2)
      jy=stack(top_of_stack-1)
      jz=stack(top_of_stack)
      top_of_stack=top_of_stack-3
c
c--- Check the 6 neighbours
c
      do joffset=1,6
c
        iindex(1)=jx
        iindex(2)=jy
        iindex(3)=jz
c
        if(joffset.eq.1)then
          iindex(1)=iindex(1)-1
        elseif(joffset.eq.2)then
          iindex(1)=iindex(1)+1
        elseif(joffset.eq.3)then
          iindex(2)=iindex(2)-1
        elseif(joffset.eq.4)then
          iindex(2)=iindex(2)+1
        elseif(joffset.eq.5)then
          iindex(3)=iindex(3)-1
        else
          iindex(3)=iindex(3)+1
        endif
c
c---- If this neighbour is not in the map
c
        if(.not.i_box(iindex,map))then
c
c----- If we are using xtal symmetry, attempt to map it
c
          if(xr_con(map))then
            do isym=1,n_xsym(map)
              if(spin_in(iindex,map,isym))goto 200
            enddo
            goto 300
c
c----- Else skip to next neighbour
c
          else
            goto 300
          endif
        endif
c
c---- If this neighbour matches our target...
c
200     continue
        if(gp3_a(map,iindex).eq.itest)then
c
c----- Keep count
c
          ipixels_set=ipixels_set+1
c
c----- Set it to the new value
c
          call pp3_a(map,iindex,inew)
c
c----- Check to see if the stack is full
c
          if((top_of_stack-bottom_of_stack+3).gt.ispace)then
            write(6,*)'AV_FLOOD-ERR% Stack full'
            if(record)write(nrec,*)'AV_FLOOD-ERR% Stack full'
            return
          endif
c
c----- Push the coords to the stack
c
          top_of_stack=top_of_stack+3
          stack(top_of_stack-2)=iindex(1)
          stack(top_of_stack-1)=iindex(2)
          stack(top_of_stack)=iindex(3)
c
        endif
c
c---- Next neighbour
c
300     continue
      enddo
c
c--- Next pixel off stack
c
      if(top_of_stack.gt.bottom_of_stack)goto 100
c
c--- Write some stats
c
      write(6,*)ipixels_set,' pixels set to ',inew
      if(record)write(nrec,*)ipixels_set,' pixels set to ',inew
c
c--- Should be done by now, so deallocate stack space
c
      if(.not.p2dal(stack(bottom_of_stack)))then
        write(6,*)'AV_FLOOD_ERR% Memory deallocate fails'
        if(record)write(nrec,*)'AV_FLOOD_ERR% Memory deallocate fails'
      endif
c
c--- And we're done...
c
      return
      end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE flood_setup
c     =====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c Flood fill (usually an envelope). See subroutine flood for details.
c
c--- Local variables:
c
      integer     map
      character*8 imap,g2chru
      integer     nopt,keyword,intval,ispace,j
      integer*2   inew
      real        realval,ostart(3)
      logical     lass_get
c
c--- Get input map if not already assigned
c
      if(.not.lass_get(imap,'MAPI'))then
        imap=G2CHRU('AV_FLOOD:INPUT_MAP',1,8,0)
        call lass_set(imap,'MAPI')
      endif
c
c--- Report
c
      if(verbose)then
        write(6,100)imap
        if(record)write(nrec,100)imap
100     format(/' FLOOD-FILL',/,
     &          ' =========='
     &         /' Input map:       ',a)
      else
        write(6,200)imap
        if(record)write(nrec,200)imap
200     format(' FLOOD-FILL  I/P map:',a)
      endif
c
c--- Get slot number for input map
c
      call av_ass(map,imap)
c
c--- The input map should exist by now !
c
      if(.not.defined(map)) then
        write(6,300)
        if(record)write(nrec,300)
300     format(' %FLOOD_SETUP-ERR: Input map empty')
        return
      endif
c
c--- Keyworded input
c---
c--- FILL - value to set filled points to
c--- FROM - initial xyz in Angstroms for fill
c--- NPIX - way of increasing stack space for big envs
c
c--- Defaults
c
      inew=100
      ostart(1)=1.0
      ostart(2)=1.0
      ostart(3)=1.0
      ispace=10000
c
400   nopt=keyword('AV_FLOOD:','FILL!FROM!NPIX!GO  !BACK',0)
c
c--- FILL
c
      if(nopt.eq.1)then
        inew=intval('AV_FLOOD:ENTER_FILL_VALUE',-32000,32000,0)
      endif
c
c--- FROM
c
      if(nopt.eq.2)then
        do j=1,3
          ostart(j)=realval('AV_FLOOD:ENTER_INIT_COORD',
     &            -100000.0,100000.0,0)
        enddo
      endif
c
c--- NPIX
c
      if(nopt.eq.3)then
        ispace=intval('AV_FLOOD:ENTER_NPIX',10000,715800000,0)
      endif
c
c--- BACK
c
      if(nopt.eq.5) goto 500
c
c--- GO
c
      if(nopt.ne.4) goto 400
c
c--- Do it!
c
      call flood(map,inew,ostart,ispace)
c
500   return
      end
c
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE fold_av(map1,map2,env1,env2,i_check,o_check)
c     =======================================================
      IMPLICIT NONE
c
c--- Fold map1 into map2, where map2 is a crystallographic asymmetric unit
c---
c---  map1    =  input map number
c---  map2    =  output map number
c---  env1    =  input envelope number
c---  env2    =  output envelope number
c---  i_check =  use input envelope?
c---  o_check =  use output envelope?
c
c--- Global variables
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- Passed variables:
c
      integer   map1,map2,env1,env2
      logical   i_check,o_check
c
c--- Local variables
c
      integer   index1(3),new_index(3)
      integer   index2(3),newer_index(3),temp_index(3)
      real      x(3),x1(3)
      real      aindex2(3)
      integer   j
      integer*2 irho,iienv
      external  gp3,gp3_a
      integer*2 gp3,gp3_a
      integer   ngadd,j_x,j_y,j_z,j_out
      integer   iz,jz,iy,jy,ix,jx,jjz,jjy,jjx,n_pix_store
c
c--- Defaults
c
      iienv=1000
c
c--- The folding process is very simple: we are driven by the
c--- input map. For each pixel in that map (accessed in a sequential
c--- fashion, but since map is bricked we move slowly through space)
c--- we use the CSOs to try and map that pixel into the output map
c
c--- Set up pointers for the map that drives the process
c
      call init_map(map1)
c
c--- Initialize R-factors, Corr coefs, etc
c
      ngadd=0
      n_pix_store=0
      call i_stat
c
c--- If we're using an input envelope we also need to set pointers
c
      if(i_check)call init_map(env1)
c
c
c--- Off to work...
c--- ==============
c
c--- loop over bricks
c
      do iz=1,n_brick(3,map1)
c
       jz=(iz-1)*brick(3,map1)
c
       do iy=1,n_brick(2,map1)
c
        jy=(iy-1)*brick(2,map1)
c
        do ix=1,n_brick(1,map1)
c
         jx=(ix-1)*brick(1,map1)
c
c---- Now loop within each brick
c
         do jjz=jz+1,jz+brick(3,map1)
c
          if(jjz.le.nx(3,map1))then
c
           do jjy=jy+1,jy+brick(2,map1)
c
            if(jjy.le.nx(2,map1))then
c
             do jjx=jx+1,jx+brick(1,map1)
c               
              if(jjx.le.nx(1,map1))then
c
c----- Check that this pixel is in a permitted region of the input map
c
               if(i_check)then
                if(gp3(env1,jjx,jjy,jjz).eq.0)goto 200
               endif
c
c----- Put index in an array
c
               index1(1)=jjx
               index1(2)=jjy
               index1(3)=jjz
c
c----- Get the orthogonalised co-ordinates of this pixel
c

               call j3_to_f(index1,x1,map1)
               call f_to_o(x1,x,map1)
c
c----- Get the co-ordinates of this pixel in map2 space
c----- Why? If they are on the same grid , its the same space?
c
               call o_to_j3(x,aindex2,index2,map2)
c
c----- Get density from map1
c
               irho=gp3_a(map1,index1)
c
c----- Loop over the crystallographic symmetry operators
c
               do j=1,n_xsym(map2)
c
c------ Apply this CSO
c
                call spin_out(map2,index2,new_index,j)
c
c------ Set initial x index ... down two unit cell edges
c
                newer_index(1)=new_index(1)-2*nunit(1,map2)
c
c------ Loop over unit cells adjacent in x
c
                do j_x= -1,1
c
c------- Move up one x
c
                 newer_index(1)=newer_index(1)+nunit(1,map2)
c
c------- Check whether within output map
c
                 if((newer_index(1).ge.nstart(1,map2)).and.
     &                 (newer_index(1).le.nend(1,map2)))then
c
c-------- Set initial y index
c
                  newer_index(2)=new_index(2)-2*nunit(2,map2)
c
c-------- Loop over unit cells adjacent in y
c
                  do j_y= -1,1
c
c--------- Move up one y
c
                   newer_index(2)=newer_index(2)+nunit(2,map2)
c
c--------- Check whether within output map
c
                   if((newer_index(2).ge.nstart(2,map2)).and.
     &                   (newer_index(2).le.nend(2,map2)))then
c
c---------- Set initial z index
c
                    newer_index(3)=new_index(3)-2*nunit(3,map2)
c
c---------- Loop over unit cells adjacent in z
c
                    do j_z= -1,1
c
c----------- Move up one z
c
                     newer_index(3)=newer_index(3)+nunit(3,map2)
c
c----------- Check whether within output map
c
                     if((newer_index(3).ge.nstart(3,map2)).and.
     &                      (newer_index(3).le.nend(3,map2)))then
c
c------------ In-frame, so revert to local indices
c------------ I don't really see what's going on here...
c
                      do j_out=1,3
                       temp_index(j_out)=newer_index(j_out)
     &                                      +norg(j_out,map2)
                      enddo
c
c------------ Check that new index is within the output envelope
c------------ iienv is never changed from 1000 if no output env
c
                      if(o_check)iienv=gp3_a(env2,temp_index)
c
c------------ If OK, put the pixel from map1 here in map2
c
                      if(iienv.ne.0)then
                       call pp3_a(map2,temp_index,irho)
                       n_pix_store=n_pix_store+1
                      endif
c
c------------ Next ouput map pixel please
c
                     endif
                    enddo
                   endif
                  enddo
                 endif
                enddo
               enddo
c
c------ Keep track of what we have done
c
               ngadd=ngadd+1
c
c------ And next input map pixel please
c
200            continue
              endif
             enddo
            endif
           enddo
          endif
         enddo
        enddo
       enddo
      enddo
c
c--- Report on what happened
c
c999   continue
      if(verbose)then
        write(6,998)ngadd,n_pix_store
        if(record)write(nrec,998)ngadd,n_pix_store
998     format(/,' PACKING DOWN FINISHED',/,
     &           ' =====================',/,
     &             i9,' pixels were filtered from i/p map',/,
     &             i9,' pixels filtered through to o/p map')
      else
        write(6,997)ngadd,n_pix_store
        if(record)write(nrec,997)ngadd,n_pix_store
997     format(' PACKING DOWN FINISHED',/,
     &    i10,' pixels were taken from i/p map',i10,
     &    ' pixels filtered through to o/p map')
      endif
c
c--- Pipe MAPO into MAPI
c
      if(pipe)then
        call pip('MAPO','MAPI',.true.)
        call pip('    ','MAPO',.false.)
      endif
c
c--- Bye!
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      SUBROUTINE fold_set
c     ===================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c     =====================
c
c----local variables:
c
      integer     map1,map2,env1,env2,nopt,keyword,j,intval,ics_local
      real        xrvol
      character*8 G2CHRU,imap,omap,ienv,oenv
      logical     i_check,o_check,quit,same,same_grid,lass_get
c
c--- Check if MAPI has an assignment - if not then get one
c--- Afterwards imap contains the input map name
c
      if(.not.lass_get(imap,'MAPI'))then
        imap=G2CHRU('AV_FOLD:INPUT_MAP',1,8,0)
        call lass_set(imap,'MAPI')
      endif
c
c--- Check if MAPO has an assignment - if not then get one
c--- Afterwards omap contains the output map name
c
      if(.not.lass_get(omap,'MAPO'))then
        omap=G2CHRU('AV_FOLD:OUTPUT_MAP',1,8,0)
        call lass_set(omap,'MAPO')
      endif
c
c--- Check if ENVI has an assignment - if not then get one
c--- Afterwards ienv contains the input envelope name
c
      if(.not.lass_get(ienv,'ENVI'))then
        ienv=G2CHRU('AV_FOLD:INPUT_ENV',1,8,0)
        call lass_set(ienv,'ENVI')
      endif
c
c--- Check if ENVO has an assignment - if not then get one
c--- Afterwards envo contains the output envelope name
c
      if(.not.lass_get(oenv,'ENVO'))then
        oenv=G2CHRU('AV_FOLD:OUTPUT_ENV',1,8,0)
        call lass_set(oenv,'ENVO')
      endif
c
c--- Report on the state of play
c
      if(verbose)then
        write(6,110) imap,omap,ienv,oenv
        if(record)write(nrec,110) imap,omap,ienv,oenv
110     format(/' FOLD ELECTRON DENSITY',/,
     &          ' ====================='
     &         /' Input map:       ',a,
     &         /' Output map:      ',a,
     &         /' Input envelope:  ',a,
     &         /' Output envelope: ',a)
      else
        write(6,111) imap,omap,ienv,oenv
        if(record)write(nrec,111) imap,omap,ienv,oenv
111     format(' FOLD ELECTRON DENSITY',/,
     &    ' I/P map:',a,'  O/P map:',a,
     &    ' I/P env:',a,'  O/P env:',a)
      endif
c
c--- Turn envelope checking off if envelope names are 'OFF'
c
      i_check=ienv.ne.'OFF'
      o_check=oenv.ne.'OFF'
c
c
c
      if(verbose)then
        if(i_check) write(6,120)
        if(.not.i_check)write(6,130)
        if(o_check) write(6,140)
        if(.not.o_check)write(6,150)
c
        if(record) then
          if(i_check) write(nrec,120)
          if(.not.i_check)write(nrec,130)
          if(o_check) write(nrec,140)
          if(.not.o_check)write(nrec,150)
        endif
      endif
c
120   format(' Input map envelope filtered')
130   format(' Input map not filtered')
140   format(' Output map envelope filtered')
150   format(' Output map not filtered')
c
c
c--- Get the slot number of the input map
c
      call av_ass(map1,imap)
c
c--- The input map should exist by now - if not then its time to leave
c
      if(.not.defined(map1))then
        write(6,160)
160     format(' %FOLD_SET-ERR: Input map empty')
        return
      endif
c
c--- Check status of the output map - if it exists we can jump to the
c--- envelope setup
c
      call av_ass(map2,omap)
      if(defined(map2))then
        write(6,170)
        if(record)write(nrec,170)
170     format(
     &    ' **WARNING** o/p map defined; MUST be a cryst asym unit')
        goto 100
      endif
c
c--- If it doesn't exist, we need to set it up
c
      if(.not.terse)then
        write(6,180)
        if(record)write(nrec,180)
180     format(' The output map limits must define an asymmetric unit')
      endif
c
c--- Set up output map header information
c
      call copy_head(map1,map2)
c
c--- Check if we want to redefine the output grid (eg coarser than input)
c
20    nopt=keyword('AV_FOLD_UPDATE','UPDA!RESE!BACK!QLIM!GO  ',0)
c
c---- UPDA - change map limits for output map
c
      if(nopt.eq.1) then
        do j=1,3
          nstart(j,map2)=intval('AV_AVERAGE:NSTART',-3000,3000,0)
          nend(j,map2)  =intval('AV_AVERAGE:NEND',-3000,3000,0)
          nunit(j,map2) =intval('AV_AVERAGE:NUNIT',-3000,3000,0)
          nx(j,map2)=nend(j,map2)-nstart(j,map2)+1
          norg(j,map2)=1-nstart(j,map2)
          if(.not.terse)then
            write(6,200)nstart(j,map2),nend(j,map2),
     &                  nx(j,map2),nunit(j,map2)
            if(record)
     &          write(nrec,200)nstart(j,map2),nend(j,map2),
     &                  nx(j,map2),nunit(j,map2)
200         format(' ',4i5)
          endif
        enddo
        npix(map2)=nx(1,map2)*nx(2,map2)*nx(3,map2)
        write(6,210)npix(map2)
        if(record)write(nrec,210)npix(map2)
210     format(' There will be',i10,' pixels in the output map')
      endif
c
c---- RESE - reset axis ordering
c
      if(nopt.eq.2)then
        do j=1,3
          iuvw(j,map2)=intval('AV_AVERAGE:ENTER_AXIS',1,3,0)
        enddo
      endif
c
c---- BACK - run away!
c
      if(NOPT.eq.3)then
        quit=.true.
        return
      endif
c
c---- QLIM - print the current map limits
c
      if(NOPT.eq.4)then
        write(6,214)(nstart(j,map2),nend(j,map2),nunit(j,map2),
     &      nx(j,map2),norg(j,map2),j=1,3)
214     format( 3(' START:',i4,' END:',i4,' GRID:',i4,
     &              ' POINTS:',i4,' ORIGIN:',i4,/) )
      endif
c
c---- Wait for the GO
c
      if(nopt.ne.5)goto 20
c
c--- Set up cell data
c
      call XRfrac(XRcell(1,map2),XRtr(1,1,map2),XRintr(1,1,map2),
     &    XRvol,.true.)
c
c--- Set up symmetry data
c
      call XR_symsy(XRcell(1,map2),lgrp(map2),ics_local,
     &                  n_xsym(map2),XR_sym(1,map2) )
c
c--- And integerized CSO's:
c
      call XR_isym(map2)
c
c--- Set up pointers indicating which symm op to try first
c
      do j=1,maxsym
        last_sym(j,map2)=1
      enddo
c
c--- Set up bricking data
c
      call brickit(map2)
c
c--- Allocate space for output map
c
      call alloc(map2)
c
c--- Set up allocations for envelopes
c
100   if(i_check)call av_ass(env1,ienv)
      if(o_check)call av_ass(env2,oenv)
c
c--- Check that the grid matches on the input and output maps
c--- Time to abandon ship if it doesn't
c
      same=same_grid(map1,map2)
      if(.not.same)then
        write(6,49)
        if(record)write(nrec,49)
49      format(' %FOLD_SET-ERR: MAPO and MAPI not on same grid')
        goto 700
      endif
c
c--- The same for the input map and envelope
c
      if(i_check)then
        same=same_grid(map1,env1)
        if(.not.same)then
          write(6,50)
          if(record)write(nrec,50)
50        format(' %FOLD_SET-ERR: ENVI and MAPI not on same grid')
          goto 700
        endif
      endif
c
c--- The same for the input map and output envelope
c
      if(o_check)then
        same=same_grid(map1,env2)
        if(.not.same)then
          write(6,51)
          if(record)write(nrec,51)
51        format(' %FOLD_SET-ERR: ENVO and MAPI not on same grid')
          goto 700
        endif
      endif
c
c--- And actually do it 
c---  map1    =  input map number
c---  map2    =  output map number
c---  env1    =  input envelope number
c---  env2    =  output envelope number
c---  i_check =  use input envelope?
c---  o_check =  use output envelope?
c
      call fold_av(map1,map2,env1,env2,i_check,o_check)
c
c--- Its all over by now
c
700   return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE frac_prot
c     ==================
      IMPLICIT NONE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      INCLUDE 'average.fcm'
c     =====================
c
c----local variables:
      integer   map1,map2,map3,env1,env2,nout,np,nc
      integer   iindex(3),new_index(3)
      real      aindex(3),x(3),x1(3)
      logical   i_check,o_check,use_in,no_trans,trans,check
      logical   iiinout,iinout,quit
      integer*2 irho,fsol
      external  gp3,gp3_a
      integer*2 gp3,gp3_a
c
      integer   nsym1,nsym2,itype,ngadd,rna,prot
      integer   iz,jz,iy,jy,ix,jx,jjz,jjy,jjx,j
c
c      equivalence  (i_temp,ch_temp)
c
      use_in=.true.
      check=.true.
      iiinout=.false.
c
      call frac_prot_setup(map1,map2,map3,env1,env2,nsym1,nsym2,i_check,
     &                o_check,rna,prot,no_trans,trans,ITYPE,quit)
c
      if(quit)then
        quit=.false.
        return
      endif
c
c
c
c
c---- the averaging process is very simple, we are driven by the
c---- output map, for each pixel in that map (accessed in a sequential
c---- fashion, but since map is bricked we move slowly through space)
c---- we pick up all the symm related pixels in the input map and take
c---- the average.
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
      call init_map(map2)
c
c---- initialize r-factors, corr coefs, etc
c
      ngadd=0
      np=0
      nc=0
      nout=0
c
      call i_stat
c
      if(o_check)call init_map(env2)
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
      do iz=1,n_brick(3,map2)
c
        jz=(iz-1)*brick(3,map2)
c
        do iy=1,n_brick(2,map2)
c
          jy=(iy-1)*brick(2,map2)
c
          do ix=1,n_brick(1,map2)
c
            jx=(ix-1)*brick(1,map2)
c
c---- now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map2)
c
              if(jjz.le.nx(3,map2))then
c
                do jjy=jy+1,jy+brick(2,map2)
c
                  if(jjy.le.nx(2,map2))then
c
                    do jjx=jx+1,jx+brick(1,map2)
c               
                      if(jjx.le.nx(1,map2))then
c
c---- check that this pixel is in required region of the output map
c
      if(o_check)then
        if(gp3(env2,jjx,jjy,jjz).ne.1)goto 200
      endif
c
      iindex(1)=jjx
      iindex(2)=jjy
      iindex(3)=jjz
c
      call j3_to_f(iindex,x1,map2)
      call f_to_o(x1,x,map2)
c
c
      in_temp=0
      sum_temp=zero
c
c
c---- loop over all the required NCS operations:
c
      do j=nsym1,nsym2
c
c---- get the appropriate postion for this copy:
c
        if(trans)    call spin (x,x1,ops(1,1,j),vecs(1,j))
c
        if(no_trans) call vmtply(ops(1,1,j),x,x1)
c
        call o_to_j3(x1,aindex,iindex,map1)
c
c---- if required check that this is in a permitted region of the Input map:
c
        if(i_check)then
c
c---- this code allows for the envelope to be limited to the cryst. aysm.
c---- portion of the map (irrespective of what area the actual maps cover)
c
          iiinout=iinout(env1,iindex(1),iindex(2),iindex(3)
     &                            ,new_index,check,j)
c
c          (use the envelope pixel that is closest to non-integral grid point)
c
          if(.not.iiinout.or.gp3_a(env1,iindex).eq.0)then
            nout=nout+1
c---- pixels outside envelope for i/p counted
            goto 300
          endif
c
        endif
c
c---- now get the interpolated electron density:
c
c--- if not using xtal symm....much quicker
        if(xr_con(map1))then
c
          call g_p_i_i_x(map1,irho,aindex,iindex,itype,j,use_in
     &                                    ,np,nc)
          call g_p_i_i_x(map3,fsol,aindex,iindex,itype,j,use_in
     &                                    ,np,nc)
c
c                            ( itype specifies interp: 1, 8 or 11 point. )
c
        else
          call g_p_i_i(map1,irho,aindex,iindex,itype,use_in,np,nc)
          call g_p_i_i(map3,fsol,aindex,iindex,itype,use_in,np,nc)
        endif
c
c---- if the pixel was found ok add it into the sums:
c
        if(.not.use_in) then
          nout=nout+1
        else
          irho=prot*(irho-rna+(rna*fsol/1000))/(prot-rna)
          sum_temp=sum_temp+irho
          in_temp=in_temp+1
          ir_temp(in_temp)=irho
          ngadd=ngadd+1
        endif
300     continue
      enddo
c
c---- now we have gathered all the contributions.
c---- a_stat evaluates mean rho and accumualtes stats for this output pixel:
c
      call a_stat
c
      irho=nint(sum_temp)
c
c---- and now finally store the pixel in the Output map:
c
      call pp3(map2,jjx,jjy,jjz,irho)
c
200   continue
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c----phew, averaged, lets see how we did ...
c
c999   continue
      if(verbose)then
        write(6,998)ngadd
        if(record)write(nrec,998)ngadd
998     format(' RNA/PROTEIN SEPARATION FINISHED',/,
     &         ' ===============================',/,
     &      i9,' pixels were inspected')
      else
        write(6,997)ngadd
        if(record)write(nrec,997)ngadd
997     format(' RNA/PROT FINISHED, ',i10,' pixels were inspected')
      endif
c
      write(6,1000)np,nout,nc
      if(record)write(nrec,1000)np,nout,nc
1000  format(' # interpolated:',i9,
     &       '  outside i/p map:',i9,
     &       '  not interpolated:',i9)
c
c---- print some stats
c
      call p_stat
c
      if(pipe)then
        call pip('MAPO','MAPI',.true.)
        call pip('    ','MAPO',.false.)
      endif
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE frac_prot_setup(map1,map2,map3,env1,env2,nsym1,nsym2,
     &  i_check,o_check,rna,prot,no_trans,trans,ITYPE,quit)
c     ============================================================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c     =====================
c
c----local variables:
      integer     map1,map2,map3,env1,env2
      character*8 imap,smap,omap,ienv,oenv
      logical     i_check,o_check,no_trans,trans,use_out
c
      character*8 G2CHRU,ch_temp
      real        xrvol
      integer     nopt,keyword,j,intval,ics_local
      integer     itype,nsym1,nsym2,rna,prot
      logical     lass_get,quit
      character*13 mess(0:1)
      data mess/'on all pairs ','from the mean'/
c
c
      no_trans=.false.
      trans=.true.
      use_out=.false.
      correl=.false.
c
c
      if(.not.lass_get(ch_temp,'MAPI'))then
c
        ch_temp=G2CHRU('AV_AVER:INPUT_MAP',1,8,0)
c
        call lass_set(ch_temp,'MAPI')
      endif
      imap=ch_temp
c
      if(.not.lass_get(ch_temp,'MAPS'))then
c
        ch_temp=G2CHRU('AV_AVER:FSOL_MAP',1,8,0)
c
        call lass_set(ch_temp,'MAPS')
      endif
      smap=ch_temp
c
c
      if(.not.lass_get(ch_temp,'MAPO'))then
c
        ch_temp=G2CHRU('AV_AVER:OUTPUT_MAP',1,8,0)
c
        call lass_set(ch_temp,'MAPO')
      endif
      omap=ch_temp
c
c
      if(.not.lass_get(ch_temp,'ENVI'))then
c
        ch_temp=G2CHRU('AV_AVER:INPUT_ENV',1,8,0)
c
        call lass_set(ch_temp,'ENVI')
      endif
      ienv=ch_temp
c
c
      if(.not.lass_get(ch_temp,'ENVO'))then
c
        ch_temp=G2CHRU('AV_AVER:OUTPUT_ENV',1,8,0)
c
        call lass_set(ch_temp,'ENVO')
      endif
      oenv=ch_temp
c
c
      if(verbose)then
        write(6,110) imap,omap,ienv,oenv
        if(record)write(nrec,110) imap,smap,omap,ienv,oenv
110     format(/' SEPARATE RNA AND PROTEIN DENSITY',/,
     &       ' ========================'
     &      /' Input map:       ',a,
     &      /' Solvent map:     ',a,
     &      /' Output map:      ',a,
     &      /' Input envelope:  ',a,
     &      /' Output envelope: ',a)
      else
        write(6,111) imap,omap,ienv,oenv
        if(record)write(nrec,111) imap,smap,omap,ienv,oenv
111     format(' SEPARATE RNA AND PROTEIN DENSITY',/,
     &    ' I/P map:',a,'  SOL map:',a,'  O/P map:',a,'  I/P env:',a,
     &    '  O/P env:',a)
      endif
C
      i_check=ienv.ne.'OFF'
      o_check=oenv.ne.'OFF'
c
      if(verbose)then
        if(i_check) write(6,120)
        if(.not.i_check)write(6,130)
        if(o_check) write(6,140)
        if(.not.o_check)write(6,150)
c
        if(record) then
          if(i_check) write(nrec,120)
          if(.not.i_check)write(nrec,130)
          if(o_check) write(nrec,140)
          if(.not.o_check)write(nrec,150)
        endif
      endif
c
120   format(' Input map envelope filtered')
130   format(' Input map not filtered')
140   format(' Output map envelope filtered')
150   format(' Output map not filtered')
c
c
c       get slot numbers for all maps
c
      call av_ass(map1,imap)
c----map 1 should exist by now !
      if(.not.defined(map1)) then
        write(6,160)
        if(record)write(nrec,160)
160     format(' %AVERAGE-ERR: Input map empty')
        return
      endif
      call av_ass(map3,smap)
c----map 3 should exist by now !
      if(.not.defined(map3)) then
        write(6,161)
        if(record)write(nrec,161)
161     format(' %AVERAGE-ERR: Solvent map empty')
        return
      endif
c---- check status of output file
      call av_ass(map2,omap)
      if(.not.terse)then
        if(defined(map2))write(6,170)
        if(.not.defined(map2))write(6,180)
c
        if(record)then
          if(defined(map2))write(nrec,170)
          if(.not.defined(map2))write(nrec,180)
        endif
      endif
c
170   format(' Output map exists: will overwrite')
180   format(' Make output from i/p template')
c
c---- set up output map header information
c
      if(.not.defined(map2))then
        call copy_head(map1,map2)
c
c---- check if we want to redefine the output grid (eg coarser than input)
c
      quit=.false.
c
20    nopt=keyword
     &  ('AV_AVERAGE_UPDATE','UPDA!RESE!NOTR!BACK!QLIM!GO  ',0)
c
      if(nopt.eq.1) then
        do j=1,3
          nstart(j,map2)=intval('AV_AVERAGE:NSTART',-10000,10000,0)
          nend(j,map2)  =intval('AV_AVERAGE:NEND',-10000,10000,0)
          nunit(j,map2) =intval('AV_AVERAGE:NUNIT',-10000,10000,0)
          nx(j,map2)=nend(j,map2)-nstart(j,map2)+1
          norg(j,map2)=1-nstart(j,map2)
          if(.not.terse)then
            write(6,200)nstart(j,map2),nend(j,map2),
     &        nx(j,map2),nunit(j,map2)
            if(record)
     &        write(nrec,200)nstart(j,map2),nend(j,map2),
     &        nx(j,map2),nunit(j,map2)
200         format(' ',4i5)
          endif
        enddo
        npix(map2)=nx(1,map2)*nx(2,map2)*nx(3,map2)
        write(6,210)npix(map2)
        if(record)write(nrec,210)npix(map2)
210     format(' There will be',i10,' pixels in the output map')
      endif
c
      if(nopt.eq.2)then
        do j=1,3
          iuvw(j,map2)=intval('AV_AVERAGE:ENTER_AXIS',1,3,0)
        enddo
      endif
c
      if(NOPT.EQ.3) then
        no_trans=.true.
        trans=.false.
      endif
c
      if(NOPT.eq.4)then
        quit=.true.
        return
      endif
c
      if(NOPT.eq.5)then
        write(6,214)(nstart(j,map2),nend(j,map2),nunit(j,map2),
     &    nx(j,map2),norg(j,map2),j=1,3)
214     format( 3(' START:',i5,' END:',i5,' GRID:',i5,
     &            ' POINTS:',i5,' ORIGIN:',i5,/) )
      endif
c
      if(nopt.ne.6)goto 20
c
c---- set up cell data
c
        call XRfrac(XRcell(1,map2),XRtr(1,1,map2),XRintr(1,1,map2),
     &              XRvol,.true.)
c
c---- set up symmetry data
c
        call XR_symsy(XRcell(1,map2),lgrp(map2),ics_local,
     &                n_xsym(map2),XR_sym(1,map2) )
c
c---- and integerized CSO's:
c
        call XR_isym(map2)
c
c---- also set up pointers indicating which symm op to try first
c
        do j=1,maxsym
          last_sym(j,map2)=1
        enddo
c
c
c---- set up bricking data
c
        call brickit(map2)
c
c---- allocate space for output map
c
        call alloc(map2)
c
      endif
c
c---- set up allocations for envelopes
c
      if(i_check) call av_ass(env1,ienv)
      if(o_check) call av_ass(env2,oenv)
c
c
c
c---- list symmetry operators to be used
c
      itype=11
      istat=1
      nsym1=1
      nsym2=nsym
      prot=0
      rna=1
c
10    nopt=keyword
     &  ('AV_AVERAGE_SYMM','ALL !STAR!END !GO  !INTP!CORQ!CORA!'//
     &           'RNA !PROT',0)
c
      if(nopt.eq.2)nsym1=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.eq.3)nsym2=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.eq.5)itype=intval('AV_AVERAGE_INTP',1,64,0)
      if(nopt.eq.6)istat=1
      if(nopt.eq.7)istat=0
      if(NOPT.eq.8)rna=intval('AV_RNA_LEVEL',-32000,32000,0)
      if(nopt.eq.9)prot=intval('AV_PROT_LEVEL',-32000,32000,0)
c
      if(nopt.ne.4)goto 10
c
      if( (nsym1.gt.nsym2).or.(nsym2.gt.nsym) )then
        write(6,230)
        if(record)write(nrec,230)
230     format(' %AVERAGE-ERR: In choice of symmetry operators ')
      endif
c
      if(verbose)then
        do j=nsym1,nsym2
          call prtsym(6,j)
        enddo
      endif
C----
c
      if(itype.ne.11.and.itype.ne.8.and.itype.ne.1.and.itype.ne.64)
     &  itype=11
      write(6,240)itype,mess(istat)
      if(record)write(nrec,240)itype,mess(istat)
240   format(i3,' point interpolation. Correlation coeff calculated ',a)
c
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE frac_sol
c     ====================
      IMPLICIT NONE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Knowledge of solvent density is used to calculate the fractional
c solvent content and non-solvent density of each pixel, and the corrected
c value is put back.
c
      INCLUDE 'average.fcm'
c     =====================
c
c----local variables:
c
      integer   maps(maxsym),map2,map3,envs(maxsym),env2,nsym1,nsym2
      logical   i_checks(maxsym),o_check,quit
c
      integer   j,iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
      external  gp3
      integer*2 gp3,irho,fsol2,poth2
c
      integer   nobs,nout,sum_nobs,jj,iin(maxsym),nover,nunder,
     &          sum_dev_sq_c(maxsym)
      real*8    mean_fsol,sig_fsol,mean_poth,sig_poth,sumx,sumxx,sumxy,
     &          sumy,denom,sum_dev_sq(maxsym),rms_dev,solv(maxsym),
     &          fsol,max_fsol,min_fsol,sum_fsol,sum_fsol_sq,
     &            poth,max_poth,min_poth,sum_poth,sum_poth_sq,
     &            irho_fs(maxsym),av_nobs,irho_r
c
c---- call setup routine
c
      call frac_sol_setup(maps,map2,map3,envs,env2,nsym1,nsym2,i_checks,
     &                o_check,quit)
c
      if(quit)then
        quit=.false.
        return
      endif
c
c---- the averaging process is very simple, we are driven by the
c---- output map, for each pixel in that map (accessed in a sequential
c---- fashion, but since map is bricked we move slowly through space)
c---- we pick up all the symm related pixels in the input map and take
c---- the average.
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
      call init_map(map2)
      if(o_check)call init_map(env2)
c
c---- initialize stats
c
      max_fsol=-32000
      min_fsol=32000
      sum_fsol=0
      sum_fsol_sq=0
      max_poth=-32000
      min_poth=32000
      sum_poth=0
      sum_poth_sq=0
      nobs=0
      nout=0
      sum_nobs=0
      do j=1,maxsym
        sum_dev_sq(j)=0.0d0
        sum_dev_sq_c(j)=0
      enddo
      nover=0
      nunder=0
c
c---- get some calculations out of the way
c
      do j=nsym1,nsym2
        solv(j)=(float(sollev(maps(j)))-moffset(maps(j)))
     &    /mscale(maps(j))
      enddo
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
      do iz=1,n_brick(3,map2)
c
        jz=(iz-1)*brick(3,map2)
c
        do iy=1,n_brick(2,map2)
c
          jy=(iy-1)*brick(2,map2)
c
          do ix=1,n_brick(1,map2)
c
            jx=(ix-1)*brick(1,map2)
c
c---- now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map2)
c
              if(jjz.le.nx(3,map2))then
c
                do jjy=jy+1,jy+brick(2,map2)
c
                  if(jjy.le.nx(2,map2))then
c
                    do jjx=jx+1,jx+brick(1,map2)
c               
                      if(jjx.le.nx(1,map2))then
c
c---- check that this pixel is in required region of the output map
c
      if(o_check)then
        if(gp3(env2,jjx,jjy,jjz).ne.1)goto 200
      endif
c
c---- loop over all the required maps:
c
      do j=nsym1,nsym2
c
c---- if required check that this is in a permitted region of the input map.
c
        if(i_checks(j))then
          if(gp3(envs(j),jjx,jjy,jjz).eq.0)goto 300
        endif
c
c---- now get the electron density.
c
        irho=gp3(maps(j),jjx,jjy,jjz)
c
c---- add it into the storage arrays
c
        nobs=nobs+1
        irho_fs(nobs)=(float(irho)-moffset(maps(j)))/mscale(maps(j))
        iin(nobs)=j
c
300     continue
      enddo
c
c---- now we have gathered all the contributions...
c
      if(nobs.gt.1)then
c
        sumx=0
        sumxx=0
        sumxy=0
        sumy=0
c
        do jj=1,nobs
          sumx=sumx+solv(iin(jj))
          sumxx=sumxx+(solv(iin(jj))**2)
          sumxy=sumxy+(solv(iin(jj))*irho_fs(jj))
          sumy=sumy+irho_fs(jj)
        enddo
c
c---- perform calculation.
c
        denom=(nobs*sumxx)-(sumx**2)
        fsol=((nobs*sumxy)-(sumx*sumy))/denom
        poth=((sumy*sumxx)-(sumx*sumxy))/denom
c
c---- stats
c
        if(fsol.gt.max_fsol)max_fsol=fsol
        if(fsol.lt.min_fsol)min_fsol=fsol
c
c---- reality checks
c
        if(fsol.gt.1)then
          fsol=1.0d0
          poth=0.0d0
          nover=nover+1
        else if(fsol.lt.0)then
          fsol=0.0d0
          poth=sumy/nobs
          nunder=nunder+1
        endif
c
c---- stats
c
        sum_fsol=sum_fsol+fsol
        sum_fsol_sq=sum_fsol_sq+fsol**2
c
        if(poth.gt.max_poth)max_poth=poth
        if(poth.lt.min_poth)min_poth=poth
        sum_poth=sum_poth+poth
        sum_poth_sq=sum_poth_sq+poth**2
c
c---- write out the corrected values
c
        do jj=1,nobs
          irho_r=(fsol*solv(iin(jj))+poth)
          irho=nint((mscale(iin(jj))*irho_r)+moffset(iin(jj)))
          call pp3(maps(iin(jj)),jjx,jjy,jjz,irho)
c
c----- stats
c
          sum_dev_sq(iin(jj))=sum_dev_sq(iin(jj))+
     &      (irho_fs(jj)-(irho_r))**2
          sum_dev_sq_c(iin(jj))=sum_dev_sq_c(iin(jj))+1
c
        enddo
c
        nout=nout+1
        sum_nobs=sum_nobs+nobs
        nobs=0
c
        fsol2=nint(fsol*10000)
        poth2=nint(poth)
        call pp3(map2,jjx,jjy,jjz,fsol2)
        call pp3(map3,jjx,jjy,jjz,poth2)
c
      endif
c
c---- and now end all the loops:
c
200   continue
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c---- phew, finished. Lets see how we did...
c
      mean_fsol=sum_fsol/nout
      mean_poth=sum_poth/nout
c
      sig_fsol=sqrt((sum_fsol_sq/nout)-mean_fsol**2)
      sig_poth=sqrt((sum_poth_sq/nout)-mean_poth**2)
c
      av_nobs=sum_nobs/nout
c
      write(6,10)nout,sum_nobs,av_nobs,nover,nunder,max_fsol,
     &  min_fsol,mean_fsol,sig_fsol,max_poth,min_poth,mean_poth,sig_poth
c
      if(record)write(nrec,10)nout,sum_nobs,av_nobs,nover,nunder,
     &  max_fsol,min_fsol,mean_fsol,sig_fsol,max_poth,min_poth,
     &  mean_poth,sig_poth
c
10    format(
     &  ' FRACTIONAL SOLVENT CALCULATION STATISTICS',/
     &  ' Total # of ouput pixels:',i10,/
     &  ' Total # of input pixels:',i10,' (',f8.3,
     &  ' input pixels per output pixel)',/
     &  ' # Pixels > 100% solvent',i10,' # Pixels < 0% solvent',i10,/
     &  ' Fractional solvent content map:',/
     &  '   Max:',f8.3,'  Min:',f8.3,'  Mean:',f8.3,'  SD:',f8.3,/
     &  ' Non-solvent density map:',/
     &  '   Max:',f12.3,'  Min:',f12.3,'  Mean:',f12.3,'  SD:',f12.3,/)
c
      do jj=1,maxsym
        if(sum_dev_sq_c(jj).gt.0)then
          rms_dev=sqrt(sum_dev_sq(jj)/sum_dev_sq_c(jj))
          write(6,20)jj,sum_dev_sq_c(jj),jj,rms_dev
          if(record)write(nrec,20)jj,sum_dev_sq_c(jj),jj,rms_dev
20        format(
     &      ' Total # of output pixels for map ',i2,':',i10,/,
     &      ' RMS Deviation from fitted values for map ',i2,':',f20.5)
        endif
      enddo
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE frac_sol_setup(maps,map2,map3,envs,env2,nsym1,nsym2,
     &  i_checks,o_check,quit)
c     ===============================================================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c     =====================
c
c----local variables:
c
      integer         maps(maxsym),map2,map3,envs(maxsym),env2,
     &                nopt,keyword,j,intval,ics_local,nsym1,nsym2
      integer         keylist,n_read
      real            xrvol
      character*8     ompf,ompo,oenv,G2CHRU,ch_temp
      logical         i_checks(maxsym),o_check,lass_get,quit
c
c---- initialise logicals
c
      quit=.false.
c
c---- get number of maps
c
      n_read=intval('AV_FSOL:NUMBER_OF_MAPS',1,MMAPS,0)
c
c---- get map and envelope labels
c
      do j=1,n_read
        maps(j)=keylist('AV_FSOL:MAP_LABEL',LABELS,NMAPS,0)
c
        envs(j)=0
        if('ENV'.eq.G2CHRU('AV_FSOL:ENV/NOE',0,3,0))
     &    envs(j)=keylist('AV_FSOL:ENV_LABEL',LABELS,NMAPS,0)
        i_checks(j)=envs(j).gt.0
      enddo
c
      write(6,1003)'Number of operators to be used:',n_read
      if(record)write(nrec,1003)'Number of operators to be used:',n_read
c
      do j=1,n_read
        if(i_checks(j))write(6,1001)j,maps(j),envs(j)
        if(.not.i_checks(j))write(6,1002)j,maps(j)
        if(record)then
          if(i_checks(j))write(nrec,1001)j,maps(j),envs(j)
          if(.not.i_checks(j))write(nrec,1002)j,maps(j)
        endif
      enddo
c
1001  format(' Operator',i4,' will use map',i3,' and Input env',i3)
1002  format(' Operator',i4,' will use map',i3,' and no Input env')
1003  format(' ',a,i4)
c
c---- get output maps
c
      if(.not.lass_get(ch_temp,'MAPF'))then
        ch_temp=G2CHRU('AV_FSOL:OUTPUT_FSOL_MAP',1,8,0)
        call lass_set(ch_temp,'MAPF')
      endif
      ompf=ch_temp
c
      if(.not.lass_get(ch_temp,'MAPO'))then
        ch_temp=G2CHRU('AV_FSOL:OUTPUT_POTHER_MAP',1,8,0)
        call lass_set(ch_temp,'MAPO')
      endif
      ompo=ch_temp
c
c---- get output envelopes
c
      if(.not.lass_get(ch_temp,'ENVO'))then
        ch_temp=G2CHRU('AV_FSOL:OUTPUT_ENV',1,8,0)
        call lass_set(ch_temp,'ENVO')
      endif
      oenv=ch_temp
      o_check=oenv.ne.'OFF'
c
c---- confirm in writing
c
      if(verbose)then
        write(6,110) ompf,ompo,oenv
        if(record)write(nrec,110) ompf,ompo,oenv
      else
        write(6,111) ompf,ompo,oenv
        if(record)write(nrec,111) ompf,ompo,oenv
      endif
c
      if(verbose)then
        if(o_check) write(6,140)
        if(.not.o_check)write(6,150)
        if(record) then
          if(o_check) write(nrec,140)
          if(.not.o_check)write(nrec,150)
        endif
      endif
c
110   format(/' FRACTIONAL SOLVENT CALCULATION',/,
     &        ' ==============================',/,
     &        ' Fraction solvent output map:      ',a,/,
     &        ' Other density output map:      ',a,/,
     &        ' Output envelope: ',a)
111   format(' FRACTIONAL SOLVENT CALCULATION',
     &  ' F_Sol O/P map:',a,'  P_Oth O/P map:',a,'  O/P env:',a)
140   format(' Output maps envelope filtered')
150   format(' Output maps not filtered')
c
c---- check status of output files
c
      call av_ass(map2,ompf)
      call av_ass(map3,ompo)
      if(.not.terse)then
        if(defined(map2).and.defined(map3))write(6,170)
        if((.not.defined(map2)).or.(.not.defined(map3)))write(6,180)
        if(record)then
          if(defined(map2).and.defined(map3))write(nrec,170)
          if((.not.defined(map2)).or.(.not.defined(map3)))
     &      write(nrec,180)
        endif
      endif
c
170   format(' F_Sol Output Map and P_Oth Output Map exist: ',
     &  'will overwrite')
180   format(' Make F_Sol and P_Oth Output Maps from i/p template')
c
c---- set up output map header information
c
      call copy_head(maps(1),map2)
      mscale(map2)=1000.0
      call copy_head(map2,map3)
c
c---- set up cell data
c
      call XRfrac(XRcell(1,map2),XRtr(1,1,map2),XRintr(1,1,map2),
     &            XRvol,.true.)
      call XRfrac(XRcell(1,map3),XRtr(1,1,map3),XRintr(1,1,map3),
     &            XRvol,.true.)
c
c---- set up symmetry data
c
      call XR_symsy(XRcell(1,map2),lgrp(map2),ics_local,
     &              n_xsym(map2),XR_sym(1,map2) )
      call XR_symsy(XRcell(1,map3),lgrp(map3),ics_local,
     &              n_xsym(map3),XR_sym(1,map3) )
c
c---- and integerized CSO's:
c
      call XR_isym(map2)
      call XR_isym(map3)
c
c
c---- set up bricking data
c
      call brickit(map2)
      call brickit(map3)
c
c---- allocate space for output map
c
      call alloc(map2)
      call alloc(map3)
c
c---- set up allocations for envelopes
c
      if(o_check) call av_ass(env2,oenv)
c
c---- list symmetry operators to be used
c
      nsym1=1
      nsym2=n_read
c
10    nopt=keyword
     &  ('AV_FSOL_SYMM','ALL !STAR!END !BACK!GO  ',0)
c
      if(nopt.eq.2)nsym1=intval('AV_FSOL_SYMM',1,n_read,0)
      if(nopt.eq.3)nsym2=intval('AV_FSOL_SYMM',1,n_read,0)
      if(nopt.eq.4)quit=.true.
      if((nopt.ne.4).and.(nopt.ne.5))goto 10
c
      write(6,221)nsym1,nsym2
      if(record)write(nrec,221)nsym1,nsym2
221   format(' Form average density from i/p map using maps',i4,' to',
     &  i4)
c
      if( (nsym1.gt.nsym2).or.(nsym2.gt.n_read) )then
        write(6,230)
        if(record)write(nrec,230)
230     format(' %AVERAGE-ERR: In choice of maps ')
      endif
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE g_p_i_i(map,irho,aj3,j3,itype,in_map,np,nc)
c       ==================================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c---- interpolation, 64, 11, 8 or 1 point !
c
c Nordman's 11 (or 8) point interpolation
c
        real aj3(3),uu(3),u,v,w,p000,p100,p010,p001,p110,p011,p101,p111
        real pn00,p0n0,p00n,sum_r,a,b,c,d,p8,p11,constant
        integer j3(3),ind(3),off(3),np,nc,j,map,itype
        equivalence (uu(1),u),(uu(2),v),(uu(3),w)
        logical in_map
        integer*2 irho
        real del
c
        parameter (constant=0.4)
c
        external gp3_a,gp3
        integer*2 gp3_a,gp3
c
c
        in_map=.true.
c
        do j=1,3
          if( (j3(j).gt.nx(j,map)) .or. (j3(j).lt.1) )then
            irho=0
            in_map=.false.
            return
          endif
          uu(j)=(aj3(j)-j3(j))
          if(uu(j).lt.(0.0)) then
                uu(j)=abs(uu(j))
                ind(j)=-1
            else
                ind(j)=1
            endif
c
        enddo
c
        sum_r=uu(1)+uu(2)+uu(3)
        if(sum_r.lt.0.1.or.itype.eq.1)then
          irho = gp3_a(map,j3)
          nc=nc+1
          return
        endif
c
        do j=1,3
          off(j)=ind(j)+j3(j)
c---- if on edge take nearest pixel
          if( (off(j).gt.nx(j,map)) .or. (off(j).lt.1)  )then
            irho = gp3_a(map,j3)
            nc=nc+1
            return
          endif
        enddo
c
        p000 = gp3(map,j3(1) ,j3(2) ,j3(3) )
        p100 = gp3(map,off(1),j3(2) ,j3(3) )
        p010 = gp3(map,j3(1) ,off(2),j3(3) )
        p001 = gp3(map,j3(1) ,j3(2) ,off(3))
        p110 = gp3(map,off(1),off(2),j3(3) )
        p011 = gp3(map,j3(1) ,off(2),off(3))
        p101 = gp3(map,off(1),j3(2) ,off(3))
        p111 = gp3(map,off(1),off(2),off(3))
c
c
        a=p100-p000
        b=p010-p000
        c=p110-p010
        d=p101-p001
c
        p8 = p000+u*(a+w*(-a+d)+v*((c-a)+w*( a-c-d-p011+p111)))
     &            + v*(b+w*(-p001+p011-b))+w*(-p000+p001)
c
c
        if (itype.eq.11) then
          do j=1,3
            off(j)=-ind(j)+j3(j)
            if( (off(j).gt.nx(j,map)) .or. (off(j).lt.1)  )then
              p11=p8
              goto 8
            endif
          enddo
          pn00 = gp3(map,off(1),j3(2),j3(3))
          p0n0 = gp3(map,j3(1),off(2),j3(3))
          p00n = gp3(map,j3(1),j3(2),off(3))
c
          del = (p000-0.5*p100-0.5*pn00)*(u-u**2) +
     &       (p000-0.5*p010-0.5*p0n0)*(v-v**2) +
     &       (p000-0.5*p001-0.5*p00n)*(w-w**2)
c
c
          p11 = p8 + constant*del
c
8       continue
        else
c
          p11 = p8
c
        endif
c
        irho = nint(p11)
c
        np=np+1
c
c

        return
c
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      SUBROUTINE g_p_i_i_x(map,irho,aj3,j3,itype,nxt,in_map,np,nc)
c       ========================================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c---- interpolation, 64, 11, 8 or 1 point !
c
c Nordman's 11 (or 8) point interpolation
c
        real aj3(3),uu(3),u,v,w,p000,p100,p010,p001,p110,p011,p101,p111
        real pn00,p0n0,p00n,sum_r,a,b,c,d,p8,p11,constant
        integer j,map,itype
        integer j3(3),ind(3),off(3),new_index(3),nxt,np,nc
        equivalence (uu(1),u),(uu(2),v),(uu(3),w)
        logical iiinout,check,in_map
        integer*2 irho
        real del
c
        parameter (constant=0.4)
c
        external gp3_a,gp3,iinout
        logical  iinout
        integer*2 gp3_a,gp3
c
        iiinout=.false.
        in_map=.true.
        irho=0
c
c----could setup check in eg. av_setup....searches all csops
        check=.true.
c--- initialise counters
c       
c
c
c----here we set up uu as the fractional indices for interpolation.
c----they are always positive in Nordmann's algorithm so use ind as
c----the sign.
c
        do j=1,3
          uu(j)=(aj3(j)-j3(j))
          if(uu(j).le.(0.0)) then
                uu(j)=abs(uu(j))
                ind(j)=-1
            else
                ind(j)=1
            endif
c
        enddo
c
        sum_r=uu(1)+uu(2)+uu(3)
c
          iiinout=iinout(map,j3(1),j3(2),j3(3),new_index,check,nxt)
          if(iiinout)then
            irho = gp3_a(map,new_index)
          else
            in_map=.false.
            return
          endif
c
c----if we are nearly integral don't interpolate
c
        if(sum_r.lt.0.1.or.itype.eq.1)then
         nc=nc+1
         return
        endif
c
        np=np+1

c
        do j=1,3
          off(j)=ind(j)+j3(j)
        enddo
c
        iiinout=iinout(map,j3(1),j3(2),j3(3),new_index,check,nxt)
        if(.not.iiinout)return
        p000 = gp3_a(map,new_index)
c
        iiinout=iinout(map,off(1),j3(2),j3(3),new_index,check,nxt)
        if(.not.iiinout)return
        p100 = gp3_a(map,new_index )
c
        iiinout=iinout(map,j3(1),off(2),j3(3),new_index,check,nxt)
        if(.not.iiinout)return
        p010 = gp3_a(map,new_index )
c
        iiinout=iinout(map,j3(1),j3(2),off(3),new_index,check,nxt)
        if(.not.iiinout)return
        p001 = gp3_a(map,new_index)
c
        iiinout=iinout(map,off(1),off(2),j3(3),new_index,check,nxt)
        if(.not.iiinout)return
        p110 = gp3_a(map,new_index )
c
        iiinout=iinout(map,j3(1),off(2),off(3),new_index,check,nxt)
        if(.not.iiinout)return
        p011 = gp3_a(map,new_index)
c
        iiinout=iinout(map,off(1),j3(2),off(3),new_index,check,nxt)
        if(.not.iiinout)return
        p101 = gp3_a(map,new_index)
c
        iiinout=iinout(map,off(1),off(2),off(3),new_index,check,nxt)
        if(.not.iiinout)return
        p111 = gp3_a(map,new_index)
c
        a=p100-p000
        b=p010-p000
        c=p110-p010
        d=p101-p001
c
        p8 = p000+u*(a+w*(-a+d)+v*((c-a)+w*( a-c-d-p011+p111)))
     &            + v*(b+w*(-p001+p011-b))+w*(-p000+p001)
c
c
         irho=nint(p8)
c
        if (itype.eq.11) then
          do j=1,3
            off(j)=-ind(j)+j3(j)
          enddo
c
          iiinout=iinout(map,off(1),j3(2),j3(3),new_index,check,nxt)
          if(.not.iiinout)return
          pn00 = gp3_a(map,new_index)
c
          iiinout=iinout(map,j3(1),off(2),j3(3),new_index,check,nxt)
          if(.not.iiinout)return
          p0n0 = gp3_a(map,new_index)
c
          iiinout=iinout(map,j3(1),j3(2),off(3),new_index,check,nxt)
          if(.not.iiinout)return
          p00n = gp3_a(map,new_index)
c
          del = (p000-0.5*p100-0.5*pn00)*(u-u**2) +
     &       (p000-0.5*p010-0.5*p0n0)*(v-v**2) +
     &       (p000-0.5*p001-0.5*p00n)*(w-w**2)
c
c
          p11 = p8 + constant*del
          irho = nint(p11)
c
        endif
c
c
c
c
        return
c
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE gcut_setup(map1,env1,i_check)
c       ========================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c----local variables:
        integer     map1,env1
        character*8 G2CHRU,ch_temp,imap,ienv
        logical     lass_get,i_check
c
c
        if(.not.lass_get(ch_temp,'MAPI'))then
c
          ch_temp=G2CHRU('AV_CUT:INPUT_MAP',1,8,0)
c
          call lass_set(ch_temp,'MAPI')
        endif
          imap=ch_temp
c
c
        if(.not.lass_get(ch_temp,'ENVI'))then
c
        ch_temp=G2CHRU('AV_CUT:INPUT_ENV',1,8,0)
c
          call lass_set(ch_temp,'ENVI')
        endif
        ienv=ch_temp
c
c
        if(verbose)then
          write(6,110) imap,ienv
          if(record)write(nrec,110) imap,ienv
110       format(/' CUT ELECTRON DENSITY',/,
     &           ' ===================='
     &      /' Input map:       ',a
     &      /' Input envelope:  ',a)
          else
          write(6,111) imap,ienv
          if(record)write(nrec,111) imap,ienv
111       format(' CUT ELECTRON DENSITY',
     &    '  I/P map:',a,'  I/P env:',a)
          endif
C
        i_check=ienv.ne.'OFF'
c
c
        if(verbose)then
          if(i_check) write(6,120)
          if(.not.i_check)write(6,130)
c
          if(record)then
            if(i_check) write(nrec,120)
            if(.not.i_check)write(nrec,130)
            endif
          endif
c
120     format(' Input map envelope filtered')
130     format(' Input map not filtered')
c
c
c       get slot numbers for all maps
c
        call av_ass(map1,imap)
c----map 1 should exist by now !
        if(.not.defined(map1)) then
          write(6,140)
          if(record)write(nrec,140)
140       format(' %GCUT_SETUP-ERR: Input map empty')
          return
          endif
c---- check status of output file
c
c---- set up allocations for envelopes
c
        if(i_check) call av_ass(env1,ienv)
c
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE gen
c       ====================
      IMPLICIT NONE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Following on from the fractional solvent content calculation, the non-
c solvent density is corrected to 100% occupancy.
c
      INCLUDE 'average.fcm'
c       =====================
c
c----local variables:
c
        integer   maps(maxsym),map2,envs(maxsym),env2
        real      frac
        logical   i_checks(maxsym),o_check,quit
c
        integer   iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
c        logical   use_in
        external  gp3
        integer*2 gp3,n_s,fs,fs_max,fs_min,ns_max,ns_min,
     &            new_ns_max,new_ns_min
        integer   ntested,nout,nover
        real      fs1,nout_pc
c
c---- call setup routine
c
        call gen_setup(maps,map2,envs,env2,i_checks,o_check,frac,quit)
c
        if(quit)then
         quit=.false.
         return
        endif
c
c---- the averaging process is very simple, we are driven by the
c---- output map, for each pixel in that map (accessed in a sequential
c---- fashion, but since map is bricked we move slowly through space)
c---- we pick up all the symm related pixels in the input map and take
c---- the average.
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
        call init_map(map2)
        if(o_check)call init_map(env2)
c
c---- initialize logicals
c
c        use_in=.true.
c
c---- initialize stats
c
        fs_max=-32000
        fs_min=32000
        ns_max=-32000
        ns_min=32000
        new_ns_max=-32000
        new_ns_min=32000
        ntested=0
        nout=0
        nover=0
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
        do iz=1,n_brick(3,map2)
c
          jz=(iz-1)*brick(3,map2)
c
          do iy=1,n_brick(2,map2)
c
            jy=(iy-1)*brick(2,map2)
c
            do ix=1,n_brick(1,map2)
c
              jx=(ix-1)*brick(1,map2)
c
c---- now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map2)
c
              if(jjz.le.nx(3,map2))then
c
                do jjy=jy+1,jy+brick(2,map2)
c
                if(jjy.le.nx(2,map2))then
c
                  do jjx=jx+1,jx+brick(1,map2)
c               
                  if(jjx.le.nx(1,map2))then
c
c---- check that this pixel is in required region of the output map
c
                    if(o_check)then
                      if(gp3(env2,jjx,jjy,jjz).ne.1)goto 200
                    endif
c
c---- loop over all the required maps:
c
          n_s=gp3(maps(1),jjx,jjy,jjz)
          fs=gp3(maps(2),jjx,jjy,jjz)
          fs1=(fs-moffset(maps(2)))/mscale(maps(2))
          ntested=ntested+1
c
c---- deal only with those passing the fraction cut
c
          if(fs1.lt.frac)then
c
c---- keep track of mins and maxes
c
            if(n_s.lt.ns_min)ns_min=n_s
            if(n_s.gt.ns_max)ns_max=n_s
            if(fs.lt.fs_min)fs_min=fs
            if(fs.gt.fs_max)fs_max=fs
c
c---- perform calculation.
c
            n_s=nint(n_s/(1-fs1))
c
c---- and more mins and maxes
c
            if(n_s.lt.new_ns_min)new_ns_min=n_s
            if(n_s.gt.new_ns_max)new_ns_max=n_s
c
c---- and keep count
c
            nout=nout+1
c
c---- and put the fixed pixel in the output map
c
            call pp3(map2,jjx,jjy,jjz,n_s)
c
c---- keep count of those failing the cut
c
          else
            nover=nover+1
          endif
c
c---- and now end all the loops:
c
200               continue
                  endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
        enddo
c
c---- phew, finished. Lets see how we did...
c
        nout_pc=100*float(nout)/float(ntested)
c
        write(6,10)ntested,nout,nout_pc,nover,fs_max,fs_min,ns_max,
     &  ns_min,new_ns_max,new_ns_min
        if(record)write(nrec,10)ntested,nout,nout_pc,nover,fs_max,
     &  fs_min,ns_max,ns_min,new_ns_max,new_ns_min
c
10      format(
     &  ' FRACTIONAL SOLVENT CALCULATION STATISTICS',/
     &  ' Total # of input pixels:',i10,/
     &  ' Total # of ouput pixels:',i10,' =',f6.2,' % ',/
     &  ' Total # of pixels failing fraction cut:',i10,/
     &  ' Fractional solvent content: Max',i8,' Min',i8,/
     &  ' Non-solvent density in    : Max',i8,' Min',i8,/
     &  ' Non-solvent density out   : Max',i8,' Min',i8,/)

c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE gen_setup(maps,map2,envs,env2,i_checks,o_check,frac,
     &  quit)
c     ===============================================================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c     =====================
c
c----local variables:
c
      integer         maps(maxsym),map2,envs(maxsym),env2,keylist,
     &                nopt,keyword,ics_local
      real            realval,Xrvol,frac
      character*8     omap,oenv,G2CHRU,ch_temp
      logical         i_checks(maxsym),o_check,lass_get,quit
c
c---- initialise logicals
c
      quit=.false.
c
c---- get map and envelope labels
c
      maps(1)=keylist('AV_GEN:NON_SOLVENT_MAP_LABEL',LABELS,NMAPS,0)
c
      envs(1)=0
      if('ENV'.eq.G2CHRU('AV_GEN:ENV/NOE',0,3,0))
     &  envs(1)=keylist('AV_GEN:ENV_LABEL',LABELS,NMAPS,0)
      i_checks(1)=envs(1).gt.0
c
      maps(2)=keylist('AV_GEN:FRAC_SOLVENT_MAP_LABEL',LABELS,NMAPS,0)
c
      envs(2)=0
      if('ENV'.eq.G2CHRU('AV_GEN:ENV/NOE',0,3,0))
     &    envs(2)=keylist('AV_GEN:ENV_LABEL',LABELS,NMAPS,0)
      i_checks(2)=envs(2).gt.0
c
      if(i_checks(1))write(6,1001)maps(1),envs(1)
      if(.not.i_checks(1))write(6,1002)maps(1)
      if(record)then
          if(i_checks(1))write(6,1001)maps(1),envs(1)
          if(.not.i_checks(1))write(6,1002)maps(1)
      endif
c
      if(i_checks(2))write(6,1003)maps(2),envs(2)
      if(.not.i_checks(2))write(6,1004)maps(2)
      if(record)then
          if(i_checks(2))write(6,1003)maps(2),envs(2)
          if(.not.i_checks(2))write(6,1004)maps(2)
      endif
c
1001  format(' Non-Solvent density map is map',i3,' with Input env',i3)
1002  format(' Non-Solvent density map is map',i3,' with no Input env')
1003  format(' Solvent content map is map',i3,' with Input env',i3)
1004  format(' Solvent content map is map',i3,' with no Input env')
c
c---- get output map
c
      if(.not.lass_get(ch_temp,'MAPO'))then
        ch_temp=G2CHRU('AV_GEN:OUTPUT_MAP',1,8,0)
        call lass_set(ch_temp,'MAPO')
      endif
      omap=ch_temp
c
c---- get output envelope
c
      if(.not.lass_get(ch_temp,'ENVO'))then
        ch_temp=G2CHRU('AV_GEN:OUTPUT_ENV',1,8,0)
        call lass_set(ch_temp,'ENVO')
      endif
      oenv=ch_temp
      o_check=oenv.ne.'OFF'
c
c---- confirm in writing
c
      if(verbose)then
        write(6,110) omap,oenv
        if(record)write(nrec,110) omap,oenv
      else
        write(6,111) omap,oenv
        if(record)write(nrec,111) omap,oenv
      endif
c
      if(verbose)then
        if(o_check) write(6,140)
        if(.not.o_check)write(6,150)
        if(record) then
          if(o_check) write(nrec,140)
          if(.not.o_check)write(nrec,150)
        endif
      endif
c
110   format(/' NON-SOLVENT DENSITY SCALING',/,
     &        ' ==============================',/,
     &        ' Output map:      ',a,/,
     &        ' Output envelope: ',a)
111   format(' FRACTIONAL SOLVENT CALCULATION',
     &  ' O/P map:',a,'  O/P env:',a)
140   format(' Output maps envelope filtered')
150   format(' Output maps not filtered')
c
c---- check status of output files
c
      call av_ass(map2,omap)
      if(.not.terse)then
        if(defined(map2))write(6,170)
        if(.not.defined(map2))write(6,180)
        if(record)then
          if(defined(map2))write(nrec,170)
          if(.not.defined(map2))write(nrec,180)
        endif
      endif
c
170   format(' Output Map exists: will overwrite')
180   format(' Make Output Map from i/p template')
c
c---- set up output map header information
c
      call copy_head(maps(1),map2)
c
c---- set up cell data
c
      call XRfrac(XRcell(1,map2),XRtr(1,1,map2),XRintr(1,1,map2),
     &            XRvol,.true.)
c
c---- set up symmetry data
c
      call XR_symsy(XRcell(1,map2),lgrp(map2),ics_local,
     &              n_xsym(map2),XR_sym(1,map2) )
c
c---- and integerized CSO's:
c
      call XR_isym(map2)
c
c---- set up bricking data
c
      call brickit(map2)
c
c---- allocate space for output map
c
      call alloc(map2)
c
c---- set up allocations for envelopes
c
      if(o_check) call av_ass(env2,oenv)
c
c---- set up fraction cut
c
      frac=1.0
c
10    nopt=keyword
     &  ('AV_GEN','FRAC!BACK!GO  ',0)
c
      if(nopt.eq.1)frac=realval('AV_GEN:FRAC',0.0,1.0,0)
      if(nopt.eq.2)quit=.true.
      if((nopt.ne.2).and.(nopt.ne.3))goto 10
c
      write(6,221)frac
      if(record)write(nrec,221)frac
221   format(' Scale up non-solvent density where fraction solvent < ',
     &  f8.6)
c
      return
      end
c


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE go_cut
c       =================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c--- local variables
c
        integer   map1,env1
        logical   i_check
        integer*2 irho
        external  gp3,gp3_a
        integer*2 gp3,gp3_a
c
        integer   i_hi,i_lo
        integer   iz,jz,iy,jy,ix,jx,jjz,jjy,jjx,ngadd
c
c
c
        call gcut_setup(map1,env1,i_check)
c
c
c
c
c---- the averaging process is very simple, we are driven by the
c---- output map, for each pixel in that map (accessed in a sequential
c---- fashion, but since map is bricked we move slowly through space)
c---- flatten sets all pixels above and below certain values to zero.
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
        call init_map(map1)
c
c---- initialize r-factors, corr coefs, etc
c
        ngadd=0
        i_hi=0
        i_lo=0
c
        call f_stat
c
        if(i_check)call init_map(env1)
c
c
        call init_map(map1)
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
        do iz=1,n_brick(3,map1)
c
          jz=(iz-1)*brick(3,map1)
c
          do iy=1,n_brick(2,map1)
c
            jy=(iy-1)*brick(2,map1)
c
            do ix=1,n_brick(1,map1)
c
              jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                if(jjy.le.nx(2,map1))then
c
                  do jjx=jx+1,jx+brick(1,map1)
c               
                  if(jjx.le.nx(1,map1))then
c
c---- check that this pixel is in required region of the output map
c
                    if(i_check)then
                      if(gp3(env1,jjx,jjy,jjz).eq.0)goto 200
                    endif
c
                        irho=gp3(map1,jjx,jjy,jjz)
c
                      ngadd=ngadd+1
c
                    call f_a_stat(irho)
c
c
c
c
                    if(gp3(map1,jjx,jjy,jjz).gt.ncut_u)then
                        call pp3(map1,jjx,jjy,jjz,ncut_u)
                        i_hi=i_hi+1
                    endif
c
c
                    if(gp3(map1,jjx,jjy,jjz).lt.ncut_l)then
                        call pp3(map1,jjx,jjy,jjz,ncut_l)
                        i_lo=i_lo+1
                    endif
c
c
c---- remem to put av_temp,etc in common
c
200               continue
                  endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
        enddo
c
c
c
c
c999     continue
        if(verbose)then
          write(6,998)ngadd
          if(record)write(nrec,998)ngadd
998       format(/,' CUTTING FINISHED,',/,
     &            ' ================',/,
     &             i9,' pixels were inspected')
          else
          write(6,997)ngadd
          if(record)write(nrec,997)ngadd
997       format(' CUTTING FINISHED,',i10,' pixels were inspected')
          endif
          write(6,996)i_hi,i_lo
          if(record)write(nrec,996)i_hi,i_lo
996       format(' Number of pixels set to upper & lower limits:',2i10)
c       
c
c
c---- print some stats
c
        call f_p_stat(.true.)
c
        return
        end
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION gp3(map,jj1,jj2,jj3)
C       =============================
      IMPLICIT NONE
C
C This version fiddled with by Robert 20/7/93 to try to speed up averaging
C of large maps. This routine assumes that the brick size is 8. This is only
C a first, minimal attempt at optimisation - this function should really be
C in-lined with the redundant calculations removed. The brick size allows
C more efficient bit twiddling, avoiding divisions

      INCLUDE 'average.fcm'
C       =====================

        integer*2 gp3
        integer map,jj1,jj2,jj3

        integer dx, dy, dz
        integer brkoff

C The first thing to do is calculate an offset within an 8x8x8 brick from
C jj1, jj2 and jj3

        dx = IAND(jj1-1,7)
        dy = IAND(jj2-1,7)
        dz = IAND(jj3-1,7)

C The brick offset is a bit more of a pain since it is not easily divided
C into powers of two, however we can avoid the division size subtracting
C the index mod 8 leaves the brick number times 8, so we divide the total
C brick size by 8 leaving 64...

        brkoff = 64 * (
     &          ( (jj3-1-dz) * n_brick(2,map)
     &          + (jj2-1-dy) ) * n_brick(1,map)
     &          + (jj1-1-dx) )

C Now just look up the answer

        gp3 = scratch(ns(map) + brkoff + 64*dz + 8*dy + dx + 1)

        return
        end

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION gp3_a(map,j3)
C       ======================
      IMPLICIT NONE

C This version fiddled with by Robert 22/7/93 the same way as GP3

      INCLUDE 'average.fcm'
C       =====================

        integer*2 gp3_a
        integer map,j3(3)

        integer dx, dy, dz, brkoff

        dx = IAND(j3(1)-1,7)
        dy = IAND(j3(2)-1,7)
        dz = IAND(j3(3)-1,7)

        brkoff = 64 * (
     &          ( (j3(3)-1-dz) * n_brick(2,map)
     &          + (j3(2)-1-dy) ) * n_brick(1,map)
     &          + (j3(1)-1-dx) )

        gp3_a = scratch(ns(map) + brkoff + 64*dz + 8*dy + dx + 1)

        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE grow
c     ===============
      IMPLICIT NONE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c This subroutine goes through the non-asymmetric portion of a labelled
c envelope (the bits marked with symmetry operators greater than one) and
c ensures that all pixels that will be required for the 11-point interpolation
c involved in the unpacking are included in the asymmetric portion of the
c envelope.
c
c It goes through the input map (envelope) and gets the pixel value. If the
c pixel value is greater than one, it uses the symmetry operator denoted by
c the pixel value to identify all the symmetry related pixels that will be
c involved in 11-point interpolation. The value in these pixels is set to one
c regardless of their original content.
c
c JMD 9/19/2001
c
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- local variables:
c
      integer     map1,iindex(3)
      real        x(3),x1(3),aindex(3)
      integer*2   i_env
      external    gp3,gp3_a
      integer*2   gp3,gp3_a,one
c
      real        uu(3),sum_r
      integer     n2,nstep,j
      integer     iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
      integer     ind(3),off(3)
c
c--- Initialisations
c
      n2=2
      nstep=1
      one=1
c
c--- Get map
c
      call grow_setup(map1)
c
c--- Set up pointers
c
      call init_map(map1)
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
      do iz=1,n_brick(3,map1)
c
        jz=(iz-1)*brick(3,map1)
c
        do iy=1,n_brick(2,map1)
c
          jy=(iy-1)*brick(2,map1)
c
          do ix=1,n_brick(1,map1)
c
            jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                  if(jjy.le.nx(2,map1))then
c
                    do jjx=jx+1,jx+brick(1,map1)
c               
                      if(jjx.le.nx(1,map1))then
c
c
      i_env=gp3(map1,jjx,jjy,jjz)
      if(i_env.le.1)goto 200
c
      iindex(1)=jjx
      iindex(2)=jjy
      iindex(3)=jjz
      call j3_to_f(iindex,x1,map1)
      call f_to_o(x1,x,map1)
c
      call spin_inv(x,x1,ops_inv(1,1,i_env),vecs(1,i_env))
      call o_to_j3(x1,aindex,iindex,map1)
c
c aindex contains spun point in real numbers 
c index contains spun point in integers, ie closest pixel
c difference between the two tells which way to go for interpolation
c
      do j=1,3
        uu(j)=aindex(j)-iindex(j)
        if(uu(j).lt.0.0)then
          uu(j)=abs(uu(j))
          ind(j)=-1
        else
          ind(j)=1
        endif
      enddo
c
      sum_r=uu(1)+uu(2)+uu(3)
      if(sum_r.lt.0.1)then
        call pp3_a(map1,iindex,1)
        goto 200
      endif
c
      do j=1,3
        off(j)=ind(j)+iindex(j)
        if((off(j).gt.nx(j,map1)).or.(off(j).lt.1))then
          call pp3_a(map1,iindex,one)
          goto 200
        endif
      enddo
c
      call pp3(map1,iindex(1),iindex(2),iindex(3),one)
      call pp3(map1,off(1),iindex(2),iindex(3),one)
      call pp3(map1,iindex(1),off(2),iindex(3),one)
      call pp3(map1,iindex(1),iindex(2),off(3),one)
      call pp3(map1,off(1),off(2),iindex(3),one)
      call pp3(map1,iindex(1),off(2),off(3),one)
      call pp3(map1,off(1),iindex(2),off(3),one)
      call pp3(map1,off(1),off(2),off(3),one)
c
      do j=1,3
        off(j)=-ind(j)+iindex(j)
        if((off(j).gt.nx(j,map1)).or.(off(j).lt.1))then
          goto 200
        endif
      enddo
c
      call pp3(map1,off(1),iindex(2),iindex(3),one)
      call pp3(map1,iindex(1),off(2),iindex(3),one)
      call pp3(map1,iindex(1),iindex(2),off(3),one)
c
200   continue
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c
c--- All over
c
      return
      end

c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE grow_setup(map1)
c     ===========================
      IMPLICIT NONE
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- local variables
c
      integer     map1
      character*8 G2CHRU,ch_temp,imap
      logical     lass_get
c
c
c
      if(.not.lass_get(ch_temp,'MAPI'))then
        ch_temp=G2CHRU('AV_GROW:INPUT_MAP',1,8,0)
        call lass_set(ch_temp,'MAPI')
      endif
      imap=ch_temp
c
      if(verbose)then
        write(6,110) imap
        if(record)write(nrec,110) imap
110     format(/' GROW LABELLED ENVELOPE',/,
     &          ' ======================'
     &         /' Input map  :  ',a)
      else
        write(6,112) imap
        if(record)write(nrec,112) imap
112     format(' GROW LABELLED ENVELOPE  I/P map:',a)
      endif
c
c----get slot numbers for all maps
c
      call av_ass(map1,imap)
c
c----map 1 should exist by now !
c
      if(.not.defined(map1)) then
        write(6,140)
        if(record)write(nrec,140)
140     format(' %GROW_SETUP-ERR: Input map empty')
        return
      endif
c
      return
      end

C=============================================================================

      SUBROUTINE GSORTH(MATRIX, ERRCOD)
C     =================================
      IMPLICIT NONE

C Perform a Gram-Schmidt orthogonalisaton by columns. If the input matrix is
C singular then it returns an error code. Code effectively stolen from the
C routine by Morten Kjeldgaard as used in the o2mol program of Janet Smith.

C MATRIX (DP) The 3x3 rotation matrix
C ERRCOD (I)  Returned as -1 is the matrix is singular, 0 otherwise

      DOUBLE PRECISION MATRIX(3,3)
      INTEGER ERRCOD

      INTEGER I, J, K
      DOUBLE PRECISION DOT, DNORM
      INTEGER N
      PARAMETER (N=3)

      ERRCOD = 0
      DO I = 1, N-1

        DNORM = 0.0
        DO J = 1, N
          DNORM = DNORM + MATRIX(I,J)**2
        ENDDO

        IF (DNORM.LT.1.0E-12) THEN
          ERRCOD = -1
          RETURN
        ENDIF

        DO J = I+1, N
          DOT = 0.0
          DO K = 1, N
            DOT = DOT + MATRIX(I,K) * MATRIX(J,K)
          ENDDO
          DO K = 1, N
            MATRIX(J,K) = MATRIX(J,K) - MATRIX(I,K) * DOT / DNORM
          ENDDO
        ENDDO

      ENDDO

C Normalise the matrix by columns

      DO J = 1, N
        DOT = 0.0
        DO I = 1, N
          DOT = DOT + MATRIX(J,I)**2
        ENDDO
        DOT = SQRT(DOT)
        DO I = 1, N
          MATRIX(J,I) = MATRIX(J,I) / DOT
        ENDDO
      ENDDO

      RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE gulp_map(lun,line,bit,lsize,nmap,speed)
C     ======================================================
      IMPLICIT NONE

C Read formatted section from unit lun into slot nmap
c
      INCLUDE 'average.fcm'
c     =====================
c
      integer order(3), istore,lsize,i,lun,j,k,nmap,isec
      integer*2 line(lsize)
      byte bit(lsize)
      character*4 speed
c
      istore=0
      do i=1,nx(iuvw(3,nmap),nmap)
        order(iuvw(3,nmap))=i
        if(speed.eq.'fast')then
          read(lun)isec
        else
          read(lun,6001)isec
        endif
6001    format(/7x,i8/)
        if(.not.terse)then
          if( i.eq.2.and.i.ne.nx(iuvw(3,nmap),nmap) )then
            write(6,*)'.....................'
            if(record)write(nrec,*)'.....................'
          endif
          if( i.le.1.or.i.eq.nx(iuvw(3,nmap),nmap) )then
            write(6,*)'Section',isec,' read.'
            if(record)write(nrec,*)'Section',isec,' read.'
          endif
        endif
c
        do j=1,nx(iuvw(2,nmap),nmap)
          order(iuvw(2,nmap))=j
          if(speed.eq.'slow')then
            read(lun,6010)(line(k),k=1,lsize)
          else
            call pmud(lun,bit,lsize)
          endif
6010      format(20i4)
c    
c       Store density recentered to zero and NOT scaled up by two
c       to fit in +/- 32000
c       (stored as +9999 to -999 formatted)
c
          if(speed.eq.'slow')then
            do k=1,lsize
              line(k)= ( float(line(k))-moffset(nmap) )
              order(iuvw(1,nmap))=k
               call pp3(nmap,order(1),order(2),order(3),line(k))
            enddo
            istore=istore+lsize
          else
            do k=1,lsize
              line(k)= bit(k)
              order(iuvw(1,nmap))=k
              call pp3(nmap,order(1),order(2),order(3),line(k))
            enddo
            istore=istore+lsize
          endif
c
        enddo
c
      enddo
c
      write(6,6020)istore
      if(record)write(nrec,6020)istore
6020  format(' Input completed,',i10,' pixels read and stored')
c
      return
      end
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE gulp_mapb(lun,rho,isize,line,lsize,nmap,scale_r)
C     =========================================================
      IMPLICIT NONE

C Read formatted section from unit lun into slot nmap
c
      INCLUDE 'average.fcm'
c     =====================
c
      integer   order(3),istore,lsize,i,lun,j,k,nmap,isec,isize,ier,
     &          inrho
      integer*2 line(lsize)
      real      rho(isize),scale_r
c
      istore=0
      isec=0
      do i=1,nx(iuvw(3,nmap),nmap)
        isec=i
        order(iuvw(3,nmap))=i
        call mgulp(lun,rho,ier)
        inrho=1
        if(.not.terse)then
          if( i.eq.2.and.i.ne.nx(iuvw(3,nmap),nmap) )then
            write(6,*)'.....................'
            if(record)write(nrec,*)'.....................'
          endif
          if( i.le.1.or.i.eq.nx(iuvw(3,nmap),nmap) )then
            write(6,*)'Section #',isec,' read.'
            if(record)write(nrec,*)'Section #',isec,' read.'
          endif
        endif
c
        do j=1,nx(iuvw(2,nmap),nmap)
          order(iuvw(2,nmap))=j
c    
c       Store density recentered to zero and NOT scaled up by two
c       to fit in +/- 32000
c       (stored as +9999 to -999 formatted)
          do k=1,lsize
            line(k)= nint(rho(inrho)*scale_r)
            inrho=inrho+1
            order(iuvw(1,nmap))=k
            call pp3(nmap,order(1),order(2),order(3),line(k))
          enddo
          istore=istore+lsize
c
        enddo
c
      enddo
c
      write(6,6020)istore
      if(record)write(nrec,6020)istore
6020  format(' Input completed,',i10,' pixels read and stored')
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE header(flag)
c     =======================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      character*20 date
      character*18 time
      integer      ilen
      logical flag
c
c
      call P2TIME(time,ilen)
      call P2DATE(date,ilen)
      if(flag)   write(6,10),time,date,maxsym,mmaps
      if(record) write(nrec,10),time,date,maxsym,mmaps
10    format(' =====================================================',/
     &       ' General Averaging Program,    V6.5.0,     02-Apr-2010',/,
     &       ' =====================================================',/
     &       ' Run at    ',a,'  on  ',a,/
     &       ' Max NCS ',i5,'               Max number of maps ',i5,/
     &       ' =====================================================')
      return
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE hist_set(map1,map2)
c     ==============================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c     =====================
c
c---- local variables:
c
      integer         map1,map2,ics_local,j
      real            xrvol
      character*8     imap,oenv,G2CHRU,ch_temp
      logical         lass_get
c
c---- get input map name and number
c
      if(.not.lass_get(ch_temp,'MAPI'))then
        ch_temp=G2CHRU('AV_HIST:INPUT_MAP',1,8,0)
        call lass_set(ch_temp,'MAPI')
      endif
      imap=ch_temp
c
c---- get output envelope name and number
c
      if(.not.lass_get(ch_temp,'ENVO'))then
        ch_temp=G2CHRU('AV_HIST:WANG_ENV',1,8,0)
        call lass_set(ch_temp,'ENVO')
      endif
      oenv=ch_temp
c
c---- for the record
c
      if(verbose)then
        write(6,110) imap,oenv
        if(record)write(nrec,110) imap,oenv
110     format(/' ENVELOPE FROM HISTOGRAM',/,
     &          ' ======================='
     &         /' Input map:       ',a,
     &         /' Output envelope: ',a)
      else
        write(6,111) imap,oenv
        if(record)write(nrec,111) imap,oenv
111     format(' ENVELOPE FROM HISTOGRAM',
     &    '  I/P map:',a,'  O/P env:',a)
      endif
c
c---- get slot number for map1, and drop out if it doesn't exist
c
      call av_ass(map1,imap)
      if(.not.defined(map1)) then
        write(6,120)
        if(record)write(nrec,120)
120     format(' %HIST_SET-ERR: Input map empty')
        return
      endif
c
c---- get slot number for map2, and overwrite if it does exist
c
      call av_ass(map2,oenv)
      if(.not.terse)then
        if(defined(map2))write(6,*)'output env exists: will overwrite'
        if(.not.defined(map2))write(6,*)'make output from i/p template'
c
        if(record)then
          if(defined(map2))write(nrec,*)
     &      'output env exists: will overwrite'
          if(.not.defined(map2))write(nrec,*)
     &      'make output from i/p template'
        endif
      endif
c
c---- set up output map header information
c
      if(.not.defined(map2))then
        call copy_head(map1,map2)
c
c---- set up cell data
c
        call XRfrac(XRcell(1,map2),XRtr(1,1,map2),XRintr(1,1,map2),
     &    XRvol,.true.)
c
c---- set up symmetry data
c
        call XR_symsy(XRcell(1,map2),lgrp(map2),ics_local,
     &    n_xsym(map2),XR_sym(1,map2) )
c
c---- and integerized CSO's:
c
        call XR_isym(map2)
c
c---- also set up pointers indicating which symm op to try first
c
        do j=1,maxsym
          last_sym(j,map2)=1
        enddo
c
c---- set up bricking data
c
        call brickit(map2)
c
c---- allocate space for output map
c
        call alloc(map2)
c
      endif
c
c---- reset scale and offset for an envelope (in case we write it out)
c
      if(.not.terse)then
        write(6,*)
     &    'NB scale and offset reset to 1 and 0 for this envelope'
        if(record)write(nrec,*)
     &    'NB scale and offset reset to 1 and 0 for this envelope'
      endif
c
      mscale(map2)=1.0
      moffset(map2)=0.0
      rholim(1,map2)=0.0
      rholim(2,map2)=1.0
      rholim(3,map2)=0.0
      rholim(4,map2)=0.0
c
c---- and return
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE histo
c     ================
      IMPLICIT NONE
c
c---- do histogram analysis of map to produce an envelope from a map
c
      INCLUDE 'average.fcm'
c      =====================
c
c---- local variables:
       integer     map1,map2
       integer*2   irho,ieset,solv,notsolv,ir
       external    gp3,gp3_a
       integer*2   gp3,gp3_a
c
       integer     hist(-32000:32000),chist(-32001:32000),ngadd
       integer     bin(0:20),nspix,j,jj,jj_old,ncut,n_set
       real        realval,fra_sol,shell
       integer     pmin,pmax,iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
c
c---- set local constants:
c
       solv=0
       notsolv=1
c
c---- initialise counters
c
       chist(-32001)=0
       do j=-32000,32000
         hist(j)=0
         chist(j)=0
       enddo
       do j=0,20
         bin(j)=0
       enddo
c
c---- get map numbers
c
       call hist_set(map1,map2)
c
c---- give some info on output
c
       write(6,1000)solv,notsolv
       if(record)write(nrec,1000)solv,notsolv
1000   format(' Solvent pixels will be set to:',i5,'  other pixels to:'
     &         ,i5)
c
c---- get target fraction solvent
c
       fra_sol=realval('Enter solvent fraction:',0.0,1.0,0)
       if(record)write(nrec,1005)fra_sol
1005   format(' Fractional solvent content will be set to:',f8.3)
c
c
c---- now loop thorugh input map once to calculate histogram
c---- initialize pointers for map1
c
       call init_map(map1)
c
c---- initialize r-factors, corr coefs, etc
c
       ngadd=0
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
       do iz=1,n_brick(3,map1)
c
         jz=(iz-1)*brick(3,map1)
c
         do iy=1,n_brick(2,map1)
c
           jy=(iy-1)*brick(2,map1)
c
           do ix=1,n_brick(1,map1)
c
             jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
             do jjz=jz+1,jz+brick(3,map1)
c
             if(jjz.le.nx(3,map1))then
c
               do jjy=jy+1,jy+brick(2,map1)
c
               if(jjy.le.nx(2,map1))then
c
                 do jjx=jx+1,jx+brick(1,map1)
c                    
                 if(jjx.le.nx(1,map1))then
c
                  ir=gp3(map1,jjx,jjy,jjz)
                  if(ir.gt.32000)ir=32000
                  if(ir.lt.-32000)ir=-32000
c
                  hist(ir)=hist(ir)+1
                  ngadd=ngadd+1
c
                  endif
                 enddo
               endif
               enddo
             endif
             enddo
           enddo
         enddo
       enddo
c
c--- set up a running total and find first and last occupied ranges
c
       pmin=32000
       pmax=-32000
       do j=-32000,32000
         chist(j)=chist(j-1)+hist(j)
         if(hist(j).gt.0)pmax=j
         if(hist(j).gt.0.and.j.lt.pmin)pmin=j
       enddo
c
c--- calculate number of pixels to be set as solvent
c
       nspix=ngadd*fra_sol
       write(6,1010)nspix
       if(record)write(nrec,1010)nspix
1010   format(' Target for number of solvent pixels',i10)
c
c--- find cutoff
c
       do j=pmin,pmax
         if (chist(j).le.nspix)ncut=j
       enddo
       write(6,1020)ncut
       if(record)write(nrec,1020)ncut
1020   format(' The cut-off is defined at rho =',i6)
c
c--- Show the histogram
c
       if(.not.terse)then
c
c---- Find bin size
c
         shell=float(pmax+1-pmin)/20.0
         write(6,*) 'Thickness of ranges for analysis:', shell
         if(record)write(nrec,*)
     &              'Thickness of ranges for analysis:', shell
c
c---- Get start of each bin - jj runs from 0 to 20
c
         jj_old=-1
         do j=pmin,pmax+1
           jj=int((j-pmin)/shell)
           if(jj.ne.jj_old)then
             bin(jj)=j
             jj_old=jj
           endif
         enddo
c
c---- Give statistics in bins - histogram!
c---- previous max_rho = bin(j)-1
c---- min_rho = bin(j)
c---- max_rho = bin(j+1) - 1
c---- this_bin = chist(max_rho) - chist(previous max_rho)
c---- cumulative = chist(max_rho)
c
         write(6,1100)
         if(record)write(nrec,1100)
         do j=0,19
           write(6,1110)
     &       bin(j),chist(bin(j+1)-1)-chist(bin(j)-1),chist(bin(j+1)-1)
           if(record)write(nrec,1110)
     &       bin(j),chist(bin(j+1)-1)-chist(bin(j)-1),chist(bin(j+1)-1)
          enddo
       endif
1100    format(' ====================================================')
1110   format(' Min rho',i6,'    # pixels',i8,'    sum #',i10)
c
c--- go through map1 again, this time creating the envelope as map2
c--- initialize pointers for map1
c
       call init_map(map1)
c
c--- off to work...
c--- =============
c
c--- loop over bricks
c
       do iz=1,n_brick(3,map1)
c
         jz=(iz-1)*brick(3,map1)
c
         do iy=1,n_brick(2,map1)
c
           jy=(iy-1)*brick(2,map1)
c
           do ix=1,n_brick(1,map1)
c
             jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
             do jjz=jz+1,jz+brick(3,map1)
c
             if(jjz.le.nx(3,map1))then
c
               do jjy=jy+1,jy+brick(2,map1)
c
               if(jjy.le.nx(2,map1))then
c
                 do jjx=jx+1,jx+brick(1,map1)
c                    
                 if(jjx.le.nx(1,map1))then
c
                  irho= gp3(map1,jjx,jjy,jjz)
c
                  ieset=solv
c
                  if(irho.gt.ncut) then
                    ieset=notsolv
                    n_set=n_set+1
                  endif
                  call pp3(map2,jjx,jjy,jjz,ieset)
c
                  endif
                 enddo
               endif
               enddo
             endif
             enddo
           enddo
         enddo
       enddo
c
c--- and write out some stats
c
       write(6,1200)ngadd-n_set,float(ngadd-n_set)/float(ngadd),n_set
       if(record)write(nrec,1200)
     &              ngadd-n_set,float(ngadd-n_set)/float(ngadd),n_set
1200   format(' WANG TYPE ENVELOPE CALCULATED,',i10,
     &        ' pixels were set to solvent'/,
     &        ' (fractional solvent content = ',f5.3,', ',i10,
     &         ' non-solvent pixels)')
c
c--- and its all over
c
       return
       end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      FUNCTION i_box(iindex,map)
c     =========================
      IMPLICIT NONE
c
c  determine whether pixel INDEX is within the available portion of map
c  MAP.
c  I_BOX returns .true. if it is.
c
c
c
      INCLUDE 'average.fcm'
c       =====================
c
      integer iindex(3),map
      logical i_box
c
      i_box=.true.
c
      if(iindex(1).lt.1 .or. iindex(2).lt.1 .or. iindex(3).lt.1
     &    ) then
        i_box=.false.
        return
      endif
c
      if( iindex(1).gt.nx(1,map) .or. iindex(2).gt.nx(2,map) .or.
     &    iindex(3).gt.nx(3,map) ) i_box=.false.
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE i_stat
c       =================
      IMPLICIT NONE
c
c--- These subrouines accumulate statistics on the averaging process.
c--- Note that the statistics may be more time consuming than the averaging
c--- and for this reason two types of correlation coefficients are provided
c--- (i) a full pair-wise analysis  -- cost proportional to  (# NCS ops)**2
c--- (ii) an analysis w.r.t. to mean -- "    "           "   (# NCS ops)**1
c
      INCLUDE 'average.fcm'
c       =====================
c
       integer n_uni,n_tot,j,k
       real*8  r_max,r_min,r_top,r_bot,r1_top
       real*8  r1_bot,cc_t,cc_b,rms,rms_del
       real*8  r_r,del,r_r_s
       logical no_printout
c
       save cc_t,cc_b,r_top,r_bot,r1_top,r1_bot,rms_del,n_tot,
     &      rms,n_uni,r_min,r_max
c
c--- Initialise counters
c
       cc=zero
       r=zero
       r1=zero
       r_max=-32000.0D0
       r_min=32000.0D0
       r_top=zero
       r_bot=zero
       r1_top=zero
       r1_bot=zero
       cc_t=zero
       cc_b=zero
       rms=zero
       rms_del=zero
       n_uni=0
       n_tot=0
       do j=1,maxsym
         h_temp(j)=0
       enddo
       return
c
c--- Accumulate statistics
c
       entry a_stat
c      ============
c
c--- One observation, check rho min/max but no averaging:
c
       if(in_temp.eq.1) then
         h_temp(in_temp)=h_temp(in_temp)+1
c
         if(sum_temp.gt.r_max)r_max=sum_temp
         if(sum_temp.lt.r_min)r_min=sum_temp
c
         return
       endif
c
c--- No observations, set rho to zero:
c
       if(in_temp.le.0)then
         sum_temp=zero
       else
c
c--- More than one obsevation - average and record stats:
c
c--- Accumulate histogram of multiplicity
c
         h_temp(in_temp)=h_temp(in_temp)+1
c
c--- Average
c
         sum_temp=sum_temp/float(in_temp)
c
c--- Check min/max on averaged rho:
c
         if(sum_temp.gt.r_max) r_max=sum_temp
         if(sum_temp.lt.r_min) r_min=sum_temp
c
c--- Running totals of input (n_tot) and output (n_uni) pixels
c
         n_uni=n_uni+1
         n_tot=n_tot+in_temp
c
c--- Accumulate <rho>**2 for calculation of RMS deviation in output map
c
         rms=rms+sum_temp**2
c
c--- Initialise cc_top, cc_bot and rms for internal loop over obs
c
         cc_top=zero
         cc_bot=zero
         rms_delta=zero
c
c--- Loop over observations for stats:
c
         do j=1,in_temp
c
c---- Accumulate r_top and r_bot for R factor calculation
c
           r_r=float(ir_temp(j))
           del=abs(sum_temp-r_r)
           r_top=r_top+del
           r_bot=r_bot+abs(r_r)
c
c---- Accumulate r1_top and r1_bot for R where <rho> > rho max/4
c---- NB This is dodgy as rho max isn't definitely known till the end
c
           if(sum_temp.gt.r_max/4.0)then
             r1_top=r1_top+del
             r1_bot=r1_bot+abs(r_r)
           endif
c
c---- Accumulate rms deviation during averaging
c
           rms_delta=rms_delta+del**2
c
c---- Accumulate cc_top and cc_bot for correlation coefficient calculation
c---- Choice of two correlation coefficients:
c----  istat.gt.0 => wrt mean
c----  istat.le.0 => pairwise
c
           if(istat.gt.0)then
c
c---- Wrt mean:
c
             cc_top=cc_top+r_r*sum_temp
             cc_bot=cc_bot+r_r**2
c
           else
c
c---- Pairwise:
c
             if(j.lt.in_temp)then
               r_r_s=r_r**2
               do k=j+1,in_temp
                 cc_top=cc_top+r_r*float(ir_temp(k))
                 cc_bot=cc_bot+r_r_s
               enddo
             endif
c
           endif
c
         enddo
c
c--- Accumulate cc_t and cc_b for running total over whole map
c
         cc_t = cc_t + cc_top
         cc_b = cc_b + cc_bot
c
c--- Accumulate rmsd during averaging for running total over whole map
c
         rms_del = rms_del + rms_delta
c
       endif
c
       return
c
c--- Calculate statistics - g_stat calculates, p_stat calculates and
c--- prints
c
       entry g_stat
c       ============
c
       no_printout=.true.
       goto 10
c
       entry p_stat
c       ============
c
       no_printout=.false.
c
10     continue
c
       cc=0.0
       if(cc_b.gt.zero)  cc=cc_t/cc_b
c
       r=0.0
       if(r_bot.ne.zero) r=r_top/r_bot
c
       r1=0.0
       if(r1_bot.ne.zero) r1=r1_top/r1_bot
       if(n_tot.gt.0.and.rms_del.gt.zero) rms_del=sqrt(rms_del/n_tot)
       if(n_uni.gt.0.and.rms.gt.zero)     rms    =sqrt(rms    /n_uni)
c
       if(no_printout) return
c
c--- Display results
c
       if(.not.terse)then
         write(6,100)
         if(record)write(nrec,100)
       endif
c
       if(istat.eq.0)write(6,110)
     &      n_tot,n_uni,r_max,r_min,rms,rms_del,r,r1,cc
       if(istat.gt.0)write(6,120)
     &      n_tot,n_uni,r_max,r_min,rms,rms_del,r,r1,cc
c
       if(record)then
         if(istat.eq.0)write(nrec,110)
     &      n_tot,n_uni,r_max,r_min,rms,rms_del,r,r1,cc
         if(istat.gt.0)write(nrec,120)
     &      n_tot,n_uni,r_max,r_min,rms,rms_del,r,r1,cc
       endif
c
100    format(' AVERAGING STATISTICS')
110    format(
     &  ' Total # of observations:',i10,';',i9,
     &  ' pixels had 2 or more observations',/
     &  ' Max density after averaging:',f8.1,
     &  '   Min density after averaging:',f8.1,/
     &  ' Rms density after averaging:',f8.1,
     &  '   Rms dev.   during averaging:',f8.1,/
     &  ' Real space R-fact (all rho):',f8.3,
     &  '   Real space R (rho>rhomax/4):',f8.3,/
     &  ' Corr.coeff. between all pairs:',f8.3)
120    format(
     &  ' Total # of observations:',i10,';',i9,
     &  ' pixels had 2 or more observations',/
     &  ' Max density after averaging:',f8.1,
     &  '   Min density after averaging:',f8.1,/
     &  ' Rms density after averaging:',f8.1,
     &  '   Rms dev.   during averaging:',f8.1,/
     &  ' Real space R-fact (all rho):',f8.3,
     &  '   Real space R (rho>rhomax/4):',f8.3,/
     &  ' Corr coeff (each  cf. mean):',f8.3)
c
       if(.not.terse)then
         write(6,130)
         if(record)write(nrec,130)
       endif
       do j=1,maxsym
         if(h_temp(j).gt.0)then
           write(6,140)h_temp(j),j
           if(record)write(nrec,140)h_temp(j),j
           h_temp(j)=0
         endif
       enddo
c
130    format(' HISTOGRAM OF AVERAGING MULTIPLICITY')
140    format(i9,' pixels were reconstructed from',
     &        i4,' observations')
c
       return
       end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      SUBROUTINE idealise_mats
c       ========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c subroutine to take non-ideal proper operators and 
c force them to be ideal.  Most subroutines written
c by R.E. 
c
        integer intval,nncs,j,x,y,ierr,nfold,nsym1,nsym2
        DOUBLE PRECISION dbROT(3,3,maxsym)
        DOUBLE PRECISION dbTRAN(3,maxsym)


        nncs=0
        ierr=0

        if (nsym.lt.1) then
          write(6,10)
            if(record) write(nrec,10)
          return
        endif
10      format('AV_IDEALIZE_ERR: No Symmetry Operators defined')
c
        nfold=intval('AV_IDEALIZE_SYM:ENTER_SYMMETRY',1,nsym,0)
c
        nsym1=intval('AV_IDEALIZE_SYM:ENTER_1_NUMBER',1,nsym,0)
        nsym2=intval('AV_IDEALIZE_SYM:ENTER_2_NUMBER',nsym1,nsym,0)
c
        if(nfold.ne.(nsym2-nsym1+1)) then
          write(6,15)
            if(record) write(nrec,15)
          return
        endif
15      format('AV_IDEALIZE_ERR: Symmetry and operators unmatched')


          write(6,20) nsym1,nsym2
          if(record) write(nrec,20) nsym1,nsym2
20      format(' Symmetry operators to idealize',i4,' to ',i4)
c
c
        if(verbose)then
          do j=nsym1,nsym2
            call prtsym(6,j)
          enddo
        endif
c
c
c-- at mo copy real ops,ops_inv,vecs into double precision 
c-- do we want to make these double precision ??
c
c-- would use subroutine copy by default but it is not double precision
c
c-- also set up the n-foldness of the symmetry, jj

        do j=nsym1,nsym2
        nncs=nncs+1
        do x=1,3
          dbTRAN(x,nncs) = vecs(x,j)
        do y=1,3
          dbrot(x,y,nncs) = ops(x,y,j)
        enddo
        enddo
        enddo

        call nfoldncs(nncs, dbrot,dbtran, dbrot, dbtran, ierr)
c
        if(ierr.gt.1) then
          write(6,30)
          if(record) write(nrec,30)
30        format(' AV_IDEALIZE_ERR: Rotation matrix is singular')
          return
        endif
c
        if(ierr.eq.1) then
          write(6,40)
          if(record) write(nrec,40)
40        format(' AV_IDEALIZE_ERR: First matrix must be identity')
          return
        endif

        if(ierr.eq.-2) then
          write(6,40)
          if(record) write(nrec,50)
50        format(' AV_IDEALIZE_ERR: Problem sorting rotation angles')
          return
        endif

c
c Now write back new operators to ops defined as reals
c
        do j=nsym1,nsym2
        do x=1,3
          vecs(x,j) =  dbTRAN(x,j-nsym1+1)
        do y=1,3
          ops(x,y,j) = dbrot(y,x,j-nsym1+1)
        enddo
        enddo
        enddo

c
c
          write(6,60)
          if(record) write(nrec,60)
60        format(' AV_IDEALIZE: Symmetry operators idealised')

        if(verbose)then
          do j=nsym1,nsym2
            call prtsym(6,j)
          enddo
        endif
c


        RETURN
        END


c=============================================================================
c
c all of rob esnoufs code above here
C **********************************
c
C=============================================================================
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      FUNCTION iinout(map,ind1,ind2,ind3,new_index,check,jncs)
c     ========================================================
      IMPLICIT NONE
c
C***********************************************************************
c  This funtion returns value .true. if the point INDEX can be folded
c  into the available portion of the map.
c  NEW_INDEX contains the folded coordinates.
c  CHECK should be sent as .true. in order to search all the
c    crystallographic symmetry operators (CSO) to find a folding strategy.
c  JNCS is the number of the current NCS operator (NCSO) (to save time we
c    keep track for each NCSO which is the currently succesful
c    CSO for folding).
c
c
      INCLUDE 'average.fcm'
c       =====================
c
      integer  iindex(3),new_index(3),now_cryst
      logical  iinout,check,i_box,spin_in
      external i_box,spin_in
      integer ind1,ind2,ind3,j,map,jncs,j3
c
      iindex(1)=ind1
      iindex(2)=ind2
      iindex(3)=ind3
c
      do j=1,3
        new_index(j)=iindex(j)
      enddo
c
c---- if not folding check trivial:
c
      if(.not.check) then
        iinout=i_box(new_index,map)
        return
      endif
c
c---- try folding.
c
c---- first of all check the last successful CSO for this NCSO 
c
      now_cryst=last_sym(jncs,map)
      iinout=spin_in(new_index,map,now_cryst)
c
c---- if it works we are done:
c
      if(iinout) return
c
c----no, march through all the opeators except that last one.
c
      do j=1,n_xsym(map)
c
        if(j.ne.now_cryst) then
c
          do j3=1,3
            new_index(j3)=iindex(j3)
          enddo
c
          iinout=spin_in(new_index,map,j)
c
c---- any good ?
c
          if(iinout) then
c
c---- ok, store this as best bet for next pixel and return
c
            last_sym(jncs,map)=j
            return
          endif
c
        endif
c
      enddo
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE ilist_sym
c       ====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        integer intval,nsym1,nsym2,j
        if (nsym.lt.1) then
          write(6,10)
            if(record) write(nrec,10)
          return
        endif
10      format(' No Symmetry Operators defined')
c
        nsym1=intval('AV_LIST_INV_SYM:ENTER_1_NUMBER',1,nsym,0)
        nsym2=intval('AV_LIST_INV_SYM:ENTER_2_NUMBER',nsym1,nsym,0)
        do j=nsym1,nsym2
        call iprtsym(6,j)
        enddo
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c JMD!PORT: This function is not used?
c
c      FUNCTION index(j,n)
c       ===================
c      IMPLICIT NONE
c        integer index,j,n
c
c        index=(j-1)/n+1
c        return
c        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE init
c       ===============
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
        integer j
c
        ncut_l=0
        ncut_u=32000
        nsym=0
        mset=1
        nmaps=0
        nxtpix=1
        do j=1,mmaps
          xr_con(j)=.false.
        enddo
        cg_check=.false.
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE init_map(map)
c       ========================
      IMPLICIT NONE
c
c       initialize pointers to start of map for sequential access
c
c
      INCLUDE 'average.fcm'
c       =====================
        integer map
c
        pix_now(map)=0
        word_now(map)=ns(map)-1
        return
        end
c   )
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE INVERS(A,B,DET)
c     ==========================
      IMPLICIT NONE
C***** A    - 3 X 3 MATRIX                                       (GIVEN)
C***** B    - 3 X 3 MATRIX IS CALCULATED INVERSE MATRIX OF A    (RESULT)
C***** DET  - CALCULATED DETERMINANT OF GIVEN MATRIX A          (RESULT)
      INTEGER   I,J,K
        real a(3,3),b(3,3), det


      DO 1 I=1,3
      J=I+1
      IF (J.GT.3)J=J-3
      K=I+2
      IF (K.GT.3)K=K-3
      B(I,1)=SNGL(DPROD(A(2,J),A(3,K))-DPROD(A(2,K),A(3,J)))
      B(I,2)=SNGL(DPROD(A(3,J),A(1,K))-DPROD(A(3,K),A(1,J)))
1     B(I,3)=SNGL(DPROD(A(1,J),A(2,K))-DPROD(A(1,K),A(2,J)))
      DET=SNGL(DPROD(A(1,1),B(1,1))+DPROD(A(1,2),B(2,1))
     &          +DPROD(A(1,3),B(3,1)))
      IF (DET.EQ.(0.0))GO TO 3
      DO 2 I=1,3
      DO 2 J=1,3
2     B(I,J)=B(I,J)/DET
3       RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE INVERS_D(A,B,DET)
c     ==========================
      IMPLICIT NONE
C***** A    - 3 X 3 MATRIX                                       (GIVEN)
C***** B    - 3 X 3 MATRIX IS CALCULATED INVERSE MATRIX OF A    (RESULT)
C***** DET  - CALCULATED DETERMINANT OF GIVEN MATRIX A          (RESULT)
      INTEGER   I,J,K
        double precision a(3,3),b(3,3), det


      DO 1 I=1,3
      J=I+1
      IF (J.GT.3)J=J-3
      K=I+2
      IF (K.GT.3)K=K-3
      B(I,1)=A(2,J)*A(3,K)-A(2,K)*A(3,J)
      B(I,2)=A(3,J)*A(1,K)-A(3,K)*A(1,J)
1     B(I,3)=A(1,J)*A(2,K)-A(1,K)*A(2,J)
      DET=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      IF (DET.EQ.(0.0))GO TO 3
      DO 2 I=1,3
      DO 2 J=1,3
2     B(I,J)=B(I,J)/DET
3       RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE INVRT(R1,RINV1,NODER)
c     ================================
      IMPLICIT NONE
c
c
        integer noder,i,j,iinv,jinv,i1,j1
        real r,r1,rwork,rinv,det,det4,rinv1
        external det
      DIMENSION R1(NODER,NODER),RWORK(3,3),RINV1(NODER,NODER)
      DIMENSION R(4,4),RINV(4,4)
c
      DO 1 I=1,4
      R(I,4)=0.
1     R(4,I)=0.
      R(4,4)=1.
      DO 2 I=1,NODER
      DO 2 J=1,NODER
2     R(J,I)=R1(J,I)
C                                                                               
C                                                                               
      DO 200 IINV=1,4
      DO 200 JINV=1,4
      I1=0
      DO 100 I=1,4
      IF(I.EQ.IINV) GO TO 100
      I1=I1+1
      J1=0
      DO 99 J=1,4
      IF(J.EQ.JINV) GO TO 99
      J1=J1+1
      RWORK(I1,J1)=R(I,J)
99    CONTINUE
100   continue
      RINV(JINV,IINV)=(-1)**(IINV+JINV)*DET(RWORK)
200   CONTINUE
      DET4=R(1,1)*RINV(1,1)+
     & R(1,2)*RINV(2,1)+R(1,3)*RINV(3,1)+R(1,4)*RINV(4,1)
      IF(DET4.EQ.0.)WRITE(5,103)DET4
      IF(DET4.EQ.0.) DET4=1.
103   FORMAT('   DISASTER  -  MATRIX DET = 0 - RESET TO 1')
      DO 201 I=1,4
      DO 201 J=1,4
201   RINV(J,I)=RINV(J,I)/DET4
C      WRITE(5,102)
C      DO 202 I=1,4
C102   FORMAT(/, 15X,'  CALLING MATRIX',
C     & 15X,'RECIPROCAL MATRIX ')
C      WRITE(5,101)(R(I,J),J=1,4),(RINV(I,K),K=1,4)
C202     CONTINUE
C101   FORMAT(4F10.4,5X,4F10.6)
      DO 3 I=1,NODER
      DO 3 J=1,NODER
3     RINV1(J,I)=RINV(J,I)
      RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE iprtsym(nun,n)
c       ==========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c This was broken - inverse translation component is not just negative of
c initial translation. Should be negative of inverse matrix * initial
c translation:
c
c [R']( [R] + T ) + T' = [I]
c [R'][R] + [R']T + T' = [I]
c [R'][R] = [I], thus [R']T + T' = 0
c T' = -[R']T
c
      integer j,k,n,nun
      real trans(3),inv_trans(3),zero_vec(3)
c
c Initialise translation component and zero vector for spinning
c
      do j=1,3
        trans(j)=vecs(j,n)
        zero_vec(j)=0.0e0
      enddo
c
c Use spin_inv to calculate [R']T
c NB As we are using a rotation centred on the origin, we could just as easily
c call spin instead of spin_inv
c
      call spin_inv(trans,inv_trans,ops_inv(1,1,n),zero_vec)
c
c Write out the answer. Hopefully it will now be right!
c
      write (nun,100)
     &  n,((ops_inv(k,j,n),j=1,3),k=1,3), (-inv_trans(j),j=1,3)
      if(record)write (nrec,100)
     &  n,((ops_inv(k,j,n),j=1,3),k=1,3), (-inv_trans(j),j=1,3)
100   format (15x,' Symm op n.o. ',I4,3(/,10x,3f10.5)/,10x,3f10.5)
c
c
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      FUNCTION irscale(rt,off,scal)
c       =============================
      IMPLICIT NONE
        integer irscale
        real rt,scal,off
c
        irscale=nint(rt*scal+off)
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE j1_to_j3(j1,j3,ngrid)
c       ================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        integer j1,j3(3),ngrid(3),jj1,jj2
c
        jj1=j1/ngrid(1)
        jj2=jj1/ngrid(2)
        j3(3)=jj2+1
        j3(2)=(jj1-jj2*ngrid(2)+1)
        j3(1)=(j1-jj1*ngrid(1))
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE j3_to_f(j3,x,map)
c       ============================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c
        integer j3(3),j,map
        real x(3)
        do j=1,3
        x(j)=float(j3(j)-norg(j,map))/float(nunit(j,map))
        enddo
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE label_env
c     ====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c----  new option 28th dec 93 j.g to set up a labelled from an input
c----  envelope defining 1 protomer/mol.  The envelope is overwritten
c----  and has incorporated the labelled envelopes (labelled with the ncs
c----  operator) for the ncs molecules.  This is then used for the simple
c----  averageing...using for the averaging onto mol1 the envelope labelled
c----  1...and then the other labels are used to interpolate the averaged
c----  density for mol1 back onto mols 2,3,4...etc
c--- local variables
c
      integer     map2,env2
      integer     iindex(3),jjx,jjy,jjz,iz,jz,iy,jy,ix,jx
      integer     nsym1,nsym2,ngadd,np,nc,nout,j_try,itype
      real        aindex(3),x(3),x1(3)
      character*8 G2CHRU,ch_temp,imap,ienv
      integer     nopt,keyword,j,intval
      integer*2   irho
      external    gp3,gp3_a
      integer*2   gp3,gp3_a
      logical     lass_get,save_op,use_in,i_check
c
c      equivalence  (i_temp,ch_temp)
c
      save_op = .false.
      i_check = .false.
      use_in = .true.
c
c
      if(.not.lass_get(ch_temp,'MAPI'))then
c
        ch_temp=G2CHRU('AV_LABEL_ENV:INPUT_MAP',1,8,0)
c
        call lass_set(ch_temp,'MAPI')
      endif
      imap=ch_temp
c
      if(.not.lass_get(ch_temp,'ENVI'))then
c
        ch_temp=G2CHRU('AV_LABEL_ENV:INPUT_ENV',1,8,0)
c
        call lass_set(ch_temp,'ENVI')
      endif
      ienv=ch_temp
      i_check=ienv.ne.'OFF'
c
c
      if(verbose)then
        write(6,110) imap,ienv
        if(record)write(nrec,110) imap,ienv
110     format(/' CREATE LABELLED ENVELOPE',/,
     &          ' ========================'
     &         /' Input map:       ',a,
     &         /' Input env:       ',a)
      else
        write(6,111) imap,ienv
        if(record)write(nrec,111) imap,ienv
111     format(' CREATE LABELLED ENVELOPE',
     &         '  I/P map:',a,'  I/P env:',a)
      endif
c
c       get slot numbers for all maps
c
      if(i_check)call av_ass(env2,ienv)
c
      call av_ass(map2,imap)
c----map 1 should exist by now !
      if(.not.defined(map2)) then
        write(6,160)
160     format(' %LABEL_ENV-ERR: Input map empty')
        return
      endif
c
c---- list symmetry operators to be used
c
      itype=1
      nsym1=1
      nsym2=nsym
c
10    nopt=keyword
     &  ('LABEL_ENV_SYMM','ALL !STAR!END !GO  ',0)
c
      if(nopt.eq.2)nsym1=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.eq.3)nsym2=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.ne.4)goto 10
c
      write(6,220)nsym1,nsym2
      if(record)write(nrec,220)nsym1,nsym2
220   format(' Will use symmetry operators from ',
     &             i4,' to',i4)
      if( (nsym1.gt.nsym2).or.(nsym2.gt.nsym) )then
        write(6,230)
        if(record)write(nrec,230)
230     format(' %LABEL_ENV-ERR: In choice of sym ops')
      endif
c
c---- warn if xtal symmetry switched on....very very dangerous
c
      if(xr_con(map2))then
        write(6,240)
        if(record)write(nrec,240)
240     format(
     &    ' %WARNING-ERR: **CRYSTALLOGRAPHIC SYMMETRY** switched on')
      endif
c
C---- incorporate bit from averaging....
c---- we will loop thro map....if pixel value = 1 goto continue
c---- therfore if pixel=0 loop thro inverse ops from 2 to nsym 
c---- and check to see if inv-ncs related pixel eq 1 (ie in env)
c---- then set pixel to be inv_ops no.
c
c---- firstly set up pointers for the maps that drive the process (backwards!)
c
      call init_map(map2)
      call init_map(env2)
c
c---- initialize r-factors, corr coefs, etc
c
      ngadd=0
      np=0
      nc=0
      nout=0
      save_op=.false.
      j_try=0
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
      do iz=1,n_brick(3,map2)
c
        jz=(iz-1)*brick(3,map2)
c
        do iy=1,n_brick(2,map2)
c
          jy=(iy-1)*brick(2,map2)
c
          do ix=1,n_brick(1,map2)
c
            jx=(ix-1)*brick(1,map2)
c
c---- now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map2)
c
              if(jjz.le.nx(3,map2))then
c
                do jjy=jy+1,jy+brick(2,map2)
c
                  if(jjy.le.nx(2,map2))then
c
                    do jjx=jx+1,jx+brick(1,map2)
c               
                      if(jjx.le.nx(1,map2))then
c
c--- check to see whether we are in input env
c
      if((gp3(env2,jjx,jjy,jjz).eq.0).and.i_check)goto 200
c
c--- check to see whether we already in env for mol1
c--- if so ignore the pixel
      if(gp3(map2,jjx,jjy,jjz).ne.0) then
        nout=nout+1
        goto 200
      endif
c
c---- really no of pixels in env...its confusing i know
c
      iindex(1)=jjx
      iindex(2)=jjy
      iindex(3)=jjz
c
      call j3_to_f(iindex,x1,map2)
      call f_to_o(x1,x,map2)
c
c
c       
c---- first of all try the operator that worked last time
c
      if(save_op)then
        j=j_try
c
c---- get the appropriate postion for this copy:
c
        call spin_inv(x,x1,ops_inv(1,1,j),vecs(1,j))
c
        call o_to_j3(x1,aindex,iindex,map2)
c
c---- now get the interpolated electron density:
c
c--- if not using xtal symm....much quicker
        if(xr_con(map2))then
c
          call g_p_i_i_x(map2,irho,aindex,iindex,itype,j,use_in
     &                                    ,np,nc)
c
c   (  itype specifies interp: 1, 8 or 11 point. )
c
        else
          call g_p_i_i(map2,irho,aindex,iindex,itype,use_in,np,nc)
        endif
c
        if(irho.eq.1) then
          irho=j
c---- and now finally store labelled envelope pixel.
c
          call pp3(map2,jjx,jjy,jjz,irho)
          ngadd=ngadd+1
          goto 200
        else
          save_op=.false.
        endif
c
      endif
c
c
c---- loop over the NCS operators from 2 to final one:
c
      do j=nsym1,nsym2
c
c---- get the appropriate postion for this copy:
c
        call spin_inv(x,x1,ops_inv(1,1,j),vecs(1,j))
c
        call o_to_j3(x1,aindex,iindex,map2)
c
c---- now get the interpolated electron density:
c
c--- if not using xtal symm....much quicker
        if(xr_con(map2))then
c
          call g_p_i_i_x(map2,irho,aindex,iindex,itype,j,use_in
     &                                    ,np,nc)
c
c       ( itype specifies interp: 1, 8 or 11 point. )
c
        else
          call g_p_i_i(map2,irho,aindex,iindex,itype,use_in,np,nc)
        endif
c
        if(irho.eq.1) then
          irho=j
c---- and now finally store labelled envelope pixel.
c
          call pp3(map2,jjx,jjy,jjz,irho)
          j_try=j
          save_op=.true.
          ngadd=ngadd+1
          goto 200
        endif
c
      enddo
c
c
c
200   continue
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c----phew, averaged, lets see how we did ...
c
c999   continue
      if(verbose)then
        write(6,998)ngadd,nout
        if(record)write(nrec,998)ngadd,nout
998     format(' LABELLED ENVELOPE EXPANSION FINISHED',/,
     &         ' ====================================',/,
     &      i9,' pixels were reconstructed from ', i9 ,' pixels')
      else
        write(6,997)ngadd,nout
        if(record)write(nrec,997)ngadd,nout
997     format(' Env produced:',
     &      i9,' pixels expanded from', i9 ,' pixels')
      endif
c
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE lass_add
c       ===================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        character*8 G2CHRU,label,stream,label_d,stream_d
        integer j
c
        label  =G2CHRU('AV_ASSIGN:MAP_LABEL',1,8,0)
        stream =G2CHRU('AV_ASSIGN:LOGICAL_STREAM',1,8,0)
        goto 999
c
        entry lass_set(label_d,stream_d)
c       ================================
c
        label=label_d
        stream=stream_d
c
999     continue
c
        do j=1,n_stream
          if(stream.eq.streams(j)) then
            labelsa(j)=label
            if(a_made(j))then
              write(6,10)label,stream
              if(record)write(nrec,10)label,stream
10            format(' ',a,' assigned to stream '
     &              ,a,'(replacement)')
            else
              write(6,15)label,stream
              if(record)write(nrec,15)label,stream
15            format(' ',a,' assigned to stream '
     &              ,a,'(new)')
            endif
            a_made(j)=.true.
            return
          endif
        enddo
c
        do j=1,n_stream
          if(streams(j).eq.'    ')then
            streams(j)=stream
            a_made(j)=.true.
            labelsa(j)=label
            write(6,20)label,stream
            if(record)write(nrec,10)label,stream
20          format(' ',a,' assigned to stream '
     &              ,a,'(new)')
            return
          endif
        enddo
c
        write(6,30)
        if(record)write(nrec,30)
30      format(' %LASS_ADD-ERR: Assign fails, no free space')
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      FUNCTION lass_get(label,stream)
c       ===============================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        character*(*) label,stream
        logical lass_get
        integer j
c
        lass_get=.false.
c
        do j=1,n_stream
          if(streams(j).eq.stream) then
            if(a_made(j))then
              label=labelsa(j)
              lass_get=.true.
              return
            else
              lass_get=.false.
              return
            endif
          endif
        enddo
c
        write(6,10)stream
        if(record)write(nrec,10)stream
10      format(' %LASS_GET-ERR: Can''t find internal stream ',a)
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE lass_init
c       ====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        integer j
c
        do j=1,n_stream
          streams(j)='    '
          labelsa(j)='OFF '
          a_made(j)=.false.
        enddo
        streams(1)='MAPI'
        streams(2)='MAPO'
        streams(3)='ENVI'
        streams(4)='ENVO'
        streams(5)='ENVM'
        a_made(3)=.true.
        a_made(4)=.true.
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE lass_q
c       =================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        integer j
c
        write(6,10)
        if(record)write(nrec,10)
10      format(' LOGICAL STREAM ASSIGMENTS:')
c
        do j=1,n_stream
          if(a_made(j))then
            write(6,20)labelsa(j),streams(j)
            if(record)write(nrec,20)labelsa(j),streams(j)
20          format(' Map ',a,' assigned to stream ',a)
          endif
        enddo
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE list_map(flag)
c       =========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c
        external    gp3_a
        integer*2   gp3_a
        character*8 map
        integer     iform, x(3,3), uvw(3), nmap
        integer     line(500), local_ind(3)
        real        scale_r
        character*1 axes(3)
        logical     m_done,flag
c
        character*8 G2CHRU,ch_temp
        real        realval
        integer     j,nopt,keyword,intval,k,ll,l,lll
        logical     lass_get
c
c
        axes(1)='x'
        axes(2)='y'
        axes(3)='z'
c
        m_done=.false.
        scale_r=0.05
        iform=4
        do j=1,3
          uvw(j)=j
        enddo
        do j=1,3
          x(3,j) = 8
        enddo
c
c---- we might have a logical assignment setup for this stream (MAPI)
c
        if(lass_get(ch_temp,'MAPI'))then
          m_done=.true.
          map=ch_temp
          call av_ass(nmap,map)
          if(.not.defined(nmap)) then
            write(6,*)'%LIST_MAP-ERR: Cannot print undefined map'
            if(record)write(nrec,*)
     &                '%LIST_MAP-ERR: Cannot print undefined map'
            else
            do j=1,3
              x(1,j) = 1
              x(2,j) = nx(j,nmap)
            enddo
            endif
          endif
c
        if(flag)then
          x(3,1)=2
          x(3,2)=2
          x(3,3)=4
        goto 99
        endif
c
10      nopt=keyword('AV_L_MAP','MAP !SCAL!STAR!END !INTE!GO  !BACK',0)
c
c
        if(nopt.eq.1)then
          ch_temp=G2CHRU('AV_L_MAP',1,8,0)
          map=ch_temp
          m_done=.true.
        call av_ass(nmap,map)
        if(.not.defined(nmap)) then
          write(6,*)'%LIST_MAP-ERR: Cannot print undefined map'
          if(record)write(nrec,*)
     &              '%LIST_MAP-ERR: Cannot print undefined map'
          goto 10
          endif
          do j=1,3
            x(1,j) = 1
            x(2,j) = nx(j,nmap)
          enddo
        endif
        if(nopt.eq.2)scale_r=realval('AV_L_MAP:ENTER_CUT',
     &                    -10000.0,10000.0,0)
c
        if(nopt.eq.3)then
        if(.not.m_done)then
          write(6,*)'%LIST_MAP-ERR: define map first'
          if(record)write(nrec,*)
     &              '%LIST_MAP-ERR: define map first'
        endif
          do j=1,3
            x(1,j)=intval('AV_L_MAP:ENTER_START',1,nx(j,nmap),0)
          enddo
        endif
c
        if(nopt.eq.4)then
        if(.not.m_done)then
          write(6,*)'%LIST_MAP-ERR: define map first'
          if(record)write(nrec,*)
     &              '%LIST_MAP-ERR: define map first'
        endif
          do j=1,3
            x(2,j)=intval('AV_L_MAP:ENTER_END',x(1,j),nx(j,nmap),0)
          enddo
        endif
c
        if(nopt.eq.5)then
          do j=1,3
            x(3,j)=intval('AV_L_MAP:ENTER_INT',1,100,0)
          enddo
        endif
c
        if(nopt.eq.7)return
c
        if(nopt.ne.6) goto 10
c
99      if(.not.m_done)then
          write(6,*)'%LIST_MAP-ERR: No map defined for printing'
          if(record)write(nrec,*)
     &              '%LIST_MAP-ERR: No map defined for printing'
        goto 10
        endif
c
        do j=1,3
          x(1,j) = 1
          x(2,j) = nx(j,nmap)
        enddo
c
          if(.not.terse)then
            write(6,20)scale_r,iform,x,uvw
            if(record)write(nrec,20)scale_r,iform,x,uvw
20          format(' PRINT MAP, scaled by:',e10.2,'  format  I',I1,
     &      /,' start end step(for x y and z)', 3(/(3I5)),/,' uvw:',3I4)
          else
            write(6,22)scale_r,iform
            if(record)write(nrec,22)scale_r,iform
22          format(' PRINT MAP, scaled by:',e10.2,'  format  I',I1)
          endif
          do j=x(1, uvw(3)), x(2, uvw(3)), x(3, uvw(3))
            write(6,30) axes(uvw(3)),(nstart(uvw(3),nmap)+j-1),
     &                  axes(uvw(1))
            if(record)write(nrec,30)
     &               axes(uvw(3)),(nstart(uvw(3),nmap)+j-1),axes(uvw(1))
30          format(/,'  section  ', a,'=', I4,/,'   ----------->  ',a)
            do k=x(1, uvw(2)), x(2, uvw(2)), x(3, uvw(2))
              ll=0
              do l=x(1, uvw(1)), x(2, uvw(1)), x(3, uvw(1))
                ll=ll+1
                  local_ind(uvw(1))=l
                  local_ind(uvw(2))=k
                  local_ind(uvw(3))=j
                line(ll)=
     &            gp3_a(nmap,local_ind)*scale_r
              enddo
              write(6,40)(line(lll),lll=1,ll)
              if(record)write(nrec,40)(line(lll),lll=1,ll)
40            format(' ',20I4)
            enddo
          enddo
c       endif
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE list_mask(flag)
c       ==========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c
        external    gp3_a
        integer*2   gp3_a,local
        character*8 map
        integer     iform,x(3,3),uvw(3),nmap,i_temp
        integer     local_ind(3),j,nopt,keyword,intval,k,ll,l,lll
        real        cut
        character*1 axes(3)
        logical     m_done,flag
        character*1 line(500)
        character*1 lookup(60)
        character*8 G2CHRU,ch_temp
        real        realval
        logical     lass_get
c
        equivalence  (i_temp,ch_temp)
        data lookup/'*','2','3','4','5','6','7','8','9','0'
     &             ,'a','b','c','d','e','f','g','h','i','j'
     &             ,'k','l','m','n','o','p','q','r','s','t'
     &             ,'u','v','w','x','y','z','A','B','C','D'
     &             ,'E','F','G','H','I','J','K','L','M','N'
     &             ,'O','P','Q','R','S','T','U','V','W','X'/
c
c
        axes(1)='x'
        axes(2)='y'
        axes(3)='z'
c
        m_done=.false.
        cut=0.0
        iform=4
        do j=1,3
          uvw(j)=j
        enddo
c
        do j=1,3
          x(3,j) = 2
        enddo
        x(3,uvw(3))=8
c
c---- we might have a logical assignment setup for this stream (MAPI)
c
        if(lass_get(ch_temp,'MAPI'))then
          m_done=.true.
          map=ch_temp
          call av_ass(nmap,map)
          if(.not.defined(nmap)) then
            write(6,*)'%LIST_MAP-ERR: Cannot print undefined map'
            if(record)write(nrec,*)
     &                '%LIST_MAP-ERR: Cannot print undefined map'
            else
            do j=1,3
              x(1,j) = 1
              x(2,j) = nx(j,nmap)
            enddo
            endif
          endif
c
        if(flag)then
          x(3,1)=1
          x(3,2)=1
          x(3,3)=2
        goto 99
        endif
c
10      nopt=keyword('AV_L_MASK','MAP !SCAL!STAR!END !INTE!GO  !BACK',0)
c
        if(nopt.eq.1)then
          ch_temp=G2CHRU('AV_L_MASK',1,8,0)
          map=ch_temp
          m_done=.true.
        call av_ass(nmap,map)
        if(.not.defined(nmap)) then
          write(6,*)'%LIST_MASK-ERR: Cannot print undefined map'
          if(record)write(nrec,*)
     &              '%LIST_MASK-ERR: Cannot print undefined map'
          goto 10
          endif
          do j=1,3
            x(1,j) = 1
            x(2,j) = nx(j,nmap)
          enddo
        endif
c
        if(nopt.eq.2)cut=realval('AV_L_MASK:ENTER_CUT',
     &                    -10000.0,10000.0,0)
c
        if(nopt.eq.3)then
        if(.not.m_done)then
          write(6,*)'%LIST_MASK-ERR: define map first'
          if(record)write(nrec,*)
     &              '%LIST_MASK-ERR: define map first'
        endif
          do j=1,3
            x(1,j)=intval('AV_L_MASK:ENTER_START',1,nx(j,nmap),0)
          enddo
        endif
c
        if(nopt.eq.4)then
        if(.not.m_done)then
          write(6,*)'%LIST_MASK-ERR: define map first'
          if(record)write(nrec,*)
     &              '%LIST_MASK-ERR: define map first'
        endif
          do j=1,3
            x(2,j)=intval('AV_L_MASK:ENTER_END',x(1,j),nx(j,nmap),0)
          enddo
        endif
c
        if(nopt.eq.5)then
          do j=1,3
            x(3,j)=intval('AV_L_MASK:ENTER_INT',1,100,0)
          enddo
        endif
c
        if(nopt.eq.7)return
c
        if(nopt.ne.6) goto 10
c
99      if(.not.m_done)then
          write(6,*)'%LIST_MASK-ERR: No map defined for printing'
          if(record)write(nrec,*)
     &              '%LIST_MASK-ERR: No map defined for printing'
          goto 10
        endif
c
c
        do j=1,3
          x(1,j) = 1
          x(2,j) = nx(j,nmap)
        enddo
c
          if(.not.terse)then
            write(6,20)cut,iform,x,uvw
            if(record)write(nrec,20)cut,iform,x,uvw
20          format(' PRINT MAP, dots <=:',e10.2,'  format  I',I1,
     &      /,' start end step(for x y and z)', 3(/(3I5)),/,' uvw:',3I4)
          else
            write(6,22)cut,iform
            if(record)write(nrec,22)cut,iform
22          format(' PRINT MAP, dots <=:',e10.2,'  format  I',I1)
          endif
c
          do j=x(1, uvw(3)), x(2, uvw(3)), x(3, uvw(3))
            write(6,30) axes(uvw(3)),(nstart(uvw(3),nmap)+j-1),
     &                  axes(uvw(1))
            if(record)write(nrec,30)
     &             axes(uvw(3)),(nstart(uvw(3),nmap)+j-1),axes(uvw(1))
30          format(/,'  section  ', a,'=', I4,/,'   ----------->  ',a)
            do k=x(1, uvw(2)), x(2, uvw(2)), x(3, uvw(2))
              ll=0
              do l=x(1, uvw(1)), x(2, uvw(1)), x(3, uvw(1))
                ll=ll+1
                  local_ind(uvw(1))=l
                  local_ind(uvw(2))=k
                  local_ind(uvw(3))=j
                line(ll)='.'
                local=gp3_a(nmap,local_ind)
                if(local.gt.cut)then
                  if(local.gt.60)local=local-60
                  if(local.gt.60)local=local-60
                  line(ll)='#'
                  if(local.le.60)line(ll)=lookup(local)
                endif
              enddo
              write(6,40)(line(lll),lll=1,ll)
              if(record)write(nrec,40)(line(lll),lll=1,ll)
40            format(' ',132a1)
            enddo
          enddo
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE list_sym
c       ===================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        integer intval,nsym1,nsym2,j
c
c
        if (nsym.lt.1) then
          write(6,10)
            if(record) write(nrec,10)
          return
        endif
10      format(' No Symmetry Operators defined ')
c
        nsym1=intval('AV_LIST_SYM:ENTER_1_NUMBER',1,nsym,0)
        nsym2=intval('AV_LIST_SYM:ENTER_2_NUMBER',nsym1,nsym,0)
        do j=nsym1,nsym2
        call prtsym(6,j)
        enddo
        return
        end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE loopover(map1,env1,i_check,apply,av_dens2)
c     =====================================================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c     =====================
c
c----local variables:
c
      integer   map1,env1,ngadd,iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
      logical   i_check,apply
      integer*2 irho
      external  gp3,gp3_a
      integer*2 gp3,gp3_a
      real      r_max,r_min,av_dens2
c
c---- the averaging process is very simple, we are driven by the
c---- output map, for each pixel in that map (accessed in a sequential
c---- fashion, but since map is bricked we move slowly through space)
c---- we pick up all the symm related pixels in the input map and take
c---- the average.
c
c---- firstly set up pointers for the maps that drive the process (forwards!)
c
      call init_map(map1)
c
c
      ngadd=0
c
      call r_stat
c
      if(i_check)call init_map(env1)
c
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
      do iz=1,n_brick(3,map1)
c
        jz=(iz-1)*brick(3,map1)
c
        do iy=1,n_brick(2,map1)
c
          jy=(iy-1)*brick(2,map1)
c
          do ix=1,n_brick(1,map1)
c
            jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                  if(jjy.le.nx(2,map1))then
c
                    do jjx=jx+1,jx+brick(1,map1)
c               
                      if(jjx.le.nx(1,map1))then
c
c---- check that this pixel is in required region of the output map
c
c
      if(.not.apply)then
c--- at the mo we will accumulate statistics from the map that is within
c--- the envelope.  therefore checking to see whether
c--- we are within an envelope.  BUT we then apply scale the map over
c--- the complete map....seems to make more physical sense.
c
        if(i_check)then
          if(gp3(env1,jjx,jjy,jjz).eq.0)goto 200
        endif

        irho=gp3(map1,jjx,jjy,jjz)
        call r_a_stat(irho)
      endif
c
      if(apply) then
        irho=gp3(map1,jjx,jjy,jjz)
        irho=( (irho-av_dens2)*sca_map) + shi_map +
     &                   av_dens2
        call pp3(map1,jjx,jjy,jjz,irho)
        call r_a_stat(irho)
      endif
c
      ngadd=ngadd+1

c
200   continue
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c----phew, averaged, lets see how we did ...
c
      if(.not.apply)then
c990     continue
        if(verbose)then
          write(6,991)ngadd
          if(record)write(nrec,991)ngadd
991       format(' SCALING MAP STATISTICS ',/,
     &           ' ======================',/,
     &       i10,' pixels were used from i/p map')
        else
          write(6,992)ngadd
          if(record)write(nrec,992)ngadd
992       format(' SCALING MAP STATISTICS,',
     &       i10,' pixels were used from i/p map')
        endif
      endif

      if(apply)then
c999     continue
        if(verbose)then
          write(6,998)ngadd
          if(record)write(nrec,998)ngadd
998       format(' MAP SCALING APPLIED ',/,
     &           ' ===================',/,
     &       i10,' pixels were used from i/p map')
        else
          write(6,996)ngadd
          if(record)write(nrec,996)ngadd
996       format(' MAP SCALING APPLIED,',
     &       i10,' pixels were used from i/p map')
        endif
      endif
c
c
      call r_p_stat(r_max,r_min)
c
      if(apply)sollev(map1)=nint(((sollev(map1)-av_dens2)
     &       *sca_map)+shi_map+av_dens2)
c
c---- print a summary of the statistics
c
      write(6,10)r_max,r_min,av_temp,sigm,sollev(map1)
c
      if(record) then
        write(nrec,10)r_max,r_min,av_temp,sigm,sollev(map1)
      endif
c
10    format(' New statistics on map pixels will replace old RHOLIM'/,
     &       ' Max rho',f8.1,
     &       '    min rho',f8.1,
     &       '    mean rho',f8.1,
     &       '    sigma',f8.1,
     &       '    solvent',i8)
c
c
c----av_temp is in common block----
c
      if(apply)then
        rholim(1,map1)=r_min
        rholim(2,map1)=r_max
        rholim(3,map1)=av_temp
        rholim(4,map1)=sigm
      endif
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE make_env
c     ===================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c---- local arrays, note we allow upto max_o spheres and max_o planes
c---- reset max_o for more complex systems
c
C RME: 23/4/2010: Made new J8 variable to be INTEGER*8 pointer for scratch
C RME: 23/4/2010: Shoudl be used for all access to scratch()
      integer*8   j8
C RME: 23/4/2010: End of change
      integer     max_o,nopt,keyword,j,intval,jplan,jcirc,inner
      parameter   (max_o=100)
c
      integer*2   val(max_o),iinout(max_o),plan_v(max_o),gp3_a
      integer     map1,map2,iindex(3),nsym1,nsym2,jjj,nopt2
      character*8 map,env
      real        cen(3,max_o),rad(max_o),plan(3,2,max_o),del_p(2)
      real        x1(3),x(3)
c
      real local_prot(3)
      integer local_v
c
      character*8 G2CHRU,ch_temp
      real        realval,xrvol,local_cen,local_rad,local_val
      real        ar_sq,del,a_rad,pad_d
      integer     a_pad(3)
      integer     n_reset,ii,jj,kk,iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
      integer     ics_local,logic,ncirc,nplan,i,nplan_old
      logical     lass_get,usei,plan_or,prot_f,ival_f
c
c--- variables assoc with atom
      integer      iflag, iii(3), ind_n(3)
      real         f_ind_n(3), local_at(3,natom)
      integer*2    vala
      real         sqdist, dist
      character*80 fil
c
      local_v=0
c
      if(.not.lass_get(ch_temp,'MAPI'))then
c
        ch_temp=G2CHRU('AV_MAKE_ENV:INPUT_MAP',1,8,0)
c
        call lass_set(ch_temp,'MAPI')
      endif
      map=ch_temp
c
c
      if(.not.lass_get(ch_temp,'ENVM'))then
c
        ch_temp=G2CHRU('AV_MAKE_ENV:OUTPUT_ENV',1,8,0)
c
        call lass_set(ch_temp,'ENVM')
      endif
      env=ch_temp
c
c
      write(6,20)env,map
      if(record)write(nrec,20)env,map
20    format(' MAKE ENVELOPE, LOGICAL NAME:',a,' map:',a,
     &  ' will act as a template.')
c
c
      call av_ass(map1, map)
c
      call av_ass(map2, env)
c
c
c       use this map as a template for envelope
c       
      call copy_head(map1,map2)
c
c---- reset scale and offset for an envelope (in case we write it out)
c
      if(.not.terse)then
        write(6,*)
     &    'NB scale and offset reset to 1 and 0 for this envelope'
        if(record)then
          write(nrec,*)
     &      'NB scale and offset reset to 1 and 0 for this envelope'
        endif
      endif
c
      mscale(map2)=1.0
      moffset(map2)=0.0
      rholim(1,map2)=0.0
      rholim(2,map2)=1.0
      rholim(3,map2)=0.0
      rholim(4,map2)=0.0
      usei=.false.
c
c---- check if we want to redefine the output grid (eg coarser than input)
c
30    nopt=keyword('AV_AVERAGE_UPDATE','UPDA!RESE!GO  !BACK!USEI',0)
c
      if(nopt.eq.1) then
        do j=1,3
          nstart(j,map2)=intval('AV_MAKE_ENV:NSTART',-3000,3000,0)
          nend(j,map2)  =intval('AV_MAKE_ENV:NEND',-3000,3000,0)
          nunit(j,map2) =intval('AV_MAKE_ENV:NUNIT',-3000,3000,0)
          nx(j,map2)=nend(j,map2)-nstart(j,map2)+1
          norg(j,map2)=1-nstart(j,map2)
          if(.not.terse)then
            write(6,200)nstart(j,map2),nend(j,map2),
     &        nx(j,map2),nunit(j,map2)
            if(record)
     &        write(nrec,200)nstart(j,map2),nend(j,map2),
     &          nx(j,map2),nunit(j,map2)
200         format(' ',4i5)
          endif
        enddo
C
        npix(map2)=nx(1,map2)*nx(2,map2)*nx(3,map2)
C
        write(6,210)npix(map2)
        if(record)write(nrec,210)npix(map2)
210     format(' There will be',i9,' pixels in the output map')
      endif
c
      if(nopt.eq.2)then
        do j=1,3
          iuvw(j,map2)=intval('AV_AVERAGE:ENTER_AXIS',1,3,0)
        enddo
      endif
c
      if(nopt.eq.4)return
c
      if(nopt.eq.5)then
        usei=.true.
        write(6,215)
        if(record)write(nrec,215)
215     format(' The input map will be used as a starting envelope,',
     &    ' which will be further elaborated')
      endif
c
      if(nopt.ne.3)goto 30
c
c
      call XRfrac(XRcell(1,map2),XRtr(1,1,map2),XRintr(1,1,map2),
     &  XRvol,.true.)
c
c---- set up symmetry data
c
      call XR_symsy(XRcell(1,map2),lgrp(map2),ics_local,
     &  n_xsym(map2),XR_sym(1,map2) )
      if(.not.terse)then
        write(6,220)n_xsym(map2)
        if(record)write(nrec,220)n_xsym(map2)
220     format(' Number of symmetry ops for this space group=',i4)
      endif
c      write(6,*)' xtal sym op one',(XR_sym(j,map2),j=1,12)
c
c---- and integerized CSO's:
c
      call XR_isym(map2)
c
c
c---- also set up pointers indicating which symm op to try first
c
      do j=1,maxsym
        last_sym(j,map2)=1
      enddo
c
c
c---- set up bricking data
c
      call brickit(map2)

c
      call alloc(map2)
c
c---- all space allocated
c
c
c----set starting value
c
      if(.not.usei)then
        logic=intval('AV_MAKE_ENV:ENTER_START_VAL',-32000,32000,0)
        write(6,230)logic
        if(record)write(nrec,230)logic
230     format(' Set all of envelope to value:',i5)
        do j8=ns(map2),ne(map2)
          scratch(j8)=logic
        enddo
        if(.not.terse)then
          write(6,232)
          if(record)write(nrec,232)
232       format(' Envelope initialized, now define objects')
        endif
      endif
c
c
      ncirc=0
      nplan=0
      nsym1=1
      nsym2=nsym
      plan_or=.false.
c
c
40    nopt=keyword('AV_MAKE_ENV:OBJECTS','SPHE!PLAN!AUTO!GO  !BACK'//
     &  '!ATOM!AUTP!P_OR!P_AN',0)
c
c
      if(nopt.eq.1) then
        if(ncirc.gt.max_o)then
          write(6,260)max_o
          if(record)write(nrec,260)max_o
260       format(' **WARNING** The number of spheres cannot
     &                    be set any higher, limit is',i6)
          goto 40
        endif
        ncirc=ncirc+1
        do j=1,3
          cen(j,ncirc)
     &      =realval('AV_MAKE_ENV:ENTER_CEN',-100000.0,100000.0,0)
        enddo
        rad(ncirc)=realval('AV_MAKE_ENV:ENTER_RAD',0.000001,100000.0,0)
        val(ncirc)=intval('AV_MAKE_ENV:ENTER_NEW_VAL',-100000,100000,0)
        iinout(ncirc)=
     &    intval('AV_MAKE_ENV:ENTER_INOUT',0,1,0)
c
        write(6,240)
     &    (cen(j,ncirc),j=1,3),rad(ncirc),val(ncirc),iinout(ncirc)
        if(record)write(nrec,240)
     &    (cen(j,ncirc),j=1,3),rad(ncirc),val(ncirc),iinout(ncirc)
240     format(' Centre:',3f8.1,'   radius:',f8.1,
     &    '   value:',i3,'   inout:',i3)
        rad(ncirc)=rad(ncirc)**2
      endif
c
c
      if(nopt.eq.2) then
        if(nplan.gt.max_o)then
          write(6,261)max_o
          if(record)write(nrec,261)max_o
261       format( ' **Warning** The number of planes cannot',
     &                  ' be set any higher, limit is',i6)
          goto 40
        endif
        nplan=nplan+1
        do i=1,2
          do j=1,3
            plan(j,i,nplan)=realval('AV_MAKE_ENV_COORDS:',
     &                                -100000.,100000.,0)
          enddo
        enddo
        plan_v(nplan)=intval('AV_MAKE_ENV_VALUE:',-32000,32000,0)
        write(6,263)
     &    ((plan(j,i,nplan),j=1,3),i=1,2),plan_v(nplan)
        if(record)write(nrec,263)
     &    ((plan(j,i,nplan),j=1,3),i=1,2),plan_v(nplan)
263     format(' Planes:(endpoints, value)',2(' ',3f8.1),i4)
c
      endif
c
c
c JMD!PORT: THis looks improperly implemented
      if(nopt.eq.3)then
c
c---- generate contact planes automatically using crystallographic symmetry
c---- to find the neighbours
c
        do j=1,3
          local_cen=realval('AV_MAKE_ENV:CENTRE',-100000.0,100000.0,0)
        enddo
        local_rad=realval('AV_MAKE_ENV:RADIUS',0.0,100000.0,0)
        local_val=intval('AV_MAKE_ENV:VALUE',-32000,32000,0)
c
        nplan_old=nplan
c JMD!PORT: Why is this call commented out?
c        call set_auto(local_val,local_cen,local_rad,nplan,map2)
        nplan_old=nplan-nplan_old
        write(6,270)nplan_old
        if(record)write(nrec,270)nplan_old
270     format(' The number of planes set up with AUTO option was:',i6)
      endif
c
c
      if(nopt.eq.7) then
c
        pad_d = 1.0
        prot_f=.false.
        ival_f=.false.
c
280     nopt2=keyword
     &    ('AV_AVERAGE_SYMM','PROT!VAL !ALL !STAR!END !PAD!GO  ',0)
c
        if(nopt2.eq.1)then
          do j=1,3
            local_prot(j)=realval
     &        ('AV_MAKE_ENV:PROT_XYZ',-100000.0,100000.0,0)
          enddo
          prot_f=.true.
        endif
c
        if(nopt2.eq.2)then
          local_v=intval('AV_MAKE_ENV_VALUE:',-32000,32000,0)
          ival_f=.true.
        endif
c
c
        if(nopt2.eq.4)nsym1=intval('AV_AVERAGE_SYMM',1,nsym,0)
        if(nopt2.eq.5)nsym2=intval('AV_AVERAGE_SYMM',1,nsym,0)
        if(nopt2.eq.6)pad_d=realval
     &    ('AV_MAKE_ENV:PAD_VAL',0.0,100000.0,0)
        if(nopt2.ne.7)goto 280
c
c
        if(.not.prot_f) then
          write(96,278)
          if(record)write(nrec,278)
          return
        endif
278     format(' %AVERAGE-ERR: PROT_XYZ not defined')
c
        if(.not.ival_f) then
          write(96,279)
          if(record)write(nrec,279)
          return
        endif
279     format(' %AVERAGE-ERR: ENV_IVAL not defined')

        write(6,281)nsym1,nsym2
        if(record)write(nrec,281)nsym1,nsym2
281     format(' Using ncs ops',
     &             i4,' to',i4,' to generate protomer env')
c
        if( (nsym1.gt.nsym2).or.(nsym2.gt.nsym) )then
          write(6,282)
          if(record)write(nrec,282)
282       format(' %AVERAGE-ERR: In choice of symmetry operators ')
          return
        endif
c
c
        do jjj = nsym1,nsym2
c
          if(nplan.gt.max_o)then
            write(6,283)max_o
            if(record)write(nrec,283)max_o
283         format( ' **Warning** The number of planes cannot',
     &                  '  be set any higher, limit is',i6)
            goto 40
          endif
c
          nplan=nplan+1
c
          call spin(local_prot,x1,ops(1,1,jjj),vecs(1,jjj))
c
          plan_v(nplan)=local_v
c
          dist = 0.0
          do j=1,3
            dist = dist + ( x1(j) - local_prot(j) )**2
          enddo
          dist = sqrt(dist)
c
          do j=1,3
            x1(j) = x1(j) + pad_d*(x1(j)-local_prot(j))/dist
          enddo
c
          do ii=1,3
            plan(ii,1,nplan) = x1(ii)
            plan(ii,2,nplan) = local_prot(ii)
          enddo
c
          if(verbose) then
            write(6,284)
     &        ((plan(j,i,nplan),j=1,3),i=1,2),plan_v(nplan)
            if(record)write(nrec,284)
     &        ((plan(j,i,nplan),j=1,3),i=1,2),plan_v(nplan)
          endif
284       format(' Planes:(endpoints, value)',2(' ',3f8.1),i5)
c
        enddo
c
      endif
c
      if(nopt.eq.8)then
        plan_or=.true.
        write(6,482)
        if(record)write(nrec,482)
482     format(' Will set pixels outside any plane to specified value')
      endif
c
      if(nopt.eq.9)then
        plan_or=.false.
        write(6,483)
        if(record)write(nrec,483)
483     format(' Will set pixels inside all planes to specified value')
      endif
c
      if(nopt.eq.5) return
c
      if(nopt.eq.6) then
c
        call read_pdb(fil)
c
        local_rad=3.5
        vala=1
        iflag=1
400     nopt=keyword('AV_MAKE_ENV:','CALP!MAIN!ALL !GO  ',0)
c
        if(nopt.eq.1) then
          a_rad=realval('AV_MAKE_ENV:RADIUS',0.0,100000.0,0)
          vala=intval('AV_MAKE_ENV:ENTER_NEW_VAL',-32000,32000,0)
          iflag=1
          write(6,500) a_rad
          if(record)write(nrec,500) a_rad
500       format(' Envelope calculation using C-alphas, with radius:',
     &      f8.1)
        endif
c
        if(nopt.eq.2) then
          a_rad=realval('AV_MAKE_ENV:RADIUS',0.0,100000.0,0)
          vala=intval('AV_MAKE_ENV:ENTER_NEW_VAL',-32000,32000,0)
          iflag=2
          write(6,510) a_rad
          if(record)write(nrec,510) a_rad
510       format(' Envelope calculation using main chain, with radius:',
     &      f8.1)
        endif
c
        if(nopt.eq.3) then
          a_rad=realval('AV_MAKE_ENV:RADIUS',0.0,100000.0,0)
          vala=intval('AV_MAKE_ENV:ENTER_NEW_VAL',-32000,32000,0)
          iflag=3
          write(6,520) a_rad
          if(record)write(nrec,520) a_rad
520       format(' Envelope calculation using all atoms, with radius:',
     &      f8.1)
        endif
c
        if(nopt.ne.4)goto 400
c
        ar_sq=a_rad**2
c
c- Bug-fix - boxing to save time broke down for non-orthonganol cell axes
c- Bodged 24/6/2002 JMD so that it looks in pixel space instead of
c- orthoganol space, and doubles the box dimensions (=8 times the pixels
c- to search over!) for non-orthogonal axes. Doing the job properly requires
c- finding the vectors perpendicular to the ab, ac and bc planes, moving
c- +- a_rad along those from the centre and mapping those points into pixel
c- space. NB each vector at a time rather than just the corners of the box!
c
        a_pad(1)=nint(float(nunit(1,map2))*a_rad/xrcell(1,map2))
        a_pad(2)=nint(float(nunit(2,map2))*a_rad/xrcell(2,map2))
        a_pad(3)=nint(float(nunit(3,map2))*a_rad/xrcell(3,map2))

        if((xrcell(4,map2).ne.90.0).and.
     &     (xrcell(5,map2).ne.90.0).and.
     &     (xrcell(6,map2).ne.90.0))then
          a_pad(1)=a_pad(1)*2
          a_pad(2)=a_pad(2)*2
          a_pad(3)=a_pad(3)*2
        endif
c
        n_reset=0
c
        do i=1,p_na
c
          if(iflag.eq.1) then
c***remember to declare local_at
            if(p_anam(i).eq.' CA  ')then
              local_at(1,i)=p_xyz(1,i)
              local_at(2,i)=p_xyz(2,i)
              local_at(3,i)=p_xyz(3,i)
            endif
          else if(iflag.eq.2)then
            if( (p_anam(i).eq.' CA  ') .or. (p_anam(i).eq.' C   ') .or.
     &         (p_anam(i).eq.' O   ') .or. (p_anam(i).eq.' N   ') )then
              local_at(1,i)=p_xyz(1,i)
              local_at(2,i)=p_xyz(2,i)
              local_at(3,i)=p_xyz(3,i)
            endif
          else
            local_at(1,i)=p_xyz(1,i)
            local_at(2,i)=p_xyz(2,i)
            local_at(3,i)=p_xyz(3,i)
          endif
c
c- Map into pixel space
c
          call o_to_j3(local_at(1,i),f_ind_n,ind_n,map2)

          do ii=ind_n(1)-a_pad(1), ind_n(1)+a_pad(1)
            do jj=ind_n(2)-a_pad(2), ind_n(2)+a_pad(2)
              do kk=ind_n(3)-a_pad(3), ind_n(3)+a_pad(3)
c
c---chech whether pixel is within a_rad
                iii(1)=ii
                iii(2)=jj
                iii(3)=kk
                call j3_to_f(iii,x1,map2)
                call f_to_o(x1,x,map2)
c
                sqdist=( x(1)-local_at(1,i) )**2 +
     &                 ( x(2)-local_at(2,i) )**2 +
     &                 ( x(3)-local_at(3,i) )**2
c
c
                if(sqdist.lt.ar_sq)then
                  call pp3(map2,ii,jj,kk,vala)
                  n_reset=n_reset+1
                endif
              enddo
            enddo
          enddo
c
        enddo
        write(6,*)'ENVELOPE CALCULATED,',n_reset,' pixels were reset'
        if(record)then
          write(nrec,*)
     &      'ENVELOPE CALCULATED,',n_reset,' pixels were reset'
        endif
c
        return
      endif
c
      if(nopt.ne.4) goto 40
c
      n_reset=0
c
c---- now skip thro' resetting
c---- loop over bricks
c
      do iz=1,n_brick(3,map2)
c
        jz=(iz-1)*brick(3,map2)
c
        do iy=1,n_brick(2,map2)
c
          jy=(iy-1)*brick(2,map2)
c
          do ix=1,n_brick(1,map2)
c
            jx=(ix-1)*brick(1,map2)
c
c---- now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map2)
c
              if(jjz.le.nx(3,map2))then
c
                do jjy=jy+1,jy+brick(2,map2)
c
                  if(jjy.le.nx(2,map2))then
c
                    do jjx=jx+1,jx+brick(1,map2)
c               
                      if(jjx.le.nx(1,map2))then
c
c---- ok this is a valid pixel, now test all logical operators...
c
      iindex(1)=jjx
      iindex(2)=jjy
      iindex(3)=jjz
c
      call j3_to_f(iindex,x1,map2)
      call f_to_o(x1,x,map2)
c
c---- test planes
c
c
c---- Can save some time here. If the current value is equal to the value it
c---- may be set to, we don't need to bother. JMD
c
      if(nplan.le.0)goto 800
      if(plan_v(1).eq.gp3_a(map2,iindex))goto 800
c
      if(.not.plan_or)then
c
        do jplan=1,nplan
          do j=1,2
            del_p(j)=0.0
            do i=1,3
              del_p(j)=del_p(j)+(plan(i,j,jplan)-x(i))**2
            enddo
          enddo
c
c---- More time saving - once its failed a test we can ditch out
c
          if(del_p(1).lt.del_p(2))goto 600
        enddo
        call pp3(map2,jjx,jjy,jjz,plan_v(1))
        n_reset=n_reset+1
c
600     continue
c
      else
c
        do jplan=1,nplan
          do j=1,2
            del_p(j)=0.0
            do i=1,3
              del_p(j)=del_p(j)+(plan(i,j,jplan)-x(i))**2
            enddo
          enddo
c
c---- More time saving - once its passed a test we can ditch out
c
          if(del_p(1).lt.del_p(2))then
            call pp3(map2,jjx,jjy,jjz,plan_v(1))
            n_reset=n_reset+1
            goto 700
          endif
        enddo
c
700     continue
c
      endif
c---- end plane timesaver
800   continue
c
c---- test spheres
c
      do jcirc=1,ncirc
        del=0.0
        do inner=1,3
          del=del+(cen(inner,jcirc)-x(inner))**2
        enddo
        if(del.le.rad(jcirc).and.iinout(jcirc).eq.0)then
          call pp3(map2,jjx,jjy,jjz,val(jcirc))
          n_reset=n_reset+1
        endif
        if(del.gt.rad(jcirc).and.iinout(jcirc).eq.1)then
          call pp3(map2,jjx,jjy,jjz,val(jcirc))
          n_reset=n_reset+1
        endif
      enddo
c
c
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c
      write(6,*)'Envelope calculated,',n_reset,' pixels were reset'
      if(record)then
        write(nrec,*)'Envelope calculated,',n_reset,
     &     ' pixels were reset'
      endif
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE map_atom
c     ===================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- local variables
c
      integer   map1,env1
      logical   i_check
      integer*2 irho
      external  gp3,gp3_a
      integer*2 gp3,gp3_a
c
      real      aocc,b
      integer   lun1,nres_local
      integer   map_cutoff,map_step(3)
      integer   iz,jz,iy,jy,ix,jx,jjz,jjy,jjx,ngadd
c
      real x1(3), x(3)
      integer iii(3)

c
c      equivalence  (i_temp,ch_temp)
c
c
c
      call map_atom_setup(map1,env1,i_check,map_cutoff,map_step,lun1)
c
c
c---- firstly set up pointers for the maps that drive the process 
c
      call init_map(map1)
c
c---- initialize r-factors, corr coefs, etc
c
      ngadd=0
      nres_local=0
c
      call f_stat
c
      if(i_check)call init_map(env1)
c
c
      call init_map(map1)
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
      do iz=1,n_brick(3,map1)
c
        jz=(iz-1)*brick(3,map1)
c
        do iy=1,n_brick(2,map1)
c
          jy=(iy-1)*brick(2,map1)
c
          do ix=1,n_brick(1,map1)
c
            jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
            do jjz=jz+1,jz+brick(3,map1),map_step(3)
c
            if(jjz.le.nx(3,map1))then
c
              do jjy=jy+1,jy+brick(2,map1),map_step(2)
c
              if(jjy.le.nx(2,map1))then
c
                do jjx=jx+1,jx+brick(1,map1),map_step(1)
c               
                if(jjx.le.nx(1,map1))then
c
c---- check that this pixel is in required region of the output map
c
      if(i_check)then
        if(gp3(env1,jjx,jjy,jjz).eq.0)goto 200
      endif
c
      irho=gp3(map1,jjx,jjy,jjz)
c
      ngadd=ngadd+1
c
      call f_a_stat(irho)
c
c
c
c
      if(gp3(map1,jjx,jjy,jjz).gt.map_cutoff)then
c
        iii(1)=jjx
        iii(2)=jjy
        iii(3)=jjz
c
        call j3_to_f(iii,x1,map1)
        call f_to_o(x1,x,map1)
c
        aocc=1.0
        b=20.0
        nres_local=nres_local+1
        write(lun1,1111)nres_local,x,aocc,b
1111    format('ATOM         OH2 WAT ',I5,4X,3F8.3,2f6.2)


      endif
c
200                     continue
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
c
c
c
c
c999   continue
      if(verbose)then
        write(6,998)ngadd,nres_local
        if(record)write(nrec,998)ngadd,nres_local
998     format(/,' MAP TO ATOM FINISHED,',/,
     &           ' ====================',/,
     &           i9,' pixels were inspected',/,
     &           i9,' atoms written')
      else
        write(6,997)ngadd,nres_local
        if(record)write(nrec,997)ngadd,nres_local
997     format(' MAP TO ATOM FINISHED,',i9,' pixels inspected',i9,
     &      ' atoms written')
      endif
c
c---- print some stats
c
      call f_p_stat(.true.)
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE map_atom_setup(map1,env1,i_check,map_cutoff,map_step,
     &  lun1)
c     ==================================================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c--- local variables
c
        integer      map1,env1,i,intval
        character*8  imap,ienv
c
        character*8  G2CHRU,ch_temp
        character*80 fil,charval
        integer      lun1
        logical      lass_get,i_check
        integer      map_cutoff,map_step(3)
c
c
        if(.not.lass_get(ch_temp,'MAPI'))then
c
          ch_temp=G2CHRU('AV_M2AT:INPUT_MAP',1,8,0)
c
          call lass_set(ch_temp,'MAPI')
        endif
          imap=ch_temp
c
c
        if(.not.lass_get(ch_temp,'ENVI'))then
c
        ch_temp=G2CHRU('AV_M2AT:INPUT_ENV',1,8,0)
c
          call lass_set(ch_temp,'ENVI')
        endif
        ienv=ch_temp
c
c
        if(verbose)then
          write(6,110) imap,ienv
          if(record)write(nrec,110) imap,ienv
110       format(/' CONVERT ELECTRON DENSITY TO ATOMS',/,
     &            ' ================================='
     &      /' Input map:       ',a
     &      /' Input envelope:  ',a)
          else
          write(6,111) imap,ienv
          if(record)write(nrec,111) imap,ienv
111       format(' MAP TO ATOMS',
     &    '  I/P map:',a,'  I/P env:',a)
          endif
C
        i_check=ienv.ne.'OFF'
c
c
        if(verbose)then
          if(i_check) write(6,120)
          if(.not.i_check)write(6,130)
c
          if(record)then
            if(i_check) write(nrec,120)
            if(.not.i_check)write(nrec,130)
            endif
          endif
c
120     format(' Input map envelope filtered')
130     format(' Input map not filtered')
c
c
c       get slot numbers for all maps
c
        call av_ass(map1,imap)
c----map 1 should exist by now !
        if(.not.defined(map1)) then
          write(6,140)
          if(record)write(nrec,140)
140       format(' %M2AT_SETUP-ERR: Input map empty')
          return
          endif
c---- check status of output file
c
c---- set up allocations for envelopes
c
        if(i_check) call av_ass(env1,ienv)
c
        lun1=1
c
        fil=charval('AV_MAP_ATOM: ENTER_FILE',1,80,0)
          if(index(fil,'.').eq.0) then
            i=index(fil,' ')
            fil=fil(1:i-1)//'.pdb'
            endif
c
        if(.not.terse)then
          write(6,40)fil(1:60), lun1
          if(record)write(nrec,40)fil(1:60), lun1
40        format(' Pdbfile:', a60,/,' opened on unit=',i3)
          else
          write(6,41)fil(1:60)
          if(record)write(nrec,41)fil(1:60)
41        format(' Pdbfile:', a60)
          endif
c
        close (unit=lun1)
        open  (unit=lun1,file=fil,status='new')
c
        map_cutoff=intval('AV_MAP_ATOM: ENTER_CUTOFF',-99999,99999,0)
c
        map_step(1)=intval('AV_MAP_ATOM: ENTER_STEP_X',1,1000,0)
        map_step(2)=intval('AV_MAP_ATOM: ENTER_STEP_Y',1,1000,0)
        map_step(3)=intval('AV_MAP_ATOM: ENTER_STEP_Z',1,1000,0)
c
        write(6,50) map_cutoff, map_step
        if(record)write(nrec,50) map_cutoff, map_step
50      format(' Map cutoff and steps in x, y and z:',4I7)
c
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE map_scale
c       ====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        real realval
c
                sca_map=realval('AV_MAP_SCALE:',-100000.0,100000.0,0.0)
c
        if(record)write(nrec,*)'MAP_SCALE changed to ',sca_map
c
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE map_shift
c       ====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        real realval
c
                shi_map=realval('AV_MAP_SHIFT:',-100000.0,100000.0,0.0)
c
        write(6,*)'MAP_SHIFT changed to ',shi_map
        if(record)write(nrec,*)'MAP_SHIFT changed to ',shi_map
c
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE map_shift_sol
c     ========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      integer     map
      character*8 g2chru,lmap
c
      lmap=g2chru('AV_MAP_SHIFT:ENTER_MAP_NAME',1,8,0)
      call av_ass(map,lmap)
      shi_map=sollev(map)
c
      write(6,*)'MAP_SHIFT changed to ',shi_map
      if(record)write(nrec,*)'MAP_SHIFT changed to ',shi_map
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE map_shift_true
c       ====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        real            realval
        character*8     G2CHRU,ch_temp
        integer         map1
        logical         lass_get
c
        if(.not.lass_get(ch_temp,'MAPI'))then
          ch_temp=G2CHRU('AV_SIG_CUT:INPUT_MAP',1,8,0)
          call lass_set(ch_temp,'MAPI')
        endif
        call av_ass(map1,ch_temp)
c
        shi_map=realval('AV_MAP_SHIFT:',-100000.0,100000.0,0.0)
        shi_map=shi_map*mscale(map1)
c
        write(6,*)'MAP_SHIFT changed to ',shi_map
        if(record)write(nrec,*)'MAP_SHIFT changed to ',shi_map
c
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE MATM33(A,B,C)
c       ========================
      IMPLICIT NONE
c
c*******************************************************************************
        integer i,j,k
        real a,b,c
c
        DIMENSION A(3,3),B(3,3),C(3,3)
        DO 100 I=1,3
        DO 100 J=1,3
        C(I,J)=0.0
        DO 100 K=1,3
        C(I,J)=C(I,J)+A(I,K)*B(K,J)
  100   CONTINUE
        RETURN
        END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE MATMUL (A,B,C)
c     =========================
      IMPLICIT NONE
c
c     MULTIPLY TWO 3X3 MATRICES (A X B = C)
c
      INTEGER   I,J,K
      REAL      A,B,C
      REAL*8    SUM_R
      DIMENSION A(3,3),B(3,3),C(3,3)
c
      DO  I=1,3
        DO J=1,3
          SUM_R=0.0D0
          DO K=1,3
            SUM_R=SUM_R+DPROD( A(I,K), B(K,J) )
          ENDDO
          C(I,J)=SNGL(SUM_R)
        ENDDO
      ENDDO
C
      RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE matmul1(a,b,c,v,cg,vt)
c     =================================
      IMPLICIT NONE
c
C---- MULTIPLY TWO 3X3 MATRICES (A X B = C)
c---- returns a rotation matrix and new translation vector having applied
c---- a centre of gravity recentring.
c---- [R]=[ROT][OPS]
c----  t'= [ROT]t-[ROT]cg+cg
c
        real a,b,c,vt,v,cg,sum_r
        integer i,j,k
      dimension a(3,3),b(3,3),c(3,3),vt(3),cg(3),v(3)
c
      do  I=1,3
        do j=1,3
          sum_r=0.0
           do k=1,3
             sum_r=sum_r+a(i,k)*b(k,j)
           enddo
             c(i,j)=sum_r
        enddo
      enddo
c
        vt(1)=(v(1)*A(1,1)+V(2)*A(1,2)+V(3)*A(1,3))-
     &        (cg(1)*A(1,1)+cg(2)*A(1,2)+cg(3)*A(1,3))+cg(1)
        vt(2)=(V(1)*A(2,1)+V(2)*A(2,2)+V(3)*A(2,3))-
     &        (cg(1)*A(2,1)+cg(2)*A(2,2)+cg(3)*A(2,3))+cg(2)
        vt(3)=(V(1)*A(3,1)+v(2)*a(3,2)+v(3)*a(3,3))-
     &        (cg(1)*A(3,1)+cg(2)*A(3,2)+cg(3)*A(3,3))+cg(3)
      RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE MATMUL_D (A,B,C)
c     =========================
      IMPLICIT NONE
c
C     MULTIPLY TWO 3X3 MATRICES (A X B = C)
c
        integer i,j,k
        double precision a,b,c,sum_r,local
      DIMENSION A(3,3),B(3,3),C(3,3),local(3,3)
c
      DO  I=1,3
        do j=1,3
          sum_r=0.0
           do k=1,3
             sum_r=sum_r+A(i,k)*B(k,j)
           enddo
             local(i,j)=sum_r
        enddo
      enddo
C
      do i=1,3
        do j=1,3
          c(i,j)=local(i,j)
        enddo
        enddo
c
        RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE MATROT(RROT,RT1,RT2,RT3,RAXIS,MODE)
c     ==============================================
      IMPLICIT NONE
c
c
c--- routine cleaved from xplor
c--- uses d.p within routine
c
C
C Routine computes Eulerian angles (MODE="EULE"), Lattman angles
C (MODE="LATT"), spherical polar angles (MODE="SPHE"), or rotation 
C axis and angle (MODE="AXIS") from unitary matrix ROT.  The angular 
C definitions are described in routine ROTMAT.  
C
C Restrictions:
C    Matrix ROT has to be unitary, i.e., det(ROT)=+1 otherwise a fatal
C    warning will be issued. 
C
C Conventions:
C   In Eulerian angle mode T2 can be forced without restriction of generality 
C   to be located in the interval 0<= T2 <= pi (this is a consequence of 
C   the identity operation t1,t2,t3 -> pi+t1,-t2,pi+t3 in Eulerian angle space).
C   T1 is forced into 0 <= t1 < 2*pi and T3 is forced into 
C   0<=T3< 2*pi.  If T2 is 0 or PI we set T3 to zero without restriction
C   of generality.  
C
C   Lattman angles are forced into the intervals (0 <= theta- <= 2*pi,
C   0<= theta2 <= pi, 0 <= theta+ < 4*pi).  
C
C   For spherical polar angles we force psi (=t1) into 0<= t1 < pi and
C   phi (=t2) into 0<= t2 < pi, kappa (=t3) into 0 <= t3 < 2*pi
C   without restriction of generality.  If kappa is equal to zero,
C   all angles will be zero.  If psi is equal to zero, phi will be set to zero. 
C
C   Spherical polar angles are used to compute the rotation axis in AXIS mode.
C   Without restriction of generality, kappa is forced into 0 <= kappa <= pi.
C
C Input:
C    MODE specifies angle mode
C    ROT(3,3) contains the rotation matrix.  The matrix is defines
C    the rotation according to r'(i)=sum_j ROT(i,j)*r(j)
C
C Ouput: 
C    T1,T2,T3 are theta1 (z), theta2 (x'), theta3 (z') for MODE="EULE"
C    T1,T2,T3 are theta+, theta2, theta- for MODE="LATT"
C    T1,T2,T3 are psi (incl. vs. y), phi (azimuthal), kappa for MODE="SPHE"
C             and MODE="AXIS"
C    T3, AXIS(3) are kappa and a 3-D vector specifying the axis for MODE="AXIS"
C Note: all rotations are counter-clockwise
C
C Author: Axel T. Brunger
C =======================
C
C I/O
      INCLUDE 'average.fcm'
C
      REAL RROT(3,3),RT1,RT2,RT3,RAXIS(3),rrot2(3,3),pi
C
C
      DOUBLE PRECISION ROT(3,3), T1, T2, T3, AXIS(3)
      CHARACTER*4 MODE
C local
      DOUBLE PRECISION S1, S2, S3, C1, C2, C3, DET, C22, C12, TP,TM
      DOUBLE PRECISION ROT2(3,3)
      LOGICAL COND
      INTEGER I, J
C parameter
      DOUBLE PRECISION RAD, ONE, TWO, SMALL, R180, R360, SMALLR
      PARAMETER (RAD=0.017453292D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (SMALL=1.0D-2, R180=180.0D0, R360=360.0D0)
      PARAMETER (SMALLR=1.0D-8)
      PARAMETER (PI=3.141593)
C
      AXIS(1)=ZERO
      AXIS(2)=ZERO
      AXIS(3)=ZERO
C
C
C---- CONVERT TO REAL TO DO
      DO I=1,3
        DO J=1,3
          ROT(I,J)=RROT(I,J)
        ENDDO
      ENDDO
C
C
C test to make sure that the matrix is unitary
      DET=(ROT(1,1)*ROT(2,2)-ROT(1,2)*ROT(2,1))*ROT(3,3)
     &   +(ROT(2,1)*ROT(3,2)-ROT(2,2)*ROT(3,1))*ROT(1,3)
     &   +(ROT(3,1)*ROT(1,2)-ROT(3,2)*ROT(1,1))*ROT(2,3)
      IF (ABS(DET-ONE).GT.SMALL) THEN
        write(6,999)
        if(record)write(nrec,999)
999     format(' %MATROT-ERR: Matrix is not unitary')
        return
      ELSE
C==================================================================
        IF (MODE.EQ.'EULE'.OR.MODE.EQ.'LATT') THEN
C
C first, let's get t2:
          C2=dmax1(-ONE,dmin1(ONE,ROT(3,3)))
          T2=ACOS(C2)/RAD
C
C second, let's compute SIN(T2) (OE positive)
          S2=SQRT(dmax1(ZERO,ONE-ROT(3,3)**2))
C
C if COS(T2) not equal +-1 then T1 and T3 are uniquely determined
          IF (ABS(ABS(C2)-ONE).GT.SMALLR) THEN
            S3=ROT(1,3)/S2
            C3=ROT(2,3)/S2
            C3=dmax1(-ONE,dmin1(ONE,C3))
            T3=ACOS(C3)
            IF (S3.LT.ZERO) T3=TWO*PI-T3
            T3=T3/RAD
            S1=ROT(3,1)/S2
            C1=-ROT(3,2)/S2
            C1=dmax1(-ONE,dmin1(ONE,C1))
            T1=ACOS(C1)
            IF (S1.LT.ZERO) T1=TWO*PI-T1
            T1=T1/RAD
          ELSE
C
C without restriction T3 can be set to zero
            T3=ZERO
            C1=ROT(1,1)
            S1=ROT(1,2)
            C1=dmax1(-ONE,dmin1(ONE,C1))
            T1=ACOS(C1)
            IF (S1.LT.ZERO) T1=TWO*PI-T1
            T1=T1/RAD
          END IF
C
          IF (MODE.EQ.'LATT') THEN
C
C in "Lattman" mode compute theta+, theta- from theta1, theta3
            TP=T1+T3
            TM=T1-T3
            T1=TP
            T3=TM
C
C we know that 0 <= theta1 < 2*pi and 0 <= theta3 < 2*pi.  We have to
C project this into the Lattman space asymmetric unit (0<= theta+ < 4*pi,
C 0 <= theta- < 2*pi).  
            IF (T3.LT.ZERO) THEN
              IF (T1.LT.R360) THEN
                T1=T1+R360
                T3=T3+R360
              ELSE
                T1=T1-R360
                T3=T3+R360
              END IF
            END IF
C
          END IF
C==========================================================
        ELSE IF (MODE.EQ.'SPHE'.OR.MODE.EQ.'AXIS') THEN
C
C first lets get COS (kappa) 
          C3=(ROT(1,1)+ROT(2,2)+ROT(3,3)-ONE)/TWO
          C3=dmax1(-ONE,dmin1(ONE,C3))
C
C special case COS(kappa)=1
          IF (ABS(C3-ONE).LT.SMALLR) THEN
            T1=ZERO
            T2=ZERO
            T3=ZERO
          ELSE
C
C determine COS^2(psi)
            C12=dmax1(ZERO,(ROT(2,2)-C3)/(ONE-C3))
C
C special case COS(psi)=+-1
            IF (ABS(C12-ONE).LT.SMALLR) THEN
              T1=ZERO
              T2=ZERO
              T3=ACOS(C3)
              IF (ROT(3,1).LT.ZERO) T3=TWO*PI-T3
              T3=T3/RAD
            ELSE
C
C determine COS^2(phi)
              C22=dmax1(ZERO, (ROT(1,1)-C3)/((ONE-C3)*(ONE-C12)) )
C
C special case COS(phi)=+-1
              IF (ABS(C22-ONE).LT.SMALLR) THEN
                C1=SIGN(ONE,(ROT(1,2)+ROT(2,1))/(ONE-C3)  )*SQRT(C12)
                C1=dmax1(-ONE,dmin1(ONE,C1))
                T1=ACOS(C1)/RAD
                T2=ZERO
                T3=ACOS(C3)
                IF (ROT(2,3)-ROT(3,2).LT.ZERO) T3=TWO*PI-T3
                T3=T3/RAD
              ELSE
C
C determine sign of COS(psi) and then determine psi 
                C1=SIGN(ONE,(-ROT(2,3)-ROT(3,2))/(ONE-C3)  )*SQRT(C12)
                C1=dmax1(-ONE,dmin1(ONE,C1))
                T1=ACOS(C1)/RAD
C
C determine sign of COS(phi) and then determine phi
                C2=SIGN(ONE,(-ROT(1,3)-ROT(3,1))/(ONE-C3)  )*SQRT(C22)
                C2=dmax1(-ONE,dmin1(ONE,C2))
                T2=ACOS(C2)/RAD
C
C determine kappa
                T3=ACOS(C3)
                IF (ABS(C2-ONE).LT.SMALL) THEN
                  IF (ROT(2,3)-ROT(3,2).LT.ZERO) T3=TWO*PI-T3
                ELSE
                  IF (ROT(1,2)-ROT(2,1).GT.ZERO) T3=TWO*PI-T3
                END IF
                T3=T3/RAD
              END IF
            END IF
          END IF
C
          IF (MODE.EQ.'AXIS') THEN
C
            AXIS(2)=COS(T1*RAD)
            AXIS(1)=COS(T2*RAD)*SQRT(dmax1(ZERO,ONE-AXIS(2)**2))
            AXIS(3)=-SQRT(dmax1(ZERO,ONE-AXIS(1)**2-AXIS(2)**2))
C
C without restriction of generality we force T3 in 0 <= T3 <= pi
            IF (T3.GT.R180) THEN
              T3=TWO*R180-T3
              AXIS(1)=-AXIS(1)
              AXIS(2)=-AXIS(2)
              AXIS(3)=-AXIS(3)
            END IF
          END IF
C
        END IF
C
        rt1=t1
        rt2=t2
        rt3=t3
c
        raxis(1)=axis(1)
        raxis(2)=axis(2)
        raxis(3)=axis(3)
c
C now back-compute the matrix as an internal consistency check
        CALL ROTMAT(RROT2,RT1,RT2,RT3,RAXIS,MODE)
        DO I=1,3
          DO J=1,3
            ROT2(I,J)=RROT2(I,J)
          ENDDO
        ENDDO
c
        COND=
     &       ABS(ROT(1,1)-ROT2(1,1)).GT.SMALL
     &   .OR.ABS(ROT(1,2)-ROT2(1,2)).GT.SMALL
     &   .OR.ABS(ROT(1,3)-ROT2(1,3)).GT.SMALL
     &   .OR.ABS(ROT(2,1)-ROT2(2,1)).GT.SMALL
     &   .OR.ABS(ROT(2,2)-ROT2(2,2)).GT.SMALL
     &   .OR.ABS(ROT(2,3)-ROT2(2,3)).GT.SMALL
     &   .OR.ABS(ROT(3,1)-ROT2(3,1)).GT.SMALL
     &   .OR.ABS(ROT(3,2)-ROT2(3,2)).GT.SMALL
     &   .OR.ABS(ROT(3,3)-ROT2(3,3)).GT.SMALL
        IF (COND) THEN
          WRITE(6,'(A,3G12.4)')
     &      ' %ROTMAT-ERR: Inconsistent, T=',T1, T2, T3
          WRITE(6,'(/A,3(/3F12.6))')
     &      ' ROT =',((ROT(I,J),J=1,3),I=1,3)
          WRITE(6,'(/A,3(/3F12.6))')
     &      ' ROT2 =',((ROT2(I,J),J=1,3),I=1,3)
          IF (RECORD) THEN
            WRITE(nrec,'(A,3G12.4)')
     &        ' %ROTMAT-ERR: Inconsistent, T=',T1, T2, T3
            WRITE(nrec,'(/A,3(/3F12.6))')
     &        ' ROT =',((ROT(I,J),J=1,3),I=1,3)
            WRITE(nrec,'(/A,3(/3F12.6))')
     &        ' ROT2 =',((ROT2(I,J),J=1,3),I=1,3)
          ENDIF
        END IF
C
      END IF
c
C
      RETURN
      END
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE max_peak(map1,max_rho,jstart)
c      =========================================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c      =====================
c
c--- Local variables:
c
       integer   map1,ngadd,iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
       integer*2 irho,max_rho
       external  gp3,gp3_a
       integer*2 gp3,gp3_a
       integer   jstart(3)
c
c--- Set up pointers for the maps that drive the process
c
       call init_map(map1)
c
c--- Zero counters
c
       ngadd=0
       max_rho=-32000
c
c
c
c--- Off to work...
c--- =============
c
c--- Loop over bricks
c
       do iz=1,n_brick(3,map1)
c
         jz=(iz-1)*brick(3,map1)
c
         do iy=1,n_brick(2,map1)
c
           jy=(iy-1)*brick(2,map1)
c
           do ix=1,n_brick(1,map1)
c
             jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
             do jjz=jz+1,jz+brick(3,map1)
c
             if(jjz.le.nx(3,map1))then
c
               do jjy=jy+1,jy+brick(2,map1)
c
               if(jjy.le.nx(2,map1))then
c
                 do jjx=jx+1,jx+brick(1,map1)
c                    
                 if(jjx.le.nx(1,map1))then
cc
                   irho=gp3(map1,jjx,jjy,jjz)
c
                   if(irho.gt.max_rho)then
                       max_rho=irho
                       jstart(1)=jjx
                       jstart(2)=jjy
                       jstart(3)=jjz
                       endif
c
c200                continue
                 endif
                 enddo
               endif
               enddo
             endif
             enddo
           enddo
         enddo
       enddo
c
c--- OK
c
c300    continue
c
c--- And it is all over...
c
       return
       end
C
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE mod_csy
c       ==================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        character*8 G2CHRU,ch_temp
c
c----update cell information
c
        integer mapn
        character*8 mapl
c
        ch_temp=G2CHRU('AV_MCSY:ENTER_MAP_NAME',1,8,0)
c
        mapl=ch_temp
c
        call av_ass(mapn,mapl)
c
        write(6,*)' UPDATE CRYSTALLOGRAPHIC SYMM SWITCH FOR MAP:',mapn
        if(record)write(nrec,*)
     &            ' UPDATE CRYSTALLOGRAPHIC SYMM SWITCH for map:',mapn
c
100     ch_temp=G2CHRU('AV_MCSY:ENTER_ON_or_OFF',1,8,0)
c
        if(ch_temp.eq.'ON')then
           write(6,*)' cryst symm switched  ON'
           if(record)write(nrec,*)' cryst symm switched  ON'
           xr_con(mapn)=.true.
c
           else
c
             if(ch_temp.eq.'OFF')then
               write(6,*)' cryst symm switched  OFF'
               if(record)write(nrec,*)' cryst symm switched  OFF'
               xr_con(mapn)=.false.
             else
                goto 100
             endif
          endif
c
c
c
        return
        end
C 
c
c==============================================================================
c
c all of rob esnoufs code below here
C **********************************
c
C=============================================================================

      SUBROUTINE NFOLDNCS(NNCS, ROT, TRAN, ROTI, TRANI, IERR)
      IMPLICIT NONE

C Idealise an NNCS-fold NCS rotation from the set of transformations
C The first operator should be the identity, the rest non-singular
C rotation matrices. The idealised rotations and translations are returned
C in ROTI and TRANI, which may be the same arrays as ROT and TRAN, in which
C case the input operators are overwritten.

C IERR > 1: That rotation matrix is singular
C IERR = 1: First operator was not identity
C IERR = 0: No error
C IERR = -1: Need to increase array sizes in subroutine
C IERR = -2: Error when sorting rotation angles

      INTEGER NNCS
      DOUBLE PRECISION ROT(3,3,NNCS)
      DOUBLE PRECISION TRAN(3,NNCS)
      DOUBLE PRECISION ROTI(3,3,NNCS)
      DOUBLE PRECISION TRANI(3,NNCS)
      INTEGER IERR

C Parameters for subroutine

      DOUBLE PRECISION RTOD
      PARAMETER (RTOD = 57.295779512)
      DOUBLE PRECISION DTOR
      PARAMETER (DTOR = 1.0/RTOD)
      INTEGER MAXNCS
      PARAMETER (MAXNCS=24)

C Comparision function for angle sort routine

      INTEGER COMP
      EXTERNAL COMP

C Work arrays... must be big enough for the highest level of NCS

      DOUBLE PRECISION VEC(3,MAXNCS)
      DOUBLE PRECISION ANG(MAXNCS), ANG2(MAXNCS)
      DOUBLE PRECISION TRAN2(3,MAXNCS)
      DOUBLE PRECISION AVVEC(3), ORIG(3)

      INTEGER NCS
      INTEGER IX(MAXNCS), STACK(MAXNCS)
      DOUBLE PRECISION THETA, DOT, SN, CS

      COMMON /ANGS/ ANG
      SAVE /ANGS/

C If the order of the NCS is too big then return with IERR = -1

      IF (NNCS.GT.MAXNCS) THEN
        IERR = -1
        RETURN
      ENDIF

C If the first operator isn't the identity the return with IERR = 1

      IF (ROT(1,1,1).NE.1.0 .OR. ROT(2,2,1).NE.1.0 .OR.
     &    ROT(3,3,1).NE.1.0 .OR.
     &    ROT(2,1,1).NE.0.0 .OR. ROT(3,1,1).NE.0.0 .OR.
     &    ROT(3,2,1).NE.0.0 .OR.
     &    ROT(1,2,1).NE.0.0 .OR. ROT(1,3,1).NE.0.0 .OR.
     &    ROT(2,3,1).NE.0.0 .OR.
     &    TRAN(1,1).NE.0.0 .OR. TRAN(2,1).NE.0.0 .OR.
     &    TRAN(3,1).NE.0.0 ) THEN
        IERR = 1
        RETURN
      ENDIF

C Orthogonalise the rotation matrices in case they are slightly approximate

      DO NCS = 2, NNCS
        CALL GSORTH(ROT(1,1,NCS),IERR)
        IF (IERR.LT.0) THEN
          IERR = NCS
          RETURN
        ENDIF
      ENDDO

C Convert the rotation matrices into a rotation around an axis. (Which
C in our convention passes through the origin.)

      DO NCS = 2, NNCS
        CALL ROTAXIS(ROT(1,1,NCS), VEC(1,NCS), ANG(NCS))
      ENDDO

C Select signs for all the axes so that the they are all within 90 degrees
C of the first axis. ie if the dot product is negative then swap the signs
C of the direction vectors and the rotation angles.

      DO NCS = 3, NNCS
        IF (VEC(1,2)*VEC(1,NCS) + VEC(2,2)*VEC(2,NCS) +
     &      VEC(3,2)*VEC(3,NCS).LT.0.0) THEN
          VEC(1,NCS) = -VEC(1,NCS)
          VEC(2,NCS) = -VEC(2,NCS)
          VEC(3,NCS) = -VEC(3,NCS)
          ANG(NCS) = -ANG(NCS)
        ENDIF
      ENDDO

C Get the average direction vector... this could be weighted according to
C amount of rotation, but isn't

      AVVEC(1) = 0.0
      AVVEC(2) = 0.0
      AVVEC(3) = 0.0
      DO NCS = 2, NNCS
        AVVEC(1) = AVVEC(1) + VEC(1,NCS)
        AVVEC(2) = AVVEC(2) + VEC(2,NCS)
        AVVEC(3) = AVVEC(3) + VEC(3,NCS)
      ENDDO
      CALL NORM(AVVEC, AVVEC)

C Make sure that a rotation of approximately 180 degrees is +180ish

      DO NCS = 2, NNCS
        IF (ANG(NCS).LT.(-360.0/NNCS)*(((NNCS+1)/2)-0.5)*DTOR ) THEN
          ANG(NCS) = ANG(NCS) + 360.0*DTOR
        ENDIF
      ENDDO

C Now sort the angles into order so that we can match the desired angles
C against the actual ones

      DO NCS = 1, NNCS
        IX(NCS) = NCS
      ENDDO
      CALL QSORTP(IX, NNCS, STACK, MAXNCS, COMP, IERR)
      IF (IERR.NE.0) THEN
        IERR = -2
        RETURN
      ENDIF

C Now store the ideaelised angles with the appropriate approximate one

      DO NCS = 1, NNCS
        THETA = (360.0/NNCS)*(NCS-((NNCS+1)/2))
        ANG2(IX(NCS)) = THETA*DTOR
      ENDDO

C Project the translation vectors into the plane perpendicular to the averaged
C rotation axis and passing through the origin

      DO NCS = 2, NNCS
        DOT = TRAN(1,NCS)*AVVEC(1) + TRAN(2,NCS)*AVVEC(2) +
     &        TRAN(3,NCS)*AVVEC(3)
        TRAN2(1,NCS) = TRAN(1,NCS) - AVVEC(1) * DOT
        TRAN2(2,NCS) = TRAN(2,NCS) - AVVEC(2) * DOT
        TRAN2(3,NCS) = TRAN(3,NCS) - AVVEC(3) * DOT
      ENDDO

C Average (all) the translation vectors to get the rotation centre

      ORIG(1) = 0.0
      ORIG(2) = 0.0
      ORIG(3) = 0.0
      DO NCS = 2, NNCS
        ORIG(1) = ORIG(1) + TRAN2(1,NCS)
        ORIG(2) = ORIG(2) + TRAN2(2,NCS)
        ORIG(3) = ORIG(3) + TRAN2(3,NCS)
      ENDDO
      ORIG(1) = ORIG(1)/NNCS
      ORIG(2) = ORIG(2)/NNCS
      ORIG(3) = ORIG(3)/NNCS

C Summarize what we have found out

      DO NCS = 1, NNCS
        WRITE(6,*)
        WRITE(6,'(A,I2)') 'Operator:', NCS
        WRITE(6,*)
        WRITE(6,'(2X,3F9.5,F16.5)') ROT(1,1,NCS), ROT(2,1,NCS),
     &      ROT(3,1,NCS), TRAN(1,NCS)
        WRITE(6,'(2X,3F9.5,F16.5)') ROT(1,2,NCS), ROT(2,2,NCS),
     &      ROT(3,2,NCS), TRAN(2,NCS)
        WRITE(6,'(2X,3F9.5,F16.5)') ROT(1,3,NCS), ROT(2,3,NCS),
     &      ROT(3,3,NCS), TRAN(3,NCS)
        WRITE(6,*)
        WRITE(6,'(A,3F9.5)') 'Direction vector',
     &                          VEC(1,NCS), VEC(2,NCS), VEC(3,NCS)
        WRITE(6,'(A,F8.2)') 'Rotation angle', ANG(NCS)*RTOD
      ENDDO

C Display the angles between the average vector and the individual vectors

      WRITE(6,*)
      WRITE(6,'(A,3F9.5)') 'Averaged vector:', AVVEC(1), AVVEC(2),
     &    AVVEC(3)
      WRITE(6,*)
      DO NCS = 2, NNCS
        THETA = ACOS(AVVEC(1)*VEC(1,NCS) + AVVEC(2)*VEC(2,NCS) +
     &            AVVEC(3)*VEC(3,NCS))
        WRITE(6,'(A,I2,A,F7.2,A)') 'Change to rotation axis', NCS,
     &            ' is', THETA*RTOD, ' degrees'
      ENDDO

      WRITE(6,*)
      DO NCS = 1, NNCS
        WRITE(6,'(A,I2,A,F7.2,A,F7.2,A)') 'Rotation', NCS,
     &      ' changed from',
     &      ANG(NCS)*RTOD, ' to', ANG2(NCS)*RTOD, ' degrees'
      ENDDO

      WRITE(6,*)
      WRITE(6,'(A,3F9.3)') 'Origin:', ORIG(1), ORIG(2), ORIG(3)

C And finally create the idealised transformations and return

      DO NCS = 1, NNCS
        SN = SIN(ANG2(NCS))
        CS = COS(ANG2(NCS))
        ROTI(1,1,NCS) = CS + AVVEC(1)**2 * (1.0-CS)
        ROTI(2,1,NCS) = -AVVEC(3)*SN + AVVEC(1)*AVVEC(2)*(1.0-CS)
        ROTI(3,1,NCS) =  AVVEC(2)*SN + AVVEC(1)*AVVEC(3)*(1.0-CS)
        ROTI(1,2,NCS) =  AVVEC(3)*SN + AVVEC(1)*AVVEC(2)*(1.0-CS)
        ROTI(2,2,NCS) = CS + AVVEC(2)**2 * (1.0-CS)
        ROTI(3,2,NCS) = -AVVEC(1)*SN + AVVEC(2)*AVVEC(3)*(1.0-CS)
        ROTI(1,3,NCS) = -AVVEC(2)*SN + AVVEC(1)*AVVEC(3)*(1.0-CS)
        ROTI(2,3,NCS) =  AVVEC(1)*SN + AVVEC(2)*AVVEC(3)*(1.0-CS)
        ROTI(3,3,NCS) = CS + AVVEC(3)**2 * (1.0-CS)
        TRANI(1,NCS) = - ROTI(1,1,NCS)*ORIG(1) - ROTI(2,1,NCS)*ORIG(2)
     &                 - ROTI(3,1,NCS)*ORIG(3) + ORIG(1)
        TRANI(2,NCS) = - ROTI(1,2,NCS)*ORIG(1) - ROTI(2,2,NCS)*ORIG(2)
     &                 - ROTI(3,2,NCS)*ORIG(3) + ORIG(2)
        TRANI(3,NCS) = - ROTI(1,3,NCS)*ORIG(1) - ROTI(2,3,NCS)*ORIG(2)
     &                 - ROTI(3,3,NCS)*ORIG(3) + ORIG(3)
      ENDDO
      IERR = 0

      RETURN
      END

C=============================================================================

      SUBROUTINE NORM(VEC1, VEC2)

      IMPLICIT NONE
C Normalise the second vector into the first

        DOUBLE PRECISION VEC1(3), VEC2(3)

        DOUBLE PRECISION LEN

        LEN = VEC2(1)**2 + VEC2(2)**2 + VEC2(3)**2
        IF (LEN.LE.0.0) THEN
          VEC1(1) = 0.0
          VEC1(2) = 0.0
          VEC1(3) = 0.0
        ELSE
          LEN = SQRT(LEN)
          VEC1(1) = VEC2(1) / LEN
          VEC1(2) = VEC2(2) / LEN
          VEC1(3) = VEC2(3) / LEN
        ENDIF

        RETURN
        END
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE normal(j,r_rot,r_vec)
c       ================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c
      real    coor(3,5), xa(3,5), xb(3,5)
      real    r_rot(3,3), r_vec(3)
      real    rmse
      integer i, j, l, k
c
      call rigset
c
c--- vectors to apply ops tp
c
      coor(1,1)=0.
      coor(2,1)=0.
      coor(3,1)=0.
c
      coor(1,2)=1.
      coor(2,2)=0.
      coor(3,2)=0.
c
      coor(1,3)=0.
      coor(2,3)=1.
      coor(3,3)=0.
c
      coor(1,4)=0.
      coor(2,4)=0.
      coor(3,4)=1.
c
      coor(1,5)=1.
      coor(2,5)=1.
      coor(3,5)=1.
c       
c       
      do i=1,5
        call vmtply(save_r(1,1,j),coor(1,i),xb(1,i))
        call vmtply(ops(1,1,j),coor(1,i),xa(1,i))
        call rigadd(xa(1,i),xb(1,i),1.0)
      enddo
c               
      call rigsol(r_rot,r_vec,rmse,.true.)
c
      write(6,990)j
      if(record)write(nrec,990)j
990   format(' For sym ops ',i4)
      write(6,1000)  rmse,((r_rot(l,k),k=1,3),l=1,3), r_vec
      if(record)write(nrec,1000)
     &                 rmse,((r_rot(l,k),k=1,3),l=1,3), r_vec
1000  format(' After Kabsch routine,  RMS error =',f10.4,/,
     &           ' matrix:',t20,3f10.5,2(/,t20,3f10.5),/,
     &           ' vector:',t20,3f10.2)
c
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE o_to_j3(xo,aj3,j3,map)
c       =================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        real    xo(3),xlocal(3),aj3(3)
        integer map,j3(3)
c
        call vmtply(xrtr(1,1,map),xo,xlocal)
        call f_to_j3(xlocal,aj3,j3,map)
c
        return
        end
c



c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      SUBROUTINE ortho_mats
c       =====================
      IMPLICIT NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      INCLUDE 'average.fcm'
c       =====================
c
c
c subroutine to call SUBROUTINE GSORTH(MATRIX, ERRCOD). 
c MATRIX is a double precision 3x3 matrix. ERRCOD integer error code.
c and orthogonalise matrix to ensure pure rotational components
c
c
        integer intval,nncs,j,x,y,ierr,nsym1,nsym2
        DOUBLE PRECISION dbROT(3,3,maxsym)
c
c
        nncs=0
        ierr=0
c
        if (nsym.lt.1) then
          write(6,10)
            if(record) write(nrec,10)
          return
        endif
10      format('AV_OMAT_ERR: No Symmetry Operators defined')
c
c
        nsym1=intval('AV_OMAT_SYM:ENTER_1_NUMBER',1,nsym,0)
        nsym2=intval('AV_OMAT_SYM:ENTER_2_NUMBER',nsym1,nsym,0)
c
c         
c
          write(6,20) nsym1,nsym2
          if(record) write(nrec,20) nsym1,nsym2
20      format(' Symmetry operators to orthogonalise',i4,' to ',i4)
c
c
        if(verbose)then
          do j=nsym1,nsym2
            call prtsym(6,j)
          enddo
        endif
c
c
c-- at mo copy real ops into double precision 
c-- do we want to make these double precision ??
c
c-- would use subroutine copy by default but it is not double precision
c


        do j=nsym1,nsym2
        nncs=nncs+1
        do x=1,3
        do y=1,3
          dbrot(x,y,nncs) = ops(x,y,j)
        enddo
        enddo
        enddo
c
c-- could check that nncs = nsym2 - nsym1 + 1 
c
c
        do j=1,nncs
c       
        call GSORTH(dbrot(1,1,j), ierr)
c
        enddo
c
        if(ierr.ne.0) then
          write(6,30)
          if(record) write(nrec,30)
30        format(' AV_OMAT_ERR: Rotation matrix is singular')
          return
        endif
c

c
c Now write back new operators to ops defined as reals
c
        do j=nsym1,nsym2
        do x=1,3
        do y=1,3
          ops(x,y,j) = dbrot(x,y,j-nsym1+1)
        enddo
        enddo
        enddo

c
c
          write(6,60)
          if(record) write(nrec,60)
60        format(' AV_OMAT: Symmetry matrices orthogonalised')

        if(verbose)then
          do j=nsym1,nsym2
            call prtsym(6,j)
          enddo
        endif
c


        RETURN
        END
C 
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE outxpl
C     =================
      IMPLICIT NONE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
C       To create pdb files ready for output to xplor...using columns
c       73 to 76 as identifiers.
c       Will appply all ncs operators to generate up all molecules 
c       within the asymmetric unit....this must include the unitary 
c       matrix.  Will label A,B,C...initially to 13 elements
C
      INCLUDE 'average.fcm'
c     ====================
C
c   we use R. ESNOUF'S common block descriptions for pdb files, his coding
c   for reading pdb files, and to add a terminal oxygen...thanks.
c
c*******************************************************************************
      integer  lun2,nsym1,nsym2,nopt,keyword,j,i,l
      integer intval
      integer*4 iflag
c
      real*4  coor1(3,natom)
c
      logical inv
      character*80  file1
      character*80  fil
      character*4   flag1
      character*1   flag2(4)
c
      equivalence(flag1,flag2)
c
C========================================================================
C
c----lun is input stream for coords
      lun2=2
      iflag=0
      inv=.false.
      write(6,100)
      if(record)write(nrec,100)
100   format('CREATE PDB FILES SUITABLE FOR X-PLOR')
c
c
c---check whether any symm ops
c
      if(nsym.lt.1)then
        write(6,*)'%OUTXPL-ERR: No symm operators'
        if(record)write(nrec,*)'%OUTXPL-ERR: No symm operators'
        return
      endif
c
c
      call read_pdb(fil)
c
      call addoxt
      write(6,*) 'Terminal Oxygen added, no of atoms is',p_na
      if(record)write(nrec,*) 'Terminal Oxygen added, no of atoms is',
     &  p_na
c
      iflag=intval('ADD_CHAIN_IDENTIFIER:',1000,9000,0)
c
      nsym1=1
      nsym2=nsym
c
75    nopt=keyword(' SYMM_OPS','ALL !STAR!END !INV !GO  ',0)
c
      if(nopt.eq.2)nsym1=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.eq.3)nsym2=intval('AV_AVERAGE_SYMM',1,nsym,0)
      if(nopt.eq.4)inv=.true.
c
      if(nopt.ne.5)goto 75
c
      if(inv)then
        write(6,*)'Inverse operations applied'
        if(record)write(nrec,*)'Inverse operations applied'
      endif
c
      do j=nsym1,nsym2
c
        iflag=iflag+1
        write(flag1,'(i4)') iflag
        flag2(1)='_'
c
        if (index(fil,'.').eq.0) then
          i=index(fil,' ')
          file1=fil(1:i-1)//flag1//'.pdb'
        else
          i=index(fil,'.')
          file1=fil(1:i-1)//flag1//'.pdb'
        endif
c
c
c
        write(6,200)file1
        if(record)write(nrec,200)file1
200     format(' Creating PDB file:',a)
        close(unit=lun2)
c JMD!PORT:     open(unit=lun2,file=file1,status=cstatus,carriagecontrol='list')
        open(unit=lun2,file=file1,status=cstatus)
c
        if(.not.inv)call spin_p( p_xyz(1,1),coor1(1,1),ops(1,1,j),
     &    vecs(1,j),p_na,inv)
c 
c
        if(inv)call spin_p( p_xyz(1,1),coor1(1,1),ops_inv(1,1,j),
     &    vecs(1,j),p_na,inv)
c 
        do l=1,p_na
          write(lun2,1111)p_rectyp(l),p_anum(l),p_anam(l),p_rnam(l),
     &      p_rnum(l),(coor1(i,l),i=1,3),p_occ(l),
     &      p_bfac(l),iflag
c
1111      format(a6,a5,x,a5,a4,a5,4x,3f8.3,2f6.2,6x,i4)
c       
        enddo
        write(lun2,1112)
1112    format('END')
C
      enddo
C
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE peak(map, inew, icut, jstart, my_ispace)
c     ================================================
      IMPLICIT NONE
c
c Always called from peak_setup
c
c This is a simple flood-fill algorithm:
c
c  Allocate ispace*3 space for a stack of indices
c  Use the passed orthogonal coordinates to get the first index
c  Get the threshold value from the calling routine (icut)
c  Set the first point to the passed new value (inew)
c  Push its index to the stack
c
c  Pull an index off the stack
c  Get values from its 6 actual neighbours (not the diagonals)
c  If value matches test, set it to the new value and push its
c    coordinates to the stack
c  Loop until there is nothing left on the stack
c
c Should make some effort to follow xtal symmetry when desired
c
c Modified (read blatantly stolen) from Robert Esnouf's Volumes Program
c JMD Jul2000 - updeated for peak analysis, DIS April 2002
c
      INCLUDE 'average.fcm'
c
      integer    map,my_ispace,ispace
      integer*2  inew, icut
      integer    jstart(3)
c
      integer    ipixels_set,stack(2),iindex(3),jx,jy,jz,joffset,isym
      integer*8  bottom_of_stack,top_of_stack,irhosum
      integer*2  itest,gp3_a
      real       x_1(3),x_2(3)
      logical    i_box,spin_in
c
c--- Allocate stack space
c
      ispace=my_ispace*3
      bottom_of_stack=p2al(stack(1),stack(2),ispace)
      if(bottom_of_stack.eq.0)then
        write(6,*)'AV_PEAK-ERR% Allocate failed for stack - ispace =',
     &    ispace
        if(record)write(6,*)'AV_PEAK-ERR% Allocate failed for stack'
      endif
c
c--- Find the starting gridpoint
c
      iindex(1)=jstart(1)
      iindex(2)=jstart(2)
      iindex(3)=jstart(3)
c---- Convert to Angstrom
      call j3_to_f(iindex,x_1,map)
      call f_to_o(x_1,x_2,map)
c
c--- Get the test value from the starting gridpoint
c
      itest=gp3_a(map,iindex)
      write(6,*)'Peak coords (Ang, grid)',x_2,iindex
      if(record)write(nrec,*)'Peak coords (Ang, grid)',x_2,iindex
      write(6,*)'  Peak value = ',itest
      if(record)write(nrec,*)'  Peak value = ',itest
c
c--- If the test value equals the new value we're in the poo, so run away
c
      if(itest.eq.inew)then
        write(6,*)'AV_PEAK-ERR% Test value equals fill value',
     &   ' - stopping'
        if(record)write(nrec,*)'AV_PEAK-ERR% Test value equals',
     &   ' fill value - stopping'
        return
      endif
c
c--- Set it to the new value
c
      call pp3_a(map,iindex,inew)
      ipixels_set=1
      irhosum=itest
c
c--- Push it to the stack
c
      top_of_stack=bottom_of_stack+3
      stack(top_of_stack-2)=iindex(1)
      stack(top_of_stack-1)=iindex(2)
      stack(top_of_stack)=iindex(3)
c
c--- Now the main loop....
c
c--- Pop the next pixel off the stack and check its neighbours
c
100   jx=stack(top_of_stack-2)
      jy=stack(top_of_stack-1)
      jz=stack(top_of_stack)
      top_of_stack=top_of_stack-3
c
c--- Check the 6 neighbours
c
      do joffset=1,6
c
        iindex(1)=jx
        iindex(2)=jy
        iindex(3)=jz
c
        if(joffset.eq.1)then
          iindex(1)=iindex(1)-1
        elseif(joffset.eq.2)then
          iindex(1)=iindex(1)+1
        elseif(joffset.eq.3)then
          iindex(2)=iindex(2)-1
        elseif(joffset.eq.4)then
          iindex(2)=iindex(2)+1
        elseif(joffset.eq.5)then
          iindex(3)=iindex(3)-1
        else
          iindex(3)=iindex(3)+1
        endif
c
c---- If this neighbour is not in the map
c
        if(.not.i_box(iindex,map))then
c
c----- If we are using xtal symmetry, attempt to map it
c
          if(xr_con(map))then
            do isym=1,n_xsym(map)
              if(spin_in(iindex,map,isym))goto 200
            enddo
            goto 300
c
c----- Else skip to next neighbour
c
          else
            goto 300
          endif
        endif
c
c---- If this neighbour matches our target...
c
200     continue
         itest=gp3_a(map,iindex)
         if(itest.ge.icut)then
c        if(gp3_a(map,index).eq.itest)then
c
c----- Keep count
c
          ipixels_set=ipixels_set+1
          irhosum=irhosum+itest
c
c----- Set it to the new value
c
          call pp3_a(map,iindex,inew)
c
c----- Check to see if the stack is full
c
          if((top_of_stack-bottom_of_stack+3).gt.ispace)then
            write(6,*)'AV_PEAK-ERR% Stack full'
            if(record)write(nrec,*)'AV_PEAK-ERR% Stack full'
            return
          endif
c
c----- Push the coords to the stack
c
          top_of_stack=top_of_stack+3
          stack(top_of_stack-2)=iindex(1)
          stack(top_of_stack-1)=iindex(2)
          stack(top_of_stack)=iindex(3)
c
        endif
c
c---- Next neighbour
c
300     continue
      enddo
c
c--- Next pixel off stack
c
      if(top_of_stack.gt.bottom_of_stack)goto 100
c
c--- Write some stats
c
      write(6,*)'    ',ipixels_set,' pixels set to ',inew
      if(record)write(nrec,*)'    ',ipixels_set,' pixels set to ',inew
      write(6,*)'    Total sum of densities=',irhosum
      if(record)write(nrec,*)'    Total sum of densities=',irhosum
c
c--- Should be done by now, so deallocate stack space
c
      if(.not.p2dal(stack(bottom_of_stack)))then
        write(6,*)'AV_PEAK_ERR% Memory deallocate fails'
        if(record)write(nrec,*)'AV_PEAK_ERR% Memory deallocate fails'
      endif
c
c--- And we're done...
c
      return
      end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE peak_setup
c     =====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c Flood fill (usually an envelope). See subroutine flood for details.
c
c--- Local variables:
c
      integer     map
      character*8 imap,g2chru
      integer     nopt,keyword,intval,ispace
      integer*2   inew,icut,max_rho
      integer     jstart(3)
      logical     lass_get
c
c--- Get input map if not already assigned
c
      if(.not.lass_get(imap,'MAPI'))then
        imap=G2CHRU('AV_PEAK:INPUT_MAP',1,8,0)
        call lass_set(imap,'MAPI')
      endif
c
c--- Report
c
      if(verbose)then
        write(6,100)imap
        if(record)write(nrec,100)imap
100     format(/' PEAK_ANALYSIS',/,
     &          ' ============='
     &         /' Input map:       ',a)
      else
        write(6,200)imap
        if(record)write(nrec,200)imap
200     format(' PEAK_ANALYSIS  I/P map:',a)
      endif
c
c--- Get slot number for input map
c
      call av_ass(map,imap)
c
c--- The input map should exist by now !
c
      if(.not.defined(map)) then
        write(6,300)
        if(record)write(nrec,300)
300     format(' %PEAK_SETUP-ERR: Input map empty')
        return
      endif
c
c--- Keyworded input
c---
c--- FILL - value to set filled points to
c--- THRE - threshold 
c--- NPIX - way of increasing stack space for big envs
c
c--- Defaults
c
      inew=0
      ispace=100000
c
400   nopt=keyword('AV_PEAK:','NEW !THRE!NPIX!GO  !BACK',0)
c
c--- FILL VALUE
c
      if(nopt.eq.1)then
        inew=intval('AV_PEAK:ENTER_NEW_VALUE',-32000,32000,0)
      endif
c
c--- THRESHOLD
c
      if(nopt.eq.2)then
         icut=intval('AV_PEAK:ENTER_THESHOLD',
     &            -32000,32000,0)
      endif
c
c--- NPIX
c
      if(nopt.eq.3)then
        ispace=intval('AV_PEAK:ENTER_NPIX',10000,715800000,0)
      endif
c
c--- BACK
c
      if(nopt.eq.5) goto 500
c
c--- GO
c
      if(nopt.ne.4) goto 400
c
c--- Do it!
c

c--- process peaks in decreasing order
450   call max_peak(map,max_rho,jstart)
      if(max_rho.lt.icut)goto 500
      call peak(map,inew,icut,jstart,ispace)
      goto 450
c---
500   return
      end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE pip(s1,s2,flag)
c       ==========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        integer j
        logical flag
        character*(*) s1,s2
        character*8   label
c
        do j=1,n_stream
          if(s1.eq.streams(j))label=labelsa(j)
        enddo
c
        do j=1,n_stream
          if(s2.eq.streams(j))then
            labelsa(j)=label
            a_made(j) =flag
            if(flag)then
              write(6,10)label,streams(j)
              if(record)write(nrec,10)label,streams(j)
10            format(' ',a,' piped to stream ',a)
            endif
            if(.not.flag)then
              write(6,20)streams(j)
              if(record)write(nrec,20)streams(j)
20            format(' Stream: ',a,' disabled')
            endif
            return
          endif
        enddo
c
        write(6,*)'%PIP-ERR: Broken pipe!'
        if(record)write(nrec,*)'%PIP-ERR: Broken pipe!'
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE pipe_set(on)
c       =======================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        logical on
c
        pipe=on
c
        if(pipe)then
          write(6,10)
          if(record)write(nrec,10)
10        format(' PIPE option enabled')
        endif
c
c
        if(.not.pipe)then
          write(6,20)
          if(record)write(nrec,20)
20        format(' PIPE option disabled')
        endif
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE pmud(lun,line,num)
c       ================================
      IMPLICIT NONE
        integer lun,num
        byte line(num)
        read(lun)line
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE point_exp
c     ====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c---- routine to expand ncs operators using appropriate xtallographic
c---- operators (ie not unitary and those with translational 
c---- components to them
c
c--- local variables
C
      LOGICAL point_ok
c
      character*8 G2CHRU,ch_temp
      integer     mapn
      character*8 mapl
      integer     n_ok
      integer     ioff,i_sym,l,ibase,i,j_sym,j_nsym,nsym_now,j,k
c
      integer*2   mat(3,3)
      real        rmat(3,3),back(3,3),vback(3),rinv(3,3),det
c
c
c      equivalence  (i_temp,ch_temp)
c
      point_ok=.false.
      n_ok=0
c
      write(6,*)'EXPAND NCS OPS USING CRYSTALLOGRAPHIC POINT OPS'
      if(record)write(nrec,*)
     &          'EXPAND NCS OPS USING CRYSTALLOGRAPHIC POINT OPS'
c
c----update cell information
      ch_temp=G2CHRU('AV_PEXP:ENTER_MAP_NAME',1,8,0)
      mapl=ch_temp
      call av_ass(mapn,mapl)
c
      if(.not.terse)then
        if(.not.defined(mapn)) write(6,180)
        if(record)then
          if(.not.defined(mapn)) write(nrec,180)
        endif
      endif
180   format(' The map does not exist, abort plane expansion')
c
      if(.not.defined(mapn)) goto 900
      if(nsym.lt.1) then
        write(6,*) 'No NCS ops, abort plane expansion'
        if(record) write(nrec,*) 'No NCS ops, abort plane expansion'
        goto 900
      endif
c
c---- now loop over xtal no of sym ops for this map
c
      do j_sym = 1,n_xsym(mapn)
c
c--- now set up pointers to pick up operator from array...stored as
c--- single array sized as 12*no of xtal symm operator
        ioff=(j_sym-1)*12
        do i_sym=1,3
          l=i_sym*4
          ibase=l-4
c--- loop thro translations...if not eq to 0....drop out...
c--- this means translation component which we dont want
          if(xr_sym(l+ioff,mapn).ne.0) goto 800
c--- now loop thro xtal rot matrix
          do i=1,3
            mat(i,i_sym)= xr_sym(ioff+ibase+i,mapn)
            rmat(i,i_sym)=float(mat(i,i_sym))
          enddo
        enddo
c---- check if unitary.....
        if( (mat(1,1).eq.1).and.(mat(2,2).eq.1)
     &              .and.(mat(3,3).eq.1) ) goto 800
        point_ok=.true.
        n_ok=n_ok+1
        goto 400
c---- if matrix ok loop over ncs ops and multiply with current xtal matrix
c---- ie not unitary.....if trans component should have dropped out....
400     if(point_ok) then
c
c--- loop of ncs symms
          do j_nsym = 1 ,nsym
c---- here we are multiplying [xtal].[NCSO] and [XTAL].vec.....
            call matmul(rmat,ops(1,1,j_nsym),back(1,1))
            call vmtply(rmat,vecs(1,j_nsym),vback(1))
            nsym_now= nsym*n_ok + j_nsym
            if(.not.terse)then
              write(6,*) ' Creating ncs ops ', nsym_now
              if(record)write(nrec,*) ' Creating ncs ops ', nsym_now
            endif
c---and creating new ncs ops
            call invers(back,rinv,det)
            call copy(9,back,ops(1,1,nsym_now))
            call copy(9,rinv,ops_inv(1,1,nsym_now))
            call copy(3,vback,vecs(1,nsym_now))
            if(.not.terse)then
              write(6,100)nsym_now,((back(j,k),k=1,3),j=1,3),
     &          (vback(i),i=1,3)
              if(record)write(nrec,100)nsym_now,((back(j,k),k=1,3),
     &          j=1,3),(vback(i),i=1,3)
100           format(10x,i6,' Symm ops now in use',/(3(10x,3f11.5,/),
     &             10x,3f11.2))
            endif
          enddo
        endif
        point_ok=.false.
c
c JMD!PORT: Changed to remove compiler warning - check if correct!
c800   enddo          
800     continue
      enddo
c
      write(6,200) nsym_now,n_ok
      if(record)write(nrec,200)
200     format(5x,i6,' Symm ops now in use, generated from',i4,'
     &    xtal ops')
c
      if(n_ok.gt.0)nsym=nsym_now
c
900   return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE POLAR(R,PHI,PSI,CHI)
c     ===============================
      IMPLICIT NONE
c
c*******************************************************************************
        real pi,piovr2,r,phi,psi,chi,zilch
c
      DIMENSION R(3,3)
      PI = 3.1415926535898
      PIOVR2 = PI/2
      ZILCH = 0.0000001
      IF(ABS(R(1,2)-R(2,1)).GT.ZILCH.OR
     &  .ABS(R(2,3)-R(3,2)).GT.ZILCH) GO TO 40
      IF(ABS(R(1,3)-R(3,1)).LT.ZILCH) GO TO 5
C 
      PHI = 0
      PSI = 0
      CHI = ATAN2(R(3,1),R(1,1))
      RETURN
C 
    5 IF(R(1,2).GT.ZILCH.OR.R(2,3).GT.ZILCH) GO TO 10
      PHI = 0
      PSI = 0
      CHI = 0
      IF(R(1,1).GT.0.AND.R(2,2).GT.0.AND.R(3,3).GT.0) RETURN
      CHI = PI
      IF(R(2,2).GT.0) RETURN
      PSI = PIOVR2
      IF(R(3,3).GT.0) RETURN
      PHI = PIOVR2
      RETURN
C 
   10 PHI = ATAN2(-R(2,3),R(1,2))
      CHI = PI
      IF(R(1,2).GT.ZILCH.AND.R(2,3).GT.ZILCH) GO TO 20
      PSI = ATAN2((R(1,2)-R(2,3)),(R(2,2)+1.0))
      RETURN
   20 IF(ABS(PHI).LT.0.10.OR.ABS(PHI).GT.3.00) GO TO 30
      PSI = ATAN2(R(3,1)/SIN(PHI),-R(1,2))
      RETURN
   30 PSI = ATAN2(R(3,1)/COS(PHI),R(2,3))
      RETURN
C 
   40 PHI = ATAN2((R(2,1)-R(1,2)),(R(2,3)-R(3,2)))
      IF(ABS(PHI).LT.0.10.OR.ABS(PHI).GT.3.00) GO TO 50
      PSI = ATAN2(((R(2,1)-R(1,2))/SIN(PHI)),(R(3,1)-R(1,3)))
      GO TO 60
   50 PSI = ATAN2(((R(2,3)-R(3,2))/COS(PHI)),(R(3,1)-R(1,3)))
   60 IF(ABS(PSI).GT.0.10.AND.ABS(PSI).LT.3.00) GO TO 80
      IF(ABS(PHI).LT.0.10.OR.ABS(PHI).GT.3.00) GO TO 70
      CHI = ATAN2(((R(2,1)-R(1,2))/(SIN(PSI)*SIN(PHI))),
     &            (R(1,1)+R(2,2)+R(3,3)-1.0))
      RETURN
   70 CHI = ATAN2(((R(2,3)-R(3,2))/(SIN(PSI)*COS(PHI))),
     &            (R(1,1)+R(2,2)+R(3,3)-1.0))
      RETURN
   80 CHI = ATAN2(((R(3,1)-R(1,3))/COS(PSI)),(R(1,1)+R(2,2)+R(3,3)-1.))
      RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE pp3(map,j1,j2,j3,irho)
c       =================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        external brickel_3
        integer  brickel_3,map
        integer j_local(3),j1,j2,j3,nbel
        integer*2 irho
c
        j_local(1)=j1
        j_local(2)=j2
        j_local(3)=j3
        nbel=brickel_3(j_local,map)
        scratch(ns(map)+nbel)=irho
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE pp3_a(map,j1,irho)
c       =============================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        external brickel_3
        integer  brickel_3,map
        integer  j1(3),nbel
        integer*2 irho
c
        nbel=brickel_3(j1,map)
        scratch(ns(map)+nbel)=irho
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE pri_cut
c       ==================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        write(6,10) ncut_u,ncut_l
        if(record)write(nrec,10) ncut_u,ncut_l
10      format(' CUT_OFF VALUES',/,
     &     ' Upper cut off for density:',i5,'   lower cut off:',i5)
c
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE prop(flag)
c       =====================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c       =====================
c
        logical flag
c
c
        proper=flag
c
        if (flag) then
          write(6,30)
          if(record)write(nrec,30)
30        format(' NON-CRYSTALLOGRAPHIC SYMMETRY IS PROPER')
         else
          write(6,40)
          if(record)write(nrec,40)
40        format(' NON-CRYSTALLOGRAPHIC SYMMETRY IS IMPROPER')
          endif
c
        return
        end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE prtsym(nun,n)
c       ==========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
        integer nun,n,k,j
c
        write (nun,100) n,((ops(k,j,n),j=1,3),k=1,3), (vecs(j,n),j=1,3)
        if(record)write (nrec,100)
     &                n,((ops(k,j,n),j=1,3),k=1,3), (vecs(j,n),j=1,3)
100     format (15x,' Symm op n.o. ',I4,3(/,10x,3f10.5)/,10x,3f10.5)
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      SUBROUTINE QKFIT (UMAT, RTSUM, R, ENTRY_EV )
c     =========================================
      IMPLICIT NONE
c 
      DOUBLE PRECISION UMAT(3,3),RTSUM
      DIMENSION  R(3,3)
      LOGICAL ENTRY_EV
c
C----                                                  
C----   THE 'EIGENVALUE ONLY' ENTRY WAS                
C----   ADAPTED FROM PROGRAM BY A.D. MCLACHAN 7/78     
C----                       
      DOUBLE PRECISION USQMAT(3,3),ROOT(3),A(3,3),B(3,3)
      DOUBLE PRECISION UTR(6),AAM,BAM,FAM,GAM,HAM
      DOUBLE PRECISION CAM,cc,dd,qq,cos3th,argsq,theta,diff
      DOUBLE PRECISION digav,EPS,third,forthr,half,one
      DOUBLE PRECISION two,three
      integer i,ia,isig,j,k
      double precision pi,du11,du21,du31,detu,b1,b2,b13,b23,b33
      double precision t,rt,s
      real r
C----                                                   
      EQUIVALENCE (AAM,USQMAT(1,1)),(BAM,USQMAT(2,2)),
     & (CAM,USQMAT(3,3))
      EQUIVALENCE (FAM,USQMAT(2,3)),(GAM,USQMAT(1,3))
     &,(HAM,USQMAT(1,2))
      EQUIVALENCE( A(1,1), USQMAT(1,1) ), (UTR(1), B(1,1))
C----                                                   
      DATA EPS/1.0D-8/
      DATA PI/3.14159265358979/
      DATA ONE/1.0/,TWO/2.0/,THREE/3.0/
      DATA HALF/0.5/,THIRD/0.33333333333/,FORTHR/1.33333333333/
      DATA USQMAT(2,1),USQMAT(3,1),USQMAT(3,2)/3*0.0/
      ISIG=1
C----                                                       
C----    IF ENTRY IS .TRUE. GET OUT THE ROTATION MATRIX     
C----                                                       
      IF (ENTRY_EV) GO TO 200
C----                                                       
C----   CALC DET OF UMAT                                    
C----                                                       
      DU11=(UMAT(2,2)*UMAT(3,3))-(UMAT(2,3)*UMAT(3,2))
      DU21=(UMAT(2,3)*UMAT(3,1))-(UMAT(2,1)*UMAT(3,3))
      DU31=(UMAT(2,1)*UMAT(3,2))-(UMAT(2,2)*UMAT(3,1))
      DETU=(UMAT(1,1)*DU11)+(UMAT(1,2)*DU21)+(UMAT(1,3)*DU31)
c      WRITE(6,999) 'detu =', DETU                                      
c999   FORMAT(/(3(a,F12.5)))                                      
      IF(DETU.LT. 0.0) ISIG=-1
C----                                                        
C----   FORM USQMAT AS POSITIVE SEMI DEFINITE MATRIX         
C----                                                        
      DO 110 J=1,3
      DO 105 I=1,J
      USQMAT(I,J)=(UMAT(1,I)*UMAT(1,J))+(UMAT(2,I)*UMAT(2,J))+
     &   (UMAT(3,I)*UMAT(3,J))
105   CONTINUE
110   CONTINUE
C----                                                         
C----   REDUCE AVG OF DIAGONAL TERMS TO ZERO                  
C----                                                         
      DIGAV=(AAM+BAM+CAM)*THIRD
      AAM=AAM-DIGAV
      BAM=BAM-DIGAV
      CAM=CAM-DIGAV
C----                                                         
C----   SETUP COEFFS OF SECULAR EQUATION OF MATRIX WITH TRACE ZERO
C----                                                             
      CC=(FAM*FAM)+(GAM*GAM)+(HAM*HAM)-(AAM*BAM)-(BAM*CAM)-(CAM*AAM)
      DD=(AAM*BAM*CAM)+TWO*(FAM*GAM*HAM)-(AAM*FAM*FAM)
     &   -(BAM*GAM*GAM)-(CAM*HAM*HAM)
C----                                                               
C----   THE SECULAR EQN IS Y**3-CC*Y-DD=0  AND DD IS DET(USQMAT)    
C----   REDUCE THIS TO THE FORM COS**3-(3/4)COS-                    
C----   (1/4)COS3THETA = 0                                          
C----   WITH SOLUTIONS COSTHETA.  SO Y=QQ*COSTHETA                  
C----                                                               
      IF(CC.LE.EPS) GO TO 115
      QQ=SQRT(FORTHR*CC)
      COS3TH=(THREE*DD)/(CC*QQ)
      IF(ABS(COS3TH).GT.ONE) COS3TH=SIGN(ONE,COS3TH)
C----                                                               
C----   FUNCTION ARCOS                                              
C----                                                               
      IF(COS3TH.NE.0.0) GOTO 1200
       THETA= 1.570796327
      GOTO 1400
1200  ARGSQ=COS3TH*COS3TH
      THETA=ATAN(SQRT(1.0-ARGSQ)/COS3TH)
      IF(COS3TH.LT.0.0) THETA=PI-THETA
1400  CONTINUE
C----                                                               
C----   ROOTS IN ORDER OF SIZE GO 1,2,3 1 LARGEST                   
C----                                                               
      THETA=THETA*THIRD
      ROOT(1)=QQ*COS(THETA)
      DIFF=HALF*SQRT(THREE*(QQ*QQ-ROOT(1)*ROOT(1)))
      ROOT(2)=-ROOT(1)*HALF+DIFF
      ROOT(3)=-ROOT(1)*HALF-DIFF
      GO TO 120
115   CONTINUE
C----                                                               
C----   SPECIAL FOR TRIPLY DEGENERATE                               
C----                                                               
      ROOT(1)=0.0
      ROOT(2)=0.0
      ROOT(3)=0.0
120   CONTINUE
C----   ADD ON DIGAV AND TAKE SQRT                                  
      DO 125 J=1,3
      RT=ROOT(J)+DIGAV
      IF(RT.LT.EPS) RT=0.0
      ROOT(J)=SQRT(RT)
125   CONTINUE
c%      WRITE(6,999) 'root', ROOT                                             
C----   IF DETU IS <0 CHANGE SIGN OF ROOT(3)                        
      IF(ISIG.EQ.-1) ROOT(3)=-ROOT(3)
      RTSUM=ROOT(1)+ROOT(2)+ROOT(3)
c%      WRITE(6,999) 'rtsum=', RTSUM                                            
      RETURN
C----                                                               
C----   THIS IS THE FANCY PART                                      
C----                                                               
200   CONTINUE
C----                                                               
C----   FORM USQ = (UT).U    (IN UPPER TRIANGULAR SYMMETRIC STORAGE MODE)
C----                                                          
      DO 220 I=1,3
      DO 220 J=I,3
      T = 0.0
      DO 210 K=1,3
      T = T + UMAT(K,I)*UMAT(K,J)
  210 CONTINUE
      IA = I + (J*J-J)/2
      UTR(IA) = T
  220 CONTINUE
c%      WRITE(6,999) 'utr=', UTR                                         
C----                                                          
C----   CALCULATE EIGENVALUES AND VECTORS                      
C----                                                          
      CALL SJ33EV (UTR, A, 3, 0)
c%      WRITE(6,999) 'utr=', UTR                                         
C----                                                          
      ROOT(1) = UTR(1)
      ROOT(2) = UTR(3)
      ROOT(3) = UTR(6)
c%      WRITE(6,999) 'root=', ROOT                                        
C----                                                          
C----   SET A3 = A1 CROSS A2                                   
C----   ROOTS ARE IN ORDER R(1) >= R(2) >= R(3) >= 0           
C----                                                          
      A(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      A(2,3) = A(3,1)*A(1,2) - A(1,1)*A(3,2)
      A(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)
c%      WRITE(6,*) 'a(1,3)=', a(1,3)
c%      WRITE(6,*) 'a(2,3)=', a(2,3)
c%      WRITE(6,*) 'a(3,3)=', a(3,3)
C----                                                          
C----   VECTOR SET B=U.A                                       
C----                                                          
      DO 240 I=1,3
      DO 240 J=1,3
      T = 0.0
      DO 230 K=1,3
  230 T = T + UMAT(J,K)*A(K,I)
      B(J,I) = T
  240 CONTINUE
C----                                                          
C----    NORMALIZE B1 AND B2 AND CALCULATE B3 = B1 CROSS B2    
C----                                                          
      B1 = SQRT( B(1,1)*B(1,1) + B(2,1)*B(2,1) + B(3,1)*B(3,1) )
      B2 = SQRT( B(1,2)*B(1,2) + B(2,2)*B(2,2) + B(3,2)*B(3,2) )
c%      WRITE(6,*) 'B1 IS', B1
      DO 250 I=1,3
      B(I,1) = B(I,1)/B1
  250 B(I,2) = B(I,2)/B2
C----                                                           
C----    CHECK FOR LEFT HANDED ROTATION                         
C----                                                           
      B13 = B(2,1)*B(3,2) - B(3,1)*B(2,2)
      B23 = B(3,1)*B(1,2) - B(1,1)*B(3,2)
      B33 = B(1,1)*B(2,2) - B(2,1)*B(1,2)
C----                                                           
      S = B13*B(1,3) + B23*B(2,3) + B33*B(3,3)
      IF (S .LT. 0.0) ISIG = -1
      B(1,3) = B13
      B(2,3) = B23
      B(3,3) = B33
c---      WRITE(6,999) B                                            
C----                                                           
C----   CALCULATE ROTATION MATRIX R                             
C----                                                           
      DO 270 I=1,3
      DO 270 J=1,3
      T = 0.0
      DO 260 K=1,3
  260 T = T + B(I,K)*A(J,K)
      R(I,J) = SNGL(T)
  270 CONTINUE
C----                                                           
C----   RMS ERROR                                               
C----                                                           
      DO 280 I=1,3
      IF (ROOT(I) .LT. 0.0) ROOT(I) = 0.0
      ROOT(I) = SQRT( ROOT(I) )
  280 CONTINUE
C----                                                           
C----   CHANGE SIGN OF EVAL #3 IF LEFT HANDED                   
C----                                                           
      IF (ISIG .LT. 0) ROOT(3)=-ROOT(3)
      RTSUM = ROOT(3) + ROOT(2) + ROOT(1)
      RETURN
      END
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE qmem
c       ===============
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        write(6,10)
        if(record)write(nrec,10)
10      format(' Version 5.0 onwards uses dynamic memory allocation')
        return
        end

C=============================================================================

      SUBROUTINE QSORTP(PARRAY, TOT, STACK, MAXSTK, COMPAR, ERRCOD)
      IMPLICIT NONE

C Quick sorting routine for ordering the different rotation matrices


        INTEGER TOT, MAXSTK, COMPAR, ERRCOD
        INTEGER PARRAY (TOT), STACK (MAXSTK)
        EXTERNAL COMPAR

        INTEGER S, L, R, I, J, X, PSWAP

C Just skip if nothing to sort

        ERRCOD = 0
        IF (TOT .LE. 1) RETURN

C Initialize stack

        S = 2
        IF (S .GT. MAXSTK) GOTO 900
        STACK (S - 1) = 1
        STACK (S)     = TOT

C Simulated repeat-until loop 1; pop stack until empty

100     L = STACK (S - 1)
        R = STACK (S)
        S = S - 2

C Simulated repeat-until loop 2; split interval until nothing left

200       I = L
          J = R

C Choose (arbitrarily) value in middle of interval as comparand

          X = (L + R) / 2

C Simulated repeat-until loop 3; partition until nothing left

300         IF (COMPAR (PARRAY (I), PARRAY (X)) .LT. 0) THEN
              I = I + 1
              GOTO 300
            ENDIF
310         IF (COMPAR (PARRAY (X), PARRAY (J)) .LT. 0) THEN
              J = J - 1
              GOTO 310
            ENDIF

C Swap array values if on correct side of each other in interval;
C beware of what happens to the pointer to the comparand slot

            IF (I .LE. J) THEN
              PSWAP      = PARRAY (I)
              PARRAY (I) = PARRAY (J)
              PARRAY (J) = PSWAP
              IF (X .EQ. I) THEN
                X = J
              ELSE IF (X .EQ. J) THEN
                X = I
              END IF
              I = I + 1
              J = J - 1
            ENDIF

C End of simulated repeat-until loop 3

          IF (I .LE. J) GOTO 300

C Recursion; process shorter interval first, and push interval to stack

          IF (J - L .LT. R - I) THEN

            IF (I .LT. R) THEN
              S = S + 2
              IF (S .GT. MAXSTK) GOTO 900
              STACK (S - 1) = I
              STACK (S)     = R
            ENDIF
            R = J

          ELSE

            IF (L .LT. J) THEN
              S = S + 2
              IF (S .GT. MAXSTK) GOTO 900
              STACK (S - 1) = L
              STACK (S)     = J
            ENDIF
            L = I

          ENDIF

C End of simulated repeat-until loop 2

          IF (L .LT. R) GOTO 200

C End of simulated repeat-until loop 1

        IF (S .GT. 0) GOTO 100

        RETURN

900     ERRCOD = 1

        RETURN
        END
c
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE qsym
c       ===============
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        integer n
c
        if(nsym.gt.0)then
          write(6,10)
          if(record)write(nrec,10)
10        format(' Currently active NCS operators:')
          do n=1,nsym
            call prtsym(6,n)
          enddo
c
        else
          write(6,*)'No NCS operators defined at present'
          if(record)write(nrec,*)'No NCS operators defined at present'
        endif
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE r_stat
c      =================
      IMPLICIT NONE
c
c--- These subroutines accumulate statistics on the map during a RESEt.
c--- Max, min, mean and sigma are calculated.
c
      INCLUDE 'average.fcm'
c      =====================
c
c--- Local variables
c
       integer   num_pixels
       integer*2 irho,irho_max,irho_min
       real*4    firho,sum_r,sum_sq,sav_tempf,ssign,ssig
       real*4    rr_max,rr_min
c
c--- Need to save running totals and min/maxes
c
       save  num_pixels,sum_r,sum_sq,irho_max,irho_min
c
c--- Initialise counters
c
       av_temp=zero
       sigm=zero
       sum_r=zero
       sum_sq=zero
       num_pixels=0
       irho_max=-32000
       irho_min=32000
c
       return
c
c--- Accumulate statistics
c
       entry r_a_stat(irho)
c      ====================
c
       firho=float(irho)
c
c--- Running totals
c
       num_pixels=num_pixels+1
       sum_r=sum_r+firho
       sum_sq=sum_sq+(firho**2)
c
c      write(6,*)'num_pixels=',num_pixels
c
c
c--- Rho min/max:
c
       if(irho.gt.irho_max)irho_max=irho
       if(irho.lt.irho_min)irho_min=irho
c
       return
c
c--- Calculate final stats
c
       entry r_p_stat(rr_max,rr_min)
c      =============================
c
c--- Sanity check
c
       if(num_pixels.le.0) then
         write(6,10)
         if(record)write(nrec,10)
         return
       endif
c
10     format(' %R_STAT-ERR: No pixels to analyse')
c
c--- Calculate statistics - av_temp, sigm, rr_max, rr_min
c--- NB av_temp and sigm are in one of the average.fcm common blocks
c
       av_temp=sum_r/float(num_pixels)
c
       sav_tempf=float(num_pixels)*(av_temp**2 )
       ssign=sav_tempf + sum_sq - (2.0 * av_temp * sum_r)
       ssig=ssign/float(num_pixels)
       if(ssig.gt.0.0)sigm=sqrt(ssig)
c
       rr_max=float(irho_max)
       rr_min=float(irho_min)
c
       return
       end
c
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE read_headers(lun,title,iuvw,mxyz,nx,norg,nstart,nend
     &      ,cell,lspgrp,rholim,scale_r,plus,lmode,terse,verbose,speed)
C     ===============================================================
      IMPLICIT NONE
c
C  Formatted, read keyworded header, fixed format
c
c
c*****************************************************************************
        common/av_rec /record,nrec
c*****************************************************************************
        logical record,terse,verbose
        integer nrec
c
c
c
      integer      lun
      character*80 title
      character*4  speed
      integer      iuvw(3),mxyz(3),nx(3),norg(3),nstart(3),nend(3)
     &            ,lspgrp,lmode
      real         cell(6),scale_r,plus,rholim(4)
      integer      i,j,limits(2,3)
      character*1  xyz(3),axes(3)
c
c
      data xyz/'X','Y','Z'/
c
c
      if(verbose)then
        write(6,5000)
        if(record)write(nrec,5000)
 5000   format(' HEADER INFORMATION AS READ FROM INPUT MAP',/
     &          ' =========================================')
        else
        if(.not.terse)then
          write(6,5010)
          if(record)write(nrec,5010)
 5010     format(' HEADER INFORMATION AS READ FROM INPUT MAP')
          endif
        endif
c
        if(speed.eq.'slow')then
c
      read(lun,6000)
      read(lun,6010) title
      read(lun,6020) axes
      read(lun,6030) mxyz
      read(lun,6040) limits
      read(lun,6050) lspgrp
      read(lun,6060) lmode
      read(lun,6070) cell
      read(lun,6080) rholim
      read(lun,6090) scale_r,plus
      read(lun,6095)
c
        else
c
      read(lun) title,axes,
     & mxyz, limits, lspgrp, lmode, cell, rholim, scale_r, plus
c
        endif
c
 6000 format('MAPEXCHANGE HEADER')
      if(.not.terse)then
        write(6,6100)
        if(record)write(nrec,6100)
        endif
6100  format(' MAPEXCHANGE HEADER')
 6010 format(/A)
      if(.not.terse)then
        write(6,6110) title
        if(record)write(nrec,6110) title
 6110   format(' TITLE'/1X,A)
        endif
 6020 format('AXIS    ',3(7X,A1))
      if(.not.terse)then
        write(6,6120) axes
        if(record)write(nrec,6120) axes
 6120   format(' AXIS    ',3(7X,A1))
        endif
 6030 format('GRID    ',3I8)
      if(.not.terse)then
        write(6,6130) mxyz
        if(record)write(nrec,6130) mxyz
 6130   format(' GRID    ',3I8)
        endif
 6040 format('XYZLIM  ',6I8)
      if(.not.terse)then
        write(6,6140) limits
        if(record)write(nrec,6140) limits
 6140   format(' XYZLIM  ',6I8)
        endif
 6050 format('SPACEGROUP',6X,I8)
      if(.not.terse)then
        write(6,6150) lspgrp
        if(record)write(nrec,6150) lspgrp
 6150   format(' SPACEGROUP',6X,I8)
      endif
 6060 format('MODE    ',I8)
      if(.not.terse)then
        write(6,6160) lmode
        if(record)write(nrec,6160) lmode
 6160   format(' MODE    ',I8)
        endif
 6070 format('CELL    ',6F10.3)
      if(.not.terse)then
        write(6,6170) cell
        if(record)write(nrec,6170) cell
 6170   format(' CELL    ',6F10.3)
        endif
 6080 format('RHOLIM  ',4G16.6)
      if(.not.terse)then
        write(6,6180) rholim
        if(record)write(nrec,6180) rholim
 6180   format(' RHOLIM  ',4G16.6)
        endif
 6090 format('PACK    ',2G16.6)
      if(.not.terse)then
        write(6,6190) scale_r,plus
        if(record)write(nrec,6190) scale_r,plus
 6190   format(' PACK    ',2G16.6)
        endif
 6095 format('END HEADER')
      if(.not.terse)then
        write(6,6195)
        if(record)write(nrec,6195)
 6195   format(' END HEADER')
        endif
c
C Get axis order
      do 10, i=1,3
         do 20, j=1,3
            if (axes(i) .eq. xyz(j)) then
               iuvw(i) = j
            endif
 20      continue
 10   continue
c
c       get limits
        do i=1,3
          nstart(i)=limits(1,i)
          nend(i)=limits(2,i)
          nx(i)=nend(i)-nstart(i)+1
          norg(i)=1-nstart(i)
        enddo
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE read_map1(speed)
c     ===========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c
      character*4     speed
      character*8     label,label2
      integer         nmap,ics_local
      character*80    fil,file2
      integer         lun,j,nline
      integer*8       nxtbyt
      byte            abyte(2)
c
      character*80    charval
      character*8     G2CHRU
      character*8     ch_temp
      real            xrvol
      logical         lass_get
c
c      equivalence     (i_temp,ch_temp)
c
c
      if(.not.lass_get(ch_temp,'MAPI'))then
c
        ch_temp=G2CHRU('AV_READ_MAP:ENTER_LABEL',1,8,0)
        write (6,20) ch_temp
        if(record)write (nrec,20) ch_temp
20      format(' label of map to be read = ', a)
        call lass_set(ch_temp,'MAPI')
      endif
      label=ch_temp

      write(6,10)ch_temp
      if(record)write(nrec,10)ch_temp
10    format(' READ MAP  label:',a)

      fil=charval('AV_READ_MAP:ENTER_FILE',1,80,0)
      lun=1
      close (unit=lun)
c JMD!PORT: if(speed.eq.'fast')open (unit=lun,file=file,status='old',readonly
      if(speed.eq.'fast')open (unit=lun,file=fil,status='old'
     &                          ,form='unformatted')
      if(.not.terse)then
        write(6,40)fil(1:60), lun
        if(record)write(nrec,40)fil(1:60), lun
40      format(' Mapfile:', a60,/,' opened on unit=',i3)
      else
        write(6,41)fil(1:60)
        if(record)write(nrec,41)fil(1:60)
41      format(' Mapfile:', a60)
      endif
c
      goto 777
c
      entry read_map(speed,label2,file2)
c       ===================================
c
      label=label2
      fil=file2
c JMD!PORT: open (unit=lun,file=file,status='old',readonly
      open (unit=lun,file=fil,status='old'
     &                          ,form='formatted')
c
777   continue
      call av_ass(nmap, label)
      if(defined(nmap))then
        write(6,*)
     &    '%READ_MAP-ERR: Map defined already, abort read'
        if(record)write(nrec,*)
     &    '%READ_MAP-ERR: Map defined already, abort read'
        goto 999
      endif
c
      call read_headers(lun,titles(nmap),iuvw(1,nmap),nunit(1,nmap),
     &  nx(1,nmap),norg(1,nmap),nstart(1,nmap),nend(1,nmap)
     &  ,XRcell(1,nmap),lgrp(nmap),
     &  rholim(1,nmap),mscale(nmap),moffset(nmap),mtype(nmap),
     &  terse,verbose,speed)
c
      if(.not.terse)then
        write(6,*) 'Scale, offset:',mscale(nmap),moffset(nmap)
        if(record)write(nrec,*)
     &    'Scale, offset:',mscale(nmap),moffset(nmap)
      endif
c
      call XRfrac(XRcell(1,nmap),XRtr(1,1,nmap),XRintr(1,1,nmap),
     &             XRvol,.true.)
c
c---- set up symmetry data
c
      call XR_symsy(XRcell(1,nmap),lgrp(nmap),ics_local,
     &                  n_xsym(nmap),XR_sym(1,nmap) )
      if(verbose)then
        write(6,*)'Number of symm ops for this space group=',
     &                  n_xsym(nmap)
        if(record)write(nrec,*)
     &    'Number of symm ops for this space group=',n_xsym(nmap)
      endif
c
c---- and integerized CSO's:
c
      call XR_isym(nmap)
c
c---- also set up pointers indicating which symm op to try first
c
      do j=1,maxsym
        last_sym(j,nmap)=1
      enddo
c
      npix(nmap)=nx(1,nmap)*nx(2,nmap)*nx(3,nmap)
c
c---- set up bricking data
c
      call brickit(nmap)
c
c---- allocate space in scratch array
c
      call alloc(nmap)
c
c How much scratch space do we need for reading?
c
      nline=nx(iuvw(1,nmap),nmap)
c
c Allocate some integer*2 scratch space
c
      nxtpix=P2AL(scratch(1),scratch(2),nline)
      if(nxtpix.eq.0) then
        write(6,*)'%READ_MAP-ERR: Memory allocate fails'
        if(record)write(nrec,*)'%READ_MAP-ERR: Memory allocate fails'
        goto 999
      endif
c
c Allocate some byte scratch space
c
      nxtbyt=P2AL(abyte(1),abyte(2),nline)
      if(nxtbyt.eq.0) then
        write(6,*)'%READ_MAP-ERR: Memory allocate fails'
        if(record)write(nrec,*)'%READ_MAP-ERR: Memory allocate fails'
        goto 999
      endif
c
c Actually read the map
c
      call gulp_map(lun,scratch(nxtpix),abyte(nxtbyt),nline,nmap,speed)
c
c Finished reading so close the unit
c
      close(unit=lun)
c
c Deallocate scratch space
c
      if(.not.P2DAL(scratch(nxtpix)))then
        write(6,*)'%READ_MAP-ERR: Memory deallocate fails'
        if(record)write(nrec,*)'%READ_MAP-ERR: Memory deallocate fails'
        goto 999
      endif
      if(.not.P2DAL(abyte(nxtbyt)))then
        write(6,*)'%READ_MAP-ERR: Memory deallocate fails'
        if(record)write(nrec,*)'%READ_MAP-ERR: Memory deallocate fails'
        goto 999
      endif
c
c Done!
c
      return
c
c On error...
c
999   write(6,9999)
      if(record)write(nrec,9999)
9999  format(' %READ_MAP-ERR: **ERROR: Abandon read** ')
      return
      end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE read_map_check(fixscale)
c       ===================================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c
        character*8  label
        character*80 fil,line
        integer      lun
c
        character*80 charval
        character*8  G2CHRU,ch_temp
        real         realval,usescale
        logical      lass_get,fixscale
c
c
c
        usescale=0.0
c
        if(.not.lass_get(ch_temp,'MAPI'))then
c
          ch_temp=G2CHRU('AV_READ_MAP:ENTER_LABEL',1,8,0)
          write (6,20) ch_temp
          if(record)write (nrec,20) ch_temp
20        format(' label of map to be read = ', a)
          call lass_set(ch_temp,'MAPI')
        endif
        label=ch_temp

        write(6,10)ch_temp
        if(record)write(nrec,10)ch_temp
10      format(' READ MAP  label:',a)

        fil=charval('AV_READ_MAP:ENTER_FILE',1,80,0)
        write(6,*)fil
        lun=1
        close (unit=lun)
c JMD!PORT: gfortran f77 open has no readonly attribute
c open (unit=lun,file=file,status='old',form='formatted',readonly)     
        open (unit=lun,file=fil,status='old',form='formatted')
c
c - Jiffy to apply defined scale factor on read in (mapb only!)
c
        if(fixscale)then
          usescale=realval('AV_READ_MAP:ENTER_SCALE',0.0,1000000.0,0)
        endif
c
          read(lun,6001,err=100) line
6001      format(a)
          if(index(line,'MAPEXCHANGE') .gt.0)then
            close(unit=lun)
            if(.not.terse)then
              write(6,*)'Formatted file'
              if(record)write(nrec,*)'Formatted file'
            endif
            call read_map('slow',label,fil)
          else
            close(unit=lun)
            if(.not.terse)then
              write(6,*)'Unformatted file'
              if(record)write(nrec,*)'Unformatted file'
            endif
            call read_mapb(label,fil,fixscale,usescale)
          endif
          return
c
100       close(unit=lun)
          if(.not.terse)then
            write(6,*)'Unformatted file'
            if(record)write(nrec,*)'Unformatted file'
          endif
          call read_mapb(label,fil,fixscale,usescale)
          return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE read_mapb(label,fil,fixscale,scale_r)
c     ================================
      IMPLICIT NONE
c
c Purpose: Read CCP4 binary maps
c
c Declare common parameters and variables
c 
      INCLUDE 'average.fcm'
c     =====================
c
c Declare passed variables
c
      character*8     label
      character*80    fil
      logical         fixscale
      real            scale_r
c
c Declare local variables
c
      integer         nmap,ics_local,lun,i,j,limits(2,3)
      integer         nline,nsize,ierr,iprint
      integer         nsec,nw1,nu1,nu2,nv1,nv2
      integer*8       npoint,np1
      real            xrvol,stub(2),rhmin,rhmax
c
c Initialise some variables
c
      ierr=1
      iprint=1
      if(terse)iprint=0
      lun=1
c
c Check we are storing the map in an empty slot
c
      call av_ass(nmap,label)
      if(defined(nmap))then
        write(6,*)'%READ_MAP-ERR: Map defined already, abort read'
        if(record)write(nrec,*)
     &    '%READ_MAP-ERR: Map defined already, abort read'
        goto 999
      endif
c
c  Open file on unit lun and read map headers
c
      call mrdhds(lun,fil(1:80),titles(nmap),nsec,iuvw(1,nmap),
     &  nunit(1,nmap),nw1,nu1,nu2,nv1,nv2,XRcell(1,nmap),lgrp(nmap),
     &  mtype(nmap),rholim(1,nmap),rholim(2,nmap),rholim(3,nmap),
     &  rholim(4,nmap),ierr,iprint)
      if(ierr.eq.-1)then
        print*,'%READ_MAP-ERR: mrdhds returned an error in read_mapb'
        return
      endif
c
c Work out scale factor to use dynamic range of integer*2
c
      if(.not.fixscale)then
        rhmin=abs(rholim(1,nmap))
        rhmax=abs(rholim(2,nmap))
        if(rhmin.gt.rhmax)rhmax=rhmin
        scale_r=20000.0/amax1( abs(rhmax) ,abs(rhmin) )
      endif
c
c Set limits on xyz
c
      limits(1,iuvw(1,nmap)) = nu1
      limits(2,iuvw(1,nmap)) = nu2
      limits(1,iuvw(2,nmap)) = nv1
      limits(2,iuvw(2,nmap)) = nv2
      limits(1,iuvw(3,nmap)) = nw1
      limits(2,iuvw(3,nmap)) = nw1+nsec-1
c
c Get limits
c
      do i=1,3
        nstart(i,nmap)=limits(1,i)
        nend(i,nmap)=limits(2,i)
        nx(i,nmap)=nend(i,nmap)-nstart(i,nmap)+1
        norg(i,nmap)=1-nstart(i,nmap)
      enddo
c
c Set up map descriptors
c
      mscale(nmap)=scale_r
      moffset(nmap)=0.0
c
      if(.not.terse)then
        write(6,*) 'Scale, offset:',mscale(nmap),moffset(nmap)
        if(record)write(nrec,*)
     &             'Scale, offset:',mscale(nmap),moffset(nmap)
      endif
c
      call XRfrac(XRcell(1,nmap),XRtr(1,1,nmap),XRintr(1,1,nmap),
     &             XRvol,.true.)
c
c---- set up symmetry data
c
      call XR_symsy(XRcell(1,nmap),lgrp(nmap),ics_local,
     &                  n_xsym(nmap),XR_sym(1,nmap) )
      if(verbose)then
        write(6,*)'Number of symm ops for this space group=',
     &                  n_xsym(nmap)
        if(record)write(nrec,*)
     &    'Number of symm ops for this space group=',n_xsym(nmap)
      endif
c
c---- and integerized CSO's:
c
      call XR_isym(nmap)
c
c---- also set up pointers indicating which symm op to try first
c
      do j=1,maxsym
        last_sym(j,nmap)=1
      enddo
c
      npix(nmap)=nx(1,nmap)*nx(2,nmap)*nx(3,nmap)
c
c---- set up bricking data
c
      call brickit(nmap)
c
c Allocate space in scratch array, for map
c
      call alloc(nmap)
c
c How much scratch space do we need for reading?
c
      nline=nx(iuvw(1,nmap),nmap)
      nsize=(nu2-nu1+1)*(nv2-nv1+1)
c
c Allocate some real scratch space
c
      npoint=P2AL(stub(1),stub(2),nsize)
      if(npoint.eq.0) then
        write(6,*)'%READ_MAPB-ERR: Memory allocate fails'
        if(record)write(nrec,*)'%READ_MAPB-ERR: Memory allocate fails'
        goto 999
      else if(verbose)then
        write(6,*)nsize,
     &   ' real elements dynamically allocated for section'
        if(record)write(nrec,*)nsize,
     &   ' real elements dynamically allocated for section'
      endif
c
c Allocate some integer*2 scratch space
c
      np1=P2AL(scratch(1),scratch(2),nline)
      if(np1.eq.0) then
        write(6,*)'%READ_MAPB-ERR: Memory allocate fails'
        if(record)write(nrec,*)'%READ_MAPB-ERR: Memory allocate fails'
        goto 999
      else if(verbose)then
        write(6,*)nline,
     &    ' integer*2 elements dynamically allocated for section'
        if(record)write(nrec,*)nline,
     &    ' integer*2 elements dynamically allocated for section'
      endif
c
c Actually read map
c
      call gulp_mapb(lun,stub(npoint),nsize,scratch(np1),nline,nmap,
     &  scale_r)
c
c CCP4 binary map close, ok to do this for read:
c
c     close(unit=lun)
c
c RME: 23/4/2010
c No, this is absolutely not OK, use the proper call!
c
	call mrclos(lun)
c
c Deallocate scratch space
c
      if(.not.P2DAL(stub(npoint)))then
        write(6,*)'%READ_MAPB-ERR: Memory deallocate fails'
        if(record)write(nrec,*)'%READ_MAPB-ERR: Memory deallocate fails'
        goto 999
      endif
      if(.not.P2DAL(scratch(np1)))then
        write(6,*)'%READ_MAPB-ERR: Memory deallocate fails'
        if(record)write(nrec,*)'%READ_MAPB-ERR: Memory deallocate fails'
        goto 999
      endif
c
c Done !
c
      return
c
c On error...
c
999   write(6,9999)
      if(record)write(nrec,9999)
9999  format(' %READ_MAPB-ERR: **ERROR: Abandon read** ')
      return
c
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE read_pdb(fil)
c       =========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c
c
        integer lun1,i
        CHARACTER*80  line
        CHARACTER*80  charval
        character*6   start
        character*80  fil
        integer*4  opna
        character*5   ornum
        equivalence(line,start)
c
        lun1=1
        P_NA=0
        p_nr=0
        ornum=' '
        opna=0
c
c
c
        write(6,10)
        if(record)write(nrec,10)
10      format(' READ ATOMIC COORDINATES (PDB FORMAT)')
        fil=charval('AV_READ_PDB:ENTER_FILE',1,80,0)
        if (index(fil,'.').eq.0) then
          i=index(fil,' ')
          fil=fil(1:i-1)//'.pdb'
        endif
c
        if(.not.terse)then
          write(6,40)fil(1:60), lun1
          if(record)write(nrec,40)fil(1:60), lun1
40        format(' Pdbfile:', a60,/,' opened on unit=',i3)
          else
          write(6,41)fil(1:60)
          if(record)write(nrec,41)fil(1:60)
41        format(' Pdbfile:', a60)
          endif
c
        close (unit=lun1)
c JMD!PORT: open  (unit=lun1,file=file,status='old',readonly)     
        open  (unit=lun1,file=fil,status='old')
c
100     READ (lun1,21,END=99)line
21      FORMAT(a80)
        if((start.ne.'ATOM  ').and.(start.ne.'HETATM')) goto 100
c
c---- check arrays are big enough
c
        if(P_NA.gt.natom)then
           write(6,*)
     &     '%READ_PDB-ERR: Natom too small-change parameter'
           if(record)write(nrec,*)
     &     '%READ_PDB-ERR: Natom too small-change parameter'
           return
           endif
c
c---  decode line
c
        P_NA=P_NA+1
        read(line,55) p_rectyp(p_na), p_anum(p_na), p_anam(p_na),
     &                p_rnam(p_na), p_rnum(p_na), (p_xyz(i,p_na),i=1,3),
     &                p_occ(p_na), p_bfac(p_na)
55      format(a6,a5,x,a5,a4,a5,4x,3f8.3,2f6.2)
C
c
c--- check whether its a new residue.
c
        if(p_rnum(p_na).ne.ornum)then
          if(p_nr.ge.nres)then
                write(6,*)
     &          '%READ_PDB-ERR: Too many residues'
                if(record)write(nrec,*)
     &          '%READ_PDB-ERR: Too many residues'
           close(lun1)
           return
          endif
c
c---record atoms per residue
          p_nr=p_nr+1
           if(p_nr.ne.1)then
             p_ratm(p_nr-1)=p_na-opna
           endif
          ornum=p_rnum(p_na)
          opna=p_na
        endif
c
        goto100
c
c--- record atom for last residue
99      p_ratm(p_nr)=(p_na+1)-opna

        write(6,70) p_na, p_nr
        if(record)write(nrec,70) p_na, p_nr
70      format(' No of atoms read =',i6,' in',i6,' residues')
        close(unit=lun1)
C
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE reco(flag)
c       =====================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c       =====================
c
c
        character*80 fil,charval
c
        logical flag
c
        if(record) then
          write(6,30)
          write(nrec,30)
30        format(' CURRENT RECORD FILE WILL NOW BE CLOSED')
          endif
c
        record=flag
        close (unit=nrec)
c
        if(.not.record)return
c
        fil=charval('AV_RECORD:ENTER_FILE',1,80,0)
        open (unit=nrec,file=fil,status=cstatus)
        if(.not.terse)then
          write(6,40)fil(1:60), nrec
40        format(' RECORD FILE:', a60,/,' opened on unit=',i3)
          else
          write(6,41)fil(1:60)
41        format(' RECORD FILE:', a60)
          endif
          if(record)write(nrec,50)
50        format(' THIS IS A RECORD FILE PRODUCED BY GAP',/,
     &             ' =====================================')
        call header(.false.)
        return
        end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE reconst
c     ==================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c---  New option Feb 1996, JMG. Reconstruct the averaged map
c---  from the averaged protomer.  This uses a labelled envelope
c---  which should have been generated up already using lenv option.
c---  This has the associated sym operator no as the value of the pixel
c---  in the envelope.
c---  The input map is overwritten. The program will go through
c---  the defined envelope and for pixel values > 1 use the inverse
c---  ncs operator read from the input envelope to get interpolated
c---  density from the averaged part of the map.
c---
c---  if env1 = 0 we are outside the envelope
c---  if env1 = 1 we are in the envelope but don't need to do anything
c---  if env1 > 1 we are in the envelope but need to get some density
c---
c---  XTAL symmetry should NOT be used as you will almost always
c---  get into great trouble.
c
c--- Local variables
c
      integer     map1,env1,iindex(3),jjx,jjy,jjz,iz,jz,iy,jy,ix,jx,
     &            nopt,keyword,intval,i_env,itype,nsym1,nsym2,ngadd,
     &            np,nc,nout
      integer*2   gp3,irho
      real        aindex(3),x(3),x1(3)
      character*8 G2CHRU,imap,ienv
      logical     lass_get,use_in
c
c--- Defaults
c
      itype=11
      nsym1=1
      nsym2=nsym
      use_in = .true.
c
c--- Check if MAPI has an assignment - if not then get one
c--- Afterwards imap contains the input map name
c
      if(.not.lass_get(imap,'MAPI'))then
        imap=G2CHRU('AV_FOLD:INPUT_MAP',1,8,0)
        call lass_set(imap,'MAPI')
      endif
c
c--- Check if ENVI has an assignment - if not then get one
c--- Afterwards ienv contains the input envelope name
c
      if(.not.lass_get(ienv,'ENVI'))then
        ienv=G2CHRU('AV_FOLD:INPUT_ENV',1,8,0)
        call lass_set(ienv,'ENVI')
      endif
c
c--- Report on the state of play
c
      if(verbose)then
        write(6,110) imap,ienv
        if(record)write(nrec,110) imap,ienv
110     format(/' RECONSTRUCT MAP',/,
     &          ' ==============='
     &         /' Input map:       ',a,
     &         /' Input envelope:  ',a)
      else
        write(6,111) imap,ienv
        if(record)write(nrec,111) imap,ienv
111     format(' RECONSTRUCT MAP',
     &    '  I/P map:',a,'  I/P env:'a)
      endif
c
c--- Get the slot number of the input map
c
      call av_ass(map1,imap)
c
c--- The input map should exist by now - if not then run away
c
      if(.not.defined(map1)) then
        write(6,160)
160     format(' %RECONST_MAP-ERR: Input map empty')
        return
      endif
c
c--- Get the slot number for the input envelope
c
      call av_ass(env1,ienv)
c
c--- The input envelope should exist by now - if not the run away
c
      if(.not.defined(env1)) then
        write(6,170)
170     format(' %RECONST_ENV-ERR: Input env empty')
        return
      endif
c
c--- Get the interpolation value to be used
c
10    nopt=keyword('AV_UNPA','GO  !INTP',0)
c
c---- INTP - interpolation method
c
      if(nopt.eq.2)then
        itype=intval('AV_UNPA_INTP',1,64,0)
        if(itype.ne.11.and.itype.ne.8.and.itype.ne.1.and.itype.ne.64)
     &  then
          itype=11
          print*,'%AV_UNPA_ERR: Invalid INTP - INTP reset to 11'
        endif
      endif
c
c---- GO - wait fo it
c
      if(nopt.ne.1)goto 10
c
c--- Report which symm ops we are using - all of them!
c
      write(6,220)nsym1,nsym2
      if(record)write(nrec,220)nsym1,nsym2
220   format(' Will use symmetry operators from ',i4,' to',i4)
c
c--- Warn if xtal symmetry switched on....very very dangerous
c
      if(xr_con(map1))then
        write(6,240)
        if(record)write(nrec,240)
240     format(
     &    ' %WARNING-ERR: **CRYSTALLOGRAPHIC SYMMETRY** switched on')
      endif
c
c--- Report interpolation scheme
c
      write(6,250)itype
      if(record)write(nrec,250)itype
250   format(i3,' point interpolation')
c
c--- Set up pointers for the maps that drive the process
c
      call init_map(map1)
      call init_map(env1)
c
c--- Initialize counter
c
      ngadd=0
      nout=0
c
c--- Loop over all bricks in map1
c
      do iz=1,n_brick(3,map1)
       jz=(iz-1)*brick(3,map1)
c
       do iy=1,n_brick(2,map1)
        jy=(iy-1)*brick(2,map1)
c
        do ix=1,n_brick(1,map1)
         jx=(ix-1)*brick(1,map1)
c
c---- Now loop within each brick
c
         do jjz=jz+1,jz+brick(3,map1)
          if(jjz.le.nx(3,map1))then
c
           do jjy=jy+1,jy+brick(2,map1)
            if(jjy.le.nx(2,map1))then
c
             do jjx=jx+1,jx+brick(1,map1)
              if(jjx.le.nx(1,map1))then
c
c----- Get pixel value from env1 - if 0 next pixel, if 1 update count
c----- and the next pixel, if > 1 get on with it
c
               i_env=gp3(env1,jjx,jjy,jjz)
               if(i_env.eq.1)nout=nout+1
               if(i_env.le.1)goto 200
c
c----- Put the co-ordinates of this pixel in an array
c
               iindex(1)=jjx
               iindex(2)=jjy
               iindex(3)=jjz
c
c----- Orthagonalise the co-ordinates
c
               call j3_to_f(iindex,x1,map1)
               call f_to_o(x1,x,map1)
c       
c----- Use inverse operator i_env to get the appropriate copy
c
               call spin_inv(x,x1,ops_inv(1,1,i_env),vecs(1,i_env))
               call o_to_j3(x1,aindex,iindex,map1)
c
c----- Now get the interpolated electron density - much quicker if
c----- not using crystallographic symmetry
c
               if(xr_con(map1))then
                call g_p_i_i_x(map1,irho,aindex,iindex,itype,i_env,
     &                                    use_in,np,nc)
               else
                call g_p_i_i(map1,irho,aindex,iindex,itype,use_in,np,nc)
               endif
c
c----- Put interpolated density into this pixel 
c
               call pp3(map1,jjx,jjy,jjz,irho)
               ngadd=ngadd+1
c
c----- Next pixel please
c
200            continue
              endif
             enddo
            endif
           enddo
          endif
         enddo
        enddo
       enddo
      enddo
c
c--- Report on what happened
c
c999   continue
      if(verbose)then
        write(6,998)ngadd,nout
        if(record)write(nrec,998)ngadd,nout
998     format(' MAP RECONSTRUCTION FINISHED',/,
     &         ' ==========================',/,
     &      i9,' pixels were reconstructed from ', i9 ,' pixels')
      else
        write(6,997)ngadd,nout
        if(record)write(nrec,997)ngadd,nout
997     format(' Reconstructed',i9,' pixels from', i9 ,' pixels')
      endif
c
c--- Bye!
c
      return
      end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      FUNCTION rem(j3,i3,nb)
c       ======================
      IMPLICIT NONE
        integer rem,j3,i3,nb
c
        rem=j3-(i3-1)*nb
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE reset
c      ================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c      =====================
c
c--- local variables
c
       integer     map1,env1
       character*8 G2CHRU,imap,ienv
       logical     i_check,same,same_grid,lass_get
c
c--- Get input map, if we don't already have it
c
       if(.not.lass_get(imap,'MAPI'))then
         imap=G2CHRU('AV_RESET:INPUT_MAP',1,8,0)
         call lass_set(imap,'MAPI')
       endif
c
c--- Get input envelope, if we don't already have it
c
       if(.not.lass_get(ienv,'ENVI'))then
         ienv=G2CHRU('AV_RESET:INPUT_ENV',1,8,0)
         call lass_set(ienv,'ENVI')
       endif
c
c--- Report
c
       if(verbose)then
         write(6,100) imap,ienv
         if(record)write(nrec,100) imap,ienv
       else
         write(6,110) imap,ienv
         if(record)write(nrec,110) imap,ienv
       endif
c
100    format(/' RESET RHOLIM',/,
     &           ' ============'
     &      /' Input map:       ',a,
     &      /' Input envelope:  ',a)
110    format(' RESET RHOLIM',
     &    '  I/P map:',a,'  I/P env:',a)
c
c--- Check to see if we have an envelope
c
       i_check=ienv.ne.'OFF'
       if(verbose)then
         if(i_check) write(6,120)
         if(.not.i_check)write(6,130)
c
         if(record) then
           if(i_check) write(nrec,120)
           if(.not.i_check)write(nrec,130)
         endif
       endif
c
120    format(' Input map envelope filtered')
130    format(' Input map not filtered')
c
c--- Get slot numbers for all maps - map 1 should exist by now !
c
       call av_ass(map1,imap)
       if(.not.defined(map1)) then
         write(6,160)
         if(record)write(nrec,160)
         goto 700
       endif
160    format(' %RESET-ERR: Input map empty')
c
c--- Get slot number for envelope and check grid
c
       if(i_check)then
         call av_ass(env1,ienv)
         same=same_grid(map1,env1)
c
         if(.not.same)then
           write(6,170)
           if(record)write(nrec,170)
           goto 700
         endif
       endif
170    format(' %RESET-ERR: MAPI and ENVI not on same grid')
c
c--- Off to business end
c
       call reset_av(map1,env1,i_check)
c
700    return
       end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE reset_av(map1,env1,i_check)
c      ======================================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c      =====================
c
c--- Local variables:
c
       integer   map1,env1,ngadd,iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
       logical   i_check
       integer*2 irho
       external  gp3,gp3_a
       integer*2 gp3,gp3_a
       real      r_max,r_min
c
c--- Set up pointers for the maps that drive the process
c
       call init_map(map1)
       if(i_check)call init_map(env1)
c
c--- Zero counters
c
       ngadd=0
       call r_stat
c
c
c
c--- Off to work...
c--- =============
c
c--- Loop over bricks
c
       do iz=1,n_brick(3,map1)
c
         jz=(iz-1)*brick(3,map1)
c
         do iy=1,n_brick(2,map1)
c
           jy=(iy-1)*brick(2,map1)
c
           do ix=1,n_brick(1,map1)
c
             jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
             do jjz=jz+1,jz+brick(3,map1)
c
             if(jjz.le.nx(3,map1))then
c
               do jjy=jy+1,jy+brick(2,map1)
c
               if(jjy.le.nx(2,map1))then
c
                 do jjx=jx+1,jx+brick(1,map1)
c                    
                 if(jjx.le.nx(1,map1))then
c
c---- check that this pixel is in required region of the output map
c
                   if(i_check)then
                     if(gp3(env1,jjx,jjy,jjz).eq.0)goto 200
                   endif
c
                   irho=gp3(map1,jjx,jjy,jjz)
c
                   call r_a_stat(irho)
                   ngadd=ngadd+1
c
200                continue
                 endif
                 enddo
               endif
               enddo
             endif
             enddo
           enddo
         enddo
       enddo
c
c--- OK, lets see how we did ...
c
c300    continue
       if(verbose)then
         write(6,310)ngadd
         if(record)write(nrec,310)ngadd
       else
         write(6,320)ngadd
         if(record)write(nrec,320)ngadd
       endif
c
310    format(' RHOLIMS RECALCULATED',/,
     &           ' ====================',/,
     &             i9,' pixels were used from i/p map')
320    format(' RHOLIMS RECALCULATED,',i10,
     &     ' pixels were used from i/p map')
c
c--- Get stats from counting subroutine
c
       call r_p_stat(r_max,r_min)
c
c--- Print a summary of the statistics
c
       write(6,330)r_max,r_min,av_temp,sigm
       if(record) then
         write(nrec,330)r_max,r_min,av_temp,sigm
       endif
c
330    format(' New statistics on map pixels will replace old RHOLIM'/,
     &     ' Max rho',f8.1,
     &  '    min rho',f8.1,
     &  '    mean rho',f8.1,
     &  '    sigma',f8.1)
c
c--- Update header without scalefactor
c
       write(6,340)r_max/mscale(map1),r_min/mscale(map1),
     &              av_temp/mscale(map1),sigm/mscale(map1)
       if(record) then
          write(nrec,340)r_max/mscale(map1),r_min/mscale(map1),
     &              av_temp/mscale(map1),sigm/mscale(map1)
       endif
c
340    format(' Header will be updated without current scalefactor'/,
     &     ' Max rho',f17.8,
     &  '    min rho',f17.8,
     &  '    mean rho',f17.8,
     &  '    sigma',f17.8)
c
       rholim(1,map1)=r_min/mscale(map1)
       rholim(2,map1)=r_max/mscale(map1)
       rholim(3,map1)=av_temp/mscale(map1)
       rholim(4,map1)=sigm/mscale(map1)
c
c--- And its all over...
c
       return
       end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
      SUBROUTINE RIGSET
C     =================
      IMPLICIT NONE
C----                                                                         
C---- SUBROUTINE TO FIT THE COORD SET XA(3,N) TO THE SET XB(3,N)            
C---- IN THE SENSE OF:                                                      
C----        XA= R*XB +V                                                    
C---- R IS A UNITARY 3.3 RIGHT HANDED ROTATION MATRIX                       
C---- AND V IS THE OFFSET VECTOR. THIS IS AN EXACT SOLUTION                 
C---- 
C---- IF ENTRY_L IS LOGICALLY FALSE ONLY THE RMS COORDINATE ERROR             
C---- WILL BE RETURNED (BUT QUICKLY)                                    
C---- 
C---- THIS SUBROUTINE IS A COMBINATION OF MCLACHLAN'S AND KABSCH'S           
C---- TECHNIQUES. SEE                                                        
C----   KABSCH, W. ACTA CRYST A34, 827,1978                                   
C----   MCLACHAN, A.D., J. MOL. BIO. NNN, NNNN 1978                           
C----   WRITTEN BY S.J. REMINGTON 11/78.                                      
C---- 
C---- THIS SUBROUTINE USES THE IBM SSP EIGENVALUE ROUTINE 'EIGEN'           
C----                                                                         
      DIMENSION  R(3,3),V(3)
      DIMENSION  XTEMPA(3),XTEMPB(3)
      REAL WT_TEMP,R,XTEMPA,XTEMPB,RMSE,V
      DOUBLE PRECISION  XA(3),XB(3),WT
      INTEGER I,J,K
      LOGICAL ENTRY_L
      DOUBLE PRECISION CMA(3),CMB(3),UMAT(3,3),
     & XN,XASQ,XBSQ,T,XNI,RTSUM
C
C JMD Port: This subroutine appears to need a SAVE statement that was previously
C unnecessary (or it was broken)
C
      SAVE XN,XASQ,XBSQ,CMA,CMB,UMAT
C----
C----
C**** FIRST ENTRY POINT, INITIALIZE.
C----
      XN=0.0
      XASQ=0.0
      XBSQ=0.0
      DO 10 I=1,3
        CMA(I) = 0.0
        CMB(I) = 0.0
        DO 10 J=1,3
   10     UMAT(I,J) = 0.0
C----
      RETURN
C----
C----
C**** SECOND ENTRY POINT, ADD AN OBSERVATION TO THE
C**** EQUATIONS.
C----   WT IS THE WEIGHT TO BE GIVEN TO THIS RELATIONSHIP.
C----
      ENTRY  RIGADD ( XTEMPA, XTEMPB, WT_TEMP )
C----
C----   STORE INTERNALLY AS DOUBLE PRECISION
C----
      DO J=1,3
        XA(J) = XTEMPA(J)
        XB(J) = XTEMPB(J)
      ENDDO
      WT=WT_TEMP
C----    
C----   ACCUMULATE UNCORRECTED (FOR C.M.) SUMS AND SQUARES  
C----  
      XN=XN+WT
C----  
      DO 30 I=1,3
        DO 20 K=1,3
          UMAT(I,K) = UMAT(I,K) + XA(I)*XB(K)*WT
   20   CONTINUE
C----  
        T = XA(I)*WT
        CMA(I) = CMA(I) + T
        XASQ = XASQ + T*T/WT
        T = XB(I)*WT
        CMB(I) = CMB(I) + T
        XBSQ = XBSQ + T*T/WT
   30 CONTINUE
      RETURN
C----
C----
C----
C*****THIRD ENTRY POINT ... SOLVE.
C----
      ENTRY RIGSOL(  R,  V,  RMSE,  ENTRY_L  )
C----                                                                    
C----   SUBTRACT CM OFFSETS                                              
C----                                                                    
      XNI=1.0/XN
      DO 50 I=1,3
        XASQ = XASQ - CMA(I)*CMA(I)*XNI
        XBSQ = XBSQ - CMB(I)*CMB(I)*XNI
        DO 50 J=1,3
          UMAT(I,J) = (UMAT(I,J) - CMA(I)*CMB(J)*XNI)*XNI
   50   CONTINUE
C----                                                                  
C----   FIT IT                                                         
C----                                                                  
      CALL QKFIT( UMAT, RTSUM, R, ENTRY_L )
      RMSE =(XASQ + XBSQ)*XNI - 2.0*RTSUM
      IF( RMSE .LT. 0.0 ) RMSE = 0.0
      RMSE = SQRT (RMSE)
C----                                                  
C---- CALCULATE OFFSET IF ENTRY_L=.TRUE.              
C----                                                  
      IF(.NOT. ENTRY_L) RETURN
C----                                                  
      DO 70 I=1,3
        T = 0.0
        DO 60 J=1,3
   60     T = T + R(I,J)*CMB(J)
        V(I) = ( CMA(I) - T )*XNI
   70 CONTINUE
      RETURN
      END

C=============================================================================

      SUBROUTINE ROTAXIS(ROT, VEC, ANG)
      IMPLICIT NONE

C Convert a rotation matrix into a direction vector for rotation axis plus
C the rotation angle. (See ~robert/programs/rotsearch/rot.show)

        DOUBLE PRECISION ROT(3,3)
        DOUBLE PRECISION VEC(3)
        DOUBLE PRECISION CS, ANG

C The rotation angle is acos( (Tr(R)-1) / 2)

        CS = (ROT(1,1)+ROT(2,2)+ROT(3,3)-1.0) / 2.0
        IF (CS.LT.-1.0) CS = -1.0
        IF (CS.GT. 1.0) CS =  1.0
        ANG = ACOS(CS)

C The pairs of off-diagonal terms can be subtracted to give simple expressions
C for the components of the direction vector.

        IF (CS.EQ.1.0) THEN
          VEC(1) = 1.0
          VEC(2) = 0.0
          VEC(3) = 0.0
        ELSE IF (CS.EQ.-1.0) THEN
          VEC(1) = (ROT(1,1) + 1.0) / 2.0
          IF (VEC(1).LT.0.0) VEC(1) = 0.0
          VEC(1) = SQRT(VEC(1))
          VEC(2) = (ROT(2,2) + 1.0) / 2.0
          IF (VEC(2).LT.0.0) VEC(2) = 0.0
          VEC(2) = SQRT(VEC(2))
          VEC(3) = (ROT(3,3) + 1.0) / 2.0
          IF (VEC(3).LT.0.0) VEC(3) = 0.0
          VEC(3) = SQRT(VEC(3))
          IF (VEC(1).GT.0.5) THEN
            IF (ROT(1,2).LT.0.0) VEC(2) = -VEC(2)
            IF (ROT(1,3).LT.0.0) VEC(3) = -VEC(3)
          ELSE IF (VEC(2).GT.0.5) THEN
            IF (ROT(2,1).LT.0.0) VEC(1) = -VEC(1)
            IF (ROT(2,3).LT.0.0) VEC(3) = -VEC(3)
          ELSE IF (VEC(3).GT.0.5) THEN
            IF (ROT(3,1).LT.0.0) VEC(1) = -VEC(1)
            IF (ROT(3,2).LT.0.0) VEC(2) = -VEC(2)
          ENDIF
        ELSE
          VEC(1) = (ROT(2,3) - ROT(3,2)) / (2.0 * SIN(ANG))
          VEC(2) = (ROT(3,1) - ROT(1,3)) / (2.0 * SIN(ANG))
          VEC(3) = (ROT(1,2) - ROT(2,1)) / (2.0 * SIN(ANG))
        ENDIF

        CALL NORM(VEC, VEC)

        RETURN
        END
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE ROTMAT(rROT,rT1,rT2,rT3,rAXIS,MODE)
c     ==============================================
      IMPLICIT NONE
c
c--- routine cleaved from xplor....
c--- uses d.precision within routine
C
C Routine computes unitary rotation matrix ROT using Eulerian angles 
C (MODE="EULE"), Lattman angles (MODE="LATT"), spherical polar angles 
C (MODE="SPHE") or a rotation about the specified axis (MODE="AXIS").  
C
C Input: 
C    MODE specifies angle mode
C    T1,T2,T3 are theta1 (z), theta2 (x'), theta3 (z') for MODE="EULE"
C    T1,T2,T3 are theta+, theta2, theta- for MODE="LATT"
C    T1,T2,T3 are psi (incl. vs. y), phi (azimuthal angle, that is,
C             the angle between the x-axis and the projection of the 
C             axis into the x,z plane ), kappa for MODE="SPHE"
C    T3, AXIS(3) are kappa and a 3-D vector specifying the axis for MODE="AXIS"
C  Note: all rotations are counter-clockwise
C
C Output:
C    ROT(3,3) contains the rotation matrix.  Should be applied as 
C    r'(i)=sum_j ROT(i,j)*r(j)
C
C    AXIS contains the normalized input AXIS vector and T1,T2,T3 will
C    contain the spherical polar angles for MODE="SPHE"
C
C Author: Axel T. Brunger
C =======================
C
c
C I/O
      INCLUDE 'average.fcm'
c     ====================
c
        integer i,j
        real rrot(3,3),rt1,rt2,rt3,raxis(3)
      DOUBLEPRECISION ROT(3,3), T1, T2, T3, AXIS(3)
      CHARACTER*4 MODE
C local
      DOUBLEPRECISION S1, S2, S3, C1, C2, C3, C3C, S1SQ, NN, TP, TM
      DOUBLEPRECISION TT1, TT3
C parameter
      DOUBLEPRECISION RAD, ONE, TWO, SMALL
        doubleprecision rsmall
c--- parameter rsmall stored in consta.fcm
      PARAMETER (RAD=0.017453292D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (SMALL=1.0D-7)
        PARAMETER (RSMALL=1.0D-10)
C begin
C--- change REAL INTO DOUBLE PRECISION
C
        T1=RT1
        T2=RT2
        T3=RT3
C
        do j=1,3
        axis(j)=raxis(j)
        enddo
C
        DO I=1,3
                DO J=1,3
                  RROT(I,J)=ZERO
                ENDDO
        ENDDO
C
c
C=================================================================
      IF (MODE.EQ.'EULE'.OR.MODE.EQ.'LATT') THEN
C
C Computes rotation matrix corresponding to the three
C Eulerian angles theta1=T1, theta2=T2, theta3=T3.  This uses 
C the Rossmann and Blow convention, i.e., a theta1 rotation 
C around the Z axis, followed by
C a theta2 rotation around the moved X axis, and followed by a 
C theta3 rotation around the moved Z axis.  These angles are positive
C if they are anti-clockwise when looking along the relevant axis.  
C
      IF (MODE.EQ.'LATT') THEN
C In "Lattman" mode, theta+=theta1+theta3=T1, theta2=T2, 
C theta-=theta1-theta3=T3
      TP=T1
      TM=T3
      TT1=(TP+TM)/TWO
      TT3=(TP-TM)/TWO
      ELSE
      TT1=T1
      TT3=T3
      END IF
C
      S1=SIN(RAD*TT1)
      S2=SIN(RAD*T2)
      S3=SIN(RAD*TT3)
      C1=COS(RAD*TT1)
      C2=COS(RAD*T2)
      C3=COS(RAD*TT3)
      ROT(1,1)=-S1*C2*S3+C1*C3
      ROT(1,2)=C1*C2*S3+S1*C3
      ROT(1,3)=S2*S3
      ROT(2,1)=-S1*C2*C3-C1*S3
      ROT(2,2)=C1*C2*C3-S1*S3
      ROT(2,3)=S2*C3
      ROT(3,1)=S1*S2
      ROT(3,2)=-C1*S2
      ROT(3,3)=C2
C=================================================================
      ELSE IF (MODE.EQ.'SPHE'.OR.MODE.EQ.'AXIS') THEN
C
      IF (MODE.EQ.'AXIS') THEN
C
C In axis mode we obtain the psi, phi spherical polar angles from 
C the AXIS vector.
C The rotation kappa corresonds to an anticlockwise rotation (t3) about 
C the specified axis. The axis is desribed as a vector (AXIS).  
      NN=SQRT(AXIS(1)**2+AXIS(2)**2+AXIS(3)**2)
      IF (NN.LT.RSMALL) THEN
      write(6,*) 'Rotation vector has zero length - angles set to zero'
      if(record)write(nrec,*)
     &           'Rotation vector has zero length - angles set to zero'
      T1=ZERO
      T2=ZERO
      T3=ZERO
      ELSE
      AXIS(1)=AXIS(1)/NN
      AXIS(2)=AXIS(2)/NN
      AXIS(3)=AXIS(3)/NN
      T1=ACOS(AXIS(2))/RAD
      IF (AXIS(1)**2+AXIS(3)**2.LT.SMALL) THEN
      T2=ZERO
      ELSE
      T2=ACOS(dmax1(-ONE,dmin1(ONE,AXIS(1)/SQRT(AXIS(1)**2+
     &    AXIS(3)**2))))/RAD
      IF (AXIS(3).GT.ZERO) T2=-T2
      END IF
      END IF
      END IF
C
C compute rotation matrix corresponding to the three
C spherical polar angles psi=t1, phi=t2 and kappa=t3.  The rotation is 
C described by specification of the direction of an axis through
C phi (azimutal angle between the x axis and the projection of the
C rotation axis on the x-z plane) and psi (inclination versus y axis).  
C The angle kappa specifies the rotation around the specified axis. 
C The kappa angle is anti-clockwise when looking along the rotation axis.
C The phi angle is anti-clockwise when looking along y. 
      S1=SIN(RAD*T1)
      S2=SIN(RAD*T2)
      S3=SIN(RAD*T3)
      C1=COS(RAD*T1)
      C2=COS(RAD*T2)
      C3=COS(RAD*T3)
      S1SQ=S1*S1
      C3C=1.0-C3
      ROT(1,1) = C3           + S1SQ*C2*C2*C3C
      ROT(1,2) = S1*C1*C2*C3C   - S1*S2*S3
      ROT(1,3) =-S1SQ*C2*S2*C3C - C1*S3
      ROT(2,1) = S1*C1*C2*C3C   + S1*S2*S3
      ROT(2,2) = C3           + C1*C1*C3C
      ROT(2,3) =-S1*C1*S2*C3C   + S1*C2*S3
      ROT(3,1) =-S1SQ*S2*C2*C3C + C1*S3
      ROT(3,2) =-S1*C1*S2*C3C   - S1*C2*S3
      ROT(3,3) = C3           + S1SQ*S2*S2*C3C
C================================================================
      END IF
C---CONVERT FROM DO TO SINGLE
        DO I=1,3
                DO J=1,3
                  RROT(I,J)=ROT(I,J)
c                 RAXIS(J)=AXIS(J)
                ENDDO
        ENDDO
c
      RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE rtnc
c     ===============
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c     =====================
c
      double precision R_local(3,3),rinv(3,3),det
      double precision R_save(3,3),d_local(3,3)
      real    t_save(3),t(3),tinv(3),r_save_s(3,3),rinv_s(3,3),realval
      integer i,j,k,j_nsym
c
      do i=1,3
        do j=1,3
          r_local(i,j)=
     &      dble(realval('RTNC:ENTER_ROT_MAT',-100.0,100.0,0))
        enddo
      enddo
c
      do i=1,3
c JMD!PORT: sngl unnecessary
c        t(i)=sngl(realval('RTNC:ENTER_T_VECTOR',-100000.0,100000.0,0))
        t(i)=realval('RTNC:ENTER_T_VECTOR',-100000.0,100000.0,0)
      enddo
c
c---    nsym number of ncs ops to be rotated
c
c       call invers_d(r_local,rinv,det)
c
c       if(abs(abs(det)-1.0).gt.0.99)then
c         write(6,*)'**WARNING** det .ne. 1'
c         if(record)write(nrec,*)'**WARNING** det .ne. 1'
c         endif
c
      do j=1,3
        tinv(j)=-t(j)
      enddo
      call invers_d(r_local,rinv,det)
      if(abs(abs(det)-1.0).gt.0.99)then
        write(6,*)'**WARNING** det .ne. 1'
        if(record)write(nrec,*)'**WARNING** det .ne. 1'
      endif

      do j_nsym=1,nsym

        call d_copy(ops(1,1,j_nsym),d_local,9)
        call matmul_d(d_local,r_local,r_save)
c---check
        write(6,100)j_nsym,((r_save(j,k),k=1,3),j=1,3),t
c
        call matmul_d(rinv,r_save,r_save)
c---check
        write(6,100)j_nsym,((r_save(j,k),k=1,3),j=1,3),t
c
        do i=1,3
          do j=1,3
            r_save_s(i,j)=sngl(r_save(i,j))
            rinv_s(i,j)=sngl(rinv(i,j))
          enddo
        enddo

        call vmtply(r_save_s,tinv,t_save)
c
        do j=1,3
          t_save(j)=t(j)+t_save(j)
        enddo
c
        call copy(9,r_save_s,ops(1,1,j_nsym))
        call copy(9,rinv_s,ops_inv(1,1,j_nsym))
        call copy(3,t_save,vecs(1,j_nsym))
        if(.not.terse)then
          write(6,100)j_nsym,((r_save(j,k),k=1,3),j=1,3),t_save
          if(record)write(nrec,100)j_nsym,((r_save(j,k),k=1,3),j=1,3),
     &      t_save
100       format(10x,i6,' Symm op is now',/,2(3(10x,3f11.5,/),
     &           10x,3f11.2))
        endif
      enddo
c
c
      return
      end

c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c       
      FUNCTION rtrue(irs,off,scal)
c       ============================
      IMPLICIT NONE
c
        real rtrue,scal
        integer irs,off
        rtrue=(irs-off)/scal
        return
        end
c
      FUNCTION same_grid(map1,map2)
c       ============================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c---logical fuction returns true if maps on same grid
c
        logical same_grid
        integer map1,map2,j
c
        same_grid=.true.
c
        do j=1,3
         if( nunit(j,map1).ne.nunit(j,map2) ) same_grid=.false.
        enddo
c
c
        return
        end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      FUNCTION same_map(map1,map2)
c       ============================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c---logical fuction returns true if maps same size
c
        logical same_map
        integer map1,map2,j
c
        same_map=.true.
c
        do j=1,3
         if( ( nstart(j,map1).ne.nstart(j,map2) ) .or.
     &      ( nend(j,map1).ne.nend(j,map2) ) )same_map=.false.
        enddo
c
c
        return
        end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE scale_map
c       ====================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c       =====================
c
c--- local variables
c
        integer     map1,map2,env1,env2
        character*8 G2CHRU,ch_temp,imap,omap,ienv,oenv
        logical     i_check,o_check
        real        rms1,rms2,av_dens1,av_dens2
        integer     nopt,keyword,no_maps
        logical     lass_get,same_grid,same
c
c
c
        no_maps=0
        av_dens1=0.0
        av_dens2=0.0
c
c--- setup whether scale map up using defined values or 2 maps together
c
20      nopt=keyword
     &  ('AV_MAP_SCALE','1MAP!2MAP!BACK!GO  ',0)
c
        if(nopt.eq.1)then
        no_maps=1
        endif
c
        if(nopt.eq.2)then
        no_maps=2
        endif
c
        if(nopt.eq.3) goto 1000
c
        if( (nopt.ne.4).or.(no_maps.eq.0) ) goto 20
c
c
c
c---- for just applying scalin to one map
c
        if(no_maps.eq.1)then
c
            if(.not.lass_get(ch_temp,'MAPI'))then
c
              ch_temp=G2CHRU('AV_SCAL:INPUT_MAP',1,8,0)
c
              call lass_set(ch_temp,'MAPI')
            endif
            imap=ch_temp
c
c
          if(verbose)then
            write(6,110) imap
            if(record)write(nrec,110) imap
110         format(/' SCALE A MAP',/,
     &       ' ============'
     &      /' Input map:       ',a)
            else
            write(6,112) imap
          if(record)write(nrec,112) imap
112       format(' SCALE A MAP',
     &    '  I/P map:',a)
          endif
C
          i_check=.false.
          env1=0
c
c
c       get slot numbers for all maps
c
            call av_ass(map1,imap)
c----map 1 should exist by now !
             if(.not.defined(map1)) then
               write(6,160)
               if(record)write(nrec,160)
160            format(' %SCALE_MAP-ERR: Input map empty')
              goto 1000
             endif
c
c
c---- set up allocations for envelopes
c
c--- av_temp sigm declared in commomn block...use f_stat for compiling
c--- statistics.
c
                write(6,170) shi_map, sca_map
                 if(record)write(nrec,170)
170              format(' Apply a shift of ',f8.1,' and
     &                 scalefactor of ',f6.1,' to map.' )

           call loopover(map1,env1,i_check,.true.,av_dens1)

c---  this ends scalin for 1 map option so return to base level

         goto 1000
        endif
c
c
c
c
        if(no_maps.eq.2)then
c
              if(.not.lass_get(ch_temp,'MAPI'))then
c
                 ch_temp=G2CHRU('AV_SCAL:INPUT_MAP',1,8,0)
c
                  call lass_set(ch_temp,'MAPI')
              endif
                 imap=ch_temp
c
c
              if(.not.lass_get(ch_temp,'MAPO'))then
c
                 ch_temp=G2CHRU('AV_SCAL:OUTPUT_MAP',1,8,0)
c
                  call lass_set(ch_temp,'MAPO')
              endif
                  omap=ch_temp
c
c
              if(.not.lass_get(ch_temp,'ENVI'))then
c
                 ch_temp=G2CHRU('AV_SCAL:INPUT_ENV',1,8,0)
c
                  call lass_set(ch_temp,'ENVI')
              endif
                  ienv=ch_temp
c
c
              if(.not.lass_get(ch_temp,'ENVO'))then
c
                 ch_temp=G2CHRU('AV_SCAL:OUTPUT_ENV',1,8,0)
c
                 call lass_set(ch_temp,'ENVO')
              endif
                  oenv=ch_temp
c
c
              if(verbose)then
                write(6,210) imap,omap,ienv,oenv
                if(record)write(nrec,210) imap,omap,ienv,oenv
210             format(/' SCALE TWO MAPS',/,
     &           ' =============='
     &          /' Input map:       ',a,
     &          /' Output map:      ',a,
     &          /' Input envelope:  ',a,
     &          /' Output envelope: ',a)
                else
                write(6,212) imap,omap,ienv,oenv
              if(record)write(nrec,212) imap,omap,ienv,oenv
212           format(' SCALE TWO MAPS',/,
     &        ' I/P map:',a,'  O/P map:',a,
     &        '  I/P env:',a,'  O/P env:',a)
              endif
C
              i_check=ienv.ne.'OFF'
              o_check=oenv.ne.'OFF'
c
c
              if(verbose)then
                if(i_check) write(6,220)
                if(.not.i_check)write(6,230)
                if(o_check) write(6,240)
                if(.not.o_check)write(6,250)
c
                if(record) then
                  if(i_check) write(nrec,220)
                  if(.not.i_check)write(nrec,230)
                  if(o_check) write(nrec,240)
                  if(.not.o_check)write(nrec,250)
                  endif
                endif
c
220            format(' Input map envelope filtered')
230            format(' Input map not filtered')
240            format(' Output map envelope filtered')
250            format(' Output map not filtered')
c
c
c       get slot numbers for all maps
c
              call av_ass(map1,imap)
c----map 1 should exist by now !
               if(.not.defined(map1)) then
                 write(6,260)
260              format(' %SCALE_MAP-ERR: Input map1 empty')
                  if(record) write(nrec,260)
                 return
               endif
c---- check status of output file
              call av_ass(map2,omap)
c----map 2 should exist by now !
              if(.not.defined(map2)) then
                write(6,270)
                  if(record) write(nrec,270)
270             format(' %SCALE_MAP-ERR: Input map2 empty')
                return
              endif
c
c---- set up allocations for envelopes
c
                 if(i_check) call av_ass(env1,ienv)
                 if(o_check) call av_ass(env2,oenv)
c
c
c
               if(i_check)then
                 same=same_grid(map1,env1)
c
                  if(.not.same)then
                    write(6,280)
                    if(record)write(nrec,280)
280      format(' %SCALE_MAP-ERR: Envi and mapi not on same grid')
                    goto 1000
                  endif
               endif
c
               if(o_check)then
                 same=same_grid(map2,env2)
c
                 if(.not.same)then
                  write(6,281)
                  if(record)write(nrec,281)
281           format(
     &        ' %SCALE_MAP-ERR: Envo and mapo not on same grid')
                  goto 1000
                 endif
               endif
c
c--- accumulate stats for map1...use f_stat
              call loopover(map1,env1,i_check,.false.,av_dens2)
                av_dens1=av_temp
                rms1=sigm
c
c--- accumulate stats for map2...use f_stat
              call loopover(map2,env2,o_check,.false.,av_dens2)
                av_dens2=av_temp
                rms2=sigm
c
                 shi_map = av_dens1 - av_dens2
                 sca_map = rms1/rms2
c
                write(6,290) shi_map, sca_map
                 if(record)write(nrec,290) shi_map, sca_map
290              format(' Apply a shift of ',f6.1,' and
     & then scalefactor of ', f6.1 ,' to map2.' )

c---apply shi and sca this time..logical true
c
                 call loopover(map2,env2,o_check,.true.,av_dens2)
                 mscale(map2)=mscale(map1)
                 moffset(map2)=moffset(map1)
        endif
c
1000    return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE set_cut_sol
c     ======================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      integer     map
      character*8 g2chru,lmap
c
      lmap=g2chru('AV_SET_CUT:ENTER_MAP_NAME',1,8,0)
      call av_ass(map,lmap)
      ncut_l=sollev(map)
      ncut_u=sollev(map)
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE set_cutl
c       ===================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
        integer intval
c
        ncut_l=intval('AV_LOWER_CUT:ENTER_CUTOFF',-32000,32000,0)
c
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE set_cutu
c       ===================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
        integer intval
c
        ncut_u=intval('AV_UPPER_CUT:ENTER_CUTOFF',-32000,32000,0)
c
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE set_map
c       ==================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        write(6,10)
        if(record) write(nrec,10)
        return
10      format(' SETM is not yet implemented ')

c       mset=intval('AV_SET_MAP:ENTER_ENV_VALUE',-32000,32000,0)
c
c
c       return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE set_sollev
c     ===============
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      character*8 mapl,G2CHRU
      integer mapn,intval,nopt,keyword
      real    n_solv,realval
c
c----update solvent level information
c
      mapl=G2CHRU('AV_SOLLEV:ENTER_MAP_NAME',1,8,0)
c
      call av_ass(mapn,mapl)
c
      write(6,9)sollev(mapn)
      if(record)write(nrec,9)sollev(mapn)
9     format(' Current solvent level :',i8)
c
      write(6,10)mapn
      if(record)write(nrec,10)mapn
10    format(' Update solvent level information for map:',i4)

600   nopt=keyword('AV_SOLLEV:','NOSC!SCAL!GO  ',0)
c
      if(nopt.eq.1)then
        sollev(mapn)=intval('AV_SOLLEV:ENTER_SOLV_LEVEL',-32000,32000,0)
      endif
c
      if(nopt.eq.2)then
        n_solv=realval('AV_SOLLEV:ENTER_SOLV_LEVEL',-1000000,1000000,0)
        sollev(mapn)=nint(n_solv*mscale(mapn))
      endif
c
      if(nopt.ne.3) goto 600
c
      write(6,8)sollev(mapn)
      if(record)write(nrec,8)sollev(mapn)
8     format(' New solvent level :',i8)
c

c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE sgrp
c     ===============
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      character*8 mapl,G2CHRU,ch_temp
      integer     intval,mapn
c
c----update cell information
c
      ch_temp=G2CHRU('AV_SPCGRP:ENTER_MAP_NAME',1,8,0)
c
      mapl=ch_temp
c
      call av_ass(mapn,mapl)
c
      write(6,10)mapn
      if(record)write(nrec,10)mapn
10    format(' Update spacegroup information for map:',i4)

      lgrp(mapn)=intval('AV_SPCGRP:ENTER_SPACEGRPNUM',1,300,0)

c
      return
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE sig_cut
c       ==================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c-- Local variables
c
        integer     map1,env1
        character*8 imap,ienv
        integer*2   irho
        external    gp3,gp3_a
        integer*2   gp3,gp3_a
c
        character*8 G2CHRU,ch_temp
        real        realval
        integer*2   ss_max,ss_min,abs_irho
        integer     ngadd,nladd,nadd,ngmin,nlmin,nzero,nbig
        integer     iz,jz,iy,jy,ix,jx,jjz,jjy,jjx
c
        logical lass_get,i_check
c
c-- Get input map and envelope
c
        if(.not.lass_get(ch_temp,'MAPI'))then
          ch_temp=G2CHRU('AV_SIG_CUT:INPUT_MAP',1,8,0)
          call lass_set(ch_temp,'MAPI')
        endif
        imap=ch_temp
c
        if(.not.lass_get(ch_temp,'ENVI'))then
          ch_temp=G2CHRU('AV_SIG_CUT:INPUT_ENV',1,8,0)
          call lass_set(ch_temp,'ENVI')
        endif
        ienv=ch_temp
        i_check=ienv.ne.'OFF'
c
c-- Report on map and envelope
c
        if(verbose)then
          write(6,110) imap,ienv
          if(record)write(nrec,110) imap,ienv
        else
          write(6,112) imap,ienv
          if(record)write(nrec,112) imap,ienv
        endif
110     format(/' ABS() AND SIG CUT ELECTRON DENSITY',/,
     &          ' =================================='
     &         /' Input map:       ',a
     &         /' Input envelope:  ',a)
112     format(' ABS & SIG CUT ELECTRON DENSITY',
     &         '  I/P map:',a'  I/P env:',a)
c
        if(verbose)then
          if(i_check) write(6,120)
          if(.not.i_check)write(6,130)
          if(record)then
            if(i_check) write(nrec,120)
            if(.not.i_check)write(nrec,130)
          endif
        endif
120     format(' Input map envelope filtered')
130     format(' Input map not filtered')
c
c-- Get slot numbers for all maps
c-- NB: map 1 should exist by now !
c
        call av_ass(map1,imap)
        if(.not.defined(map1)) then
          write(6,140)
          if(record)write(nrec,140)
          return
        endif
140     format(' %SCUT-ERR: Input map empty')
c
c-- set up allocations for envelopes
c
        if(i_check) call av_ass(env1,ienv)
c
c-- Warning for the fainthearted
c
        if(.not.terse)then
          if(.not.i_check)then
            write(6,*)'**WARNING** input map in peril'
          endif
        endif
c
c-- Get cut thresholds
c
        sig_min=realval('AV_SIG_CUT:MIN NUM_SIG',0.0,1000.0,0)
        sig_max=realval('AV_SIG_CUT:MAX NUM_SIG',sig_min,1000.0,0)
c
        write(6,150)sig_min,sig_max
        if(record)write(nrec,150)sig_min,sig_max
150     format(' If abs value <',f6.1
     &   ,' sig, set to zero. If abs value >',f6.1,' sig, truncate')
c
c-- Set up pointers for the maps that drive the process (backwards!)
c
        call init_map(map1)
        if(i_check)call init_map(env1)
c
c-- Initialize r-factors, corr coefs, etc
c
        ngadd=0
        ngmin=0
        nlmin=0
        nzero=0
        nadd= 0
        nladd=0
        nbig= 0
c
        call f_stat
c
c-- First loop - find sigma
c
c
c-- off to work...
c-- =============
c
c-- Loop over bricks
c
        do iz=1,n_brick(3,map1)
c
          jz=(iz-1)*brick(3,map1)
c
          do iy=1,n_brick(2,map1)
c
            jy=(iy-1)*brick(2,map1)
c
            do ix=1,n_brick(1,map1)
c
              jx=(ix-1)*brick(1,map1)
c
c--- Now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                if(jjy.le.nx(2,map1))then
c
                  do jjx=jx+1,jx+brick(1,map1)
c                      
                  if(jjx.le.nx(1,map1))then
c
c---- Check that this pixel is in required region of the output map
c
                    if(i_check)then
                      if(gp3(env1,jjx,jjy,jjz).ne.0)then
                        irho=gp3(map1,jjx,jjy,jjz)
                        ngadd=ngadd+1
                        call f_a_stat(irho)
                      endif
                    else
                      irho=gp3(map1,jjx,jjy,jjz)
                      ngadd=ngadd+1
                      call f_a_stat(irho)
                    endif
c
                   endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
         enddo
c
c-- Print stats on that pass
c
        call f_p_stat(.true.)
c
c-- Calculate thresholds
c
        ss_min=nint(sigm*sig_min)
        ss_max=nint(sigm*sig_max)
        write(6,160)ss_min,ss_max
        if(record)write(nrec,160)ss_min,ss_max
160     format(' The lower and upper cuts on abs(rho) are:',2i10)
c
c-- Now work through the map taking the absolute density and setting pixels
c-- below the lower threshold to zero and capping pixels above the higher
c-- threshold.
c
        call init_map(map1)
        if(i_check)call init_map(env1)
c
c-- Off to work...
c-- =============
c
c-- Loop over bricks
c
        do iz=1,n_brick(3,map1)
c
          jz=(iz-1)*brick(3,map1)
c
          do iy=1,n_brick(2,map1)
c
            jy=(iy-1)*brick(2,map1)
c
            do ix=1,n_brick(1,map1)
c
              jx=(ix-1)*brick(1,map1)
c
c--- Now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                if(jjy.le.nx(2,map1))then
c
                  do jjx=jx+1,jx+brick(1,map1)
c                      
                  if(jjx.le.nx(1,map1))then
c
c---- Check that this pixel is in required region of the output map
c
                    if(i_check)then
                      if(gp3(env1,jjx,jjy,jjz).eq.0)goto 400
                    endif
c
c---- Get this pixel density
c
                    irho=gp3(map1,jjx,jjy,jjz)
c
c---- Take absolute value
c
                    abs_irho=abs(irho)
c
c---- Deal with pixels below lower threshold
c
                    if(abs_irho.lt.ss_min)then
                      if(irho.ge.0)ngmin=ngmin+1
                      if(irho.lt.0)nlmin=nlmin+1
                      nzero=nzero+1
                      abs_irho=0
c
c---- Deal with pixels above higher threshold
c
                    else if(abs_irho.gt.ss_max) then
                      if(irho.gt.ss_max)nadd=nadd+1
                      if(irho.lt.-ss_max)nladd=nladd+1
                      nbig=nbig+1
                      abs_irho=ss_max
                    endif
c
c---- Deal with remaining pixels
c
                    call pp3(map1,jjx,jjy,jjz,abs_irho)
c
400                 endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
         enddo
c
c-- And we're done
c
c999     continue
        if(verbose)then
          write(6,998)ngadd
          if(record)write(nrec,998)ngadd
        else
          write(6,996)ngadd
          if(record)write(nrec,996)ngadd
        endif
996     format(' SIG CUTTING FINISHED,',i10,' pixels were inspected')
998     format(' SIG CUTTING FINISHED,',/,
     &         ' ====================',/,
     &           i9,' pixels were inspected')
c
        write(6,1000) nzero, ngmin, nlmin
        if(record)write(nrec,1000)nzero,ngmin,nlmin
c
        write(6,1010) nbig, ss_max, nadd, nladd
        if(record)write(nrec,1010) nbig, ss_max,  nadd, nladd
c
1000    format(i10,' pixels set to 0 (of which',i9,
     &  ' were pos and',i9,' neg)')
1010    format(i10,' pixels set to',i8,' (of which',i9,
     &  ' were pos and',i9,' neg)')
c
        if(verbose)then
          if(.not.i_check)then
            write(6,*)'(the whole map was processed)'
            if(record)write(nrec,*)'(the whole map was processed)'
          endif
        endif
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE SJ33EV(A,R,N,MV)
c     ===========================
      IMPLICIT NONE
C
C
        DIMENSION A(1),R(1)
C
        INTEGER MV,IQ,N,J,I,IJ,IA,IND,L,M,MQ,LQ,LM,LL,MM,ILQ,IMQ,IM
        INTEGER IL,ILR,IMR,JQ,K
C
C----   SUBROUTINE EIGEN                                        
C----                                                           
C----      PURPOSE                                              
C----         COMPUTE EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C----          MATRIX                                                 
C----                                                                 
C----       USAGE                                                     
C----          CALL EIGEN(A,R,N,MV)                                   
C----                                                                 
C----       DESCRIPTION OF PARAMETERS                                 
C----          A - ORIGINAL MATRIX (SYMMETRIC), DESTROYED IN COMPUTATION.
C----              RESULTANT EIGENVALUES ARE DEVELOPED IN DIAGONAL OF    
C----              MATRIX A IN DESCENDING ORDER.                         
C---- P.S. THE EQUATION FOR THE INDICES OF A IS A(I,J) -> A(I+(J*J-J)/2))
C----          R - RESULTANT MATRIX OF EIGENVECTORS (STORED COLUMNWISE,  
C----              IN SAME SEQUENCE AS EIGENVALUES)                      
C----          N - ORDER OF MATRICES A AND R                             
C----          MV- INPUT CODE                                            
C----                  0   COMPUTE EIGENVALUES AND EIGENVECTORS          
C----                  1   COMPUTE EIGENVALUES ONLY (R NEED NOT BE       
C----                      DIMENSIONED BUT MUST STILL APPEAR IN CALLING  
C----                      SEQUENCE)                                     
C----                                                                    
C----       REMARKS                                                      
C----          ORIGINAL MATRIX A MUST BE REAL SYMMETRIC (STORAGE MODE=1) 
C----          MATRIX A CANNOT BE IN THE SAME LOCATION AS MATRIX R       
C----                                                                    
C----       SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                
C----          NONE                                                      
C----                                                                    
C----       METHOD                                                       
C----          DIAGONALIZATION METHOD ORIGINATED BY JACOBI AND ADAPTED   
C----          BY VON NEUMANN FOR LARGE COMPUTERS AS FOUND IN 'MATHEMATICAL
C----          METHODS FOR DIGITAL COMPUTERS', EDITED BY A. RALSTON AND    
C----          H.S. WILF, JOHN WILEY AND SONS, NEW YORK, 1962, CHAPTER 7   
C----       IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  
C----       C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      
C----       STATEMENT WHICH FOLLOWS.                                       
C----                                                                      
            DOUBLE PRECISION A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,
     &                 COSX2,SINCS,RANGE_D
C----                                                                      
C----       THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    
C----       APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      
C----       ROUTINE.                                                       
C----                                                                      
C----       THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      
C----       CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTS
C----       40, 68, 75, AND 78 MUST BE CHANGED TO DSQRT.  ABS IN STATEMENT 
C----       62 MUST BE CHANGED TO DABS. THE CONSTANT IN STATEMENT 5 SHOULD 
C----       BE CHANGED TO 1.0D-12.                                         
C----                                                                      
C----       ...............................................................
C----                                                                      
C----       GENERATE IDENTITY MATRIX                                       
C----                                                                      
      RANGE_D=1.0D-12
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C     IF(MV-1) 10,25,10
      IF(MV.EQ.1) GOTO 25
      IQ=-N
      DO 20 J=1,N
      IQ=IQ+N
      DO 20 I=1,N
      IJ=IQ+I
      R(IJ)=0.0
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C     IF(I-J) 20,15,20
      IF(I.NE.J) GOTO 20
      R(IJ)=1.0
   20 CONTINUE
C----                                                           
C----       COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANORMX)  
C----                                                           
   25 ANORM=0.0
      DO 35 I=1,N
      DO 35 J=I,N
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C     IF(I-J) 30,35,30
      IF(I.EQ.J) GOTO 35
      IA=I+(J*J-J)/2
      ANORM=ANORM+A(IA)*A(IA)
   35 CONTINUE
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C     IF(ANORM) 165,165,40
      IF(ANORM.LE.0.0) GOTO 165
C? 40 ANORM=1.414D0*DSQRT(ANORM)                                 
      ANORM=1.41421356*DSQRT(ANORM)
      ANRMX=ANORM*RANGE_D/FLOAT(N)
C----                                                           
C----       INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR    
C----                                                           
      IND=0
      THR=ANORM
   45 THR=THR/FLOAT(N)
   50 L=1
   55 M=L+1
C----                                                           
C----       COMPUTE SIN AND COS                                 
C----                                                           
   60 MQ=(M*M-M)/2
      LQ=(L*L-L)/2
      LM=L+MQ
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C     IF( DABS(A(LM))-THR) 130,65,65
      IF( DABS(A(LM))-THR .LT. 0.0) GOTO 130
      IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5*(A(LL)-A(MM))
      Y=-A(LM)/ DSQRT(A(LM)*A(LM)+X*X)
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C     IF(X) 70,75,75
      IF(X.GE.0.0) GOTO 75
      Y=-Y
   75 SINX=Y/ DSQRT(2.0*(1.0+( DSQRT(1.0-Y*Y))))
      SINX2=SINX*SINX
      COSX= DSQRT(1.0-SINX2)
      COSX2=COSX*COSX
      SINCS =SINX*COSX
C----                                                           
C----       ROTATE L AND M COLUMNS                              
C----                                                           
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=(I*I-I)/2
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C     IF(I-L) 80,115,80
      IF(I.EQ.L) GOTO 115
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C  80 IF(I-M) 85,115,90
      IF(I.EQ.M) GOTO 115
      IF(I.GT.M) GOTO 90
      IM=I+MQ
      GO TO 95
   90 IM=M+IQ
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C  95 IF(I-L) 100,105,105
   95 IF(I.GE.L) GOTO 105
      IL=I+LQ
      GO TO 110
  105 IL=L+IQ
  110 X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C 115 IF(MV-1) 120,125,120
  115 IF(MV.EQ.1) GOTO 125
      ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
  125 CONTINUE
      X=2.0*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C----                                                           
C----       TESTS FOR COMPLETION                                
C----                                                           
C----       TEST FOR M = LAST COLUMN                            
C----                                                           
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C 130 IF(M-N) 135,140,135
  130 IF(M.EQ.N) GOTO 140
      M=M+1
      GO TO 60
C----                                                           
C----       TEST FOR L = SECOND FROM LAST COLUMN                
C----                                                           
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C 140 IF(L-(N-1)) 145,150,145
  140 IF(L.EQ.(N-1)) GOTO 150
      L=L+1
      GO TO 55
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C 150 IF(IND-1) 160,155,160
  150 IF(IND.NE.1) GOTO 160
      IND=0
      GO TO 50
C----                                                           
C----       COMPARE THRESHOLD WITH FINAL NORM                   
C----                                                           
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C 160 IF(THR-ANRMX) 165,165,45
  160 IF(THR-ANRMX .GT. 0.0) GOTO 45
C----                                                           
C----       SORT EIGENVALUES AND EIGENVECTORS                   
C----                                                           
  165 IQ=-N
      DO 185 I=1,N
      IQ=IQ+N
      LL=I+(I*I-I)/2
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=J+(J*J-J)/2
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C     IF(A(LL)-A(MM)) 170,185,185
      IF(A(LL).GE.A(MM)) GOTO 185
      X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
C Obsolete arithmetic IFs reprogrammed, RME 22/4/2010
C     IF(MV-1) 175,185,175
      IF(MV.EQ.1) GOTO 185
      DO 180 K=1,N
      ILR=IQ+K
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
  180 R(IMR)=X
  185 CONTINUE
      RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE spew_map(lun,line,lsize,nmap,speed,bytsize)
c     ==============================================
      IMPLICIT NONE
c
c Purpose: Write formatted section to unit lun from slot nmap
c
c
c Declare common variables and parameters
c
      INCLUDE 'average.fcm'
c     =====================
c
c Declare passed variables
c
      character*4     speed
      integer         lsize,nmap,lun,n_err
      integer*2       line(lsize)
      byte            bytsize(lsize)
c
c Declare internal variables
c
      integer         order(3),istore,i,isec,j,k
c
c Declare external variables
c
      external        gp3
      integer*2       gp3
c
c Initialise counters
c
      istore=0
      n_err=0
c
      do i=1,nx(iuvw(3,nmap),nmap)
        order(iuvw(3,nmap))=i
        isec=nstart(iuvw(3,nmap),nmap)+i-1
        if(speed.eq.'slow')then
          write(lun,6001)isec
        else
          write(lun)isec
        endif
        if(.not.terse)then
          if( i.eq.2.and.i.ne.nx(iuvw(3,nmap),nmap) )then
            write(6,*)'.....................'
            if(record)write(nrec,*)'.....................'
          endif
          if( i.le.1.or.i.eq.nx(iuvw(3,nmap),nmap) )then
            write(6,*)'section',isec,' being written.'
            if(record)write(nrec,*)'section',isec,' being written.'
          endif
        endif
c    
c       write density decentered from zero and NOT scaled down by two
c       (stored as +9999 to -999 formatted)
c
        do j=1,nx(iuvw(2,nmap),nmap)
          order(iuvw(2,nmap))=j
c
          if(speed.eq.'slow')then
            do k=1,lsize
c
              order(iuvw(1,nmap))=k
              line(k) = gp3(nmap,order(1),order(2),order(3))
              line(k) = nint( float(line(k)) + moffset(nmap) )
              if(line(k).gt.9999) then
                line(k)=9999
                n_err=n_err+1
                if(n_err.lt.3)then
                  write(6,*) '%SPEW_MAP-ERR: Pixel value gt 9999'
                  if(record)write(nrec,*)
     &              '%SPEW_MAP-ERR: Pixel value gt 9999'
                endif
              endif
              if(line(k).lt.-999) then
                line(k)=-999
                n_err=n_err+1
                if(n_err.lt.3)then
                  write(6,*) '%SPEW_MAP-ERR: Pixel value lt -999'
                  if(record)write(nrec,*)
     &              '%SPEW_MAP-ERR: Pixel value lt -999'
                endif
              endif
              istore=istore+1
            enddo
c
c If speed = 'fast', convert integer*2 into byte
c
          else
            do k=1,lsize
              order(iuvw(1,nmap))=k
              line(k) = gp3(nmap,order(1),order(2),order(3))
              bytsize(k)=line(k)
              istore=istore+1
            enddo
          endif
c
c Write the line to unit lun
c
          if(speed.eq.'slow')then
            write(lun,6010)(line(k),k=1,lsize)
          else
            write(lun)bytsize
          endif
c
        enddo
c
      enddo
c
      write(6,6020)istore
      if(record)write(nrec,6020)istore
c
      if(n_err.gt.0) then
        write(6,6030)n_err
        if(record)write(nrec,6030)n_err
      endif
c
c Format statements
c
6001  format(/7x,i8/)
6010  format(20i4)
6020  format(' Output completed,',i10,' pixels written')
6030  format(' %SPEW_MAP-ERR:',i10,' pixels outside dynamic range')
c
      return
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE spew_mapb(lun,rho,isize,line,lsize,nmap,scale_flag)
C     ===================================================
      IMPLICIT NONE
c
C   Write formatted section to unit lun from slot nmap
c
c
      INCLUDE 'average.fcm'
c      =====================
c
      integer order(3), istore,lsize,i,nmap,isec,lun,j,k,isize
      integer inrho
      integer*2 line(lsize)
      real rho(isize)
      logical scale_flag
c
      external gp3
      integer*2 gp3
c
c--- Zero counter
c
      istore=0
c
c--- Loop over all pixels on slow axis
c
      do i=1,nx(iuvw(3,nmap),nmap)
c
c--- Set order(slow_axis)=slow_axis position
c
        order(iuvw(3,nmap))=i
c
c--- Eye candy...
c
        isec=nstart(iuvw(3,nmap),nmap)+i-1
        if(.not.terse)then
          if( i.eq.2.and.i.ne.nx(iuvw(3,nmap),nmap) )then
            write(6,*)'.....................'
            if(record)write(nrec,*)'.....................'
          endif
          if( i.le.1.or.i.eq.nx(iuvw(3,nmap),nmap) )then
            write(6,*)'section',isec,' being written.'
            if(record)write(nrec,*)'section',isec,' being written.'
          endif
        endif
c
c--- Reset counter
c
        inrho=1
c
c----  Loop over all pixels in medium axis
c
        do j=1,nx(iuvw(2,nmap),nmap)
c
c---- Set order(medium_axis)=medium_axis position
c
          order(iuvw(2,nmap))=j
c
c----- Loop over all pixels in fast axis
c
          do k=1,nx(iuvw(1,nmap),nmap)
c
c----- Set order(fast_axis)=fast_axis position
c
            order(iuvw(1,nmap))=k
c
c----- Get pixel value - not entirely sure why this needs goes into an array
c
            line(k) = gp3(nmap,order(1),order(2),order(3))
c
c----- Increment counter
c
            istore=istore+1
c
c----- Are we scaling?
c
            if(scale_flag)then
              rho(inrho) = float(line(k))/mscale(nmap)
            else
              rho(inrho) = float(line(k))
            endif
c
c----- Increment marker
c
            inrho=inrho+1
c
c----- End fast loop
c
          enddo
c
c---- End medium loop
c
        enddo
c
c---- Write section to map
c
        call mspew(lun,rho)
c
c--- End slow loop
c
      enddo
c
c--- Print stats on writing
c
      write(6,6020)istore
      if(record)write(nrec,6020)istore
6020  format(' Output completed,',i10,'  pixels written')
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE spin(x,xo,c,v)
C       =========================
      IMPLICIT NONE
C----
C---- spin point x by matrix c centered at position v
C---- output is in xo
C----
c JMD!PORT: j should be int?
c        real x(3),xo(3),c(3,3),v(3),xi(3),j
        real x(3),xo(3),c(3,3),v(3),xi(3)
        integer j
C----
        do j=1,3
          xi(j)=x(j)
        enddo
c
        xo(1)=xi(1)*C(1,1)+xi(2)*C(1,2)+xi(3)*C(1,3)
        xo(2)=xi(1)*C(2,1)+xi(2)*C(2,2)+xi(3)*c(2,3)
        xo(3)=xi(1)*C(3,1)+xi(2)*C(3,2)+xi(3)*C(3,3)
C----
        do j=1,3
          xo(j)=xo(j)+v(j)
        enddo
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      FUNCTION spin_in(iindex,map,i_xrcs)
c     ==================================
      IMPLICIT NONE
C
c  This function returns .true. if the point INDEX can be folded into
c  the available bit of map MAP by crystallographic operator stored in
c  I_XRCS.
c
c
      INCLUDE 'average.fcm'
c       =====================
c
      integer iindex(3),map,t_index(3),i_xrcs,l
      integer ioff,i,i4,ibase
      logical spin_in,s_plus,e_in
c
      spin_in = .false.
c
c---- first convert to true indices
c
        t_index(1) = iindex(1) - norg(1,map)
        t_index(2) = iindex(2) - norg(2,map)
        t_index(3) = iindex(3) - norg(3,map)
c
c---- now apply crystallographic translation/rotation
c
        ioff=( i_xrcs-1 )*12
c
      do i=1,3
        i4=i*4
        ibase=i4-4
        l=xr_insym(ioff+i4,map)
     &     + t_index(1) * xr_insym(ibase+1+ioff,map)
     &     + t_index(2) * xr_insym(ibase+2+ioff,map)
     &     + t_index(3) * xr_insym(ibase+3+ioff,map)
c
c---- now put into cell
c
        l=l+n10(i,map)
        l=mod(l,nunit(i,map))
c
c---- if it fails on this axis, return 
c---- Now, beware when the start and end limits straddle the xtal cell.
c
          s_plus=nstart(i,map).ge.0
          e_in=nend(i,map).le.nunit(i,map)
c
          if(s_plus.or.e_in)then
            if(s_plus.and.e_in)then
              if(l.lt.nstart(i,map))return
              if(l.gt.nend(i,map)  )return
            else
              if(.not.s_plus)then
                if(l.gt.nend(i,map))then
                  l=l-nunit(i,map)
                  if(l.lt.nstart(i,map))return
                endif
              else
                if(l.gt.nunit(1,map))then
                  if(l.lt.nstart(i,map))then
                    l=l+nunit(i,map)
                    if(l.gt.nend(i,map))return
                  endif
                endif
              endif
            endif
          endif
        iindex(i)=l
      enddo
c
c---- ok in-frame, so revert to local indices
c
        iindex(1)=iindex(1)+norg(1,map)
        iindex(2)=iindex(2)+norg(2,map)
        iindex(3)=iindex(3)+norg(3,map)
c
c
      spin_in=.true.
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE spin_inv(x,xo,c,v)
C       =============================
      IMPLICIT NONE
C----
C---- spin point x by inverse matrix c but applying neg position v
C---- output is in xo
C---- used in simple averaging..and labelled envelope stuff
c
c JMD!PORT: j should be int?
c        real x(3),xo(3),c(3,3),v(3),xi(3),j
        real x(3),xo(3),c(3,3),v(3),xi(3)
        integer j
C----
        do j=1,3
          xi(j)=x(j)-v(j)
        enddo
c
        xo(1)=xi(1)*C(1,1)+xi(2)*C(1,2)+xi(3)*C(1,3)
        xo(2)=xi(1)*C(2,1)+xi(2)*C(2,2)+xi(3)*c(2,3)
        xo(3)=xi(1)*C(3,1)+xi(2)*C(3,2)+xi(3)*C(3,3)
C----
      return
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE spin_out(map,iindex,new_index,i_xrcs)
c     ===============================================
      IMPLICIT NONE
C
      INCLUDE 'average.fcm'
c       =====================
c
      integer iindex(3),map,t_index(3),i_xrcs,l,new_index(3)
      integer index1(3),j,ioff,i,i4,ibase
c
c
c---- first convert to true indices
c
      do j=1,3
        new_index(j)=0
        t_index(j) = iindex(j) - norg(j,map)
      enddo
c
c---- now apply crystallographic translation/rotation
c       l=0.0
c
        ioff=( i_xrcs-1 )*12
c
        do i=1,3
        i4=i*4
        ibase=i4-4
        l=xr_insym(ioff+i4,map)

          do j=1,3
            l = l + t_index(j) * xr_insym(ibase+j+ioff,map)
          enddo
c
c---- now put into cell
c
        l=l+n10(i,map)
        l=mod(l,nunit(i,map))
c
c
c
        index1(i)=l
       enddo
c
         do j=1,3
                new_index(j)=index1(j)
         enddo
c
c
      return
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE spin_p(x,xo,c,v,n,inv)
c       =================================
      IMPLICIT NONE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      INCLUDE 'average.fcm'
c       ====================
c
        real x(3,natom), xo(3,natom)
        real c(3,3), v(3)
        integer i,ii
        integer*4 n
        logical inv
c
c
        if(.not.inv)then
c
          do i=1,n
           xo(1,i)=x(1,i)*c(1,1)+x(2,i)*c(1,2)+x(3,i)*c(1,3) +v(1)
           xo(2,i)=x(1,i)*c(2,1)+x(2,i)*c(2,2)+x(3,i)*c(2,3) +v(2)
           xo(3,i)=x(1,i)*c(3,1)+x(2,i)*c(3,2)+x(3,i)*c(3,3) +v(3)
          enddo
c
        else
c
          do i=1,n
             do ii=1,3
                x(ii,i)=x(ii,i) - v(ii)
             enddo
           xo(1,i)=x(1,i)*c(1,1)+x(2,i)*c(1,2)+x(3,i)*c(1,3)
           xo(2,i)=x(1,i)*c(2,1)+x(2,i)*c(2,2)+x(3,i)*c(2,3)
           xo(3,i)=x(1,i)*c(3,1)+x(2,i)*c(3,2)+x(3,i)*c(3,3)
          enddo
c
        endif
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE spinn(x,xo,c,v,rot)
C       ================================
      IMPLICIT NONE
C----
C---- routine to carry out rotation search at the origin using [rot]
C---- output is in xo.  
C---- returns a point.
c
c---- x'=[ops]x+t
c---- x''= [R]x'
c---- 
c
        real x(3),xo(3),c(3,3),v(3),xi(3),rot(3,3)
C----
c
          xi(1)=x(1)*C(1,1)+   x(2)*C(1,2)+    x(3)*C(1,3)+v(1)
          xi(2)=x(1)*C(2,1)+   x(2)*C(2,2)+    x(3)*c(2,3)+v(2)
          xi(3)=x(1)*C(3,1)+   x(2)*C(3,2)+    x(3)*C(3,3)+v(3)
C----
c
          xo(1)=xi(1)*rot(1,1)+xi(2)*rot(1,2)+xi(3)*rot(1,3)
          xo(2)=xi(1)*rot(2,1)+xi(2)*rot(2,2)+xi(3)*rot(2,3)
          xo(3)=xi(1)*rot(3,1)+xi(2)*rot(3,2)+xi(3)*rot(3,3)
c
      return
      end
c*******************************************************************************
c
      SUBROUTINE spinr(x,xo,c,v,cg,rot)
C       =================================
      IMPLICIT NONE
C----
C---- routine to carry out rotation search at the origin using [rot]
C---- output is in xo.  cg is the vector approx of the centre of gravity
c---- of the object to be refined.
C---- returns a point.
c---- x'=[ops]x+t
c---- x''=x'-cg
c---- x'''= [R]x''
c---- x''''=x'''+cg
c
        real x(3),xo(3),c(3,3),v(3),xi(3),rot(3,3),cg(3)
        integer j
C----
c
        xi(1)=x(1)*C(1,1)+x(2)*C(1,2)+x(3)*C(1,3)
        xi(2)=x(1)*C(2,1)+x(2)*C(2,2)+x(3)*c(2,3)
        xi(3)=x(1)*C(3,1)+x(2)*C(3,2)+x(3)*C(3,3)
C----
        do j=1,3
          xi(j)=xi(j)+v(j)-cg(j)
        enddo
c
          xo(1)=xi(1)*rot(1,1)+xi(2)*rot(1,2)+xi(3)*rot(1,3)
          xo(2)=xi(1)*rot(2,1)+xi(2)*rot(2,2)+xi(3)*rot(2,3)
          xo(3)=xi(1)*rot(3,1)+xi(2)*rot(3,2)+xi(3)*rot(3,3)
c
        do j=1,3
          xo(j)=xo(j)+cg(j)
        enddo
      return
      end
c   (
      FUNCTION SPLINE (DX,Y0,Y1,Z0,Z1)
      IMPLICIT NONE
c
        real r6,c1,z0,c2,z1,c3,y0,c4,y1,hx,dx,spline
c
        R6 = 1./6.
        C1 = R6*Z0
        C2 = R6*Z1
        C3 = Y0-C1
        C4 = Y1-C2
        HX = 1.-DX
        SPLINE = (C1*HX**2+C3)*HX + (C2*DX**2+C4)*DX
        RETURN
        END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE av_stat
c       ===============
      IMPLICIT NONE
c
c option for unix use...set status to unknown
c
      INCLUDE 'average.fcm'
c       =====================
c
        integer nopt,keyword

10      nopt=keyword('AV_STAT','NEW !UNKN',0)
c
        if(nopt.eq.1) then
          cstatus='new    '
          goto 20
        endif
        if(nopt.eq.2) then
          cstatus='unknown'
          goto 20
        endif
        if(nopt.ne.1.or.nopt.ne.2) goto 10
c
20      return
        end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE sterefi
c     ==================
      IMPLICIT NONE
C
C Basically cycles over increments in all eulerian theta and trans calcing
C new symm ops etc and looking for improvements in corr coeff.  Options are
C avaliable to refine either translational components or rotational.  New 
C symmetry ops are saved if corr coeff is improved.
C 
C
      INCLUDE 'average.fcm'
c     =====================
C begin
      real rot(3,3),ops_t(3,3,maxsym)
      real v(3),vecs_t(3,maxsym),t(3),vt(3,maxsym)
      real raxis(3),target_r,target_top
      real r_rot(3,3), r_vec(3)
      integer map1,map2,env1,env2,lun,nps,npe,ncycle
      logical better1,better2,rota,transa
      logical i_check,o_check,use_out,no_trans,trans
      real t1,t2,t3,small,dstep,tstep,dstop,tstop,ddelta
      real tdelta,dwhatever,twhatever,realval,del
      integer nstep,nsym1,nsym2,itype,j,nopt,keyword,intval
      integer ll,nj,i,mm,itarget
c
      character*6 par(6)
c
c
      character*80 newfile,charval
      logical fil, r_check,quit
      character*13 mess(0:1)
      data mess/'on all pairs ','from the mean'/
      data par/'rot x ','rot y ','rot z ','transx','transy'
     &          ,'transz'/
c
      target_r=0
      target_top=0
c
C default values
c
c--- smallest worthwhile improvement in corr coef:
      fil=.false.
      small = 0.001
      rota=.true.
      transa=.true.
      r_check=.false.
      NSTEP=10
      DSTEP=4.0
      tstep=1.0
c---  could have tstep for trans part                          
      DSTOP=0.25
      tstop=0.1
      DDELTA=DSTEP
      tdelta=tstep
      itype=11
      istat=1
      nsym1=1
      nsym2=nsym
c
      itarget = 1
c
C--- SET THESE TO ANYTHING....SO 1/2, 1/3 WHATEVER
      DWHATEVER=2.0
      TWHATEVER=2.0
c
c
      call av_setup(map1,map2,env1,env2,nsym1,nsym2,i_check,o_check,
     &                use_out,no_trans,trans,ITYPE,quit)
c
      if(quit)then
        quit=.false.
        return
      endif
c
c----copying sym ops to backups
c
      do j=nsym1,nsym2
        CALL COPY(9,OPS(1,1,j),OPS_back(1,1,j))
        CALL COPY(9,OPS_INV(1,1,j),OPS_INVback(1,1,j))
        CALL COPY(3,VECS(1,j),VECS_back(1,j))
        CALL COPY(9,OPS(1,1,j),save_r(1,1,j))
        CALL COPY(3,VECS(1,j),save_v(1,j))
      enddo
c
c
c
cc
c
20    nopt=keyword('AV_REFI'
     &  ,'NCYC!ROTA!TRAN!NORO!NOTR!GO!WRIT!BACK!TARG',0)
c
c
c
      if(nopt.eq.1)nstep=intval('AV_REFI:NO OF CYCLES',1,100,0)
c
      if(nopt.eq.2)then
        write(6,35)
        if(record)write(nrec,35)
35      format(' Will do rotational search')
        rota=.true.
        dstep=realval('AV_ROTA:INIT_STEP VALUE',0.001,360.0,0)
        dstop=realval('AV_ROTA:FINAL_STEP VALUE',0.001,dstep,0)
        DDELTA=DSTEP
c      
c----could set dwhatever to 2, 3 etc....
      endif
c
      if(nopt.eq.3)then
        write(6,45)
        if(record)write(nrec,45)
45      format(' Will do translational search')
        transa=.true.
        tstep=realval('AV_TRAN:INIT_STEP VALUE',0.001,1000.0,0)
        tstop=realval('AV_TRAN:FINAL_STEP VALUE',0.001,tstep,0)
        tdelta=tstep
      endif
c
      if(nopt.eq.4) then
        write(6,55)
        if(record)write(nrec,55)
55      format(' Don"t perform rotational search')
        rota=.false.
      endif
c
      if(nopt.eq.5) then
        write(6,65)
        if(record)write(nrec,65)
65      format(' Don"t perform translational search')
        transa=.false.
      endif
c
c
      if(nopt.eq.7) then
        write(6,*)' New operators will be written to a file'
        if(record)write(nrec,*)
     &    ' New operators will be written to a file'
c
        newfile=charval('AV_FILE:ENTER_FILENAME',1,80,0)
        fil=.true.
        lun=1
        close(unit=lun)
        open (unit=lun,file=newfile,status=cstatus)
c
      endif
c
      if(nopt.eq.8)return
c
      if(nopt.eq.9)then
        itarget=keyword('AV_REFI:TARGET_FUNCTION','RFAC!JOIN!CORR',0)-2
        if(itarget.eq.-1)then
          write(6,*)' R-factor used as target function'
          if(record)write(nrec,*)
     &      ' R-factor used as target function'
        elseif (itarget.eq.0)then
          write(6,*)' ( CC + 1 - R ) used as target function'
          if(record)write(nrec,*)
     &      ' ( CC + 1 - R ) used as target function'
        elseif (itarget.eq.1)then
          write(6,*)' Correlation coeff used as target function'
          if(record)write(nrec,*)
     &      ' Correlation coeff used as target function'
        endif
c
      endif
c
      if (nopt.ne.6) goto 20
c
c
      if(verbose)then
        write(6,75)
        if(record)write(nrec,75)
75      format(' PERFORM STEEPEST DESCENT REFINEMENT',/,
     &         ' ===================================')
      else
        if(.not.terse)then
          write(6,76)
          if(record)write(nrec,76)
76        format(' PERFORM STEEPEST DESCENT REFINEMENT')
        endif
      endif
      if(rota)then
        write(6,*)'Rotational refinement to be done'
        if(record)write(nrec,*)'Rotational refinement to be done'
      endif
      if(transa)then
        write(6,*)'Translational refinement to be done'
        if(record)write(nrec,*)'Translational refinement to be done'
      endif
c
      write(6,85)nsym1,nsym2
      if(record)write(nrec,85)nsym1,nsym2
85    format(' Will use symmetry operators from ',i4,' to ',i4)
      if( (nsym1.gt.nsym2).or.(nsym2.gt.nsym) )
     &  write(6,*)'%STEREFI-ERR: In choice of symm ops'
      if( ( (nsym1.gt.nsym2).or.(nsym2.gt.nsym) ) .and. record )
     &  write(nrec,*)'%STEREFI-ERR: In choice of symm ops'
c
      if(verbose)then
        do j=nsym1,nsym2
          call prtsym(6,j)
        enddo
      endif
C
      if(itype.ne.11.and.itype.ne.8.and.itype.ne.1.and.itype.ne.64)
     &  itype=11
      write(6,95)itype,mess(istat)
      if(record)write(nrec,95)itype,mess(istat)
95    format(i3,' point interpolation. Correlation coeff calculated ',a)
c
c
      NPS=1
      NPE=6
      if(.not.ROTA)  NPS=4
      IF(.NOT.TRANSA)NPE=3
      IF(NPE.LT.NPS) THEN
        WRITE(6,*)'%STEREFI-ERR: Nothing to refine'
        if(record)WRITE(nrec,*)'%STEREFI-ERR: Nothing to refine'
        RETURN
      endif
c
      if(TRANSA) then
        IF(NO_TRANS)write(6,*)'Trans refinement overrides NOTR'
        IF(NO_TRANS.AND.RECORD)
     &    write(NREC,*)'Trans refinement overrides NOTR'
        TRANS=.true.
        NO_TRANS=.not.TRANS
      endif
c
C
C calculate and output initial correlation
C
      call subaver(map1,map2,env1,env2,nsym1,nsym2,
     &  i_check,o_check,use_out,no_trans,ITYPE,rot(1,1),r_check)
c
c
      call p_stat
c
      if(itarget.eq.-1)target_top =  (1-r)
      if(itarget.eq.1) target_top = cc
      if(itarget.eq.0) target_top =  cc + (1-r)

c
      do j=nsym1,nsym2
c
        if(.not.terse)then
          write(6,105)j
          if(record)write(nrec,105)j
105       format(' Analysis of NCS operator number',i4,
     &    ' before refinement')
          call matrot(ops(1,1,j),t1,t2,t3,rAXIS,'EULE')
          write(6,120)
     &    'Corresponding Eulerian angles (th1,th2,th3)         '
     &    , t1, t2, t3
          if(record)write(nrec,120)
     &    'Corresponding Eulerian angles (th1,th2,th3)         '
     &    , t1, t2, t3
c
          call matrot(ops(1,1,j),t1,t2,t3,rAXIS,'LATT')
          write(6,120)
     &    'Corresponding Lattman angles (th+,th2,th-)          '
     &    , t1, t2, t3
          if(record)write(nrec,120)
     &    'Corresponding Lattman angles (th+,th2,th-)          '
     &    , t1, t2, t3
c
          call matrot(ops(1,1,j),t1,t2,t3,rAXIS,'SPHE')
          write(6,120)
     &    'Corresponding spherical polar angles (psi,phi,kappa)'
     &    , t1, t2, t3
          if(record)write(nrec,120)
     &    'Corresponding spherical polar angles (psi,phi,kappa)'
     &    , t1, t2, t3
c
          call matrot(ops(1,1,j),t1,t2,t3,rAXIS,'AXIS')
          write(6,130)
     &    'Corresponding rotation angle', t3,' about axis', raxis(1)
     &              , raxis(2), raxis(3)
          if(record)write(nrec,130)
     &    'Corresponding rotation angle', t3,' about axis', raxis(1)
     &              , raxis(2), raxis(3)
        endif
      enddo
c
c cycle over increments in angles to try to improve correlation
c do nstep full cycles unless delta already v. small if so skip round 
c
      DO NCYCLE=1,NSTEP
c
c---- if we have run out of steam on the rotational refinement, switch it off.
c     
        IF(DDELTA.LT.DSTOP.and.TRANSA)THEN
          NPS=4
          NPE=6
        ENDIF
c
c---- ditto for translations.
c     
        IF (TDELTA.LT.TSTOP.and.ROTA)THEN
          NPS=1
          NPE=3
        ENDIF
c
c---- is it all over ?
c
        IF(DDELTA.LT.DSTOP.AND.TDELTA.LT.TSTOP) GOTO 300
        IF(DDELTA.LT.DSTOP.AND.(.not.TRANSA))   GOTO 300
        IF(TDELTA.LT.TSTOP.AND.(.not.ROTA))     GOTO 300
c
c---- TRY 6 PARAMETERS USUALLY.... both + and - delta
c
        DO LL=-1,1,2
c
c---- first part applies to rotational steps:
c      
          DO NJ=NPS,NPE
            r_check=.false.
            DO MM=1,3
              T(MM)=0.0
              V(MM)=0.0
            ENDDO
c
            do j=nsym1,nsym2
              call COPY(9,OPS_back(1,1,J),OPS(1,1,J))
              CALL COPY(3,VECS_back(1,J),VECS(1,J))
            enddo
c
c---- set up test angles, Eulerian, theta one to three..
c
            IF (ROTA) THEN
              IF (NJ.eq.2) THEN
                T(1)=-90.0
                T(2)=FLOAT(LL)*DDELTA
                T(3)=90.0
              ENDIF
c               
              IF (NJ.eq.1) T(2)=FLOAT(LL)*DDELTA
              IF (NJ.eq.3) T(3)=FLOAT(LL)*DDELTA
c
              IF (NJ.LT.4) THEN
                r_check=.true.
                do i=1,3
                  do j=1,3
                    rot(i,j)=0.0
                  enddo
                enddo
c
c---- calc rot matrix from eulerian angles
c
                CALL ROTMAT(ROT,T(1),T(2),T(3),rAXIS,'EULE')
c       
c---- need to multiply rot mat with the nonxtal symm mat 
c---- so call matmul to multiply all the nonxtal oper with the rot mat
c---- for delta thetas.
c
              ENDIF
            ENDIF
c
c---- add vector increment to vect of symm ops
c---- loop over all symm ops...from nsym1 to nsym2      
c
            IF (TRANSA) THEN
              IF (NJ.GT.3) V(NJ-3)=FLOAt(LL)*TDELTA
c               
              IF (NJ.GT.3) THEN
                R_CHECK=.FALSE.
                DO J=NSYM1,NSYM2
                  VECS(NJ-3,J)=VECS_BACK(NJ-3,J) + V(NJ-3)
                ENDDO
              ENDIF
            ENDIF
c
c---- call subaver with the rot mat and trans vec and will
c---- hopefully give back corrcoeff....using symops ops(1,1,j) vecs(1,j)
c---- from j = nsym1 to nsym2
c
            call subaver(map1,map2,env1,env2,nsym1,nsym2,
     &  i_check,o_check,use_out,no_trans,ITYPE,rot(1,1),r_check)
c
            call g_stat
c
c
            del=tdelta
            if(nj.lt.4)del=ddelta
            del=del*ll
c
c
            if(itarget.eq.-1)target_r =  (1-r)
            if(itarget.eq.1) target_r = cc
            if(itarget.eq.0) target_r =  cc + (1-r)
c
            IF(target_r.GT.(target_TOP+small)) then
              target_TOP=target_r
c
              better1=.false.
              better2=.false.
c
              write(6,100), 'Cycle',ncycle,par(nj),del,' corr coeff=',
     &                  cc,' r_fac=',r,' r_top=',r1,' **'
              if(record)write(nrec,100)
     &                  , 'Cycle',ncycle,par(nj),del,' corr coeff=',cc,
     &                  ' r_fac=',r,' r_top=',r1,' **'
100           format(1x,a,i3,2x,a,f7.3,1x,a,f6.3,a,f5.3,a,f5.3,a)
c          
              IF (NJ.LT.4) BETTER1=.TRUE.
              IF (NJ.GT.3) BETTER2=.TRUE.
c       
              do j=nsym1,nsym2
                if(better1)then
                  if(cg_check)then
                    vt(1,j)=0.0
                    vt(2,j)=0.0
                    vt(3,j)=0.0
c               
                    call matmul1(rot(1,1),ops(1,1,j),ops_t(1,1,j),
     &                  vecs_back(1,j),cg(1),vt(1,j))
                    CALL COPY(3,Vt(1,J),VECS_T(1,J))
                  else
                    call matmul(rot(1,1),ops(1,1,j),ops_t(1,1,j))
                  endif
                endif
c
                if(better2)CALL COPY(3,VECS(1,J),VECS_T(1,J))
              enddo
c       
            ELSE
c
              if(.not.terse)then
                write(6,100),'Cycle',ncycle,par(nj),del, ' corr coeff=',
     &                  cc,' r_fac=',r,' r_top=',r1,' '
                if(record)write(nrec,100)
     &                  ,'Cycle',ncycle,par(nj),del, ' corr coeff=',cc
     &                  ,' r_fac=',r,' r_top=',r1,' '
              endif
c
            ENDIF
c
          ENDDO
c
c---- THATS END OF THAT PARAMETER
c
        ENDDO
c
c---- thats the end of this cycle, note we intend only to use the best of
c---- the trans OR rot shift, even if we search BOTH
c     
        do j=nsym1,nsym2
          IF(BETTER1) then
            CALL COPY(9,OPS_T(1,1,j),OPS_BACK(1,1,j))
            CALL COPY(3,VECS_T(1,j),VECS_BACK(1,j))
          endif
          IF(BETTER2)CALL COPY(3,VECS_T(1,j),VECS_BACK(1,j))
        enddo
c
        IF( .NOT.(BETTER1.OR.BETTER2) )THEN
c
          DDELTA=DDELTA/DWHATEVER
          TDELTA=TDELTA/TWHATEVER
          if(.not.terse)then
            write(6,110), 'Cycle',ncycle,'  FAILED to improve CC'
     &      ,',  *halve step size(s)*'
            if(record)write(nrec,110)
     &                , 'Cycle',ncycle,'  FAILED to improve CC'
     &      ,',  *halve step size(s)*'
110         format(1x,a,i3,a,a)
          endif
c
        ELSE
c
          if(.not.terse)then
            write(6,110),'Cycle',ncycle,'  *IMPROVED* the CC',
     &        ',  *maintain step size unaltered*'
            if(record)write(nrec,110),'Cycle',ncycle,
     &        '  *IMPROVED* the CC',
     &        ',  *maintain step size unaltered*'
          endif
c         
        ENDIF
c
        better1=.false.
        better2=.false.
c
      ENDDO
c
c---- finish off with stats on the final position
c
300   do j=nsym1,nsym2
       call COPY(9,OPS_back(1,1,j),OPS(1,1,j))
       call COPY(3,VECS_back(1,j),VECS(1,j))
      enddo
c
c
C
      r_check=.false.
      call subaver(map1,map2,env1,env2,nsym1,nsym2,
     &  i_check,o_check,use_out,no_trans,ITYPE,rot(1,1),r_check)
c
C
      CALL P_STAT
c
c----justy print symss to look at ops_back
C
      do j=nsym1,nsym2
        call normal(j,r_rot,r_vec)
        call matmul(r_rot,save_r(1,1,j),ops(1,1,j))
c
c        write (6,3367) j,((r_rot(k,i),i=1,3),k=1,3)
c        write (6,3367) j,((save_r(k,i,j),i=1,3),k=1,3)
c        write (6,3367) j,((ops(k,i,j),i=1,3),k=1,3)
c3367    format (15x,' Symm op n.o. ',I4,3(/,10x,3f10.5)/,10x,3f10.5)
c
        do i=1,3
          vecs(i,j)=vecs(i,j) + r_vec(i)
        enddo
      enddo
c
      if(.not.terse)then
        write(6,*) 'Symmetry operators after refinement'
        if(record)write(nrec,*) 'Symmetry operators after refinement'
      endif
      do j=nsym1,nsym2
        call prtsym(6,j)
        if(fil)call wrtsym(lun,j)
      enddo
      close(unit=lun)
c
      do j=nsym1,nsym2
        call domain(j,save_r(1,1,j),ops(1,1,j))
      enddo
c
c
c
      do j=nsym1,nsym2
        if(.not.terse)then
          write(6,115)j
          if(record)write(nrec,105)j
115       format(' Analysis of updated NCS operator number',i4)
          call matrot(ops(1,1,j),t1,t2,t3,rAXIS,'EULE')
          write(6,120)
     &    'Corresp Eulerian angles (th1,th2,th3)         '
     &    , t1, t2, t3
          if(record)write(nrec,120)
     &    'Corresp Eulerian angles (th1,th2,th3)         '
     &   , t1, t2, t3
c
          call matrot(ops(1,1,j),t1,t2,t3,rAXIS,'LATT')
          write(6,120)
     &    'Corresp Lattman angles (th+,th2,th-)          '
     &    , t1, t2, t3
          if(record)write(nrec,120)
     &    'Corresponding Lattman angles (th+,th2,th-)          '
     &    , t1, t2, t3
c
          call matrot(ops(1,1,j),t1,t2,t3,rAXIS,'SPHE')
          write(6,120)
     &    'Corresp spherical polar angles (psi,phi,kappa)'
     &    , t1, t2, t3
          if(record)write(nrec,120)
     &    'Corresp spherical polar angles (psi,phi,kappa)'
     &    , t1, t2, t3
c
          call matrot(ops(1,1,j),t1,t2,t3,rAXIS,'AXIS')
          write(6,130)
     &    'Corresp rotation angle', t3,' about axis', raxis(1)
     &              , raxis(2), raxis(3)
          if(record)write(nrec,130)
     &    'Corresp rotation angle', t3,' about axis', raxis(1)
     &              , raxis(2), raxis(3)
        endif
      enddo
120   format (x,a,x, 3(f9.4))
130   format (x,a,f9.4,a,x,3(f9.4)/)
c
c
      r_check=.false.
      call subaver(map1,map2,env1,env2,nsym1,nsym2,
     &  i_check,o_check,use_out,no_trans,ITYPE,rot(1,1),r_check)
      CALL P_STAT


      RETURN
      END
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c      SUBROUTINE subaver(map1,map2,env1,env2,nsym1,nsym2,
c     &    i_check,o_check,use_out,no_trans,ITYPE,rot,r_check)
c     ================================================================
c
c Used to be located here...
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE tidy
c     ===============
      IMPLICIT NONE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c---- perform sharpening of envelope
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- local variables:
c
C RME: 23/4/2010: Made index to scratch space explicitly integer*8
      integer*8 j8
C
      integer     map1,env1,new_i(3),nfrac,nsweep,n2,nstep,np_cut,ns_cut
      integer     nsump,nopt,keyword,intval,j,k,ngadd,ns_sum,nns,np,
     &            nsums
      integer     map2,iz,jz,iy,jy,ix,jx,jjz,jjy,jjx,jjz_z,jjy_y,jjx_x
      integer*8   i8temp,offset
      integer*2   irho,nirho
      integer*2   gp3,gp3_a
      external    gp3,gp3_a
      real        psolv
      logical     check,i_check,iiinout,iinout
      character*8 ch_temp
c
c--- initialise variables
c
      check=.false.
      iiinout=.false.
      nfrac=1
      nsweep=1
      n2=2
      nstep=1
      np_cut=4
      ns_cut=4
c
c--- get map numbers and envelope flag
c
      call tidy_setup(map1,env1,i_check)
c
c--- control cards for this routine
c
10    nopt=keyword('AV_TIDY','NSWP!PCUT!SCUT!YSYM!NSYM!GO  !BACK',0)
c
c NSWP
c
      if(nopt.eq.1)nsweep=intval('AV_TIDY:NO_OF_SWEEPS',1,30000,0)
c
c PCUT
c
      if(nopt.eq.2)np_cut =intval('AV_TIDY:NP_CUT',1,6,0)
c
c SCUT
c
      if(nopt.eq.3)ns_cut =intval('AV_TIDY:NS_CUT',1,6,0)
c
c YSYM
c
      if(nopt.eq.4)then
        write(6,900)
        if(record)write(nrec,900)
900     format(' Crystallographic symmetry switched on. Beware!')
        check=.true.
      endif
c
c NSYM
c
      if(nopt.eq.5)then
        write(6,901)
        if(record)write(nrec,901)
901     format(' Crystallographic symmetry switched off.')
        check=.false.
      endif
c
c BACK
c
      if(nopt.eq.7) return
c
c GO
c
      if(nopt.ne.6) goto 10
c
c--- get our dummy map
c
      ch_temp="        "
      call av_ass(map2,ch_temp)
      call copy_head(map1,map2)
      call brickit(map2)
      call alloc(map2)
      offset=ns(map2)-ns(map1)
      do j8=ns(map2),ne(map2)
        scratch(j8)=scratch(j8-offset)
      enddo
c
c---- outer loop - do it twice to mimic earlier behaviour
c
      do j=1,nsweep
        do k=1,2
c
c---- firstly set up pointers for the input map and envelope
c
        call init_map(map1)
        if(i_check)call init_map(env1)
c
c---- initialize r-factors, corr coefs, etc
c
        ngadd=0
        ns_sum=0
        nns=0
        np=0
c
c---- off to work...
c---- =============
c
c---- loop over bricks
c
        do iz=1,n_brick(3,map1)
c
          jz=(iz-1)*brick(3,map1)
c
          do iy=1,n_brick(2,map1)
c
            jy=(iy-1)*brick(2,map1)
c
            do ix=1,n_brick(1,map1)
c
              jx=(ix-1)*brick(1,map1)
c
c---- now loop within each brick
c
              do jjz=jz+1,jz+brick(3,map1)
c
              if(jjz.le.nx(3,map1))then
c
                do jjy=jy+1,jy+brick(2,map1)
c
                if(jjy.le.nx(2,map1))then
c
                  do jjx=jx+1,jx+brick(1,map1)
c               
                  if(jjx.le.nx(1,map1))then
c
c--- envelope check
c
                    if(i_check)then
                      if(gp3(env1,jjx,jjy,jjz).eq.0)goto 200
                    endif
c
c--- stats
c
                    ngadd=ngadd+1
                    nsums=0
                    nsump=0
c
c--- get points adjacent in z
c
                    do jjz_z=jjz-nstep,jjz+nstep,n2
                      iiinout=iinout(map1,jjx,jjy,jjz_z,new_i,check,1)
                      if(.not.iiinout)goto 201
                      irho=gp3_a(map1,new_i)
                      if(irho.eq.0)then
                        nsums=nsums+1
                      else
                        nsump=nsump+1
                      endif
                    enddo
c
c--- get points adjacent in y
c
201                 do jjy_y=jjy-nstep,jjy+nstep,n2
                      iiinout=iinout(map1,jjx,jjy_y,jjz,new_i,check,1)
                      if(.not.iiinout)goto 202
                      irho=gp3_a(map1,new_i)
                      if(irho.eq.0)then
                        nsums=nsums+1
                      else
                        nsump=nsump+1
                      endif
                    enddo
c
c--- get points adjacent in x
c
202                 do jjx_x=jjx-nstep,jjx+nstep,n2
                      iiinout=iinout(map1,jjx_x,jjy,jjz,new_i,check,1)
                      if(.not.iiinout)goto 203
                      irho=gp3_a(map1,new_i)
                      if(irho.eq.0)then
                        nsums=nsums+1
                      else
                        nsump=nsump+1
                      endif
                    enddo
c
c--- get this point
c
203                 nirho=gp3(map1,jjx,jjy,jjz)
c
c--- does it need swapping to protein?
c
                    if((nsump.ge.np_cut).and.(nirho.eq.0))then
                      nirho=1
                      np=np+1
                    endif
c
c--- does it need swapping to solvent?
c
                    if((nsums.ge.ns_cut).and.(nirho.ne.0))then
                      nirho=0
                      nns=nns+1
                    endif
c
c--- stick the new value in the output map
c
                    call pp3(map2,jjx,jjy,jjz,nirho)
c
c--- keep track of solvent content
c
                    if(nirho.eq.0) ns_sum=ns_sum+1
c
c--- and round we go again
c
200               continue
                  endif
                  enddo
                endif
                enddo
              endif
              enddo
            enddo
          enddo
        enddo
c
c
c--- print stats on the map after this pass
c
        psolv=(float(ns_sum)/float(ngadd))*100.0
        write(6,205)ngadd,np,nns,psolv
        if(record)then
          write(nrec,205) ngadd, np, nns, psolv
        endif
205     format(' Of',i10,' pixels,',i8,' set to protein,'
     &        ,i8,' to solvent,',f5.1,'%solv')
c
c--- copy new map over original map
c
        i8temp=ns(map1)
        ns(map1)=ns(map2)
        ns(map2)=i8temp
c
        i8temp=ne(map1)
        ne(map1)=ne(map2)
        ne(map2)=i8temp
c
c--- and second pass
c
        enddo
c
c--- and next step
c
      enddo
c
c--- delete temporary map - use parser deallocation routine
c
      if(.not.P2DAL(scratch(ns(map2)))) then
        write(6,*)' %DEL_MAP-ERROR: memory deallocation fails'
        if(record)write(nrec,*)
     &    ' %DEL_MAP-ERROR: memory deallocation fails'
      endif
c
c--- set first free pointer
c
      nmaps=nmaps-1
c
c--- and its all over
c
      return
      end
C
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE tidy_setup(map1,env1,i_check)
c       ========================================
      IMPLICIT NONE
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      INCLUDE 'average.fcm'
c       =====================
c
c--- local variables
c
        logical     i_check,lass_get
        integer     map1,env1
        character*8 g2chru,ch_temp,imap,ienv
c
c--- set up input map if not already assigned
c
        if(.not.lass_get(ch_temp,'MAPI'))then
          ch_temp=G2CHRU('AV_TIDY:INPUT_MAP',1,8,0)
          call lass_set(ch_temp,'MAPI')
        endif
        imap=ch_temp
c
c--- set up input envelope if not already assigned
c
        if(.not.lass_get(ch_temp,'ENVI'))then
          ch_temp=G2CHRU('AV_TIDY:INPUT_ENV',1,8,0)
          call lass_set(ch_temp,'ENVI')
        endif
        ienv=ch_temp
c
c--- report what we've got
c
        if(verbose)then
          write(6,100)imap,ienv
          if(record)write(nrec,100)imap,ienv
        else
          write(6,110)imap,ienv
          if(record)write(nrec,110)imap,ienv
        endif
100     format(/' TIDY ENVELOPE',/,
     &          ' =============',/,
     &          ' Input map:       ',a,/,
     &          ' Input envelope:  ',a)
110     format(' TIDY ENVELOPE I/P map:',a,'  I/P env:',a)
c
c--- turn off input envelope checking if its called OFF
c
        i_check=ienv.ne.'OFF'
c
c--- report if using input envelope
c
        if(verbose)then
          if(i_check)write(6,120)
          if(.not.i_check)write(6,130)
          if(record)then
            if(i_check) write(nrec,120)
            if(.not.i_check)write(nrec,130)
          endif
        endif
120     format(' Input map envelope filtered')
130     format(' Input map not filtered')
c
c--- get slot number for input map
c
        call av_ass(map1,imap)
c
c--- map 1 should exist by now !
c
        if(.not.defined(map1)) then
          write(6,140)
          if(record)write(nrec,140)
140       format(' %TIDY_SETUP-ERR: Input map empty')
          return
        endif
c
c--- get slot number for input envelope
c
        if(i_check)call av_ass(env1,ienv)
c
c--- and thats it
c
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE time_get
c       ===================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c
        character*20 date
        character*18 time
        integer      ilen
c
c
        call P2TIME(time,ilen)
        call P2DATE(date,ilen)
        write(6,10),time,date
        if(record)write(nrec,10),time,date
10      format(' TIME: it is now ',a,' on ',a)
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE title
c       ================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
        character*72 i_title
        character*72 charval
c
        i_title=charval('AV_TITLE',1,72,0)
        write(6,20)i_title
        if(record)write(nrec,20)i_title
20      format(' ',a72,/,' ',72('='))
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE tosser
C     =================
      IMPLICIT NONE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                                      
C       TOSSER ...  PRODUCE MATRICES FOR NCS AVERAGING   
c               
C      Altered from eyj & dis toss_bits that checked 3-fold matrices 
c      in retrospect from .pdb coords.
C      Now the several groups of atoms related by ncs are read in all
C      in one pdb file......with each atom flagged to a particular group
c      by a number in the "10's" column of residue number. 
c       ie. for trimer 1, res nos go from 11 to 13....and similarly
c      for trimer 2, res nos go from 21 to 23.
C                                                                          
C       june-92  dis/jmg                         
C   ********************************************************************
C
C
      INCLUDE 'average.fcm'
c     =====================
C
      integer      lun,maxres,nstep,mincop,maxcop,icon,nnres
      integer      nc,jstep,maxr
      real         rinv(3,3),tinv(3),det,rtoss,ttoss,rmse,wt_cent
      REAL         X(3),B,WT
      CHARACTER*4  ACHECK,start,atom,big,startup
      CHARACTER*1  ALAB(8),line(80)
      character*80 long_line
      character*6  chain,chain_old
      character*4  atom_type
      equivalence (line,start)
      equivalence (line,long_line)
      common /scrtch/rtoss(3,3),ttoss(3)
c
      logical      entry_l,tossed,cen
c
      character*80 charval
      real         centre(3)
      integer      nread,maxres1,ncop,num,jmax,kov,i,j,j1,j2
c      equivalence  (i_temp,ch_temp)
      equivalence  (alab,atom_type)
      integer      k
      character*80 fil
      CHARACTER*4  toss,cent,atus,step,atuse
      data atom/'ATOM'/,big/'BIG '/
      data toss/'TOSS'/,cent/'CENT'/,atus/'ATUS'/,step/'STEP'/
c

C========================================================================
C
c----lun is input stream for coords
      lun=1
      maxr=0
      nc=0
      maxres1=0
      jstep=100
      chain_old='@@@@@@'
      nstep=10
      mincop=9999
      maxcop=-9999
      tossed=.false.
c
c
c
c ------read in pdb format for sub1
c
c
      write(6,20)
      if(record)write(nrec,20)
20    format(' AV_TOSS> ** Evaluate NCS matrices from a PDB file **')
      fil=charval('AV_READ_PDB:ENTER_FILE',1,80,0)
      if(.not.terse)then
        write(6,40)fil(1:60), lun
        if(record)write(nrec,40)fil(1:60), lun
40      format(' AV_TOSS> Pdbfile:', a60,/,'          opened on unit=',
     &         i3)
      else
        write(6,41)fil(1:60)
        if(record)write(nrec,41)fil(1:60)
41      format(' AV_TOSS> Pdbfile:', a60)
      endif
c
      close (unit=lun)
c JMD!PORT: open  (unit=lun,file=file,status='old',readonly)     
      open  (unit=lun,file=fil,status='old')
c
c
100   READ (lun,21,END=99)line
21    FORMAT(80a1)
c251   format(' ',80a1)
      icon=4
c
c---- now the header records can contain keyworded commands - devious
c
c----convert to upper case (this is a routine in Robert Esnouf's GETWORD)
c
      startup=start
      call P2UCAS(startup)
c----BIG
      if(startup.eq.big)then
        write(6,252)
        if(record)write(nrec,252)
252     format(' AV_TOSS> BIG selected: up to 99 copies allowed')
        nstep=100
      endif
c----TOSS
      if(startup.eq.toss)then
        write(6,262)
        if(record)write(nrec,262)
262     format(' AV_TOSS> Coordinate file uses GAP/XPLOR chain ids')
        tossed=.true.
      endif
c----CENT
      if(startup.eq.cent)then
        wt_cent=10.0
        read(long_line(5:80),*)centre,wt_cent
        cen=.true.
        write(6,272)CENTRE,wt_cent
        if(record)write(nrec,272)CENTRE,wt_cent
272     format(' AV_TOSS> CENTre of object will be restrained',/,
     & '          to be at:',3f9.2,', with weight:',f8.1)
      endif
c----ATUS
      if(startup.eq.atus)then
        read(long_line(5:80),273)atuse
273     format(a)
        write(6,282)atuse
        if(record)write(nrec,282)atuse
282     format(' AV_TOSS> Pick out one atom per residue, label:',a)
      endif
c----STEP
      if(startup.eq.step)then
      read(long_line(5:80),*)jstep
        write(6,292)jstep
        if(record)write(nrec,292)jstep
292     format(' AV_TOSS> Use only every',i5,'th residue')
      endif
c----
      if(startup.ne.atom)goto 100
c
c-------decode line
c
      read(long_line,55)acheck,alab,nnres,x,wt,b,chain
55    format(a4,8x,4a1,1x,4a1,1x,i4,4x,3f8.3,2f6.2,6x,a6)
      nread=nread+1
C
      if(.not.tossed)MAXRES=nNRES/nstep
      if(.not.tossed)NCOP=nNRES-(nstep*MAXRES)
C
      if(tossed) then
        if(chain.ne.chain_old)then
          chain_old=chain
          nc=nc+1
          if(record)write(NREC,62)chain,nc,maxr
          write(6,62)chain,nc,maxr
62        format(' NEW CHAIN:',a,'is number',i5,
     &          ' (prev. chain used',i5,' guide points)')
          maxr=0
        endif
C        write(6,622)atom_type,atuse
C622     format('ATOMTYPE in FILE:',a,'###ATOMTYPE TARGET:',a,'###')
        if(atom_type.ne.atuse) goto 100
        if(nnres/jstep*jstep.ne.nnres)goto 100
        maxr=maxr+1
        maxres=maxr
        ncop=nc
      endif
c
      maxres1=maxres1+1
      if(maxres1.gt.natom)then
        write(6,*)'%TOSSER-ERR: NATOM too small-change parameter'
        if(record)write(nrec,*)
     &            '%TOSSER-ERR: NATOM too small-change parameter'
        return
      endif
c
      if(ncop.le.mincop)mincop=ncop
      if(ncop.gt.maxcop)maxcop=ncop
      if(.not.terse)then
        if(nread.le.2)then
          WRITE(6,65)MAXRES,NCOP,X
          IF(RECORD)WRITE(NREC,65)MAXRES,NCOP,X
65        FORMAT(' ATOM READ, ATOM #:',I4,'  COPY',I4,' POS:',3F8.2)
        endif
      endif
c
      do num=1,3
        coords(num,maxres,ncop)=x(num)
      enddo
      goto100
99    write(6,919) NREAD,MINCOP,MAXCOP
      if(record)write(NREC,919) NREAD,MINCOP,MAXCOP
919   format(' Number of guide points scanned:',i6,
     &       '  Min NCS op:',i4,'  Max NCS op:',i4)
      close(unit=lun)
C
c
c7777  continue
c
c----coordinates are assumed to be orthogonal and in A.
c
      jmax=1
      if(proper)jmax=MAXCOP
      do kov=0, MAXCOP-1
        call rigset
c----firstly restrain centre if required
        if(cen) call rigadd(centre,centre,10.0)
        do i=1,maxres
          do j=1,jmax
c
            j1=j
            j2=j+kov
            if(j2.gt.MAXCOP) j2=j2-MAXCOP
c
c----unit weights for least squares
c
            call rigadd(coords(1,i,j2), coords(1,i,j1),1.0)
c
          enddo
        enddo
c
        entry_l=.true.
        call rigsol(rtoss,ttoss,rmse,entry_l)
        write(6,32) 1+kov
        if(record)write(nrec,32) 1+kov
32      format(' Symm element 1 onto symm element',i4,' matrix')
        write(6,320)rmse
        if(record)write(nrec,320)rmse
320     format(' RMS error from toss:',f10.2)
c
        if(verbose)then
          WRITE(6,31)((RTOSS(I,J),j=1,3),i=1,3),TTOSS
          if(record)WRITE(nrec,31)((RTOSS(I,J),j=1,3),i=1,3),TTOSS
31        FORMAT(' Symmetry operators for input to GAP ',
     &       /,3(T20,3F10.5,/),(T20,3F10.5))
        endif
c

        nsym=nsym+1
c
        call invers(rtoss,rinv,det)
c
        if(abs(abs(det)-1.0).gt.0.99)then
          write(6,*)'%TOSSER-ERR: Determ. .ne.1'
          if(record)write(nrec,*)'%TOSSER-ERR: Determ. .ne.1'
        endif
c
        do j=1,3
          tinv(j)=-ttoss(j)
        enddo
        call copy(9,rtoss,ops(1,1,nsym))
        call copy(9,rinv,ops_inv(1,1,nsym))
        call copy(3,ttoss,vecs(1,nsym))
c
        write(6,101)nsym,((rtoss(k,j),j=1,3),k=1,3), (ttoss(j),j=1,3)
     &                  ,  ((rinv(k,j),j=1,3),k=1,3), (tinv(j),j=1,3)
        if(record)write(nrec,101)nsym,
     &     ((rtoss(k,j),j=1,3),k=1,3), (ttoss(j),j=1,3)
     &  ,  ((rinv(k,j),j=1,3),k=1,3), (tinv(j),j=1,3)
101     format(' 'i6,' Symm ops now in use',/,2(3(/,10x,3f10.5)/,10x,
     &   3f10.5,/))
c
      enddo
C
      return
      end
c
c*******************************************************************************
c
      SUBROUTINE VADD(V1,V2,V3)
c       =========================
      IMPLICIT NONE
c*******************************************************************************
c
C Subroutine to add vectors V2 and V3, storing the result in V1
c
        REAL V1(3),V2(3),V3(3)
c
        V1(1)=V2(1)+V3(1)
        V1(2)=V2(2)+V3(2)
        V1(3)=V2(3)+V3(3)
c
        RETURN
        END

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      FUNCTION VDOT(V1,V2)
c       ====================
      IMPLICIT NONE
c*******************************************************************************
c
C Function to calculate the dot product of vectors V1 and V2
c
        REAL V1(3),V2(3),vdot
c
c       write(6,*) ' in sub vdot'
        VDOT=V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
c
        RETURN
        END
c
c*******************************************************************************
c
      FUNCTION VLEN(V1)
c       =================
      IMPLICIT NONE
c*******************************************************************************
c
C Function to return the length of vector V1
c
        REAL V1(3),vlen
c
        VLEN=SQRT(V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))
c
        RETURN
        END

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE vmtply(c,x,xo)
C       =========================
      IMPLICIT NONE
C----
C---- MULTIPLIES VECTOR X by MATRIX C, result to xo
C----
        real x(3),xo(3),c(3,3)
C----
        xo(1)=x(1)*C(1,1)+x(2)*C(1,2)+x(3)*C(1,3)
        xo(2)=x(1)*C(2,1)+x(2)*C(2,2)+x(3)*c(2,3)
        xo(3)=x(1)*C(3,1)+x(2)*C(3,2)+x(3)*C(3,3)
C----

      return
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE VSUB(V1,V2,V3)
c       =========================
      IMPLICIT NONE
c*******************************************************************************
c
C Subroutine to subtract vectors V2 and V3, storing the result in V1
c
        REAL V1(3),V2(3),V3(3)
c
        V1(1)=V2(1)-V3(1)
        V1(2)=V2(2)-V3(2)
        V1(3)=V2(3)-V3(3)
c
        RETURN
        END
c*******************************************************************************
c
      SUBROUTINE VUNIT(V1,V2)
c       =======================
      IMPLICIT NONE
c*******************************************************************************
c
C Subroutine to calculate the unit vector in the direction of vector V2 and
C store the result in V1
c
        REAL V1(3),V2(3),vlen,a
        external vlen
c
        A=VLEN(V2)
        V1(1)=V2(1)/A
        V1(2)=V2(2)/A
        V1(3)=V2(3)/A
c
        RETURN
        END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE write_headers(lun,title,iuvw,mxyz,nstart,nend
     &       ,cell,lspgrp,rholim,scale_r,plus,lmode,terse,verbose,speed)
C     ================================================================
      IMPLICIT NONE
c
c
c*****************************************************************************
      common/av_rec /record,nrec
c*****************************************************************************
      character*4 speed
      logical record,terse,verbose
      integer nrec
c
c
c
c
c  Formatted, write keyworded header
c
      integer      lun
      character*80 title
      integer      iuvw(3),mxyz(3),nstart(3),nend(3)
     &            ,lspgrp,lmode
      real         cell(6),scale_r,plus,rholim(4)
c
      integer      i,limits(2,3)
      character*1  xyz(3)
      data xyz/'X','Y','Z'/
c
c     setup limits
c
      do i=1,3
        limits(1,i)=nstart(i)
        limits(2,i)=nend(i)
      enddo
c
      if(verbose)then
        write(6,5000)
        if(record)write(nrec,5000)
 5000   format(' HEADER INFORMATION AS WRITTEN TO OUTPUT MAP',/
     &          ' ============================================')
      else
        if(.not.terse)then
          write(6,5010)
          if(record)write(nrec,5010)
 5010     format(' HEADER INFORMATION AS WRITTEN TO OUTPUT MAP')
        endif
      endif
c
      if(speed.eq.'slow')then
c
        write(lun,6000)
 6000   format('MAPEXCHANGE HEADER')
        write(lun,6010) title
 6010   format(/A)
        write(lun,6020) xyz(iuvw(1)),xyz(iuvw(2)),xyz(iuvw(3))
 6020   format('AXIS    ',3(7X,A1))
        write(lun,6030) mxyz
 6030   format('GRID    ',3I8)
        write(lun,6040) limits
 6040   format('XYZLIM  ',6I8)
        write(lun,6050) lspgrp
 6050   format('SPACEGROUP',6X,I8)
        write(lun,6060) lmode
 6060   format('MODE    ',I8)
        write(lun,6070) cell
 6070   format('CELL    ',6F10.3)
        write(lun,6080) rholim
 6080   format('RHOLIM  ',4G16.6)
        write(lun,6090) scale_r,plus
 6090   format('PACK    ',2G16.6)
        write(lun,6095)
 6095   format('END HEADER')
c
      else
c
        write(lun) title,
     &    xyz(iuvw(1)),xyz(iuvw(2)),xyz(iuvw(3)),
     &    mxyz, limits, lspgrp, lmode, cell, rholim, scale_r, plus
c
      endif
c
      if(.not.terse)then
        write(6,6100)
        if(record)write(nrec,6100)
 6100   format(' MAPEXCHANGE HEADER')
        write(6,6110) title
        if(record)write(nrec,6110) title
 6110   format(' TITLE'/1X,A)
        write(6,6120) xyz(iuvw(1)),xyz(iuvw(2)),xyz(iuvw(3))
        if(record)write(nrec,6120) xyz(iuvw(1)),xyz(iuvw(2)),
     &    xyz(iuvw(3))
 6120   format(' AXIS    ',3(7X,A1))
        write(6,6130) mxyz
        if(record)write(nrec,6130) mxyz
 6130   format(' GRID    ',3I8)
        write(6,6140) limits
        if(record)write(nrec,6140) limits
 6140   format(' XYZLIM  ',6I8)
        write(6,6150) lspgrp
        if(record)write(nrec,6150) lspgrp
 6150   format(' SPACEGROUP',6X,I8)
        write(6,6160) lmode
        if(record)write(nrec,6160) lmode
 6160   format(' MODE    ',I8)
        write(6,6170) cell
        if(record)write(nrec,6170) cell
 6170   format(' CELL    ',6F10.3)
        write(6,6180) rholim
        if(record)write(nrec,6180) rholim
 6180   format(' RHOLIM  ',4G16.6)
        write(6,6190) scale_r,plus
        if(record)write(nrec,6190) scale_r,plus
 6190   format(' PACK    ',2G16.6)
        write(6,6195)
        if(record)write(nrec,6195)
 6195   format(' END HEADER')
        endif
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE write_map(speed)
c     ===========================
      IMPLICIT NONE
c
c Purpose: Write out an ascii map (speed='slow') or an envelope (speed='fast')
c
c Declare common parameters and variables
c
      INCLUDE 'average.fcm'
c     =====================
c
c Declare passed variables
c
      character*4     speed
c
c Declare local variables
c
      character*80    charval,fil
      character*8     G2CHRU,ch_temp,label
      integer         lun,nmap,nline
      integer*8       nxtbyte
      logical         lass_get
      byte            abyte(2)
c
c Say we got here
c
      write(6,10)speed
      if(record)write(nrec,10)speed
c
c Do we know which map?
c
      if(.not.lass_get(ch_temp,'MAPI'))then
        ch_temp=G2CHRU('AV_WRITE_MAP:ENTER_LABEL',1,8,0)
        write (6,20) ch_temp
        if(record)write (nrec,20) ch_temp
        call lass_set(ch_temp,'MAPI')
      endif
      label=ch_temp
c
c Get file name and open file on unit 1
c
      fil=charval('AV_WRITE_MAP:FILE',1,80,0)
      lun=1
      close (unit=lun)
c JMD!PORT: gfortran f77 open has no carriagecontrol attribute
c     if(speed.eq.'slow')open(unit=lun,file=file,status=cstatus,
c    &  carriagecontrol='list')     
      if(speed.eq.'slow')open(unit=lun,file=fil,status=cstatus)
      if(speed.eq.'fast')open(unit=lun,file=fil,status=cstatus,
     &  form='unformatted')
      if(.not.terse)then
        write(6,40)fil(1:60), lun
        if(record)write(nrec,40)fil(1:60), lun
      else
        write(6,41)fil(1:60)
        if(record)write(nrec,41)fil(1:60)
      endif
c
c Turn the map name into a number
c
      call av_ass(nmap, label)
      if(.not.defined(nmap))then
        write(6,*)'%WRITE_MAP-ERR: Map not defined abort write'
        if(record)write(nrec,*)
     &    '%WRITE_MAP-ERR: Map not defined abort write'
        goto 999
      endif
c
c Write out the map header information
c
      call write_headers(lun,titles(nmap),iuvw(1,nmap),nunit(1,nmap),
     &  nstart(1,nmap),nend(1,nmap)
     &  ,XRcell(1,nmap),lgrp(nmap),
     &  rholim(1,nmap),mscale(nmap),moffset(nmap),mtype(nmap)
     &  ,terse,verbose,speed)
c
c How much scratch space do we need?
c
      nline=nx(iuvw(1,nmap),nmap)
c
c Allocate some integer*2 scratch space
c 
      nxtpix=P2AL(scratch(1),scratch(2),nline)
      if(nxtpix.eq.0) then
        write(6,*)'%WRITE_MAP-ERR: Memory allocate fails'
        if(record)write(nrec,*)'%WRITE_MAP-ERR: Memory allocate fails'
        goto 999
      endif
c
c Allocate some byte scratch space
c 
      nxtbyte=P2AL(abyte(1),abyte(2),nline)
      if(nxtbyte.eq.0) then
        write(6,*)'%WRITE_MAP-ERR: Memory allocate fails'
        if(record)write(nrec,*)'%WRITE_MAP-ERR: Memory allocate fails'
        goto 999
      endif
c
c Call actual writing routine
c
      call spew_map(lun,scratch(nxtpix),nline,nmap,speed,abyte(nxtbyte))
c
c Finished writing so close lun
c
      close(unit=lun)
c
c Deallocate scratch space used
c
      if(.not.P2DAL(scratch(nxtpix)))then
        write(6,*)'%WRITE_MAP-ERR: Memory deallocate fails'
        if(record)write(nrec,*)'%WRITE_MAP-ERR: Memory deallocate fails'
        goto 999
      endif
      if(.not.P2DAL(abyte(nxtbyte)))then
        write(6,*)'%WRITE_MAP-ERR: Memory deallocate fails'
        if(record)write(nrec,*)'%WRITE_MAP-ERR: Memory deallocate fails'
        goto 999
      endif
c
c Done!
c
      return
c
c On error...
c
999   write(6,*)'%WRITE_MAP-ERR: **ERROR: abandon write**'
      if(record)write(nrec,*)'%WRITE_MAP-ERR: **ERROR: abandon write**'
c
c Format statements
c
10    format(' WRITE MAP: ',a4)
20    format(' label for map to be written = ', a)
40    format(' Mapfile:', a60,/,' opened on unit=',i3)
41    format(' Mapfile:', a60)
c
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE write_mapb(scale_flag)
c     =====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
c----- write out a binary CCP4 map
c
      character*80 charval
      character*8  G2CHRU
      character*8  ch_temp
      real         stub(2)
      integer      nmap,nline,limits(2,3),nsize
      integer*8    np1,npoint
      character*8  label
      character*80 fil
      logical      scale_flag,lass_get
c
      integer      lun,i,nu1,nu2,nv1,nv2,nw1,nsec
      character*80 symop_lib
c
      if(record)write(nrec,10)
10    format(' WRITE MAP: CCP4 BINARY')
      if(.not.lass_get(ch_temp,'MAPI'))then
c
        ch_temp=G2CHRU('AV_WRITE_MAP:ENTER_LABEL',1,8,0)
        write (6,20) ch_temp
        if(record)write (nrec,20) ch_temp
20      format(' label for map to be written = ', a)
        call lass_set(ch_temp,'MAPI')
      endif
      label=ch_temp
c
      fil=charval('AV_WRITE_MAP:FILE',1,80,0)
      lun=1
      close (unit=lun)
      if(.not.terse)then
        write(6,40)fil(1:60), lun
        if(record)write(nrec,40)fil(1:60), lun
40      format(' Mapfile:', a60,/,' opened on unit=',i3)
      else
        write(6,41)fil(1:60)
        if(record)write(nrec,41)fil(1:60)
41      format(' Mapfile:', a60)
      endif
c
      call av_ass(nmap, label)
      if(.not.defined(nmap))then
        write(6,*)'%WRITE_MAP-ERR: Map not defined abort write'
        if(record)write(nrec,*)
     &            '%WRITE_MAP-ERR: Map not defined abort write'
        goto 999
      endif
c
c--- set limits on xyz
c
      do i=1,3
        limits(1,i) = nstart(i,nmap)
        limits(2,i) = nend(i,nmap)
      enddo
      nu1 = limits(1,iuvw(1,nmap))
      nu2 = limits(2,iuvw(1,nmap))
      nv1 = limits(1,iuvw(2,nmap))
      nv2 = limits(2,iuvw(2,nmap))
      nw1 = limits(1,iuvw(3,nmap))
      nsec=limits(2,iuvw(3,nmap))-limits(1,iuvw(3,nmap))+1
c
      call mwrhdl(lun,fil,titles(nmap),nsec,iuvw(1,nmap),nunit(1,nmap),
     &  nw1,nu1,nu2,nv1,nv2,Xrcell(1,nmap),lgrp(nmap),mtype(nmap))
c
c Get path to symop.lib from the environment variable CLIBD
c
      call ugtenv('CLIBD',symop_lib)
      symop_lib=symop_lib(1:index(symop_lib,' ')-1)//'/symop.lib'
c
c Write symops to unit lun
c
c JMD!PORT: open(4,file=symop_lib,status='OLD',readonly)
      open(4,file=symop_lib,status='OLD')
      call msyput(4,lgrp(nmap),lun)
      close(unit=4)
c
      nline=nx(iuvw(1,nmap),nmap)
      nsize= (nu2-nu1+1)*(nv2-nv1+1)
      npoint=P2AL(stub(1),stub(2),nsize)
      np1=P2AL(scratch(1),scratch(2),nline)
c---- check there is enough scratch space 
      if(npoint.eq.0) then
        write(6,*)'%WRITE_MAP-ERR: Memory allocate fails'
        if(record)
     &    write(nrec,*)'%WRITE_MAP-ERR: Memory allocate fails'
        goto 999
      endif
      if(np1.eq.0) then
        write(6,*)'%WRITE_MAP-ERR: Memory allocate fails'
        if(record)
     &    write(nrec,*)'%WRITE_MAP-ERR: Memory allocate fails'
        goto 999
      endif
c
      call spew_mapb(lun,stub(npoint),nsize,scratch(np1),nline,nmap,
     &     scale_flag)
c
      call mwclose(lun)
c
      if(.not.P2DAL(stub(npoint)))then
        write(6,*)'%WRITE_MAP-ERR: Memory deallocate fails'
        if(record)
     &    write(nrec,*)'%WRITE_MAP-ERR: Memory deallocate fails'
        goto 999
      endif
      if(.not.P2DAL(scratch(np1)))then
        write(6,*)'%WRITE_MAP-ERR: Memory deallocate fails'
        if(record)
     &    write(nrec,*)'%WRITE_MAP-ERR: Memory deallocate fails'
        goto 999
      endif
      return
c
999   write(6,*)'%WRITE_MAP-ERR: **ERROR: abandon write**'
      if(record)write(nrec,*)'%WRITE_MAP-ERR: **ERROR: abandon write**'
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE write_mean
c     ====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      character*80    newfile,charval
      character*8     G2CHRU,ch_temp
      integer         lun,map1
      logical         lass_get
c
      if(.not.lass_get(ch_temp,'MAPI'))then
        ch_temp=G2CHRU('AV_SIG_CUT:INPUT_MAP',1,8,0)
        call lass_set(ch_temp,'MAPI')
      endif
      call av_ass(map1,ch_temp)
c
      newfile=charval('AV_WRIT_MEAN:ENTER_FILENAME',1,80,0)
      lun=1
      close(unit=lun)
      open (unit=lun,file=newfile,status=cstatus)
c
      write(lun,10)rholim(3,map1)
10    format(f21.10)
c
      close(unit=lun)
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE write_sym
c     ====================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      character*80 newfile,charval
      integer intval,lun,nsym1,nsym2,j
c
      if (nsym.lt.1) then
        write(6,10)
        if(record) write(nrec,10)
        return
      endif
10    format(' No Symmetry Operators defined')
c
      newfile=charval('AV_WRIT_SYM:ENTER_FILENAME',1,80,0)
      lun=1
      close(unit=lun)
      open (unit=lun,file=newfile,status=cstatus)
c
      nsym1=intval('AV_WRIT_SYM:ENTER_1_NUMBER',1,nsym,0)
      nsym2=intval('AV_WRIT_SYM:ENTER_2_NUMBER',nsym1,nsym,0)
c
      write(6,20)' SYMM OPS WILL BE WRITTEN TO FILE, filename:'
     &  ,newfile
      if(record)write(nrec,20)
     &  ' SYMM OPS WILL BE WRITTEN TO FILE, filename:'
     &  ,newfile
c
20    format(a,/,a)
      do j=nsym1,nsym2
        call wrtsym(lun,j)
      enddo
c
      close(unit=lun)
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE write_sym_inv
c     ========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
c
      character*80 newfile,charval
      integer intval,lun,nsym1,nsym2,j
c
      if (nsym.lt.1) then
        write(6,10)
          if(record) write(nrec,10)
        return
      endif
10    format(' No Symmetry Operators defined')
c
      newfile=charval('AV_WRIT_SYM:ENTER_FILENAME',1,80,0)
      lun=1
      close(unit=lun)
      open (unit=lun,file=newfile,status=cstatus)
c
      nsym1=intval('AV_WRIT_INV_SYM:ENTER_1_NUMBER',1,nsym,0)
      nsym2=intval('AV_WRIT_INV_SYM:ENTER_2_NUMBER',nsym1,nsym,0)
c
      write(6,20)' INVERSE OPS WILL BE WRITTEN TO FILE, filename:'
     &                     ,newfile
      if(record)write(nrec,20)
     &  ' INVERSE OPS WILL BE WRITTEN TO FILE, filename:',newfile
20    format(a,/,a)
c
      do j=nsym1,nsym2
        call wrtsym_inv(lun,j)
      enddo
c
      close(unit=lun)
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE wrtsym(lun,n)
c     ========================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c     =====================
      integer n,k,j,lun
c
      if(.not.terse)then
        if(.not.xplor)then
          write(6,90)n
          if(record)write(nrec,90)n
90        format(' Store operator',i4,' on file in GAP format')
        else
          write(6,95)n
          if(record)write(nrec,95)n
95        format(' Store operator',i4,' on file in XPLOR format')
        endif
      endif
      if(.not.xplor)
     &  write(lun,100) ((ops(k,j,n),j=1,3),k=1,3), (vecs(j,n),j=1,3)
100     format (' ',/,5x,'ASYM ',3(/,10x,3f10.5)/,/,15x,3f10.4,/)
      if(xplor)
     &  write(lun,110) ((ops(k,j,n),j=1,3),k=1,3), (vecs(j,n),j=1,3)
110       format (/,' xncsrel',/
     &    ,'     matrix=',3(/,'                  (',3f11.5, ')' ) /
     &    ,'     translation= (',3f11.4,')    end')
c
      return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE wrtsym_inv(lun,n)
c       ============================
      IMPLICIT NONE
c
      INCLUDE 'average.fcm'
c       =====================
c
c This was broken - inverse translation component is not just negative of
c initial translation. Should be negative of inverse matrix * initial
c translation:
c
c [R']( [R] + T ) + T' = [I]
c [R'][R] + [R']T + T' = [I]
c [R'][R] = [I], thus [R']T + T' = 0
c T' = -[R']T
c
      integer j,k,lun,n
      real trans(3),inv_trans(3),zero_vec(3)
c
c Initialise translation component and zero vector for spinning
c
      do j=1,3
        trans(j)=vecs(j,n)
        zero_vec(j)=0.0e0
      enddo
c
c Use spin_inv to calculate [R']T
c NB As we are using a rotation centred on the origin, we could just as easily
c call spin instead of spin_inv
c
      call spin_inv(trans,inv_trans,ops_inv(1,1,n),zero_vec)
c
c Write out the answer. Hopefully it will now be right!
c
      if(.not.terse)then
        write(6,90)n
        if(record)write(nrec,90)n
      endif
90    format(' Store operator',i4,' on file')
c
      write(lun,100)((ops_inv(k,j,n),j=1,3),k=1,3),(-inv_trans(j),j=1,3)
100   format (' ',/,5x,'ASYM ',3(/,10x,3f10.5)/,/,15x,3f10.4,/)
c
c
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE XR_isym(map)
C     =======================
      IMPLICIT NONE
c
c  make a set of integerized crystallographic symm matrices for map MAP.
c
c
      INCLUDE 'average.fcm'
c       =====================
c
      integer map,j,ioff,i,l,ibase,j_sym,jpoint
c
      do j_sym=1,n_xsym(map)
c
c---- set up pointers
c
        ioff=(j_sym-1)*12
c
        do i=1,3
          l=i*4
          ibase=l-4
c
c---- firstly deal with the translation part.
c
          XR_insym(l+ioff,map)=(XR_sym(l+ioff,map))*nunit(i,map)/12
c
c---- now the rotational parts
c
          do j=1,3
            jpoint=j+ibase+ioff
            XR_insym(jpoint,map) = (XR_sym(jpoint,map))
          enddo
c
        enddo
c
      enddo
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE XR_symsy(CELL,IGROUP,ICS,NEQ,ML)
C     ===========================================
      IMPLICIT NONE
c
C***********************************************************************
C
C
C  OBTAIN EQUIVALENT POSITIONS FROM SPACE GROUP NUMBER.
C  ( ONLY THE 65 ENANTIOMORPHIC SPACE GROUPS ARE IMPLEMENTED )
C  THE SPACE GROUP NUMBER FOR YOUR CRYSTAL CAN BE LOOKED UP IN THE
C    INTERNATIONAL TABLES FOR X-RAY CRYSTALLOGRAPHY VOL.I
C
C              W. KABSCH      5-1986
C
C***********************************************************************
C  CELL  - UNIT CELL PARAMETERS IN ANGSTROEM AND DEGREES         (GIVEN)
C IGROUP - SPACE GROUP NUMBER OF ONE OF THE 65 ENANTIOMORPHIC    (GIVEN)
C          GROUPS AS OBTAINED FROM THE INTERNATIONAL TABLES.
C  ICS   - NUMBER SPECIFYING CRYSTAL SYSTEM                     (RESULT)
C          0 : UNIT CELL PARAMETERS ARE INCONSISTENT WITH
C              SPACE GROUP
C          1 : TRICLINIC
C          2 : MONOCLINIC FIRST SETTING
C          3 : MONOCLINIC SECOND SETTING
C          4 : ORTHORHOMBIC
C          5 : TETRAGONAL
C          6 : TRIGONAL
C          7 : HEXAGONAL
C          8 : CUBIC
C  NEQ   - NUMBER OF EQUIVALENT POSITIONS IN SPACE GROUP.       (RESULT)
C          A VALUE OF 0 FOR NEQ INDICATES AN ILLEGAL SPACE
C          GROUP NUMBER. OTHERWISE NEQ ASSUMES VALUES BETWEEN
C          1 AND 96.
C   ML   - INTEGER ARRAY OF LENGTH AT LEAST 12*NEQ SPECIFIED IN (RESULT)
C          THE MAIN PROGRAM, TO DESCRIBE SYMMETRY. EACH OPERATOR
C          IS REPRESENTED BY 12 CONSECUTIVE NUMBERS AS FOLLOWS:
C          IF X',Y',Z' IS AN EQUIVALENT POSITION TO X,Y,Z, THEN
C              X'=X*ML(1)+Y*ML( 2)+Z*ML( 3)+ML( 4)/12.0
C              Y'=X*ML(5)+Y*ML( 6)+Z*ML( 7)+ML( 8)/12.0
C              Z'=X*ML(9)+Y*ML(10)+Z*ML(11)+ML(12)/12.0
C          THE NEXT OPERATOR OCCUPIES POSITIONS ML(13)...ML(24),
C          AND SO ON...
C          AS MENTIONED, THE LARGEST NUMBER OF EQUIVALENT
C          POSITIONS IS 96, FOR TWO CUBIC SPACE GROUPS. HENCE,
C          THE CALLING PROGRAM MUST DECLARE   
C          INTEGER   ML(1152)     
C          TO HANDLE ALL POSSIBLE CASES !!!
C          TO SAVE MEMORY SPACE, ML,NEQ,ICS,IGROUP MAY ALL BE
C          DECLARED AS INTEGER*2 IN THE CALLING PROGRAM. BUT MAKE
C          SHURE YOU ALSO CHANGE THE DECLARATION HERE IN THIS
C          SUBROUTINE !!!
C***********************************************************************
      REAL        CELL(6)
      INTEGER     NEQ,ICS,IGROUP
      INTEGER*2   ML(*)
      INTEGER*2   IG,I,J,L,M,N,MM,NN,IADR,ITADR,IMA,ITA
      INTEGER*2   MAT(288),T(84),GROUP(70),NM(70),NT(70),NA(70),MATN(78)
      INTEGER*2   TRANS(902),TMONO(16),TORTHO(64),TETRA(128),TRIG(81),
     &            THEX(108),TCUB(504)
      INTEGER*2   T1(1),T3F(2),T4F(2),T5F(4),T3S(2),T4S(2),T5S(4)
      INTEGER*2   T16(4),T17(4),T18(4),T19(4),T20(8),T21(8),T22(16),
     &            T23(8),T24(8)
      INTEGER*2   T75(4),T76(4),T77(4),T78(4),T79(8),T80(8),
     &            T89(8),T90(8),T91(8),T92(8),T93(8),T94(8),T95(8),
     &            T96(8),T97(16),T98(16)
      INTEGER*2   T143(3),T144(3),T145(3),T146H(9),T146R(3),T149(6),
     &            T150(6),T151(6),T152(6),T153(6),T154(6),T155H(18),
     &            T155R(6)
      INTEGER*2   T168(6),T169(6),T170(6),T171(6),T172(6),T173(6),
     &            T177(12),T178(12),T179(12),T180(12),T181(12),T182(12)
      INTEGER*2   T195(12),T196(48),T197(24),T198(12),T199(24),
     &            T207(24),T208(24),T209(96),T210(96),T211(48),
     &            T212(24),T213(24),T214(48)
      EQUIVALENCE ( TMONO( 1),T3F(1)),( TMONO( 3),T4F(1)),
     &            ( TMONO( 5),T5F(1)),( TMONO( 9),T3S(1)),
     &            ( TMONO(11),T4S(1)),( TMONO(13),T5S(1))
      EQUIVALENCE (TORTHO( 1),T16(1)),(TORTHO( 5),T17(1)),
     &            (TORTHO( 9),T18(1)),(TORTHO(13),T19(1)),
     &            (TORTHO(17),T20(1)),(TORTHO(25),T21(1)),
     &            (TORTHO(33),T22(1)),(TORTHO(49),T23(1)),
     &            (TORTHO(57),T24(1))
      EQUIVALENCE          (TETRA( 1),T75(1)),(TETRA( 5),T76(1)),
     &  (TETRA( 9),T77(1)),(TETRA(13),T78(1)),(TETRA(17),T79(1)),
     &  (TETRA(25),T80(1)),(TETRA(33),T89(1)),(TETRA(41),T90(1)),
     &  (TETRA(49),T91(1)),(TETRA(57),T92(1)),(TETRA(65),T93(1)),
     &  (TETRA(73),T94(1)),(TETRA(81),T95(1)),(TETRA(89),T96(1)),
     &  (TETRA(97),T97(1)),(TETRA(113),T98(1))
      EQUIVALENCE (TRIG( 1), T143(1)),(TRIG( 4), T144(1)),
     &            (TRIG( 7), T145(1)),(TRIG(10),T146H(1)),
     &            (TRIG(19),T146R(1)),(TRIG(22), T149(1)),
     &            (TRIG(28), T150(1)),(TRIG(34), T151(1)),
     &            (TRIG(40), T152(1)),(TRIG(46), T153(1)),
     &            (TRIG(52), T154(1)),(TRIG(58),T155H(1)),
     &            (TRIG(76),T155R(1))
      EQUIVALENCE (THEX( 1),T168(1)),(THEX( 7),T169(1)),
     &            (THEX(13),T170(1)),(THEX(19),T171(1)),
     &            (THEX(25),T172(1)),(THEX(31),T173(1)),
     &            (THEX(37),T177(1)),(THEX(49),T178(1)),
     &            (THEX(61),T179(1)),(THEX(73),T180(1)),
     &            (THEX(85),T181(1)),(THEX(97),T182(1))
      EQUIVALENCE (TCUB(  1),T195(1)),
     &            (TCUB( 13),T196(1)),(TCUB( 61),T197(1)),
     &            (TCUB( 85),T198(1)),(TCUB( 97),T199(1)),
     &            (TCUB(121),T207(1)),(TCUB(145),T208(1)),
     &            (TCUB(169),T209(1)),(TCUB(265),T210(1)),
     &            (TCUB(361),T211(1)),(TCUB(409),T212(1)),
     &            (TCUB(433),T213(1)),(TCUB(457),T214(1))
      EQUIVALENCE (TRANS(  1),    T1(1)),(TRANS(  2),TMONO(1)),
     &            (TRANS( 18),TORTHO(1)),(TRANS( 82),TETRA(1)),
     &            (TRANS(210),  TRIG(1)),(TRANS(291), THEX(1)),
     &            (TRANS(399),  TCUB(1))
      DATA T1   / 1/
      DATA T3F  / 1,1/
      DATA T4F  / 1,2/
      DATA T5F  / 1,1,5,5/
      DATA T3S  / 1,1/
      DATA T4S  / 1,3/
      DATA T5S  / 1,1,6,6/
      DATA T16  / 4*1/
      DATA T17  / 1,2,1,2/
      DATA T18  / 1,1,6,6/
      DATA T19  / 1,5,6,7/
      DATA T20  / 1,2,1,2,6,8,6,8/
      DATA T21  / 4*1,4*6/
      DATA T22  / 4*1,4*7,4*5,4*6/
      DATA T23  / 4*1,4*8/
      DATA T24  / 1,5,6,7,8,3,2,4/
      DATA T75  / 4*1/
      DATA T76  / 1,2,10,9/
      DATA T77  / 1,1,2,2/
      DATA T78  / 1,2,9,10/
      DATA T79  / 4*1,4*8/
      DATA T80  / 1,1,11,11,8,8,12,12/
      DATA T89  / 8*1/
      DATA T90  / 4*1,4*6/
      DATA T91  / 1,2,9,10,10,9,1,2/
      DATA T92  / 1,2,2,1,14,13,13,14/
      DATA T93  / 1,1,4*2,1,1/
      DATA T94  / 4*1,4*8/
      DATA T95  / 1,2,10,9,9,10,1,2/
      DATA T96  / 1,2,2,1,13,14,14,13/
      DATA T97  / 8*1,8*8/
      DATA T98  / 4*1,4*11,4*8,4*12/
      DATA T143 / 3*1/
      DATA T144 / 1,15,16/
      DATA T145 / 1,16,15/
      DATA T146H/ 3*1,3*17,3*18/
      DATA T146R/ 3*1/
      DATA T149 / 6*1/
      DATA T150 / 6*1/
      DATA T151 / 1,1,15,15,16,16/
      DATA T152 / 1,1,15,15,16,16/
      DATA T153 / 1,1,16,16,15,15/
      DATA T154 / 1,1,16,16,15,15/
      DATA T155H/ 6*1,6*17,6*18/
      DATA T155R/ 6*1/
      DATA T168 / 6*1/
      DATA T169 / 1,15,16,2,19,20/
      DATA T170 / 1,16,15,2,20,19/
      DATA T171 / 1,16,15,1,16,15/
      DATA T172 / 1,15,16,1,15,16/
      DATA T173 / 3*1,3*2/
      DATA T177 / 12*1/
      DATA T178 / 1,1,2,2,15,15,19,19,16,16,20,20/
      DATA T179 / 1,1,2,2,16,16,20,20,15,15,19,19/
      DATA T180 / 4*1,4*16,4*15/
      DATA T181 / 4*1,4*15,4*16/
      DATA T182 / 2*1,2*2,2*1,2*2,2*1,2*2/
      DATA T195 / 12*1/
      DATA T196 / 12*1,12*7,12*5,12*6/
      DATA T197 / 12*1,12*8/
      DATA T198 / 3*1,3*6,3*7,3*5/
      DATA T199 / 3*1,3*6,3*7,3*5,3*8,3*2,3*4,3*3/
      DATA T207 / 24*1/
      DATA T208 / 12*1,12*8/
      DATA T209 / 24*1,24*7,24*5,24*6/
      DATA T210 / 12*1,12*21,12*7,12*22,12*5,12*23,12*6,12*24/
      DATA T211 / 24*1,24*8/
      DATA T212 / 3*1,3*6,3*7,3*5,3*24,3*23,3*21,3*22/
      DATA T213 / 3*1,3*6,3*7,3*5,3*27,3*26,3*25,3*28/
      DATA T214 / 3*1,3*6,3*7,3*5,3*24,3*23,3*21,3*22,
     &            3*8,3*2,3*4,3*3,3*27,3*26,3*25,3*28/
      DATA MATN / 1,3, 1,2,4,3, 1,2,5,6, 1,2,7,8,5,6,3,4, 1,11,12,
     &            7,17,18, 1,13,9,14,10,7, 1,8,9,15,10,16,
     &            1,9,10,2,19,20, 1,16,2,14,9,8,7,19,10,15,20,13,
     &            1,11,12,4,21,22,3,23,24,2,25,26,
     &            5,29,30,6,27,28,7,17,18,8,31,32/
      DATA MAT/ 1,3*0,1,3*0,1,  -1,3*0,-1,3*0,1, -1,3*0,1,3*0,-1,
     & 1,3*0,-1,3*0,-1,  0,1,0,-1,4*0,1,  0,-1,0,1,4*0,1,
     & 0,-1,0,-1,4*0,-1,  0,1,0,1,4*0,-1,  0,-1,0,1,-1,3*0,1,
     & -1,1,0,-1,4*0,1,  0,0,1,1,3*0,1,0,  0,1,3*0,1,1,2*0,
     & 1,0,0,1,-1,3*0,-1,  -1,1,0,0,1,3*0,-1,  -1,0,0,-1,1,3*0,-1,
     & 1,-1,0,0,-1,3*0,-1,  0,0,-1,0,-1,0,-1,0,0,  -1,4*0,-1,0,-1,0,
     & 0,1,0,-1,1,3*0,1,  1,-1,0,1,4*0,1,  0,0,1,-1,3*0,-1,0,
     & 0,1,3*0,-1,-1,0,0,  0,0,-1,1,3*0,-1,0,  0,-1,3*0,1,-1,0,0,
     & 0,0,-1,-1,3*0,1,0,  0,-1,3*0,-1,1,0,0,  -1,4*0,1,0,1,0,
     & 0,0,-1,0,1,0,1,0,0,  1,4*0,-1,0,1,0,  0,0,1,0,-1,0,1,0,0,
     & 1,4*0,1,0,-1,0,  0,0,1,0,1,0,-1,0,0/
      DATA   T/ 0,0,0,  0,0,6,  0,6,0,  6,0,0,  6,0,6,  6,6,0,
     &          0,6,6,  6,6,6,  0,0,3,  0,0,9,  0,6,3,  6,0,9,
     &          6,6,3,  6,6,9,  0,0,4,  0,0,8,  4,8,8,  8,4,4,
     &         0,0,10,  0,0,2,  3,3,3,  3,9,9,  9,3,9,  9,9,3,
     &          9,9,9,  3,9,3,  3,3,9,  9,3,3/
      DATA GROUP/
     &  1,  3,4,5,515,516,517,    16,17,18,19,20,21,22,23,24,
     & 75,76,77,78,79,80,         89,90,91,92,93,94,95,96,97,98,
     & 143,144,145,146,658,       149,150,151,152,153,154,155,667,
     & 168,169,170,171,172,173,   177,178,179,180,181,182,
     & 195,196,197,198,199,       207,208,209,210,211,212,213,214/
      DATA NM/  1, 6*2,  15*4, 10*8, 5*3, 14*6, 11*12, 8*24/
      DATA NT/  3*1,2,2*1,2,4*1,2*2,4,2*2,4*1,2*2,8*1,2*2,3*1,
     &  3,7*1,3,14*1,4,2,1,2,2*1,2*4,2,2*1,2/
      DATA NA/  0,3*2,3*0,9*2,6*6,10*10,4*36,54,24,30,24,30,24,
     &  2*30,18,6*36,6*42,13*54/
C***********************************************************************
      NEQ=0
      ICS=0
      IG=IGROUP
      DO 10 J=1,6
      IF (CELL(J).LE.0.0)RETURN
10    CONTINUE
C..... DETERMINE CRYSTAL SYSTEM NUMBER FROM SPACEGROUP NUMBER
C..... AND UNIT CELL PARAMETERS
      IF (IG.EQ.1)ICS=1
      IF ((IG.GE.3).AND.(IG.LE.5))THEN
         IF ((CELL(4).EQ.90.0).AND.(CELL(5).EQ.90.0))ICS=2
         IF ((CELL(4).EQ.90.0).AND.(CELL(6).EQ.90.0))ICS=3
         IF (ICS.EQ.3)IG=IG+512
      ENDIF
      IF ((IG.GE.16).AND.(IG.LE.24).AND.(CELL(4).EQ.90.0).AND.
     &    (CELL(5).EQ.90.0).AND.(CELL(6).EQ.90.0))ICS=4
      IF ((IG.GE.75).AND.(IG.LE.98).AND.(CELL(1).EQ.CELL(2)).AND.
     &(CELL(4).EQ.90.0).AND.(CELL(5).EQ.90.0).AND.(CELL(6).EQ.90.0))
     & ICS=5
      IF ((IG.GE.143).AND.(IG.LE.155).AND.(CELL(1).EQ.CELL(2)))THEN
        IF ((CELL(2).EQ.CELL(3)).AND.(CELL(4).EQ.CELL(5)).AND.
     &      (CELL(5).EQ.CELL(6)))ICS=6
        IF ((CELL(4).EQ.90.0).AND.(CELL(5).EQ.90.0).AND.
     &      (CELL(6).EQ.120.0))  ICS=7
        IF ((ICS.EQ.6).AND.((IG.EQ.146).OR.(IG.EQ.155)))IG=IG+512
      ENDIF
      IF ((IG.GE.168).AND.(IG.LE.182).AND.(CELL(1).EQ.CELL(2)).
     &   AND.(CELL(4).EQ.90.0).AND.(CELL(5).EQ.90.0).AND.
     &       (CELL(6).EQ.120.0)) ICS=7
      IF ((IG.GE.195).AND.(IG.LE.214).AND.(CELL(1).EQ.CELL(2))
     &   .AND.(CELL(2).EQ.CELL(3)).AND.(CELL(4).EQ.90.0).AND.
     &        (CELL(5).EQ.90.0).AND.(CELL(6).EQ.90.0))ICS=8
C..... OBTAIN EQUIVALENT POSITIONS
      ITADR=0
      DO 20 J=1,70
        M=NT(J)
        N=NM(J)
        L=NA(J)
        IF (IG.EQ.GROUP(J))GO TO 30
20    ITADR=ITADR+N*M
      RETURN
30    NEQ=N*M
      IADR=0
      DO 60 MM=1,M
      DO 60 NN=1,N
      ITADR=ITADR+1
      IMA=9*(MATN(L+NN)-1)
      ITA=3*(TRANS(ITADR)-1)
      DO 50 J=1,3
      DO 40 I=1,3
      IMA=IMA+1
40    ML(IADR+I)=MAT(IMA)
      ITA=ITA+1
      ML(IADR+4)=T(ITA)
50    IADR=IADR+4
60    CONTINUE
      RETURN
      END
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE xrfrac(rxrcell,rxrtr,rxrintr,rxrvol,print_l)
c     =====================================================
      IMPLICIT NONE
c
c Computes orthogonal to fractional transformation matrix lXRTR and
c its inverse.  Similar to routine ORTHO.FOR in PROTEIN program
c package (W. Steigemann, MPI Biochemie, FRG).
c
c New coordinates are obtained by r'(i)=sum_j matrix(i,j)*r(j)
c where matrix is lXRTR and lXRINTR for orthogonal to fractional and
c fractional to orthogonal, respectively.
c
c The convention to setup the matrices is as follows: x same direction
c as a, y is in (a,b) plane. 
c
c Author: Axel T. Brunger 
c ======================= 
c
c input/output
c
      INCLUDE 'average.fcm'
c     =====================
c
c--- Passed variables
c---  NB rxrcell and print are the only inputs...
c
c--- rxrcell real(9)   contains matrix describing lattice parameters
c--- rxtr    real(3,3) contains matrix to go from o to f
c--- rxintr  real(3,3) contains matrix to go from f to o
c--- rxvol   real      contains unit cell volume
c--- print   logical   flag for verbose output
c
      real rxrcell(9), rxrtr(3,3), rxrintr(3,3), rxrvol
      logical print_l
c
c--- Local variables
c
      integer i, j,k
      doubleprecision lxrcell(9), lxrtr(3,3), lxrintr(3,3), lxrvol
      doubleprecision cabg(3), sabg(3), cabgs(3), abcs(3), sabgs1
c
c--- Parameters
c
      doubleprecision rad, one, two
      parameter (one=1.0d0, two=2.0d0, rad=0.017453292d0)
c
c--- Begin
c
c
c--- Copy lattice matrix to double precision
c
      do j=1,6
        lxrcell(j)=rxrcell(j)
      enddo
c
c--- Zero working matrices
c
      do i=1,3
        do j=1,3
          lxrtr(i,j)=zero
          lxrintr(i,j)=zero
        enddo
      enddo
c
c
c
      do i=1,3
        cabg(i)=cos(lxrcell(i+3)*rad)
        sabg(i)=sin(lxrcell(i+3)*rad)
      enddo
      cabgs(1)=(cabg(2)*cabg(3)-cabg(1))/(sabg(2)*sabg(3))
      cabgs(2)=(cabg(3)*cabg(1)-cabg(2))/(sabg(3)*sabg(1))
      cabgs(3)=(cabg(1)*cabg(2)-cabg(3))/(sabg(1)*sabg(2))
      lxrvol=lxrcell(1)*lxrcell(2)*lxrcell(3)*
     &                   sqrt(one+two*cabg(1)*cabg(2)*cabg(3)
     &                  -cabg(1)**2-cabg(2)**2-cabg(3)**2)
      abcs(1)=lxrcell(2)*lxrcell(3)*sabg(1)/lxrvol
      abcs(2)=lxrcell(1)*lxrcell(3)*sabg(2)/lxrvol
      abcs(3)=lxrcell(1)*lxrcell(2)*sabg(3)/lxrvol
      sabgs1=sqrt(one-cabgs(1)**2)
c
c--- Cartesian to fractional conversion
c
      lxrtr(1,1)=one/lxrcell(1)
      lxrtr(1,2)=-cabg(3)/(sabg(3)*lxrcell(1))
      lxrtr(1,3)=-(cabg(3)*sabg(2)*cabgs(1)+cabg(2)*sabg(3))/
     &           (sabg(2)*sabgs1*sabg(3)*lxrcell(1))
      lxrtr(2,2)=one/(sabg(3)*lxrcell(2))
      lxrtr(2,3)=cabgs(1)/(sabgs1*sabg(3)*lxrcell(2))
      lxrtr(3,3)=one/(sabg(2)*sabgs1*lxrcell(3))
c
c--- Fractional to cartesian
c
      lxrintr(1,1)= lxrcell(1)
      lxrintr(1,2)= cabg(3)*lxrcell(2)
      lxrintr(1,3)= cabg(2)*lxrcell(3)
      lxrintr(2,2)= sabg(3)*lxrcell(2)
      lxrintr(2,3)=-sabg(2)*cabgs(1)*lxrcell(3)
      lxrintr(3,3)=sabg(2)*sabgs1*lxrcell(3)
c                
c--- Compute norm of a*, b*, c*
c
      lxrcell(7)=sqrt(lxrtr(1,1)**2 + lxrtr(1,2)**2 +lxrtr(1,3)**2)
      lxrcell(8)=sqrt(lxrtr(2,1)**2 + lxrtr(2,2)**2 +lxrtr(2,3)**2)
      lxrcell(9)=sqrt(lxrtr(3,1)**2 + lxrtr(3,2)**2 +lxrtr(3,3)**2)
c
c--- Display results
c
      if(print_l.and.verbose)then
        write(6,10)(lxrcell(j),j=1,6),((lxrtr(j,k),k=1,3),j=1,3)
c
c     ,((lxrintr(j,k),k=1,3),j=1,3)
c
        if(record)write(nrec,10)(lxrcell(j),j=1,6)
     &      ,((lxrtr(j,k),k=1,3),j=1,3)
c
c      ,((lxrintr(j,k),k=1,3),j=1,3)
c
      endif
10    format(' xrcell:',t10,6f11.2,/,' xrtr:',t10,3f11.5,2(/,t10,
     &          3f11.5) )
c
c     1 ,' xrintr:',3(10x,3(f11.5)/))
c
c--- Copy back into real matrices
c
      do j=1,3
        do k=1,3
          rxrtr(j,k)=lxrtr(j,k)
          rxrintr(j,k)=lxrintr(j,k)
        enddo
      enddo
      rxrvol=sngl(lxrvol)
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE zero_sym
c       ===================
      IMPLICIT NONE
c
c
      INCLUDE 'average.fcm'
c       =====================
c
        nsym=0
        write(6,10)'NCS OPERATORS INITIALIZED'
        if(record)write(nrec,10)'NCS OPERATORS INITIALIZED'
10      format(' ',a)
        return
        end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
