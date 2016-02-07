C***********************************************************************
C Parameters used by the P2 parser. These are the machine-independent
C parameters.
C
C Robert Esnouf, 4/7/93 - 31/12/93
C
C***********************************************************************

      INTEGER PUNMIN, PUNMAX
      INTEGER IXSPEC, IXMBAS, IXMSIZ, IXBAS, IXSIZ, IXNUM, IXABAS

C Range of unit numbers possibly used by the P2. This range may be
C modified by the machine dependent parameters in M2.

      PARAMETER (PUNMIN=1)
      PARAMETER (PUNMAX=99)

C Parameters used by memory allocation routines. These are the indexes
C to the housekeeping variables stored with each memory allocation.
C They are the special value for the machine specific part, the actual
C base address of the actual allocation, the actual size of allocation
C in bytes, the base address of the array used by the user, the size of
C each element, the number of elements to allocate and the address
C which needs to be supplied for deallocation (just a checksum really)

      PARAMETER (IXSPEC=0)
      PARAMETER (IXMBAS=1)
      PARAMETER (IXMSIZ=2)
      PARAMETER (IXBAS =3)
      PARAMETER (IXSIZ =4)
      PARAMETER (IXNUM =5)
      PARAMETER (IXABAS=6)

C***********************************************************************
C Common blocks used by the P2 routines. This is the main set of common
C blocks controlling most aspects of the parser.
C
C Robert Esnouf, 4/7/93 - 24/7/93
C
C***********************************************************************

      COMMON /P2/ TALK, CHATTY, ISTERM, IPUNIT, OPUNIT, UNMIN, UNMAX
      COMMON /P2_2/ OLDCPU, CASEIN, UPPERC, LOWERC, FRECL
      COMMON /P2_SL/ LMCTYP, LOSTYP, LOSVER, LMCNAM, LUSERN, LMACH
      COMMON /P2_S/ MCTYP, OSTYP, OSVER, MCNAM, USERN, MACH
      COMMON /P2_SL2/ LSTTIM, LSTDAT, LARROW
      COMMON /P2_S2/ STTIM, STDAT, PARROW
      COMMON /P2_C/ FNAME, FUNIT
      COMMON /P2_MEM/ IBKEEP, SIZINT, BKEEPS, CURALC, MAXALC

      SAVE /P2/, /P2_2/, /P2_SL/, /P2_S/, /P2_SL2/, /P2_S2/, /P2_C/, /P2_MEM/

C Flags indicating whether the parser should produce many run-time
C messages.

      LOGICAL TALK, CHATTY

C Flag indicating whether the base level input is from an interactive
C terminal

      LOGICAL ISTERM

C The units used as standard input and output by P2. These are just
C set to be the values in the machine-dependent parameters allowing
C P2 not to have to include M2.INC

      INTEGER IPUNIT
      INTEGER OPUNIT

C The limits on the available units for the machine on which the
C parser is running. These are set from the machine-independent
C parameters at initialisation time and not changed

      INTEGER UNMIN
      INTEGER UNMAX

C CPU time at last P2PCPU call, this allows the user an easy way to
C time routines coarsely

      INTEGER OLDCPU

C Flags controlling whether strings are compared case insensitively,
C and whether input is uppercased, lowercased or not changed.

      LOGICAL CASEIN, UPPERC, LOWERC

C Strings describing the machine and username

      INTEGER LMCTYP, LOSTYP, LOSVER, LMCNAM, LUSERN, LMACH
      CHARACTER*32 MCTYP, OSTYP, OSVER, MCNAM, USERN, MACH

C Strings to store the starting time and date

      INTEGER LSTTIM, LSTDAT
      CHARACTER*8 STTIM
      CHARACTER*19 STDAT

C Array storing the arrow to use constructing a prompt

      INTEGER LARROW
      CHARACTER*16 PARROW

C Array storing status of possible FORTRAN units

      CHARACTER*10 FUNIT(PUNMIN:PUNMAX)
      INTEGER FRECL(PUNMIN:PUNMAX)

C Array storing names of currently open FORTRAN files

      CHARACTER*80 FNAME(PUNMIN:PUNMAX)

C Memory allocation bits. IBKEEP is a small array to show how big an
C integer is on the machine and used for boundary alignment. SIZINT is
C the size of an integer as calculated by P2MEMI. BKEEPS is the start
C address of the IBKEEP array. The final two are the current and
C maximum allocations of memory

      INTEGER (KIND=8) IBKEEP(2)
      INTEGER (KIND=8) SIZINT
      INTEGER (KIND=8) BKEEPS
      INTEGER (KIND=8) CURALC
      INTEGER (KIND=8) MAXALC
