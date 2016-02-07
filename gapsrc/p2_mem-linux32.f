C***********************************************************************
C These are the machine independent routines for allocating and deall-
C ocating blocks of memory for the P2 routines. They also keep track of
C what has been allocated so that deallocation can be masked from the
C user as much as possible. One routine is needed for allocation of non-
C character data and one for character data. Similarly, two routines are
C required for deallocation. The first two arguments to the allocation
C calls are the first two elements of an array (thus declare arrays 2
C elements long). From these it can calculate the size of an element and
C the start of the array in memory. This allows a single routine to work
C for a variety of data types. For deallocation just the first element
C is needed along with the index returned by the call. This way of doing
C things also allows multi-dimensional arrays to be passed as arguments
C provided the last dimension is the one dimensioned to 2. For the
C allocates a return value of 0 for the index indicates failure, since
C that indexes the already existing variable. Deallocates retrun TRUE
C or FALSE. Rehashed for 64-bit addressing in 1997
C
C NOTE for Fortran 2003, the INTEGER*8 is no longer defined, so for
C compilation in this way the INTEGER (KIND=8) statement is substituted.
C This, technically, makes this routine machine dependent.
C
C Robert Esnouf, 23/7/93 - 2/4/97
C
C***********************************************************************

C***********************************************************************
C 32-bit routine to allocate memory for non-character arrays
C RME 2/4/97 - 2/4/97

      INTEGER FUNCTION P2AL32(ITEM1, ITEM2, NITEMS)
      INTEGER ITEM1, ITEM2, NITEMS

      INCLUDE 'p2.inc'

      INTEGER (KIND=8) P2AL, IX
      INTEGER (KIND=8) M2VLOC, ADRESS, ISIZE
      LOGICAL FLAG, P2MEMD

C Try and allocate memory using the 64-bit routines

      IX = P2AL(ITEM1, ITEM2, NITEMS)
      P2AL32 = IX

C See if the returned index can be expressed in 32 bits. If not then
C deallocate the block of memory and give a suitable error message

      IF (P2AL32.NE.IX) THEN
        ADRESS = M2VLOC(ITEM1)
        ISIZE = M2VLOC(ITEM2) - ADRESS
        ADRESS = ADRESS + ISIZE*(IX-1)
        FLAG = P2MEMD(ADRESS)
        WRITE(OPUNIT, 1000)
        P2AL32 = 1
      ENDIF

      RETURN

1000  FORMAT('P2AL32: 64-bit address beyond range of 32-bit index')

      END

C***********************************************************************
C 32-bit routine to allocate memory for character arrays
C RME 2/4/97 - 2/4/97

      INTEGER FUNCTION P2CAL32(ITEM1, ITEM2, NITEMS)
      CHARACTER*(*) ITEM1, ITEM2
      INTEGER NITEMS

      INCLUDE 'p2.inc'

      INTEGER (KIND=8) P2CAL, IX
      INTEGER (KIND=8) M2CLOC, ADRESS, ISIZE
      LOGICAL FLAG, P2MEMD

C Try and allocate memory using the 64-bit routines

      IX = P2CAL(ITEM1, ITEM2, NITEMS)
      P2CAL32 = IX

C See if the returned index can be expressed in 32 bits. If not then
C deallocate the block of memory and give a suitable error message

      IF (P2CAL32.NE.IX) THEN
        ADRESS = M2CLOC(ITEM1)
        ISIZE = M2CLOC(ITEM2) - ADRESS
        ADRESS = ADRESS + ISIZE*(IX-1)
        FLAG = P2MEMD(ADRESS)
        WRITE(OPUNIT, 1000)
        P2CAL32 = 1
      ENDIF

      RETURN

1000  FORMAT('P2CAL32: 64-bit address beyond range of 32-bit index')

      END

C***********************************************************************
C 64-bit routine to allocate memory for a non-character array. The first
C two arguments are the reference arguments for the array to be allocated
C to and the last argument is the number of elements for the array. Returns
C and index of 0 if the allocation failed
C RME 2/4/97 - 2/4/97

      INTEGER (KIND=8) FUNCTION P2AL(ITEM1, ITEM2, NITEMS)
      INTEGER ITEM1, ITEM2, NITEMS

      INTEGER (KIND=8) BASEAD, P2MEMA, M2VLOC
      INTEGER (KIND=8) ISIZE

C Find out the size of each element if the array

      BASEAD = M2VLOC(ITEM1)
      ISIZE = M2VLOC(ITEM2) - BASEAD

C Then pass the start address unit size and number of items to the
C variable-general routine

      P2AL = P2MEMA(BASEAD, ISIZE, NITEMS)

      RETURN

      END

C***********************************************************************
C 64-bit routine to allocate memory for a character array. The first
C two arguments are the reference arguments for the array to be allocated
C to and the last argument is the number of elements for the array. Returns
C and index of 0 if the allocation failed
C RME 2/4/97 - 2/4/97

      INTEGER (KIND=8) FUNCTION P2CAL(ITEM1, ITEM2, NITEMS)
      CHARACTER*(*) ITEM1, ITEM2
      INTEGER NITEMS

      INTEGER (KIND=8) BASEAD, P2MEMA, M2CLOC
      INTEGER (KIND=8) ISIZE

C Find out the size of each element if the array

      BASEAD = M2CLOC(ITEM1)
      ISIZE = M2CLOC(ITEM2) - BASEAD

C Then pass the start address unit size and number of items to the
C variable-general routine

      P2CAL = P2MEMA(BASEAD, ISIZE, NITEMS)

      RETURN

      END

C***********************************************************************
C Routine to deallocate memory for non-character arrays. It needs to be
C told the first element of the reference array and the index returned by
C the allocation call. Returns TRUE if deallocation was successful
C RME 2/4/97 - 2/4/97

      LOGICAL FUNCTION P2DAL(ITEM)
      INTEGER ITEM

      INTEGER (KIND=8) ADRESS, M2VLOC
      LOGICAL P2MEMD

C Find the address of the variable to be released and deallocate that block

      ADRESS = M2VLOC(ITEM)
      P2DAL = P2MEMD(ADRESS)

      RETURN

      END

C***********************************************************************
C Routine to deallocate memory for non-character arrays. It needs to be
C told the first element of the reference array and the index returned by
C the allocation call. Returns TRUE if deallocation was successful
C RME 2/4/97 - 2/4/97

      LOGICAL FUNCTION P2CDAL(ITEM)
      CHARACTER*(*) ITEM

      INTEGER (KIND=8) ADRESS, M2CLOC
      LOGICAL P2MEMD

C Find the address of the variable to be released and deallocate that block

      ADRESS = M2CLOC(ITEM)
      P2CDAL = P2MEMD(ADRESS)

      RETURN
      END

C***********************************************************************
C Initialisation routine for the memory allocation stuff. All it needs
C to do is find out the size of a long integer for the book-keeping
C variables.
C RME 24/7/93 - 24/7/93
C***********************************************************************

      SUBROUTINE P2MEMI

      INCLUDE 'p2.inc'

      INTEGER (KIND=8) M2VLOC

C Set up constants for linking to book-keeping variables

      BKEEPS = M2VLOC(IBKEEP(1))
      SIZINT = M2VLOC(IBKEEP(2)) - BKEEPS

C Set the counters for the total size allocated

      MAXALC = 0
      CURALC = 0

      RETURN
      END

C***********************************************************************
C General routine to allocate memory. The base address size of an
C element and number of elements have been calculated before entry. The
C amount of memory actually allocated is big enough to allow for align-
C ment of the required size on to the base address - possibly a problem
C with multi-dimensional arrays if each element is large. This routine
C is not designed for users to call - they should stick to P2AL and
C P2CAL. Modified for 64-bit addressing. Returns 0 if the allocation failed
C RME 23/7/93 - 2/4/97
C***********************************************************************

      INTEGER (KIND=8) FUNCTION P2MEMA(BASEAD, ISIZE, NITEMS)
      INTEGER (KIND=8) BASEAD, ISIZE
      INTEGER NITEMS

      INCLUDE 'p2.inc'

      INTEGER IS, IE, STATUS
      INTEGER (KIND=8) NBYTES, START, ASTART, AIX, BSTART, BIX
      LOGICAL M2MEMA
      CHARACTER*16 LINE

C Check we should allocate something

      IF (NITEMS.LE.0) THEN
        IF (TALK) WRITE(OPUNIT, 1000)
        P2MEMA = 0
        RETURN
      ENDIF

C Check the size of each item

      IF (ISIZE.LE.0) THEN
        WRITE(OPUNIT, 1010)
        P2MEMA = 0
        RETURN
      ENDIF

C Calculate how many bytes are needed allowing for all the boundary
C alignment and book-keeping variables

      NBYTES = NITEMS*ISIZE + 7*SIZINT + ISIZE + SIZINT

C Now try to allocate that much memory, note that NBYTES may be changed
C by the routine if a larger size is allocated (maybe)

      IF (.NOT.M2MEMA(NBYTES, START, STATUS)) THEN
        WRITE(OPUNIT, 1020)
        WRITE(LINE, '(I16)') STATUS
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1120) LINE(IS:IE)
        WRITE(LINE, '(I16)') NBYTES
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1130) LINE(IS:IE)
        WRITE(LINE, '(Z16)') START
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1140) LINE(IS:IE)
        WRITE(LINE,'(F16.1)') CURALC / 1024.0
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1210) LINE(IS:IE)
        WRITE(LINE,'(F16.1)') MAXALC / 1024.0
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1220) LINE(IS:IE)
        P2MEMA = 0
        RETURN
      ENDIF

C Find the address for the returned index leaving room for
C book-keeping. The address is calculated for part of the
C first element then rounded down to a tidy boundary.
C Finally this address is calculated for book-keeping

      ASTART = START + 8*SIZINT + ISIZE - 1
      IF (ASTART.GE.BASEAD) THEN
        AIX = 1 + (ASTART - BASEAD) / ISIZE
      ELSE
        AIX = - ((BASEAD - ASTART - 1) / ISIZE)
      ENDIF
      ASTART = BASEAD + (AIX-1) * ISIZE

C We have the memory, now work out the start index for the
C book-keeping long integers relative to the new start so we can
C find them again

      BSTART = ASTART - 7*SIZINT
      IF (BSTART.GE.BKEEPS) THEN
        BIX = 1 + (BSTART - BKEEPS) / SIZINT
      ELSE
        BIX = - ((BKEEPS - BSTART - 1) / SIZINT)
      ENDIF

C Now save the information

      CALL P2STBK(IBKEEP(BIX), STATUS, NITEMS, START, NBYTES, BASEAD,
     &            ISIZE, ASTART)

C If P2 is talkative or chatty then summarise the allocate in
C potentially exhausting detail

      CURALC = CURALC + NBYTES
      MAXALC = MAX(CURALC, MAXALC)

      IF (TALK) THEN
        IF (CHATTY) THEN
          WRITE(LINE,'(I16)') ISIZE
          CALL P2SLEN(LINE, IS, IE)
          WRITE(OPUNIT, 1100) LINE(IS:IE)
          WRITE(LINE,'(I16)') NITEMS
          CALL P2SLEN(LINE, IS, IE)
          WRITE(OPUNIT, 1110) LINE(IS:IE)
          WRITE(LINE,'(I16)') STATUS
          CALL P2SLEN(LINE, IS, IE)
          WRITE(OPUNIT, 1120) LINE(IS:IE)
          WRITE(LINE,'(I16)') NBYTES
          CALL P2SLEN(LINE, IS, IE)
          WRITE(OPUNIT, 1130) LINE(IS:IE)
          WRITE(LINE,'(Z16)') START
          CALL P2SLEN(LINE, IS, IE)
          WRITE(OPUNIT, 1140) LINE(IS:IE)
          WRITE(LINE,'(Z16)') BASEAD
          CALL P2SLEN(LINE, IS, IE)
          WRITE(OPUNIT, 1150) LINE(IS:IE)
          WRITE(LINE,'(Z16)') ASTART
          CALL P2SLEN(LINE, IS, IE)
          WRITE(OPUNIT, 1160) LINE(IS:IE)
          WRITE(LINE,'(I16)') AIX
          CALL P2SLEN(LINE, IS, IE)
          WRITE(OPUNIT, 1170) LINE(IS:IE)
        ENDIF
        WRITE(LINE,'(F16.1)') NBYTES / 1024.0
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1200) LINE(IS:IE)
        WRITE(LINE,'(F16.1)') CURALC / 1024.0
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1210) LINE(IS:IE)
        WRITE(LINE,'(F16.1)') MAXALC / 1024.0
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1220) LINE(IS:IE)
      ENDIF

C All book-keeping done safely so return with the index to use

      P2MEMA = AIX
      RETURN

1000  FORMAT('P2MEMA: Zero (or fewer) items requested for allocation')
1010  FORMAT('P2MEMA: Size of item specified was zero or negative')
1020  FORMAT('P2MEMA: Memory allocation call failed')

1100  FORMAT('P2MEMA: Size for each item in array: ',A)
1110  FORMAT('P2MEMA: Number of items to allocate: ',A)
1120  FORMAT('P2MEMA: Return code from allocation: ',A)
1130  FORMAT('P2MEMA: Actual byte total allocated: ',A)
1140  FORMAT('P2MEMA: Start address of allocation: ',A)
1150  FORMAT('P2MEMA: Base address of dummy array: ',A)
1160  FORMAT('P2MEMA: Start address for the array: ',A)
1170  FORMAT('P2MEMA: Starting index of the array: ',A)

1200  FORMAT('P2MEMA: Amount of memory allocated was: ',A,' Kbytes')
1210  FORMAT('P2MEMA: Total current memory allocated: ',A,' Kbytes')
1220  FORMAT('P2MEMA: Maximum total memory allocated: ',A,' Kbytes')

      END

C***********************************************************************
C General routine to deallocate memory. The address passed should be
C the address of the array element returned as the index by P2MEMA.
C P2MEMD then looks before that element to find the information about
C the allocate and then if that checks out it does a deallocate. The
C return value is TRUE if the deallocation is successful. Modified for
C 64-bit addressing
C RME 24/7/93 - 2/4/97
C***********************************************************************

      LOGICAL FUNCTION P2MEMD(ADRESS)
      INTEGER (KIND=8) ADRESS

      INCLUDE 'p2.inc'

      INTEGER (KIND=8) START, BASEAD, NBYTES, ISIZE
      INTEGER (KIND=8) ASTART, BSTART, BIX
      INTEGER IS, IE, NITEMS, STATUS
      LOGICAL M2MEMD
      CHARACTER*16 LINE

C Find the book-keeping variables

      BSTART = ADRESS - 7*SIZINT
      IF (BSTART.GE.BKEEPS) THEN
        BIX = 1 + (BSTART - BKEEPS) / SIZINT
      ELSE
        BIX = - ((BKEEPS - BSTART - 1) / SIZINT)
      ENDIF

      CALL P2GTBK(IBKEEP(BIX), STATUS, NITEMS, START, NBYTES, BASEAD,
     &            ISIZE, ASTART, 7)

C Show the book-keeping variables if CHATTY whether or not they make
C sense

      IF (CHATTY) THEN
        WRITE(LINE,'(I16)') ISIZE
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1100) LINE(IS:IE)
        WRITE(LINE,'(I16)') NITEMS
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1110) LINE(IS:IE)
        WRITE(LINE,'(I16)') STATUS
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1120) LINE(IS:IE)
        WRITE(LINE,'(I16)') NBYTES
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1130) LINE(IS:IE)
        WRITE(LINE,'(Z16)') START
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1140) LINE(IS:IE)
        WRITE(LINE,'(Z16)') BASEAD
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1150) LINE(IS:IE)
        WRITE(LINE,'(Z16)') ASTART
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1160) LINE(IS:IE)
        WRITE(LINE,'(Z16)') ADRESS
        CALL P2SLEN(LINE, IS, IE)
        WRITE(OPUNIT, 1170) LINE(IS:IE)
      ENDIF

C See if they make sense as a set

      P2MEMD = .FALSE.
      IF (ADRESS.NE.ASTART) THEN

        WRITE(OPUNIT, 1000)

      ELSE IF (ISIZE.LE.0) THEN

        WRITE(OPUNIT, 1010)

      ELSE IF (MOD(ABS(ASTART-BASEAD),ISIZE).NE.0) THEN

        WRITE(OPUNIT, 1020)

      ELSE IF (NBYTES .LT. ISIZE*(NITEMS+1) + 8*SIZINT - 2) THEN

        WRITE(OPUNIT, 1030)

      ELSE IF ((ASTART-START)-7*SIZINT.LT.0 .OR.
     &         (START-ASTART)+NBYTES-ISIZE*NITEMS.LT.0) THEN

        WRITE(OPUNIT, 1040)

      ELSE IF (STATUS.LE.0 .OR. STATUS.GT.8) THEN

C An unusual status is used to attempt to catch blocks previously
C deallocated (see just below - even if dealloc fails status will
C have been nobbled since we can't write to it after deallocation

        WRITE(OPUNIT, 1050)

      ELSE

        IBKEEP(BIX+IXSPEC) = 0
        IF (M2MEMD(NBYTES, START, STATUS)) THEN
          CURALC = CURALC - NBYTES
          IF (TALK) THEN
            WRITE(LINE,'(F16.1)') NBYTES / 1024.0
            CALL P2SLEN(LINE, IS, IE)
            WRITE(OPUNIT, 1200) LINE(IS:IE)
            WRITE(LINE,'(F16.1)') CURALC / 1024.0
            CALL P2SLEN(LINE, IS, IE)
            WRITE(OPUNIT, 1210) LINE(IS:IE)
            WRITE(LINE,'(F16.1)') MAXALC / 1024.0
            CALL P2SLEN(LINE, IS, IE)
            WRITE(OPUNIT, 1220) LINE(IS:IE)
          ENDIF
          P2MEMD = .TRUE.
        ELSE
          WRITE(OPUNIT, 1060)
        ENDIF

      ENDIF

      RETURN

1000  FORMAT('P2MEMD: Address for deallocation invalid')
1010  FORMAT('P2MEMD: Size of array elements is zero or negative')
1020  FORMAT('P2MEMD: Calculated index for deallocation non-integral')
1030  FORMAT('P2MEMD: Calculated deallocation size too large')
1040  FORMAT('P2MEMD: Calculated deallocation addresses invalid')
1050  FORMAT('P2MEMD: Memory block has already been deallocated')
1060  FORMAT('P2MEMD: Memory deallocation call failed')

1100  FORMAT('P2MEMD: Size for each item in array: ',A)
1110  FORMAT('P2MEMD: Number of items to allocate: ',A)
1120  FORMAT('P2MEMD: Return code from allocation: ',A)
1130  FORMAT('P2MEMD: Actual byte total allocated: ',A)
1140  FORMAT('P2MEMD: Start address of allocation: ',A)
1150  FORMAT('P2MEMD: Base address of dummy array: ',A)
1160  FORMAT('P2MEMD: Start address for the array: ',A)
1170  FORMAT('P2MEMD: Address asked to deallocate: ',A)

1200  FORMAT('P2MEMD:   Amount of memory deallocated: ',A,' Kbytes')
1210  FORMAT('P2MEMD: Total current memory allocated: ',A,' Kbytes')
1220  FORMAT('P2MEMD: Maximum total memory allocated: ',A,' Kbytes')

      END

C***********************************************************************
C Added routine to set book-keeping variables to avoid segmentation
C errors on Linux machines
C RME 21/4/10 - 21/4/10
C Now save the information

      SUBROUTINE P2STBK(IBK, STATUS, NITEMS, START, NBYTES, BASEAD,
     &                  ISIZE, ASTART)
      INTEGER (KIND=8) IBK(7)
      INTEGER (KIND=8) START, NBYTES, BASEAD, ISIZE, ASTART
      INTEGER STATUS, NITEMS
	
      INCLUDE 'p2_mem.inc'

      IBK(IXSPEC+1) = STATUS
      IBK(IXMBAS+1) = START
      IBK(IXMSIZ+1) = NBYTES
      IBK(IXBAS+1) = BASEAD
      IBK(IXSIZ+1) = ISIZE
      IBK(IXNUM+1) = NITEMS
      IBK(IXABAS+1) = ASTART

      RETURN
      END

C***********************************************************************
C Added routine to access book-keeping variables to avoid segmentation
C errors on Linux machines
C RME 21/4/10 - 21/4/10

      SUBROUTINE P2GTBK(IBK, STATUS, NITEMS, START, NBYTES, BASEAD,
     &                  ISIZE, ASTART, IPAR)
      INTEGER (KIND=8) IBK(IPAR)
      INTEGER (KIND=8) START, NBYTES, BASEAD, ISIZE, ASTART
      INTEGER STATUS, NITEMS, IPAR
	
      INCLUDE 'p2_mem.inc'

      STATUS = IBK(IXSPEC+1)
      START = IBK(IXMBAS+1)
      NBYTES = IBK(IXMSIZ+1)
      BASEAD = IBK(IXBAS+1)
      ISIZE = IBK(IXSIZ+1)
      NITEMS = IBK(IXNUM+1)
      ASTART = IBK(IXABAS+1)

      RETURN
      END
