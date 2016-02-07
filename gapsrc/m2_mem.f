C***********************************************************************
C Return the address at which the specified item can be found. Used for
C the memory allocation routines. (Non-character variables)
C RME 24/7/93 - 24/7/93
C***********************************************************************

@IF LINUX64 | LINUX32
      INTEGER (KIND=8) FUNCTION M2VLOC(ITEM)
@ELSE
      INTEGER*8 FUNCTION M2VLOC(ITEM)
@ENDIF
      INTEGER ITEM

@IF VMS | CONVEX | HP
      M2VLOC = %LOC(ITEM)
@ELSEIF SGI
      EXTERNAL C2LOC

      CALL C2LOC(M2VLOC, ITEM)
@ELSEIF LINUX64 | LINUX32
      M2VLOC = LOC(ITEM)
@ELSE
      M2VLOC = 0
@ENDIF

      RETURN
      END

C***********************************************************************
C Return the address at which the specified item can be found. Used for
C the memory allocation routines. Character variables ONLY!
C RME 24/7/93 - 24/7/93
C***********************************************************************

@IF LINUX64 | LINUX32
      INTEGER (KIND=8) FUNCTION M2CLOC(ITEM)
@ELSE
      INTEGER*8 FUNCTION M2CLOC(ITEM)
@ENDIF
      CHARACTER*(*) ITEM

@IF VMS | CONVEX | HP
      M2CLOC = %LOC(ITEM)
@ELSEIF SGI
      EXTERNAL C2CLOC

      CALL C2LOC(M2CLOC, ITEM)
@ELSEIF LINUX64 | LINUX32
      M2CLOC = LOC(ITEM)
@ELSE
      M2CLOC = 0
@ENDIF

      RETURN
      END

C***********************************************************************
C Routine to allocate a block of memory of at least NBYTES bytes. If the
C amount allocated is differrent then NBYTES is altered. START is the
C address for the allocation and STATUS can be used by the routine as a
C machine dependent variable which will be passed to the deallocate
C routine.
C RME 24/7/93 - 15/11/93
C***********************************************************************

      LOGICAL FUNCTION M2MEMA(NBYTES, START, STATUS)
@IF LINUX64 | LINUX32
      INTEGER (KIND=8) NBYTES, START
@ELSE
      INTEGER*8 NBYTES, START
@ENDIF
      INTEGER STATUS

      INCLUDE 'p2.inc'
@IF VMS
      INCLUDE 'm2_vms.inc'

      INTEGER ISTAT, SYS$EXPREG, LIB$GET_VM
      INTEGER PAGLET, RETADR(2)
      INTEGER PGLSIZ
      PARAMETER (PGLSIZ=512)
@ELSEIF HP

      EXTERNAL MALLOC_
@ELSEIF SGI | LINUX64 | LINUX32

      EXTERNAL C2MEMA
@ENDIF

@IF VMS
C On a VAX and Alpha we use different routines depending on the amount
C of memory requested. This cutoff is 512 bytes on the VAX and 8192
C bytes on the Alpha

      IF (NBYTES.GE.MEMLIM) THEN
        STATUS = 1
        PAGLET = (NBYTES-1)/PGLSIZ + 1
        ISTAT = SYS$EXPREG(%VAL(PAGLET), RETADR, %VAL(0), %VAL(0))
        IF (.NOT.ISTAT) THEN
          CALL LIB$SIGNAL(%VAL(ISTAT))
          GOTO 99
        ENDIF
        START = RETADR(1)
        NBYTES = RETADR(2) - RETADR(1) + 1
      ELSE
        STATUS = 2
        ISTAT = LIB$GET_VM(NBYTES, START)
        IF (.NOT.ISTAT) THEN
          CALL LIB$SIGNAL(%VAL(ISTAT))
          GOTO 99
        ENDIF
      ENDIF

      M2MEMA = .TRUE.
@ELSEIF HP
C On the HP we use the U77 MALLOC routine

      STATUS = 1
      M2MEMA = .TRUE.
      CALL MALLOC_(NBYTES, START)
      IF (START.EQ.0) GOTO 99
@ELSEIF SGI | LINUX64 | LINUX32
C On the SG we use the an interface to the C malloc routine

      STATUS = 1
      M2MEMA = .TRUE.
      CALL C2MEMA(NBYTES, START)
      IF (START.EQ.0) GOTO 99
@ELSE
      STATUS = 0
      NBYTES = 0
      START = 0
      M2MEMA = .FALSE.
@ENDIF

      RETURN

@IF VMS | HP | SGI | LINUX64 | LINUX32
C Return FALSE if the allocation failed with NBYTES and START = 0

99    NBYTES = 0
      START = 0
      M2MEMA = .FALSE.
      RETURN
@ENDIF
      END

C***********************************************************************
C Routine to deallocate a block of memory of at least NBYTES bytes which
C start from START. The variable STATUS was the value returned by the
C allocation call, and can be used for any machine dependent purpose. If
C the deallocate fails then return FALSE.
C RME 24/7/93 - 24/7/93
C***********************************************************************

      LOGICAL FUNCTION M2MEMD(NBYTES, START, STATUS)
@IF LINUX64 | LINUX32
      INTEGER (KIND=8) NBYTES, START
@ELSE
      INTEGER*8 NBYTES, START
@ENDIF
      INTEGER STATUS

      INCLUDE 'p2.inc'
@IF VMS

      INTEGER ISTAT, SYS$DELTVA, LIB$FREE_VM
      INTEGER DELADR(2), RETADR(2)
@ELSEIF HP

      EXTERNAL FREE_
@ELSEIF SGI | LINUX64 | LINUX32

      EXTERNAL C2MEMD
@ENDIF

@IF VMS
C On a VAX and Alpha we used different routines depending on the amount
C of memory to be allocated. The routine used is given by STATUS, so
C use the complementary routine for each case.

      IF (STATUS.EQ.1) THEN
        DELADR(1) = START
        DELADR(2) = START + NBYTES - 1
        ISTAT = SYS$DELTVA(DELADR, RETADR, %VAL(0))
        IF (.NOT.ISTAT) THEN
          CALL LIB$SIGNAL(%VAL(ISTAT))
          GOTO 99
        ENDIF
        IF (RETADR(1).NE.DELADR(1) .OR. RETADR(2).NE.DELADR(2)) THEN
          WRITE(OPUNIT,1000) RETADR, DELADR
          GOTO 99
        ENDIF
      ELSE IF (STATUS.EQ.2) THEN
        ISTAT = LIB$FREE_VM(NBYTES, START)
        IF (.NOT.ISTAT) THEN
          CALL LIB$SIGNAL(%VAL(ISTAT))
          GOTO 99
        ENDIF
      ELSE
        GOTO 99
      ENDIF

      M2MEMD = .TRUE.
@ELSEIF HP
      CALL FREE_(START)
      M2MEMD = .TRUE.
@ELSEIF SGI | LINUX64 | LINUX32
      NBYTES = NBYTES
      STATUS = STATUS
      CALL C2MEMD(START)
      M2MEMD = .TRUE.
@ELSE
      M2MEMD = .FALSE.
@ENDIF

      RETURN

@IF VMS
C Return FALSE if the deallocation failed

99    M2MEMD = .FALSE.
      RETURN

1000  FORMAT('M2MEMD: *ERROR* Please tell Robert: ',2Z12,' vs. ',
     +      2Z12)
@ENDIF
      END
