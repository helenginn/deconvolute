C***********************************************************************
C Return the address at which the specified item can be found. Used for
C the memory allocation routines. (Non-character variables)
C RME 24/7/93 - 24/7/93
C***********************************************************************

      INTEGER (KIND=8) FUNCTION M2VLOC(ITEM)
      INTEGER ITEM

      M2VLOC = LOC(ITEM)

      RETURN
      END

C***********************************************************************
C Return the address at which the specified item can be found. Used for
C the memory allocation routines. Character variables ONLY!
C RME 24/7/93 - 24/7/93
C***********************************************************************

      INTEGER (KIND=8) FUNCTION M2CLOC(ITEM)
      CHARACTER*(*) ITEM

      M2CLOC = LOC(ITEM)

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
      INTEGER (KIND=8) NBYTES, START
      INTEGER STATUS

      INCLUDE 'p2.inc'

      EXTERNAL C2MEMA

C On the SG we use the an interface to the C malloc routine

      STATUS = 1
      M2MEMA = .TRUE.
      CALL C2MEMA(NBYTES, START)
      IF (START.EQ.0) GOTO 99

      RETURN

C Return FALSE if the allocation failed with NBYTES and START = 0

99    NBYTES = 0
      START = 0
      M2MEMA = .FALSE.
      RETURN
      END

C***********************************************************************
C Routine to deallocate a block of memory of at least NBYTES bytes which
C start from START. The variable STATUS was the value returned by the
C allocation call, and can be used for any machine dependent purpose. If
C the deallocate fails then return FALSE.
C RME 24/7/93 - 24/7/93
C***********************************************************************

      LOGICAL FUNCTION M2MEMD(NBYTES, START, STATUS)
      INTEGER (KIND=8) NBYTES, START
      INTEGER STATUS

      INCLUDE 'p2.inc'

      EXTERNAL C2MEMD

      NBYTES = NBYTES
      STATUS = STATUS
      CALL C2MEMD(START)
      M2MEMD = .TRUE.

      RETURN

      END
