C***********************************************************************
C Parameters used by the G2 routines, these are used to fool GETWORD14
C into running.
C
C Robert Esnouf, 4/7/93 - 4/7/93
C
C***********************************************************************

      INTEGER G2MINU, G2MAXU, NOUNIT

C Minimum possible unit number for GETWORD to be using

      PARAMETER (G2MINU=1)

C Maximum possible unit number for GETWORD to be using

      PARAMETER (G2MAXU=99)

C Value to show that a GETWORD unit is not currently attached to a real
C FORTRAN unit

      PARAMETER (NOUNIT=-1)

C***********************************************************************
C Common blocks used by the G2 routines, inorder to fool GETWORD14 into
C running.
C
C Robert Esnouf, 4/7/93 - 4/7/93
C
C***********************************************************************

      COMMON /G2/ G2UNIT
      SAVE /G2/

C List mapping GETWORD unit numbers to real FORTRAN unit numbers

      INTEGER G2UNIT(G2MINU:G2MAXU)

