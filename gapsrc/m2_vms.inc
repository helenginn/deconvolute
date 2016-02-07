C***********************************************************************
C Parameters and variables used by the M2 routines for the VMS version
C of the parser. This makes the main include file M2.INC machine independent
C
C This is for the VMS version only
C
C Robert Esnouf, 4/7/93 - 1/4/97
C
C***********************************************************************

C Under VMS we use two units with different carriage control. This
C unit is for the FORTRAN formatted carriage control SYS$OUTPUT for
C prompt generation.

      INTEGER OUTP2
      PARAMETER (OUTP2=99)

C For VMS memory allocation we use two allocation routines. The size
C above which we use the system service to do the allocate is governed by
C the integer MEMLIM

      INTEGER MEMLIM
      COMMON /M2_COM/ MEMLIM
      SAVE /M2_COM/
