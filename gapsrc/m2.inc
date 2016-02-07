C***********************************************************************
C Parameters used by the M2 routines, these are machine-dependent
C routines, but this include file is machine independent.
C
C Robert Esnouf, 4/7/93 - 21/11/93
C
C***********************************************************************

      INTEGER MNUNIT, MXUNIT, INP, OUTP

C Minimum unit number available on machine

      PARAMETER (MNUNIT=1)

C Maximum unit number available on machine

      PARAMETER (MXUNIT=99)

C Unit used for standard input (from a terminal etc)

      PARAMETER (INP=5)

C Unit used for standard input (from a terminal etc)

      PARAMETER (OUTP=6)

