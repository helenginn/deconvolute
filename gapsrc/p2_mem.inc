C***********************************************************************
C Parameters used by the P2 parser for memory allocation
C
C Robert Esnouf, 21/4/10 - 21/4/10
C
C***********************************************************************

      INTEGER IXSPEC, IXMBAS, IXMSIZ, IXBAS, IXSIZ, IXNUM, IXABAS

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

