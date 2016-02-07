C***********************************************************************
C Execute an operating system command - obviously a machine dependent
C bit of code! If the string is being passed on to C then it has to be
C '\0' terminated.
C RME 7/7/93 - 16/7/93
C***********************************************************************

      SUBROUTINE M2SHEL(COMAND)
      CHARACTER*(*) COMAND

      INCLUDE 'm2.inc'

      INTEGER ISTAT, ILEN
      INTEGER C2SYSTEM

      ILEN = LEN(COMAND)
      IF (COMAND(ILEN:ILEN).NE.' ') THEN
        WRITE(OUTP,1005)
      ELSE
10      ILEN = ILEN - 1
        IF (COMAND(ILEN:ILEN).EQ.' ') GOTO 10
        COMAND(ILEN+1:ILEN+1) = CHAR(0)
        ISTAT = C2SYSTEM(COMAND)
        IF (ISTAT.EQ.-1) WRITE(OUTP,1000)
      ENDIF

      RETURN

1000  FORMAT('M2SHEL: Error executing shell command')
1005  FORMAT('M2SHEL: Shell command too long to terminate')
      END

C***********************************************************************
C Replace the system EXIT call with something crude but safe
C RME 19-4-2010
C***********************************************************************

      SUBROUTINE M2EXIT

      STOP "Program terminated"
      RETURN
      END

C***********************************************************************
C Get the string associated with a symbol (VMS) or with an environment
C variable (UNIX). If string is not found ILEN is returned as 0.
C RME 13/1/94 - 13/1/94
C***********************************************************************

      SUBROUTINE M2SYMB(VNAME, WORD, ILEN)
      CHARACTER*(*) VNAME
      CHARACTER*(*) WORD
      INTEGER ILEN

      INTEGER IS, IE

C Get the length of the variable name

      CALL P2SLEN(VNAME, IS, IE)
      WORD = ' '
      ILEN = 0

      RETURN
      END

C***********************************************************************
C Initialise the system-specific variables and start off the timer. The
C system-specific things cannot change during a run - eg: machine name,
C type and username.
C RME 14/7/93 - 19/11/93
C***********************************************************************

      SUBROUTINE M2INIT

      INCLUDE 'm2.inc'
      INCLUDE 'p2.inc'


C Get the machine type, detailed machine name, OS type, detailed OS version
C node name and username

      MACH = 'Linux Box'
      LMACH = 9

      CALL C2MDAT(OSTYP,LOSTYP,OSVER,LOSVER,MCNAM,LMCNAM,MCTYP,LMCTYP)
      IF (LOSTYP.EQ.0) THEN
        WRITE(OUTP,1000) 'OS type'
        WRITE(OUTP,1000) 'OS version'
        WRITE(OUTP,1000) 'host name'
        WRITE(OUTP,1000) 'machine type'
        OSTYP = 'unknown OS'
        LOSTYP = 10
        OSVER = 'unknown OS version'
        LOSVER = 18
        MCNAM = 'unknown host'
        LMCNAM = 12
        MCTYP = 'unknown hardware'
        LMCTYP = 16
      ELSE
        IF (LOSTYP.LT.LEN(OSTYP)) OSTYP(LOSTYP+1:) = ' '
        IF (LOSVER.LT.LEN(OSVER)) OSVER(LOSVER+1:) = ' '
        IF (LMCNAM.LT.LEN(MCNAM)) MCNAM(LMCNAM+1:) = ' '
        IF (LMCTYP.LT.LEN(MCTYP)) MCTYP(LMCTYP+1:) = ' '
      ENDIF

      CALL C2USER(USERN, LUSERN)
      IF (LUSERN.EQ.0) THEN
        WRITE(OUTP,1000) 'user name'
        USERN = 'unknown user'
        LUSERN = 12
      ELSE
        IF (LUSERN.LT.LEN(USERN)) USERN(LUSERN+1:) = ' '
      ENDIF

      RETURN

1000  FORMAT('M2INIT: Error trying to find out ',A)
      END

C***********************************************************************
C Return the amount of CPU time used since the initialisation call in
C 100ths of a second.
C RME 14/7/93 - 14/7/93
C***********************************************************************

      INTEGER FUNCTION M2CPU()

      REAL COMBT, ETIME, TARRAY(2)

C Time is returned real in seconds

      COMBT = ETIME(TARRAY)
      M2CPU = INT(COMBT * 100.0)

      RETURN
      END

C***********************************************************************
C Return the time in the format 12:34:56. It assumes that the passed
C string is 8 characters, 'cos it should only be called from P2TIME,
C which I know does this
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE M2TIME(STRING)
      CHARACTER*(*) STRING

      INTEGER HR, MN, SC
      EXTERNAL C2TIME

      CALL C2TIME(HR, MN, SC)
      WRITE(STRING,1000) HR, MN, SC
      IF (STRING(1:1).EQ.' ') STRING(1:1) = '0'
      IF (STRING(4:4).EQ.' ') STRING(4:4) = '0'
      IF (STRING(7:7).EQ.' ') STRING(7:7) = '0'
1000  FORMAT(I2,':',I2,':',I2)

      RETURN
      END

C***********************************************************************
C Return the date as a series of three integers. This is then used by
C the P2 date routines to construct formatted date strings. The year is
C returned as a 4 digit number - whether the value will be right come
C year 2000 I don't know...!
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE M2IDAT(DAY, MONTH, YEAR)
      INTEGER DAY, MONTH, YEAR

      CALL C2IDATE(DAY, MONTH, YEAR)

      RETURN
      END

