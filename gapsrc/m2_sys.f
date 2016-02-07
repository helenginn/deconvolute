C***********************************************************************
C Execute an operating system command - obviously a machine dependent
C bit of code! If the string is being passed on to C then it has to be
C '\0' terminated.
C RME 7/7/93 - 16/7/93
C***********************************************************************

      SUBROUTINE M2SHEL(COMAND)
      CHARACTER*(*) COMAND

@IF !VMS
      INCLUDE 'm2.inc'

@ENDIF
@IF VMS
      INTEGER ISTAT
      LOGICAL P2ICMP
      INTEGER LIB$SPAWN
@ELSEIF CONVEX | HP | SGI
      INTEGER ISTAT
      INTEGER SYSTEM
@ELSEIF LINUX64 | LINUX32
      INTEGER ISTAT, ILEN
      INTEGER C2SYSTEM
@ENDIF

@IF VMS
C If command is '#SPAWN' then don't spawn two shells

      IF (P2ICMP(COMAND,'SPAWN')) THEN
        ISTAT = LIB$SPAWN()
      ELSE
        ISTAT = LIB$SPAWN(COMAND)
      ENDIF

      IF (.NOT.ISTAT) THEN
        CALL LIB$SIGNAL(%VAL(ISTAT))
      ENDIF
@ELSEIF CONVEX | HP | SGI
      ISTAT = SYSTEM(COMAND)
      IF (ISTAT.EQ.127) THEN
        WRITE(OUTP,1000)
      ENDIF
@ELSEIF LINUX64 | LINUX32
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
@ENDIF

      RETURN
@IF CONVEX | HP | SGI | LINUX64 | LINUX32

1000  FORMAT('M2SHEL: Error executing shell command')
@ENDIF
@IF LINUX64 | LINUX32
1005  FORMAT('M2SHEL: Shell command too long to terminate')
@ENDIF
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
@IF VMS
      INTEGER ISTAT, LIB$GET_SYMBOL
      INCLUDE '($LIBDEF)'
@ENDIF

C Get the length of the variable name

      CALL P2SLEN(VNAME, IS, IE)
@IF VMS
      ISTAT = LIB$GET_SYMBOL(VNAME(IS:IE), WORD, ILEN)
      IF (.NOT.ISTAT) THEN
        IF (ISTAT.EQ.LIB$_NOSUCHSYM .OR. ISTAT.EQ.LIB$_INVSYMNAM) THEN
          WORD = ' '
          ILEN = 0
        ELSE
          CALL LIB$SIGNAL(%VAL(ISTAT))
        ENDIF
      ENDIF
@ELSE
      WORD = ' '
      ILEN = 0
@ENDIF

      RETURN
      END

C***********************************************************************
C Initialise the system-specific variables and start off the timer. The
C system-specific things cannot change during a run - eg: machine name,
C type and username.
C RME 14/7/93 - 19/11/93
C***********************************************************************

      SUBROUTINE M2INIT

@IF VMS
      INCLUDE '($SYIDEF)'
      INCLUDE '($JPIDEF)'
      INCLUDE 'm2_vms.inc'
@ENDIF
      INCLUDE 'm2.inc'
      INCLUDE 'p2.inc'

@IF VMS
      INTEGER ISTAT, IS
      INTEGER LIB$INIT_TIMER, LIB$GETSYI, LIB$GETJPI
@ELSEIF CONVEX
      INTEGER GETLOG, GETUID, HOSTNM
@ELSEIF HP
      INTEGER HOSTNM_
@ELSEIF HP
      INTEGER I
@ENDIF

@IF VMS
C Initialise the timer if necessary (not on Unix)

      ISTAT = LIB$INIT_TIMER()
      IF (.NOT.ISTAT) THEN
        CALL LIB$SIGNAL(%VAL(ISTAT))
        WRITE(OUTP,1010)
      ENDIF      

@ENDIF
C Get the machine type, detailed machine name, OS type, detailed OS version
C node name and username

@IF VMS
      ISTAT = LIB$GETSYI(SYI$_CPU, IS)
      IF (.NOT.ISTAT) THEN
        CALL LIB$SIGNAL(%VAL(ISTAT))
        WRITE(OUTP,1000) 'machine'
        MACH = 'DEC UNKNOWN'
        LMACH = 11
        MEMLIM = 8192
      ELSE IF (IS.LT.128) THEN
        MACH = 'DEC VAX'
        LMACH = 7
        MEMLIM = 512
      ELSE
        MACH = 'DEC ALPHA'
        LMACH = 9
        MEMLIM = 8192
      ENDIF

      OSTYP = 'VMS'
      LOSTYP = 3

      ISTAT = LIB$GETSYI(SYI$_HW_NAME,,MCTYP,LMCTYP)
      IF (.NOT.ISTAT) THEN
        CALL LIB$SIGNAL(%VAL(ISTAT))
        WRITE(OUTP,1000) 'machine type'
        MCTYP = 'Unknown VMS machine'
        LMCTYP = 19
      ENDIF      
      CALL P2SLEN(MCTYP, IS, LMCTYP)
      IF (IS.GT.1) THEN
        MCTYP = MCTYP(IS:)
        LMCTYP = LMCTYP - IS + 1
      ENDIF

      ISTAT = LIB$GETSYI(SYI$_VERSION,,OSVER,LOSVER)
      IF (.NOT.ISTAT) THEN
        CALL LIB$SIGNAL(%VAL(ISTAT))
        WRITE(OUTP,1000) 'OS version'
        OSVER = 'Unknown OS version'
        LOSVER = 18
      ENDIF      
      CALL P2SLEN(OSVER, IS, LOSVER)
      IF (IS.GT.1) THEN
        OSVER = OSVER(IS:)
        LOSVER = LOSVER - IS + 1
      ENDIF

      ISTAT = LIB$GETSYI(SYI$_NODENAME,,MCNAM,LMCNAM)
      IF (.NOT.ISTAT) THEN
        CALL LIB$SIGNAL(%VAL(ISTAT))
        WRITE(OUTP,1000) 'node name'
        MCNAM = 'UNKNOWN'
        LMCNAM = 7
      ENDIF      
      CALL P2SLEN(MCNAM, IS, LMCNAM)
      IF (IS.GT.1) THEN
        MCNAM = MCNAM(IS:)
        LMCNAM = LMCNAM - IS + 1
      ENDIF

      ISTAT = LIB$GETJPI(JPI$_USERNAME,,,,USERN,LUSERN)
      IF (.NOT.ISTAT) THEN
        CALL LIB$SIGNAL(%VAL(ISTAT))
        WRITE(OUTP,1000) 'user name'
        USERN = 'ANONYMOUS'
        LUSERN = 9
      ENDIF      
      CALL P2SLEN(USERN, IS, LUSERN)
      IF (IS.GT.1) THEN
        USERN = USERN(IS:)
        LUSERN = LUSERN - IS + 1
      ENDIF

@ELSEIF CONVEX
      MACH = 'CONVEX C-SERIES'
      LMACH = 15

      OSTYP = 'UNIX'
      LOSTYP = 4

      WRITE(OUTP,1000) 'machine type'
      MCTYP = 'Convex C2xx'
      LMCTYP = 11

      WRITE(OUTP,1000) 'OS version'
      OSVER = 'Unknown OS version'
      LOSVER = 18

      ISTAT = HOSTNM(MCNAM)
      IF (ISTAT.NE.0) THEN
        WRITE(OUTP,1000) 'host name'
        MCNAM = 'unknown'
        LMCNAM = 7
      ENDIF
      CALL P2SLEN(MCNAM, IS, LMCNAM)
      IF (IS.GT.1) THEN
        MCNAM = MCNAM(IS:)
        LMCNAM = LMCNAM - IS + 1
      ENDIF

      CALL GETLOG(USERN)
      IF (USERN.EQ.' ') THEN
        IS = GETUID()
        WRITE(USERN,'(I10)') IS
      ENDIF
      CALL P2SLEN(USERN, IS, LUSERN)
      IF (IS.GT.1) THEN
        USERN = USERN(IS:)
        LUSERN = LUSERN - IS + 1
      ENDIF

@ELSEIF HP
      MACH = 'HP PA-RISC'
      LMACH = 10

      OSTYP = 'UNIX'
      LOSTYP = 4

      CALL C2MDAT(OSVER, LOSVER, MCTYP, LMCTYP)
      IF (LOSVER.EQ.0) THEN
        WRITE(OUTP,1000) 'machine type'
        WRITE(OUTP,1000) 'OS version'
        MCTYP = 'Unknown HP machine'
        LMCTYP = 18
        OSVER = 'Unknown OS version'
        LOSVER = 18
      ELSE
        IF (LOSVER.LT.LEN(OSVER)) OSVER(LOSVER+1:) = ' '
        IF (LMCTYP.LT.LEN(MCTYP)) MCTYP(LMCTYP+1:) = ' '
      ENDIF

      ISTAT = HOSTNM_(MCNAM)
      IF (ISTAT.NE.0) THEN
        WRITE(OUTP,1000) 'host name'
        MCNAM = 'unknown'
        LMCNAM = 7
      ENDIF
      CALL P2SLEN(MCNAM, IS, LMCNAM)
      IF (IS.GT.1) THEN
        MCNAM = MCNAM(IS:)
        LMCNAM = LMCNAM - IS + 1
      ENDIF

      CALL C2USER(USERN, LUSERN)
      IF (LUSERN.EQ.0) THEN
        WRITE(OUTP,1000) 'user name'
        USERN = 'anonymous'
        LUSERN = 9
      ELSE
        IF (LUSERN.LT.LEN(USERN)) USERN(LUSERN+1:) = ' '
      ENDIF

@ELSEIF SGI
      MACH = 'SGI MIPS-RISC'
      LMACH = 13

      OSTYP = 'UNIX'
      LOSTYP = 4

      CALL C2MDAT(OSVER, LOSVER, MCNAM, LMCNAM, MCTYP, LMCTYP)
      IF (LOSVER.EQ.0) THEN
        WRITE(OUTP,1000) 'machine type'
        WRITE(OUTP,1000) 'host name'
        WRITE(OUTP,1000) 'OS version'
        MCTYP = 'Unknown SGI machine'
        LMCTYP = 18
        USERN = 'unknown host'
        LUSERN = 12
        OSVER = 'Unknown OS version'
        LOSVER = 18
      ELSE
        IF (LOSVER.LT.LEN(OSVER)) OSVER(LOSVER+1:) = ' '
        IF (LMCNAM.LT.LEN(MCNAM)) MCNAM(LMCNAM+1:) = ' '
        IF (LMCTYP.LT.LEN(MCTYP)) MCTYP(LMCTYP+1:) = ' '
      ENDIF

      CALL C2USER(USERN, LUSERN)
      IF (LUSERN.EQ.0) THEN
        WRITE(OUTP,1000) 'user name'
        USERN = 'anonymous'
        LUSERN = 9
      ELSE
        IF (LUSERN.LT.LEN(USERN)) USERN(LUSERN+1:) = ' '
      ENDIF

@ELSEIF LINUX64 | LINUX32
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

@ELSE
      MACH = 'UNSUPPORTED'
      LMACH = 11

      WRITE(OUTP,1000) 'OS type'
      OSTYP = 'UNKNOWN'
      LOSTYP = 7

      WRITE(OUTP,1000) 'machine type'
      MCTYP = 'Unknown machine'
      LMCTYP = 15

      WRITE(OUTP,1000) 'OS version'
      OSVER = 'Unknown OS version'
      LOSVER = 18

      WRITE(OUTP,1000) 'host name'
      MCNAM = 'unknown host'
      LMCNAM = 12

      WRITE(OUTP,1000) 'user name'
      USERN = 'anonymous'
      LUSERN = 9

@ENDIF
      RETURN

1000  FORMAT('M2INIT: Error trying to find out ',A)
@IF VMS
1010  FORMAT('M2INIT: Error trying to initialise timer')
@ENDIF
      END

C***********************************************************************
C Return the amount of CPU time used since the initialisation call in
C 100ths of a second.
C RME 14/7/93 - 14/7/93
C***********************************************************************

      INTEGER FUNCTION M2CPU()

@IF VMS
      INTEGER ICODE, LIB$STAT_TIMER, ISTAT
@ELSEIF CONVEX | SGI | LINUX64 | LINUX32
      REAL COMBT, ETIME, TARRAY(2)
@ELSEIF HP
      REAL COMBT, ETIME_, TARRAY(2)
@ENDIF

@IF VMS
C Just use the LIB$ call for the result

      ICODE = 2
      ISTAT = LIB$STAT_TIMER(ICODE, M2CPU)
      IF (.NOT.ISTAT) THEN
        CALL LIB$SIGNAL(%VAL(ISTAT))
        M2CPU = 0
      ENDIF      
@ELSEIF CONVEX
C Time is returned real in seconds to an accuracy of 1/60th second

      COMBT = ETIME(TARRAY)      
      M2CPU = INT(COMBT * 100.0)
@ELSEIF HP
C Time is returned real in seconds to the accuracy of the architecture
C this can be measured using sysconf and for HP735 was 0.01s

      COMBT = ETIME_(TARRAY)      
      M2CPU = INT(COMBT * 100.0)
@ELSEIF SGI | LINUX64 | LINUX32
C Time is returned real in seconds

      COMBT = ETIME(TARRAY)      
      M2CPU = INT(COMBT * 100.0)
@ELSE
C Cannot determine a time for unknown machine

      M2CPU = 0
@ENDIF

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

@IF VMS | CONVEX | HP | SGI 
      CALL TIME(STRING)
@ELSEIF LINUX64 | LINUX32
      INTEGER HR, MN, SC
      EXTERNAL C2TIME

      CALL C2TIME(HR, MN, SC)
      WRITE(STRING,1000) HR, MN, SC
      IF (STRING(1:1).EQ.' ') STRING(1:1) = '0'
      IF (STRING(4:4).EQ.' ') STRING(4:4) = '0'
      IF (STRING(7:7).EQ.' ') STRING(7:7) = '0'
1000  FORMAT(I2,':',I2,':',I2)
@ELSE
      STRING = '00:00:00'
@ENDIF

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

@IF CONVEX
      INTEGER IARRAY(3)

@ELSEIF HP
      EXTERNAL IDATE_
      INTEGER IARRAY(3)

@ENDIF
@IF VMS | SGI
      CALL IDATE(MONTH, DAY, YEAR)
      IF (YEAR.LT.90) THEN
        YEAR = 2000 + YEAR
      ELSE
        YEAR = 1900 + YEAR
      ENDIF
@ELSEIF CONVEX
      CALL IDATE(IARRAY)
      DAY = IARRAY(1)
      MONTH = IARRAY(2)
      YEAR = IARRAY(3)
@ELSEIF HP
      CALL IDATE_(IARRAY)
      DAY = IARRAY(1)
      MONTH = IARRAY(2)
      YEAR = IARRAY(3)
@ELSEIF LINUX64 | LINUX32
      CALL C2IDATE(DAY, MONTH, YEAR)
@ELSE
      DAY = 1
      MONTH = 1
      YEAR = 1900
@ENDIF

      RETURN
      END

