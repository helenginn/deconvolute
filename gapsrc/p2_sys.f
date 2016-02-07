C***********************************************************************
C System calls for the new parser P2. These largely just interface to
C M2 calls and print the answer or return a string to the caller
C
C Robert Esnouf, 4/7/93 - 12/11/93
C
C***********************************************************************

C***********************************************************************
C Produce a quick simple timestamp for a program. It gives machine name
C username date and CPU time. Sophisticated programs will probably do
C their own thing using the more basic calls.
C
C Robert Esnouf, 4/7/93 - 12/11/93
C
C***********************************************************************

      SUBROUTINE P2STMP

      INCLUDE 'p2.inc'

      INTEGER ICENT, M2CPU, ILEN
      CHARACTER*28 CPUTIM
      CHARACTER*19 DATE
      CHARACTER*8 TIME

      WRITE(OPUNIT,1000)
      WRITE(OPUNIT,1005)

C Write the details of the machine

      WRITE(OPUNIT,1010) USERN(:LUSERN), MCNAM(:LMCNAM)
      WRITE(OPUNIT,1020) MCTYP(:LMCTYP)
      WRITE(OPUNIT,1030) OSTYP(:LOSTYP), OSVER(:LOSVER)
      WRITE(OPUNIT,1005)

C Write the timing details including CPU usage but don't bother if usage
C is less than 0.1 secs, since this is probably a startup call

      ICENT = M2CPU()

      IF (ICENT.GT.10) THEN
        WRITE(OPUNIT,1040) STTIM(:LSTTIM), STDAT(:LSTDAT)
      ENDIF

      CALL P2TIME(TIME, ILEN)
      CALL P2DATE(DATE, ILEN)
      WRITE(OPUNIT,1050) TIME, DATE(:ILEN)

      IF (ICENT.GT.10) THEN
        CALL P2FTIM(ICENT, CPUTIM, ILEN)
        WRITE(OPUNIT,1060) CPUTIM(:ILEN)
      ENDIF
      WRITE(OPUNIT,1005)

      WRITE(OPUNIT,1000)
      RETURN

1000  FORMAT('*****************************************************')
1005  FORMAT(' ')
1010  FORMAT('   Program run for user ',A,' on ',A)
1020  FORMAT('   Type of machine: ',A)
1030  FORMAT('   Operating system: ',A,' ',A)
1040  FORMAT('   Job was started ',A,' on ',A)
1050  FORMAT('   Current time is ',A,' on ',A)
1060  FORMAT('   Total CPU time used: ',A)
      END

C***********************************************************************
C Return the current amount of CPU time used and the amount of CPU time
C used in the last call to this or P2PCPU or since the program started
C RME 20/7/93 - 20/7/93
C***********************************************************************

      SUBROUTINE P2CPU(ICURR, ISINCE)
      INTEGER ICURR
      INTEGER ISINCE

      INCLUDE 'p2.inc'

      INTEGER M2CPU

      ICURR = M2CPU()
      ISINCE = ICURR - OLDCPU
      OLDCPU = ICURR

      RETURN
      END

C***********************************************************************
C Calculate the CPU time used and the time since the previous call to
C this routine, since calling P2CPU or since the program started.
C RME 20/7/93 - 20/7/93
C***********************************************************************

      SUBROUTINE P2PCPU

      INCLUDE 'p2.inc'

      INTEGER ICENT, M2CPU, ILEN
      CHARACTER*28 CPUTIM

C Now construct a string of total CPU usage

      ICENT = M2CPU()
      CALL P2FTIM(ICENT, CPUTIM, ILEN)
      WRITE(OPUNIT,1000)
      WRITE(OPUNIT,1010)
      WRITE(OPUNIT,1020) CPUTIM(:ILEN)

C Then a string of change in CPU time

      CALL P2FTIM(ICENT-OLDCPU, CPUTIM, ILEN)
      WRITE(OPUNIT,1030) CPUTIM(:ILEN)
      WRITE(OPUNIT,1010)
      WRITE(OPUNIT,1000)
      OLDCPU = ICENT

      RETURN

1000  FORMAT(' ')
1010  FORMAT('*********************************************')
1020  FORMAT('   CPU usage currently: ',A)
1030  FORMAT('   CPU since last call: ',A)
      END

C***********************************************************************
C Return the actual time in the format 12:34:56. All returned strings
C are thus 8 characters long so ILEN = 8 or the length of the input
C string. This is just an interface to the M2 call to call the operating
C system for us.
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE P2TIME(STRING, ILEN)
      CHARACTER*(*) STRING
      INTEGER ILEN

      CHARACTER*8 TSTR

      CALL M2TIME(TSTR)
      ILEN = MIN(LEN(STRING), 8)
      STRING = TSTR(:ILEN)

      RETURN
      END

C***********************************************************************
C Return the date in the long format. The maximum length of the date
C string is 19 characters. An example of the format '14th July 1993'.
C Note that 5 short of ILEN just gives the day and month. Uses the
C M2IDAT routine to find the day month and year.
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE P2DATE(STRING, ILEN)
      CHARACTER*(*) STRING
      INTEGER ILEN

      CHARACTER*19 DSTR
      CHARACTER*9 MONTH(12)
      INTEGER IDAY, IMON, IYEAR

      DATA MONTH /'January  ','February ','March    ','April    ',
     +            'May      ','June     ','July     ','August   ',
     +            'September','October  ','November ','December '/

      CALL M2IDAT(IDAY, IMON, IYEAR)
      DSTR = ' '
      WRITE(DSTR(1:2), '(I2)') IDAY
      IF (IDAY.EQ.1 .OR. IDAY.EQ.21 .OR. IDAY.EQ.31) THEN
        DSTR(3:4) = 'st'
      ELSE IF (IDAY.EQ.2 .OR. IDAY.EQ.22) THEN
        DSTR(3:4) = 'nd'
      ELSE IF (IDAY.EQ.3 .OR. IDAY.EQ.23) THEN
        DSTR(3:4) = 'rd'
      ELSE
        DSTR(3:4) = 'th'
      ENDIF
      DSTR(6:14) = MONTH(IMON)
      WRITE(DSTR(16:19), '(I4)') IYEAR
      CALL P2SSET(DSTR, ILEN)

      ILEN = MIN(LEN(STRING), ILEN)
      STRING = DSTR(:ILEN)

      RETURN
      END

C***********************************************************************
C Return the date in the form 14/7/93. The maximum length the string can
C be is 8. The length is returned in ILEN. This is, of course, the
C British way of giving the date. Uses the M2IDAT function to get each
C bit of the date.
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE P2DAT(STRING, ILEN)
      CHARACTER*(*) STRING
      INTEGER ILEN

      CHARACTER*8 DSTR
      INTEGER IS, IE
      INTEGER IDAY, IMON, IYEAR

      CALL M2IDAT(IDAY, IMON, IYEAR)
      IF (IMON.LT.10) THEN
        WRITE(DSTR,1000) IDAY, IMON, MOD(IYEAR,100)
      ELSE
        WRITE(DSTR,1010) IDAY, IMON, MOD(IYEAR,100)
      ENDIF

      CALL P2SLEN(DSTR, IS, IE)
      IF (IS.GT.1) DSTR = DSTR(IS:)

      ILEN = MIN(LEN(STRING), IE+IS-1)
      STRING = DSTR(:ILEN)

      RETURN

1000  FORMAT(I2,'/',I1,'/',I2,' ')
1010  FORMAT(I2,'/',I2,'/',I2)
      END

C***********************************************************************
C Return the simple machine type string. The actual length of the string
C is returned in ILEN, but is truncated if the string is too short.
C This function which should be used to determine if machine
C specific code needs to be skipped.
C
C RME 16/7/93 - 16/7/93
C***********************************************************************

      SUBROUTINE P2MACH(STRING, ILEN)
      CHARACTER*(*) STRING
      INTEGER ILEN

      INCLUDE 'p2.inc'

      ILEN = MIN(LEN(STRING),LMACH)
      STRING = MACH(:ILEN)

      RETURN
      END

C***********************************************************************
C Return the machine type string. The actual length of the string is
C returned in ILEN, but is truncated if the string is too short.
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE P2SYST(STRING, ILEN)
      CHARACTER*(*) STRING
      INTEGER ILEN

      INCLUDE 'p2.inc'

      ILEN = MIN(LEN(STRING),LMCTYP)
      STRING = MCTYP(:ILEN)

      RETURN
      END

C***********************************************************************
C Return the machine name string. The actual length of the string is
C returned in ILEN, but is truncated if the string is too short.
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE P2SYSN(STRING, ILEN)
      CHARACTER*(*) STRING
      INTEGER ILEN

      INCLUDE 'p2.inc'

      ILEN = MIN(LEN(STRING),LMCNAM)
      STRING = MCNAM(:ILEN)

      RETURN
      END

C***********************************************************************
C Return the operating system string. The actual length of the string is
C returned in ILEN, but is truncated if the string is too short. The
C possible types are 'VMS' and 'UNIX'
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE P2SYSO(STRING, ILEN)
      CHARACTER*(*) STRING
      INTEGER ILEN

      INCLUDE 'p2.inc'

      ILEN = MIN(LEN(STRING),LOSTYP)
      STRING = OSTYP(:ILEN)

      RETURN
      END

C***********************************************************************
C Return the OS version string. The actual length of the string is
C returned in ILEN, but is truncated if the string is too short.
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE P2SYSV(STRING, ILEN)
      CHARACTER*(*) STRING
      INTEGER ILEN

      INCLUDE 'p2.inc'

      ILEN = MIN(LEN(STRING),LOSVER)
      STRING = OSVER(:ILEN)

      RETURN
      END

C***********************************************************************
C Return the username string. The actual length of the string is
C returned in ILEN, but is truncated if the string is too short.
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE P2SYSU(STRING, ILEN)
      CHARACTER*(*) STRING
      INTEGER ILEN

      INCLUDE 'p2.inc'

      ILEN = MIN(LEN(STRING),LUSERN)
      STRING = USERN(:ILEN)

      RETURN
      END

C***********************************************************************
C Execute an operating system command. This call is simply a front end
C to the machine-dependent routine which does the job.
C RME 17/7/93 - 17/7/93
C***********************************************************************

      SUBROUTINE P2SHEL(COMAND)
      CHARACTER*(*) COMAND

      CALL M2SHEL(COMAND)

      RETURN
      END

