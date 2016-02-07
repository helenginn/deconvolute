C***********************************************************************
C This function opens the standard I/O units and sets the output to be
C line buffered if necessary. It also checks whether the standard input
C is a terminal or not. If it is then at the bottom level we prompt for
C new lines of data and don't echo them.
C RME 7/7/93 - 13/7/93
C***********************************************************************

      LOGICAL FUNCTION M2STD()

      INCLUDE 'm2.inc'
      INCLUDE 'p2.inc'

      LOGICAL ISATTY

C Set the machine independent input and output unit variables

      IPUNIT = INP
      OPUNIT = OUTP

C Work out the range of FORTRAN units available, and check there
C are some

      UNMIN = MAX(PUNMIN, MNUNIT)
      UNMAX = MIN(PUNMAX, MXUNIT)
      IF (UNMIN.GT.UNMAX) GOTO 99

C Open the standard output file and make it line-buffered. Usually the
C output unit is opened properly on entry. The line buffering is made
C using a call to a C routine

      CALL C2LBUF

C See where input on stdin is actually coming from

      ISTERM = ISATTY(INP)

      M2STD = .TRUE.
      RETURN

C There was an error performing the standard unit initialisation

99    M2STD = .FALSE.
      RETURN

      END

C***********************************************************************
C Open a file on a specified FORTRAN unit. The filename and characteris-
C tics can be found in the p2.inc common blocks.
C RME 7/7/93 - 21/11/93
C***********************************************************************

      LOGICAL FUNCTION M2OPEN(UNIT)
      INTEGER UNIT

      INCLUDE 'm2.inc'
      INCLUDE 'p2.inc'

      INTEGER RL
      CHARACTER*7 STAT
      CHARACTER*10 ACC
      CHARACTER*11 FORM

      STAT = 'UNKNOWN'
      ACC = 'SEQUENTIAL'
      FORM = 'FORMATTED'

      IF (FUNIT(UNIT)(3:3).EQ.'R') THEN
        STAT='OLD'
      ELSE IF (FUNIT(UNIT)(4:4).EQ.'W') THEN
        IF (FUNIT(UNIT)(7:7).EQ.'?') THEN
          STAT = 'UNKNOWN'
        ELSE
          STAT = 'NEW'
        ENDIF
      ELSE IF (FUNIT(UNIT)(5:5).EQ.'A') THEN
        STAT='OLD'
        ACC='APPEND'
      ELSE IF (FUNIT(UNIT)(7:7).EQ.'?') THEN
        STAT = 'UNKNOWN'
      ELSE IF (FUNIT(UNIT)(6:6).EQ.'S') THEN
        STAT='SCRATCH'
      ENDIF

      IF (FUNIT(UNIT)(8:8).EQ.'F') THEN
        FORM = 'FORMATTED'
      ELSE IF (FUNIT(UNIT)(9:9).EQ.'U') THEN
        FORM = 'UNFORMATTED'
      ENDIF

      IF (FUNIT(UNIT)(10:10).EQ.'D') THEN
        RL = FRECL(UNIT)
        IF (FUNIT(UNIT)(9:9).EQ.'U') RL = (RL+3)/4
        OPEN(UNIT=UNIT,FILE=FNAME(UNIT),STATUS=STAT,ACCESS='DIRECT',
     +       FORM=FORM,RECL=RL,ERR=99)
      ELSE IF (FUNIT(UNIT)(3:3).EQ.'R') THEN
        OPEN(UNIT=UNIT,FILE=FNAME(UNIT),STATUS=STAT,ACCESS=ACC,
     +       FORM=FORM,ERR=99)
      ELSE
        OPEN(UNIT=UNIT,FILE=FNAME(UNIT),STATUS=STAT,ACCESS=ACC,
     +       FORM=FORM,ERR=99)
      ENDIF

C OPEN went OK so return true, otherwise return false

      M2OPEN = .TRUE.
      RETURN

99    M2OPEN = .FALSE.
      RETURN

      END

C***********************************************************************
C Close a file on the specified unit number, the filename and the file
C characteristics can be found in the P2 common blocks. The CHRS string
C controls keeping, 'K', or deleting, 'D', of files
C RME 8/7/93 - 21/11/93
C***********************************************************************

      LOGICAL FUNCTION M2CLOS(UNIT, CHRS)
      INTEGER UNIT
      CHARACTER*1 CHRS

      CHRS = CHRS
      CLOSE(UNIT, ERR=99)
      M2CLOS = .TRUE.
      RETURN

99    M2CLOS = .FALSE.
      RETURN

      END

C***********************************************************************
C Produce a prompt on the standard output unit. The new line is suppres-
C -sed using the $ format.
C RME 13/7/93 - 13/7/93
C***********************************************************************

      SUBROUTINE M2PROM(PROMPT)
      CHARACTER*(*) PROMPT

      INCLUDE 'm2.inc'
      INCLUDE 'p2.inc'

      WRITE(OUTP,1000,ADVANCE='NO') PROMPT, PARROW(:LARROW)
      RETURN

1000  FORMAT(A,A)
      END

C***********************************************************************
C Read a line from the specified unit. It uses the Q format to find out
C how many characters were read
C RME 9/7/93 - 9/7/93
C***********************************************************************

      SUBROUTINE M2READ(UNIT,LINE,ILEN)
      INTEGER UNIT
      CHARACTER*(*) LINE
      INTEGER ILEN
      INTEGER MAXLEN

      READ(UNIT,'(A)',ERR=99,END=10) LINE
      IF (LINE.EQ.' ') THEN
        ILEN = 0
      ELSE
        MAXLEN  = LEN(LINE)
20      IF (LINE(MAXLEN/2+1:).EQ.' ') THEN
          MAXLEN = MAXLEN/2
          GOTO 20
        ENDIF
        DO ILEN = MAXLEN, 1, -1
          IF (LINE(ILEN:ILEN).NE.' ') GOTO 21
        ENDDO
        ILEN = 0
21      CONTINUE
      ENDIF
      RETURN

C For end of file return ILEN = -1

10    ILEN = -1
      RETURN

C If there is an error reading the file then return -2

99    ILEN = -2
      RETURN

      END

C***********************************************************************
C Check the status of the specified FORTRAN unit, this routine has to
C fill the FNAME and FUNIT arrays for a unit which is not already open
C RME 8/7/93 - 21/11/93
C***********************************************************************

      SUBROUTINE M2FCHK(UNIT)
      INTEGER UNIT

      INCLUDE 'm2.inc'
      INCLUDE 'p2.inc'

      INTEGER RL
      LOGICAL NAMED, OPENED
      CHARACTER*80 NAME
      CHARACTER*12 ACC, FORM

C Assume that if P2 thinks the file is open, and the file is open then
C the information is correct - this may not always be true if random
C CLOSE statements are used, however we can not find out all the info
C we want about an opened file.

      INQUIRE(UNIT=UNIT,ACCESS=ACC,FORM=FORM,NAME=NAME,NAMED=NAMED,
     +        OPENED=OPENED,RECL=RL)
      IF (OPENED .AND. FUNIT(UNIT)(1:1).NE.'-') RETURN

C If the file is not opened then clear its information and return. If
C the unit is opened then try to fill in information about it - this
C will not be complete since the ISINP and ISSCR flags cannot be found.
C Line buffering is not implemented on the VMS version

      IF (.NOT.OPENED) THEN
        FUNIT(UNIT) = '-'
	FRECL(UNIT) = 0
        FNAME(UNIT) = ' '
        RETURN
      ENDIF

C Test if unit is a standard one or not

      IF (UNIT.EQ.OUTP) THEN
        FUNIT(UNIT) = ' * W   F  '
	FRECL(UNIT) = 0
        FNAME(UNIT) = 'Standard FORTRAN output unit'
        RETURN
      ELSE IF (UNIT.EQ.INP) THEN
        FUNIT(UNIT) = ' *R    F  '
	FRECL(UNIT) = 0
        FNAME(UNIT) = 'Standard FORTRAN input unit'
        RETURN
      ENDIF

C Test if file is formatted or not

      IF (FORM.EQ.'FORMATTED') THEN
        FUNIT(UNIT) = '     ? F  '
      ELSE
        FUNIT(UNIT) = '     ?  U '
      ENDIF

C Test for direct access

      IF (ACC.EQ.'DIRECT') THEN
        IF (FORM.EQ.'UNFORMATTED') RL = 4*RL
        FUNIT(UNIT)(10:10) = 'D'
        FRECL(UNIT) = RL
      ELSE
        FRECL(UNIT) = 0
      ENDIF

C Get the filename if it is known

      IF (NAMED) THEN
        FNAME(UNIT) = NAME
      ELSE
        FNAME(UNIT) = 'Cannot determine a filename'
      ENDIF

      RETURN
      END

