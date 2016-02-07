C Program to create the machine specific files for the parser. It
C reads the standard input and writes to the standard output the code for
C the selected target machine. This program should be eventually replaced
C by the C preprocessor.
C
C makevers MACHINE < file.in > file.out
C
C It looks for @IF...@ENDIF combinations and selects the
C correct bits of code. This is not the world's most robust bit
C of code - don't expect it to cope with sloppy syntax!

	INTEGER		NMACS
	PARAMETER	(NMACS=8)

	CHARACTER*8	ARG1, ARG2, MLIST(NMACS)
	CHARACTER*8	MSPEC(NMACS)
	INTEGER		NSTR, I

	INTEGER		command_argument_count
	LOGICAL		DEBUG

	DATA MLIST /'VAX     ', 'ALPHA   ', 'CONVEX  ', 'HP      ','SGI     ',
     +		    'LINUX64 ', 'LINUX32 ', 'UNSUPP  '/

	NSTR = 0
	DEBUG = .FALSE.
	IF (command_argument_count().EQ.1) THEN

	  CALL get_command_argument(1,ARG1)
	  CALL UPCASE(ARG1)
	  DO I = 1, NMACS
	    IF (ARG1.EQ.MLIST(I)) THEN
	      NSTR = 1
	      MSPEC(1) = MLIST(I)
	    ENDIF
	  ENDDO

	ELSE IF (command_argument_count().EQ.2) THEN

	  CALL get_command_argument(1,ARG1)
	  CALL UPCASE(ARG1)
	  CALL get_command_argument(2,ARG2)
	  CALL UPCASE(ARG2)
	  IF (ARG1.EQ.'DEBUG') THEN
	    DEBUG = .TRUE.
	    ARG1 = ARG2
	  ELSE IF (ARG2.EQ.'DEBUG') THEN
	    DEBUG = .TRUE.
	  ELSE
	    ARG1 = '%%%%%%'
	  ENDIF

	  DO I = 1, NMACS
	    IF (ARG1.EQ.MLIST(I)) THEN
	      NSTR = 1
	      MSPEC(1) = MLIST(I)
	    ENDIF
	  ENDDO

	ENDIF

C Syntax was incorrect so report this and stop

	IF (NSTR.EQ.0) THEN
	  WRITE(0,*)
	  WRITE(0,'(A)') 'Syntax is: makevers MACHINE < file.in > file.out'
	  WRITE(0,'(A)') '           where MACHINE is one of "sgi", "hp",'//
     &				     ' "alpha", "linux64" or "linux32"'
	  WRITE(0,*)
	  CALL ERROREXIT
	ENDIF

C Set up related name list

	IF (MSPEC(1).EQ.'VAX' .OR. MSPEC(1).EQ.'ALPHA') THEN
	  NSTR = NSTR + 1
	  MSPEC(NSTR) = 'VMS   '
	ELSE
	  NSTR = NSTR + 1
	  MSPEC(NSTR) = 'UNIX  '
	ENDIF

	CALL PREPFILE(5, 6, 0, NSTR, MSPEC, DEBUG)

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C This is the subroutine that actually parses the input file
C

	SUBROUTINE PREPFILE(UN1, UN2, UERR, NSTR, MSPEC, DEBUG)
	INTEGER UN1, UN2, UERR, NSTR
	CHARACTER*(*) MSPEC(NSTR)
	LOGICAL DEBUG

	INTEGER MXNEST
	PARAMETER (MXNEST=10)

	LOGICAL COPY(0:MXNEST), EVALIF, ELIF(0:MXNEST)
	LOGICAL DELSE(0:MXNEST), WRFLAG
	INTEGER LLEN, NESTLEV, LCIN(0:MXNEST), LCNOW
	CHARACTER*80	LINE, ULINE
	CHARACTER*80    SPC


	SPC = '..'

	COPY(0) = .TRUE.
	LCIN(0) = 0

	WRFLAG = .FALSE.
	LCNOW = 0
	NESTLEV = 0
10	READ(UN1,'(A)',END=99) LINE
	IF (LINE.EQ.' ') GOTO 11
	DO LLEN = LEN(LINE), 1 , -1
	  IF (LINE(LLEN:LLEN).NE.' ') GOTO 12
	ENDDO
11	LLEN = 0
12	LCNOW = LCNOW + 1

C Check for @IFs and deal with them

	IF (LINE(:3).EQ.'@IF') THEN
	  WRFLAG = .FALSE.

C Check we can cope with another level of nesting.

	  IF (NESTLEV.EQ.MXNEST) THEN
	    WRITE(UERR,1000) LCNOW
	    CALL ERROREXIT
	  ENDIF

	  NESTLEV = NESTLEV + 1
	  LCIN(NESTLEV) = LCNOW
	  IF (DEBUG) WRITE(UERR,1100) SPC(:2*NESTLEV), LINE(4:LLEN)

C Condition should not be evaluated if a current @IF was false.
C Otherwise print and evaluate the condition

	  IF (COPY(NESTLEV-1).EQV..FALSE.) THEN
	    COPY(NESTLEV)=.FALSE.
	  ELSE
	    ULINE = LINE
	    COPY(NESTLEV) = EVALIF(ULINE(4:), NSTR, MSPEC)
	  ENDIF
	  ELIF(NESTLEV) = COPY(NESTLEV)
	  DELSE(NESTLEV) = .FALSE.

C Check for @ELSEIFs and deal with them

	ELSE IF (LINE(:7).EQ.'@ELSEIF') THEN
	  WRFLAG = .FALSE.

C If a previous condition has been printed then we never check the @elseif
C - uses the ELIF logical array
C Cannot have a @elseif without a @if
C An @elseif cannot follow a @else for the same @if 

	  IF (DEBUG) WRITE(UERR,1110) SPC(:2*NESTLEV), LINE(8:LLEN)
	  IF (NESTLEV.EQ.0) THEN
	    WRITE(UERR,1010) LCNOW
	    CALL ERROREXIT
	  ELSE IF (DELSE(NESTLEV)) THEN
	    WRITE(UERR,1020) LCIN(NESTLEV), LCNOW
	    CALL ERROREXIT
	  ELSE IF (COPY(NESTLEV-1).EQV..FALSE.) THEN
	    COPY(NESTLEV)=.FALSE.
	  ELSE IF (ELIF(NESTLEV)) THEN
	    COPY(NESTLEV)=.FALSE.
	  ELSE
	    ULINE = LINE
	    COPY(NESTLEV) = EVALIF(ULINE(8:), NSTR, MSPEC)
	    ELIF(NESTLEV) = COPY(NESTLEV)
	  ENDIF

C Check for @ELSEs and deal with them
C Cannot have a @else without a @if
C Only one @else allowed for each @IF

	ELSE IF (LINE(:5).EQ.'@ELSE') THEN
	  WRFLAG = .FALSE.

	  IF (DEBUG) WRITE(UERR,1120) SPC(:2*NESTLEV)
	  IF (NESTLEV.EQ.0) THEN
	    WRITE(UERR,1030) LCNOW
	    CALL ERROREXIT
	  ELSE IF (DELSE(NESTLEV)) THEN
	    WRITE(UERR,1040) LCIN(NESTLEV), LCNOW
	    CALL ERROREXIT
	  ELSE
	    DELSE(NESTLEV) = .TRUE.
	    COPY(NESTLEV) = COPY(NESTLEV-1).AND..NOT.ELIF(NESTLEV)
	  ENDIF

C Now transfer the line if wanted

	ELSE IF (LINE(:6).NE.'@ENDIF') THEN

	  IF (COPY(NESTLEV)) THEN
	    IF (LLEN.EQ.0) THEN
	      WRITE(UN2,*)
	    ELSE IF (LINE(1:1).NE.'@') THEN
	      WRITE(UN2,'(A)') LINE(:LLEN)
	      IF (.NOT.WRFLAG) THEN
		IF (DEBUG) WRITE(UERR,1130) SPC(:2*NESTLEV+2)
		WRFLAG = .TRUE.
	      ENDIF
	    ENDIF
	  ELSE
	    IF (.NOT.WRFLAG) THEN
	      IF (DEBUG) WRITE(UERR,1140) SPC(:2*NESTLEV+2)
	      WRFLAG = .TRUE.
	    ENDIF
	  ENDIF

C Check for @ENDIFs and deal with them

	ELSE
	  WRFLAG = .FALSE.

	  IF (DEBUG) WRITE(UERR,1150) SPC(:2*NESTLEV)
	  IF (NESTLEV.EQ.0) THEN
	    WRITE(UERR,1050) LCNOW
	    CALL ERROREXIT
	  ENDIF
	  NESTLEV = NESTLEV - 1

	ENDIF

C Now read the next line of input

	GOTO 10
99	RETURN

C Format statements for error messages and debug output

1000	FORMAT('Nesting of conditions too deep (line',I4,')')
1010	FORMAT('@ELSEIF found outside @IF...@ENDIF (line',I4,')')
1020	FORMAT('@ELSEIF follows a @ELSE statment (lines',I4,' -',I4,')')
1030	FORMAT('@ELSE found outside @IF...@ENDIF (line',I4,')')
1040	FORMAT('@IF has multiple @ELSE statments (lines',I4,' -',I4,')')
1050	FORMAT('@ENDIF without an @IF statement (line',I4,')')

1100	FORMAT(A,'@IF condition: ',A)
1110	FORMAT(A,'@ELSEIF condition: ',A)
1120	FORMAT(A,'@ELSE')
1130	FORMAT(A,'Writing to output file...')
1140	FORMAT(A,'Skipping lines from input file...')
1150	FORMAT(A,'@ENDIF')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Evaluate the IF or ELSEIF string returning true or false. Also need
C to pass the machine name array

	LOGICAL FUNCTION EVALIF(STRING, NSTR, MSPEC)
	CHARACTER*(*) STRING
	INTEGER NSTR
	CHARACTER*(*) MSPEC(NSTR)

	INTEGER I, LMSP, IX
	LOGICAL EVALCOND

C Uppercase the line

	CALL UPCASE(STRING)

C Substitute '+' for each string matching the name list

	DO I = 1, NSTR
	  DO LMSP = 1, LEN(MSPEC(I))
	    IF (MSPEC(I)(LMSP:LMSP).EQ.' ') GOTO 10
	  ENDDO
10	  LMSP = LMSP - 1
20	  IX = INDEX(STRING,MSPEC(I)(:LMSP))
	  IF (IX.GT.0) THEN
	    STRING(IX:IX+LMSP-1) = '+'
	    GOTO 20
	  ENDIF
	ENDDO

C Sets of characters remaining other than ' +-!()&|' should now be set to '-'

	DO I = 1, LEN(STRING)
	  IX = INDEX(' +-!()&|', STRING(I:I))
	  IF (IX.EQ.0) THEN
	    IF (I.EQ.1) THEN
	      STRING(I:I) = '-'
	    ELSE IF (STRING(I-1:I-1).EQ.'-') THEN
	      STRING(I-1:I) = ' -'
	    ELSE
	      STRING(I:I) = '-'
	    ENDIF
	  ENDIF
	ENDDO

C Now evaluate the condition statement

	EVALIF = EVALCOND(STRING)

	RETURN

	END

C Convert string to upper case

	SUBROUTINE UPCASE(STRING)
	CHARACTER*(*) STRING

	INTEGER I, IC

	DO I = 1, LEN(STRING)
	  IC = ICHAR(STRING(I:I))
	  IF (IC.GE.ICHAR('a') .AND. IC.LE.ICHAR('z')) THEN
	    STRING(I:I) = CHAR(ICHAR('A')+IC-ICHAR('a'))
	  ENDIF
	ENDDO

	RETURN
	END

C Convert string to lower case

	SUBROUTINE LOCASE(STRING)
	CHARACTER*(*) STRING

	INTEGER I, IC

	DO I = 1, LEN(STRING)
	  IC = ICHAR(STRING(I:I))
	  IF (IC.GE.ICHAR('A') .AND. IC.LE.ICHAR('Z')) THEN
	    STRING(I:I) = CHAR(ICHAR('a')+IC-ICHAR('A'))
	  ENDIF
	ENDDO

	RETURN
	END

C Remove all spaces from the supplied string (or at least push them
C to the end of the string)

	SUBROUTINE REMSP(STRING)
	CHARACTER*(*) STRING

	INTEGER I

	DO I = 1, LEN(STRING)
10	  IF (STRING(I:).EQ.' ') RETURN
	  IF (STRING(I:I).EQ.' ') THEN
	    STRING(I:) = STRING(I+1:)
	    GOTO 10
	  ENDIF
	ENDDO

	RETURN
	END

C Subroutine to exit the program, as this is no longer a standard...

	SUBROUTINE ERROREXIT

	STOP 'Program "makevers" could not parse input file'

	RETURN
	END

C Evaluate a logical condition returning .TRUE. or .FALSE.
C Lists of adjacent TRUEs or FALSEs are ORed together

	LOGICAL FUNCTION EVALCOND(STRING)
	CHARACTER*(*) STRING

	INTEGER NCONDS
	PARAMETER (NCONDS=20)

	CHARACTER*3 COND(NCONDS)
	CHARACTER*1 RES(NCONDS)
	INTEGER	LCOND(NCONDS)
	INTEGER I, IX

	DATA COND /'(+)','(-)',
     +		   '!+','!-',
     +		   '+&+','+&-','-&+','-&-',
     +		   '+%+','+%-','-%+','-%-',
     +		   '+|+','+|-','-|+','-|-',
     +		   '++','+-','-+','--'/

	DATA RES  /'+','-',
     +		   '-','+',
     +		   '+','-','-','-',
     +		   '-','+','+','-',
     +		   '+','+','+','-',
     +		   '+','+','+','-'/

	DATA LCOND/3,3,
     +		   2,2,
     +		   3,3,3,3,
     +		   3,3,3,3,
     +		   3,3,3,3,
     +		   2,2,2,2/

C	WRITE(0,1000) STRING

C Line may contain spaces so remove them all

	CALL REMSP(STRING)
	IF (STRING.EQ.' ') GOTO 98
	GOTO 40

C Now repeatedly replace sub-strings of the logical construction
C until no further changes are made

30	CONTINUE
C	WRITE(0,1005) STRING
40	DO I = 1, NCONDS
	  IX = INDEX(STRING,COND(I))
	  IF (IX.NE.0) THEN
	    STRING(IX:) = RES(I)//STRING(IX+LCOND(I):)
	    GOTO 30
	  ENDIF
	ENDDO

C Now test the end result - only '+' is taken as true

	IF (STRING.EQ.'+') GOTO 99

98	EVALCOND = .FALSE.
C	WRITE(0,1010)
	RETURN

99	EVALCOND = .TRUE.
C	WRITE(0,1020)
	RETURN

C1000	FORMAT(' @IF Logically: ',A64)
C1005	FORMAT(' Will simplify: ',A64)
C1010	FORMAT(' @IF Result --> FALSE')
C1020	FORMAT(' @IF Result --> TRUE')
	END
