C***********************************************************************
C Interface routines between GETWORD 1.4 and P2. These are simply to
C prolong the lifetime of GETWORD and the programs which use it! It
C should not be machine-dependent.
C
C Robert Esnouf, 2/7/93 - 16/7/93
C
C***********************************************************************

C***********************************************************************
C Initialise the parser, opening the standard unit numbers. UNIT is the
C unit number that GETWORD thinks it is using, though G2 knows better
C and maps it to unit IPUNIT (always 5?).
C RME 2/7/93 - 2/7/93
C***********************************************************************

      LOGICAL FUNCTION G2INIT(UNIT)
      INTEGER UNIT

      INCLUDE 'p2.inc'
      INCLUDE 'g2.inc'

      LOGICAL P2INIT
      INTEGER I

      DO I = G2MINU, G2MAXU
        G2UNIT(I) = NOUNIT
      ENDDO

      IF (P2INIT()) THEN
        G2UNIT(UNIT) = IPUNIT
        UNIT = OPUNIT
        IF (CHATTY) WRITE(OPUNIT,1000)
        IF (CHATTY .AND. ISTERM) WRITE(OPUNIT,1010)
        IF (CHATTY .AND. .NOT.ISTERM) WRITE(OPUNIT,1020)
        G2INIT = .TRUE.
      ELSE
        G2INIT = .FALSE.
      ENDIF

      RETURN

1000  FORMAT('G2INIT: P2 Initialised successfully')
1010  FORMAT('G2INIT: Standard input is interactive terminal')
1020  FORMAT('G2INIT: Standard input is non-interactive')
      END

C***********************************************************************
C Open a file other than one of the standard input files, this is used
C by GETWORD for all calls other than the initialising call. Again UNIT
C is the number that GETWORD thinks it is using although G2 uses the
C P2OPEN function to find a suitable unit number and sets up mapping
C between the two. The file is assumed to exist and be formatted.
C RME 2/7/93 - 2/7/93
C***********************************************************************

      LOGICAL FUNCTION G2OPEN(UNIT,FNAM)
      INTEGER UNIT
      CHARACTER*(*) FNAM

      INCLUDE 'p2.inc'
      INCLUDE 'g2.inc'

      INTEGER IS, IE, P2UNIT
      LOGICAL P2OPEN

C If the file unit is already being used then close it first

      IF (G2UNIT(UNIT).NE.NOUNIT) CALL P2CLOS(G2UNIT(UNIT), ' ')

      IF (P2OPEN(P2UNIT, FNAM, 'RF', 0)) THEN
        G2UNIT(UNIT) = P2UNIT
        G2OPEN = .TRUE.
        IF (TALK) THEN
          CALL P2SLEN(FNAM, IS, IE)
          WRITE(OPUNIT,1000) FNAM(IS:IE), P2UNIT, UNIT
        ENDIF
      ELSE
        G2OPEN = .FALSE.
      ENDIF

      RETURN

1000  FORMAT('G2OPEN: File ',A,' opened using FORTRAN unit ',I2,
     +' as GETWORD unit ',I2)
      END

C***********************************************************************
C Close a file after the EOF has been read - a slight enhancement over
C the original GETWORD.
C RME 16/7/93 - 16/7/93
C***********************************************************************

      SUBROUTINE G2CLOS(UNIT)
      INTEGER UNIT

      INCLUDE 'p2.inc'
      INCLUDE 'g2.inc'

      CHARACTER*80 FNAM
      INTEGER IS, IE

C Don't try to close it if it's not open

      IF (G2UNIT(UNIT).NE.NOUNIT) THEN
        FNAM = FNAME(G2UNIT(UNIT))
        CALL P2CLOS(G2UNIT(UNIT), ' ')
        IF (TALK) THEN
          CALL P2SLEN(FNAM, IS, IE)
          WRITE(OPUNIT,1000) FNAM(IS:IE), G2UNIT(UNIT), UNIT
        ENDIF
        G2UNIT(UNIT) = NOUNIT
      ENDIF

      RETURN

1000  FORMAT('G2CLOS: File ',A,' closed on FORTRAN unit ',I2,
     +' and GETWORD unit ',I2)
      END

C***********************************************************************
C Read a line from the input file that GETWORD thinks is attached to
C UNIT. It can use the string PROMPT as a prompt, leaving P2 to decide
C whether a prompt is needed or if the prompt and line should be echoed
C to the output file. If ECHO is not zero then no prompt will be used
C if the input is not from an interactive terminal. The line is returned
C completely unparsed in string LINE, and the number of characters read
C into the string is ILEN. If ILEN is zero then a blank line was read
C and if ILEN is -1 then the end of file was reached.
C RME 2/7/93 - 2/7/93
C***********************************************************************

      SUBROUTINE G2READ(UNIT,PROMPT,ECHO,LINE,ILEN)
      INTEGER UNIT
      CHARACTER*(*) PROMPT
      INTEGER ECHO
      CHARACTER*(*) LINE
      INTEGER ILEN

      INCLUDE 'p2.inc'
      INCLUDE 'g2.inc'

C If the input is from a terminal then produce a prompt and read a line

      IF (G2UNIT(UNIT).EQ.IPUNIT .AND. ISTERM) THEN
        CALL M2PROM(PROMPT)
        CALL M2READ(IPUNIT, LINE, ILEN)

C If we're not at a terminal and there is echo then read a line and echo
C it

      ELSE
        CALL M2READ(G2UNIT(UNIT), LINE, ILEN)
        IF (ECHO.EQ.0 .AND. ILEN.GT.0) THEN
          WRITE(OPUNIT, 1000) PROMPT, LINE(:ILEN)
        ELSE IF (ECHO.EQ.0 .AND. ILEN.EQ.0) THEN
          WRITE(OPUNIT, 1000) PROMPT, ' '
        ENDIF
      ENDIF

      RETURN

1000  FORMAT(A,'> ',A)
      END

