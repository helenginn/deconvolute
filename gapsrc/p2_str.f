C***********************************************************************
C Useful string functions from the P2 parser. Don't get too excited
C but its better than rewriting them every time.
C
C Robert Esnouf, 4/7/93 - 22/11/93
C
C***********************************************************************

C***********************************************************************
C Convert the case of a string depending on the values of the logical
C variables UPPERC and LOWERC. If UPPERC is .TRUE. then the string is
C uppercased otherwise if LOWERC is .TRUE. then it is lowercased other-
C wise the string is returned unchanged (the default behaviour)
C RME 22/11/93 - 22/11/93
C***********************************************************************

      SUBROUTINE P2CNV(STRING)
      CHARACTER*(*) STRING

      INCLUDE 'p2.inc'

      IF (UPPERC) THEN
        CALL P2UCAS(STRING)
      ELSE IF (LOWERC) THEN
        CALL P2LCAS(STRING)
      ENDIF

      RETURN
      END

C***********************************************************************
C Uppercase a string in situ. Assuming no crappy character sets so that
C all letters are sequential
C RME 21/11/93 - 21/11/93
C***********************************************************************

      SUBROUTINE P2UCAS(STRING)
      CHARACTER*(*) STRING

      INTEGER I, IOFF

      IOFF = ICHAR('A') - ICHAR('a')
      DO I = 1, LEN(STRING)
        IF (STRING(I:I).GE.'a' .AND. STRING(I:I).LE.'z') THEN
          STRING(I:I) = CHAR(ICHAR(STRING(I:I)) + IOFF)
        ENDIF
      ENDDO

      RETURN
      END

C***********************************************************************
C Lowercase a string in situ. Assuming no crappy character sets so that
C all letters are sequential
C RME 21/11/93 - 21/11/93
C***********************************************************************

      SUBROUTINE P2LCAS(STRING)
      CHARACTER*(*) STRING

      INTEGER I, IOFF

      IOFF = ICHAR('a') - ICHAR('A')
      DO I = 1, LEN(STRING)
        IF (STRING(I:I).GE.'A' .AND. STRING(I:I).LE.'Z') THEN
          STRING(I:I) = CHAR(ICHAR(STRING(I:I)) + IOFF)
        ENDIF
      ENDDO

      RETURN
      END

C***********************************************************************
C Compare two strings either case sensitively or case insensitively
C depending on the logical variable CASEIN. Does not alter the strings.
C RME 22/11/93 - 22/11/93
C***********************************************************************

      LOGICAL FUNCTION P2CMP(STR1, STR2)
      CHARACTER*(*) STR1
      CHARACTER*(*) STR2

      INCLUDE 'p2.inc'

      LOGICAL P2ICMP

      IF (CASEIN) THEN
        P2CMP = P2ICMP(STR1, STR2)
      ELSE
        P2CMP = (STR1.EQ.STR2)
      ENDIF

      RETURN
      END

C***********************************************************************
C Compare two strings either case insensitively. Does not alter the
C strings passed to it. For switchable case comparisons use the P2CMP
C call.
C RME 22/11/93 - 22/11/93
C***********************************************************************

      LOGICAL FUNCTION P2ICMP(STR1, STR2)
      CHARACTER*(*) STR1
      CHARACTER*(*) STR2

      INTEGER I, IOFF, I1, I2

C Get the lengths of the strings

      I1 = LEN(STR1)
      I2 = LEN(STR2)
      P2ICMP = .FALSE.
      IOFF = ICHAR('a') - ICHAR('A')

C Check that if one string is longer than the other it is filled only
C with blanks

      IF (I1.GT.I2) THEN
        IF (STR1(I2+1:).NE.' ') RETURN
      ELSE IF (I1.LT.I2) THEN
        IF (STR2(I1+1:).NE.' ') RETURN
      ENDIF

C Now compare the main parts of the strings

      DO I = 1, MIN(I1, I2)
        IF (STR1(I:I).NE.STR2(I:I)) THEN
          IF (STR1(I:I).GE.'A' .AND. STR1(I:I).LE.'Z') THEN
            IF (ICHAR(STR2(I:I))-ICHAR(STR1(I:I)) .NE. IOFF) RETURN
          ELSE IF (STR1(I:I).GE.'a' .AND. STR1(I:I).LE.'z') THEN
            IF (ICHAR(STR1(I:I))-ICHAR(STR2(I:I)) .NE. IOFF) RETURN
          ELSE
            RETURN
          ENDIF
        ENDIF
      ENDDO

      P2ICMP = .TRUE.
      RETURN
      END

C***********************************************************************
C Find out the length of the used part of the string in the supplied
C string.
C RME 5/7/93 - 5/7/93
C***********************************************************************

      SUBROUTINE P2SLEN(STRING, IS, IE)
      CHARACTER*(*) STRING
      INTEGER IS
      INTEGER IE

      INTEGER I

      IE = 1
      IF (STRING.EQ.' ') THEN
        IS = 1
      ELSE
        IS = 0
        DO I = 1, LEN(STRING)
          IF (IS.EQ.0 .AND. STRING(I:I).NE.' ') IS = I
          IF (IS.NE.0 .AND. STRING(I:I).NE.' ') IE = I
        ENDDO
      ENDIF

      RETURN
      END

C***********************************************************************
C Crunch multiple spaces out of a string (as well as leading ones) and
C return the actual length of the string - slow but makes pretty format-
C ting easier. If the string is empty the ILEN is set to 1.
C RME 14/7/93 - 14/7/93
C***********************************************************************

      SUBROUTINE P2SSET(STRING, ILEN)
      CHARACTER*(*) STRING
      INTEGER ILEN

      INTEGER I

C Check for an empty string and return immediately

      IF (STRING.EQ.' ' .OR. LEN(STRING).EQ.1) THEN
        ILEN = 1
      ELSE
10      IF (STRING(1:1).EQ.' ') THEN
          STRING = STRING(2:)
          GOTO 10
        ENDIF
        I = 2
20      IF (STRING(I-1:I).EQ.'  ') THEN
          STRING(I-1:) = STRING(I:)
        ELSE
          I = I + 1
        ENDIF
        IF (I.GT.LEN(STRING)) GOTO 30
        IF (STRING(I:).NE.' ') GOTO 20
30      ILEN = I - 1
      ENDIF

      RETURN
      END

C***********************************************************************
C Format a value in hundredths of a second in to hours, minutes and
C seconds neatly. This is mainly used for formatting CPU usage strings.
C RME 20/7/93 - 20/7/93
C***********************************************************************

      SUBROUTINE P2FTIM(ICENT, STRING, ILEN)
      INTEGER ICENT
      CHARACTER*(*) STRING
      INTEGER ILEN

      INTEGER IHRS, IMINS, IHSECS
      CHARACTER*28 STIM

      IHRS = ICENT/360000
      IMINS = (ICENT - 360000*IHRS)/6000
      IHSECS = ICENT - 6000*IMINS - 360000*IHRS
      STIM = '      hrs    mins       secs'
      WRITE(STIM(1:5),'(I5)') IHRS
      WRITE(STIM(11:12),'(I2)') IMINS
      WRITE(STIM(19:23),'(F5.2)') IHSECS/100.0
      IF (IHRS.EQ.0) STIM(1:9) = '         '
      IF (IHRS.EQ.1) STIM(9:9) = ' '
      IF (IMINS.EQ.0) STIM(10:17) = '        '
      IF (IMINS.EQ.1) STIM(17:17) = ' '
      CALL P2SSET(STIM, ILEN)

      ILEN = MIN(LEN(STRING),LEN(STIM))
      STRING = STIM(:ILEN)

      RETURN
      END

C***********************************************************************
C Function to return which keyword has been typed. This is more
C complex than it might seem at first since we don't know the
C length of the word or of any option, the case comparisons should
C be case insensitive (unless requested), commands may be abbreviated
C and one command may be a substring of another. Waaah! Return one of:
C the option number for a command, 0 for no match or -1 for ambiguity
C RME 22/11/93 - 22/11/93
C***********************************************************************

      INTEGER FUNCTION P2NOPT(OPTS, WORD)
      CHARACTER*(*) OPTS
      CHARACTER*(*) WORD

      INCLUDE 'p2.inc'

      INTEGER I1, I2, IPS, IPE, CP, LO, LW, IOFF, OPNUM
      INTEGER IS, IE
      INTEGER LMATCH
      LOGICAL EXACT, AMBIG

      IOFF = ICHAR('a') - ICHAR('A')

C Set the best match details to nothing

      LMATCH = 0
      P2NOPT = 0
      AMBIG = .FALSE.
      EXACT = .FALSE.

C Get the length of the options string and the word to be matched

      CALL P2SLEN(WORD, IS, IE)
      LW = IE - IS + 1
      LO = LEN(OPTS)

      CP = 1
      IPS = 1
      OPNUM = 1

10    IF (CP.LE.LO .AND. OPTS(CP:CP).NE.'!') THEN
        CP = CP + 1
        GOTO 10
      ENDIF

C Find the length of the command in OPTS (trailing spaces don't count)

      IPE = CP - 1
20    IF (IPE.GT.IPS .AND. OPTS(IPE:IPE).EQ.' ') THEN
        IPE = IPE - 1
        GOTO 20
      ENDIF

C Blank commands ('...!!...' or '...!   !...') will match anything,
C but with very low priority, otherwise find how many characters match
C the input word. Now match against the word. We will need to know if
C the match is a perfect one or not.

      IF (IPS.LE.IPE .AND. OPTS(IPS:IPE).NE.' ') THEN

        DO I2 = 1, MIN(LW, IPE-IPS+1)
          I1 = IPS + I2 - 1
          IF (OPTS(I1:I1).NE.WORD(I2:I2)) THEN
            IF (.NOT.CASEIN) THEN
              GOTO 30
            ELSE IF (OPTS(I1:I1).GE.'A' .AND. OPTS(I1:I1).LE.'Z') THEN
              IF(ICHAR(WORD(I2:I2))-ICHAR(OPTS(I1:I1)).NE.IOFF)GOTO 30
            ELSE IF (OPTS(I1:I1).GE.'a' .AND. OPTS(I1:I1).LE.'z') THEN
              IF(ICHAR(OPTS(I1:I1))-ICHAR(WORD(I2:I2)).NE.IOFF)GOTO 30
            ELSE
              GOTO 30
            ENDIF
          ENDIF
        ENDDO

C Have matched as well as possible and I2 is one more than the number
C of characters matched. If we have a perfect match then set the
C EXACT logical to TRUE - but both the length and contents have to
C match for that. We can't have unmatched characters left in BOTH the
C option and word strings - if that's the case (eg: "RMAP" against
C "RMERGE") then the option doesn't count.

30      I2 = I2 - 1
        IF (I2.EQ.LW .OR. I2.EQ.IPE-IPS+1) THEN
          IF (I2.GT.LMATCH) THEN
            P2NOPT = OPNUM
            LMATCH = I2
            AMBIG = .FALSE.
            EXACT = (I2.EQ.LW .AND. I2.EQ.IPE-IPS+1)
          ELSE IF (I2.EQ.LMATCH) THEN
            IF (EXACT) THEN
              IF (I2.EQ.LW .AND. I2.EQ.IPE-IPS+1) THEN
                AMBIG = .TRUE.
              ENDIF
            ELSE
              IF (I2.EQ.LW .AND. I2.EQ.IPE-IPS+1) THEN
                P2NOPT = OPNUM
                AMBIG = .FALSE.
                EXACT = .TRUE.
              ELSE IF (P2NOPT.NE.0) THEN
                AMBIG = .TRUE.
              ENDIF
            ENDIF
          ENDIF
        ELSE IF (LMATCH.NE.0 .OR. P2NOPT.EQ.0) THEN
          IF (I2.GE.LMATCH .AND. .NOT.EXACT) THEN
            LMATCH = I2
            P2NOPT = 0
          ENDIF
        ENDIF

C The options string contains a blank option. Match this unless anything else
C matches at all.

      ELSE IF (P2NOPT.EQ.0) THEN

        LMATCH = 0
        P2NOPT = OPNUM

      ENDIF

C Now see if we're at the end of the OPTS string

      IF (CP.LE.LO) THEN
        OPNUM = OPNUM + 1
        IPS = CP + 1
        CP = IPS
        GOTO 10
      ENDIF

C Finally (!) we can return a result. P2NOPT is OK unless we've found
C an ambiguous case

      IF (AMBIG) P2NOPT = -1

      RETURN
      END

C***********************************************************************
C Subroutine to display the options actually contained in an options
C string. Extra spaces on the end of options are ignored and null
C options (spaces only or empty strings are ignored). Fields of output
C are 8 characters-wide with at least 1 space at the right hand end. FW
C is the length to allow for each field in the output. OUTP is the
C string to use to construct the output lines and FOUT is the format
C to use for printing output.
C RME 23/11/93 - 23/11/93
C***********************************************************************

      SUBROUTINE P2DOPT(OPTS, OUTP, FOUT, FW)
      CHARACTER*(*) OPTS
      CHARACTER*(*) OUTP
      CHARACTER*(*) FOUT
      INTEGER FW

      INCLUDE 'p2.inc'

      INTEGER IPS, IPE, CP, LO
      INTEGER LTOT, LUSED, LREQ

C Get the length of the output line

      LTOT = LEN(OUTP)
      LUSED = 0
      OUTP = ' '

C Get the length of the options string

      CP = 1
      IPS = 1
      LO = LEN(OPTS)

10    IF (CP.LE.LO .AND. OPTS(CP:CP).NE.'!') THEN
        CP = CP + 1
        GOTO 10
      ENDIF

C Find the length of the command in OPTS (trailing spaces don't count)

      IPE = CP - 1
20    IF (IPE.GT.IPS .AND. OPTS(IPE:IPE).EQ.' ') THEN
        IPE = IPE - 1
        GOTO 20
      ENDIF

C Blank commands are ignored, otherwise try to add output field to the
C current list. We need LREQ spaces to be left in OUTP. If that would make
C the line too long then print out what we've got

      IF (OPTS(IPS:IPE).NE.' ') THEN
        LREQ = MIN(IPE-IPS+1, LTOT)
        IF (LUSED+LREQ .GT. LTOT) THEN
          WRITE(OPUNIT,FOUT) OUTP(:LUSED)
          OUTP = ' '
          LUSED = 0
        ENDIF
        OUTP(LUSED+1:LUSED+LREQ) = OPTS(IPS:IPE)
        LUSED = MIN(((LUSED+LREQ+FW)/FW)*FW, LTOT)
      ENDIF

C Now see if we're at the end of the OPTS string

      IF (CP.LT.LO) THEN
        IPS = CP + 1
        CP = IPS
        GOTO 10
      ENDIF

C Print out any remaining options string

      IF (LUSED.GT.0) WRITE(OPUNIT,FOUT) OUTP(:LUSED)

      RETURN
      END

C***********************************************************************
C Subroutine to display the options actually contained in an options
C array. Extra spaces on the end of options are ignored and null
C options (spaces only or empty strings are ignored). FW
C is the length to allow for each field in the output. OUTP is the
C string to use to construct the output lines and FOUT is the format
C to use for printing output.
C RME 14/1/94 - 14/1/94
C***********************************************************************

      SUBROUTINE P2DARY(OPTS, NOPTS, OUTP, FOUT, FW)
      INTEGER NOPTS
      CHARACTER*(*) OPTS(NOPTS)
      CHARACTER*(*) OUTP
      CHARACTER*(*) FOUT
      INTEGER FW

      INCLUDE 'p2.inc'

      INTEGER I, IS, IE
      INTEGER LTOT, LUSED, LREQ

C Get the length of the output line

      LTOT = LEN(OUTP)
      LUSED = 0
      OUTP = ' '

C Go through the options string

      DO I = 1, NOPTS

C Find the length of the command in OPTS array (trailing spaces
C don't count). Blank commands are ignored.

        IF (OPTS(I).NE.' ') THEN
          CALL P2SLEN(OPTS(I), IS, IE)

C Try to add output field to the current list. We need LREQ spaces
C to be left in OUTP. If that would make the line too long then print
C out what we've got

          LREQ = MIN(IE-IS+1, LTOT)
          IF (LUSED+LREQ .GT. LTOT) THEN
            WRITE(OPUNIT,FOUT) OUTP(:LUSED)
            OUTP = ' '
            LUSED = 0
          ENDIF
          OUTP(LUSED+1:LUSED+LREQ) = OPTS(I)(IS:IE)
          LUSED = MIN(((LUSED+LREQ+FW)/FW)*FW, LTOT)

        ENDIF
      ENDDO

C Print out any remaining options string

      IF (LUSED.GT.0) WRITE(OPUNIT,FOUT) OUTP(:LUSED)

      RETURN
      END

