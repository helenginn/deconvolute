C***********************************************************************
C File open and close functions for the new parser P2. These routines
C give basic file open and close calls which are machine independent
C that then make calls suitable for a particular machine. These are the
C reworked versions
C
C Robert Esnouf, 19/11/93 - 22/11/93
C
C***********************************************************************

C***********************************************************************
C Open a file with the specified characteristics. The unit number is
C assigned by this routine, so if this routine is used to do all file
C opens then there will be no problems with units clashing. UNIT is the
C selected unit, FNAME is the filename to use and CHRS are the file
C characteristics to use. If the filename is blank then the 'VMS'-type
C filename FOR005.DAT will be used. The units are used in reverse order as
C I think this will reduce the chance of an accidental clash with a
C program which is not using the parser fully.
C RME 5/7/93 - 21/11/93
C***********************************************************************

      LOGICAL FUNCTION P2OPEN(UNIT, FNAM, CHRS, RECL)
      INTEGER UNIT, RECL
      CHARACTER*(*) FNAM
      CHARACTER*(*) CHRS

      INCLUDE 'p2.inc'

      CHARACTER*10 TCHRS
      INTEGER I, IS, IE
      LOGICAL M2OPEN, P2OCHK

C Sort out the characteristics string first of all. If no options are
C specified then use UNKNOWN and FORMATTED. If there is an error then
C use what options one can guess. If the options are actually incompat-
C ible then return, signalling an error

      CALL P2SLEN(CHRS, IS, IE)
      IF (CHATTY) WRITE(OPUNIT,1000) CHRS(IS:IE)
      IF (.NOT.P2OCHK(CHRS(IS:IE), TCHRS)) THEN
        IF (.NOT.CHATTY) WRITE(OPUNIT,1000) CHRS(IS:IE)
        WRITE(OPUNIT,1010)
        P2OPEN = .FALSE.
        RETURN
      ELSE IF (TALK) THEN
        WRITE(OPUNIT, 1020) TCHRS
      ENDIF

C Find an available unit for the file. To avoid the chance of a clash
C with files not using the P2OPEN command, the status of the free
C unit is checked before opening

      DO I = UNMAX, UNMIN, -1
        IF (FUNIT(I)(1:1).EQ.'-') THEN
          CALL M2FCHK(I)
          IF (FUNIT(I)(1:1).EQ.'-') GOTO 10
        ENDIF
      ENDDO

C Cannot find a free unit, so this needs an error message

      WRITE(OPUNIT, 1030)
      P2OPEN = .FALSE.
      RETURN

C Now we have a safe unit number construct a file name if necessary

10    UNIT = I
      FUNIT(UNIT) = TCHRS
      FRECL(UNIT) = RECL
      IF (FNAM.EQ.' ') THEN
        FNAME(UNIT) = 'FOR000.DAT'
        IF (UNIT.LT.10) THEN
          WRITE(FNAME(UNIT)(6:6),'(I1)') UNIT
        ELSE IF (UNIT.LT.100) THEN
          WRITE(FNAME(UNIT)(5:6),'(I2)') UNIT
        ELSE
          WRITE(FNAME(UNIT)(4:6),'(I3)') UNIT
        ENDIF
      ELSE
        FNAME(UNIT) = FNAM
      ENDIF

C Now actually try to open the file

      IF (M2OPEN(UNIT)) THEN
        P2OPEN = .TRUE.
        IF (TALK) THEN
          CALL P2SLEN(FNAME(UNIT), IS, IE)
          WRITE(OPUNIT, 1040) FNAME(UNIT)(IS:IE), UNIT
        ENDIF
      ELSE
        P2OPEN = .FALSE.
        FNAME(UNIT) = ' '
        FUNIT(UNIT) = '-'
      ENDIF

      IF (CHATTY) CALL P2FLST

      RETURN

1000  FORMAT('P2OPEN: Options sent to open statement: ',A)
1010  FORMAT('P2OPEN: Incompatible modes, could not open file')
1020  FORMAT('P2OPEN: Options used on open statement: ',A)
1030  FORMAT('P2OPEN: Could not find a free FORTRAN unit for OPENing')
1040  FORMAT('P2OPEN: Opened file ',A,' on FORTRAN unit ',I2)

      END

C***********************************************************************
C Close a file. If the unit number is 0 then the file with the name
C matching FNAM is closed otherwise the file on unit number UNIT is
C closed.
C RME 5/7/93 - 21/11/93
C***********************************************************************

      SUBROUTINE P2CLOS(UNIT, FNAM)
      INTEGER UNIT
      CHARACTER*(*) FNAM

      CALL P2CLS2(UNIT, FNAM, ' ')
      RETURN
      END

C***********************************************************************
C Close a file. If the unit number is 0 then the file with the name
C matching FNAM is closed otherwise the file on unit number UNIT is
C closed. The file is deleted on CLOSE whether or not it is SCRATCH.
C RME 5/7/93 - 21/11/93
C***********************************************************************

      SUBROUTINE P2CLSD(UNIT, FNAM)
      INTEGER UNIT
      CHARACTER*(*) FNAM

      CALL P2CLS2(UNIT, FNAM, 'D')
      RETURN
      END

C***********************************************************************
C Close a file. If the unit number is 0 then the file with the name
C matching FNAM is closed otherwise the file on unit number UNIT is
C closed. The file is kept on CLOSE whether or not it is SCRATCH.
C RME 5/7/93 - 21/11/93
C***********************************************************************

      SUBROUTINE P2CLSK(UNIT, FNAM)
      INTEGER UNIT
      CHARACTER*(*) FNAM

      CALL P2CLS2(UNIT, FNAM, 'K')
      RETURN
      END

C***********************************************************************
C Close a file. If the unit number is 0 then the file with the name
C matching FNAM is closed otherwise the file on unit number UNIT is
C closed. The characteristics string controls the type of close: DELETE
C KEEP or unspecified. It should only be called from the previous three
C routines so no error checking is performed
C RME 5/7/93 - 21/11/93
C***********************************************************************

      SUBROUTINE P2CLS2(UNIT, FNAM, CHRS)
      INTEGER UNIT
      CHARACTER*(*) FNAM
      CHARACTER*1 CHRS

      INCLUDE 'p2.inc'

      INTEGER IU, IS, IE
      LOGICAL M2CLOS

C For both UNIT to be 0 and FNAM to be ' ' is an error

      IF (UNIT.EQ.0 .AND. FNAM.EQ.' ') THEN
        WRITE(OPUNIT, 1000)
        RETURN
      ENDIF

C If the unit number is 0 then try to find the unit associated with
C the filename

      IF (UNIT.EQ.0) THEN
        DO IU = UNMIN, UNMAX
          IF (FNAME(IU).EQ.FNAM) GOTO 10
        ENDDO

C Did not match the filename, so cannot close

        CALL P2SLEN(FNAM, IS, IE)
        WRITE(OPUNIT, 1010) FNAM(IS:IE)
	RETURN
      ELSE
        IU = UNIT
      ENDIF

C Check that the unit number is not out of range

10    IF (IU.LT.UNMIN .OR. IU.GT.UNMAX) THEN
        WRITE(OPUNIT, 1020)

C Check we are not trying to close one of the standard FORTRAN units

      ELSE IF (FUNIT(IU)(2:2).EQ.'*') THEN
        WRITE(OPUNIT, 1030)

C Check that the unit is OPEN in the first place

      ELSE IF (FUNIT(IU)(1:1).EQ.'-') THEN
        WRITE(OPUNIT, 1040)

C Close the desired unit number

      ELSE
        CALL P2SLEN(FNAME(IU), IS, IE)
        IF (M2CLOS(IU, CHRS)) THEN
          IF (TALK) WRITE(OPUNIT, 1050) FNAME(IU)(IS:IE), IU
          FNAME(IU) = ' '
          FUNIT(IU) = '-'
        ELSE
          WRITE(OPUNIT, 1060) FNAME(IU)(IS:IE), IU
        ENDIF

C End of close routine

      ENDIF
      RETURN

1000  FORMAT('P2CLOS: No unit number or filename specified for CLOSE')
1010  FORMAT('P2CLOS: Could not find unit for ',A,' so cannot CLOSE')
1020  FORMAT('P2CLOS: Unit number out of range for CLOSE')
1030  FORMAT('P2CLOS: Cannot CLOSE standard FORTRAN I/O units')
1040  FORMAT('P2CLOS: No OPEN file to CLOSE on FORTRAN unit ',I2)
1050  FORMAT('P2CLOS: CLOSEd file ',A,' on FORTRAN unit ',I2)
1060  FORMAT('P2CLOS: Could not CLOSE file ',A,' on FORTRAN unit ',I2)
      END

C***********************************************************************
C Check the options on an open statement returning .TRUE. if the options
C string can be interpreted and .FALSE. if the options specified are
C actually incompatible. The 'parsed' options string is returned in the
C second argument. If any of the specified options are duplicated or
C cannot be interpreted then a warning message is produced.
C RME 21/11/93 - 21/11/93
C***********************************************************************

      LOGICAL FUNCTION P2OCHK(CHRS, TCHRS)
      CHARACTER*(*) CHRS
      CHARACTER*(*) TCHRS

      INCLUDE 'p2.inc'

      INTEGER I, IX
      CHARACTER*10 REFCHR
      CHARACTER*20 INCOMP
      CHARACTER*32 TCHRS2

C Set up reference list of options in correct order and list of
C incompatible pairs of options

      P2OCHK = .TRUE.
      REFCHR = '-*RWAS?FUD'
      INCOMP = 'RWRAR?RSWAWSA?AS?SFU'

C If no options are specified then use UNKNOWN and FORMATTED

      IF (CHRS.EQ.' ') THEN
        IF (TALK) WRITE(OPUNIT,1000) 'No options specified'
        TCHRS = '      ?F       '

C If the string is too long then use the default case and give an
C error message

      ELSE IF (LEN(CHRS) .GT. LEN(TCHRS2)) THEN
        WRITE(OPUNIT,1000) 'Options string too long'
        TCHRS = '      ?F       '

C Copy options string to a temporary string and convert to upper case

      ELSE

        TCHRS = ' '
        TCHRS2 = CHRS
        CALL P2UCAS(TCHRS2)

C Check for incompatibilities in the options

        DO I = 1, LEN(INCOMP), 2
          IF (INDEX(TCHRS2,INCOMP(I:I)).NE.0 .AND.
     +        INDEX(TCHRS2,INCOMP(I+1:I+1)).NE.0) THEN
            WRITE(OPUNIT, 1010) INCOMP(I:I), INCOMP(I+1:I+1)
            P2OCHK = .FALSE.
          ENDIF
        ENDDO
        IF (.NOT.P2OCHK) RETURN

C Check for the options excluding direct access and construct TCHRS

        DO I = 3, 9
          IX = INDEX(TCHRS2, REFCHR(I:I))
          IF (IX.NE.0) THEN
            TCHRS2(IX:IX) = ' '
            TCHRS(I:I) = REFCHR(I:I)
          ENDIF
        ENDDO

C Sort out direct access option. I only allow DIRECT access files for files
C opened status=UNKNOWN.

        IX = INDEX(TCHRS2, 'D')
        IF (IX.GT.0) THEN
          TCHRS(3:7) = '   ? '
          TCHRS(10:10) = 'D'
          TCHRS2(IX:IX) = ' '
        ENDIF

C If there is anything left in the TCHRS2 string then there was an
C error or multiple definition

        IF (TCHRS2.NE.' ') WRITE(OPUNIT, 1020)

C We managed to get an options string so the open can continue

      ENDIF
      RETURN

1000  FORMAT('P2OPEN: ',A,', opening file UNKNOWN and FORMATTED')
1010  FORMAT('P2OPEN: Options ''',A,''' and ''',A,''' incompatible')
1020  FORMAT('P2OPEN: Multiple declaration or invalid specification',
     +       ' of options')
      END

C***********************************************************************
C Test the status of all possible FORTRAN units available to the parser.
C This routine simply updates the FNAME and FUNIT arrays. It should only
C be used if the program is quite sure things have got badly screwed up
C since all parser routines check whether units are available before
C using them.
C RME 22/11/93 - 22/11/93
C***********************************************************************

      SUBROUTINE P2INQA

      INCLUDE 'p2.inc'

      INTEGER I

C Loop over all units checking them

      DO I = UNMIN, UNMAX
        CALL M2FCHK(I)
      ENDDO

      RETURN
      END

C***********************************************************************
C Return the options associated with a particular output unit. The
C routine checks the status of the file just in case the information
C in the file access arrays is incorrect. The filename and character-
C istics are returned in FNAM (max 80 chars), CHRS (max 10 chars) and
C RECL, the record length for DIRECT access files.
C RME 22/11/93 - 22/11/93
C***********************************************************************

      SUBROUTINE P2INQ(UNIT, FNAM, CHRS, RECL)
      INTEGER UNIT, RECL
      CHARACTER*(*) FNAM
      CHARACTER*(*) CHRS

      INCLUDE 'p2.inc'

C Check that the unit number is not out of range

      IF (UNIT.LT.UNMIN .OR. UNIT.GT.UNMAX) THEN
        WRITE(OPUNIT, 1000)
      ELSE
        CALL M2FCHK(UNIT)
        FNAM = FNAME(UNIT)
        CHRS = FUNIT(UNIT)
        RECL = FRECL(UNIT)
      ENDIF

      RETURN
1000  FORMAT('P2INQ: Unit number specified was out of range')
      END

C***********************************************************************
C Print out a summary of all the currently open files. Can be called
C from the main program but will eventually be a P2 command, possibly
C .PRINTFILES
C RME 22/11/93 - 22/11/93
C***********************************************************************

      SUBROUTINE P2FLST

      INCLUDE 'p2.inc'

      INTEGER I, IS, IE

C Print out a heading

      WRITE(OPUNIT,1000)
      WRITE(OPUNIT,1010)
      WRITE(OPUNIT,1020)

C Print out the unit, status string and name of open files

      DO I = UNMIN, UNMAX
        IF (FUNIT(I)(1:1).NE.'-') THEN
          CALL P2SLEN(FNAME(I), IS, IE)
          IF (FUNIT(I)(10:10).EQ.'D') THEN
            WRITE(OPUNIT,1030) I, FUNIT(I), FRECL(I), FNAME(I)(IS:IE)
          ELSE
            WRITE(OPUNIT,1035) I, FUNIT(I), FNAME(I)(IS:IE)
          ENDIF
        ENDIF
      ENDDO

C Print an explanation of the status string to finish off

      WRITE(OPUNIT,1020)
      WRITE(OPUNIT,1040)
      WRITE(OPUNIT,1050)
      WRITE(OPUNIT,1060)
      WRITE(OPUNIT,1070)
      WRITE(OPUNIT,1080)
      WRITE(OPUNIT,1000)

      RETURN

1000  FORMAT(' ')
1010  FORMAT('Unit|  OPEN options  |  Filename')
1020  FORMAT('----+----------------+---------------------------------')
1030  FORMAT(I3,' |',A10,I5,' | ',A)
1035  FORMAT(I3,' |',A10,5X,' | ',A)
1040  FORMAT('Key: ''*'' Standard IO unit,   ',
     +            '''R'' OPENed for reading,')
1050  FORMAT('     ''W'' OPENed for writing, ',
     +            '''A'' OPENed for write (append),')
1060  FORMAT('     ''?'' OPENed as UNKNOWN,  ',
     +            '''S'' OPENed as SCRATCH,')
1070  FORMAT('     ''F'' Formatted file,     ',
     +            '''U'' Unformatted file,')
1080  FORMAT('     ''D'' Direct access file (followed ',
     +       'by record length)')
      END
