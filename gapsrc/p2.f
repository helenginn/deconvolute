C***********************************************************************
C The new version parser, P2, by RME. It replaces the ailing GETWORD
C which is not particularly machine-independent. It is proposed to make
C this firstly act as a way of getting GETWORD to run tidily on all
C machines and then to replace it all together giving added functional-
C ity for both parsing and implementing normal machine-dependent
C actions.
C
C Robert Esnouf, 4/7/93 - 22/11/93
C
C***********************************************************************

C***********************************************************************
C Initialise the P2 parser. The first thing it does is construct arrays
C to control the allocation of FORTRAN units by calls to P2OPEN
C RME 5/7/93 - 24/7/93
C***********************************************************************

      LOGICAL FUNCTION P2INIT()

      INCLUDE 'p2.inc'

      INTEGER I
      LOGICAL M2STD

C Start the parser in terse mode (actually in chatty mode for debug)
C Select no case conversion and case insensitive comparisons by default

      CALL P2TERS
C     CALL P2CHAT
      CALL P2STCI
      CALL P2NCNV

C Initialise the FORTRAN unit numbers, and set up the values for the
C standard I/O units

      DO I = PUNMIN, PUNMAX
        FUNIT(I) = '-         '
        FRECL(I) = 0
        FNAME(I) = ' '
      ENDDO
      IF (M2STD()) THEN
        FUNIT(IPUNIT) = ' *R    F  '
        FNAME(IPUNIT) = 'Standard FORTRAN input unit'
        FUNIT(OPUNIT) = ' * W   F  '
        FNAME(OPUNIT) = 'Standard FORTRAN output unit'
        P2INIT = .TRUE.
      ELSE
        P2INIT = .FALSE.
      ENDIF

C Set up the variables for memory allocation

      CALL P2MEMI

C Find out system information, storing them in relevant strings and save
C the starting time and date also set the initial CPU used time to zero

      CALL M2INIT
      OLDCPU = 0
      CALL P2TIME(STTIM, LSTTIM)
      CALL P2DATE(STDAT, LSTDAT)
      IF (CHATTY) CALL P2STMP

C Set the default arrow type for constructing prompts

      CALL P2AROW('> ')

      RETURN
      END

C***********************************************************************
C Set the parser to operate in terse mode - no messages at all other
C than fatal error messages. This is the default mode for the parser.
C RME 5/7/93 - 5/7/93
C***********************************************************************

      SUBROUTINE P2TERS

      INCLUDE 'p2.inc'

      CHATTY = .FALSE.
      TALK = .FALSE.

      RETURN
      END

C***********************************************************************
C Set the parser to operate in 'talk' mode - occasional messages
C for important events.
C RME 5/7/93 - 5/7/93
C***********************************************************************

      SUBROUTINE P2TALK

      INCLUDE 'p2.inc'

      CHATTY = .FALSE.
      TALK = .TRUE.

      RETURN
      END

C***********************************************************************
C Set the parser to operate in verbose mode - loads of messages produced
C RME 5/7/93 - 5/7/93
C***********************************************************************

      SUBROUTINE P2CHAT

      INCLUDE 'p2.inc'

      CHATTY = .TRUE.
      TALK = .TRUE.

      RETURN
      END

C***********************************************************************
C Set the parser compare strings in a case insensitive fashion (default)
C RME 22/11/93 - 22/11/93
C***********************************************************************

      SUBROUTINE P2STCI

      INCLUDE 'p2.inc'

      CASEIN = .TRUE.

      RETURN
      END

C***********************************************************************
C Set the parser compare strings in a case sensitive fashion
C RME 22/11/93 - 22/11/93
C***********************************************************************

      SUBROUTINE P2STCS

      INCLUDE 'p2.inc'

      CASEIN = .FALSE.

      RETURN
      END

C***********************************************************************
C Set the parser not to convert the case of strings (default)
C RME 22/11/93 - 22/11/93
C***********************************************************************

      SUBROUTINE P2NCNV

      INCLUDE 'p2.inc'

      UPPERC = .FALSE.
      LOWERC = .FALSE.

      RETURN
      END

C***********************************************************************
C Set the parser to convert the case of strings to uppercase
C RME 22/11/93 - 22/11/93
C***********************************************************************

      SUBROUTINE P2UCNV

      INCLUDE 'p2.inc'

      UPPERC = .TRUE.
      LOWERC = .FALSE.

      RETURN
      END

C***********************************************************************
C Set the parser to convert the case of strings to lowercase
C RME 22/11/93 - 22/11/93
C***********************************************************************

      SUBROUTINE P2LCNV

      INCLUDE 'p2.inc'

      UPPERC = .FALSE.
      LOWERC = .TRUE.

      RETURN
      END

C***********************************************************************
C Set the prompt arrow. It stores it in the string PARROW up to a maxi-
C mum of 16 characters. This string should contain a space if desired.
C RME 18/7/93 - 18/7/93
C***********************************************************************

      SUBROUTINE P2AROW(ARROW)
      CHARACTER*(*) ARROW

      INCLUDE 'p2.inc'

      LARROW = MIN(LEN(PARROW),LEN(ARROW))
      PARROW = ARROW(:LARROW)

      RETURN
      END

