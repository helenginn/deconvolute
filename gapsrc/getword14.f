CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									       C
C	GETWORD.FOR Version 1.4 - A relatively robust input parsing routine    C
C									       C
C	Robert Esnouf, LMB Oxford, 20 December 1988			       C
C			  Updated: 26 January 1990			       C
C			Vers. 1.4: 2 July 1993 - the beginning of the end      C
C			         - 13 January 1994                             C
C									       C
C	Five recommended entry points are supported, to use other points       C
C	consult the comments blocks at the start of the relevant routines.     C
C									       C
C	1) CALL GETOPEN(' ',1)						       C
C									       C
C	   Subroutine to start the ball rolling by opening the initial file    C
C	   and initialising everything.					       C
C									       C
C	2) NOPT=KEYWORD(level,options,help)				       C
C	   level	CHAR Program level used in prompt generation etc.      C
C	   options	CHAR Contains the allowed options for input	       C
C	   help		INT  If not zero then typing "HELP" gives information  C
C									       C
C	   Integer function to return a number of the keyword typed in the     C
C	   option string. Options are separated by '!'s and can be of any      C
C	   length. eg: 'READ!WRITE!COPY!CLEAR!LIST!END'			       C
C									       C
C	2a)NOPT=KEYLIST(level,optarray,nopts,help)			       C
C	   level	CHAR Program level used in prompt generation etc.      C
C	   optarray	CHAR Contains the allowed options for input	       C
C	   nopts        INT  Number of options in character array	       C
C	   help		INT  If not zero then typing "HELP" gives information  C
C									       C
C	   Integer function to return a number of the keyword typed in the     C
C	   option array.						       C
C									       C
C	3) I=INTVAL(level,min,max,help)					       C
C	   level	CHAR Program level used in prompt generation etc.      C
C	   min		INT  Minimum allowed value			       C
C	   max		INT  Maximum allowed value			       C
C	   help		INT  If not zero then typing "HELP" gives information  C
C									       C
C	   Integer function to get an integer value within the specified       C
C	   range. Non-integer values and character strings are rejected. The   C
C	   help variable is not 0 if special help is available.		       C
C									       C
C	4) X=REALVAL(level,rmin,rmax,help)				       C
C	   level	CHAR Program level used in prompt generation etc.      C
C	   rmin		REAL Minimum allowed value			       C
C	   rmax		REAL Maximum allowed value			       C
C	   help		INT  If not zero then typing "HELP" gives information  C
C									       C
C	   Real function to get an real number value within the specified      C
C	   range. Character strings are rejected, though input is standard     C
C	   form is acceptable. If the help variable is not 0 if special help   C
C	   is available.						       C
C									       C
C	5) STR=CHARVAL(level,min,max,help)				       C
C	   STR=G2CHRL(level,min,max,help)				       C
C	   STR=G2CHRU(level,min,max,help)				       C
C	   level	CHAR Program level used in prompt generation etc.      C
C	   min		INT  Shortest permissible string		       C
C	   max		INT  Longest permissible string			       C
C	   help		INT  If not zero then typing "HELP" gives information  C
C									       C
C	   Character function to get string of characters within a specified   C
C	   length range. If the string needs to contain spaces then it should  C
C	   be enclosed in quotation marks. Quotation marks can be included in  C
C	   strings by using quote images (Pairs of quotation marks). If the    C
C	   help variable is not zero then special help is available.	       C
C	   Alternative forms return the string either upper or lowercased.     C
C									       C
C									       C
C	Further to these calls, a string of words to ignore can be defined by  C
C	CALL G2SIGN(string), where (string) is just a list of these words      C
C	separated by '!'s. For example CALL G2SIGN('!TO!FROM!THE!') would      C
C	cause the parser to seek another word from the input if any of those   C
C	four words were typed. Note that exclamation marks are also placed at  C
C	both ends of the string. When GETOPEN is called at the beginning of a  C
C	program it does CALL G2SIGN('!'), so no words are ignored. CALL        C
C	P2UCNV is used to force conversion of input to upper case. By default  C
C	no case conversion is performed. However, by default all comparisons   C
C	are case insensitive.  The P2 calls P2STCI and P2STCS control this.    C
C	CALL P2NCNV disables case conversion if it has been enabled.           C
C	The production of prompts is supressed (in a nested fashion) by typing C
C	#NOPROMPT or #NOECHO and enabled (in a nested fashion) by typing #ECHO C
C	or #PROMPT. Typing #SILENT causes unit 6 to be reassigned to NL:, thus C
C	supressing output. This does not affect the parser which writes to the C
C	unit passed to GETOPEN. Typing #VERBOSE behaves as you would expect by C
C	assigning unit 6 to be SYS$OUTPUT. This implementation is only simple  C
C	and may not be suitable for all situations. By default the parser runs C
C	in VERBOSE. If #INCLUDE is typed then the parser gets another word     C
C	from the input stream and tries to open a file with that name. This    C
C	process allows ready prepared files to be included at any place in a   C
C	program. When the file is finished, further input is got from the      C
C	original file. #INCLUDEs may be nested up to 29 files deep, and under  C
C	VMS the same file may be opened more than once. As an alternative to   C
C	#INCLUDE, the user may just type @filename, which is intended for use  C
C	with 'trusted' files as a #NOECHO is forced before the file is opened  C
C	and after the file closes a #ECHO is executed. If a command beginning  C
C	with a '#' is entered other than the ones above it is assumed to be an C
C	operating system command and so passed to it with a call to LIB$SPAWN. C
C	The ECHO/NOECHO nesting can be overridden by a call to GETECHO, which  C
C	starts echo. If the calls above do not satisfy requirements then a     C
C	CALL GETPARSE(level,word,ilen) can be used to develop your own style   C
C	of input, whilst retaining the basic functions of the parser. The      C
C	arguments are:							       C
C									       C
C			level	CHAR Used for prompts as above		       C
C			word	CHAR Word is returned in this variable after   C
C				     parsing is complete		       C
C			ilen	INT  Contains the length of the actual word    C
C				     returned in the 'word' argument	       C
C									       C
C	Variables may be assigned using the form $(varname)=(string). If the   C
C	string is surrounded by quotation marks then these are removed before  C
C	processing. (string) may contain variable references, ie $(varname)    C
C	and these will be substituted recursively (if the variable is defined) C
C	or left untouched. If a variable is referred to as $$(varname) then it C
C	is substituted when referenced rather than at definition. The contents C
C	of a variable can be examined by typing $(varname)= and a variable     C
C	can be undefined by typing $(varname)="" without any (string). All     C
C	variable references mentioned above can have a '$' sign trailing the   C
C	variable name, and this should be used when the intention would	       C
C	otherwise be obscure or impossible for the parser to interpret. For    C
C	example:							       C
C			$FILE="PROTEIN"					       C
C			$FILE1=$FILE1.PDB is ambiguous, however		       C
C			$FILE1=$FILE$1.PDB makes the message interpretable     C
C                                                                              C
C	When using files for input it may be that it is desirable to ask for   C
C	the value of a variable at run time. This is achieved by typing one    C
C	of the following:-						       C
C				$(var)=?				       C
C				$(var)=?prompt				       C
C				$(var)=?"lengthy prompt"		       C
C									       C
C	In the first case, if a prompt is required then the name of the        C
C	variable is used. In these commands, the input is taken from the next  C
C	word which appears in the previous input level, unless we are at the   C
C	top level, in which case the next word is used. If no word is ready    C
C	for input then a prompt is generated.				       C
C	The variables $DATE and $TIME are supported, and providing an easy,    C
C	machine independent method for getting the date and time. The routines C
C	that get this information are also available to the user. These are:   C
C									       C
C		CALL P2DATE(word,ilen)   Returns date in the character string  C
C					 (word) and the actual length of the   C
C					 string in (ilen).		       C
C		CALL P2TIME(word,ilen)   Returns time in the character string  C
C					 (word) and the actual length of the   C
C					 string in (ilen).		       C
C									       C
C	It would be tedious to have to code calls to frequently used routines  C
C	every time they may be required. An example of this would be calls to  C
C	read PDB files or list the names of proteins currently in memory. In   C
C	other cases, it would be desirable to have access to a subroutine from C
C	any input level, for example STOP to stop the program in a controlled  C
C	fashion checking that there is no data that needs saving. To allow for C
C	these possibilities a set of subroutines have been added to the parser C
C	version 1.3 to deal with utilities. These commands are checked for at  C
C	every input. They are initiated by CALL SETUTIL(string). Where the     C
C	argument is the commands in the same format as KEYWORD. By default     C
C	CALL GETOPEN sets the string to 'UTIL!STOP'. When one of these words   C
C	is detected, a call is made to a routine UTILITIES(string) which is    C
C	supplied by the programmer, although a minimal default one is part of  C
C	the parser. (string) is just the 4 letter name of the utility required C
C	and so this subroutine just needs to perform a series of statements    C
C	of the form:							       C
C									       C
C			IF (UTIL.EQ.'READ') CALL READPDB		       C
C									       C
C	The user then has to write the necessary subroutines to perform the    C
C	utility. The functions 2 - 5 above are available within the utility    C
C	though they have to be called as KEYWORD2, INTVAL2 etc to stop the     C
C	possibility of recursive calling. These functions behave exactly the   C
C	same way as above except that the utility commands are not checked for C
C	by the input. To list the commands currently available as utilities    C
C	then CALL UTIL will do this in much the same format as HELP. As the    C
C	utilities are checked for before 'HELP' then defining HELP as one of   C
C	the utilities will allow the programmer to provide their own style of  C
C	help within a program. Finally, CALL GETUTIL(word,ilen) will return    C
C	will return the current string of utilities in the string (word) with  C
C	length (ilen). This can be used to temporarily disable some utulities  C
C	or to add to the current list even if the current list is unknown.     C
C									       C
C	As an update, version 1.31 allows for the '=' sign to be treated as    C
C	a space. This is switched on by CALL EQUALS and off by CALL NOEQUALS.  C
C	By default, GETOPEN calls NOEQUALS. The purpose of this is to allow    C
C	for routines that set parameters to look like SET WEIGHT=2.0 END       C
C	which looks nice. Routines to set parameters probably ought to be      C
C	utility routines so that they are available at all times and the use   C
C	of '=' as a space will not affect any other part of the program.       C
C									       C
C	To use the help facility, the programmer must include SUBROUTINE       C
C	G2HELP(help,opt) in his/her program. The arguments are the integer     C
C	variable passed to one of the functions 2 - 5 and the second is the    C
C	number of the command of KEYWORD calls. This subroutine is called only C
C	if 'help' is ever passed to one of the functions with a non zero value C
C	but it should always be present in at least its minimal form in order  C
C	to prevent LINK warnings. A suitable form for the help function is     C
C	shown below. Note that if G2HELP needs to use the parser at all it     C
C	must use call to G2PRS1(level,word,ilen) otherwise FORTRAN will        C
C	probably object to being called recursively!			       C
C									       C
C	SUBROUTINE G2HELP(HELP,OPT)					       C
C	INTEGER		HELP,OPT					       C
C									       C
C	WRITE(6,*)							       C
C									       C
C	IF (HELP.EQ.1) THEN						       C
C	  WRITE(6,*) 'Lots of useful help information can make programs'       C
C	  WRITE(6,*) 'so much easier to use.'				       C
C	ELSE IF (HELP.EQ.2.AND.OPT.EQ.1) THEN				       C
C	  WRITE(6,*) 'If you want more help then Robert''s the man to ask'     C
C	ELSE IF (HELP.EQ.2.AND.OPT.EQ.2) THEN				       C
C	  WRITE(6,*) 'You can put as many different help sections as you'      C
C	  WRITE(6,*) 'like. Even have an interactive help a bit like VMS'      C
C	  WRITE(6,*) 'and better still the parser can be used to handle'       C
C	  WRITE(6,*) 'the input!'					       C
C	ELSE								       C
C	  WRITE(6,*) 'Remember to trap for unsupported values'		       C
C	ENDIF								       C
C									       C
C	WRITE(6,*)							       C
C									       C
C	RETURN								       C
C	END								       C
C									       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine to open files on the appropriate units for input - if unit is
C positive then this is the start of a new set of files, internal calls usually
C pass zero. When a file is opened using an '@', -1 is passed. This has the same
C effect as passing zero except that the string '#ECHO ' is appended onto the
C beginning of the current line. This is to make the '@filename' form execute
C the file without echoing anything to the screen. When called with a positive
C value of I, then the opened file is on unit (I+1) and unit I is used for
C output of messages by the parser. This allows the user to reassign unit 6
C to NL: or otherwise and not affect the parser. The idea intended is to have
C subroutines which can run in SILENT or VERBOSE modes.

	SUBROUTINE GETOPEN(FILENAME, I)
	CHARACTER*(*)	FILENAME
	INTEGER		I

	INCLUDE 'gw14.inc'

	LOGICAL		G2OPEN, G2INIT
        INTEGER		NUN, J

C Save the remaining line from the old file and then set the new input unit. If
C the nesting gets too deep then don't include it.

	IF (I.LE.0) THEN
	  LINEBAK(NINP-NST+1) = LINE
	  NINP = NINP + 1
	  IF (NINP.GT.99 .OR. NINP-NST.GT.29) THEN
	    WRITE(NOUT, 1000)
	    NINP = NINP - 1
	    RETURN
	  ENDIF
	  ECHO(NINP-NST+1) = ECHO(NINP-NST) - I
	  IF (.NOT.G2OPEN(NINP,FILENAME)) GOTO 10
	  LINE = ' '
	ELSE

C This is the initialisation call, it does not use either the filename
C of the unit number, relying simply on P2INIT to know what to do

	  CALL P2UCNV
	  CALL G2ECHO
	  CALL G2NOEQ
	  CALL G2PSHI
	  CALL G2SIGN('!')
	  CALL G2SUT('UTIL!STOP')
	  NVAR = 0
	  NUSED = 0
	  NINP = 1
	  NST = 1
          NUN = NINP
	  IF (.NOT.G2INIT(NUN)) GOTO 10
          NOUT = NUN
	  LINE = ' '
	ENDIF

	RETURN

C File could not be opened, so report the error and carry on without including
C the file

10	NINP = NINP - 1
	DO J = 1, 40
	  IF (FILENAME(J:).EQ.' ') GOTO 20
	ENDDO
20	WRITE(NOUT, 1010) FILENAME(1:J-1)
	RETURN

1000	FORMAT('GETOPEN: *** Nesting of files too deep ***')
1010	FORMAT('GETOPEN: *** Cannot open ',A,' ***')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Integer function to get an option from a passed option list. Invalid values
C are trapped and an error message produced incorporating the first string
C passed which indicates the level from which the subroutine was called. The
C options are passed in the second string as commands separated by '!'s. The
C number of significant characters is arbitrary and abbreviations are
C accepted.

	INTEGER FUNCTION KEYWORD(LEVEL, OPTIONS, HELP)
	CHARACTER*(*)	LEVEL, OPTIONS
	INTEGER		HELP

	INCLUDE 'gw14.inc'

	INTEGER		G2KWRD

	INHELP = .FALSE.
10	UTDO = .FALSE.
	UTON = .TRUE.
	KEYWORD = G2KWRD(LEVEL, OPTIONS, HELP)
	UTON = .FALSE.
	IF (UTDO) THEN
	  UTDO = .FALSE.
	  CALL G2UTLS(UTILITY)
	  GOTO 10
	ENDIF

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Integer function to get an option from an array. Very similar to the KEYWORD
C option below. See there for details. Exact (case insensitive) matches only in
C this current version.

	INTEGER FUNCTION KEYLIST(LEVEL, OPTLIST, NOPTS, HELP)
	INTEGER		NOPTS,HELP
	CHARACTER*(*)	LEVEL, OPTLIST(NOPTS)

	INCLUDE 'gw14.inc'

	INTEGER		G2KLST

	INHELP = .FALSE.
10	UTDO = .FALSE.
	UTON = .TRUE.
	KEYLIST = G2KLST(LEVEL, OPTLIST, NOPTS, HELP)
	UTON = .FALSE.
	IF (UTDO) THEN
	  UTDO = .FALSE.
	  CALL G2UTLS(UTILITY)
	  GOTO 10
	ENDIF

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Integer function to extract an integer value between specified ranges from
C the next word in the input stream.

	INTEGER FUNCTION INTVAL(LEVEL, MIN, MAX, HELP)
	CHARACTER*(*)	LEVEL
	INTEGER		MIN, MAX, HELP

	INCLUDE 'gw14.inc'

	INTEGER		G2IVAL

10	UTDO = .FALSE.
	UTON = .TRUE.
	INTVAL = G2IVAL(LEVEL, MIN, MAX, HELP)
	UTON = .FALSE.
	IF (UTDO) THEN
	  UTDO = .FALSE.
	  CALL G2UTLS(UTILITY)
	  GOTO 10
	ENDIF

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Real function to extract an real value between specified ranges from the next
C word in the input stream.

	REAL FUNCTION REALVAL(LEVEL, RMIN, RMAX, HELP)
	CHARACTER*(*)	LEVEL
	REAL		RMIN, RMAX
	INTEGER		HELP

	INCLUDE 'gw14.inc'

	REAL		G2RVAL

10	UTDO = .FALSE.
	UTON = .TRUE.
	REALVAL = G2RVAL(LEVEL, RMIN, RMAX, HELP)
	UTON = .FALSE.
	IF (UTDO) THEN
	  UTDO = .FALSE.
	  CALL G2UTLS(UTILITY)
	  GOTO 10
	ENDIF

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Character function to take the next word from the input stream as a character
C string whose length is between the specified limits.

	FUNCTION CHARVAL(LEVEL, MIN, MAX, HELP)
	CHARACTER*(*)	CHARVAL,LEVEL
	INTEGER		MIN, MAX, HELP

	INCLUDE 'gw14.inc'

	CHARACTER*256	G2CVAL

10	UTDO = .FALSE.
	UTON = .TRUE.
	CHARVAL = G2CVAL(LEVEL, MIN, MAX, HELP)
	UTON = .FALSE.
	IF (UTDO) THEN
	  UTDO = .FALSE.
	  CALL G2UTLS(UTILITY)
	  GOTO 10
	ENDIF

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Same call and KEYWORD but doesn't act on utilities

	INTEGER FUNCTION G2KWRD(LEVEL, OPTIONS, HELP)
	CHARACTER*(*)	LEVEL,OPTIONS
	INTEGER		HELP

	INCLUDE 'gw14.inc'

	CHARACTER*80	WORD
	INTEGER		ILEN, IX

	CHARACTER*75	OUTP
	LOGICAL		P2CMP
	INTEGER		P2NOPT

C Check to see if we're returning from an interrupted HELP section

        G2KWRD = 0
C RME 19-4-2010 Modified for F2003
C	IF (INHELP) GOTO 30
	IF (INHELP) GOTO 12

C Get a word and look for it in the options string, if the end of file is found
C then stop the program with a useful message

10	CALL G2PRS1(LEVEL, WORD, ILEN)
	IF (UTDO) RETURN

C Check for the end of the input stream

	IF (ILEN.EQ.0) THEN
	  WRITE(NOUT, 1000) LEVEL
	  CALL M2EXIT
	ENDIF

C If HELP is typed and a help option is available then give the help. The help
C is implemented in a subroutine G2HELP, where the passed parameter is just
C the value of the HELP variable. It is intended that G2HELP(help,opt) is
C implemented in the user program as a subroutine which just does a series
C of write statements to unit 6 depending on the values of (help) and (opt).
C If HELP is typed without a help option being available, or a '?' is typed then
C give a command summary for the available level - ie just a list of the words
C in OPTIONS.

C RME 19-4-2010 Modified for F2003
C	IF (P2CMP(WORD(:ILEN),'HELP') .OR. WORD(:ILEN).EQ.'?') THEN
12      IF (P2CMP(WORD(:ILEN),'HELP') .OR. WORD(:ILEN).EQ.'?'
     +      .OR. INHELP) THEN
C RME 19-4-2010 Added for F2003
          IF (INHELP) GOTO 30
15	  WRITE(NOUT, *)
	  WRITE(NOUT, 1020)
	  WRITE(NOUT, *)
	  CALL P2DOPT(OPTIONS, OUTP, '(4X,A)', 8)
	  WRITE(NOUT, *)
	  IF (WORD(:ILEN).EQ.'?' .OR. HELP.EQ.0) GOTO 10

C Now get the option if help is available

30	  INHELP = .TRUE.
	  CALL G2PRS1('Enter command or type Q', WORD, ILEN)
	  IF (UTDO) RETURN
	  INHELP = .FALSE.
	  IF (ILEN.EQ.0) THEN
	    WRITE(NOUT, 1000) LEVEL
	    CALL M2EXIT
	  ENDIF
	  IF (P2CMP(WORD(:ILEN),'Q')) THEN
	    WRITE(NOUT, *)
	    GOTO 10
	  ENDIF
	  IF (WORD(:ILEN).EQ.'?') THEN
	    WORD = 'HELP'
	    ILEN = 4
	  ENDIF
	  IF (P2CMP(WORD(:ILEN),'HELP')) GOTO 15

C Now try and match the word - accepting shortened forms providing that they
C are not ambiguous. Whizzo version (RME 22/11/92)

	  IX = P2NOPT(OPTIONS, WORD(:ILEN))
	  IF (IX.EQ.-1) THEN
	    WRITE(NOUT,1005) 'HELP', WORD(:ILEN)
	    GOTO 30
	  ENDIF

	  IF (IX.EQ.0) THEN
	    WRITE(NOUT,1010) 'HELP', WORD(:ILEN)
	    GOTO 30
	  ELSE
	    CALL G2HELP(HELP,IX)
	    GOTO 30
	  ENDIF
	ENDIF

C Now try and match the word - accepting shortened forms providing that they
C are not ambiguous. Whizzo version (RME 22/11/92)

	IX = P2NOPT(OPTIONS, WORD(:ILEN))
	IF (IX.EQ.-1) THEN
	  WRITE(NOUT,1005) LEVEL, WORD(:ILEN)
	  GOTO 10
	ENDIF

C If the option is not found then IX=0, so print a warning message and get a new
C word, otherwise return IX

	IF (IX.EQ.0) THEN
	  WRITE(NOUT, 1010) LEVEL, WORD(1:ILEN)
	  GOTO 10
	ELSE
	  G2KWRD = IX
	  RETURN
	ENDIF

1000	FORMAT(A,': *** End of input -- PROGRAM TERMINATED ***')
1005	FORMAT(A,': *** Keyword "',A,'" ambiguous ***')
1010	FORMAT(A,': *** Keyword "',A,'" ignored ***')
1020	FORMAT('Available commands at this level:-')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Same call as KEYLIST but doesn't act on utilities

	INTEGER FUNCTION G2KLST(LEVEL, OPTLIST, NOPTS, HELP)
	INTEGER		NOPTS, HELP
	CHARACTER*(*)	LEVEL, OPTLIST(NOPTS)

	INCLUDE 'gw14.inc'

	CHARACTER*80	WORD
	INTEGER		ILEN, IX, I

	CHARACTER*75	OUTP
	LOGICAL		P2CMP

C Check to see if we're returning from an interrupted HELP section

        G2KLST = 0
C RME 19-4-2010 Modified for F2003
C	IF (INHELP) GOTO 30
	IF (INHELP) GOTO 12

C Get a word and look for it in the options string, if the end of file is found
C then stop the program with a useful message

10	CALL G2PRS1(LEVEL, WORD, ILEN)
	IF (UTDO) RETURN

C Check for the end of the input stream

	IF (ILEN.EQ.0) THEN
	  WRITE(NOUT, 1000) LEVEL
	  CALL M2EXIT
	ENDIF

C If HELP is typed and a help option is available then give the help. The help
C is implemented in a subroutine G2HELP, where the passed parameter is just
C the value of the HELP variable. It is intended that G2HELP(help,opt) is
C implemented in the user program as a subroutine which just does a series
C of write statements to unit 6 depending on the values of (help) and (opt).
C If HELP is typed without a help option being available, or a '?' is typed then
C give a command summary for the available level - ie just a list of the words
C in OPTIONS.

C RME 19-4-2010 Modified for F2003
C	IF (P2CMP(WORD(:ILEN),'HELP') .OR. WORD(:ILEN).EQ.'?') THEN
12      IF (P2CMP(WORD(:ILEN),'HELP') .OR. WORD(:ILEN).EQ.'?'
     +      .OR. INHELP) THEN
C RME 19-4-2010 Added for F2003
          IF (INHELP) GOTO 30
15	  WRITE(NOUT, *)
	  WRITE(NOUT, 1020)
	  WRITE(NOUT, *)
	  CALL P2DARY(OPTLIST, NOPTS, OUTP, '(4X,A)', 8)
	  WRITE(NOUT, *)
	  IF (WORD(:ILEN).EQ.'?' .OR. HELP.EQ.0) GOTO 10

C Now get the option if help is available

30	  INHELP = .TRUE.
	  CALL G2PRS1('Enter command or type Q', WORD, ILEN)
	  IF (UTDO) RETURN
	  INHELP = .FALSE.
	  IF (ILEN.EQ.0) THEN
	    WRITE(NOUT, 1000) LEVEL
	    CALL M2EXIT
	  ENDIF
	  IF (P2CMP(WORD(:ILEN),'Q')) THEN
	    WRITE(NOUT, *)
	    GOTO 10
	  ENDIF
	  IF (WORD(:ILEN).EQ.'?') THEN
	    WORD = 'HELP'
	    ILEN = 4
	  ENDIF
	  IF (P2CMP(WORD(:ILEN),'HELP')) GOTO 15

C Now try and match the word

	  IX = 0
	  DO I = 1, NOPTS
	    IF (P2CMP(OPTLIST(I),WORD(:ILEN))) THEN
	      IX = I
	      GOTO 45
	    ENDIF
	  ENDDO

45	  IF (IX.EQ.-1) THEN
	    WRITE(NOUT, 1005) 'HELP', WORD(:ILEN)
	    GOTO 30
	  ENDIF

	  IF (IX.EQ.0) THEN
	    WRITE(NOUT, 1010) 'HELP', WORD(:ILEN)
	    GOTO 30
	  ELSE
	    CALL G2HELP(HELP,IX)
	    GOTO 30
	  ENDIF
	ENDIF

C Now try and match the word

	IX = 0
	DO I = 1, NOPTS
	  IF (P2CMP(OPTLIST(I),WORD(:ILEN))) THEN
	    IX = I
	    GOTO 40
	  ENDIF
	ENDDO

40	IF (IX.EQ.-1) THEN
	  WRITE(NOUT,1005) LEVEL, WORD(:ILEN)
	  GOTO 10
	ENDIF

C If the option is not found then IX=0, so print a warning message and get a new
C word, otherwise return IX

	IF (IX.EQ.0) THEN
	  WRITE(NOUT, 1010) LEVEL, WORD(:ILEN)
	  GOTO 10
	ELSE
	  G2KLST = IX
	  RETURN
	ENDIF

1000	FORMAT(A,': *** End of input -- PROGRAM TERMINATED ***')
1005	FORMAT(A,': *** Keyword "',A,'" ambiguous ***')
1010	FORMAT(A,': *** Keyword "',A,'" ignored ***')
1020	FORMAT('Available commands at this level:-')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Same call as INTVAL but doesn't act on utilities

	INTEGER FUNCTION G2IVAL(LEVEL, MIN, MAX, HELP)
	CHARACTER*(*)	LEVEL
	INTEGER		MIN,MAX,HELP

	INCLUDE 'gw14.inc'

	CHARACTER*80	WORD,LLIM,ULIM
	INTEGER		ILEN,ILL,IUL
	REAL		RVAL

	G2IVAL = 0
10	CALL G2PRS1(LEVEL, WORD, ILEN)
	IF (UTDO) RETURN
	IF (ILEN.EQ.0) THEN
	  WRITE(NOUT, 1000) LEVEL
	  CALL M2EXIT
	ENDIF

C Check for HELP or command summary ('?') being asked for

	IF (WORD.EQ.'HELP' .AND. HELP.NE.0) THEN
	  CALL G2HELP(HELP, 0)
	  GOTO 10
	ENDIF

	IF (WORD.EQ.'HELP' .OR. WORD.EQ.'?') THEN
	  WRITE(NOUT, *)
	  LLIM = ' '
	  WRITE(LLIM, *) MIN
	  CALL G2NORM(LLIM, ILL)
	  ULIM = ' '
	  WRITE(ULIM, *) MAX
	  CALL G2NORM(ULIM, IUL)
	  WRITE(NOUT, 1010) LLIM(:ILL), ULIM(:IUL)
	  WRITE(NOUT, *)
	  GOTO 10
	ENDIF

C Now check the number is integral and decode it into an integer

	READ(WORD, *, ERR=20) G2IVAL
	READ(WORD, *, ERR=20) RVAL
	IF (RVAL.NE.FLOAT(INT(RVAL))) THEN
	  WRITE(NOUT, 1020) LEVEL, WORD(:ILEN)
	  GOTO 10
	ENDIF

C Check the limits

	IF (G2IVAL.LT.MIN) THEN
	  WRITE(NOUT, 1030) LEVEL, WORD(:ILEN)
	ELSE IF (G2IVAL.GT.MAX) THEN
	  WRITE(NOUT, 1040) LEVEL, WORD(:ILEN)
	ELSE
	  RETURN
	ENDIF
	GOTO 10

C String can not be interpreted as an integer

20	WRITE(NOUT, 1050) LEVEL, WORD(:ILEN)
	GOTO 10

1000	FORMAT(A,': *** End of input -- PROGRAM TERMINATED ***')
1010	FORMAT('Enter an integer between ',A,' and ',A)
1020	FORMAT(A,': *** ',A,' is non integral ***')
1030	FORMAT(A,': *** ',A,' is too small ***')
1040	FORMAT(A,': *** ',A,' is too large ***')
1050	FORMAT(A,': *** Cannot read ',A,' as an integer ***')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Same call as REALVAL but doesn't act on utilities

	REAL FUNCTION G2RVAL(LEVEL, RMIN, RMAX, HELP)
	CHARACTER*(*)	LEVEL
	REAL		RMIN,RMAX
	INTEGER		HELP

	INCLUDE 'gw14.inc'

	CHARACTER*80	WORD, LLIM, ULIM
	INTEGER		ILEN,ILL,IUL

	G2RVAL = 0.0
10	CALL G2PRS1(LEVEL, WORD, ILEN)
	IF (UTDO) RETURN
	IF (ILEN.EQ.0) THEN
	  WRITE(NOUT, 1000) LEVEL
	  CALL M2EXIT
	ENDIF

C Check for HELP or command summary ('?') being asked for

	IF (WORD.EQ.'HELP' .AND. HELP.NE.0) THEN
	  CALL G2HELP(HELP, 0)
	  GOTO 10
	ENDIF

	IF (WORD.EQ.'HELP' .OR. WORD.EQ.'?') THEN
	  WRITE(NOUT, *)
	  LLIM = ' '
	  IF (RMIN.EQ.0.0) THEN
	    LLIM = '0.0'
	  ELSE
	    WRITE(LLIM, *) RMIN
	  ENDIF
	  CALL G2NORM(LLIM, ILL)
	  ULIM = ' '
	  IF (RMAX.EQ.0.0) THEN
	    ULIM = '0.0'
	  ELSE
	    WRITE(ULIM, *) RMAX
	  ENDIF
	  CALL G2NORM(ULIM, IUL)
	  WRITE(NOUT, 1010) LLIM(:ILL), ULIM(:IUL)
	  WRITE(NOUT, *)
	  GOTO 10
	ENDIF

C Now check word can be read as a number and decode it

	READ(WORD,*,ERR=20) G2RVAL

C Check the limits

	IF (G2RVAL.LT.RMIN) THEN
	  WRITE(NOUT, 1030) LEVEL, WORD(:ILEN)
	ELSE IF (G2RVAL.GT.RMAX) THEN
	  WRITE(NOUT,1 040) LEVEL, WORD(:ILEN)
	ELSE
	  RETURN
	ENDIF
	GOTO 10

C String can not be interpreted as a number

20	WRITE(NOUT, 1050) LEVEL, WORD(:ILEN)
	GOTO 10

1000	FORMAT(A,': *** End of input -- PROGRAM TERMINATED ***')
1010	FORMAT('Enter a real number between ',A,' and ',A)
1030	FORMAT(A,': *** ',A,' is too small ***')
1040	FORMAT(A,': *** ',A,' is too large ***')
1050	FORMAT(A,': *** Cannot read ',A,' as a real number ***')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Same call as CHARVAL but doesn't act on utilities

	FUNCTION G2CVAL(LEVEL, MIN, MAX, HELP)
	CHARACTER*(*)	G2CVAL, LEVEL
	INTEGER		MIN, MAX, HELP

	INCLUDE 'gw14.inc'

	CHARACTER*256	WORD, LLIM, ULIM
	INTEGER		ILEN, ILL, IUL

	G2CVAL = ' '
10	CALL G2PRS1(LEVEL, WORD, ILEN)
	IF (UTDO) RETURN
	IF (ILEN.EQ.0) THEN
	  WRITE(NOUT, 1000) LEVEL
	  CALL M2EXIT
	ENDIF

C Check for HELP or command summary ('?') being asked for

	IF (WORD.EQ.'HELP' .AND. HELP.NE.0) THEN
	  CALL G2HELP(HELP, 0)
	  GOTO 10
	ENDIF

	IF (WORD.EQ.'HELP' .OR. WORD.EQ.'?') THEN
	  WRITE(NOUT, *)
	  LLIM = ' '
	  WRITE(LLIM, *) MIN
	  CALL G2NORM(LLIM, ILL)
	  ULIM = ' '
	  WRITE(ULIM, *) MAX
	  CALL G2NORM(ULIM, IUL)
	  WRITE(NOUT, 1010) LLIM(:ILL), ULIM(:IUL)
	  WRITE(NOUT, *)
	  GOTO 10
	ENDIF

C Check the length of the string

	IF (ILEN.LT.MIN) THEN
	  WRITE(NOUT, 1030) LEVEL, WORD(:ILEN)
	ELSE IF (ILEN.GT.MAX) THEN
	  WRITE(NOUT, 1040) LEVEL, WORD(:ILEN)
	ELSE
	  G2CVAL = WORD(:ILEN)
	  RETURN
	ENDIF
	GOTO 10

1000	FORMAT(A,': *** End of input -- PROGRAM TERMINATED ***')
1010	FORMAT('Enter a string of between ',A,' and ',A,' characters')
1030	FORMAT(A,': *** ',A,' is too short ***')
1040	FORMAT(A,': *** ',A,' is too long ***')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Character function to take the next word from the input stream as a character
C string whose length is between the specified limits. Resulting string is
C lowercase

	FUNCTION G2CHRL(LEVEL, MIN, MAX, HELP)
	CHARACTER*(*)	G2CHRL, LEVEL
	INTEGER		MIN, MAX, HELP

	CHARACTER*256	CHARVAL	

	G2CHRL = CHARVAL(LEVEL, MIN, MAX, HELP)
	CALL P2LCAS(G2CHRL)

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Character function to take the next word from the input stream as a character
C string whose length is between the specified limits. Resulting string is
C uppercase

	FUNCTION G2CHRU(LEVEL, MIN, MAX, HELP)
	CHARACTER*(*)	G2CHRU, LEVEL
	INTEGER		MIN, MAX, HELP

	CHARACTER*256	CHARVAL	

	G2CHRU = CHARVAL(LEVEL, MIN, MAX, HELP)
	CALL P2UCAS(G2CHRU)

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Utility subroutine to take a string containing a number and left justify it.
C It also calculates the number of significant characters in the string.

	SUBROUTINE G2NORM(WORD, ILEN)
	CHARACTER*(*)	WORD
	INTEGER		ILEN

	INCLUDE 'gw14.inc'

	INTEGER		I, IX, ILEN2
	CHARACTER*80	TPWORD

        ILEN2 = 0

C Left justify the string

	DO I = 1, LEN(WORD)
	  IF (WORD(I:I).NE.' ') GOTO 10
	ENDDO
10	WORD = WORD(I:)

C Now find its length

	DO I = 1, LEN(WORD)
	  IF (WORD(I:I).EQ.' ') GOTO 20
	ENDDO
20	ILEN = I - 1

C Find if its a number in exponent form and find the length of the exponent

	TPWORD = ' '
	IX = INDEX(WORD(1:ILEN), 'E')
	IF (IX.NE.0) THEN
	  TPWORD = WORD(IX:ILEN)
	  ILEN2 = ILEN - IX + 1
	  ILEN = IX - 1
	ENDIF

C Now see if there is a decimal point, and if so then remove spurious 0s

	IF (INDEX(WORD(:ILEN),'.').NE.0) THEN
	  DO I = ILEN, 1, -1
	    IF (WORD(I:I).NE.'0') GOTO 30
	  ENDDO
30	  IF (WORD(I:I).EQ.'.') THEN
	    WORD(I+1:I+1) = '0'
	    I = I + 1
	  ENDIF
	  ILEN = I
	ENDIF

C Now put the exponent back on if necessary

	IF (IX.NE.0) THEN
	  WORD = WORD(:ILEN)//TPWORD(:ILEN2)
	  ILEN = ILEN + ILEN2
	ENDIF

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine initialise the word saving mechanism of the parser
 
	SUBROUTINE G2PSHI

	INCLUDE 'gw14.inc'

	IXPSH = 0
	INEXT = IXPSH
	PSHFST = .TRUE.
	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine to tell the next call to GETPARSE to return the previous word on
C the stack. It may be called up to 30 times in a row to repeat the previous
C n words which are stoed on a rotating LIFO stack.
 
	SUBROUTINE G2PUSH

	INCLUDE 'gw14.inc'

C Force the next call to GETPARSE to look one further back in the stack

	IF (INEXT.EQ.IXPSH+1 .OR. (INEXT.EQ.1.AND.IXPSH.EQ.30) .OR.
     &      (INEXT.EQ.1.AND.PSHFST)) THEN
	  WRITE(NOUT, 1000)
	ELSE IF (INEXT.GT.1) THEN
	  INEXT = INEXT - 1
	ELSE
	  INEXT = 30
	ENDIF

	RETURN

1000	FORMAT('G2PUSH: *** Stack limit reached, PUSH ignored ***')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine to get the next word from the input. The word has been fully
C parsed at this stage: files have been substituted and variables have been
C expanded. All the parsing is done by the routine G2PRS2 below - this
C routine just serves to create a push back functionality of 30 words
 
	SUBROUTINE G2PRS1(LEVEL, WORD, ILEN)
	CHARACTER*(*)	LEVEL, WORD
	INTEGER		ILEN

	INCLUDE 'gw14.inc'

C If there are pushed back words then get the next one off the stack in a LIFO
C order

	IF (IXPSH.EQ.INEXT) THEN
	  CALL G2PRS2(LEVEL, WORD, ILEN)
	  IXPSH = IXPSH + 1
	  IF (IXPSH.GT.30) THEN
            IXPSH = 1
	    PSHFST = .FALSE.
	  ENDIF
	  INEXT = IXPSH
	  PSHWRD(IXPSH) = WORD
	  PSHLEN(IXPSH) = ILEN
	ELSE
	  INEXT = INEXT + 1
	  IF (INEXT.GT.30) INEXT = 1
	  WORD = PSHWRD(INEXT)
	  ILEN = PSHLEN(INEXT)
	ENDIF

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine to parse fully an input file word - #INCLUDE causes the next word
C to be treated as a filename and includes that file. Files can be nested 30
C deep. By default, the word is converted to uppercase - this is switched off
C by a call to GETNOCONV and restored by GETCONV. Prompts may also be enabled
C and disabled in this routine by the commands #PROMPT and #NOPROMPT which may
C be deeply nested. (The alternative commands #ECHO and #NOECHO may also be
C used)

	SUBROUTINE G2PRS2(LEVEL, WORD, ILEN)
	CHARACTER*(*)	LEVEL, WORD
	INTEGER		ILEN

	INCLUDE 'gw14.inc'

	INTEGER		ILEN1, VLEN, VVLEN
	INTEGER		I, ISTART, IDOLL
        LOGICAL		G2VCHR, P2CMP
	CHARACTER*256	TPWORD, VARNAM, VARVAL
	CHARACTER*1	C1

C Get the next word with the low level call

10	CALL G2WORD(LEVEL, WORD, ILEN)

C Variable names start with a '$' symbol, and may be placed anywhere in the
C word. The part of the string following the '$' is then extracted and passed
C to the G2GVAR routine, which returns a string for substitution in the original
C word. If the variable is not found then no substitution is made so that '$'s
C should be able to be used in strings without hassle. The first stage is to see
C if there is a dollar in the string.

	ISTART = 1
11	IDOLL = INDEX(WORD(ISTART:ILEN), '$')
	IF (IDOLL.NE.0) THEN
	  IDOLL = ISTART + IDOLL - 1

C Get the variable name its length

	  VARNAM = WORD(IDOLL:)
	  DO I = 2, 80
	    C1 = VARNAM(I:I)
	    IF (.NOT.G2VCHR(C1)) GOTO 12
	  ENDDO
12	  IF (C1.EQ.'$') THEN
	    WORD(IDOLL+I-1:ILEN) = WORD(IDOLL+I:ILEN)//' '
	    ILEN = ILEN - 1
	    VARNAM = WORD(IDOLL:)
	    C1 = VARNAM(I:I)
	  ENDIF
	  VLEN = I - 1

C Variables are defined by following the variable name with an equals sign. Also
C the '$' must be the first character of the string

	  IF (IDOLL.EQ.1 .AND. C1.EQ.'=') THEN
	    CALL G2SVAR(WORD, ILEN)
	    GOTO 10
	  ENDIF

C Go and look for the variable name in the current list

	  CALL G2GVAR(VARNAM(:VLEN), VARVAL, VVLEN)

C If the variable is found then make the substitution

	  IF (VVLEN.NE.0) THEN
	    IF (IDOLL.EQ.1) THEN
	      WORD = VARVAL(:VVLEN)//WORD(VLEN+1:ILEN)
	    ELSE
	      WORD = WORD(:IDOLL-1)//VARVAL(:VVLEN)//WORD(IDOLL+VLEN:ILEN)
	    ENDIF
	    ILEN = ILEN + VVLEN - VLEN
	    ISTART = IDOLL
	  ELSE
	    ISTART = IDOLL + 1
	  ENDIF

C And now go and check if there is another variable to be substituted

	  GOTO 11

	ENDIF

C Now process special directives in the input file. If the word is '#INCLUDE'
C then the next word is taken as a file name and the file is opened. Up to 30
C files may be opened in a nested fashion. Also supported is the form @filename
C for an included file. This is intended for trusted files, and incorporates a
C #NOECHO before doing the include and a #ECHO after the included file has
C completed executing. There should be no space between the @ and the filename.
C If the commands #PROMPT, #ECHO, #NOPROMPT or #NOECHO are typed then update
C the counter showing the degree of nesting of these calls.

	IF (WORD(1:1).EQ.'#') THEN
	  IF (P2CMP(WORD(2:), 'INCLUDE')) THEN
	    CALL G2WORD('Filename', WORD, ILEN)
	    CALL GETOPEN(WORD, 0)
	  ELSE IF (P2CMP(WORD(2:),'ECHO') .OR. P2CMP(WORD(2:),'PROMPT')) THEN
	    IF (ECHO(NINP-NST+1).GT.0) ECHO(NINP-NST+1) = ECHO(NINP-NST+1) - 1
	  ELSE IF (P2CMP(WORD(2:),'NOECHO').OR.P2CMP(WORD(2:),'NOPROMPT')) THEN
	    ECHO(NINP-NST+1) = ECHO(NINP-NST+1) + 1
	  ELSE IF (P2CMP(WORD(2:),'TIME')) THEN
	    CALL P2PCPU
	  ELSE IF (P2CMP(WORD(2:),'SILENT')) THEN
	    CALL P2TERS
	  ELSE IF (P2CMP(WORD(2:),'TERSE')) THEN
	    CALL P2TALK
	  ELSE IF (P2CMP(WORD(2:),'VERBOSE')) THEN
	    CALL P2CHAT
	  ELSE IF (P2CMP(WORD(2:),'FILES')) THEN
            CALL P2INQA
	    CALL P2FLST
	  ELSE
	    CALL P2SHEL(WORD(2:))
	  ENDIF
	  GOTO 10
	ELSE IF (WORD(1:1).EQ.'@') THEN
	  CALL GETOPEN(WORD(2:), -1)
	  GOTO 10
	ENDIF

C And case convert into a temporary string

	ILEN1 = ILEN
        TPWORD = WORD
        CALL P2UCAS(TPWORD(:ILEN1))

C Check if the word is in the IGNORE list, and if so then go right back
C to the beginning of the routine and do it again!

	IF (INDEX(IGNORE(1:IGLEN),'!'//TPWORD(1:ILEN1)//'!').NE.0) GOTO 10

C Check if the word is one of the utility commands, only if we're at an input
C that allows it

	IF (UTON) THEN
	  IF (INDEX(UTILS(:UTLEN), TPWORD(:4)) .NE. 0) THEN
	    UTDO = .TRUE.
	    UTILITY = TPWORD(:4)
	    RETURN
	  ENDIF
	ENDIF

C Finally, if conversion of the return string is required, then do it

        CALL P2CNV(WORD(:ILEN))
	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine to set parser to treat '=' as a space

	SUBROUTINE G2EQ

	INCLUDE 'gw14.inc'

	EQUAL = .TRUE.
	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine to set parser to treat '=' as an ordinary character

	SUBROUTINE G2NOEQ

	INCLUDE 'gw14.inc'

	EQUAL = .FALSE.
	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine to initialise printing of prompts - overrides the nesting feature
C of the #NOECHO command

	SUBROUTINE G2ECHO

	INCLUDE 'gw14.inc'

	INTEGER		I

	DO I = 1, 30
	  ECHO(I) = 0
	ENDDO

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine to set up the string of words to be ignored on the input stream.
C The test for IGNOREing is at the end of the G2PRS2 subroutine. By default
C the ignore string is set to '!'. When called the string contains the words
C separated by '!'s.

	SUBROUTINE G2SIGN(WORD)
	CHARACTER*(*)	WORD

	INCLUDE 'gw14.inc'

C Find the length of the input string

	IGLEN = LEN(WORD)
	IGNORE = WORD

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine to get the next word from the input unit NINP and return it in
C WORD with the length in ILEN. If ILEN=0 then the end of the file has been
C reached

	SUBROUTINE G2WORD(LEVEL, WORD, ILEN)
	CHARACTER*(*)	LEVEL, WORD
	INTEGER		ILEN

	INCLUDE 'gw14.inc'

	INTEGER		I, I2, IOLD, INEW, ILIN, NPAR
	LOGICAL		QUOTE

C Clear the word and word length variable

5	WORD = ' '
	ILEN = 0

C Now clear the counter of the degree of brace nesting

	NPAR = 0

C And clear the quotation mark count

	QUOTE = .FALSE.

C Check if a new line needs to be read

	IF (LINE.NE.' ') GOTO 20

C Read in a line from the current file, outputting a prompt if ECHO=0. If the
C source file is not SYS$INPUT or TT then the line is also echoed if ECHO=0.

10	CALL G2READ(NINP, LEVEL, ECHO(NINP-NST+1), LINE, ILIN)
	IF (ILIN.EQ.-1) GOTO 40
	IF (ILIN.EQ. 0) GOTO 10
	IF (LINE.EQ.' ') GOTO 10

C Search through the line for special characters: An exclamation mark when not
C in braces causes the rest of the line to be ignored. Braces (which can be
C deeply nested) surround comments. Braces should not be used (except in pairs)
C after exclamation marks, although exclamation marks may be used freely within
C braces. NPAR is the counter of the current level of brace nesting. A pair of
C braces or (more likely) an exclamation mark may be used to spread single words
C over more than one line in the buffer. Either put an exclamation mark at the
C end of the line or an start in the first character of the next line with the
C rest of the word or put a brace as the last and first character of a line.
C Spaces may be included by enclosing phrases in quotation marks. (However these
C should not be extended across line breaks unless a 'tie' character like '!'
C is used!)

20	DO I = 1, 132
	  IF (LINE(I:I).EQ.'=' .AND. EQUAL) LINE(I:I) = ' '
	  IF (LINE(I:I).EQ.'!' .AND. NPAR.EQ.0) GOTO 10
	  IF (LINE(I:I).EQ.'{') NPAR = NPAR + 1
	  IF (NPAR.EQ.0) THEN
	    IF (LINE(I:I).EQ.'"') THEN
	      QUOTE = .NOT. QUOTE
	    ENDIF
	    IF (LINE(I:I).NE.' ' .AND. LINE(I:I).NE.CHAR(9)) THEN
	      ILEN = ILEN + 1
	      WORD(ILEN:ILEN) = LINE(I:I)
	    ELSE IF (QUOTE) THEN
	      IF (LINE(I:).EQ.' ') THEN
	        WRITE(NOUT, 1020) LEVEL
	        LINE = ' '
	        GOTO 5
	      ELSE
	        ILEN = ILEN + 1
	        WORD(ILEN:ILEN) = ' '
	      ENDIF
	    ELSE IF (ILEN.NE.0) THEN
	      LINE = LINE(I:)
	      GOTO 30
	    ELSE IF (LINE(I:).EQ.' ') THEN
	      GOTO 10
	    ENDIF
	  ELSE
	    IF (LINE(I:).EQ.' ') GOTO 10
	  ENDIF
	  IF (LINE(I:I).EQ.'}') THEN
	    NPAR = NPAR - 1
	    IF (NPAR.LT.0) NPAR = 0
	  ENDIF
	ENDDO

	WRITE(NOUT, 1030) LEVEL
	CALL M2EXIT

C Now remove quotation marks and deal with quote images

30	IOLD = 0
	DO I = ILEN, 1, -1
	  IF (WORD(I:I).EQ.'"') THEN
	    IOLD = IOLD + 1
	  ELSE IF (IOLD.GT.0) THEN
	    INEW = IOLD / 2
	    DO I2 = 1, INEW
	      WORD(I+I2:I+I2) = '"'
	    ENDDO
	    WORD(I+INEW+1:) = WORD(I+IOLD+1:)
	    ILEN = ILEN - IOLD + INEW
	    IOLD = 0
	  ENDIF
	ENDDO

	IF (IOLD.GT.0) THEN
	  INEW = IOLD / 2
	  DO I2 = 1, INEW
	    WORD(I2:I2) = '"'
	  ENDDO
	  WORD(INEW+1:) = WORD(IOLD+1:)
	    ILEN = ILEN - IOLD + INEW
	ENDIF

	RETURN

C At the end of the file return to the next one

40	IF (NINP.GT.NST) THEN
	  CALL G2CLOS(NINP)
	  NINP = NINP - 1
	  LINE = LINEBAK(NINP-NST+1)
	  GOTO 20
	ELSE
	  RETURN
	ENDIF

1020	FORMAT(A,': *** Unbalanced quotation marks ***')
1030	FORMAT(A,': *** Unknown fatal error -- PROGRAM TERMINATED ***')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Variable saving subroutine, called when a line starting $varname=... is typed
C The string is copied with substitution into a buffer so that if a reference to
C another variable is typed then this variable will be substituted when the
C variable is defined, rather than when it is referenced. Substitution can be
C stopped by typing $$varname instead of $varname. Variable names have 16
C character significance, though this number is totally arbitrary. The contents
C of a variable may be examined by typing '$varname=?'.

	SUBROUTINE G2SVAR(WORD, ILEN)
	CHARACTER*(*)	WORD
	INTEGER		ILEN

	INCLUDE 'gw14.inc'

	INTEGER		I, IEQUAL, NOLD, ILEN2, ISTART, IDOLL
	INTEGER		VLEN, VVLEN, IVLEN
	LOGICAL		G2VCHR
	CHARACTER*1	C1
	CHARACTER*80	VARNAM, VARVAL, INVAR, INVAL
	CHARACTER*132	LOLD

	IEQUAL=INDEX(WORD(:ILEN), '=')

C Check to see if the contents of the variable are to be examined, ie. the
C equals sign isn't followed by anything

	IF (ILEN.EQ.IEQUAL) THEN
	  IVLEN = IEQUAL - 1
	  IF (IVLEN.GT.16) IVLEN = 16
	  CALL G2GVAR(WORD(:IVLEN), VARVAL, VVLEN)
	  IF (VVLEN.EQ.0) THEN
	    WRITE(NOUT, 1000) WORD(:IEQUAL-1)
	  ELSE
	    WRITE(NOUT, 1010) WORD(:IEQUAL-1), VARVAL(:VVLEN)
	  ENDIF
	  RETURN
	ENDIF

C As a special case, too, allow variables to consist of several words which may
C be case converted, if the first and last characters of the variable are
C quotation marks then these are removed before continuing

	IF (WORD(IEQUAL+1:IEQUAL+1).EQ.'"' .AND. WORD(ILEN:ILEN).EQ.'"') THEN
	  WORD(ILEN:ILEN) = ' '
	  WORD(IEQUAL+1:ILEN) = WORD(IEQUAL+2:ILEN)
	  ILEN = ILEN - 2
	ENDIF

C Variable input using the '?' notation is performed here, so first check if
C we have the correct form

	IF (WORD(IEQUAL+1:IEQUAL+1).EQ.'?') THEN
	  IF (IEQUAL+1.EQ.ILEN) THEN
	    INVAR = WORD(:IEQUAL-1)
	    IVLEN = IEQUAL - 1
	  ELSE
	    INVAR = WORD(IEQUAL+2:ILEN)
	    IVLEN = ILEN - IEQUAL - 1
	    IF (INVAR(1:1).EQ.'"' .AND. INVAR(IVLEN:IVLEN).EQ.'"') THEN
	      INVAR(IVLEN:IVLEN) = ' '
	      INVAR(:IVLEN) = INVAR(2:IVLEN)
	      IVLEN = IVLEN - 2
	    ENDIF
	  ENDIF
	  IF (NINP.GT.NST) THEN
	    NOLD = NINP
	    LOLD = LINE
	    NINP = NINP - 1
	    LINE = LINEBAK(NINP-NST+1)
	    CALL G2WORD(INVAR(:IVLEN), INVAL, ILEN2)
	    LINEBAK(NINP-NST+1) = LINE
	    DO I = NINP+1, NOLD
	      LINEBAK(I-NST+1) = ' '
	    ENDDO
	    NINP = NOLD
	    LINE = LOLD
	  ELSE
	    CALL G2WORD(INVAR(:IVLEN), INVAL, ILEN2)
	  ENDIF
	  WORD(IEQUAL+1:) = INVAL(:ILEN2)
	  ILEN = IEQUAL + ILEN2
	ENDIF

C Now go through and make the variable substitutions as in the GETPARSE routine

	ISTART = IEQUAL + 1
10	IDOLL = INDEX(WORD(ISTART:ILEN), '$')
	IF (IDOLL.NE.0) THEN
	  IDOLL = ISTART + IDOLL - 1

C Check that next character isn't a '$' as then we want to inhibit variable
C substitution. Instead we remove one of the dollar signs and continue

	  IF (WORD(IDOLL+1:IDOLL+1).EQ.'$') THEN
	    WORD = WORD(:IDOLL)//WORD(IDOLL+2:ILEN)
	    ILEN = ILEN - 1
	    ISTART = IDOLL + 1
	    GOTO 10
	  ENDIF

C Get the variable name its length as this one has to be substituted

	  VARNAM = WORD(IDOLL:)
	  DO I = 2, 80
	    C1 = VARNAM(I:I)
	    IF (.NOT.G2VCHR(C1)) GOTO 12
	  ENDDO
12	  IF (C1.EQ.'$') THEN
	    WORD(IDOLL+I-1:ILEN) = WORD(IDOLL+I:ILEN)
	    ILEN = ILEN - 1
	    VARNAM = WORD(IDOLL:)
	    C1 = VARNAM(I:I)
	  ENDIF
	  VLEN = I - 1

C Go and look for the variable name in the current list

	  CALL G2GVAR(VARNAM(:VLEN), VARVAL, VVLEN)

C If the variable is found then make the substitution

	  IF (VVLEN.NE.0) THEN
	    IF (IDOLL.EQ.1) THEN
	      WORD = VARVAL(:VVLEN)//WORD(VLEN+1:ILEN)
	    ELSE
	      WORD = WORD(:IDOLL-1)//VARVAL(:VVLEN)//WORD(IDOLL+VLEN:ILEN)
	    ENDIF
	    ILEN = ILEN + VVLEN - VLEN
	    ISTART = IDOLL
	  ELSE
	    ISTART = IDOLL + 1
	  ENDIF

C And now go and check if there is another variable to be substituted

	  GOTO 10

	ENDIF

C At this stage the variable substitution has been completed and now we have to
C store it in the appropriate place, that is either to replace an existing
C variable name or to create a new one

	IVLEN = IEQUAL - 1
	IF (IVLEN.GT.16) IVLEN = 16
	CALL G2SV2(WORD(:IVLEN), WORD(IEQUAL+1:ILEN), ILEN-IEQUAL)
	RETURN

1000	FORMAT(' *** ',A,' is not defined as a variable ***')
1010	FORMAT(' *** ',A,' equals "',A,'" ***')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Variable decoding subroutine. The parameters passed are the variable name
C (including the dollar) as a string of the correct length. On return, WORD
C is set to be the value of the variable and ILEN is the length. ILEN is
C returned as 0 if the variable is not found.

	SUBROUTINE G2GVAR(VNAME, WORD, ILEN)
	CHARACTER*(*)	VNAME, WORD
	INTEGER		ILEN

	INTEGER		I, ISTART, IDOLL, VLEN, VVLEN
	CHARACTER*1	C1
	CHARACTER*80	VARNAM,VARVAL
	LOGICAL		G2VCHR

C Try and find the variable's contents

	CALL G2GV2(VNAME, WORD, ILEN)
	IF (ILEN.EQ.0) RETURN

C Now go through and make the variable substitutions as in the G2PRS2 routine

	ISTART = 1
10	IDOLL = INDEX(WORD(ISTART:ILEN), '$')
	IF (IDOLL.NE.0) THEN
	  IDOLL = ISTART + IDOLL - 1

C Get the variable name its length as this one has to be substituted

	  VARNAM = WORD(IDOLL:)
	  DO I = 2, 80
	    C1 = VARNAM(I:I)
	    IF (.NOT.G2VCHR(C1)) GOTO 12
	  ENDDO
12	  IF (C1.EQ.'$') THEN
	    WORD(IDOLL+I-1:ILEN) = WORD(IDOLL+I:ILEN)//' '
	    ILEN = ILEN - 1
	    VARNAM = WORD(IDOLL:)
	    C1 = VARNAM(I:I)
	  ENDIF
	  VLEN = I - 1

C Go and look for the variable name in the current list

	  CALL G2GV2(VARNAM(:VLEN), VARVAL, VVLEN)

C If the variable is found then make the substitution

	  IF (VVLEN.NE.0) THEN
	    IF (IDOLL.EQ.1) THEN
	      WORD = VARVAL(:VVLEN)//WORD(VLEN+1:ILEN)
	    ELSE
	      WORD = WORD(:IDOLL-1)//VARVAL(:VVLEN)//WORD(IDOLL+VLEN:ILEN)
	    ENDIF
	    ILEN = ILEN + VVLEN - VLEN
	    ISTART = IDOLL
	  ELSE
	    ISTART = IDOLL + 1
	  ENDIF

C And now go and check if there is another variable to be substituted

	  GOTO 10

	ENDIF

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine which actually saves the value of a variable, either over the top
C of an existing one or at the end of the list.

	SUBROUTINE G2SV2(VNAME, WORD, ILEN)
	CHARACTER*(*)	VNAME, WORD
	INTEGER		ILEN

	INCLUDE 'gw14.inc'

	INTEGER		I, J, IVLEN

C Don't bother trying to set variables with reserved names

	IF (VNAME.EQ.'$DATE' .OR. VNAME.EQ.'$TIME') THEN
	  WRITE(NOUT, 1000)
	  RETURN
	ENDIF

C See if variable with same name has already been defined. Is so then delete the
C current entry and then proceed

	DO I = 1, NVAR
	  IF (VNAME.EQ.VNAMES(I)) THEN

	    IVLEN = VEND(I) - VST(I) + 1
	    IF (NVAR.EQ.1) THEN
	      VARBUF = ' '
	    ELSE
	      VARBUF(VST(I):NUSED) = VARBUF(VST(I+1):NUSED)//' '
	    ENDIF
	    NVAR = NVAR - 1
	    NUSED = NUSED - IVLEN
	    DO J = I, NVAR
	      VNAMES(J) = VNAMES(J+1)
	      VST(J) = VST(J+1) - IVLEN
	      VEND(J) = VEND(J+1) - IVLEN
	    ENDDO
	    GOTO 10

	  ENDIF
	ENDDO

C If length is zero and it is a new variable name then just ignore the call

10	IF (ILEN.EQ.0) RETURN

C It is a new variable name so create it, checking that there is room for it

	IF (ILEN.GT.2048-NUSED) THEN
	  WRITE(NOUT, 1010)
	  RETURN
	ELSE IF (NVAR.EQ.100) THEN
	  WRITE(NOUT, 1020)
	  RETURN
	ENDIF

	NVAR = NVAR + 1
	VNAMES(NVAR) = VNAME
	VST(NVAR) = NUSED + 1
	VEND(NVAR) = NUSED + ILEN
	VARBUF(NUSED+1:NUSED+ILEN) = WORD(:ILEN)
	NUSED = NUSED + ILEN

	RETURN

1000	FORMAT('SETVAR: *** Variable name reserved ***')
1010	FORMAT('SETVAR: *** Not enough variable space ***')
1020	FORMAT('SETVAR: *** Too many variables defined ***')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine which actually finds the string associated with a variable name, no
C further decoding of the resulting string is performed

	SUBROUTINE G2GV2(VNAME, WORD, ILEN)
	CHARACTER*(*)	VNAME, WORD
	INTEGER		ILEN

	INCLUDE 'gw14.inc'

	INTEGER		I

C Check for the special variables $DATE and $TIME and get these if necessary

	ILEN = 0
	IF (VNAME.EQ.'$DATE') THEN
	  CALL P2DATE(WORD, ILEN)
	ELSE IF (VNAME.EQ.'$TIME') THEN
	  CALL P2TIME(WORD, ILEN)
	ELSE
	  DO I = 1, NVAR
	    IF (VNAME.EQ.VNAMES(I)) THEN
	      WORD = VARBUF(VST(I):VEND(I))
	      ILEN = VEND(I) - VST(I) + 1
	      GOTO 10
	    ENDIF
	  ENDDO
	ENDIF

C If we haven't found the variable then it may be a symbol (VMS) or
C environment variable (UNIX - not implemented) so try and get its value

	IF (ILEN.EQ.0) CALL M2SYMB(VNAME(2:), WORD, ILEN)
10	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Logical function to return TRUE if the character can be in a variable
C name

	LOGICAL FUNCTION G2VCHR(C1)
	CHARACTER*1 C1

	G2VCHR = ((C1.GE.'0'.AND.C1.LE.'9').OR.
     &		  (C1.EQ.'_').OR.
     &		  (C1.GE.'A'.AND.C1.LE.'Z').OR.
     &		  (C1.GE.'a'.AND.C1.LE.'z'))

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutine that will be used by default if a G2HELP call is made and one
C hasn't been supplied in the program. HELP is the number supplied with the
C call to the parser and OPT is the command of a KEYWORD call about which
C help is required

	SUBROUTINE G2HELP(HELP, OPT)
	INTEGER		HELP, OPT

	INCLUDE 'gw14.inc'

C Little bit of dummy code

	WRITE(NOUT,*)
	WRITE(NOUT,*) 'No HELP is available - try ''?'''
	WRITE(NOUT,*) 'Help number:', HELP
	WRITE(NOUT,*) 'Option number:', OPT
	WRITE(NOUT,*)

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Subroutines to deal with UTILITY functions for programs. These functions
C can be accessed from any input without the main program having to be aware
C of the call. This first routine sets up the names of routines which can
C be called

	SUBROUTINE G2SUT(WORD)
	CHARACTER*(*)	WORD

	INCLUDE 'gw14.inc'

	UTLEN = LEN(WORD)
	UTILS = WORD
	UTON = .FALSE.
	INHELP = .FALSE.
	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This subroutine returns the current setting UTILS so that calls can be
C temporarily altered and then restored.

	SUBROUTINE G2GUT(WORD, ILEN)
	CHARACTER*(*)	WORD
	INTEGER		ILEN

	INCLUDE 'gw14.inc'

	ILEN = UTLEN
	IF (LEN(WORD).LT.UTLEN) THEN
	  WORD = UTILS(:LEN(WORD))
	ELSE
	  WORD = UTILS
	ENDIF
	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This is the prototype main UTILITY subroutine, it real programs this should
C be expanded to handle all the possible utility calls. By default the call
C to GETOPEN sets the UTILS string to be 'UTIL!STOP' and it is probably a good
C idea to keep these two options in the string at all times.

	SUBROUTINE G2UTLS(UTILITY)
	CHARACTER*4	UTILITY

	IF (UTILITY.EQ.'UTIL') CALL G2UTIL
	IF (UTILITY.EQ.'STOP') CALL G2STOP

	RETURN

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This is a generally useful subroutine even when the user supplies their own
C set of utility routines as it supplies a list of the currently available
C UTILITIES. If the user desires then this routine may, of course, be replaced
C by an expanded one which gives an explanation of the use of the UTILITIES.

	SUBROUTINE G2UTIL

	INCLUDE 'gw14.inc'

	CHARACTER*75 OUTP

	WRITE(NOUT, *)
	WRITE(NOUT, 1000)
	WRITE(NOUT, *)
	CALL P2DOPT(UTILS, OUTP, '(4X,A)', 8)
	WRITE(NOUT, *)
	RETURN

1000	FORMAT('Currently available utilities are:-')

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Default subroutine to use when STOP is typed, if the program needs to check
C that it is OK to stop then the user may supply their own routine with the same
C name

	SUBROUTINE G2STOP

	CALL P2STMP
	CALL M2EXIT

	END
