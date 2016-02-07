C***********************************************************************
C The machine dependent routines for the new parser, P2, by RME. It
C replaces the ailing GETWORD. These routines should be kept to a
C minimum and will have to be changed (or checked) for each implement-
C ation of P2. Functions performed by these routines include file
C manipulation, system functions like memory allocation and date, CPU
C time etc as well as interfaces to system services written in C.
C
@IF VMS
C This is the VAX/Alpha VMS version
@ELSEIF CONVEX
C This is the Convex C-series version
@ELSEIF HP
C This is the HP PA-RISC version
@ELSEIF SGI
C This is the SGI MIPS-RISC version
@ELSE
C This is a bland version for an unspecified machine
@ENDIF
C
C Robert Esnouf, 7/7/93 - 24/7/93
C
C***********************************************************************

