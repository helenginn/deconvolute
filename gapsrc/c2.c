@IF !VMS
/* Machine dependent routines for the P2 parser. Following the M2.FOR
   routines, all functions begin with 'C2'. Since this language
   conversion stuff is a real pain it will be kept to a minimum */

@IF CONVEX
/* This is the Convex C-series version */
@ELSEIF HP
/* This is the HP PA-RISC version */
@ELSEIF SGI
/* This is the SGI MIPS-RISC version */
@ELSEIF LINUX64
/* This is the 64-bit Linux version */
@ELSEIF LINUX32
/* This is the 32-bit Linux version */
@ELSE
/* This is the bland, unsupported Unix version */
@ENDIF

/* Robert Esnouf, 18/7/93 - 12/11/93 */

@IF CONVEX
#include <stdio.h>
@ELSEIF HP
#include <stdio.h>
#include <sys/utsname.h>
#include <pwd.h>
@ELSEIF SGI
#include <stdio.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <pwd.h>
#include <malloc.h>
@ELSEIF LINUX64 | LINUX32
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <malloc.h>
#include <time.h>
@ENDIF

/* Make the standard output unit line buffered so that sniff will
   work nicely
*/

@IF CONVEX
void c2lbuf_()
{ setlinebuf(stdout);
}
@ELSEIF HP
void c2lbuf()
{ setvbuf(stdout, NULL, _IOLBF, 0);
}
@ELSEIF SGI
void c2lbuf_()
{ setlinebuf(stdout);
}
@ELSEIF LINUX64 | LINUX32
void c2lbuf_()
{ setvbuf(stdout, (char *) NULL, _IOLBF, 0);
}
@ENDIF

/* Get the operating system version and the machine type into a string
   for FORTRAN
*/

@IF CONVEX
void c2mdat(osver, ilen1, mctyp, ilen2)
char *osver;
int *ilen1;
char *mctyp;
int *ilen2;
{ struct utsname name;
  char cosver[2*UTSLEN+1];

  if (uname(&name) == -1)
    *ilen1 = 0;
  else
  { strcpy(cosver, name.sysname);
    strcat(cosver, " ");
    strcat(cosver, name.release);
    strcpy(osver, cosver);
    *ilen1 = strlen(name.sysname) + strlen(name.release) + 1;
    strcpy(mctyp, "HP ");
    strcat(mctyp, name.machine);
    *ilen2 = strlen(name.machine) + 3;
  }
  
}
@ELSEIF SGI
void c2mdat_(osver, ilen1, mcnam, ilen2, mctyp, ilen3)
char *osver;
int *ilen1;
char *mcnam;
int *ilen2;
char *mctyp;
int *ilen3;
{ struct utsname name;
  char cosver[2*UTSLEN+1];

  if (uname(&name) == -1)
    *ilen1 = 0;
  else
  { strcpy(cosver, name.sysname);
    strcat(cosver, " ");
    strcat(cosver, name.release);
    strcpy(osver, cosver);
    *ilen1 = strlen(name.sysname) + strlen(name.release) + 1;
    strcpy(mcnam, name.nodename);
    *ilen2 = strlen(name.nodename);
    strcpy(mctyp, "SGI ");
    strcat(mctyp, name.machine);
    *ilen3 = strlen(name.machine) + 4;
  }
  
}
@ELSEIF LINUX64 | LINUX32
#define UTSLEN 9

void c2mdat_(ostyp, ilen1, osver, ilen2, mcnam, ilen3, mctyp, ilen4)
char *ostyp;
int *ilen1;
char *osver;
int *ilen2;
char *mcnam;
int *ilen3;
char *mctyp;
int *ilen4;
{ struct utsname name;

  if (uname(&name) == -1)
    *ilen1 = 0;
  else
  { strcpy(ostyp, name.sysname);
    *ilen1 = strlen(name.sysname);
    strcpy(osver, name.release);
    *ilen2 = strlen(name.release);
    strcpy(mcnam, name.nodename);
    *ilen3 = strlen(name.nodename);
    strcat(mctyp, name.machine);
    *ilen4 = strlen(name.machine);
  }
  
}
@ENDIF

/* Get the username associated with the current process even if the job
   is running detatched from a terminal
*/

@IF CONVEX
void c2user(usernm, ilen)
char *usernm;
int *ilen;
{ struct passwd *pws;

  if ((pws = getpwuid(getuid())) == NULL)
    *ilen = 0;
  else
  { strcpy(usernm, pws->pw_name);
    *ilen = strlen(pws->pw_name);
  }
}
@ELSEIF HP
void c2user(usernm, ilen)
char *usernm;
int *ilen;
{ struct passwd *pws;

  if ((pws = getpwuid(getuid())) == NULL)
    *ilen = 0;
  else
  { strcpy(usernm, pws->pw_name);
    *ilen = strlen(pws->pw_name);
  }
}
@ELSEIF SGI
void c2user_(usernm, ilen)
char *usernm;
int *ilen;
{ struct passwd *pws;

  if ((pws = getpwuid(getuid())) == NULL)
    *ilen = 0;
  else
  { strcpy(usernm, pws->pw_name);
    *ilen = strlen(pws->pw_name);
  }
}
@ELSEIF LINUX64 | LINUX32
void c2user_(usernm, ilen)
char *usernm;
int *ilen;
{ struct passwd *pws;

  if ((pws = getpwuid(getuid())) == NULL)
    *ilen = 0;
  else
  { strcpy(usernm, pws->pw_name);
    *ilen = strlen(pws->pw_name);
  }
}
@ENDIF

/* Allocate the specified number of bytes and return the starting address
   of the allocation
*/

@IF CONVEX
void c2mema_(nbytes, start)
long long *nbytes;
long long *start;
{ *start = (unsigned long long) malloc((size_t)*nbytes);
}
@ELSEIF HP
void c2mema_(nbytes, start)
long long *nbytes;
long long *start;
{ *start = (unsigned long long) malloc((size_t)*nbytes);
}
@ELSEIF SGI
void c2mema_(nbytes, start)
long long *nbytes;
long long *start;
{ *start = (unsigned long long) malloc((size_t)*nbytes);
}
@ELSEIF LINUX64 | LINUX32
void c2mema_(nbytes, start)
long long *nbytes;
uintptr_t *start;
{ *start = (uintptr_t)malloc((size_t)*nbytes);
}
@ENDIF

/* Deallocate the memory allocated with the c2mema function
*/

@IF CONVEX
void c2memd_(start)
long long *start;
{ free((void *)*start);
}
@ELSEIF HP
void c2memd_(start)
long long *start;
{ free((void *)*start);
}
@ELSEIF SGI
void c2memd_(start)
long long *start;
{ free((void *)*start);
}
@ELSEIF LINUX64 | LINUX32
void c2memd_(start)
uintptr_t *start;
{ free((void *)*start);
}
@ENDIF

/* Return the 64-bit address of a variable because %LOC is only a 32-bit call
*/

@IF CONVEX
void c2loc_(addr, item)
long long *addr;
int item;
{ *addr = (unsigned long long) item;
}
@ELSEIF HP
void c2loc_(addr, item)
long long *addr;
int item;
{ *addr = (unsigned long long) item;
}
@ELSEIF SGI
void c2loc_(addr, item)
long long *addr;
int item;
{ *addr = (unsigned long long) item;
}
@ENDIF

/* Return the 64-bit address of a character variable because %LOC is only 
   a 32-bit call
*/

@IF CONVEX
void c2cloc_(addr, item)
long long *addr;
int item;
{ *addr = (unsigned long long) item;
}
@ELSEIF HP
void c2cloc_(addr, item)
long long *addr;
int item;
{ *addr = (unsigned long long) item;
}
@ELSEIF SGI
void c2cloc_(addr, item)
long long *addr;
int item;
{ *addr = (unsigned long long) item;
}
@ELSEIF LINUX64 | LINUX32
void c2cloc_(addr, item)
long long *addr;
int item;
{ *addr = (unsigned long long) item;
}
@ENDIF

@IF LINUX32 | LINUX64

/* Extra system functions moved from Fortran to C in Linux */

/* Execute a system command using the standard C call
*/

int c2system_(comand)
char *comand;
{ return (system(comand));
}

/* Return the current time as an 8 character string of the form 12:34:56
*/

void c2time_(hour, min, sec)
int *hour;
int *min;
int *sec;

{ time_t tnow;
  struct tm *dnow;

  tnow = time((time_t *)NULL); 
  dnow = localtime(&tnow); 

  *hour = dnow->tm_hour;
  *min = dnow->tm_min;
  *sec = dnow->tm_sec;
}

/* Return the current date as three integers DAY, MONTH, TIME (4 digits)
*/

void c2idate_(day, month, year)
int *day;
int *month;
int *year;
{ time_t tnow;
  struct tm *dnow;

  tnow = time((time_t *)NULL); 
  dnow = localtime(&tnow); 
  *day = dnow->tm_mday;
  *month = dnow->tm_mon + 1;
  *year = dnow->tm_year + 1900;
}

@ENDIF
@ENDIF
