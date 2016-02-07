/* Machine dependent routines for the P2 parser. Following the M2.FOR
   routines, all functions begin with 'C2'. Since this language
   conversion stuff is a real pain it will be kept to a minimum */

/* This is the 32-bit Linux version */

/* Robert Esnouf, 18/7/93 - 12/11/93 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
//#include <malloc.h>
#include <time.h>

/* Make the standard output unit line buffered so that sniff will
   work nicely
*/

void c2lbuf_()
{ setvbuf(stdout, (char *) NULL, _IOLBF, 0);
}

/* Get the operating system version and the machine type into a string
   for FORTRAN
*/

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

/* Get the username associated with the current process even if the job
   is running detatched from a terminal
*/

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

/* Allocate the specified number of bytes and return the starting address
   of the allocation
*/

void c2mema_(nbytes, start)
long long *nbytes;
uintptr_t *start;
{ *start = (uintptr_t)malloc((size_t)*nbytes);
}

/* Deallocate the memory allocated with the c2mema function
*/

void c2memd_(start)
uintptr_t *start;
{ free((void *)*start);
}

/* Return the 64-bit address of a variable because %LOC is only a 32-bit call
*/


/* Return the 64-bit address of a character variable because %LOC is only
   a 32-bit call
*/

void c2cloc_(addr, item)
long long *addr;
int item;
{ *addr = (unsigned long long) item;
}


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

