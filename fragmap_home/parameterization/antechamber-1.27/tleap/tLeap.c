/*
 *      File:   leap.c
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Description:
 *              This is the file that starts the LEaP parser.
 *
 *
 *      $Author: case $
 *      $Date: 2005/04/01 22:07:04 $
 *      $Log: tLeap.c,v $
 *      Revision 1.1.1.1  2005/04/01 22:07:04  case
 *      make a stand-alone antechamber
 *
 *      Revision 1.1.1.1  2005/03/17 23:15:27  case
 *      based on dock5.1 code, with post-ACS modifications to nab, etc.
 *
 *      Revision 1.1  2005/02/07 19:55:36  case
 *      copied over from the nab5 distribution
 *
 *      Revision 1.4  2004/04/14 18:10:36  case
 *      update to Amber 8 version of tleap
 *
 *      Revision 7.2  1999/06/15 03:58:06  ross
 *      teLeap compiles
 *
 *      Revision 7.1  1999/04/26 18:45:18  ross
 *      ansify
 *
 *      Revision 7.0  1997/10/08 04:21:38  case
 *      version 5.0
 *
 * Revision 6.1  97/03/29  09:14:30  09:14:30  ross (wilson ross)
 * don't fee line
 * 
 * Revision 6.0  95/03/31  20:53:17  20:53:17  ross (wilson ross)
 * "4.1 release"
 * 
 * Revision 1.6  95/03/23  01:45:46  01:45:46  ross (wilson ross)
 * 1995 copyright
 * 
 * Revision 1.5  95/03/22  04:08:15  04:08:15  ross (wilson ross)
 * romsky
 * 
 * Revision 1.3  93/04/06  00:01:38  schaf
 * Corrected ProgramName stuff
 * 
 * Revision 1.1  92/11/10  20:08:52  schaf
 * Initial revision
 * 
 * Revision 1.1  92/11/10  16:52:17  schaf
 * Initial revision
 * 
 * Revision 1.2  90/11/30  12:52:54  schaf
 * Testing
 * 
 *
 *      The main program for LEaP.
 */



#include        <ctype.h>


#include	"basics.h"

#include        "classes.h"
#include        "dictionary.h"

#include	"parser.h"
#include	"block.h"
#include        "leap.h"

#include        "getline.h"

/* GNU libc has a getline with different behavior */
#define getline readline

/*
 *******************************************************************
 *
 *	MAIN PROGRAM
 *
 *
 *      Amended: Vladimir Romanovski (1994)
 *        It was added the ``input-edit'' library facility.
 */

static char *
stripwhite (char *string)
{
  register char *s, *t;

  for (s = string; isspace (*s); s++)
    /*EMPTY*/;
    
  if (*s == '\0')
    return s;

  t = s + strlen (s) - 1;
  while (t > s && isspace (*t))
    t--;
  
  *++t = '\0';

  return s;
}


main( int argc, char *argv[] )
{
BOOL		bUseStartup, bGotCmd;
BLOCK		bCmd;
char		c;
RESULTt		rResult;


    setbuf(stdout,NULL);
    
    strcpy( GsProgramName, argv[0] );

    BasicsInitialize();

    GbGraphicalEnvironment = FALSE;

    ParseArguments( argc, argv );

    ParseInit( &rResult );

    bGotCmd = FALSE;
    bCmd = bBlockCreate();

    while ( rResult.iCommand != CQUIT ) {
      char *line=NULL;
      char *cp;

      line = getline("> ");
      cp = stripwhite (line);


      if (*cp){
	gl_histadd(cp);
	
	for(; *cp; cp++) 
	  (void) bBlockAddChar( bCmd, *cp );
	
	if ( bBlockAddChar( bCmd, '\n' ) ) {
	  ParseBlock( bCmd, &rResult );
	  BlockEmpty( bCmd );
	}
      }
      
    } while ( rResult.iCommand != CQUIT );

    ParseShutdown();

    LISTUNFREEDMEMORYTOLOGFILE();

    exit(0);
}

