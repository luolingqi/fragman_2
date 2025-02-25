/*
 *      File:   unitio.h
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
 *              Part of the UNIT object.
 *              All file input/output routines have been
 *              placed in the file 'unitio.c'
 *
 */
 
#ifndef UNITIO_H
#define UNITIO_H
 

extern BOOL     zbUnitIOLoadTables(UNIT uUnit, DATABASE db);
extern void     zUnitIOSaveTables(UNIT uUnit, DATABASE db);

extern void     zUnitIOBuildTables(UNIT uUnit, PARMLIB plParameters,
                        BOOL *bPGeneratedParameters, BOOL bPert, BOOL bCheck);
extern void     zUnitIOBuildFromTables(UNIT uUnit);
extern void     zUnitIODestroyTables(UNIT uUnit);

extern void     zUnitIOSaveAmberParmFormat(UNIT uUnit, FILE *fOut,
                        FILE *fCrd, BOOL bPolar, BOOL bPert);
extern void     zUnitIOSaveAmberParmFormat_old(UNIT uUnit, FILE *fOut,
                        FILE *fCrd, BOOL bPolar, BOOL bPert);

extern void     UnitIOSaveAmberPrep( UNIT uUnit, FILE *fOut );

#endif  /* UNITIO_H */
