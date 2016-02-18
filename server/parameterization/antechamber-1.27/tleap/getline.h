#ifndef GETLINE_H
#define GETLINE_H

/* unix systems can #define POSIX to use termios, otherwise 
 * the bsd or sysv interface will be used 
 */
/* GNU libc has a getline with different behavior */
#define getline readline

extern char	*getline(char *prompt);		/* read a line of input */
void            gl_histadd(char *buf);		/* adds entries to hist */

#undef getline

#endif /* GETLINE_H */
