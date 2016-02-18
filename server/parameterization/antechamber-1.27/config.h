#  This file should work for a standard gcc/g77 installation.
#  You will need to edit it if you do not have the GNU compilers and tools
#     installed.

SFX=
#  uncomment the following line for cygwin:
#SFX=.exe

CC=      gcc
CXX=     g++
CFLAGS=  -O2 
OCFLAGS= -O3
FC=      gfortran
FFLAGS	= -O2 -fno-automatic -finit-local-zero
LEX=     flex
RANLIB=  ranlib
AR=      ar rv 
LM=      -lm
