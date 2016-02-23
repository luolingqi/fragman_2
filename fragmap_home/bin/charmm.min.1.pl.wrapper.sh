#!/bin/bash

#this program is purely to set ulimit to unlimited before running charmm
ulimit -s unlimited
charmm.min.1.pl $@
