#!/bin/bash


FILE="CuPoQStrR100WedNlFePeAx0y0z0N4t"

oomph-convert -z $FILE*.dat
makePvd $FILE $FILE.pvd


