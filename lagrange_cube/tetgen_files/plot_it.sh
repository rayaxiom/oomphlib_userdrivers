#!/bin/bash


FILE="CuPoQStrR100WedNlFePeAx0y0z0N4t0"

oomph-convert -z $FILE*.dat
makePvd $FILE $FILE.pvd


