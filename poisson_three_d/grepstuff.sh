#!/bin/bash

grep "Processor $1:" $2 | grep "Time to generate Jacobian\|Setting up BoomerAMG\|Time for preconditioner setup\|Time to generate Trilinos matrix\|Linear solver iterations\|Time for trilinos solve itself\|Time for complete trilinos solve\|Time for linear solver\|Total time for linear solver\|Total time for Newton solver\|Time outside linear solver"


