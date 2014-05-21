#!/bin/bash

# The name of the program
PROGRAM="sq_lgr"

# Create the FILEBASE, this is the folder where the testing will be done.
THISFILE=$0 # This contains "./", which we do not want.
THISFILE=${THISFILE:2} # Gets rid of "./"
FILEBASE=${THISFILE%%.*} # Get rid of the extension (in this case, ".sh")

# Create the new folder (remove old one)
touch $FILEBASE
rm -rf $FILEBASE
mkdir $FILEBASE

# Get the current directory and the oomph-base
CURRENT_DIR=`pwd`
OOMPHROOT_DIR=$(make -s --no-print-directory print-top_builddir)

# folder of where the iteration counts will be.
RESITS_DIR="res_iterations"

# Get version of oomph-lib
cd $CURRENT_DIR
cd $OOMPHROOT_DIR
git log -1 > $CURRENT_DIR/$FILEBASE/oomphlib_revision
cd $CURRENT_DIR

# Get version of user drivers
cd $OOMPHROOT_DIR/user_drivers
git log -1 > $CURRENT_DIR/$FILEBASE/user_driver_revision
cd $CURRENT_DIR

# make the program and move it into the test folder.
make $PROGRAM
mv $PROGRAM ./$FILEBASE
cd $FILEBASE

###############################################################################
# Now we are inside the test folder (FILEBASE)

# Make the results directory. I have it in an if statement...
# Because some times when it exists, we wish to reuse it... of course not this 
# time...
if [ ! -d "$RESITS_DIR" ]; then
  mkdir $RESITS_DIR
fi

# There may be two test lists, so we declare the strings up here.
# Then we change this right before we call the gen_testxy functions.
TEST_FILEBASE=""
TEST_LIST=""

# NOTE: The F amg settings are the same as --f_solver 69

## Preconditioner parameters.
Famg_BASE="--f_solver 96"
Famg_ITER="--f_amg_iter 1"
Famg_SMITER="--f_amg_smiter 2"
Famg_SSMOOTHER="--f_amg_sim_smoo 1" # GS
Famg_CSMOOTHER=""
Famg_DAMP="--f_amg_damp -1"
Famg_STRN_SIM="--f_amg_str 0.25"
Famg_STRN_STR="--f_amg_str 0.668"
Famg_COARSE="--f_amg_coarse 1" #RS

Prec_WLu_NSLu="--w_solver 0 --ns_solver 0"
Prec_WLu_NSLscLu="--w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0"
Prec_WLu_NSLscPamgFlu="--w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 0"

Famg_sim="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SSMOOTHER $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
Famg_str="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SSMOOTHER $Famg_DAMP $Famg_STRN_STR $Famg_COARSE"
Prec_WLu_NSLscPLuFamgsim="--w_solver 0 --ns_solver 1 --p_solver 0 $Famg_sim"
Prec_WLu_NSLscPLuFamgstr="--w_solver 0 --ns_solver 1 --p_solver 0 $Famg_str"
Prec_WLu_NSLscPamgFamgsim="--w_solver 0 --ns_solver 1 --p_solver 1 $Famg_sim"
Prec_WLu_NSLscPamgFamgstr="--w_solver 0 --ns_solver 1 --p_solver 1 $Famg_str"


## Creates test lists for all prec combinations for noel = 4 to 128
function gen_tests256()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - W SuperLU, NS SuperLU
# 1 - W SuperLU, NS LSC: P SuperLU, F SuperLU
# 2 - W SuperLU, NS LSC: P AMG, F LU
# 3 - W SuperLU, NS LSC: P Lu, F AMG
# 4 - W Super LU, NS LSC: P AMG F AMG

PRECLIST="4"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0 1"
ANGLIST="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360"
RELIST="200"
NOELLIST="256"

for PREC  in $PRECLIST
do

  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do
case "$PREC" in
  0)
    PRECPARAM="$Prec_WLu_NSLu"
    ;;
  1)
    PRECPARAM="$Prec_WLu_NSLscLu"
    ;;
  2)
    PRECPARAM="$Prec_WLu_NSLscPamgFlu"
    ;;
  3)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$Prec_WLu_NSLscPLuFamgsim"
    else
      PRECPARAM="$Prec_WLu_NSLscPLuFamgstr"
    fi
    ;;
  4)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$Prec_WLu_NSLscPamgFamgsim"
    else
      PRECPARAM="$Prec_WLu_NSLscPamgFamgstr"
    fi
    ;;
esac

echo "mpirun -np 1 ./$PROGRAM --dist_prob --trilinos_solver --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST

        done
      done
    done
  done
done
} # gen_tests function


TEST_FILEBASE="tests_256_varyAng"
TEST_LIST="$TEST_FILEBASE.list"
gen_tests256
TEST_RUN="$TEST_FILEBASE.sh"
echo "#!/bin/bash" >> $TEST_RUN
cat $TEST_LIST >> $TEST_RUN

cp ./../$0 .

