This directory contains data for sample pose estimation runs with posest

Monocular:
==========

Typing

cd mono; posest_demo.exe K.txt 32D.txt 

should produce output similar to

--------------------------------------------------------------------------------------------------

Camera pose estimation using 891 image matches
PnP pose: -0.0538266 0.950572 0.781393 -1.53928 -0.821313 0.843891

Refinement using 1430 measurements, 6 variables
LM returned 6 in 6 iter, reason 2, error 0.667234 [initial 0.862387], 12/6 func/fjac evals

Estimated motion ([rv t]) [176 outliers, 19.75%]
-0.05235749 0.9534564 0.7841375 -1.541927 -0.8226765 0.8480734 

Elapsed time: 0.01 seconds, 5.00 msecs

Reprojection RMS and RMedS errors for input points: 110.167 0.563253
	25%, 50% and 75% quartiles: 0.308561 0.563253 1.58899

--------------------------------------------------------------------------------------------------

Typing
posest_demo.exe K_f0.txt 32D.txt

will also estimate the focal length and should produce output similar to

--------------------------------------------------------------------------------------------------

Camera pose estimation using 891 image matches

Refinement using 1432 measurements, 7 variables
LM returned 28 in 28 iter, reason 2, error 0.785415 [initial 0.968801], 33/28 func/fjac evals

Estimated motion & focal length [175 outliers, 19.64%]
-0.05403624 0.9542258 0.7833138 -1.541982 -0.8253089 0.9091973 3078.976 

Elapsed time: 0.06 seconds, 61.00 msecs

Reprojection RMS and RMedS errors for input points: 110.228 0.66289
	25%, 50% and 75% quartiles: 0.387375 0.66289 1.65318

--------------------------------------------------------------------------------------------------


Binocular:
==========

Typing
cd bino; binocposest_demo Pleft.txt matches32Dleft.txt Pright.txt matches32Dright.txt 

should produce output similar to

--------------------------------------------------------------------------------------------------

PnP pose: -0.157923 2.38138 1.77346 -175.537 -603.184 1186.9
Init refinement to PnP pose (in-place): 0.0001 0 0 17.7496 -34.9395 627.689
MLSL return reason 5, evals 92243, samples 250, minf 199.063
MLSL pose: 0.0045704 -0.00567251 -0.0043006 19.9922 -34.0209 625.749, error 1.6184

Refinement using 246 measurements, 6 variables
LM returned 1 in 1 iter, reason 2, error 1.6184 [initial 1.6184], 1/1 func/fjac evals
posestBinoc(): estimated motion from left frame ([rv t]) [106 outliers, 46.29%]
	-0.156919 2.37181 1.77406 -178.903 -603.97 1181.25
PnP pose: -0.161443 2.35616 1.78562 -132.784 -608.457 1177.59
Init refinement to PnP pose (in-place): 0.0001 0 0 17.3507 -32.6357 630.41
MLSL return reason 5, evals 77250, samples 250, minf 153.212
MLSL pose: -0.00996569 -0.0020437 -0.00089337 18.2696 -36.2153 628.746, error 1.38029

Refinement using 222 measurements, 6 variables
LM returned 1 in 1 iter, reason 2, error 1.38029 [initial 1.38029], 1/1 func/fjac evals
posestBinoc(): estimated motion from right frame ([rv t]) [110 outliers, 49.77%]
	-0.163473 2.3637 1.77273 -133.5 -606.421 1181.32

Binocular refinement #1 (L-R) using 468 measurements, 6 variables
posestBinoc(): estimated motion from left & right frames [216 outliers in 450 combined matches, 48.00%]
LM returned 20 in 20 iter, reason 2, error 2.99204 [initial 7.56598], 26/20 func/fjac evals
Refined L-R motion ([rv t])
	-0.157404 2.38373 1.76061 -180.299 -601.512 1198.02

Binocular refinement #2 (R-L) using 468 measurements, 6 variables
posestBinoc(): estimated motion from left & right frames [216 outliers in 450 combined matches, 48.00%]
LM returned 21 in 21 iter, reason 2, error 2.99204 [initial 10.5223], 25/21 func/fjac evals
Refined R-L motion ([rv t])
	-0.164433 2.37358 1.75406 -136.687 -603.568 1201.63

Binocularly refined motion ([rv t])
	-0.157404 2.38373 1.76061 -180.299 -601.512 1198.02

Refined motion ([rv t  scale]) [216 outliers, 48.00%]
	-0.157404 2.38373 1.76061 -180.299 -601.512 1198.02  1.000000


Elapsed time: 0.82 seconds, 825.00 msecs
Reprojection errors 25%, 50% and 75% quartiles: 1.36653 3.74435 116.657

--------------------------------------------------------------------------------------------------
