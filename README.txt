DESCRIPTION
Scripts to reproduce figures and results from Young et al. (2022) PNAS

There are two main scripts to be run that process ApRES data to obtain basal melt rates:
- A_trackBed.m
- B_basalMelt.m

The script 'A_trackBed.m' produces a time series of bed depths (relative to the vertical location of the ApRES antenna array) and tracks its change through time. This is done through observing the convergence of bursts within each daily bin. 

The script 'B_basalMelt.m' solves for basal melt rates by combining contemporaneous observations of vertical strain rate and bed depth.

Note 1: Script 'B_basalmelt.m' is currently incomplete--the newest version of this script is on an old computer halfway around the world (sorry!). 
Note 2: Scripts are dependent on the fmcw package from Stewart (2018) and the output data from Young et al. (2019). 

REFERENCES
Stewart, C.L. (2018) Ice-ocean interactions beneath the north-western Ross Ice Shelf, Antarctica (Doctoral thesis). Scott Polar Research Institute, University of Cambridge. 
Young, T.J. et al. (2019) Physical Conditions of Fast Glacier Flow: 3. Seasonally-Evolving Ice Deformation on Store Glacier, West Greenland. Journal of Geophysical Research: Earth Surface 124(1):245-267. 
Young, T.J. et al. (2022) Rapid basal melting of the Greenland Ice Sheet from surface meltwater drainage. Proceedings of the National Academy of Sciences 119 (e2116036119). 

METADATA
Author: T. J. Young, Poul Christoffersen
Institution: Scott Polar Research Institute, University of Cambridge
Associated Project: RESPONDER
Date compiled: 10 February 2022
