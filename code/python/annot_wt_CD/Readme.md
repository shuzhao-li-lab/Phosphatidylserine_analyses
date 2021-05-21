# Handle CD-exported Cpd Table in Rafi-DeepScan project
- The annotations can be from DeepScan (assuming there are no DDA runs).
- This part of script handle Compound Table exported from Compound Discoverer.

## Issues to be solved
- [] As the Compound table only contains m/z pertained to probably the most abundant ions (we don't know if the degenerate ions, or features associated to a particular empirical compound can be showed all together in the Compound Table).

## How it works
- Basically match m/z with certain ppm in both the compound table and feature table. So then the annotated table can be exported as 
  - full table 
  - only annotated (including chemical composition predictions...)

## Resoures
- `/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/code/python/annotation`