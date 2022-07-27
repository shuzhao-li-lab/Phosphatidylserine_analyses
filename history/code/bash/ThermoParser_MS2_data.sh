#!/bin/bash
I_PATH=/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/input/04292021_Rafi_RP_neg_AX_DeepScan/04292021_Rafi_RP_neg_AX_DeepScan
O_PATH=/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/input/04292021_Rafi_RP_neg_AX_DeepScan/mzML
mono /Users/gongm/Documents/packages/ThermoRawFileParser/ThermoRawFileParser.exe -d $I_PATH -o $O_PATH