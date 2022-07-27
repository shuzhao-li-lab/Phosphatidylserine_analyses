#!/bin/bash
I_PATH=/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/input/batch1/raw
O_PATH=/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/input/batch1/mzML
mono /Users/gongm/Documents/packages/ThermoRawFileParser/ThermoRawFileParser.exe -d $I_PATH -o $O_PATH