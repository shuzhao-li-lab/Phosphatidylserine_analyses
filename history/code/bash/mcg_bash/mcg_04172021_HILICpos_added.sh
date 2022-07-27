#!/bin/bash

I_PATH="/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/output/v0412_fix/HILICpos/ttest_equal/"

O_PATH="/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/output/v0412_fix/HILICpos/mcg_HILICpos/"

cd $O_PATH

FILE_NAME="ttest_res_padj_R5posvsR5neg.txt"
mummichog -f $I_PATH$FILE_NAME -c 0.05 -o mcg_HILICpos_padj_0.05_R5posvsR5neg
