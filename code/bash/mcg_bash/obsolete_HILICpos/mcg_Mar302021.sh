#!/bin/bash

I_PATH="/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/output/HILICpos_031821run/cleanup_and_stat_test/"

O_PATH="/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/output/HILICpos_031821run/mcg/mcg_032021/"


mkdir $O_PATH

cd $O_PATH

FILE_NAME="ttest_padj4mcg/ttest_res_padj_R5posvsR5neg.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_HILICpos_padj_default_R5posvsR5neg
mummichog -f $I_PATH$FILE_NAME -c 0.25 -o mcg_HILICpos_padj_025_R5posvsR5neg

FILE_NAME="ttest_padj4mcg/ttest_res_padj_R5posvsNaive.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_HILICpos_padj_default_R5posvsNaive
mummichog -f $I_PATH$FILE_NAME -c 0.25 -o mcg_HILICpos_padj_025_R5posvsNaive

FILE_NAME="ttest_padj4mcg/ttest_res_padj_R5negvsNaive.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_HILICpos_padj_default_R5negvsNaive
mummichog -f $I_PATH$FILE_NAME -c 0.25 -o mcg_HILICpos_padj_025_R5negvsNaive





FILE_NAME="ttest_RawPval4mcg/ttest_res_rawpval_R5posvsR5neg.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_HILICpos_RawPval_default_R5posvsR5neg
mummichog -f $I_PATH$FILE_NAME -c 0.05 -o mcg_HILICpos_RawPval_005_R5posvsR5neg


FILE_NAME="ttest_RawPval4mcg/ttest_res_rawpval_R5negvsNaive.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_HILICpos_RawPval_default_R5negvsNaive
mummichog -f $I_PATH$FILE_NAME -c 0.05 -o mcg_HILICpos_RawPval_005_R5negvsNaive


FILE_NAME="ttest_RawPval4mcg/ttest_res_rawpval_R5posvsNaive.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_HILICpos_RawPval_default_R5posvsNaive
mummichog -f $I_PATH$FILE_NAME -c 0.05 -o mcg_HILICpos_RawPval_005_R5posvsNaive

