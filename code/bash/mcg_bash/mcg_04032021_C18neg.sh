#!/bin/bash

I_PATH="/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/output/RPneg_0402/ttest_equal/"

O_PATH="/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/output/RPneg_0402/mcg/"
mkdir $O_PATH

cd $O_PATH

FILE_NAME="ttest_res_padj_R5posvsR5neg.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_C18neg_padj_default_R5posvsR5neg -m negative
mummichog -f $I_PATH$FILE_NAME -c 0.25 -o mcg_C18neg_padj_025_R5posvsR5neg -m negative

FILE_NAME="ttest_res_padj_R5posvsNaive.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_C18neg_padj_default_R5posvsNaive -m negative
mummichog -f $I_PATH$FILE_NAME -c 0.25 -o mcg_C18neg_padj_025_R5posvsNaive -m negative

FILE_NAME="ttest_res_padj_R5negvsNaive.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_C18neg_padj_default_R5negvsNaive -m negative
mummichog -f $I_PATH$FILE_NAME -c 0.25 -o mcg_C18neg_padj_025_R5negvsNaive -m negative





FILE_NAME="ttest_res_rawpval_R5posvsR5neg.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_C18neg_RawPval_default_R5posvsR5neg -m negative
mummichog -f $I_PATH$FILE_NAME -c 0.05 -o mcg_C18neg_RawPval_005_R5posvsR5neg -m negative


FILE_NAME="ttest_res_rawpval_R5negvsNaive.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_C18neg_RawPval_default_R5negvsNaive -m negative
mummichog -f $I_PATH$FILE_NAME -c 0.05 -o mcg_C18neg_RawPval_005_R5negvsNaive -m negative


FILE_NAME="ttest_res_rawpval_R5posvsNaive.txt"
mummichog -f $I_PATH$FILE_NAME -o mcg_C18neg_RawPval_default_R5posvsNaive -m negative
mummichog -f $I_PATH$FILE_NAME -c 0.05 -o mcg_C18neg_RawPval_005_R5posvsNaive -m negative

