## Metabolomic processing using XCMS
- mode: C18neg
- author: Minghao Gong
- origin: the notebook mainly follows Johannes's LCMS process using xcms (https://github.com/sneumann/xcms)
- version: reanalysis from 04/12/2021

## Main modifications from 04/02/2021
- Inspired from extensive discussion here: https://github.com/sneumann/xcms/issues/557
- And a lot of ideas below have been tested in the `parameter_tune4xcms` branch of https://github.com/gmhhope/Subrbier_CHIKV/settings/branches. 
- `ppm` in `cwp` set as 1 instead of 10 as the ppm here is not refering 
  - `XIC plot` will be also good to test on this idea.
- `bw` use 4-5 instead of 1.8 for this experiment. But still needs some test
  - So I will add a chunk for `bw-test`
- In `quantify` using `method = sum` so that split peak intensities can be merged.

## Inclusion list
- Mainly focusing on R5pos versus R5neg
- And incorporation of mcg results can be involved (see if I have time)

## Steps
### XCMS, QC, summarization and mummichog
1. XCMS analysis: code/R/C18neg/X1_Rafi_CD8_xcms_RPneg_wt_pdf.Rmd
- 32726
- Merging reduced 32726 chromPeaks to 25018.
- Processing 1032060 mz slices ... OK (binSize = 0.002)

2. Further QC and summarization: code/R/C18neg/X2_featureSummary_C18neg.R

3. Statistical test (here t-test): code/R/C18neg/X3_creat_ttest4mcg_C18neg.R


### Mummichog analysis and Heatmap visualization on heatmap
4. mcg analysis is performed and may use bash script: /Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/code/bash/mcg_bash

5. After mcg analysis, select the best result and use it for heatmap visualization: code/R/C18neg/pheatmap_Empd_based_on_stat.R

### Inclusion list
1. median and stat filtering: code/R/C18neg/incl_list_X1_stat_med_filtering4incl_list.Rmd

2. merge wth mcg output: code/R/C18neg/Incl_list_X2_proj_merge_filt_table_and_mcg_output.Rmd

3. peak further selection by visualization of chromPeak:code/R/C18neg/chrom_peak_further_select_post_mcg.Rmd