## Metabolomic processing using XCMS
- mode: HILICpos
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
1. XCMS analysis: code/R/HILICpos/Rafi_CD8_xcms_HILICpos_wt_pdf_v0402.Rmd
2. Further QC and summarization: code/R/HILICpos/featureSummary_HILICpos.R

3. Statistical test (here t-test): code/R/HILICpos/creat_ttest4mcg_HILICpos_equa_var_t_test.R


### Mummichog analysis and Heatmap visualization on heatmap
4. mcg analysis is performed and may use bash script: /Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/code/bash/mcg_bash

5. After mcg analysis, select the best result and use it for heatmap visualization: code/R/HILICpos/pheatmap_Empd_based_on_stat_HILICpos.R

### Inclusion list
1. median and stat filtering: code/R/HILICpos/incl_list_X1_stat_med_filtering4incl_list.Rmd

2. [ ] merge wth mcg output (may do it or not): 

3. [ ] peak further selection by visualization of chromPeak: 