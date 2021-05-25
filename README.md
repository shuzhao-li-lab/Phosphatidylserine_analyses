Rafi_Ahmed_exhaustedCD8Tcell analysis post 04022021
==============================

This branch of analysis document the analysis after finding of mislabeling of HILICpos & RPneg data. So I decide to redo analysis from the very beginning to avoid any twists across the wrong paths.

The data repository: https://www.synapse.org/#!Synapse:syn12613596/wiki/581963

## This branch  - MS2_Spectra_massbank_match_05152021
- This branch is created to match MS2 spetra using deepScan Acquire X mode with the significant list. 
- Because of the significant list previously created is based on a stringent criteria (padj < 0.05 & FC > 2).
- To maximize the search I will recreate a list with less stringent criteria (R5pos vs R5neg; padj < 0.05)
- Then I will use the Spectra massbank matching script developed in CHIKV project to do a search.
  - I have also done the first round of Compound Discoverer compound identification, though very few was matched
  - The results of HILICpos part is stored and documented in https://www.synapse.org/#!Synapse:syn12613596/wiki/610530

## The branch - annotation-with-CD-exported-tables
- This branch is to merge feature table with CD exported Compound tables.
- So that we have another annotation from Compound Discoverer (esp. mzCloud) together with public available annotations （e.g., HMDB & massbank).

## The branch - Prepare4Shuzhao-li_Seminar_052321
- A rushing Sunday preparing multiple things for Shuzhao Li's seminar

### The code
- python script to merge annotation of Compound Discoverer with xcms-feature-table: `/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/code/python/annot_wt_CD`
- mzCloud, massbank, hmdb results sum-up: `/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/code/R/HILICpos/mzCloud_mbank_hmdb_res_sumup` 
  - Importantly, this part has assumed a manual step to reorganize the spectra-exported massbank/hmdb table and the CD-exported mzCloud table. A lot needs to improve here.

### Reports and figures
- The slides are in `/Users/gongm/Documents/projects/Presentation_seminar_poster/Shuzhao_Seminar/slides`
- There are also additional figures/tables created in 
  - Compound-Disocover exported table merged with the xcms-feature-table (for convenience, here only using significant feature table): `/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/output/annotation-with-CD-exported-tables`
    - some Compound-Discoverer mirror plots in `/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/output/annotation-with-CD-exported-tables/Rafi-CD8-from-Compound-Discoverer/mirror_plot_from_CD`
  - output folder for manual reorganization of MS2 matched tables from mzCloud, mbank and hmdb: `/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/output/manual_compiled_MS2_mzCloud_mbank_hmdb_res`
- This branch is finally merged into `annotation-with-CD-exported tables` branch. 

### Summary
- this branch has two main improvments:
  - merge xcms-feature-table with Compound-Discoverer-exported table
  - reorganization of hmdb/mbank/mzCloud and render nice heatmaps, chromPeaks of those annotated features.
- In need of enhancement:
  - [ ] How to perform Spectra matching together with generating a table documenting all the matched precursor information, annotations and so on
    - The specificies can refer to the manual reorganized table (`/Users/gongm/Documents/projects/Rafi_Ahmed_exhaustedCD8Tcell/data/input/manual_compiled_MS2_mzCloud_mbank_hmdb_res/tables/HILICpos_manual_input_mzCloud_mbank_hmdb_matched_spectra_in_Rafi_Sign_list.csv`)
  - [ ] X1: figure out a way to put together hmdb/mbank table (generated in the previous step) & the mzCloud table, facilitating the X2: heatmap generations.
  - [ ] X4: XIC plot after chromPeak plotting to make sure the mz traces look good.

