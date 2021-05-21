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
- 