### DIA-NN

DIA-NN - a fast and easy to use tool for processing data-independent acquisition (DIA) proteomics data.  
DIA-NN implements deep neural networks to improve precursor ion identification.  
DIA-NN now also supports library-free search and spectral library generation.

**DIA-NN: neural networks and interference correction   
enable deep proteome coverage in high throughput**  
Vadim Demichev, Christoph B. Messner, Spyros I. Vernardis, Kathryn S. Lilley, Markus Ralser  
https://www.nature.com/articles/s41592-019-0638-x

DIA-NN encompasses all stages of DIA-MS data processing in a single program.   
As input it takes raw data files and a spectral library / FASTA database.  
A report with protein and precursor ion quantities is produced. 
  
**Download**: https://github.com/vdemichev/DiaNN/releases   
GUI manual: https://github.com/vdemichev/DiaNN/blob/master/DIA-NN%20GUI%20manual.pdf   

### New in version 1.7.0  

Considerably improved identification performance     
Tunable quantification algorithms        
It is no longer necessary to define unknown modifications featured in the spectral library      
Experimental support for SILAC and absolute quantification with spike-in standards    

### Installation

None required (for .raw, .mzML and .dia processing). Two executables are provided: DiaNN.exe (a command-line tool) and DIA-NN.exe (a GUI implemented as a wrapper for DiaNN.exe).  
For .wiff support, first download and install ProteoWizard (http://proteowizard.sourceforge.net/download.html - choose the version that supports "vendor files"), then use the installer provided (DIA-NN-Setup.msi) and specify the Proteowizard directory (e.g. C:\Program Files (x86)\ProteoWizard\ProteoWizard \[version\], where \[version\] is the ProteoWizard version number, e.g. 3.0.11537) as the installation directory.  

### Input files

Raw data files: Sciex .wiff, Thermo .raw, .mzML or .dia (format used by DIA-NN to store spectra).  
Reading Thermo .raw files requires Thermo MS File Reader (https://thermo.flexnetoperations.com/control/thmo/login?nextURL=%2Fcontrol%2Fthmo%2Fdownload%3Felement%3D6306677) to be installed. It is essential to use specifically the version by the link above (3.0 SP3).     
The .mzML files should be centroided and contain data as spectra (e.g. SWATH) and not chromatograms.  

Spectral library: comma-separated (.csv) or tab-separated (.tsv, .xlx or .txt) file.    
  
### Command-line tool usage
```
diann.exe [commands]  
```
Commands can be supplied in arbitrary order.     

#### Required commands:  
```
--f <data file> 
```
Specifies a data file to be processed. Use --f for each file to be processed. 
```
--lib <spectral library file>
```
Specifies the spectral library. Example:
```
diann.exe --f run1.mzML --f run2.mzML --lib yeast.tsv  
```
#### Auxiliary commands:  
```
--cfg <config file> 
```
Specifies a file with a set of commands.
```
--dir <folder> 
```
Specifies a folder containing raw files to be processed. All files in the folder must be in .raw, .mzML or .dia format.  
```
--threads <thread number> 
```
Specifies the number of CPU threads to be used.  
```
--convert
```
With this option DIA-NN converts .mzML files (specified using the --f command) to the .dia format. Unlike .mzML files, .dia files can be loaded quickly (seconds), so it is recommended to convert files that are going to be analysed multiple times.    
```
--ext <string>
```
Add a string to the end of each file name (specified with --f).  
```
--prefix <string>
```
Add a string to the beginning of each file name.  
```
--out <output file> 
```
Specifies the output file; by default, the output is saved to quant.tsv in the current working directory.


#### More commands: 


```
--clear-mods : Modification names specified in the spectral library will be used for annotating output. Only suitable for library-based analysis. 
--threads : Thread number. 
--fasta-search : Library-free search enabled. 
--pg-level : this determines which peptides are considered 'proteotypic' and thus affects protein FDR calculation. 
--no-batch-mode : Batch mode disabled. 
--verbose : verbose
--export-windows : export infomation, default true. 
--export-library : export infomation, default true. 
--export-decoys : export infomation, default true. 
--prosit : prosit, default true.
--vis : XICs for precursors corresponding peptides will be saved. 
--cal-info : Save calculate Info, default true. 
--compact-report : Extended report, default false. 
--no-isotopes : Isotopologue chromatograms will not be used. 
--no-ms2-range : MS2 range inference will not be performed. 
--min-peak : Minimum peak height
--no-cal-filter : Peptides with modifications that can cause interferences with isotopologues will not be filtered out for mass calibration. 
--no-nn-filter : Peptides with modifications that can cause interferences with isotopologues will be used for neural network training. 
--nn-cross-val : Neural network cross-validation will be used to tackle potential overfitting. 
--guide-classifier : A separate classifier for the guide library will be used. 
--int-removal : Number of interference removal iterations. 
--int-margin : Interference correlation margin. 
--strict-int-removal : Potentially interfering peptides with close (but not the same) elution times will also be discarded. 
--reverse-decoys : Decoys will be generated using the pseudo-reverse method. 
--force-frag-rec : Decoys will be generated only for precursors with all library fragments recognised. 
--max-rec-charge : Fragment recognition module will consider charges up. 
--max-rec-loss : Fragment recognition module will consider losses with the index.
--gen-spec-lib : A spectral library will be generated. 
--lib-gen-direct-q : When generating a spectral library, run-specific q-values will be used instead of profile q-values. 
--save-original-lib : All entries from the library provided will be saved to the newly generated library. 
--fast-wiff : Custom fast centroiding will be used when processing .wiff files (WARNING: experimental).
--dir : file path. 
--lib : library. 
--fasta : fasta file.
--fasta-filter : fasta filter. 
--ref : reference. 
--out : output file.
--out-gene : out genes.
--qvalue : qvalue
--protein-qvalue : Output will be filtered at protein-level FDR. 
--no-prot-inf : Protein inference will not be performed. 
--no-swissprot : SwissProt proteins will not be prioritised for protein inference. 
--force-swissprot: Only SwissProt proteins will be considered in library-free search. 
--species-genes : Species suffix will be added to gene annotation; this affects proteotypicity definition . 
--duplicate-proteins : Duplicate proteins in FASTA files will not be skipped. 
--out-lib : out library. 
--learn-lib : learn library. 
--out-measured-rt : When generating a spectral library without a guide library but with a training library, iRT values (Tr_recalibrated) will correspond to the measured retention times.
--library-headers : library headers. 
--output-headers : output headers. 
--mod : modification. 
--fixed-mod : modification will be considered as fixed. 
--var-mod : modification will be considered as variable. 
--ref-cal : Reference peptides will be used for calibration. 
--gen-ref : A library of reference peptides will be generated. 
--window : Scan window.
--cut-after : In silico digest will include cuts after amino acids : . 
--no-cut-before : In silico digest will not include cuts before amino acids: . 
--min-pep-len : Min peptide length. 
--max-pep-len : Max peptide length. 
--min-fr-corr : Minimum fragment profile correlation for the inclusion into the spectral library. 
--min-gen-fr : Minimum number of fragments for library generation. 
--min-pr-mz : Min precursor m/z. 
--max-pr-mz : Max precursor m/z.
--min-fr-mz : Min fragment m/z. 
--max-fr-mz : Max fragment m/z. 
--max-fr : Maximum number of fragments. 
--min-fr : Minimum number of fragments for library export. 
--min-search-fr : Minimum number of fragments required for a precursor to be searched .
--missed-cleavages : Maximum number of missed cleavages. 
--unimod4 : Cysteine carbamidomethylation enabled as a fixed modification. 
--unimod35 : Methionine oxidation enabled as a variable modification. 
--var-mods : Maximum number of variable modifications. 
--met-excision : N-terminal methionine excision enabled. 
--no-rt-window : Full range of retention times will be considered. 
--disable-rt : All RT-related scores disabled. 
--min-rt-win: Minimum acceptable RT window scale. 
--no-window-inference : Scan window inference disabled. 
--individual-windows : Scan windows will be inferred separately for different runs. 
--individual-mass-acc : Mass accuracy will be determined separately for different runs. 
--individual-reports : Reports will be generated separately for different runs (in the respective folders). 
--no-stats : Run statstics infomation, default false.
--convert : MS data files will be converted to .dia format.
--out-dir : output directory.
--remove-quant : .quant files will be removed when the analysis is finished. 
--no-quant-files : .quant files will not be saved to the disk. 
--use-rt : Existing .quant files will be used for RT profiling.
--use-quant : Existing .quant files will be used. 
--quant-only : Quantification will be performed anew using existing identification info. 
--report-only : Report will be generated using .quant files. 
--iter : Number of iterations . 
--profiling-qvalue : RT profiling q-value threshold. 
--quant-qvalue : Q-value threshold for cross-run quantification. 
--protein-quant-qvalue : Precursor Q-value threshold for protein quantification. 
--top : precursors will be used for protein quantification in each run . 
--out-lib-qvalue : Q-value threshold for spectral library generation . 
--rt-profiling : RT profiling enabled. 
--prefix : prefix added to input file names. 
--ext : extension added to input file names. 
--lc-all-scores : All scores will be used by the linear classifier (not recommended). 
--peak-center : Fixed-width center of each elution peak will be used for quantification. 
--peak-boundary : Peak boundary intensity factor. 
--standardisation-scale : Standardisation scale.
--no-ifs-removal : Interference removal from fragment elution curves disabled. 
--no-fr-selection : Cross-run selection of fragments for quantification disabled (not recommended) . 
--no-fr-exclusion : Exclusion of fragments shared between heavy and light labelled peptides from quantification disabled . 
--peak-translation : Translation of retention times between peptides within the same elution group enabled. 
--no-standardisation : Scores will not be standardised for neural network training. 
--no-nn : Neural network classifier disabled. 
--nn-iter : Neural network classifier will be used starting from the interation number. 
--nn-bagging : Neural network bagging . 
--nn-epochs : Neural network epochs number. 
--nn-learning-rate : Neural network learning rate. 
--nn-reg : Neural network regularisation. 
--nn-hidden : Number of hidden layers. 
--mass-acc-cal : Calibration mass accuracy. 
--fix-mass-acc : Force mass accuracy, default true.
--mass-acc : Global mass accuracy. 
--mass-acc-ms1 : Global MS1 accuracy.
--gen-acc : Fragmentation spectrum generator accuracy. 
--min-corr : Only peaks with correlation sum exceeding Min Ms1 correlation will be considered.
--corr-diff : Peaks with correlation sum below from maximum will not be considered. 
--peak-apex : Peaks apex height. 
--all-peaks : The number of putative elution peaks considered in library-free mode will not be reduced to decrease RAM usage. 
--norm-qvalue : Q-value threshold for cross-run normalisation. 
--norm-fraction : Global normalisation peptides fraction. 
--norm-radius : Local normalisation radius. 
--global-norm : Median-based local normalisation disabled. 
--q1-cal : Q1 calibration enabled. 
--no-calibration : Mass calibration disabled. 
--mass-cal-bins : Maximum number of mass calibration bins. 
--min-cal : Minimum number of precursors identified at 10% FDR used for calibration. 
--min-class : Minimum number of precursors identified at 10% FDR used for linear classifier training. 
--scanning-swath : All runs will be analysed as Scanning SWATH runs. 
--regular-swath : All runs will be analysed as regular SWATH runs. 
--no-q1 : Q1 scores disabled. 
--use-q1 : Q1 scores will be used for regular SWATH runs. 

```


#### Example:
```
diann.exe --threads 4 --f run1 --f run2 --lib yeast.tsv --prefix C:\Data\ --ext .mzML --out run1_2.tsv    
```
This following command is supported in development versions uploaded on 20/12/2019 or later.     
```   
--vis <N>,<Peptide 1>,<Peptide 2>,... 
```
Specifies a number of peptides for which chromatograms of length >= N will be extracted and saved as a text table. These allow for PSM visualisation using third-party software (e.g. R or Python).    

#### Example (will save XICs to report.XIC.tsv):
```
diann.exe --f run1.mzML --f run2.mzML --lib yeast.tsv --out report.tsv --vis 20,KVYPDVLYTSK,TAIEGSYIDK,DSATHELTK       
```

For a full list of supported commands see the arguments() function in /src/diann.cpp.

### Building

**Windows**: A Visual C++ solution file is provided with the source code. Changing the SDK to a Windows 10 one in the project settings might be required. Tested on Windows 7 and 10.      
**Linux** (GCC 7 or later required):		     		
```  		
git clone https://www.github.com/vdemichev/diann  		
cd diann/mstoolkit  		
make      		
```		  
Bash scripts for building without .mzML support (for this, uncomment "//#undef MSTOOLKIT" after "#ifdef LINUX") are also provided in the root directory.   		
   		
### Tutorial

This is a simple tutorial which covers the generation of a spectral library from DIA data and its use to analyse other DIA runs (using DIA-NN 1.7.1). The tutorial is only meant to illustrate how to use the GUI; for the optimal performance, spectral libraries should be created from multiple gas-phase fractionation DIA runs.   		
1. Place DIA-NN.exe (the GUI) and DiaNN.exe (the command line tool) in the same folder.  
2. Download "Fig2 HeLa-0-5h_MHRM_R01_T0.raw", "Fig2 HeLa-1h_MHRM_R01_T0.raw" and "uniprot_sprot_2014-12-11_HUMAN_ISOFORMS.fasta" from https://www.ebi.ac.uk/pride/archive/projects/PXD005573/files.  
3. (Optional) Download "ecolihumanyeast_concat_mayu_IRR_cons_openswath_64var_curated.csv" from https://www.ebi.ac.uk/pride/archive/projects/PXD002952/files. Rename the file extension to .tsv (as this is a tab-separated file).  
4. (Optional) Download http://www.peptideatlas.org/builds/human/201712/APD_Hs_all.fasta (browse http://www.peptideatlas.org/builds/ and download APD_Hs_all.fasta corresponding to the Jan 2018 Human build) to the same folder as DIA-NN.exe.    
5. Download and install MSFileReader (see above; to let DIA-NN access .raw files directly).  
6. In the "Input" panel click "Add raw data" and select "Fig2 HeLa-1h_MHRM_R01_T0.raw". Also in the "Input" panel, click "Add FASTA" and select "uniprot_sprot_2014-12-11_HUMAN_ISOFORMS.fasta". (Optional) in the "Additional options" text box ("Output" panel) type "--fasta-filter APD_Hs_all.fasta". This will instruct DIA-NN to only look for peptides known to be detected in human samples, speeding up the search significantly. In the "Library-free search" panel check "Use library-free search" and (optional) select "ecolihumanyeast_concat_mayu_IRR_cons_openswath_64var_curated.tsv" as the Training Library. The latter allows to train peptide fragmentation and retention time predictors (see Manual). To speed up analysis, set the number of missed cleavages to 0. Click "Add to pipeline".   
7. "Clear list" at the top of the "Input" panel, then click "Add raw data" and select "Fig2 HeLa-0-5h_MHRM_R01_T0.raw". Uncheck "Use library-free search" in the "Library-free search" panel. In the "Spectral library" field ("Input" panel) type "lib.tsv". Clear the "Additional options" text box, click "Add to pipeline".    
8. Click "Execute". First, DIA-NN will generate a spectral library (consisting of 42310 precursors) from the 1h run (this will take 10-30 minutes on an average desktop). Afterwards, DIA-NN will process the 0.5h HeLa run and report 31770 precursors and 3859 proteins identified at 1% FDR (this will take about a minute).     
