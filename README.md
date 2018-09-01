### DIA-NN

DIA-NN - a fast and easy to use tool for processing data independent acquisition (DIA) proteomics data.  
DIA-NN implements deep neural networks to improve precursor ion identification.  
DIA-NN now also supports library-free search and spectral library generation.

**DIA-NN: Deep neural networks substantially improve the   
identification performance of Data-independent acquisition (DIA) in proteomics**  
Vadim Demichev, Christoph B. Messner, Kathryn S. Lilley, Markus Ralser  
https://doi.org/10.1101/282699

DIA-NN encompasses all stages of DIA-MS data processing in a single program.   
As input it takes raw data files and a spectral library / FASTA database.  
A report with protein and precursor ion quantities is produced. 
  
Manual: https://github.com/vdemichev/DiaNN/blob/master/DIA-NN%20GUI%20manual.pdf   
  
### New in version 1.5

Significantly improved identification performance.  
Automatic mass correction and mass accuracy inference.  
Automatic removal of interfering peptides.  
Increased speed and greatly reduced memory usage.  
Library-free search and spectral library generation from DIA data.  
GUI wrapper.  

### Installation

None required. Two executables are provided: DiaNN.exe (a command-line tool) and DIA-NN.exe (a GUI implemented as a wrapper for DiaNN.exe).

### Input files

Raw data files: Thermo .raw, .mzML or .dia (format used by DIA-NN to store spectra).  
Reading Thermo .raw files requires Thermo MS File Reader (https://thermo.flexnetoperations.com/control/thmo/login?nextURL=%2Fcontrol%2Fthmo%2Fdownload%3Felement%3D6306677) to be installed.   
The .mzML files should be centroided and contain data as spectra (e.g. SWATH) and not chromatograms.  

Spectral library: comma-separated (.csv) or tab-separated (.tsv or .txt) file.    
  
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

#### Example:
```
diann.exe --threads 4 --f run1 --f run2 --lib yeast.tsv --prefix C:\Data\ --ext .mzML --out run1_2.tsv    
```

For a full list of supported commands see the arguments() function in /src/diann.cpp.

### Building

A Visual C++ solution file is provided with the source code for building on Windows. It should also be possible to buld DIA-NN on Linux, as all the libraries it utilises are compatible with Linux. However, so far it has only been tested on Windows (versions 7 and 10).  

### Tutorial

This is a simple tutorial which covers generation of a spectral library from DIA data and its use to analyse other DIA runs. 
1. Download DIA-NN.exe (the GUI) and DiaNN.exe (the command line tool) to the same folder. 
2. Download "Fig2 HeLa-0-5h_MHRM_R01_T0.raw", "Fig2 HeLa-1h_MHRM_R01_T0.raw" and "uniprot_sprot_2014-12-11_HUMAN_ISOFORMS.fasta" from https://www.ebi.ac.uk/pride/archive/projects/PXD005573/files. 
3. (Optional) Download "ecolihumanyeast_concat_mayu_IRR_cons_openswath_64var_curated.csv" from https://www.ebi.ac.uk/pride/archive/projects/PXD002952/files. Rename the file extension to .tsv (as this is a tab-separated file). 
4. (Optional) Download http://www.peptideatlas.org/builds/human/201712/APD_Hs_all.fasta (browse http://www.peptideatlas.org/builds/ and download APD_Hs_all.fasta corresponding to the Jan 2018 Human build) to the same folder as DIA-NN.exe. 
5. Download and install MSFileReader (see above; to let DIA-NN access .raw files directly).  
6. Run DIA-NN.exe. In the "Files" panel click "Add raw data" and select "Fig2 HeLa-1h_MHRM_R01_T0.raw". In the "Library-free search" panel check "Use library-free search / generate spectral library" and (optional) select "ecolihumanyeast_concat_mayu_IRR_cons_openswath_64var_curated.tsv" as the Training Library. The latter allows to train peptide fragmentation and retention time predictors (see Manual). In the "Main Settings" panel click "Add FASTA" and select "uniprot_sprot_2014-12-11_HUMAN_ISOFORMS.fasta". (Optional) in the "Additional options" text box type "--fasta-filter APD_Hs_all.fasta". This will instruct DIA-NN to only look for peptides known to be detected in human samples, speeding up the search significantly. 
7. Click Run. DIA-NN will perform library-free analysis of the 1h HeLa run and create a spectral library (file "lib.tsv", 41571 precursor ions). Full analysis report will be saved to "quant.tsv". This will take 10-30 minutes on an average desktop.  
8. "Clear list" in the "Files" panel, then click "Add raw data" and select "Fig2 HeLa-0-5h_MHRM_R01_T0.raw". Uncheck "Use library-free search / generate spectral library" in the "Library-free search" panel. In the "Spectral library" field ("Main Settings" panel) type "lib.tsv". Clear the "Additional options" text box and press Run. In about a minute DIA-NN will process the 0.5h HeLa run and report 27418 precursors and 3378 proteins identified at 1% FDR. 






