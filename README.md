### DIA-NN

DIA-NN - a fast and easy to use tool for processing data-independent acquisition (DIA) proteomics data.  
DIA-NN implements deep neural networks to improve precursor ion identification.  
DIA-NN now also supports library-free search and spectral library generation.

**DIA-NN: Neural networks and interference correction  	
enable deep coverage in high-throughput proteomics**  
Vadim Demichev, Christoph B. Messner, Spyros I. Vernardis, Kathryn S. Lilley, Markus Ralser  
https://doi.org/10.1101/282699

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

#### Example:
```
diann.exe --threads 4 --f run1 --f run2 --lib yeast.tsv --prefix C:\Data\ --ext .mzML --out run1_2.tsv    
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

This is a simple tutorial which covers the generation of a spectral library from DIA data and its use to analyse other DIA runs (using DIA-NN 1.7.0). The tutorial is only meant to illustrate how to use the GUI; for the optimal performance, spectral libraries should be created from multiple gas-phase fractionation DIA runs.   		
1. Place DIA-NN.exe (the GUI) and DiaNN.exe (the command line tool) in the same folder.  
2. Download "Fig2 HeLa-0-5h_MHRM_R01_T0.raw", "Fig2 HeLa-1h_MHRM_R01_T0.raw" and "uniprot_sprot_2014-12-11_HUMAN_ISOFORMS.fasta" from https://www.ebi.ac.uk/pride/archive/projects/PXD005573/files.  
3. (Optional) Download "ecolihumanyeast_concat_mayu_IRR_cons_openswath_64var_curated.csv" from https://www.ebi.ac.uk/pride/archive/projects/PXD002952/files. Rename the file extension to .tsv (as this is a tab-separated file).  
4. (Optional) Download http://www.peptideatlas.org/builds/human/201712/APD_Hs_all.fasta (browse http://www.peptideatlas.org/builds/ and download APD_Hs_all.fasta corresponding to the Jan 2018 Human build) to the same folder as DIA-NN.exe.    
5. Download and install MSFileReader (see above; to let DIA-NN access .raw files directly).  
6. In the "Input" panel click "Add raw data" and select "Fig2 HeLa-1h_MHRM_R01_T0.raw". Also in the "Input" panel, click "Add FASTA" and select "uniprot_sprot_2014-12-11_HUMAN_ISOFORMS.fasta". (Optional) in the "Additional options" text box ("Output" panel) type "--fasta-filter APD_Hs_all.fasta". This will instruct DIA-NN to only look for peptides known to be detected in human samples, speeding up the search significantly. In the "Library-free search" panel check "Use library-free search" and (optional) select "ecolihumanyeast_concat_mayu_IRR_cons_openswath_64var_curated.tsv" as the Training Library. The latter allows to train peptide fragmentation and retention time predictors (see Manual). To speed up analysis, set the number of missed cleavages to 0. Click "Add to pipeline".   
7. "Clear list" at the top of the "Input" panel, then click "Add raw data" and select "Fig2 HeLa-0-5h_MHRM_R01_T0.raw". Uncheck "Use library-free search" in the "Library-free search" panel. In the "Spectral library" field ("Input" panel) type "lib.tsv". Clear the "Additional options" text box, click "Add to pipeline".    
8. Click "Execute". First, DIA-NN will generate a spectral library (consisting of 42310 precursors) from the 1h run (this will take 10-30 minutes on an average desktop). Afterwards, DIA-NN will process the 0.5h HeLa run and report 31770 precursors and 3859 proteins identified at 1% FDR (this will take about a minute).     
