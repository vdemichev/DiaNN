
### DIA-NN

DIA-NN - a fast and easy to use tool for processing data-independent acquisition (DIA) proteomics data.  
DIA-NN implements deep neural networks to improve precursor ion identification.  
DIA-NN now also supports library-free search and spectral library generation.

Please cite   
**DIA-NN: neural networks and interference correction   
enable deep proteome coverage in high throughput**  
https://www.nature.com/articles/s41592-019-0638-x

When using DIA-NN for dia-PASEF analysis, please cite  
**High sensitivity dia-PASEF proteomics with DIA-NN and FragPipe**
https://www.biorxiv.org/content/10.1101/2021.03.08.434385v1  

DIA-NN encompasses all stages of DIA-MS data processing in a single program.   
As input it takes raw data files and a spectral library / FASTA database.  
A report with protein and precursor ion quantities is produced.

**Download**: https://github.com/vdemichev/DiaNN/releases/tag/1.7.16
(it's recommended to use the latest version - DIA-NN 1.7.16)  

**DIA-NN manual**: https://github.com/vdemichev/DiaNN/blob/master/DIA-NN%20GUI%20manual.pdf  
Please also check the commands listed below.   

**R package** with some useful functions for dealing with DIA-NN's reports: https://github.com/vdemichev/diann-rpackage

### New since version 1.7.0  

1. Deep learning-based generation of high quality spectral libraries *in silico* (Windows only). This boosts library-free search efficiency.   

2. Analysis logs (which also specify the DIA-NN version used and the search settings) are now saved automatically.  

3. Multi-tab GUI, ideal for dealing with multiple experiments in parallel.  

4. Experimental support for NIST .msp and SpectraST .sptxt libraries.   

5. MaxLFQ-based protein quantification and improved protein grouping.

6. Significantly higher numbers of proteins identified.     

7. Global (experiment-wide) precursor-level and protein-level FDR control (since 1.7.15).

8. Full dia-PASEF support (since 1.7.15).

### Installation

Download the installer (https://github.com/vdemichev/DiaNN/releases, click 'Assets') and run it. Make sure not to run the installer from a network drive.        
For .wiff support, first download and install ProteoWizard (http://proteowizard.sourceforge.net/download.html - choose the version (64-bit) that supports "vendor files"). Then copy all files with 'Clearcore' or 'Sciex' in their name (these will be .dll files) from the Proteowizard folder to the DIA-NN installation folder.   

### Input files

Raw data files: Sciex .wiff, Bruker .d, Thermo .raw, .mzML or .dia (format used by DIA-NN to store spectra).  
Reading Thermo .raw files requires Thermo MS File Reader (https://thermo.flexnetoperations.com/control/thmo/login?nextURL=%2Fcontrol%2Fthmo%2Fdownload%3Felement%3D6306677) to be installed. It is essential to use specifically the version by the link above (3.0 SP3).     
The .mzML files should be centroided and contain data as spectra (e.g. SWATH) and not chromatograms.  

Spectral library: support for comma-separated (.csv) or tab-separated (.tsv, .xls or .txt), .speclib (compact format used by DIA-NN), .sptxt (SpectraST, experimental) and .msp (NIST, experimental) files.    

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

### Building (DIA-NN 1.7.12 only)

**Windows**: A Visual C++ solution file is provided with the source code. Changing the SDK to a Windows 10 one in the project settings might be required. Tested on Windows 7 and 10. Configuration should be set to "Release" (source code for deep learning-based generation of spectral libraries (PyTorch) is not included, so this functionality is only supported in the official DIA-NN builds).      

**Linux** (GCC 7 or later required):		     		
```  		
git clone https://www.github.com/vdemichev/diann  		
cd diann/mstoolkit  		
make      		
```		  
Bash scripts for building without .mzML support (for this, uncomment "//#undef MSTOOLKIT" after "#ifdef LINUX") are also provided in the root directory.   		

### Tutorial

This is a simple tutorial which covers the generation of a spectral library from gas-phase fractionated DIA data and its use to analyse other DIA runs (using DIA-NN 1.7.11).  
1. Download all the files from https://osf.io/w5dr6/files/?view_only=00c8a68bfb824835b7fa304e31922ffa - these are 9 yeast runs on TripleTOF 6600 converted to .dia format: select "OSF Storage (Germany - Frankfurt)" and then choose "Download as zip" above (7.6 Gb total). Unpack the zip file.
2. Download the UniProt canonical *S.cerevisiae* proteome from https://www.uniprot.org/proteomes/UP000002311 (Under "Components" click "Download" and then "Go").  
3. Download the collection of peptides known to be detectable in yeast samples from http://www.peptideatlas.org/builds/: find the Yeast ("Build name" column) Mar 2013 ("Date" column) entry, right click on "APD_Sc_all.fasta" in the "Peptide Sequences" column and choose "Save link as".  
4. Launch DIA-NN. Click **Add raw data** in the **Input** panel, select the 8 gas-phase fractionation runs (with names 400-500, 495-600, ... 1095-1250). Click **Add FASTA** and select APD_Sc_all.fasta. Check **FASTA digest for library-free search / library generation** in the **Precursor ion generation** panel. Make sure **Deep learning-based spectra and RTs prediction** (same panel) as well as **Generate spectral library** (**Output** panel) are checked. Choose where to save the **Main output** file (e.g. report.tsv in the same folder). Choose where to save the **Output library** (e.g. lib.tsv in the same folder). Click **Add to pipeline** (just above the pipeline window).  

5. Copy-paste the **Output library** field to the **Spectral library** field in the **Input** panel. Uncheck **FASTA digest for library-free search / library generation**, **Deep learning-based spectra and RTs prediction** and **Generate spectral library**. Click **Clear list** next to the **Add raw data** button and then add pre_qc_1.wiff.dia (this is 23-min microflow gradient SWATH run on Sciex TripleTOF 6600). Click **Clear list** below the **Add FASTA** button, click **Add FASTA** and select the UniProt sequence database, check **Reannotate** (to add protein information from the sequence database to the spectral library). Click **Add to pipeline**.  Click **Execute** (below the pipeline window). On Intel i3-8350K, the whole pipeline (creation of the spectral library plus its use to analyse the 23-minute run) takes about 17 minutes and results in ~25000 precursors identified at 1% FDR.  



**Comments**:
1. DIA-NN can create libraries from any DIA/SWATH runs (not necessarily gas-phase fractionation runs). Sometimes even a single long-gradient run might be good enough.
2. Library-free search was used to create a spectral library. Regular spectral library-based search can also be used to generate new libraries. This is useful e.g. for refining spectra and retention times in a public library.  
3. The peptide list (APD_Sc_all) was used primarily to speed up the analysis (as this is a tutorial), it is OK to use the full sequence database instead.     

4. Each of the pipeline steps can be executed separately (select a step and click **Run**). Pipeline creation is also optional (can just click **Run** to execute the current configuration).  

### Acknowledgements
Many thanks to Christoph Messner, Spyros Vernardis, Kathryn Lilley and Markus Ralser as well as to the members of the Ralser lab (Francis Crick Institute and Charité - Universitätsmedizin Berlin) for their continuous support in improving DIA-NN and expanding its functionality.

Prosit https://www.nature.com/articles/s41592-019-0426-7 online service (part of ProteomeTools) https://www.proteomicsdb.org/prosit/ has been utilised when developing DIA-NN versions uploaded after 24/01/2020.

MaxLFQ https://doi.org/10.1074/mcp.M113.031591 protein quantification algorithm is used in DIA-NN versions uploaded after 02/04/2020.  
