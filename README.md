
<p align="center" style="margin-bottom: 0px !important;">
  <img src="https://raw.githubusercontent.com/vdemichev/DiaNN/master/GUI/DIA-NN%20art%20-%20128x128.png" width="100" height="100">
</p>
<h1 align="center" style="margin-top: -10px; font-size: 20px">DIA-NN</h1>

DIA-NN - a universal software suite for data-independent acquisition (DIA) proteomics data processing. Conceived at the University of Cambridge, UK, in the laboratory of Kathryn Lilley (Cambridge Centre for Proteomics), DIA-NN opened a new chapter in proteomics, introducing a number of algorithms which enabled reliable, robust and quantitatively accurate large-scale experiments using high-throughput methods. DIA-NN is currently being further developed in the laboratory of Vadim Demichev at the Charité (University Medicine Berlin, Germany).

DIA-NN is built on the following principles:    
- **Reliability** achieved via stringent statistical control
- **Robustness** achieved via flexible modelling of the data and automatic parameter selection
- **Reproducibility** promoted by thorough recording of all analysis steps
- **Ease of use**: high degree of automation, an analysis can be set up in several mouse clicks, no bioinformatics expertise required
- **Powerful tuning options** to enable unconventional experiments
- **Scalability and speed**: up to 1000 mass spec runs processed per hour

**Download**: https://github.com/vdemichev/DiaNN/releases/tag/1.9
(it's recommended to use the latest version - DIA-NN 1.9)  

<img src="https://raw.githubusercontent.com/vdemichev/DiaNN/master/GUI/GUI%20window.png"></br>  

Please cite:   
**DIA-NN: neural networks and interference correction   
enable deep proteome coverage in high throughput** [Nature Methods, 2020](https://www.nature.com/articles/s41592-019-0638-x)

Using DIA-NN for the analysis of post-translation modifications (PTMs), such as phosphorylation or ubiquitination: **Time-resolved in vivo ubiquitinome profiling by DIA-MS reveals USP7 targets on a proteome-wide scale** [Nature Communications, 2021](https://www.nature.com/articles/s41467-021-25454-1)

Using DIA-NN's ion mobility module for timsTOF data analysis or using DIA-NN in combination with FragPipe-generated spectral libraries: **dia-PASEF data analysis using FragPipe and DIA-NN for deep proteomics of low sample amounts** [Nature Communications, 2022](https://www.nature.com/articles/s41467-022-31492-0)

Using DIA-NN for the analysis of multiplexed samples (SILAC, mTRAQ, etc): **Increasing the throughput of sensitive proteomics by plexDIA** [Nature Biotechnology, 2022](https://www.nature.com/articles/s41587-022-01389-w)

Using DIA-NN as part of the CysQuant workflow: **CysQuant: Simultaneous quantification of cysteine oxidation and protein abundance using data dependent or independent acquisition mass spectrometry** [Redox Biology, 2023](https://doi.org/10.1016/j.redox.2023.102908)

Using DIA-NN's QuantUMS module for quantification: **QuantUMS: uncertainty minimisation enables confident quantification in proteomics** [biorxiv](https://www.biorxiv.org/content/10.1101/2023.06.20.545604v1)

Using DIA-NN to process Slice-PASEF data: **Slice-PASEF: fragmenting all ions for maximum sensitivity in proteomics** [biorxiv](https://www.biorxiv.org/content/10.1101/2022.10.31.514544v1)

**Other key papers**  
- Using DIA-NN for large-scale plasma & serum proteomics:  
[Cell Systems, 2020](https://doi.org/10.1016/j.cels.2020.05.012) and [Cell Systems, 2021](https://doi.org/10.1016/j.cels.2021.05.005)
- Ultra-fast proteomics with DIA-NN and Scanning SWATH:   
[Nature Biotechnology, 2021](https://www.nature.com/articles/s41587-021-00860-4)

**R package** with some useful functions for dealing with DIA-NN's output reports: https://github.com/vdemichev/diann-rpackage

**Visualisation** of peptide positions in the protein: https://github.com/MannLabs/alphamap (AlphaMap by Mann lab)

**Notes and discussions** on proteomics in general and the use of DIA-NN: https://github.com/vdemichev/DiaNN/discussions/categories/dia-proteomics-in-detail (this section will be further expanded).

### Table of Contents
**[Installation](#installation)**<br>
**[Getting started](#getting-started)**<br>
**[Raw data formats](#raw-data-formats)**<br>
**[Spectral library formats](#spectral-library-formats)**<br>
**[Output](#output)**<br>
**[Library-free search](#library-free-search)**<br>
**[Creation of spectral libraries](#creation-of-spectral-libraries)**<br>
**[Match-between-runs](#match-between-runs)**<br>
**[Changing default settings](#changing-default-settings)**<br>
**[Command-line tool](#command-line-tool)**<br>
**[Visualisation](#visualisation)**<br>
**[Automated pipelines](#automated-pipelines)**<br>
**[PTMs and peptidoforms](#ptms-and-peptidoforms)**<br>
**[Multiplexing using plexDIA](#multiplexing-using-plexdia)**<br>
**[GUI settings reference](#gui-settings-reference)**<br>
**[Command-line reference](#command-line-reference)**<br>
**[Main output reference](#main-output-reference)**<br>
**[Frequently asked questions (FAQ)](#frequently-asked-questions)**<br>
**[Support](#support)**<br>

### Installation

On **Windows**, download the .exe installer and run it. Make sure not to run the installer from a network drive. It is recommended to install DIA-NN into the default folder suggested by the installer.   

On **Linux**, download the .deb or .rpm package. Install it using the preferred way for your Linux distribution, e.g. using gdebi for .deb installation on Ubuntu/Debian-based distributions. DIA-NN will be installed to the /usr/diann/[version number]/ folder. The Linux system must have glibc 2.27 or later (for example, Ubuntu 18.04 or CentOS 8 and later versions are fine). Note: a Linux release of the latest DIA-NN 1.9 will follow shortly, this section will be updated once it is available.

It is also possible to run DIA-NN on Linux using **Wine** 6.8 or later.

### Getting Started

DIA mass spectrometry data can be analysed in two ways: by searching against a sequence database (library-free mode), or by using a "spectral library" - a set of known spectra and retention times for selected peptides. We discuss in detail when to use each of these approaches in the [Library-free search](#library-free-search) section. For both kinds of analyses, using DIA-NN is very simple:

1. Click **Raw** (in the **Input** pane), select your raw mass spectrometry data files. See [Raw data formats](#raw-data-formats) for information on supported formats.  
2. Click **Add FASTA**, add one or more sequence databases in UniProt format.
3. If you want to use a spectral library, click **Spectral library** and select the library. Alternatively, for library-free analysis, select **FASTA digest for library-free search/library generation** (in the **Precursor ion generation** pane).
4. Specify **Main output** file name in the **Output** pane and click **Run**.  
5. If you kept 'report.tsv' as the main output (located, by default, in the DIA-NN installation folder), it will contain the list of all precursor ions identified, along with different kinds of quantities, quality metrics and annotations. The output file report.pg_matrix.tsv will contain protein group quantities, report.gg_matrix.tsv - gene group quantities, report.pr_matrix.tsv - precursor ion quantities.

Now, the above information is sufficient for one to start using DIA-NN, it's indeed this easy! The rest of this Documentation might be helpful, but is not essential for 99% of the projects.

The above is how to run DIA-NN with default settings, and these yield optimal or almost optimal performance for most experiments. In some cases, however, it is better to adjust the settings, see [Changing default settings](#changing-default-settings) for more details.

DIA-NN also offers powerful tuning options for fancy experiments. DIA-NN is implemented as a user-friendly graphical interface that automatically invokes a command-line tool. But the user can also pass options/commands to the command-line tool directly, via the **Additional options** text box in the interface. All these options start with a double dash -- followed by the option name and, if applicable, some parameters to be set. So if you see some option/command with -- in its name mentioned in this Documentation, it means this command is meant to be typed in the **Additional options** text box.

### Raw data formats

Formats supported: Sciex .wiff, Bruker .d, Thermo .raw, .mzML and .dia (format used by DIA-NN to store spectra). Conversion from any supported format to .dia is possible. When running on Linux (native builds, not Wine), only .d, .mzML, and .dia data are supported.  

For .wiff support, download and install [ProteoWizard](http://proteowizard.sourceforge.net/download.html) - choose the version (64-bit) that supports "vendor files"). Then copy all files with 'Clearcore' or 'Sciex' in their name (these will be .dll files) from the ProteoWizard folder to the DIA-NN installation folder (the one which contains diann.exe, DIA-NN.exe and a bunch of other files).  

Reading Thermo .raw files requires [Thermo MS File Reader](https://thermo.flexnetoperations.com/control/thmo/login?nextURL=%2Fcontrol%2Fthmo%2Fdownload%3Felement%3D6306677) to be installed. It is essential to use specifically the version by the link above (3.0 SP3).

.mzML files should be centroided and contain data as spectra (e.g. SWATH/DIA) and not chromatograms.

<details>
  <summary>Technology support</summary>

- DIA and SWATH are supported  
- Acquisition schemes with overlapping windows are supported
- Gas-phase fractionation is supported
- Scanning SWATH is supported
- dia-PASEF/py-diAID is supported
- Slice-PASEF is supported (add --tims-scan to Additional options)
- midia-PASEF and Synchro-PASEF are supported (add --tims-scan to Additional options), but DIA-NN currently does not benefit from Q1 dimension deconvolution
- Orbitrap Astral is supported
- FAIMS with constant CV is supported
- FAIMS with multiple CVs is supported after splitting the runs, see [here](https://github.com/vdemichev/DiaNN/issues/67)
- BoxCar-DIA is supported, but DIA-NN has not been optimised for it
- Bruker Impact II DIA data are supported after conversion to .mzML
- multiplexing with non-isobaric tags and SILAC is supported
- MSX-DIA is not supported
</details>

<details>
  <summary>Conversion</summary>

Many mass spec formats, including those few that are not supported by DIA-NN directly, can be converted to .mzML using the MSConvertGUI application from [ProteoWizard](http://proteowizard.sourceforge.net/download.html). This works for all supported formats except Bruker .d and SCIEX Scanning SWATH - these need to be accessed by DIA-NN directly. The following MSConvert settings must be used for conversion:

<img src="https://raw.githubusercontent.com/vdemichev/DiaNN/master/GUI/MSConvert.png"></br>
</details>

### Spectral library formats

DIA-NN supports comma-separated (.csv) or tab-separated (.tsv, .xls or .txt), .speclib (compact format used by DIA-NN), .sptxt (SpectraST, experimental) and .msp (NIST, experimental) library files. Important: the library must not contain non-fragmented precursor ions as 'fragments': each fragment ion must actually be produced by the peptide backbone fragmentation.

<details>
  <summary>In detail</summary>

Libraries in the PeakView format as well as libraries produced by FragPipe, TargetedFileConverter (part of OpenMS), exported from Spectronaut (Biognosys) in the .xls format or generated by DIA-NN itself are supported “as is”.  

For .tsv/.xls/.txt libraries generated by other means, DIA-NN might require the header names to be specified (separated by commas) (for the columns it requires) using the --library-headers command. Use the * symbol instead of the name of a header to keep its recognition automatic. See below the descriptions of the respective columns (in the order the headers need to be specified).

Required columns:
- **Modified & labelled peptide sequence**
- **Precursor charge**
- **Precursor m/z**
- **Reference retention time** - arbitrary RT scale can be used
- **Fragment ion m/z**
- **Relative intensity of the fragment ion**

It is strongly recommended that columns containing the following are also present in the library:  
- **Protein IDs** - identifiers for the protein isoforms
- **Protein names**
- **Gene names**
- **Proteotypicity** - a column containing 0/1 values, depending on whether the peptide in question is 'proteotypic', that is specific to a particular protein isoform, protein name or gene
- **Decoy** - Indicates whether the peptide is a decoy. If there are decoy peptides in the library, DIA-NN uses these and does not generate its own decoys. It is strongly recommended **not** to include any decoy peptides in the library.
- **Fragment ion charge**
- **Fragment ion type** - either y or b; for x and z fragments also specify fragment type as y, and for a and c - as b
- **Fragment series number**
- **Fragment neutral loss type**
- **Q-value**
- **Elution group identifier** - if not specified, DIA-NN will infer elution groups automatically; not needed for most workflows
- **Exclude fragment indicator** - a column containing 0/1 values, with 1 meaning that the fragment ion should not be used for quantification; not needed for most workflows
- **Ion Mobility** - 1/K0 value for the precursor, arbitrary IM scale can be used

For example, a --library-headers command which specifies all column names except for the 'Decoy' column can look like this:

<nobr>--library-headers ModifiedPeptide</nobr>,PrecursorCharge,PrecursorMz,Tr_recalibrated,ProductMz,LibraryIntensity,UniprotID,ProteinName,Genes,Proteotypic,*,FragmentCharge,FragmentType,FragmentSeriesNumber,FragmentLossType,QValue,ExcludeFromAssay,IonMobility

Use --sptxt-acc to set the fragment filtering mass accuracy (in ppm) when reading .sptxt/.msp libraries.

MaxQuant msms.txt can also be used (experimental) as a spectral library in DIA-NN, although fixed modifications might not be read correctly.

DIA-NN can **convert** any library it supports into its own .tsv format. For this, click **Spectral library** (**Input** pane), select the library you want to convert, select the **Output library** file name (**Output** pane), click **Run**. If you use some exotic library format, it's a good idea to convert it to DIA-NN's .tsv and then examine the resulting library to see if the contents make sense.

All .tsv/.xls/.txt/.csv libraries are just simple tables with human-readable data, and can be explored/edited, if necessary, using Excel or (ideally) R/Python.

Importantly, when any library is being saved to this kind of text format, all numbers are rounded using certain decimal precision, meaning that they might not be exactly the same as in the original library (there might be a tiny difference). Thus, although the performance when analysing using a converted library will be comparable, the results will not match exactly.
</details>

### Output

The **Output** pane allows to specify where the output should be saved as well as the file names for the main output report and (optionally) the output spectral library. DIA-NN uses these file names to derive the names of all of its output files. Below one can find information on different types of DIA-NN output. For most workflows one only needs the main report (for analysis in R or Python - recommended) or the matrices (simplified output for MS Excel).

<details>
  <summary>Main report</summary>

A text table containing precursor and protein IDs, as well as plenty of associated information. Most column names are self-explanatory, and the full reference can be found in [Main output reference](#main-output-reference). The following keywords are used when naming columns:
- **PG** means protein group
- **GG** means gene group
- **Quantity** means non-normalised quantity
- **Normalised** means normalised quantity
- **MaxLFQ** means normalised protein quantity calculated using the [MaxLFQ](https://www.mcponline.org/article/S1535-9476(20)33310-7/fulltext) algorithm - it is strongly recommended to use these MaxLFQ quantities and not the regular quantities (also reported by DIA-NN)
- **Global** refers to a global q-value, that is calculated for the entire experiment
- **Lib** refers to the respective value saved in the spectral library, e.g. Lib.Q.Value means q-value for the respective library precursor

**Note:** since version 1.9, DIA-NN produces a report in the Apache .parquet format. This is a compressed text table format (~10x size reduction) that can be loaded in a single line of code using R 'arrow' package or Python 'pyarrow' package. Most of new functionality (introduced in DIA-NN 1.9) is only reflected in the parquet report, so it is recommended to use it instead of the legacy .tsv report in all cases, while the .tsv report is still generated only for compatibility with old analysis workflows. The generation of the legacy .tsv report can be turned off with --no-main-report. In addition to using R or Python, you can also view .parquet files with the [TAD Viewer](https://www.tadviewer.com/).

</details>

<details>
  <summary>Matrices</summary>

These contain normalised MaxLFQ quantities for protein groups ('pg_matrix'), gene groups ('gg_matrix'), unique genes ('unique_genes_matrix'; i.e. genes identified and quantified using only proteotypic, that is gene-specific, peptides) as well as normalised quantities for precursors ('pr_matrix'). They are filtered at 1% FDR, using global q-values for protein groups and both global and run-specific q-values for precursors. Additional 5% run-specific protein-level FDR filter is applied to the protein matrices, use --matrix-spec-q to adjust it. Sometimes DIA-NN will report a zero as the best estimate for a precursor or protein quantity. Such zero quantities are omitted from protein/gene matrices.
</details>

<details>
  <summary>Stats report</summary>

Contains a number of QC metrics which can be used for data filtering, e.g. to exclude failed runs or as a readout for method optimisation. Note that the number of proteins reported here corresponds to the number of unique proteins (i.e. identified with proteotypic precursors) in a given run at 1% unique protein q-value. This number can be reproduced from the main report generated using precursor FDR threshold of 100% and filtered using Protein.Q.Value <= 0.01 & Proteotypic == 1. What's counted as 'protein' here depends on the 'Protein inference' setting.    

</details>

<details>
  <summary>PDF report</summary>

A visualisation of a number of QC metrics, based on the main report as well as the stats report. The PDF report should be used only for quick preliminary assessment of the data and should not be used in publications.

</details>

<details>
  <summary>Flexible reanalysis</summary>

The **Output** pane allows to control how to handle the '.quant files'. Now, to explain what these are, let us consider how DIA-NN processes raw data. It first performs the computationally-demanding part of the processing separately for each individual run in the experiment, and saves the identifications and quantitative information to a separate .quant file. Once all the runs are processed, it collects the information from all the .quant files and performs some cross-run steps, such as global q-value calculation, protein inference, calculation of final quantities and normalisation. This allows DIA-NN to be used in a very flexible manner. For example, you can stop the processing at any moment, and then resume processing starting with the run you stopped at. Or you can remove some runs from the experiment, add some extra runs, and quickly re-run the analysis, without the need to redo the analysis of the runs already processed. All this is enabled by the **Use existing .quant files when available** option. The .quant files are saved to/read from the **Temp/.dia dir** (or the same location as the raw files, if there is no temp folder specified). When using this option, the user must ensure that the .quant files had been generated with the exact same settings as applied in the current analysis, with the exception of **Precursor FDR** (provided it is <= 5%), **Threads**, **Log level**, **MBR**, **Cross-run normalisation** and **Library generation** - these settings can be different. It is actually possible to even transfer .quant files to another computer and reuse them there - without transferring the original raw files. Important: it is strongly recommended to only reuse .quant files when both mass accuracies and the scan window are fixed to some values (non-zero), otherwise DIA-NN will perform optimisation of these yet again using the first run for which a .quant file has not been found. Further, when using MBR or creating a spectral library from DIA data with **Library generation** set to smart or full profiling, .quant files should only be reused if they have been generated in exactly the same order as the current order of raw files, that is with MBR DIA-NN currently cannot combine multiple separate analyses together.

</details>  

**Note:** the main report in .parquet format provides the full output information for any kind of downstream processing. All other output types are there to simplify the analysis when using MS Excel or similar software. The numbers of precursors and proteins reported in different types of output files might appear different due to different filtering used to generate those, please see the descriptions above. All the 'matrices' can be reproduced from the main .parquet report, if generated with precursor FDR set to 5%, using R or Python.

### Library-free search
DIA-NN has a very advanced library-free module, which is, for certain types of experiments, better than using a high quality project-specific spectral library. In general, the following makes library-free search perform better in comparison to spectral libraries (while the opposite favours spectral libraries):  
- high peptide numbers detectable per run;
- heterogeneous data (e.g. cancer tissue samples are quite heterogeneous, while replicate injections of the same sample are not);
- long chromatographic gradients as well as good separation of peptides in the ion mobility dimension;
- large dataset (although processing a large dataset in library-free mode might take time).

Please note that in 99% of cases it is essential that [MBR](#match-between-runs) is enabled for a quantitative library-free analysis. It gets activated by default when using the DIA-NN GUI.

For most experiments it does indeed make sense to try library-free search. For medium and large-scale experiments it might make sense to first try library-free analysis of a subset of the data, to see whether the performance is OK (on the whole dataset it will typically be a lot better, so no need to be too stringent here). Ourselves we also often perform a quick preliminary QC assessment of the experiment using some public library.

It is often convenient to perform library-free analysis in two steps: by first [creating](#creation-of-spectral-libraries) an in silico-predicted spectral library from the sequence database and then analysing with this library. This is the strategy that must be used in all cases except for quick preliminary analyses. Note that the pipeline functionality in DIA-NN allows to easily schedule sequences of tasks, such as creation of a predicted library followed by multiple analyses using this library.

<details>
  <summary>Comment</summary>

Note that the larger the search space (the total number of precursors considered), the more difficult it is for the analysis software to identify peptides, and the more time the search takes. DIA-NN is very good at handling very large search spaces, but even DIA-NN cannot do magic and produce as good results with a 100 million search space, as it would with a 2 million search space. So one needs to be careful about enabling all possible variable modifications at once. For example, allowing max 5 variable modifications, while having methionine oxidation, phospho and deamidation enabled simultaneously, is probably not a good idea.

Here lies an important distinction between DIA and DDA data analysis. In DDA allowing all possible variable modifications makes a lot of sense also because the search engine needs to match the spectrum to something - and if it is not matched to the correct modified peptide, it will be matched falsely. In DIA the approach is fundamentally different: the best-matching spectrum is found in the data for each precursor ion being considered (this is a very simplified view just to illustrate the concept). So not being able to identify a particular spectrum is never a problem in DIA (in fact most spectra are highly multiplexed in DIA - that is originate from multiple peptides - and only a fraction of these can be identified). And therefore it only makes sense to enable a particular variable modification if either you are specifically interested in it or if the modification is really ubiquitous.

See [PTMs and peptidoforms](#ptms-and-peptidoforms) for information on distinguishing between peptidoforms bearing different sets of modifications.
</details>

### Creation of spectral Libraries
DIA-NN can create a spectral library from any DIA dataset. This can be done in both spectral library-based and library-free modes: just select the **Generate spectral library** option in the output pane.

DIA-NN can further create an in silico-predicted spectral library out of either a sequence database (make sure **FASTA digest** is enabled) or another spectral library (often useful for public libraries): just run DIA-NN without specifying any raw files and enable the **Deep learning-based spectra, RTs and IMs prediction** option in the **Precursor ion generation** pane. The modifications currently supported by the deep learning predictor are: C(cam), M(ox), N-term acetyl, N/Q(dea), S/T/Y(phos), K(-GG), nK(mTRAQ) and nK(TMT). Of note, if the predictor module in DIA-NN does not recognise some modification, it will still carry out prediction just ignoring it. To make DIA-NN instead discard any peptides with modifications unknown to the predictor, use --skip-unknown-mods.

Spectral libraries can also be created from DDA data, and in fact offline fractionation + DDA has been the 'gold standard' way of creating libraries since the introduction of SWATH/DIA proteomics. For this we recommend using [FragPipe](https://fragpipe.nesvilab.org/), which is based on the ultra-fast and highly robust MSFragger search engine. FragPipe can further be used to create DIA-NN-compatible libraries also from DIA data, similar to DIA-NN itself.

### Match-between-runs
**MBR** is a powerful mode in DIA-NN, which is beneficial for most quantitative experiments, both with a spectral library and in library-free mode. MBR typically results in both higher average ID numbers, but also a lot better data completeness, that is a lot less missing values.

While processing any dataset, DIA-NN gathers a lot of useful information which could have been used to process the data better. And that is what is enabled by MBR. With MBR, DIA-NN first creates a spectral library from DIA data, and then re-processes the same dataset with this spectral library. The algorithmic innovation implemented in DIA-NN ensures that the FDR is stringently controlled: MBR has been validated on datasets ranging from 2 runs to over 1000 runs.

MBR should be enabled for any quantitative experiment, unless you have a very high quality project-specific spectral library, which you think (i) is likely to provide almost complete coverage of detectable peptides, that is there is no point in trying library-free search + MBR, and (ii) most of the peptides in the library are actually detectable in the DIA experiment. If only (i) is true, might be worth still trying MBR along with **Library generation** set to **IDs profiling**.

MBR should not be used for non-quantitative experiments, that is when you only want to create a spectral library, which you would then use on some other dataset.

One can manually 'imitate' MBR using a two-step approach that will result in comparable performance. First, run DIA-NN to create a spectral library from the DIA runs (the whole experiment or just its subset, which can be a lot quicker for large-scale experiments or experiments including blanks/failed runs). Then use this library to analyse the whole experiment. In either case, run DIA-NN with MBR disabled.

When using MBR (or its imitation) and relying on the main .parquet report (recommended) instead of the quantitative matrices, use the following q-value filters:
- Lib.Q.Value instead of Global.Q.Value
- When applying a filter to Q.Value that is more stringent than the FDR threshold used to generate the DIA library (e.g. Q.Value < 0.001 filter), always apply the same filter to Lib.Q.Value
- Lib.PG.Q.Value instead of Global.PG.Q.Value
- Lib.Peptidoform.Q.Value instead of Global.Peptidoform.Q.Value, when using peptidoform scoring

### Changing default settings
DIA-NN can be successfully used to process almost any experiment with default settings. In general, it is recommended to only change settings when specifically advised to do so in this Documentation (like below), for a specific experiment type, or if there is a very clear and compelling rationale for the change.

In many cases, one might want to change several parameters in the **Algorithm** pane.
- **MBR** should be enabled in most cases, see [Match-between-runs](#match-between-runs).
- **Mass accuracies**: when set to 0.0, DIA-NN determines mass tolerances automatically, based either on the first run in the experiment (default), or, if **Unrelated runs** option is selected, for each run separately. However, the automatic algorithm can be affected by the noise in the data, so even for replicate injections, say, acquired on TripleTOF 6600, it can easily yield recommended MS2 mass accuracy tolerances in the range 15ppm - 25ppm - this is perfectly OK. So what we prefer to do in most cases is run DIA-NN on several acquisitions from the experiment, with any spectral library (can choose some small one which allows for quick analysis), see what mass accuracies DIA-NN sets automatically (it prints its recommendations), and set the values to approximate averages of these. Also, often it is known already what DIA-NN parameters are optimal for particular LC-MS settings.
- **Scan window**: ideally should correspond to the approximate average number of data points per peak. Similarly to mass accuracies, can be determined by DIA-NN automatically, but we prefer to have it fixed to some average value.                   

Please also see the guidance on [Library-free search](#library-free-search), [PTMs and peptidoforms](#ptms-and-peptidoforms) and [Multiplexing using plexDIA](#multiplexing-using-plexdia), if these are relevant for your experiment.                                      

Note that once you select a particular option in the DIA-NN GUI, some other settings might get activated automatically. For example, whenever you choose to perform an in silico FASTA database digest (for library-free search), or just generate a spectral library from DIA data, [MBR](#match-between-runs) will get automatically selected too - because in 99% of cases it is beneficial.

### Command-line tool
DIA-NN is implemented as a graphical user interface (GUI), which invokes a command-line tool (diann.exe). The command-line tool can also be used separately, e.g. as part of custom automated processing pipelines. Further, even when using the GUI, one can pass options/commands to the command-line tool, in the **Additional options** text box. Some of such useful options are mentioned in this Documentation, and the full reference is provided in [Command-line reference](#command-line-reference).

When the GUI launches the command-line tool, it prints in the log window the exact set of commands it used. So in order to reproduce the behaviour observed when using the GUI (e.g. if you want to do the analysis on a Linux cluster), one can just pass exactly the same commands to the command-line tool directly.
```
diann.exe [commands]  
```
Commands are processed in the order they are supplied, and with most commands this order can be arbitrary.

For convenience, as well as for handling experiments consisting of thousands of files, some of the options/commands can be stored in a config file. For this, create a text file with any extension, say, diann_config.cfg, type in any commands supported by DIA-NN in there, and then reference this file with --cfg diann_config.cfg (in the **Additional options** text box or in the command used to invoke the diann.exe command-line tool).

### Visualisation
DIA-NN provides two visualisation options.

**Skyline**. To visualise chromatograms/spectra in Skyline, analyse your experiment with MBR and a FASTA database specified and then click the 'Skyline' button. DIA-NN will automatically launch Skyline (make sure you have Skyline/Skyline daily version 23.1.1.459 or later installed). Currently this workflow does not support multiplexing and will not work with modifications in any format other than UniMod.

**DIA-NN Viewer**. Analyse your experiment with the "XICs" checkbox checked and click the 'Viewer' button. By default "XICs" option will make DIA-NN extract chromatograms for the library fragment ions only and within 10s from the elution apex. Use --xic [N] to set the retention time window to N seconds (e.g. --xic 60 will extract chromatograms within a minute from the apex) and --xic-theoretical-fr to extract all charge 1 and 2 y/b-series fragments, including those with common neutral losses. Note that using --xic-theoretical-fr, especially in combination with large retention time window, might require a significant amount of disk space in the output folder. However the visualisation itself is effectively instantaneous, for any experiment size.

<img src="https://raw.githubusercontent.com/vdemichev/DiaNN/master/GUI/Viewer.png"></br>

**Note**: The chromatograms extracted with "XICs" are saved in Apache .parquet format (file names end with '.xic.parquet') and can be readily accessed using R or Python. This can be sometimes convenient to prepare publication-ready figures (although can do that with Skyline or DIA-NN Viewer too), or even to set up automatic custom quality control for LC-MS performance.

Peptide & modification positions within a protein can be visualised using AlphaMap by the Mann lab https://github.com/MannLabs/alphamap.

### Automated pipelines
The pipeline window within the DIA-NN GUI allows to combine multiple analysis steps into pipelines. Each pipeline step is a set of settings as displayed by the GUI. One can add such steps to the pipeline, update existing steps, remove steps, move steps up/down in the pipeline, disable/enable (by double mouse-click) certain steps within the pipeline and save/load pipelines. Further, individual pipeline steps can be copy-pasted between different GUI tabs/windows (use Copy and Paste buttons for this). We always assemble all DIA-NN runs for a particular publication in a pipeline. One can also use DIA-NN pipelines to store configuration templates.

### PTMs and peptidoforms

DIA-NN GUI features built-in workflows (**Precursor ion generation** pane) for detecting methionine oxidation, N-terminal protein acetylation, phosphorylation and ubiquitination (via the detection of remnant -GG adducts on lysines). Other modificaitons can be declared using --var-mod or --fixed-mod in **Additional options**.

Distinguishing between peptidoforms bearing different sets of modifications is a non-trivial problem in DIA: without special peptidoform scoring the effective peptidoform FDR can be in the range 5-10% for library-free analyses. DIA-NN implements a statistical target-decoy approach for peptidoform scoring, which is enabled by the **Peptidoforms** option (**Algorithm** pane) and is also activated automatically whenever a variable modification is declared, via the GUI settings or the --var-mod command. The resulting peptidoform q-values reflect DIA-NN's confidence in the correctness of the set of modifications reported for the peptide as well as the correctness of the amino acid sequence identified. These q-values, however, do not guarantee the absence of low mass shifts due to some amino acid substitutions or modifications such as deamidation (note that DDA does not guarantee this either).

Further, DIA-NN features an algorithm which reports PTM localisation confidence estimates (as posterior probabilities for correct localisation of all variable PTM sites on the peptide as well as scores for individual sites), included in the .parquet output report. The phosphosites_90 and phosphosites_99 .tsv files contain phosphosite-specific quantities, calculated using the Top 1 method (experimental), that is the highest intensity among precursors with the site localised with the specified confidence (0.9 or 0.99, respectively) is used as the phosphosite quantity in the given run. The 'top 1' algorithm is used here as it is likely the most robust against outliers and mislocalisation errors. However, whether or not this is indeed the best option needs to be investigated, which is currently challenging due to the lack of benchmarks with known ground truth.

In general, when looking for PTMs, we recommend the following:

* Essential: the variable modifications you are looking for must be specified as variable (via the GUI checkboxes or the **Additional options**) both when generating an in silico predicted library and also when analysing the raw data using any predicted or empirical library

* Settings for phosphorylation: max 3 variable modifications, max 1 missed cleavage, phosphorylation is the only variable modification specified, precursor charge range 2-3; to reduce RAM usage, make sure that the precursor mass range specified (when generating a predicted library) is not wider than the precursor mass range selected for MS/MS by the DIA method; to speed up processing when using a predicted library, first generate a DIA-based library from a subset of experiment runs (e.g. 10+ best runs) and then analyse the whole dataset using this DIA-based library with MBR disabled

* When the above succeeds, also try max 2 missed cleavages

* When looking for PTMs other than phosphorylation, in 95% of cases best to use max 1 to 3 variable modifications and max 1 missed cleavage

* When not looking for PTMs, i.e. when the goal is relative protein quantification, enabling variable modifications typically does not yield higher proteomic depth. While it usually does not hurt either, it will make the processing slower.

To the best of our knowledge, there is no published validation of the identification confidence for the detection of deamidated peptides (which are easy to confuse to heavier isotopologues, unless the mass spec has a very high resolution and a tight mass accuracy/tolerance setting is used by the search engine), even for DDA. One way to gain confidence in the identification of deamidated peptides is to check if anything is identified if the mass delta for deamidation is declared to be 1.022694, instead of the correct value 0.984016. DIA-NN does pass this test successfully on several datasets (that is no IDs are reported when specifying this 'decoy modification mass'), but we do recommend also trying such 'decoy modification mass' search on several runs from the experiment to be analysed, if looking for deamidated peptides. In each case (correct or decoy mass), --ptm-qvalues should be used to enable PTM-specific scoring for deamidation, in addition to the peptidoform scoring, and either PTM.Q.Value or Global.Q.Value/Lib.Q.Value used for filtering.

Of note, when the ultimate goal is the identification of proteins, it is largely irrelevant if a modified peptide is misidentified, by being matched to a spectrum originating from a different peptidoform. Therefore, if the purpose of the experiment is to identify/quantify specific PTMs, amino acid substitutions or distinguish proteins with high sequence identity, then the **Peptidoforms** scoring option is recommended. In all other cases peptidoform scoring is typically OK to use but not necessary, and will usually lead to a somewhat slower processing and a slight decrease in identification numbers when using MBR.

<details>
  <summary>Does DIA-NN need to recognise modifications in the spectral library?</summary>

In general, yes. However, most workflows will work without the need to recognise modifications. Although if unknown modifications are detected in the library, DIA-NN will print a warning listing those, and it is strongly recommended to declare them using --mod. Note that DIA-NN already recognises many common modifications and can also load the whole UniMod database, see the --full-unimod option.
</details>

### Multiplexing using plexDIA
In collaboration with the Slavov laboratory, we have developed plexDIA based on DIA-NN, a technology that allows to benefit from non-isobaric multiplexing (mTRAQ, dimethyl, SILAC) in combination with DIA. To analyse a plexDIA experiment, one needs an in silico predicted or empirical spectral library. DIA-NN then needs to be supplied with the following sets of commands, depending on the analysis scenario.

**Scenario 1**. The library is a regular label-free library (empirical or predicted), and multiplexing is achieved purely with isotopic labelling, i.e. without chemical labelling with tags such as mTRAQ or dimethyl. DIA-NN then needs the following options to be added to **Additional options**:
* --fixed-mod, to declare the base name for the channel labels and the associated amino acids
* --lib-fixed-mod, to in silico apply the modification declared with --fixed-mod to the library
* --channels, to declare the mass shifts for all the channels considered
* --original-mods, to prevent DIA-NN from converting the declared modifications to UniMod

Example for L/H SILAC labels on K and R:
```
--fixed-mod SILAC,0.0,KR,label
--lib-fixed-mod SILAC
--channels SILAC,L,KR,0:0; SILAC,H,KR,8.014199:10.008269
--original-mods
```
Note that in the above SILAC is declared as label, i.e. it is not supposed to change the retention time of the peptide. It is also a zero-mass label here, as it only serves to designate the amino acids that will be labelled. What the combination of --fixed-mod and --lib-fixed-mod does here is simply put (SILAC) after each K or R in the precursor id sequence, in the internal library representation used by DIA-NN. --channels then splits each library entry into two, one with masses 0 (K) and 0 (R) added upon each occurrence of K(SILAC) or R(SILAC) in the sequence, respectively, and another one with  8.014199 (K) and 10.008269 (R).

**Scenario 2**. The library is a regular label-free library (empirical or predicted), and multiplexing is achieved via chemical labelling with mTRAQ.

<u>Scenario 2: Step 1</u>. Label the library in silico with mTRAQ and run the deep learning predictor to adjust spectra/RTs/IMs. For this, run DIA-NN with the input library in the **Spectral library** field, an **Output library** specified, **Deep learning-based spectra, RTs and IMs prediction** enabled, list of raw data files empty and the following options in **Additional options**:
```
--fixed-mod mTRAQ,140.0949630177,nK
--lib-fixed-mod mTRAQ
--channels mTRAQ,0,nK,0:0; mTRAQ,4,nK,4.0070994:4.0070994;mTRAQ,8,nK,8.0141988132:8.0141988132
--original-mods
```
Use the .predicted.speclib file with the name corresponding to the **Output library** as the spectral library for the next step.

<u>Scenario 2: Step 2</u>. Run DIA-NN with the following options:
```
--fixed-mod mTRAQ,140.0949630177,nK
--channels mTRAQ,0,nK,0:0; mTRAQ,4,nK,4.0070994:4.0070994;mTRAQ,8,nK,8.0141988132:8.0141988132
--original-mods
```
Note that --lib-fixed-mod is no longer necessary as the library generated in Step 1 already contains (mTRAQ) at the N-terminus and lysines of each peptide.

**Scenario 3**. The library is a regular label-free library (empirical or predicted), and multiplexing is achieved via chemical labelling with a label other than mTRAQ. The reason this scenario is treated differently from Scenario 2 is that DIA-NN's in silico predictor has not been specifically trained for labels other than mTRAQ, and therefore an extra step to generate predictions is not necessary. Simply run DIA-NN as you would do in Scenario 1, except the --fixed-mod declaration will have a non-zero mass in this case and will not be a label. For example, for 5-channel dimethyl as described by [Thielert et al](https://doi.org/10.15252/msb.202211503):
```
‐‐fixed‐mod Dimethyl, 28.0313, nK
--lib-fixed-mod Dimethyl
‐‐channels Dimethyl,0,nK,0:0; Dimethyl,2,nK,2.0126:2.0126; Dimethyl,4,nK,4.0251:4.0251; Dimethyl,6,nK,6.0377:6.0377; Dimethyl,8,nK,8.0444:8.0444
--original-mods
```

**Scenario 4**. The library is an empirical DIA library generated by DIA-NN from a multiplexed DIA dataset. For example, this could be a library generated by DIA-NN in the first pass of MBR (and you'd like to reuse it to analyse the same or some other runs). The **Additional options** will then be the same as in Scenario 1, Scenario 2: Step 2 or Scenario 3, except (important!) --lib-fixed-mod must not be supplied.

**In all the scenarios above**, an extra option specifying the normalisation strategy must be included in **Additional options**. This can be either --channel-run-norm (pulsed SILAC, protein turnover) or -channel-spec-norm (multiplexing of independent samples).

**Output**. We recommend using the main report in .parquet format for all downstream analyses. Note that PG.Q.Value and GG.Q.Value in the main report are channel-specific, when using multiplexing. The quantities PG.MaxLFQ, Genes.MaxLFQ and Genes.MaxLFQ.Unique are only channel-specific if (i) QuantUMS is used and (ii) either the report corresponds to the second pass of MBR or MBR is not used. Alternatively, one can use the matrices (not recommended), these are precursor-level only. When using matrices, it is essential to specify --matrix-ch-qvalue, with reasonable thresholds 0.01 to 0.5. This setting will not affect the extracted MS1 matrix, which simply reports MS1 signals corresponding to each channel, whenever a precursor is identified in any of the channels - using this matrix is normally not recommended.

### GUI settings reference
<details>
  <summary>Description of selected options</summary>

**Input** pane  
* **Convert to .dia** Convert the selected raw files to DIA-NN's .dia format, for faster subsequent processing, and save them either to the same folder as the respective source raw files or to **Temp/.dia dir** (**Output** pane), if the latter is specified. Conversion is recommended for Sciex files, typically makes little difference for Thermo files and is not recommended for Bruker files.
* **Reannotate** option allows to reannotate the spectral library with protein information from the FASTA database, using the specified digest specificity
* **Contaminants** Adds common contaminants from the Cambridge Centre for Proteomics (CCP) database and automatically excludes them from quantification, see the description of the --cont-quant-exclude option. This option applies when generating a predicted spectral library from a FASTA database or analysing using such a library, if it was generated with **Contaminants** enabled.

**Precursor ion generation** pane  
* **FASTA digest** instructs DIA-NN to in silico digest the sequence database, for library-free search or to generate a spectral library in silico
* **Deep learning-based spectra, RTs and IMs prediction** instructs DIA-NN to perform deep learning-based prediction of spectra, retention times and ion mobility values. This allows not only to make in silico spectral libraries from sequence databases, but also to replace spectra/RTs/IMs in existing libraries with predicted values

**Output** pane
* **Use existing .quant files when available** reuse IDs/quantification information from a previous analysis, see [Output](#output)
* **Temp/.dia dir** specify where .quant files or converted .dia files will be saved, see [Output](#output)

**Algorithm** pane
* **Mass Accuracy** set the MS2 mass tolerance (in ppm), see [Changing default settings](#changing-default-settings)
* **Mass Accuracy MS1** set the MS1 mass tolerance (in ppm), see [Changing default settings](#changing-default-settings)
* **Scan window** sets the scan window radius to a specific value. Ideally, should be approximately equal to the average number of data points per peak, see [Changing default settings](#changing-default-settings)
* **Unrelated runs** determine mass accuracies and scan window, if automatic, independently for different runs, see [Changing default settings](#changing default settings
* **Peptidoforms** activates peptidoform confidence scoring, see [PTMs and peptidoforms](#ptms-and-peptidoforms)
* **MBR** enables MBR, should be enabled for most quantitative experiments, see [MBR](#match-between-runs)
* **No shared spectra** whether to use a spectrum centric-like algorithm to remove interfering precursors. This algorithm is particularly important when considering variable modifications and should always be enabled
* **Neural network classifier** here 'single-pass' mode is the default option and is recommended. The 'double-pass' mode might be better in some scenarios, but it is almost twice slower and it might make reported FDR values slightly less conservative. Double-pass mode must be tested against single-pass on the specific dataset, before a decision is made to use it.
* **Protein inference** this setting primarily affects proteotypicity definition, the default "Genes" is recommended for almost all applications, provided the gene-level information is actually present in the database (non-UniProt databases might lack it). When set to "Off", protein groups from the spectral library are used - this makes sense if protein inference has already been performed during library generation
* **Quantification strategy** QuantUMS (high-precision) is recommended for most scenarious, use QuantUMS (high-accuracy) for experiments where elimination of any ratio compression bias is critical
* **Cross-run normalisation** whether to use global, RT-dependent (recommended) or also signal-dependent (experimental, be very careful about it) cross-run normalisation. Normalisation can also be disabled completely using --no-norm
* **Library generation** this setting determines if and how empirical RTs/IMs and spectra are added to the newly generated library, instead of the theoretical values. IDs, RT & IM profiling is strongly recommended for almost all workflows. When analysing with a high-quality project-specific library, can switch to IDs profiling. Full profiling means always using empirical information, and might only be beneficial (in very rare cases) when having less than ~1000 peptides identified per run, and only if downstream processing is not very sensitive to a bit higher FDR.
* **Speed and RAM usage** this setting is primarily useful for library-free analyses. The first three modes will typically have little difference in terms of ID numbers, while the **Ultra-fast** mode is rather extreme: about 5x faster, but ID numbers are not as good, and the effective FDR might be somewhat higher. The setting affects only the first pass when using MBR
</details>

### Command-line reference
<details>
  <summary>Description of available options/commands</summary>

Note that some options below are strongly detrimental to performance and are only there for benchmarking purposes. So the recommendation is to only use the options which are expected to be beneficial for a particular experiment (e.g. those recommended in the present documentation) based on some clear rationale.

* **--cfg [file name]** specifies a file to load options/commands from
* **--channel-run-norm** normalisation of multiplexed samples will be performed in run-specific manner, i.e. to perform normalisation, for each precursor ion DIA-NN will sum the respective channels within each run and will normalise these sums across runs: use e.g. for protein turnover SILAC experiments
* **--channel-spec-norm** normalisation of multiplexed samples will be performed in channel-specific manner, i.e. each channel in each run is treated as a separate sample to be normalised: use to analyse experiments wherein multiplexing of independent samples is used to boost throughput
* **--channels [channel 1]; [channel 2]; ...** lists multiplexing channels, wherein each channel declaration has the form [channel] = [label group],[channel name],[sites],[mass1:mass2:...], wherein [sites] has the same syntax as for --var-mod and if N sites are listed, N masses are listed at the end of the channel declaration. The spectral library will be automatically split into multiple channels, for precursors bearing the [label group] modification. To add the latter to a label-free spectral library, can use --lib-fixed-mod, e.g. --fixed-mod SILAC,0.0,KR,label --lib-fixed-mod SILAC. See [Multiplexing using plexDIA](#multiplexing-and-plexdia) for usage examples
* **--clear-mods** makes DIA-NN 'forget' all built-in modification (PTM) names
* **--compact-report** instructs DIA-NN to provide less information in the main .tsv report
* **--cont-quant-exclude [tag]** peptides corresponding to protein sequence ids tagged with the specified tag will be excluded from normalisation as well as quantification of protein groups that do not include proteins with the tag
* **--convert** makes DIA-NN convert the mass spec files to the .dia format. The files are either saved to the same location as the input files, or in the Temp/.dia dir, if it is specified (in the GUI or using the --temp option)
* **--cut [specificty 1],[specificity 2],...** specifies cleavage specificity for the in silico digest. Cleavage sites (pairs of amino acids) are listed separated by commas, '\*' indicates any amino acid, and '!' indicates that the respective site will not be cleaved. Examples: "--cut K\*,R\*,!*P" - canonical tryptic specificity, "--cut " - digest disabled
* **--decoy-channel [channel]** specifies the decoy channel masses, wherein [channel] has the same syntax as for --channels
* **--decoys-preserve-spectrum** informs DIA-NN that decoy peptides in the library are already annotated with 'decoy' spectra
* **--dir [folder]** specifies a folder containing raw files to be processed. All files in the folder must be in .raw, .mzML or .dia format
* **--direct-quant** disable QuantUMS and use legacy DIA-NN quantification algorithms instead, also disables channel-specific protein quantification when analysing multiplexed samples
* **--dl-no-im** when using the deep learning predictor, prediction of ion mobilities will not be performed
* **--dl-no-rt** when using the deep learning predictor, prediction of retention times will not be performed
* **--duplicate-proteins** instructs DIA-NN not to skip entries in the sequence database with duplicate IDs (while by default if several entries have the same protein ID, all but the first entry will be skipped)
* **--exact-fdr** approximate FDR estimation for confident peptides based on parametric modelling will be disabled
* **--export-quant** add fragment quantities and quality information to the .parquet output report
* **--ext [string]** adds a string to the end of each file name (specified with --f)
* **--f [file name]** specifies a run to be analysed, use multiple --f commands to specify multiple runs
* **--fasta [file name]** specifies a sequence database in FASTA format (full support for UniProt proteomes), use multiple --fasta commands to specify multiple databases
* **--fasta-filter [file name]** only consider peptides matching the stripped sequences specified in the text file provided (one sequence per line), when processing a sequence database
* **--fasta-search** instructs DIA-NN to perform an in silico digest of the sequence database
* **--fixed-mod [name],[mass],[sites],[optional: 'label']** - adds the modification name to the list of recognised names and specifies the modification as fixed. Same syntax as for --var-mod.
* **--force-swissprot** only consider SwissProt (i.e. marked with '>sp|') sequences when processing a sequence database
* **--foreign-decoys** informs DIA-NN that any decoys included in the library have been generated by a tool other than this version of DIA-NN
* **--full-unimod** loads the complete UniMod modification database and disables the automatic conversion of modification names to the UniMod format
* **--gen-spec-lib** instructs DIA-NN to generate a spectral library
* **--gen-fr-restriction** annotates the library with fragment exclusion information, based on the runs being analysed (fragments least affected by interferences are selected for quantification, why the rest are excluded)
* **--global-mass-cal** disables RT-dependent mass calibration
* **--global-norm** instructs DIA-NN to use simple global normalisation instead of RT-dependent normalisation
* **--high-acc** QuantUMS settings will be otimised for maximum accuracy, i.e. to minimise any ratio compression quantitative bias
* **--ids-to-names** protein sequence ids will also be used as protein names and genes, any information on actual protein names or genes will be ignored
* **--il-eq** (experimental) when using the 'Reannotate' function, peptides will be matched to proteins while considering isoleucine and leucine equivalent
* **--im-window [x]** fixes IM extraction window to the specific value
* **--im-window-factor [x]** controls the minimum size of the IM extraction window, default is 2.0
* **--individual-mass-acc** mass accuracies, if set to automatic, will be determined independently for different runs
* **--individual-reports** a separate output report will be created for each run
* **--individual-windows** scan window, if set to automatic, will be determined independently for different runs
* **--int-removal 0** disables the removal of interfering precursors
* **--lib [file name]** specifies a spectral library. The use of multiple --lib commands (experimental) allows to load multiple libraries in .tsv format
* **--lib-fixed-mod [name]** in silico applies a modification, previously declared using --fixed-mod, to a spectral library
* **--library-headers [name 1],[name 2],...** specifies column names in the spectral library to be used, in the order described in [Spectral library formats](#spectral-library-formats). Use '*' (without quotes) instead of the column name if a particular column is irrelevant, or if DIA-NN already recognises its name
* **--mass-acc [N]** sets the MS2 mass accuracy to N ppm
* **--mass-acc-cal [N]** sets the mass accuracy used during the calibration phase of the search to N ppm (default is 100 ppm, which is adjusted automatically to lower values based on the data)
* **--mass-acc-ms1 [N]** sets the MS1 mass accuracy to N ppm
* **--matrices** output quantities matrices
* **--matrix-ch-qvalue [x]** sets the 'channel q-value' used to filter the output matrices
* **--matrix-qvalue [x]** sets the q-value used to filter the output matrices
* **--matrix-tr-qvalue [x]** sets the 'translated q-value' used to filter the output matrices
* **--matrix-spec-q** run-specific protein q-value filtering will be used, in addition to the global q-value filtering, when saving protein matrices. The ability to filter based on run-specific protein q-values, which allows to generate highly reliable data, is one of the advantages of DIA-NN
* **--max-fr** specifies the maximum number of fragments per precursors in the spectral library being saved
* **--max-pep-len [N]** sets the maximum precursor length for the in silico library generation or library-free search
* **--max-pr-charge [N]** sets the maximum precursor charge for the in silico library generation or library-free search
* **--mbr-fix-settings** when using the 'Unrelated runs' option in combination with MBR, the same settings will be used to process all runs during the second MBR pass
* **--met-excision** enables protein N-term methionine excision as variable modification for the in silico digest
* **--min-fr** specifies the minimum number of fragments per precursors in the spectral library being saved
* **--min-peak** sets the minimum peak height to consider. Must be 0.01 or greater
* **--min-pep-len [N]** sets the minimum precursor length for the in silico library generation or library-free search
* **--min-pr-charge [N]** sets the minimum precursor charge for the in silico library generation or library-free search
* **--min-pr-mz [N]** sets the minimum precursor m/z for the in silico library generation or library-free search
* **--missed-cleavages [N]** sets the maximum number of missed cleavages
* **--mod [name],[mass],[optional: 'label']** declares a modification name. Examples: "--mod UniMod:5,43.005814", "--mod SILAC-Lys8,8.014199,label"
* **--mod-no-scoring [name]** disable scoring of the modification specified (i.e. disable localisation as well as confidence scoring with --ptm-qvalues)
* **--mod-only** only consider peptides bearing modifications that are being scored, i.e. modifications listed with --var-mod
* **--no-calibration** disables mass calibration
* **--no-cut-after-mod [name]** discard peptides generated via in silico cuts after residues bearing a particular modification
* **--no-decoy-channel** disables the use of a decoy channel for channel q-value calculation
* **--no-fr-selection** the selection of fragments for quantification based on the quality assessment of the respective extracted chromatograms will be disabled
* **--no-im-window** disables IM-windowed search
* **--no-isotopes** do not extract chromatograms for heavy isotopologues
* **--no-lib-filter** the input library will be used 'as is' without discarding fragments that might be harmful for the analysis; use with caution
* **--no-main-report** do not produce the main report
* **--no-maxlfq** disables MaxLFQ for protein quantification
* **--no-norm** disables cross-run normalisation
* **--no-peptidoforms** disables automatic activation of peptidoform scoring when variable modifications are declared
* **--no-prot-inf** disables protein inference (that is protein grouping) - protein groups from the spectral library will be used instead
* **--no-quant-files** instructs DIA-NN not to save .quant files to disk and store them in memory instead
* **--no-rt-window** disables RT-windowed search
* **--no-stats** disables the generation of the stats file
* **--no-swissprot** instruct DIA-NN not to give preference for SwissProt proteins when inferring protein groups
* **--original-mods** disables the automatic conversion of known modifications to the UniMod format names
* **--out [file name]** specifies the name of the main output report. The names of all other report files will be derived from this one
* **--out-lib [file name]** specifies the name of a spectral library to be generated
* **--out-lib-copy** copies the spectral library used into the output folder
* **--out-measured-rt** instructs DIA-NN to save raw empirical retention times in the spectral library being generated, instead of saving RTs aligned to a particular scale
* **--peak-translation** instructs DIA-NN to take advantage of the co-elution of isotopologues, when identifying and quantifying precursors; automatically activated when using --channels
* **--peptidoforms** enables peptidoform confidence scoring
* **--pg-level [N]** controls the protein inference mode, with 0 - isoforms, 1 - protein names (as in UniProt), 2 - genes
* **--predict-n-frag [N]** specifies the maximum number of fragments predicted by the deep learning predictor
* **--predictor** instructs DIA-NN to perform deep learning-based prediction of spectra, retention times and ion mobility values
* **--prefix [string]** adds a string at the beginning of each file name (specified with --f) - convenient when working with automatic scripts for the generation of config files
* **--pr-filter [file name]** specify a file containing a list of precursors (same format as the Precursor.Id column in DIA-NN output), FASTA digest will be filtered to only include these precursors
* **--ptm-qvalues** enable the legacy PTM-specific q-value calculation for the scored modifications, i.e. those declared with --var-mod
* **--quant-acc [X]** sets the precision-accuracy balance for QuantUMS to X, where X must be between 0 and 1
* **--quant-fr [N]** sets the number of top fragment ions among which the fragments that will be used for quantification are chosen for the legacy (pre-QuantUMS) quantification mode. Default value is 6
* **--quick-mass-acc** (experimental) when choosing the MS2 mass accuracy setting automatically, DIA-NN will use a fast heuristical algorithm instead of IDs number optimisation
* **--quant-no-ms1** instructs QuantUMS not to use the recorded MS1 quantities directly
* **--quant-sel-runs [N]** instructs QuantUMS to train its parameters on N automatically chosen runs, to speed up training for large experiments, N here must be 6 or greater
* **--quant-train-runs [N1]:[N2]** instructs QuantUMS to train its parameters on runs with indices in the range N1 to N2 (inclusive), e.g. --quant-train-runs 0:5 will perform training on 6 runs with indices 0 to 5
* **--reanalyse** enables MBR
* **--reannotate** reannotate the spectral library with protein information from the FASTA database, using the specified digest specificity
* **--ref [file name]** (experimental) specify a special (small) spectral library which will be used exclusively for calibration - this function can speed up calibration in library-free searches; this special library must have the same RT and IM scales as the main predicted/experimental library used for the analysis; use with caution
* **--regular-swath** all runs will be analysed as if they were not Scanning SWATH runs
* **--relaxed-prot-inf** instructs DIA-NN to use a very heuristical protein inference algorithm (similar to the one used by FragPipe and many other software tools), wherein DIA-NN aims to make sure that no protein is present simultaneously in multiple protein groups. This mode (i) is recommended for method optimisation & benchmarks, (ii) might be convenient for gene set enrichment analysis and related kinds of downstream processing. However the alternative protein inference strategy of DIA-NN is more reliable for differential expression analyses (this is one of the advantages of DIA-NN). Equivalent to the 'Heuristic protein inference' GUI setting, which is activated by default since DIA-NN 1.8.1
* **--report-lib-info** adds extra library information on the precursor and its fragments to the main output report
* **--restrict-fr** some fragments will not be used for quantification, based on the value in the ExcludeFromAssay spectral library column; marking all fragments of a precursor as excluded will, together with --restrict-fr, suppress the use of this precursor for protein quantification whenever possible
* **--scanning-swath** all runs will be analysed as if they were Scanning SWATH runs
* **--semi** (experimental) when using the 'Reannotate' function, a peptide will be matched to a protein also if it could be obtained with one specific and one non-specific cut (at either of the termini)
* **--skip-unknown-mods** instructs DIA-NN to ignore peptides with modifications that are not supported by the deep learning predictor, when performing the prediction
* **--smart-profiling** enables an intelligent algorithm which determines how to extract spectra, when creating a spectral library from DIA data. The performance has so far always been observed to be comparable to IDs, RT and IM profiling
* **--species-genes** instructs DIA-NN to add the organism identifier to the gene names - useful for distinguishing genes from different species, when analysing mixed samples. Works with UniProt sequence databases.
* **--species-ids** instructs DIA-NN to add the organism identifier to the sequence ids/protein isoform ids - useful for distinguishing protein ids from different species, when analysing mixed samples. Works with UniProt sequence databases.
* **--sptxt-acc [N]** sets the fragment filtering mass accuracy (in ppm) when reading .sptxt/.msp libraries
* **--tag-to-ids [tag]** proteins that have the respective FASTA header start with tag (i.e. the string following '>' starts with tag) will have the tag incorporated in the respective protein sequence ids, names and genes
* **--temp [folder]** specifies the Temp/.dia directory
* **--threads [N]** specifies the number of CPU threads to use
* **--tims-min-int [N]** specifies the minimum peak intensity (an integer ion count) to be considered when loading timsTOF data
* **--tims-ms1-cycle [N]** merge the MS/MS spectra from N consecutive cycles, each cycle defined as the set of MS/MS scans following an MS1 scan, use this option only in combination with --tims-scan
* **--tims-scan** enables support for diagonal scanning methods on timsTOF, e.g. Slice-PASEF, midia-PASEF, Synchro-PASEF
* **--tims-skip-errors** DIA-NN will ignore errors when loading dia-PASEF data
* **--use-quant** use existing .quant files, if available
* **--verbose [N]** sets the level of detail of the log. Reasonable values are in the range 0 - 4
* **--var-mod [name],[mass],[sites],[optional: 'label']** - adds the modification name to the list of recognised names and specifies the modification as variable. [sites] can contain a list of amino acids and 'n' which codes for the N-terminus of the peptide. '\*n' indicates protein N-terminus. Examples: "--var-mod UniMod:21,79.966331,STY" - phosphorylation, "--var-mod UniMod:1,42.010565,*n" - N-terminal protein acetylation. Similar to --mod can be followed by 'label'
* **--var-mods** sets the maximum number of variable modifications
* **--xic [optional: X]** instructs DIA-NN to extract MS1/fragment chromatograms for identified precursors within X seconds from the elution apex, with X set to 10s if not provided; the chromatograms are saved in .parquet files (one per run) located in a folder that is created in the same location as the main output report; equivalent to the 'XICs' option in the GUI
* **--xic-theoretical-fr** makes DIA-NN extract chromatograms for all theoretical charge 1-2 fragment ions, including those with common neutral losses, if --xic is used
* **--window [N]** sets the scan window radius to a specific value. Ideally, should be approximately equal to the average number of MS/MS data points per peak
</details>

### Main output reference
<details>
  <summary>Description of selected columns in the main .parquet report</summary>

* **Protein.Group** - inferred proteins. See the description of the **Protein inference** [GUI setting](#gui-settings-reference) and the --relaxed-prot-inf [option](#command-line-reference).
* **Protein.Ids** - all proteins matched to the precursor in the library or, in case of library-free search, in the sequence database
* **Protein.Names** names (UniProt names) of the proteins in the Protein.Group
* **Genes** gene names corresponding to the proteins in the Protein.Group. Note that protein and gene names will only be present in the output if either (i) **Protein inference** is disabled in DIA-NN and no sequence (FASTA) database is provided - then protein & gene names from the spectral library will be used, or (ii) protein IDs in the Protein.Group can be found in the sequence database provided, or (iii) the **Reannotate** function is used. (ii) and (iii) also require that DIA-NN recognises protein & gene names as recorded in the sequence database (e.g. the UniProt format is supported)
* **Precursor.Channels.Group** channel-oblivious identifier of the precursor  
* **Q.Value** run-specific precursor q-value
* **PEP** run-specific posterior error probability for the precursor
* **Global.Q.Value** global precursor q-value
* **Protein.Q.Value** run-specific q-value for the unique protein/gene, that is protein/gene identified with proteotypic (=specific to it) peptides, not channel-specific
* **PG.Q.Value** run-specific q-value for the protein group, channel-specific
* **PG.PEP** run-specific posterior error probability for the protein group, channel-specific
* **Global.PG.Q.Value** global q-value for the protein group
* **GG.Q.Value** run-specific q-value for the gene group (identified by the value in the Genes column), channel-specific
* **Lib.Q.Value** q-value for the respective library entry, 'global' if the library was created by DIA-NN. In case of [MBR](#match-between-runs), this applies to the library created after the first MBR pass
* **Lib.PG.Q.Value** protein group q-value for the respective library entry, 'global' if the library was created by DIA-NN. In case of [MBR](#match-between-runs), this applies to the library created after the first MBR pass
* **Translated.Q.Value** the q-value a precursor would have had if it was not 'translated' between multiplexed channels
* **Channel.Q.Value** reflects the confidence that the precursor is indeed present in the respective channel
* **Ms1.Profile.Corr** correlation between MS1 and MS2 chromatograms
* **Precursor.Quantity** MS2-based quantity of the precursor (charged peptide), non-normalised
* **Precursor.Normalised** MS2-based quantity of the precursor (charged peptide), normalised, if normalisation is not disabled
* **Ms1.Area** non-normalised MS1 peak area
* **Ms1.Normalised** normalised MS1 peak area
* **Ms1.Apex.Area** MS1 peak integrated only in the vicinity of the apex
* **Normalisation.Factor** normalisation factor applied to the precursor in the specific run
* **Peptidoform.Q.Value** q-value reflecting peptidoform confidence if peptidoform scoring is used
* **Global.Peptidoform.Q.Value** q-value reflecting global peptidoform confidence if peptidoform scoring is used
* **Lib.Peptidoform.Q.Value** peptidoform q-value for the respective library entry. In case of [MBR](#match-between-runs), this applies to the library created after the first MBR pass
* **PTM.Site.Confidence** run-specific confidence in site localisation (=estimated 1 - Site.PEP, where Site.PEP is the posterior error probability for incorrect modification localisation)
* **Lib.PTM.Site.Confidence** the site localisation confidence recorded in the library. In case of [MBR](#match-between-runs), this applies to the library created after the first MBR
* **Quantity.Quality** when using QuantUMS is equal to 1.0 / (1.0 + SD), where SD is the standard deviation of the LC-MS-derived error in relative precursor quantification
* **Empirical.Quality** when using QuantUMS reflects the agreement of relative precursor quantities obtained using different quantitative features (MS1 / fragment ions)
* **Normalisation.Noise** reflects the noisiness of the RT-dependent normalisation factor calculated for the precursor
* **RT.Start** and **RT.Stop** peak boundaries
* **FWHM** estimated peak width at half-maximum; note that the accuracy of such estimates sometimes strongly depends on the DIA cycle time and sample injection amount, i.e. they can only be used to evaluate chromatographic performance in direct comparisons with similar settings, including the scan window; another caveat is that FWHM does not reflect any peak tailing
* **PG.MaxLFQ** QuantUMS or MaxLFQ normalised quantity for the protein group, channel-specific
* **Genes.MaxLFQ** QuantUMS or MaxLFQ normalised quantity for the genes group (corresponding to the 'Genes' column), channel-specific
* **Genes.MaxLFQ.Unique** as above, but calculated solely using proteotypic (gene-specific) precursors, channel-specific
* **PG.MaxLFQ.Quality** when using QuantUMS reflects the quality of PG.MaxLFQ
</details>

### Frequently asked questions
**Q: Why DIA-NN?**  
**A:** Proteomic depth, quantitative precision, reliability and speed. In particular, DIA-NN is transformative for (i) experiments acquired with fast chromatographic gradients and (ii) library-free analyses. Further, DIA-NN is the basis of such technologies as plexDIA or Slice-PASEF.

**Q: I have a regular experiment, which parts of this Documentation do I really need to read?**  
**A:** [Getting Started](#getting-started) is sufficient to start working with DIA-NN. Afterwards can also look at [Changing default settings](#changing-default-settings). If something is not working, check [Raw data formats](#raw-data-formats) and [Spectral library formats](#spectral-library-formats). If you think you might want to analyse your data without a spectral library, check [Library-free search](#library-free-search).

**Q: I am new to DIA proteomics, what papers would you suggest?**  
**A:** [Ludwig et al 2018](https://www.embopress.org/doi/full/10.15252/msb.20178126) is an excellent introduction to DIA. The [DIA-NN paper](https://doi.org/10.1038/s41592-019-0638-x) describes how neural network classifier in DIA-NN enables confident proteomics with fast methods and is an example of how to benchmark DIA software performance. Please note that the field of DIA proteomics is developing very rapidly, and things get outdated very quickly. For example, the latest version of DIA-NN is overwhelmingly better in every aspect than DIA-NN 1.6, which was featured in the initial DIA-NN paper.

**Q: How to reduce memory usage/speed up DIA-NN when analysing in library-free mode?**  
**A:** Try adjusting the **Speed and RAM usage** setting. Can even try the **Ultra-fast** mode there, it does sacrifice identification performance and it might have somewhat higher effective FDR, but it is about ~5x faster. Another option is to reduce the precursor mass range, that is search mass ranges 400-500, 500-600, 600-700, etc, separately - create a spectral library from DIA data separately for each mass range, then merge these libraries (e.g. by supplying multiple --lib commands to DIA-NN) and reanalyse the whole dataset with the merged library (without MBR).

The most important factor in library-free searches is the search space size. So here are some ways to reduce the search space and thus speed up library-free analyses:
- Only search peptides which are likely to be detectable in your experiment. For example, for plasma/serum library-free analyses can consider using peptides from the PeptideAtlas [Human Plasma Non-Glyco 2017-04](https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/buildDetails?atlas_build_id=465) build. Technically, the way to implement this is to make a FASTA database wherein each entry corresponds to a single peptide sequence, in combination with "--cut " to disable the digest. Once an in silico spectral library is generated, can **Reannotate** it with protein information. Alternatively, can convert the predicted library from .speclib to .tsv, using DIA-NN, pull protein information from PeptideAtlas and annotate using R/Python.
- Don't consider charge 1 precursors, that is use "--min-pr-charge 2" or the respective GUI option, when generating an in silico predicted library. This is particularly relevant for dia-PASEF.
- To speed up deep learning-based prediction, limit the precursor mass range to the actual mass range of your runs. This also reduces RAM usage, in particular when using multiplexing.
- When analysing a mulitplexed (plexDIA) experiment in lib-free setting, create an empirical DIA-based library by analysing only a single channel as a fixed modification (i.e. without using --channels). Then use this library with --channels (make sure to indicate correct mass shifts for the channels, with respect to the fixed modification used) and with MBR disabled.
- Do not include blanks or failed runs in your lib-free analysis. That is, if you actually want to process those, create a DIA-based library from 'good' runs (or a subset of those), and then analyse the whole experiment, including blanks, with this newly created library.
- Try --ref to speed up calibration.

**Q: How do I ensure reproducibility of my analyses?**  
**A:** DIA-NN has a strong focus on reproducibility. First, the pipeline functionality (and the ability to save pipelines) allows to conveniently store the experiment settings - with DIA-NN it can never happen that a user forgets what settings were used for a particular analysis. Second, DIA-NN's log contains complete information on the settings used and is automatically saved for each experiment, along with a separate pipeline step - which allows to easily rerun the analysis, even if the user forgot to save the whole pipeline.  

In terms of computational reproducibility, please note the following:
- DIA-NN's analysis using a spectral library will always be the same across different Windows machines.
- In silico prediction of spectra/RTs/IMs using deep learning will yield very slightly different results depending on the instruction set supported by the CPU - that's why DIA-NN prints this instruction set upon being launched. This is a consequence of computational optimisations for modern architectures as implemented in the Pytorch deep learning framework used by DIA-NN.
- DIA-NN's analysis on any Linux machine will be very slightly different from that on Windows, because certain mathematical functions are calculated slightly differently by Linux libraries. It is not guaranteed but very likely that analysis results will be the same across all reasonable Linux machines. So for reproducibility purposes in practice it is enough to just keep track of what kind of Linux distribution with what glibc version is being used.

**Q: What is the purpose of the [diann R](https://github.com/vdemichev/diann-rpackage) package?**
**A:** Ourselves we primarily use it for [MaxLFQ](https://doi.org/10.1074/mcp.M113.031591) protein quantification after we do batch correction at the precursor level (for large experiments). Of note, the recent version of the **[iq](https://cran.r-project.org/web/packages/iq/index.html)** R package actually does a better job at calculating MaxLFQ quantities: it is over an order of magnitude faster. So can also use the **iq** package.

**Q: How do I do batch correction?**    
**A:** Ideally on the precursor-level. Then can calculate protein quantities with MaxLFQ, using corrected precursor quantities. And then, if necessary, do another round of batch correction, now on the protein level. In many cases, just doing protein-level batch correction, directly using protein quantities reported by DIA-NN, will also work quite well. Finally, can just incorporate batch information as covariates into statistical models.

**Q: Do I need imputation?**  
**A:** For many downstream applications - no. Most statistical tests are fine handling missing values. When imputation is necessary (e.g. for certain kinds of machine learning, or when a protein is completely absent in some of the biological conditions), we prefer to perform it on the protein level. Note that many papers that discuss imputation methods for proteomic data benchmark them on DDA data. The fundamental difference between DDA and DIA is that in DIA if a value is missing it is a lot more likely that the respective precursor/protein is indeed quite low-abundant. Because of this, minimal value imputation, for example, performs better for DIA than it does for DDA. An important consideration with imputation is that Gaussian statistical methods (like t-test) should be used with caution on imputed data, if the latter is strongly non-Gaussian.

**Q: Can I implement incremental data processing of large-scale experiments based on DIA-NN?**  
**A:** Yes! Sometimes it is important to have an ability to analyse the data incrementally, e.g. first process one cohort of patients, then in a half a year samples from another cohort arrive, and it would be great to merge the data but keep all protein quantities for the old samples unchanged. This is fully supported by DIA-NN (with the 'legacy' quantification mode), and there are several options on how to implement this.
- Option 1. Process each batch separately, then do batch correction. This is the easiest way and might actually work quite well in some cases.
- Option 2. First, fix mass accuracies and the scan window setting. Second, use a spectral library which specifies which fragment ions should or should not be used for peptide quantification (using an ExlcudeFromAssay column), in combination with --restrict-fr and --no-fr-selection. Can generate such an annotation for the library using --gen-fr-restriction, see [Command-line reference](#command-line-reference).This will minimise the batch effects introduced by the separate processing of mass spec batches.
- Option 3. In addition to 2, disable normalisation in DIA-NN (--no-norm) and disable MaxLFQ (--no-maxlfq). Then all precursor and protein quantities will be calculated completely independently for different samples, and thus there will be no batch effects associated with data processing. One can then perform, separately for each batch, normalisation & "batch correction" with a script in R/Python, utilising the QC samples processed within each batch. For example, can replace each precursor quantity with [(precursor quantity) / (average quantity of this precursor in QC samples in the same batch)] * (average quantity of this precursor in QC samples of batch #1).

**Q: What if I do not need the quatities of the older samples to stay the same and just want to save time on processing?**  
**A:** Also possible. Just have the .quant files saved somewhere, and reanalyse the expanded dataset while reusing the old .quant files, see [Output](#output). The --cfg option, which allows to put some commands in a separate file (e.g. all the --f commands specifying the input mass spec files), might prove to be convenient here.

**Q: When designing an experiment, what number of data points per peak should I aim for?**  
**A:** We typically aim for 3-3.5 data points per peak at FWHM (full width at half maximum), which is about 6-7 data points in total. If you set **Scan window** to 0 in DIA-NN, it will print its estimate for the average "Peak width" - this is actually the estimate for the number of data points at FWHM. For example, for [LFQbench](https://www.nature.com/articles/nbt.3685) data this figure is approximately 3.0-3.1. This is enough to achieve good precision of quantification. In terms of identification performance, methods with longer cycle time and hence less data points per peak typically yield better ID numbers, and thus are preferable when the sole purpose is to generate a spectral library, which would then be used on a different dataset.

**Q: How useful are MS1 scans for DIA?**  
**A:** For what is called 'regular proteomics' - not particularly useful. However, MS1 information is crucially important for distinguishing peptidoforms/PTMs and for multiplexing. Also, MS1 scans are, in most cases, quite inexpensive time-wise, so it is almost always worth having them.

**Q: How to choose the FDR/q-value threshold?**  
**A:** Depends on the downstream analysis. In general, relaxing FDR filtering up to 5% might be reasonable in some cases. But it is only justified if it leads to a significant improvement in identification numbers. One option is to reduce the numbers of missing values by applying 1% global q-value filter and, say, 3% or 5% run-specific q-value filter. Again, only justified if the numbers of missing values are reduced significantly. Numbers of modified peptides detected sometimes strongly depend on the q-value filter threshold chosen, with 2% or 3% yielding significantly more IDs than 1% filtering. In such cases it might be worth it to also relax global q-value filter to 2% or 3%.

When it comes to protein q-value filtering, it is a common practice to only filter using global 1% q-value, to ensure at least 99% of proteins are correctly detected in at least one run in the experiment. This also guarantees low effective run-specific protein FDR, provided the spectral library used contains a lot more precursors and proteins than are actually detectable in the DIA experiment. For example, if the spectral library contains 11k proteins and 500k precursors, but only 3k proteins pass 1% global q-value filter, and only 100k precursors map to these proteins, then if the output is filtered at 1% run-specific precursor q-value and 1% global protein q-value, the effective run-specific protein FDR will also be well controlled at about 1%. However, if almost all the proteins in the library pass 1% global q-value threshold, 1% global protein q-value filtering has no effect. As a consequence, the effective run-specific protein FDR is only controlled indirectly by run-specific precursor q-value, and empirically it usually results in 5-10% effective run-specific protein FDR, given 1% run-specific precursor q-value filtering.

Is having 5% effective run-specific protein FDR a problem? In many cases it is not. Suppose a protein passes 1% run-specific protein q-value filtering in 50% of the runs, and in the other 50% it is quantified but has run-specific q-value > 1%. What are the options? One is to just use those quantities, most of them are likely to be correct. Another is to impute by replacing them via random sampling from some distribution, for example. For most applications the first option is clearly better.

However in some cases one might want to be very confident in protein identifications in specific runs. And here what helps a lot is the ability of DIA-NN to calculate run-specific protein q-values, which is absent from many alternative software options.

There is an important distinction between q-values and posterior error probabilities (PEPs). The q-values allow to select a subset of identifications while controlling the proportion of false ones. The PEPs reflect the probabilities that individual identifications are correct. It is not uncommon for identifications with q-values below 1% to have PEP values above 50%, making them unreliable for quantification. Therefore, if you require extra confidence in specific identifications, it does make sense to apply a PEP filter, respectively at the precursor (PEP) or protein (PG.PEP) level.

**Q: What is normalisation and how does it work?**  
**A:** The purpose of normalisation is to cancel out some of the undesirable effects which can occur during sample preparation and acquisition. For example, pipetting errors can result in different total amounts being injected, and thus signals of different intensity being recorded by the mass spec. The basic global normalisation addresses this. It should in general be used under the following assumption: when only considering 'biological' variability of interest, most of the peptides are not differentialy abundant, or the numbers of upregulated / downregulated (in comparison to experiment-wide averages) peptides in each sample are about the same. If this condition is not satisfied (e.g. when analysing fractions obtained with some separation technique, like SEC, for instance), then normalisation should not be used.  

DIA-NN also supports a more advanced normalisation method, which is enabled by default: RT-dependent normalisation. Sometimes, undesirable technical effects can lead to some samples having either hydrophilic or hydrophobic peptides depleted, relative to other samples. These effects are corrected by RT-dependent normalisation. Potentially it might also correct for spray stability issues.

The choice of the normalisation strategy in DIA-NN determines the way .Normalised and .MaxLFQ quantities (in the main report), as well as values in output matrices are calculated. So if normalisation is disabled, these quantities are actually non-normalised.  

**Q: How to compare the identification performance of two DIA software tools?**  
**A:** This is often tricky because FDR estimates are never precise. For example, recently we have [benchmarked](https://www.nature.com/articles/s41467-022-31492-0) (see Supplementary Figure S3 therein) DIA-NN on data acquired with different chromatographic gradients, and what we observed is that although DIA-NN is always somewhat conservative in its FDR estimates (that is when it reports FDR = 1%, the effective FDR as estimated with a two-species spectral library method is in fact lower than 1%), the accuracy of its estimates varies a bit, depending on the data. When comparing different software tools, the differences in FDR estimation accuracy will often be even larger, due to different algorithms used. Note also that certain kinds of data are very sensitive to the FDR threshold, e.g. the ID number can be, say, 20k peptides at 1% FDR, and go up to 30k at 2%. So if software A says it's 1%, real FDR is also 1%, and it reports 20k, while software B says it's 1%, but real FDR is 2%, and it still reports 20k, these two software tools are actually quite far from each other in terms of real performance.

Because of this, it is important to always independently estimate the FDR using an unbiased benchmarking procedure, which is not affected by the internal FDR estimates of the software tools. One of the ways to do that is use two-species spectral libraries, wherein peptides and proteins from one of the species act as 'decoys' (but the processing software does not know that) and are not expected to be found in the sample. Two important considerations here. First, need to make sure that the peptides from the decoy species are properly annotated as not belonging to the target species, and make sure that they absolutely cannot be present in the sample and that they are not among the common contaminants. Because of this it is not recommended, for example, to use a human library as a decoy when analysing a non-human sample. Second, it is important to apply the pi0 correction when estimating false discovery rates, that is estimate the prior probability of incorrect identification, see an [example here](https://www.nature.com/articles/s41467-022-31492-0). Of note, when estimating precursor-level FDR using this method, it is absolutely essential not to apply any kind of protein q-value filtering.

Finally, an important consideration is to make sure to compare apples to apples. For example, it only makes sense to compare global protein q-value filtering in one software with that in another, and not to compare global protein q-value filtering to run-specific one.

Importantly, there is no 'one true way' to benchmark FDR, and therefore it is OK and expected that FDR estimates using controlled benchmarks can deviate from the internal FDR estimates of the software.

**Q: How to compare the quantification performance of two DIA software tools?**   
**A:** The basic level at which this is often done is comparing the CV values. Or, better, comparing the numbers of precursors quantified with CVs less than a given threshold, e.g. 0.1 or 0.2. The difficulty here is that there are algorithms which can significantly improve CV values, but at the same time make quantification less reliable, because what is being integrated is not only the signal from the target precursor, but also signal from some (unknown) interfering precursor. One example is provided by DIA-NN's high accuracy and high precision modes: the high accuracy mode demonstrates better accuracy, while high precision mode yields lower CV values.

[LFQbench](https://www.nature.com/articles/nbt.3685)-type experiments, featuring two mixtures of different species digests, allow to assess the accuracy of quantification. So what makes sense is to require a quantification method to yield both reasonably good accuracy and reasonably good CV values. However there is currently no known way to determine which quantification method is the best, provided the methods considered are OK in terms of both accuracy and CVs. Of note, using the number of proteins detected as differentially abundant in LFQbench as a readout does not work well: most interfering signals come from the species injected in 1:1 ratio, and therefore "integrating interferences" in this case has little to no negative effect on the ability to detect differential abundance.

**Q: What happens if I process data from some related strains of my favourite organism, using a canonical spectral library?**   
**A:** If you don't have the genomes of those related strains sequenced, what might happen is that a canonical peptide in the library will get matched to a spectrum originating from a mutated peptide. For example, it can happen that PEPTIDE is matched to a spectrum originating from P**N**PTIDE. Will this happen often? Probably not. However if you wish to void this, **Peptidoforms** scoring in DIA-NN is the solution.

**Q: What if I don't know what kind of modifications might be there on my peptides and I want to make sure that the peptidoform is correctly identified?**   
**A:** **Peptidoforms** scoring in DIA-NN solves this.

**Q: Can I install multiple versions of DIA-NN on the same PC?**   
**A:** Yes! And this is one of the advantages of DIA-NN. You can keep the old versions for however long you need them, to enable reproducible reprocessing of old experiments, while also having the latest cutting-edge DIA-NN installed.

**Q: Are DIA-NN version numbers indicative of the performance?**   
**A:** No, they are just chronological. In fact, the performance improvement from 1.7.12 to 1.7.16, for example, is greater than from 1.6 to 1.7.12.

**Q: I have a cool idea, can DIA-NN do [some magic]?**   
**A:** Maybe DIA-NN already supports it (there are some major capabilities that are implemented in DIA-NN but not yet reflected in this Documentation), maybe it is realistic to implement. Please email!

### Support
Please post any questions, feedback, comments or suggestions on the [GitHub Discussion board](https://github.com/vdemichev/DiaNN/discussions) (alternatively can create a GitHub issue) or email Vadim Demichev (email: my first name dot my last name at gmail dot com). General discussions on proteomics & mass specs are also welcome on the [GitHub Discussion board](https://github.com/vdemichev/DiaNN/discussions). Can also contact me via Twitter [@DemichevLab](https://twitter.com/DemichevLab) or LinkedIn. If a query is unanswered for several days, please post a reminder (and **keep posting until answered** :)).

For questions of the type "I ran DIA-NN on some experiment and the results are surprising/strange/not as expected, why is that?", please attach, whenever possible, the full DIA-NN log recorded at Log level = 3 (can be saved using the "Save log" button).

In most cases, the reasons for any unexpected performance can be identified with the following steps:
* Check if DIA-NN has printed any errors or warnings (in the log window of the DIA-NN GUI) when it had finished analysing - if yes, read those and assume that these are the reason for the unexpected behaviour you observe, try to fix them (DIA-NN will often suggest how) and see if it works.

* Check if what DIA-NN has printed immediately after being launched makes sense, i.e. could it be a typo or some other mistake which led to DIA-NN being configured incorrectly?

* If you have changed the default settings, was this done following the recommendations in this documentation? If not, try to only change settings when recommended. For example, if you've searched M(ox) as a variable modification, then, unless you are specifically interested in biological signficance of methionine oxidations introduced in vivo, try analysing without M(ox). If the max variable modifications settings was set to some high number, try reducing it to 1.

* If you see zero or too few identifications from a particular run with a predicted library that you have generated from a sequence database, see what results you can get with an empirical spectral library instead. Can use any suitable public library or a DIA-based library generated by DIA-NN from any other runs (could be from a different platform, e.g. some public data). If still not working, try to see if it's a mass calibration problem, by searching with a wide mass window first using --mass-acc-cal 100. Further, make sure that the predicted/empirical library used reflects the background proteome of your samples, not just specific proteins of interest.

* If DIA-NN has exited unexpectedly, could it be that it ran out of memory? Memory usage is expected to be high when (i) the search space is large, e.g. for phospho or metaproteomics searches, or when allowing lots of variable modifications - see the number of library precursors reported by DIA-NN: RAM usage to store the spectral library in memory is approximately 1Gb per 1 million of precursors; (ii) analysing multiplexed data - the library size is then multiplied by [1 + the number of channels]. Try following the steps for reducing RAM usage outlined in [Frequently asked questions (FAQ)](#frequently-asked-questions). Of note, future DIA-NN versions will have significantly reduced memory usage.
