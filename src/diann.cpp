/*
Copyright 2020, Vadim Demichev

This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license,
visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

Uncomment "#define _NO_THERMORAW" in RAWReader.cpp, MSReader.cpp and MSReader.h to enable building without Thermo MSFileReader installed
*/

#define _HAS_ITERATOR_DEBUGGING 0 
#define _ITERATOR_DEBUG_LEVEL 0
#define _CRT_SECURE_NO_WARNINGS
#define EIGEN_MPL2_ONLY

#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4305)

#ifdef _MSC_VER
#define INPUT INPUT_
#include <windows.h>
#undef INPUT
#define noninline __declspec(noinline)
#else
#define noninline __attribute__((noinline))
#endif

// comment if no MKL installation available
// #define CRANIUM_USE_MKL

#ifdef __linux__
#define LINUX
#elif __unix__
#define LINUX
#endif

#define MSTOOLKIT
#define WIFFREADER
#define CPP17

#ifdef LINUX
// #undef MSTOOLKIT
#undef WIFFREADER
#undef CPP17
#if (__GNUC__ >= 7) 
#define CPP17
#endif
#endif

#include "cpu_info.h"
#include "../cranium/src/cranium.h"
#include "../eigen/Eigen/Dense"
#include "../eigen/Eigen/Sparse"
#ifdef CPP17
#include <experimental/filesystem> // requires C++17
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <string>
#include <regex>
#include <algorithm>
#include <random>
#include <atomic>
#include <thread>
#include <math.h>
#include <list>
#include <set>
#include <map>

#ifdef MSTOOLKIT
#include "MSReader.h"
#endif

#ifdef WIFFREADER
char * wiff_load_func = "?diann_wiff_load@@YAPEAXPEAD_NH@Z";
#endif

#define Min(x, y) ((x) < (y) ? (x) : (y))
#define Max(x, y) ((x) > (y) ? (x) : (y))
#define Abs(x) ((x) >= 0 ? (x) : (-(x)))
#define Sgn(x) ((x) >= 0.0 ? (1) : (-1))
#define Sqr(x) ((x) * (x))
#define Cube(x) ((x) * Sqr(x))

DStream dsout;
#define Throw(x) { dsout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << x << "\n"; exit(-1); }
#define Warning(x) { if (Verbose >= 1) dsout << "WARNING: " << x << "\n"; }

typedef std::chrono::high_resolution_clock Clock;
static auto StartTime = Clock::now();
#define Minutes(x) (floor((x) * 0.001 * (1.0/60.0)))
#define Seconds(x) (floor((x) * 0.001) - 60 * Minutes(x))
inline void Time() {
	double time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - StartTime).count();
	int min = Minutes(time);
	int sec = Seconds(time);
	dsout << "[" << min << ":";
	if (sec < 10) dsout << "0";
	dsout << sec << "] ";
}

const double E = 0.000000001;
const double INF = 10000000.0;

int Threads = 1;
int iN = 12;
int nnIter = 11;

bool LCAllScores = false;

bool Standardise = true;
float StandardisationScale = 1.0;

int QuantMode = 0; // 0 - full peak integration, 1 - fixed window
double PeakBoundary = 5.0;
bool NoIfsRemoval = false;
bool NoFragmentSelectionForQuant = false;
bool QuantFitProfiles = false;
bool RestrictFragments = false;

bool BatchMode = true;
int MinBatch = 2000;
int MaxBatches = 10000;
int Batches = 1;

const int MaxLibSize = 1000000;

bool RTProfiling = false;

int ScanFactor = 150;
int WindowRadius = 0;
double ScanScale = 2.2;
bool InferWindow = true;
bool IndividualWindows = false;
bool Calibrate = true;

const int MaxCycleLength = 2000;
double MinPeakHeight = 0.01;

int iRTTargetTopPrecursors = 50;
int iRTTopPrecursors;
int iRTRefTopPrecursors = 250;
double iRTMaxQvalue = 0.1; // only precursors with lower q-value are used for iRT estimation
double MassCalQvalue = 0.1;
double MassCalMs1Corr = 0.90;
int MassCalBinsMax = 1;
int MassCalBins = 1;
int MassCalBinsMs1 = 1;
int MinMassCalBin = 500;
bool MassCalFilter = true;

int RTSegments = 20;
int MinRTPredBin = 20;

bool ForceScanningSWATH = false;
bool ForceNormalSWATH = false;
bool UseQ1 = true;
bool ForceQ1 = false;
bool Q1Cal = false;
bool Q1CalLinear = false;
int MaxQ1Bins = 3;
int MaxIfsSpan = 4;
bool ForceQ1Apex = false;

bool RefCal = false;
int RTWindowIter = 4;
const int RTRefWindowFactor = 10;
const int RTWindowFactor = 40;
const double RTWindowLoss = 0.20;
const double RTWindowMargin = 2.0;
double MinRTWinFactor = 1.0;
int RTWinSearchMinRef = 10;
int RTWinSearchMinCal = 100;

int nnBagging = 12;
int nnEpochs = 1;
int nnHidden = 5;
double nnLearning = 0.003;
double Regularisation = E;
bool nnFilter = true;
bool nnCrossVal = false;

bool TightMassAcc = true;
bool TightMassAauxForCal = false;
double TightMassAccRatioOne = 0.45;
double TightMassAccRatioTwo = 0.2;

bool UseIsotopes = true;
int IDsInterference = 1;
double InterferenceCorrMargin = 3.0;
bool StrictIntRemoval = false;
bool MS2Range = true;
bool GuideClassifier = false;

bool RTWindowedSearch = true;
bool DisableRT = false;

const int CalibrationIter = 4;
bool ForceMassAcc = false;
bool IndividualMassAcc = false;
double CalibrationMassAccuracy = 100.0 / 1000000.0;
double GlobalMassAccuracy = 20.0 / 1000000.0;
double GlobalMassAccuracyMs1 = 20.0 / 1000000.0;
double GeneratorAccuracy = 5.0 / 1000000.0;
double MaxProfilingQvalue = 0.001;
double MaxQuantQvalue = 0.01;
double ProteinQuantQvalue = 0.01;
double ProteinIDQvalue = 0.01;
double PeakApexEvidence = 0.99; // force peak apexes to be close to the detection evidence maxima
double MinCorrScore = 0.5;
double MinMs1Corr = -INF;
double MaxCorrDiff = 2.0;
bool ForceMs1 = false;
bool MS1PeakSelection = true;
int TopN = 1;

bool TranslatePeaks = false; // translate peaks between precursors belonging to the same elution group
bool ExcludeSharedFragments = true; // exclude from quantification fragments that are shared between light and heavy labelled peptides

int MinMassDeltaCal = 20;
int MinMassDeltaCalRef = 8;

int MinCalRec = 100;
int MinCal = 1000;
int MinClassifier = 2000;

double FilterFactor = 1.5;

bool Convert = false;
bool Reannotate = false;
bool UseRTInfo = false; // use .quant files created previously for RT profiling
bool UseQuant = false; // use .quant files created previously; implies UseRTInfo
bool QuantOnly = false; // quantification will be performed anew using identification info from .quant files
bool ReportOnly = false; // generate report from quant files
bool GenRef = false; // update the .ref file with newly computed data
std::string args, lib_file, learn_lib_file, out_file = "report.tsv", prosit_file = "lib.prosit.csv", predictor_file = "lib.predicted.speclib", out_lib_file = "lib.tsv", out_dir = "", out_gene_file = "report.genes.tsv", temp_folder = "";
std::string ref_file, gen_ref_file, all_fastas;
std::vector<std::string> ms_files, fasta_files, fasta_filter_files;
std::set<std::string> failed_files;
bool FastaSearch = false;
bool GenSpecLib = false;
bool ProfileQValue = true;
bool RTLearnLib = false;
bool SaveOriginalLib = false;
bool SaveCalInfo = false;
bool iRTOutputFromLearnLib = true;
bool PredictorSaved = false;

int Verbose = 1;
bool ExportWindows = false;
bool ExportLibrary = false;
bool ExportDecoys = false;
bool GuideLibrary = false;
bool OverwriteLibraryPGs = true;
bool InSilicoRTPrediction = true;
bool ReverseDecoys = false;
bool ForceFragRec = false;
bool ExportRecFragOnly = true;
bool GenFrExclusionInfo = false;

bool ExportProsit = false;

int VisWindowRadius = 8;
std::vector<std::string> Visualise;

enum {
	loss_none, loss_H2O, loss_NH3, loss_CO,
	loss_N, loss_other
};

const double Loss[loss_other + 1] = { 0.0, 18.011113035, 17.026548, 27.994915, 0.0, 0.0 };

int MaxRecCharge = 19;
int MaxRecLoss = loss_N;

double MinMs1RangeOverlap = 0.25;
double QuadrupoleError = 1.0;

std::string CutAfter = "RK";
std::string NoCutBefore = "";

int MissedCleavages = 1;
int MinPeptideLength = 7;
int MaxPeptideLength = 30;
bool NMetExcision = false;
bool AddLosses = false;
int MaxVarMods = 1;
int MinFrAAs = 3;
double MinFrMz = 200.0, MaxFrMz = 1800.0;
double MinPrMz = 300.0, MaxPrMz = 1800.0;
double MinGenFrInt = 0.00001;
double MinRelFrHeight = 0.00001;
double MinGenFrCorr = 0.5;
double MinRareFrCorr = 0.95;
double FrCorrBSeriesPenalty = 0.1;
int MinGenFrNum = 4;
int MinOutFrNum = 4;
int MinSearchFrNum = 3;
const int MaxGenCharge = 4;

double QFilter = 1.0;
double RubbishFilter = 0.5;
int MinRubbishFactor = 5;
int MinNonRubbish = 5000;

bool InferPGs = true;
bool PGsInferred = false;

double ReportQValue = 0.01;
double ReportProteinQValue = 1.0;
int PGLevel = 2, PGLevelSet = PGLevel; // 0 - ids, 1 - names, 2 - genes
std::vector<std::string> ImplicitProteinGrouping = { "isoform IDs", "protein names", "genes" };
bool IndividualReports = false;
bool SwissProtPriority = true;
bool ForceSwissProt = false;
bool FastaProtDuplicates = false;
bool SpeciesGenes = false;

bool ExtendedReport = true;
bool RemoveQuant = false;
bool QuantInMem = false;

#define Q1 true
#define AAS false
#define LOGSC false
#define REPORT_SCORES false
#define EXTERNAL 1
#define ELUTION_PROFILE false
const int ElutionProfileRadius = 10;

enum {
	libPr, libCharge, libPrMz,
	libiRT, libFrMz, libFrI,
	libPID, libPN, libGenes, libPT, libIsDecoy, 
	libFrCharge, libFrType, libFrNumber, libFrLoss,
	libQ, libEG, libFrExc, libCols
};

std::vector<std::string> library_headers = {
	" ModifiedPeptide LabeledSequence FullUniModPeptideName \"FullUniModPeptideName\" modification_sequence ",
	" PrecursorCharge Charge \"PrecursorCharge\" prec_z ",
	" PrecursorMz \"PrecursorMz\" Q1 ",
	" iRT iRT RetentionTime NormalizedRetentionTime Tr_recalibrated \"Tr_recalibrated\" RT_detected ",
	" FragmentMz ProductMz \"ProductMz\" Q3 ",
	" RelativeIntensity RelativeFragmentIntensity RelativeFragmentIonIntensity LibraryIntensity \"LibraryIntensity\" relative_intensity ",
	" UniprotId UniProtIds UniprotID \"UniprotID\" uniprot_id ",
	" Protein Name Protein.Name ProteinName \"ProteinName\" ",
	" Genes Gene Genes \"Genes\" ",
	" IsProteotypic Is.Proteotypic Proteotypic ",
	" Decoy decoy \"Decoy\" \"decoy\" ",
	" FragmentCharge frg_z \"FragmentCharge\" ",
	" FragmentType frg_type \"FragmentType\" ",
	" FragmentNumber frg_nr FragmentSeriesNumber \"FragmentSeriesNumber\" ",
	" FragmentLossType \"FragmentLossType\" ",
	" Q.Value QValue Qvalue qvalue \"Q.Value\" \"QValue\" \"Qvalue\" \"qvalue\" ",
	" ElutionGroup ModifiedSequence ",
	" ExcludeFromAssay ExcludeFromQuantification "
};

const double proton = 1.007825035;
const double OH = 17.003288;
const double C13delta = 1.003355;

std::vector<std::pair<std::string, std::string> > UniMod;
std::vector<std::pair<std::string, int> > UniModIndex;
std::vector<int> UniModIndices;

struct MOD {
	std::string name;
	std::string aas;
	float mass = 0.0;
	int label = 0;

	MOD() {};
	MOD(const std::string &_name, float _mass, int _label = 0) { name = _name, mass = _mass; label = _label; }
	void init(const std::string &_name, const std::string &_aas, float _mass, int _label = 0) { name = _name, aas = _aas, mass = _mass; label = _label;  }
	friend bool operator < (const MOD &left, const MOD &right) { return left.name < right.name || (left.name == right.name && left.mass < right.mass); }
};
std::vector<MOD> FixedMods, VarMods;

std::vector<MOD> Modifications = {
	MOD("UniMod:4", (float)57.021464),
	MOD("Carbamidomethyl (C)", (float)57.021464),
	MOD("Carbamidomethyl", (float)57.021464),
	MOD("CAM", (float)57.021464),
	MOD("+57", (float)57.021464),
	MOD("+57.0", (float)57.021464),
	MOD("UniMod:26", (float)39.994915),
	MOD("PCm", (float)39.994915),
	MOD("UniMod:5", (float)43.005814),
	MOD("Carbamylation (KR)", (float)43.005814),
	MOD("+43", (float)43.005814),
	MOD("+43.0", (float)43.005814),
	MOD("CRM", (float)43.005814),
	MOD("UniMod:7", (float)0.984016),
	MOD("Deamidation (NQ)", (float)0.984016),
	MOD("Deamidation", (float)0.984016),
	MOD("Dea", (float)0.984016),
	MOD("+1", (float)0.984016),
	MOD("+1.0", (float)0.984016),
	MOD("UniMod:35", (float)15.994915),
	MOD("Oxidation (M)", (float)15.994915),
	MOD("Oxidation", (float)15.994915),
	MOD("Oxi", (float)15.994915),
	MOD("+16", (float)15.994915),
	MOD("+16.0", (float)15.994915),
	MOD("Oxi", (float)15.994915),
	MOD("UniMod:1", (float)42.010565),
	MOD("Acetyl (Protein N-term)", (float)42.010565),
	MOD("+42", (float)42.010565),
	MOD("+42.0", (float)42.010565),
	MOD("UniMod:255", (float)28.0313),
	MOD("AAR", (float)28.0313),
	MOD("UniMod:254", (float)26.01565),
	MOD("AAS", (float)26.01565),
	MOD("UniMod:122", (float)27.994915),
	MOD("Frm", (float)27.994915),
	MOD("UniMod:1301", (float)128.094963),
	MOD("+1K", (float)128.094963),
	MOD("UniMod:1288", (float)156.101111),
	MOD("+1R", (float)156.101111),
	MOD("UniMod:27", (float)-18.010565),
	MOD("PGE", (float)-18.010565),
	MOD("UniMod:28", (float)-17.026549),
	MOD("PGQ", (float)-17.026549),
	MOD("UniMod:526", (float)-48.003371),
	MOD("DTM", (float)-48.003371),
	MOD("UniMod:325", (float)31.989829),
	MOD("2Ox", (float)31.989829),
	MOD("UniMod:342", (float)15.010899),
	MOD("Amn", (float)15.010899),
	MOD("UniMod:1290", (float)114.042927),
	MOD("2CM", (float)114.042927),
	MOD("UniMod:359", (float)13.979265),
	MOD("PGP", (float)13.979265),
	MOD("UniMod:30", (float)21.981943),
	MOD("NaX", (float)21.981943),
	MOD("UniMod:401", (float)-2.015650),
	MOD("-2H", (float)-2.015650),
	MOD("UniMod:528", (float)14.999666),
	MOD("MDe", (float)14.999666),
	MOD("UniMod:385", (float)-17.026549),
	MOD("dAm", (float)-17.026549),
	MOD("UniMod:23", (float)-18.010565),
	MOD("Dhy", (float)-18.010565),
	MOD("UniMod:129", (float)125.896648),
	MOD("Iod", (float)125.896648),
	MOD("Phosphorylation (ST)", (float)79.966331),
	MOD("UniMod:21", (float)79.966331),
	MOD("+80", (float)79.966331),
	MOD("+80.0", (float)79.966331),
	MOD("UniMod:259", (float)8.014199, 1),
	MOD("Lys8", (float)8.014199, 1),
	MOD("UniMod:267", (float)10.008269, 1),
	MOD("Arg10", (float)10.008269, 1),
	MOD("UniMod:268", (float)6.013809, 1),
	MOD("UniMod:269", (float)10.027228, 1)
};

inline char to_lower(char c) { return c + 32; }

void init_unimod() {
	std::set<int> indices;
	auto &mods = Modifications;
	std::sort(mods.begin(), mods.end());
	UniMod.resize(mods.size());
	UniModIndex.resize(mods.size());
	for (int i = 0; i < mods.size(); i++) {
		UniMod[i].first = UniModIndex[i].first = mods[i].name;
		if (!std::memcmp(&(mods[i].name[0]), "UniMod", 6)) UniMod[i].second = mods[i].name, indices.insert(UniModIndex[i].second = i);
		else {
			double min = INF;
			int min_index = 0;
			for (int j = 0; j < mods.size(); j++) if (j != i) {
				if (std::memcmp(&(mods[j].name[0]), "UniMod", 6)) continue;
				double delta = Abs(((double)mods[i].mass) - (double)mods[j].mass);
				if (delta < min) min = delta, min_index = j;
			}
			if (min > 0.00001) {
				if (Verbose >= 1) dsout << "Cannot find a UniMod modification match for " << mods[i].name << ": " << min << " minimal mass discrepancy; using the original modificaiton name\n";
				UniMod[i].second = UniMod[i].first;
				indices.insert(UniModIndex[i].second = i);
			} else {
				UniMod[i].second = mods[min_index].name;
				indices.insert(UniModIndex[i].second = min_index);
			}
		}
	}
	UniModIndices.insert(UniModIndices.begin(), indices.begin(), indices.end());
}

std::string to_unimod(const std::string &mod) {
	for (int i = 0; i < UniMod.size(); i++) if (mod == UniMod[i].first) return UniMod[i].second;
	Warning("cannot convert to UniMod: unknown modification");
	return "";
}

int unimod_index(const std::string &mod) {
	for (int i = 0; i < UniModIndex.size(); i++) if (mod == UniModIndex[i].first) return UniModIndex[i].second;
	Warning("cannot convert to UniMod: unknown modification");
	return 0;
}

int unimod_index_number(int index) {
	return std::distance(UniModIndices.begin(), std::lower_bound(UniModIndices.begin(), UniModIndices.end(), index));
}

#ifdef WIFFREADER
std::vector<std::string> MSFormatsExt = { ".dia", ".mzml", ".raw", ".wiff" };
HMODULE wiff_dll = NULL;
bool skip_wiff = false, fast_wiff = false;
#else
std::vector<std::string> MSFormatsExt = { ".dia", ".mzml", ".raw" };
#endif

enum {
	outFile, outPG, outPID, outPNames, outGenes, outPGQ, outPGN, outGQ, outGN, outGQP, outGNP, outModSeq, outStrSeq,
	outPrId, outCharge, outQv, outPQv, outPPt, outPrQ, outPrN, outPrLR,
	outRT, outiRT, outpRT, outpiRT,
	outCols
};

std::vector<std::string> oh = { // output headers
	"File.Name",
	"Protein.Group",
	"Protein.Ids",
	"Protein.Names",
	"Genes",
	"PG.Quantity",
	"PG.Normalised",
	"Gene.Group.Quantity",
	"Gene.Group.Normalised",
	"Gene.Quantity.Unique",
	"Gene.Normalised.Unique",
	"Modified.Sequence",
	"Stripped.Sequence",
	"Precursor.Id",
	"Precursor.Charge",
	"Q.Value",
	"Protein.Q.Value",
	"Proteotypic",
	"Precursor.Quantity",
	"Precursor.Normalised",
	"Label.Ratio",
	"RT",
	"iRT",
	"Predicted.RT",
	"Predicted.iRT",
};

double NormalisationQvalue = 0.01;
double NormalisationPeptidesFraction = 0.4;
bool LocalNormalisation = true;
bool NoRTDepNorm = false;
bool NoSigDepNorm = true;
int LocNormRadius = 250;
double LocNormMax = 2.0;

#define HASH 0
#if (HASH > 0)
unsigned int hash_map[0x10000];
void init_hash() {
	std::mt19937_64 gen(1);
	for (int i = 0; i < 0x10000; i++) hash_map[i] = gen();
}
template<class T> inline unsigned int hashS(T x) { // T mist be 32-bit
	assert(sizeof(T) == 4);
	int y = *((int*)(&x));
	return hash_map[(unsigned int)(y & 0xFFFF)] ^ hash_map[(unsigned int)((y >> 16) & 0xFFFF)];
}
inline unsigned int hashA(int * data, int n) {
	unsigned int res = 0;
	for (int i = 0; i < n; i++) res ^= hash_map[(unsigned int)(data[i] & 0xFFFF)] ^ hash_map[(unsigned int)((data[i] >> 16) & 0xFFFF)];
	return res;
}
unsigned int hashD(float **dataset, int n, int m) {
	unsigned int res = 0;
	for (int i = 0; i < n; i++) res ^= hashA((int*)(dataset[i]), m);
	return res;
}
#endif

bool Predictor = false;
#ifdef PREDICTOR
namespace predictor {
	extern void init_predictor(int id_range, int threads);
	extern bool code_from_precursor(vector<long long> &code, const string &precursor, string &temp_s, vector<string> &temp_sv, const vector<pair<string, int> > &dict);

	class Predictor {
	public:
		Predictor();
		void set_instance(int instance_id);
		map<string, int> get_aa_indices();
		void predict(vector<vector<float> > &spectra, vector<pair<vector<long long>, int> > &codes, int mode, bool verbose = false);
	};
}

predictor::Predictor P;
#endif

class Parameter {
public:
	int min_iter_seek, min_iter_learn;

	Parameter(int _min_iter_seek, int _min_iter_learn) { min_iter_seek = _min_iter_seek, min_iter_learn = _min_iter_learn; }
};

int MaxF = INF, MinF = 0;
const int TopF = 6, auxF = 12;
const int nnS = 2, nnW = (2 * nnS) + 1;
const int qL = 3;
const int QSL[qL] = { 1, 3, 5 };
enum {
	pTimeCorr,
	pLocCorr, pMinCorr, pTotal, pCos, pCosCube, pMs1TimeCorr, pNFCorr, pdRT, pResCorr, pResCorrNorm, pTightCorrOne, pTightCorrTwo, pShadow, pHeavy,
	pMs1TightOne, pMs1TightTwo,
	pMs1Iso, pMs1IsoOne, pMs1IsoTwo,
	pMs1Ratio,
	pBestCorrDelta, pTotCorrSum,
	pMz, pCharge, pLength, pFrNum, 
#if AAS
	pMods,
	pAAs,
	pRT = pAAs + 20,
#else
	pRT,
#endif
	pAcc,
	pRef = pAcc + TopF,
	pSig = pRef + auxF - 1,
	pCorr = pSig + TopF,
	pShadowCorr = pCorr + auxF,
	pShape = pShadowCorr + TopF,
#if Q1
	pQLeft = pShape + nnW,
	pQRight = pQLeft + qL,
	pQPos = pQRight + qL,
	pQNFCorr = pQPos + qL,
	pQCorr = pQNFCorr + qL,
	pN = pQCorr + qL
#else
	pN = pShape + nnW
#endif
};

const int TIC_n = 250;
class RunStats {
public:
	RunStats() {}
	int precursors, proteins;
	double tot_quantity, MS1_tot, MS2_tot, MS1_acc, MS1_acc_corrected, MS2_acc, MS2_acc_corrected, RT_acc, pep_length, pep_charge, pep_mc, fwhm_scans, fwhm_RT;
	int cnt_valid;
	float TIC[TIC_n];
};
bool SaveRunStats = true;

std::vector<std::vector<float*> > training, training_classes;
std::vector<int> nn_mod, nn_shift, nn_scale, nn_cnt;
int nn_mod_n = 4, nn_scale_min = 7, nn_scale_max = 13;
float target_nn[2] = { 1.0, 0.0 };
float decoy_nn[2] = { 0.0, 1.0 };

inline bool train_on_index(int scan, int net) {
	if (!nnCrossVal) return true;
	int mod = nn_mod[net], shift = nn_shift[net], scale = nn_scale[net];
	int index = (scan - shift + scale) / scale;
	return ((index % nn_mod_n) != mod);
}

inline bool test_on_index(int scan, int net) {
	if (!nnCrossVal) return true;
	return !train_on_index(scan, net);
}

std::string location_to_file_name(const std::string& loc) {
	auto result = loc;
	for (int i = 0; i < result.size(); i++) {
		char c = result[i];
		if (c == ':' || c == '/' || c == '\\' || c == '.') result[i] = '_';
	}
	return temp_folder + result;
}

std::string get_extension(const std::string &file) {
	std::string res;
	int pos = file.find_last_of('.');
	if (pos == std::string::npos) return res;
	res = file.substr(pos);
	pos = res.length();
	if (pos > 0) while (res[pos - 1] == '\"') {
		pos--;
		if (!pos) break;
	}
	res = res.substr(0, pos);
	std::transform(res.begin(), res.end(), res.begin(), ::tolower);
	return res;
}

std::string remove_extension(const std::string &file) {
	std::string res;
	int pos = file.find_last_of('.');
	if (pos == std::string::npos) return file;
	res = file.substr(0, pos);
	return res;
}

std::string get_file_name(const std::string &file) {
	std::string res;
	int pos = file.size() - 1;
	while (pos >= 0) {
		if (file[pos] == '\\' || file[pos] == '/') { pos++; break; }
		pos--;
	}
	res = file.substr(pos);
	return res;
}

inline double logp(double x) { return (x >= E) ? log(x) : log(E); }
template <class T> inline double sum(T * x, int n) {
	double ex = 0.0;
	for (int i = 0; i < n; i++) ex += x[i];
	return ex;
}
template <class T> inline double sum(std::vector<T> &x) { return sum(&(x[0]), x.size()); }

template <class T> inline double mean(T * x, int n) {
	assert(n >= 1);
	return sum(x, n) / (double)n;
}
template <class T> inline double mean(std::vector<T> &x) {
	int n = x.size();
	assert(n >= 1);
	return sum(x) / (double)n;
}

template <class T> inline double var(T * x, int n) {
	double ex = mean(x, n), s2 = 0.0;
	for (int i = 0; i < n; i++) s2 += Sqr(x[i] - ex);
	if (n > 1) s2 /= (double)(n - 1);
	return s2;
}

template <class T> inline double var(std::vector<T> &x) { return var(&(x[0]), x.size()); }

template <class Tx, class Ty> inline double corr(Tx * x, Ty * y, int n) {
	assert(n >= 1);

	double sx = 0.0, sy = 0.0, sx2 = 0.0, sy2 = 0.0, sxy = 0.0, dn = (double)n;
	for (int i = 0; i < n; i++) {
		double xi = x[i], yi = y[i];
		sx += xi, sy += yi;
		sx2 += Sqr(xi), sy2 += Sqr(yi);
		sxy += xi * yi;
	}
	return (dn * sxy - sx * sy) / sqrt(Max(E, (dn * sx2 - Sqr(sx)) * (dn * sy2 - Sqr(sy))));
}

template <class T> inline void get_corr_matrix(double * mat, T * x, int p, int imul, int n) {
	double *s = (double*)alloca(p * sizeof(double)), *s2 = (double*)alloca(p * sizeof(double)), dn = (double)n;

	int i, j, k;
	for (i = 0; i < p; i++) {
		s[i] = s2[i] = 0.0;
		auto m_ind = i * p;
		for (j = i + 1; j < p; j++) mat[m_ind + j] = 0.0;
	}

	for (i = 0; i < p; i++) {
		auto v = x + i * imul;
		auto ps = s + i;
		auto ps2 = s2 + i;
		for (k = 0; k < n; k++) *ps += v[k], *ps2 += Sqr(v[k]);
		auto m_ind = i * p, ii = i * imul;
		for (j = i + 1; j < p; j++) {
			auto ji = j * imul;
			for (k = 0; k < n; k++) mat[m_ind + j] += x[ii + k] * x[ji + k];
		}
	}

	for (i = 0; i < p; i++) {
		auto m_ind = i * p;
		double si = s[i], s2i = s2[i], di = dn * s2i - Sqr(si);
		for (j = i + 1; j < p; j++) 
			mat[j * p + i] = mat[m_ind + j] = (dn * mat[m_ind + j] - si * s[j]) / sqrt(Max(E, di * (dn * s2[j] - Sqr(s[j]))));
	}
}

template <class Tx, class Ty> inline double grad_corr(Tx * x, Ty * y, int n) {
	assert(n >= 1);

	double sx2 = 0.0, sy2 = 0.0, sxy = 0.0;
	int j = 1;
	for (int i = 0; i < n - 1; i = j) {
		j = i + 1;
		double xi = x[j] - x[i], yi = y[j] - y[i];
		sx2 += Sqr(xi), sy2 += Sqr(yi);
		sxy += xi * yi;
	}
	return sxy / sqrt(Max(E, sx2 * sy2));
}

template <class T> inline void get_grad_corr_matrix(double * mat, T * x, int p, int imul, int n) {
	double *s2 = (double*)alloca(p * sizeof(double));

	int i, j, k;
	for (i = 0; i < p; i++) {
		s2[i] = 0.0;
		auto m_ind = i * p;
		for (j = i + 1; j < p; j++) mat[m_ind + j] = 0.0;
	}

	for (i = 0; i < p; i++) {
		auto v = x + i * imul;
		auto ps2 = s2 + i;
		for (k = 0; k < n; k++) *ps2 += Sqr(v[k]);
		auto m_ind = i * p, ii = i * imul;
		for (j = i + 1; j < p; j++) {
			auto ji = j * imul;
			for (k = 0; k < n; k++) mat[m_ind + j] += x[ii + k] * x[ji + k];
		}
	}

	for (i = 0; i < p; i++) {
		auto m_ind = i * p;
		double s2i = s2[i];
		for (j = i + 1; j < p; j++)
			mat[j * p + i] = mat[m_ind + j] = mat[m_ind + j] / sqrt(Max(E, s2i * s2[j]));
	}
}

template <class Tx, class Ty> inline double cos(Tx * x, Ty * y, int n) {
	double s2x = 0.0, s2y = 0.0, r = 0.0;
	for (int i = 0; i < n; i++) s2x += Sqr(x[i]), s2y += Sqr(y[i]), r += x[i] * y[i];
	if (s2x < E || s2y < E) return 0.0;
	return r / sqrt(s2x * s2y);
}

template <class Tx, class Ty> inline double scalar(Tx * x, Ty * y, int n) {
	double r = 0.0;
	for (int i = 0; i < n; i++) r += x[i] * y[i];
	return r;
}

template <class Tx, class Ty> inline double scalar(Tx * x, Ty * y, bool * mask, int n) {
	double r = 0.0;
	for (int i = 0; i < n; i++) if (mask[i]) r += x[i] * y[i];
	return r;
}

template <class T> inline double signal_level(T * x, int center, int n) {
	int i;
	float noise = 0.0, signal = 0.0;
	for (i = 0; i < center; i++) {
		float delta = x[i + 1] - x[i], change = Abs(delta);
		if (delta >= 0.0) signal += change;
		else noise += change;
	}
	for (i = center + 1; i < n; i++) {
		float delta = x[i] - x[i - 1], change = Abs(delta);
		if (delta >= 0.0) noise += change;
		else signal += change;
	}
	float total = signal + noise;
	if (total > E) return signal / total;
	else return 0.0;
}

template<class T> inline void smooth(T * dst, T * src, int n) {
	assert(n >= 2);

	dst[0] = (2.0 / 3.0) * src[0] + (1.0 / 3.0) * src[1];
	dst[n - 1] = (2.0 / 3.0) * src[n - 1] + (1.0 / 3.0) * src[n - 2];
	for (int i = 1; i < n - 1; i++)
		dst[i] = 0.5 * src[i] + 0.25 * (src[i - 1] + src[i + 1]);
}

template <class Tw, class Tx> inline double centroid_coo(Tw * weights, Tx * x, int n) {
	int i;
	double u, t;
	for (i = 0, u = t = 0.0; i < n; i++) u += weights[i] * x[i], t += weights[i];
	return u / Max(E, t);
}

template <class T> void solve(T * x, Eigen::MatrixXd &A, Eigen::VectorXd &b, int n) {
	Eigen::VectorXd X = A.fullPivHouseholderQr().solve(b);

	for (int i = 0; i < n; i++) x[i] = X(i);
}

// implementation of the algorithm from
// stat.wikia.com/wiki/Isotonic_regression
std::vector<int> pava_index;
std::vector<double> pava_w, pava_a;
template <class Tx, class Tfunc> std::vector<std::pair<Tx, Tfunc> > PAVA(std::vector<std::pair<Tx, Tfunc> > &data) {
	int i, j, k, l, n = data.size();
	pava_index.resize(n + 1);
	pava_w.resize(n), pava_a.resize(n);

	std::vector<std::pair<Tx, Tfunc> > result(n);
	for (i = 0; i < n; i++) result[i].first = data[i].first;

	pava_a[0] = data[0].second;
	pava_w[0] = 1.0;
	j = 0;
	pava_index[0] = 0;
	pava_index[1] = 1;

	for (i = 1; i < n; i++) {
		j++;
		pava_a[j] = data[i].second;
		pava_w[j] = 1.0;
		while (j > 0 && pava_a[j] < pava_a[j - 1]) {
			pava_a[j - 1] = (pava_w[j] * pava_a[j] + pava_w[j - 1] * pava_a[j - 1])
				/ (pava_w[j] + pava_w[j - 1]);
			pava_w[j - 1] = pava_w[j] + pava_w[j - 1];
			j--;
		}
		pava_index[j + 1] = i + 1;
	}

	for (k = 0; k <= j; k++)
		for (l = pava_index[k]; l < pava_index[k + 1]; l++)
			result[l].second = pava_a[k];

	return result;
}

std::vector<double> rt_p, rt_e, rt_s;
typedef Eigen::Triplet<double> Triplet;
std::vector<Triplet> TrL;

template <class Tx, class Tfunc> void spline(std::vector<double> &coeff, std::vector<double> &points, std::vector<std::pair<Tx, Tfunc> > &data, int segments) {
	int i, k, m, n = data.size();
	auto pava = PAVA(data);
	double w;

	rt_p.resize(n);
	rt_e.resize(n);
	m = segments;
	coeff.resize(2 * m);
	rt_s.resize(m + 1), points.resize(m);
	for (i = 0; i < n; i++) rt_p[i] = data[i].second, rt_e[i] = data[i].first;

	double split = ((double)(n - 1)) / (double)m;
	for (i = 0; i < m; i++) rt_s[i] = rt_e[Min(split * (double)i, n - 1)]; rt_s[m] = rt_e[n - 1];
	for (i = 0; i < m; i++) points[i] = 0.5 * (rt_s[i] + rt_s[i + 1]);
	for (i = 1; i < m; i++) if (points[i] - points[i - 1] < E) {
		for (k = i - 1; k < m - 1; k++) points[k] = points[k + 1];
		m--;
	}

	TrL.clear();
	for (i = k = 0; i < n; i++) {
		double x = rt_e[i];
		if (!k && x <= points[0]) {
			TrL.push_back(Triplet(i, 0, 1.0));
			TrL.push_back(Triplet(i, 1, x - points[0]));
		}
		else {
			for (; k < m && points[k] < x; k++);
			if (k == m) {
				TrL.push_back(Triplet(i, (k - 1) * 2, 1.0));
				TrL.push_back(Triplet(i, (k - 1) * 2 + 1, x - points[k - 1]));
			}
			else {
				w = points[k] - points[k - 1];
				double u = (x - points[k - 1]) / w, v = u - 1.0;
				TrL.push_back(Triplet(i, (k - 1) * 2, (1.0 + 2.0 * u) * Sqr(v)));
				TrL.push_back(Triplet(i, (k - 1) * 2 + 1, w * u * Sqr(v)));

				TrL.push_back(Triplet(i, k * 2, Sqr(u) * (1.0 - 2.0 * v)));
				TrL.push_back(Triplet(i, k * 2 + 1, w * Sqr(u) * v));
			}
		}
	}

	auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(rt_p.data(), rt_p.size());
	Eigen::SparseMatrix<double, Eigen::RowMajor> A(n, coeff.size());
	A.setFromTriplets(TrL.begin(), TrL.end());
	Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
	SQR.compute(A);
	auto X = SQR.solve(B);
	for (i = 0; i < coeff.size(); i++) coeff[i] = X(i);
}

double calc_spline(std::vector<double> &coeff, std::vector<double> &points, double x) {
	double r = 0.0;
	int m = points.size();
	if (!m) return 0.0;
	if (x <= points[0]) {
		r += coeff[0];
		r += coeff[1] * (x - points[0]);
		return r;
	} else if (x >= points[m - 1]) {
		r += coeff[(m - 1) * 2];
		r += coeff[(m - 1) * 2 + 1] * (x - points[m - 1]);
		return r;
	}
	int low = 0, k = m - 1;
	while (k > low + 1) {
		int middle = (low + k) / 2;
		if (x < points[middle]) k = middle;
		else low = middle;
	}
	double w = points[k] - points[k - 1];
	double u = (x - points[k - 1]) / w, v = u - 1.0;
	r += coeff[(k - 1) * 2] * (1.0 + 2.0 * u) * Sqr(v);
	r += coeff[(k - 1) * 2 + 1] * w * u * Sqr(v);

	r += coeff[k * 2] * Sqr(u) * (1.0 - 2.0 * v);
	r += coeff[k * 2 + 1] * w * Sqr(u) * v;
	return r;
}

std::vector<float> rt_diff, rt_diff_sorted;
std::vector<std::pair<float, float> > temp_data;
template <class Tx, class Tfunc> void map_RT(std::vector<double> &coeff, std::vector<double> &points, std::vector<std::pair<Tx, Tfunc> > &data, int segments) {
	spline(coeff, points, data, segments);
	if (data.size() < 3) return;

	rt_diff.clear();
	temp_data.clear();
	rt_diff.reserve(data.size());
	temp_data.reserve(data.size());

	int pos = 0;
	for (auto &p : data) rt_diff.push_back(Abs(p.second - calc_spline(coeff, points, p.first)));

	rt_diff_sorted = rt_diff;
	std::sort(rt_diff_sorted.begin(), rt_diff_sorted.end());
	double max = rt_diff_sorted[0.8 * (double)rt_diff_sorted.size()];
	for (int i = 0; i < data.size(); i++) if (rt_diff[i] <= max + E) temp_data.push_back(data[i]);

	spline(coeff, points, temp_data, segments);
}

template <class Tx, class Ty> inline void flip_map(std::vector<std::pair<Tx, Ty> > &data) {
	int n = data.size();
	for (int i = 0; i < n; i++) data[i] = std::pair<Ty, Tx>(data[i].second, data[i].first);
}

class Group {
public:
	mutable int index, size;
	std::vector<int> e;
	mutable std::vector<int> s;

	Group(std::vector<int> &_e, std::vector<int> &_s) { e = _e, s = _s; }
	inline friend bool operator < (const Group &left, const Group &right) { return left.e < right.e; }
};

std::vector<Group> bipartite_set_cover(std::vector<std::vector<int> > &sets, int N) { // greedy algorithm
	int i, max_size, size, M = sets.size();

	// merge identical sets into groups
	std::set<Group> merged;
	std::vector<Group> groups; // result
	std::vector<int> v(1);
	for (i = 0; i < M; i++) {
		auto &s = sets[i];
		if (!s.size()) continue;
		v[0] = i;
		auto ins = merged.insert(Group(sets[i], v));
		if (!ins.second) ins.first->s.push_back(i);
	}
	if (!merged.size()) return groups;

	// index
	i = max_size = 0;
	int S = merged.size();
	for (auto it = merged.begin(); it != merged.end(); it++) {
		it->index = i++;
		it->size = it->e.size();
		if (it->size > max_size) max_size = it->size;
	}

	// annotate elements
	std::vector<int> removed(N);
	std::vector<std::vector<std::set<Group>::iterator> > elements(N);
	for (auto it = merged.begin(); it != merged.end(); it++) for (auto &e : it->e) elements[e].push_back(it);

	// create size index
	std::vector<int> added(S);
	std::vector<std::list<std::set<Group>::iterator> > size_index(max_size + 1);
	std::vector<std::list<std::set<Group>::iterator>::iterator> index(S);
	for (auto it = merged.begin(); it != merged.end(); it++) {
		auto &s = *it;
		size = s.size;
		size_index[size].push_back(it);
		auto end = size_index[size].end();
		index[it->index] = --end;
	}

	while ((size = size_index.size() - 1) > 0) {
		auto se = size_index.end();
		auto &si = *(--se);
		auto ie = si.end();
		auto &sl = *(--ie); // iterator to one of the groups with the maximum number of elements

		int ind = sl->index;
		added[ind] = 1;
		groups.push_back(*sl);
		si.pop_back();
		auto &e = groups[groups.size() - 1].e;
		for (auto &el : e) if (!removed[el]) {
			removed[el] = 1;
			auto &gr = elements[el];
			for (auto &g : gr) if (!added[g->index]) {
				size = g->size--;
				size_index[size].erase(index[g->index]);
				if (g->size >= 1) size_index[g->size].push_back(g);
				auto end = size_index[g->size].end();
				index[g->index] = --end;
			}
		}
		while (!size_index[size_index.size() - 1].size() && size_index.size()) size_index.resize(size_index.size() - 1);
		if (!size_index.size()) break;
	}

	return groups;
}

class Peak {
public:
    float mz;
    float height;

	Peak() {}
	Peak(float _mz, float _height) {
		mz = _mz;
		height = _height;
	}
    void init(float _mz, float _height) {
        mz = _mz;
        height = _height;
    }
	friend inline bool operator < (const Peak &left, const Peak &right) { return left.mz < right.mz; }

#if (HASH > 0)
	unsigned int hash() { return hashS(mz) ^ hashS(height); }
#endif
};

const int fTypeB = 1 << 0;
const int fTypeY = 1 << 1;
const int fExclude = 1 << 6;

class Product {
public:
	float mz = 0.0;
	float height = 0.0;
	char charge = 0, type = 0, index = 0, loss = 0;

	Product() {}
	Product(float _mz, float _height, int _charge) {
		mz = _mz;
		height = _height;
		charge = _charge;
	}
	Product(float _mz, float _height, int _charge, int _type, int _index, int _loss) {
		mz = _mz;
		height = _height;
		charge = _charge;
		type = _type;
		index = _index;
		loss = _loss;
	}
	void init(float _mz, float _height, int _charge) {
		mz = _mz;
		height = _height;
		charge = _charge;
	}
	friend inline bool operator < (const Product &left, const Product &right) { return left.mz < right.mz; }

	inline int ion_code() { return (((int)type * 20 + (int)charge) * loss_N + (int)loss) * 100 + (int)index + 1; }
#if (HASH > 0)
	unsigned int hash() { return hashS(mz) ^ hashS(height); }
#endif
};

double AA[256];
int AA_index[256];
const char * AAs = "GAVLIFMPWSCTYHKRQEND";
const char * MutateAAto = "LLLVVLLLLTSSSSLLNDQE";

void init_aas() {
	for (int i = 0; i < 256; i++) AA_index[i] = 0, AA[i] = 0.0;

	AA['G'] = 57.021464;
	AA['A'] = 71.037114;
	AA['V'] = 99.068414;
	AA['L'] = 113.084064;
	AA['I'] = 113.084064;

	AA['F'] = 147.068414;
	AA['M'] = 131.040485;
	AA['P'] = 97.052764;
	AA['W'] = 186.079313;
	AA['S'] = 87.032028;

	AA['C'] = 103.009185;
	AA['T'] = 101.047679;
	AA['Y'] = 163.06332;
	AA['H'] = 137.058912;
	AA['K'] = 128.094963;

	AA['R'] = 156.101111;
	AA['Q'] = 128.058578;
	AA['E'] = 129.042593;
	AA['N'] = 114.042927;
	AA['D'] = 115.026943;

	AA['U'] = 150.95363;
	AA['X'] = 196.995499;

	for (int i = 0; i < 20; i++)
		AA_index[AAs[i]] = i;
	AA_index['U'] = AA_index['C'];
	AA_index['X'] = AA_index['M'];
}

std::vector<double> yDelta;
int yDeltaS = 2, yDeltaW = (yDeltaS * 2 + 1);
inline int y_delta_composition_index(char aa) { return AA_index[aa]; }
inline int y_delta_charge_index() { return y_delta_composition_index(AAs[19]) + 1; }
inline int y_delta_index(char aa, int shift) {
	assert(Abs(shift) <= yDeltaS);
	return y_delta_charge_index() + AA_index[aa] * yDeltaW + shift + yDeltaS + 1;
}
int yCTermD = 9, yNTermD = 8;
inline int y_cterm_index(char aa, int shift) {
	int row = (aa == 'K' ? 0 : (aa == 'R' ? 1 : 2));
	return row * yCTermD + Min(shift, yCTermD - 1) + y_delta_index(AAs[19], yDeltaS) + 1;
}
inline int y_nterm_index(int shift) { return Min(shift, yNTermD - 1) + y_cterm_index(AAs[19], yCTermD - 1) + 1; }
inline int y_delta_size() { return y_nterm_index(yNTermD - 1) + 1; }

inline double y_ratio(int i, int charge, const std::string &aas) {
	int n = aas.length(), j;
	double v = yDelta[y_delta_charge_index()] * (double)charge;
	for (j = 0; j < n; j++) v += yDelta[y_delta_composition_index(aas[j])];
	for (j = Max(i - yDeltaS, 0); j <= Min(i + yDeltaS, n - 1); j++)
		v += yDelta[y_delta_index(aas[j], j - i)];
	v += yDelta[y_nterm_index(i - 2)];
	v += yDelta[y_cterm_index(aas[n - 1], n - 1 - i)];

	return v;
}

void y_scores(std::vector<double> &s, int charge, const std::string &aas) {
	int n = aas.length(), i;
	s.resize(n);
	s[0] = s[1] = 0.0;
	for (i = 2; i < n; i++)
		s[i] = s[i - 1] + y_ratio(i, charge, aas);
}

void to_exp(std::vector<double> &s) { for (int i = 0; i < s.size(); i++) s[i] = exp(s[i]); }

std::vector<double> yLoss;
inline int y_loss_aa(char aa, bool NH3) { return AA_index[aa] + 20 * (int)NH3; }
inline int y_loss(bool NH3) { return y_loss_aa(AAs[19], true) + 1 + (int)NH3; }
inline int y_loss_size() { return y_loss(true) + 1; }

void y_loss_scores(std::vector<double> &s, const std::string &aas, bool NH3) {
	int n = aas.length(), i;
	double sum = yLoss[y_loss(NH3)];
	s.resize(n); s[0] = 0.0;
	for (i = n - 1; i >= 1; i--) {
		sum += yLoss[y_loss_aa(aas[i], NH3)];
		s[i] = sum;
	}
}

std::vector<double> InSilicoRT;
int RTTermD = 1, RTTermDScaled = 5;
inline int aa_rt_nterm(char aa, int shift) { return 1 + Min(shift, RTTermD - 1) * 20 + AA_index[aa]; }
inline int aa_pointsterm(char aa, int shift) { return aa_rt_nterm(AAs[19], RTTermD - 1) + 1 + Min(shift, RTTermD - 1) * 20 + AA_index[aa]; }
inline int aa_rt_nterm_scaled(char aa, int shift) { return aa_pointsterm(AAs[19], RTTermD - 1) + 1 + Min(shift, RTTermDScaled - 1) * 20 + AA_index[aa]; }
inline int aa_pointsterm_scaled(char aa, int shift) { return aa_rt_nterm_scaled(AAs[19], RTTermDScaled - 1) + 1 + Min(shift, RTTermDScaled - 1) * 20 + AA_index[aa]; }
inline int mod_rt_index(int mod) { return aa_pointsterm_scaled(AAs[19], RTTermDScaled - 1) + 1 + unimod_index_number(mod); }
inline int mod_rt_index(const std::string &mod) { return mod_rt_index(unimod_index(mod)); }
inline int in_silico_rt_size() { return UniModIndices.size() ? (mod_rt_index(UniModIndices[UniModIndices.size() - 1]) + 1) : (aa_pointsterm_scaled(AAs[19], RTTermDScaled - 1) + 1); }

void init_prediction() {
	yDelta.resize(y_delta_size(), 0.0);
	yLoss.resize(y_loss_size(), 0.0);
	InSilicoRT.resize(in_silico_rt_size(), 0.0);
}

inline int closing_bracket(const std::string &name, char symbol, int pos) {
	int end, par;
	char close = (symbol == '(' ? ')' : ']');
	for (end = pos + 1, par = 1; end < name.size(); end++) {
		char s = name[end];
		if (s == close) {
			par--;
			if (!par) break;
		} else if (s == symbol) par++;
	}
	return end;
}

inline int peptide_length(const std::string &name) {
	int i, n = 0;
	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i = closing_bracket(name, symbol, i);
			continue;
		}
		n++;
	}
	return n;
}

inline int missed_KR_cleavages(const std::string &name) {
	int i, n = 0, last = 0, par, end;
	for (i = 0; i < name.size() - 1; i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i = closing_bracket(name, symbol, i);
			continue;
		}
		if (symbol == 'K' || symbol == 'R') n++, last = 1;
		else last = 0;
	}
	return n - last;
}

std::string get_aas(const std::string &name) {
	int i, end, par;
	std::string result;

	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i = closing_bracket(name, symbol, i);
			continue;
		}
		result.push_back(symbol);
	}
	return result;
}

std::vector<bool> get_mod_state(const std::string &name, int length) {
	int i, j;
	std::vector<bool> result(length, false);

	for (i = j = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			result[j] = true;
			i = closing_bracket(name, symbol, i);
			continue;
		}
		j++;
	}
	return result;
}

std::vector<double> get_sequence(const std::string &name, int * no_cal = NULL) {
	int i, j;
	double add = 0.0;
	std::vector<double> result;
	if (no_cal != NULL) *no_cal = 0;

	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i++;

			int end = closing_bracket(name, symbol, i);
			if (end >= name.size()) Throw(std::string("incorrect peptide name format: ") + name);

			std::string mod = name.substr(i, end - i);
			for (j = 0; j < Modifications.size(); j++) if (Modifications[j].name == mod) {
				if (!result.size()) add += Modifications[j].mass;
				else result[result.size() - 1] += Modifications[j].mass;
				if (no_cal != NULL && Modifications[j].mass < 4.5 && Modifications[j].mass > 0.0 && Abs(Modifications[j].mass - (double)floor(Modifications[j].mass + 0.5)) < 0.1) *no_cal = 1;
				break;
			}
			if (j == Modifications.size())
				Throw(std::string("unknown modification: ") + mod);
			i = end;
			continue;
		}

		result.push_back(AA[symbol] + add);
		add = 0.0;
	}
	return result;
}

std::string to_canonical(const std::string &name) {
	int i, j;
	std::string result;

	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i++;

			int end = closing_bracket(name, symbol, i);
			if (end == std::string::npos) Throw(std::string("incorrect peptide name format: ") + name);

			std::string mod = name.substr(i, end - i);
			for (j = 0; j < UniMod.size(); j++) if (UniMod[j].first == mod) {
				result += std::string("(") + UniMod[j].second + std::string(")");
				break;
			}
			if (j == Modifications.size()) result += std::string("(") + mod + std::string(")"); // no warning here
			i = end;
			continue;
		} else result.push_back(symbol);
	}
	return result;
}
inline std::string to_canonical(const std::string &name, int charge) { return to_canonical(name) + std::to_string(charge); }

std::string to_eg(const std::string &name) {
	int i, j;
	std::string result;

	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i++;

			int end = closing_bracket(name, symbol, i);
			if (end == std::string::npos) Throw(std::string("incorrect peptide name format: ") + name);

			std::string mod = name.substr(i, end - i);
			for (j = 0; j < Modifications.size(); j++) if (Modifications[j].name == mod) {
				if (!Modifications[j].label) result += std::string("(") + UniMod[j].second + std::string(")");
				break;
			}
			if (j == Modifications.size()) result += std::string("(") + mod + std::string(")"); // no warning here
			i = end;
			continue;
		} else result.push_back(symbol);
	}
	return result;
}

std::string to_charged_eg(const std::string &name, int charge) { return to_eg(name) + std::to_string(charge); }

std::string to_prosit(const std::string &name) {
	int i, j;
	std::string result;

	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i++;

			int end = closing_bracket(name, symbol, i);
			if (end == std::string::npos) Throw(std::string("incorrect peptide name format: ") + name);

			if (i + 9 < name.size()) if (!memcmp(&(name[i]), "UniMod:35", 9)) result += "(ox)"; // methionine oxidation
			i = end;
			continue;
		} else result.push_back(symbol);
	}
	return result;
}

std::string pep_name(const std::string &pr_name) {
	int i;
	std::string result = pr_name;
	for (i = result.size() - 1; i >= 0; i--) if (result[i] < '0' || result[i] > '9') break;
	result.resize(i + 1);
	return result;
}

void count_rt_aas(std::vector<double> &count, const std::string &name) {
	int i, pos, l = peptide_length(name);
	double scale = 1.0 / (double)l;

	for (i = pos = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i++;

			int end = closing_bracket(name, symbol, i);
			if (end == std::string::npos) Throw(std::string("incorrect peptide name format: ") + name);
			count[mod_rt_index(name.substr(i, end - i))]++;
			i = end;
			continue;
		}
		count[aa_rt_nterm(symbol, pos)] += 1.0;
		count[aa_pointsterm(symbol, l - 1 - pos)] += 1.0;
		count[aa_rt_nterm_scaled(symbol, pos)] += scale;
		count[aa_pointsterm_scaled(symbol, l - 1 - pos)] += scale;
		pos++;
	}
}

double predict_irt(const std::string &name) {
	int i, pos, l = peptide_length(name);
	double iRT = InSilicoRT[0], scale = 1.0 / (double)l;

	for (i = pos = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i++;

			int end = closing_bracket(name, symbol, i);
			if (end == std::string::npos) Throw(std::string("incorrect peptide name format: ") + name);
			iRT += InSilicoRT[mod_rt_index(name.substr(i, end - i))];
			i = end;
			continue;
		}
		iRT += InSilicoRT[aa_rt_nterm(symbol, pos)];
		iRT += InSilicoRT[aa_pointsterm(symbol, l - 1 - pos)];
		iRT += InSilicoRT[aa_rt_nterm_scaled(symbol, pos)] * scale;
		iRT += InSilicoRT[aa_pointsterm_scaled(symbol, l - 1 - pos)] * scale;
		pos++;
	}
	return iRT;
}

void arguments(int argc, char *argv[]) {
	int i, start, next, end, cleared_mods = false;
	std::string prefix, ext;
	std::vector<std::string> files;

	for (i = 0; i < argc; i++) args += argv[i] + std::string(" ");
	args = std::regex_replace(args, std::regex("\n|\r"), " ");
	if (args.find("--clear-mods") != std::string::npos) Modifications.clear(), cleared_mods = true,
		dsout << "Modification names specified in the spectral library will be used for annotating output. Only suitable for library-based analysis.\n";
	start = args.find("--", 0);
	dsout.s << args << "\n\n";
	while (start != std::string::npos) {
		start += 2;
		next = args.find("--", start);
		end = (next == std::string::npos) ? args.length() : next;

		bool flag = false;
		if (!memcmp(&(args[start]), "cfg", 3)) { // parse the contents of the file
			std::string name = args.substr(start + 4, end - start - 4);
			std::ifstream in(name);
			if (in.is_open()) {
				std::stringstream stream;
				stream << in.rdbuf();
				args = stream.str() + std::string(" ")
					+ args.substr(end, std::string::npos);
				in.close();

				args = std::regex_replace(args, std::regex("\n|\r"), " ");
				start = args.find("--", 0);
				continue;
			} else dsout << "WARNING: failed to open cfg file\n";
		} else if (!memcmp(&(args[start]), "clear-mods ", 11)) {} else if (!memcmp(&(args[start]), "threads ", 8)) Threads = std::stoi(args.substr(start + 8, std::string::npos)), dsout << "Thread number set to " << Threads << "\n";
		else if (!memcmp(&(args[start]), "fasta-search ", 13)) FastaSearch = true, dsout << "Library-free search enabled\n";
		else if (!memcmp(&(args[start]), "pg-level ", 9)) PGLevelSet = PGLevel = std::stoi(args.substr(start + 9, std::string::npos)),
			dsout << "Implicit protein grouping: " << ImplicitProteinGrouping[PGLevel]
			<< "; this determines which peptides are considered 'proteotypic' and thus affects protein FDR calculation\n";
		else if (!memcmp(&(args[start]), "no-batch-mode ", 14)) BatchMode = false, dsout << "Batch mode disabled\n";
		else if (!memcmp(&(args[start]), "verbose ", 8)) Verbose = std::stoi(args.substr(start + 8, std::string::npos));
		else if (!memcmp(&(args[start]), "export-windows ", 15)) ExportWindows = true;
		else if (!memcmp(&(args[start]), "export-library ", 15)) ExportLibrary = true;
		else if (!memcmp(&(args[start]), "gen-fr-restriction ", 19)) GenFrExclusionInfo = true, 
			dsout << "The spectral library will be annotated with information on fragment exclusion from quantification based on the runs analysed\n";
		else if (!memcmp(&(args[start]), "export-decoys ", 14)) ExportDecoys = true;
		else if (!memcmp(&(args[start]), "prosit ", 7)) ExportProsit = true;
		else if (!memcmp(&(args[start]), "vis ", 4)) {
			std::string word;
			std::stringstream list(trim(args.substr(start + 4, end - start - 4)));
			for (i = 0; std::getline(list, word, ','); i++) {
				try {
					if (i == 0) VisWindowRadius = Max(1, std::stoi(word) / 2);
				} catch (std::exception& e) {}
				if (i == 0) continue;
				Visualise.push_back(get_aas(word));
			}
			dsout << "XICs for precursors corresponding to " << Visualise.size() << " peptides will be saved\n";
		} else if (!memcmp(&(args[start]), "cal-info ", 9)) SaveCalInfo = true;
		else if (!memcmp(&(args[start]), "compact-report ", 15)) ExtendedReport = false;
		else if (!memcmp(&(args[start]), "no-isotopes ", 12)) UseIsotopes = false, dsout << "Isotopologue chromatograms will not be used\n";
		else if (!memcmp(&(args[start]), "no-ms2-range ", 13)) MS2Range = false, dsout << "MS2 range inference will not be performed\n";
		else if (!memcmp(&(args[start]), "min-peak ", 9)) MinPeakHeight = std::stod(args.substr(start + 9, std::string::npos)), dsout << "Minimum peak height set to " << MinPeakHeight << "\n";
		else if (!memcmp(&(args[start]), "no-cal-filter ", 14)) MassCalFilter = false,
			dsout << "Peptides with modifications that can cause interferences with isotopologues will not be filtered out for mass calibration\n";
		else if (!memcmp(&(args[start]), "no-nn-filter ", 13)) nnFilter = false,
			dsout << "Peptides with modifications that can cause interferences with isotopologues will be used for neural network training\n";
		else if (!memcmp(&(args[start]), "nn-cross-val ", 13)) nnCrossVal = true,
			dsout << "Neural network cross-validation will be used to tackle potential overfitting\n";
		else if (!memcmp(&(args[start]), "guide-classifier ", 17)) GuideClassifier = true, dsout << "A separate classifier for the guide library will be used\n";
		else if (!memcmp(&(args[start]), "int-removal ", 12)) IDsInterference = std::stoi(args.substr(start + 12, std::string::npos)),
			dsout << "Number of interference removal iterations set to " << IDsInterference << "\n";
		else if (!memcmp(&(args[start]), "int-margin ", 11)) InterferenceCorrMargin = std::stod(args.substr(start + 11, std::string::npos)),
			dsout << "Interference correlation margin set to " << InterferenceCorrMargin << "\n";
		else if (!memcmp(&(args[start]), "strict-int-removal ", 19)) StrictIntRemoval = true, dsout << "Potentially interfering peptides with close (but not the same) elution times will also be discarded\n";
		else if (!memcmp(&(args[start]), "reverse-decoys ", 15)) ReverseDecoys = true, dsout << "Decoys will be generated using the pseudo-reverse method\n";
		else if (!memcmp(&(args[start]), "force-frag-rec ", 15)) ForceFragRec = true, dsout << "Decoys will be generated only for precursors with all library fragments recognised\n";
		else if (!memcmp(&(args[start]), "max-rec-charge ", 15)) MaxRecCharge = Max(2, std::stoi(args.substr(start + 15, std::string::npos))),
			dsout << "Fragment recognition module will consider charges up to " << MaxRecCharge << "\n";
		else if (!memcmp(&(args[start]), "max-rec-loss ", 13)) MaxRecLoss = Max(0, std::stoi(args.substr(start + 13, std::string::npos))),
			dsout << "Fragment recognition module will consider losses with the index up to " << MaxRecLoss << "\n";
		else if (!memcmp(&(args[start]), "gen-spec-lib ", 13)) GenSpecLib = true, dsout << "A spectral library will be generated\n";
		else if (!memcmp(&(args[start]), "lib-gen-direct-q ", 17)) ProfileQValue = false, dsout << "When generating a spectral library, run-specific q-values will be used instead of profile q-values\n";
		else if (!memcmp(&(args[start]), "save-original-lib ", 18)) SaveOriginalLib = true, dsout << "All entries from the library provided will be saved to the newly generated library\n";
#ifdef WIFFREADER
		else if (!memcmp(&(args[start]), "fast-wiff ", 10)) fast_wiff = true, dsout << "Custom fast centroiding will be used when processing .wiff files (WARNING: experimental).\n";
#endif
		else if (!memcmp(&(args[start]), "f ", 2)) files.push_back(trim(args.substr(start + 2, end - start - 2)));
#ifdef CPP17
		else if (!memcmp(&(args[start]), "dir ", 4)) {
			auto dir = trim(args.substr(start + 4, end - start - 4));
			std::replace(dir.begin(), dir.end(), '\\', '/');
			for (auto &file : std::experimental::filesystem::directory_iterator(dir)) {
				auto fname = file.path().string();
				auto extension = get_extension(fname);
				auto ext = std::lower_bound(MSFormatsExt.begin(), MSFormatsExt.end(), extension);
				if (*ext == extension) files.push_back(file.path().string());
			}
		}
#endif
		else if (!memcmp(&(args[start]), "lib ", 4)) lib_file = trim(args.substr(start + 4, end - start - 4));
		else if (!memcmp(&(args[start]), "temp ", 5)) temp_folder = trim(args.substr(start + 5, end - start - 5));
		else if (!memcmp(&(args[start]), "fasta ", 6)) fasta_files.push_back(trim(args.substr(start + 6, end - start - 6)));
		else if (!memcmp(&(args[start]), "fasta-filter ", 13)) fasta_filter_files.push_back(trim(args.substr(start + 13, end - start - 13)));
		else if (!memcmp(&(args[start]), "ref ", 4)) ref_file = trim(args.substr(start + 4, end - start - 4));
		else if (!memcmp(&(args[start]), "out ", 4)) out_file = trim(args.substr(start + 4, end - start - 4));
		else if (!memcmp(&(args[start]), "out-gene ", 9)) out_gene_file = trim(args.substr(start + 9, end - start - 9));
		else if (!memcmp(&(args[start]), "qvalue ", 7)) ReportQValue = std::stod(args.substr(start + 7, std::string::npos)), dsout << "Output will be filtered at " << ReportQValue << " FDR\n";
		else if (!memcmp(&(args[start]), "protein-qvalue ", 15)) ReportProteinQValue = std::stod(args.substr(start + 15, std::string::npos)),
			dsout << "Output will be filtered at " << ReportProteinQValue << " protein-level FDR\n";
		else if (!memcmp(&(args[start]), "no-prot-inf ", 12)) InferPGs = false, dsout << "Protein inference will not be performed\n";
		else if (!memcmp(&(args[start]), "no-swissprot ", 13)) SwissProtPriority = false, dsout << "SwissProt proteins will not be prioritised for protein inference\n";
		else if (!memcmp(&(args[start]), "force-swissprot ", 16)) ForceSwissProt = true, dsout << "Only SwissProt proteins will be considered in library-free search\n";
		else if (!memcmp(&(args[start]), "species-genes ", 14)) SpeciesGenes = true, dsout << "Species suffix will be added to gene annotation; this affects proteotypicity definition\n";
		else if (!memcmp(&(args[start]), "duplicate-proteins ", 19)) FastaProtDuplicates = true, dsout << "Duplicate proteins in FASTA files will not be skipped\n";
		else if (!memcmp(&(args[start]), "out-lib ", 8)) out_lib_file = trim(args.substr(start + 8, end - start - 8));
		else if (!memcmp(&(args[start]), "learn-lib ", 10)) learn_lib_file = trim(args.substr(start + 10, end - start - 10));
		else if (!memcmp(&(args[start]), "predictor ", 10)) {
#ifdef PREDICTOR
			dsout << "Deep learning will be used to generate a new in silico spectral library from peptides provided\n";
			Predictor = true;
#else
			dsout << "WARNING: DIA-NN has been compiled without support for deep learning-based in silico spectral library generation\n";
#endif
		} else if (!memcmp(&(args[start]), "out-measured-rt ", 16)) iRTOutputFromLearnLib = false,
			dsout << "When generating a spectral library without a guide library but with a training library, iRT values (Tr_recalibrated) will correspond to the measured retention times\n";
		else if (!memcmp(&(args[start]), "library-headers ", 16)) {
			std::string word;
			std::stringstream list(trim(args.substr(start + 16, end - start - 16)));
			for (i = 0; std::getline(list, word, ','); i++) {
				if (i >= libCols) {
					dsout << "WARNING: " << word << ": extra headers will be ignored\n";
					break;
				}
				if (word.find("*") == std::string::npos) library_headers[i] = std::string(" ") + word + std::string(" ");
			}
		} else if (!memcmp(&(args[start]), "output-headers ", 15)) {
			std::string word;
			std::stringstream list(trim(args.substr(start + 15, end - start - 15)));
			for (i = 0; std::getline(list, word, ','); i++) {
				if (i >= outCols) {
					dsout << "WARNING: " << word << ": extra headers will be ignored\n";
					break;
				}
				if (word.find("*") == std::string::npos) oh[i] = std::string(" ") + word + std::string(" ");
			}
		} else if (!memcmp(&(args[start]), "mod ", 4)) {
			int label = 0;
			std::string name, mass, l;
			std::stringstream list(trim(args.substr(start + 4, end - start - 4)));
			if (!std::getline(list, name, ',')) dsout << "WARNING: no modification name, modification ignored\n";
			else if (!std::getline(list, mass, ',')) dsout << "WARNING: no modification mass, modification ignored\n";
			else {
				if (std::getline(list, l, ',')) if (l.find("label") != std::string::npos) label = 1;
				Modifications.push_back(MOD(name, (float)std::stod(mass), label));
			}
			dsout << "Modification " << name << " with mass delta " << (float)std::stod(mass) << " added to the list of recognised modifications for spectral library-based search";
			if (label) dsout << " (as label)\n";
			else dsout << "\n";
		} else if (!memcmp(&(args[start]), "fixed-mod ", 10)) {
			int label = 0;
			std::string name, mass, aas, l;
			std::stringstream list(trim(args.substr(start + 10, end - start - 10)));
			if (!std::getline(list, name, ',')) dsout << "WARNING: no modification name, modification ignored\n";
			else if (!std::getline(list, mass, ',')) dsout << "WARNING: no modification mass, modification ignored\n";
			else if (!std::getline(list, aas, ',')) dsout << "WARNING: no amino acids to be modified, modification ignored\n";
			else {
				if (std::getline(list, l, ',')) if (l.size()) if (l.find("label") != std::string::npos) label = 1;
				Modifications.push_back(MOD(name, (float)std::stod(mass), label));
				FixedMods.resize(FixedMods.size() + 1);
				FixedMods[FixedMods.size() - 1].init(name, aas, (float)std::stod(mass));
				dsout << "Modification " << name << " with mass delta " << (float)std::stod(mass) << " at " << aas << " will be considered as fixed\n";
			}
		} else if (!memcmp(&(args[start]), "var-mod ", 8)) {
			int label = 0;
			std::string name, mass, aas, l;
			std::stringstream list(trim(args.substr(start + 8, end - start - 8)));
			if (!std::getline(list, name, ',')) dsout << "WARNING: no modification name, modification ignored\n";
			else if (!std::getline(list, mass, ',')) dsout << "WARNING: no modification mass, modification ignored\n";
			else if (!std::getline(list, aas, ',')) dsout << "WARNING: no amino acids to be modified, modification ignored\n";
			else {
				if (std::getline(list, l, ',')) if (l.size()) if (l.find("label") != std::string::npos) label = 1;
				Modifications.push_back(MOD(name, (float)std::stod(mass), label));
				VarMods.resize(VarMods.size() + 1);
				VarMods[VarMods.size() - 1].init(name, aas, (float)std::stod(mass));
				dsout << "Modification " << name << " with mass delta " << (float)std::stod(mass) << " at " << aas << " will be considered as variable\n";
			}
		} else if (!memcmp(&(args[start]), "ref-cal ", 8)) RefCal = true, dsout << "Reference peptides will be used for calibration\n";
		else if (!memcmp(&(args[start]), "gen-ref ", 8)) GenRef = true, gen_ref_file = trim(args.substr(start + 8, end - start - 8)),
			dsout << "A library of reference peptides will be generated\n";
		else if (!memcmp(&(args[start]), "window ", 7)) {
			WindowRadius = std::stoi(args.substr(start + 7, std::string::npos));
			if (WindowRadius <= 0) dsout << "WARNING: scan window radius should be a positive integer\n";
			else InferWindow = false, dsout << "Scan window radius set to " << WindowRadius << "\n";
		} else if (!memcmp(&(args[start]), "cut-after ", 10)) CutAfter = trim(args.substr(start + 10, end - start - 10)),
			dsout << "In silico digest will include cuts after amino acids: " << CutAfter << "\n";
		else if (!memcmp(&(args[start]), "no-cut-before ", 14)) NoCutBefore = trim(args.substr(start + 14, end - start - 14)),
			dsout << "In silico digest will not include cuts before amino acids: " << NoCutBefore << "\n";
		else if (!memcmp(&(args[start]), "min-pep-len ", 12)) MinPeptideLength = std::stoi(args.substr(start + 12, std::string::npos)),
			dsout << "Min peptide length set to " << MinPeptideLength << "\n";
		else if (!memcmp(&(args[start]), "max-pep-len ", 12)) MaxPeptideLength = std::stoi(args.substr(start + 12, std::string::npos)),
			dsout << "Max peptide length set to " << MaxPeptideLength << "\n";
		else if (!memcmp(&(args[start]), "min-fr-corr ", 12)) MinGenFrCorr = std::stod(args.substr(start + 12, std::string::npos)),
			dsout << "Minimum fragment profile correlation for the inclusion into the spectral library set to " << MinGenFrCorr << "\n";
		else if (!memcmp(&(args[start]), "min-gen-fr ", 11)) MinOutFrNum = std::stoi(args.substr(start + 11, std::string::npos)),
			dsout << "Minimum number of fragments for library generation set to " << MinOutFrNum << "\n";
		else if (!memcmp(&(args[start]), "min-pr-mz ", 10)) MinPrMz = std::stod(args.substr(start + 10, std::string::npos)),
			dsout << "Min precursor m/z set to " << MinPrMz << "\n";
		else if (!memcmp(&(args[start]), "max-pr-mz ", 10)) MaxPrMz = std::stod(args.substr(start + 10, std::string::npos)),
			dsout << "Max precursor m/z set to " << MaxPrMz << "\n";
		else if (!memcmp(&(args[start]), "min-fr-mz ", 10)) MinFrMz = std::stod(args.substr(start + 10, std::string::npos)),
			dsout << "Min fragment m/z set to " << MinFrMz << "\n";
		else if (!memcmp(&(args[start]), "max-fr-mz ", 10)) MaxFrMz = std::stod(args.substr(start + 10, std::string::npos)),
			dsout << "Max fragment m/z set to " << MaxFrMz << "\n";
		else if (!memcmp(&(args[start]), "max-fr ", 7)) MaxF = std::stoi(args.substr(start + 7, std::string::npos)),
			dsout << "Maximum number of fragments set to " << MaxF << "\n";
		else if (!memcmp(&(args[start]), "min-fr ", 7)) MinF = std::stoi(args.substr(start + 7, std::string::npos)),
			dsout << "Minimum number of fragments for library export set to " << MinF << "\n";
		else if (!memcmp(&(args[start]), "min-search-fr ", 14)) MinSearchFrNum = std::stoi(args.substr(start + 14, std::string::npos)),
			dsout << "Minimum number of fragments required for a precursor to be searched set to " << MinSearchFrNum << "\n";
		else if (!memcmp(&(args[start]), "missed-cleavages ", 17)) MissedCleavages = std::stoi(args.substr(start + 17, std::string::npos)),
			dsout << "Maximum number of missed cleavages set to " << MissedCleavages << "\n";
		else if (!memcmp(&(args[start]), "unimod4 ", 8)) {
			int n = FixedMods.size();
			FixedMods.resize(n + 1), FixedMods[n].init(std::string("UniMod:4"), std::string("C"), 57.021464), dsout << "Cysteine carbamidomethylation enabled as a fixed modification\n";
		} else if (!memcmp(&(args[start]), "unimod35 ", 9)) {
			int n = VarMods.size();
			VarMods.resize(n + 1), VarMods[n].init(std::string("UniMod:35"), std::string("M"), 15.994915), dsout << "Methionine oxidation enabled as a variable modification\n";
		} else if (!memcmp(&(args[start]), "var-mods ", 9)) MaxVarMods = std::stoi(args.substr(start + 9, std::string::npos)),
			dsout << "Maximum number of variable modifications set to " << MaxVarMods << "\n";
		else if (!memcmp(&(args[start]), "met-excision ", 6)) NMetExcision = true, dsout << "N-terminal methionine excision enabled\n";
		else if (!memcmp(&(args[start]), "no-rt-window ", 13)) RTWindowedSearch = false, dsout << "Full range of retention times will be considered\n";
		else if (!memcmp(&(args[start]), "disable-rt ", 11)) DisableRT = true, dsout << "All RT-related scores disabled\n";
		else if (!memcmp(&(args[start]), "min-rt-win ", 11)) MinRTWinFactor = std::stod(args.substr(start + 11, std::string::npos)),
			dsout << "Minimum acceptable RT window scale set to " << MinRTWinFactor << "\n";
		else if (!memcmp(&(args[start]), "no-window-inference ", 20)) InferWindow = false, dsout << "Scan window inference disabled\n";
		else if (!memcmp(&(args[start]), "individual-windows ", 19)) IndividualWindows = true, dsout << "Scan windows will be inferred separately for different runs\n";
		else if (!memcmp(&(args[start]), "individual-mass-acc ", 20)) IndividualMassAcc = true, dsout << "Mass accuracy will be determined separately for different runs\n";
		else if (!memcmp(&(args[start]), "individual-reports ", 19)) IndividualReports = true, dsout << "Reports will be generated separately for different runs (in the respective folders)\n";
		else if (!memcmp(&(args[start]), "no-stats ", 9)) SaveRunStats = false;
		else if (!memcmp(&(args[start]), "convert ", 8)) Convert = true, dsout << "MS data files will be converted to .dia format\n";
		else if (!memcmp(&(args[start]), "reannotate ", 11)) Reannotate = true, dsout << "Library precursors will be reannotated using the FASTA database\n";
		else if (!memcmp(&(args[start]), "out-dir ", 8)) out_dir = trim(args.substr(start + 8, end - start - 8));
#ifdef CPP17
		else if (!memcmp(&(args[start]), "remove-quant ", 13)) RemoveQuant = true, dsout << ".quant files will be removed when the analysis is finished\n";
#endif
		else if (!memcmp(&(args[start]), "no-quant-files ", 15)) QuantInMem = true, dsout << ".quant files will not be saved to the disk\n";
		else if (!memcmp(&(args[start]), "use-rt ", 7)) UseRTInfo = true, dsout << "Existing .quant files will be used for RT profiling\n";
		else if (!memcmp(&(args[start]), "use-quant ", 10)) UseQuant = true, dsout << "Existing .quant files will be used\n";
		else if (!memcmp(&(args[start]), "quant-only ", 11)) QuantOnly = true, dsout << "Quantification will be performed anew using existing identification info\n";
		else if (!memcmp(&(args[start]), "report-only ", 12)) ReportOnly = true, dsout << "Report will be generated using .quant files\n";
		else flag = true;

		if (!flag) {} else if (!memcmp(&(args[start]), "iter ", 5)) iN = Max(CalibrationIter + 4, std::stoi(args.substr(start + 5, std::string::npos))),
			dsout << "Number of iterations set to " << iN << "\n";
		else if (!memcmp(&(args[start]), "profiling-qvalue ", 17)) MaxProfilingQvalue = std::stod(args.substr(start + 17, std::string::npos)),
			dsout << "RT profiling q-value threshold set to " << MaxProfilingQvalue << "\n";
		else if (!memcmp(&(args[start]), "quant-qvalue ", 13)) MaxQuantQvalue = std::stod(args.substr(start + 13, std::string::npos)),
			dsout << "Q-value threshold for cross-run quantification set to " << MaxQuantQvalue << "\n";
		else if (!memcmp(&(args[start]), "protein-quant-qvalue ", 21)) ProteinQuantQvalue = std::stod(args.substr(start + 21, std::string::npos)),
			dsout << "Precursor Q-value threshold for protein quantification set to " << ProteinQuantQvalue << "\n";
		else if (!memcmp(&(args[start]), "top ", 4)) TopN = std::stoi(args.substr(start + 4, std::string::npos)),
			dsout << "Top " << TopN << " precursors will be used for protein quantification in each run\n";
		else if (!memcmp(&(args[start]), "out-lib-qvalue ", 15)) ReportQValue = std::stod(args.substr(start + 15, std::string::npos)),
			dsout << "Q-value threshold for spectral library generation set to " << ReportQValue << "\n";
		else if (!memcmp(&(args[start]), "rt-profiling ", 13)) RTProfiling = true, dsout << "RT profiling enabled\n";
		else if (!memcmp(&(args[start]), "prefix ", 7)) prefix = trim(args.substr(start + 7, end - start - 7)); // prefix added to input file names
		else if (!memcmp(&(args[start]), "ext ", 4)) ext = trim(args.substr(start + 4, end - start - 4)); // extension added to input file names
		else if (!memcmp(&(args[start]), "lc-all-scores ", 14)) LCAllScores = true, dsout << "All scores will be used by the linear classifier (not recommended)\n";
		else if (!memcmp(&(args[start]), "peak-center ", 12)) QuantMode = 1, dsout << "Fixed-width center of each elution peak will be used for quantification\n";
		else if (!memcmp(&(args[start]), "peak-boundary ", 14)) PeakBoundary = std::stod(args.substr(start + 14, std::string::npos)),
			dsout << "Peak boundary intensity factor set to " << PeakBoundary << "\n";
		else if (!memcmp(&(args[start]), "standardisation-scale ", 22)) StandardisationScale = std::stod(args.substr(start + 22, std::string::npos)),
			dsout << "Standardisation scale set to " << StandardisationScale << "\n";
		else if (!memcmp(&(args[start]), "no-ifs-removal ", 15)) NoIfsRemoval = true, dsout << "Interference removal from fragment elution curves disabled\n";
		else if (!memcmp(&(args[start]), "no-fr-selection ", 16)) NoFragmentSelectionForQuant = true, dsout << "Cross-run selection of fragments for quantification disabled (not recommended)\n";
		else if (!memcmp(&(args[start]), "restrict-fr ", 12)) RestrictFragments = true, dsout << "Certain fragments (based on the library annotation) will not be used when quantifying peptides\n";
		else if (!memcmp(&(args[start]), "no-fr-exclusion ", 16)) ExcludeSharedFragments = false, dsout << "Exclusion of fragments shared between heavy and light labelled peptides from quantification disabled\n";
		else if (!memcmp(&(args[start]), "peak-translation ", 17)) TranslatePeaks = true, dsout << "Translation of retention times between peptides within the same elution group enabled\n";
		else if (!memcmp(&(args[start]), "no-standardisation ", 19)) Standardise = false, dsout << "Scores will not be standardised for neural network training\n";
		else if (!memcmp(&(args[start]), "no-nn ", 6)) nnIter = INF, dsout << "Neural network classifier disabled\n";
		else if (!memcmp(&(args[start]), "nn-iter ", 8)) nnIter = Max(CalibrationIter + 3, std::stoi(args.substr(start + 8, std::string::npos))),
			dsout << "Neural network classifier will be used starting from the interation number " << nnIter << "\n";
		else if (!memcmp(&(args[start]), "nn-bagging ", 11)) nnBagging = std::stoi(args.substr(start + 11, std::string::npos)),
			dsout << "Neural network bagging set to " << nnBagging << "\n";
		else if (!memcmp(&(args[start]), "nn-epochs ", 10)) nnEpochs = std::stoi(args.substr(start + 10, std::string::npos)),
			dsout << "Neural network epochs number set to " << nnEpochs << "\n";
		else if (!memcmp(&(args[start]), "nn-learning-rate ", 17)) nnLearning = std::stod(args.substr(start + 17, std::string::npos)),
			dsout << "Neural network learning rate set to " << nnLearning << "\n";
		else if (!memcmp(&(args[start]), "nn-reg ", 7)) Regularisation = std::stod(args.substr(start + 7, std::string::npos)),
			dsout << "Neural network regularisation set to " << Regularisation << "\n";
		else if (!memcmp(&(args[start]), "nn-hidden ", 10)) nnHidden = Max(1, std::stoi(args.substr(start + 10, std::string::npos))),
			dsout << "Number of hidden layers set to " << nnHidden << "\n";
		else if (!memcmp(&(args[start]), "mass-acc-cal ", 13)) CalibrationMassAccuracy = std::stod(args.substr(start + 13, std::string::npos)) / 1000000.0,
			dsout << "Calibration mass accuracy set to " << CalibrationMassAccuracy << "\n";
		else if (!memcmp(&(args[start]), "fix-mass-acc ", 13)) ForceMassAcc = true;
		else if (!memcmp(&(args[start]), "mass-acc ", 9)) GlobalMassAccuracy = std::stod(args.substr(start + 9, std::string::npos)) / 1000000.0, ForceMassAcc = true;
		else if (!memcmp(&(args[start]), "mass-acc-ms1 ", 13)) GlobalMassAccuracyMs1 = std::stod(args.substr(start + 13, std::string::npos)) / 1000000.0, ForceMassAcc = true;
		else if (!memcmp(&(args[start]), "gen-acc ", 8)) GeneratorAccuracy = std::stod(args.substr(start + 8, std::string::npos)) / 1000000.0,
			dsout << "Fragmentation spectrum generator accuracy set to " << GeneratorAccuracy << "\n";
		else if (!memcmp(&(args[start]), "min-corr ", 9)) MinCorrScore = Min(std::stod(args.substr(start + 9, std::string::npos)), -1.1 + (double)TopF),
			dsout << "Only peaks with correlation sum exceeding " << MinCorrScore << " will be considered\n";
		else if (!memcmp(&(args[start]), "min-ms1-corr ", 13)) MinMs1Corr = Min(std::stod(args.substr(start + 13, std::string::npos)), 0.99),
			dsout << "Only peaks with MS1 profile correlation exceeding " << MinMs1Corr << " will be considered\n";
		else if (!memcmp(&(args[start]), "corr-diff ", 10)) MaxCorrDiff = Max(std::stod(args.substr(start + 10, std::string::npos)), E),
			dsout << "Peaks with correlation sum below " << MaxCorrDiff << " from maximum will not be considered\n";
		else if (!memcmp(&(args[start]), "peak-apex ", 10)) PeakApexEvidence = Min(std::stod(args.substr(start + 10, std::string::npos)), 0.99),
			dsout << "Peaks must have apex height which is at least " << PeakApexEvidence << " of the maximum within the m/z scan window\n";
		else if (!memcmp(&(args[start]), "norm-qvalue ", 12)) NormalisationQvalue = std::stod(args.substr(start + 11, std::string::npos)),
			dsout << "Q-value threshold for cross-run normalisation set to " << NormalisationQvalue << "\n";
		else if (!memcmp(&(args[start]), "norm-fraction ", 14)) NormalisationPeptidesFraction = std::stod(args.substr(start + 14, std::string::npos)),
			dsout << "Global normalisation peptides fraction set to " << NormalisationPeptidesFraction << "\n";
		else if (!memcmp(&(args[start]), "norm-radius ", 12)) LocNormRadius = std::stoi(args.substr(start + 12, std::string::npos)),
			dsout << "Local normalisation radius set to " << LocNormRadius << "\n";
		else if (!memcmp(&(args[start]), "global-norm ", 12)) LocalNormalisation = false, dsout << "Median-based local normalisation disabled\n";
		else if (!memcmp(&(args[start]), "no-rt-norm ", 11)) NoRTDepNorm = true, dsout << "Median-based RT-dependent local normalisation disabled\n";
		else if (!memcmp(&(args[start]), "sig-norm ", 9)) NoSigDepNorm = false, dsout << "Median-based signal-dependent local normalisation enabled\n";
#if Q1
		else if (!memcmp(&(args[start]), "q1-cal ", 7)) Q1Cal = true, dsout << "Q1 calibration enabled\n";
#endif
		else if (!memcmp(&(args[start]), "no-calibration ", 15)) Calibrate = false, dsout << "Mass calibration disabled\n";
		else if (!memcmp(&(args[start]), "mass-cal-bins ", 14)) MassCalBinsMax = std::stoi(args.substr(start + 14, std::string::npos)),
			dsout << "Maximum number of mass calibration bins set to " << MassCalBinsMax << "\n";
		else if (!memcmp(&(args[start]), "min-cal ", 8)) MinCal = std::stoi(args.substr(start + 8, std::string::npos)),
			dsout << "Minimum number of precursors identified at 10% FDR used for calibration set to " << MinCal << "\n";
		else if (!memcmp(&(args[start]), "min-class ", 10)) MinClassifier = std::stoi(args.substr(start + 10, std::string::npos)),
			dsout << "Minimum number of precursors identified at 10% FDR used for linear classifier training set to " << MinClassifier << "\n";
		else if (!memcmp(&(args[start]), "scanning-swath ", 15)) ForceScanningSWATH = true, dsout << "All runs will be analysed as Scanning SWATH runs\n";
		else if (!memcmp(&(args[start]), "regular-swath ", 14)) ForceNormalSWATH = true, dsout << "All runs will be analysed as regular SWATH runs\n";
#if Q1
		else if (!memcmp(&(args[start]), "no-q1 ", 6)) UseQ1 = false, dsout << "Q1 scores disabled\n";
		else if (!memcmp(&(args[start]), "use-q1 ", 8)) ForceQ1 = true, dsout << "Q1 scores will be used for regular SWATH runs\n";
#endif
		else dsout << "WARNING: unrecognised option [--" << trim(args.substr(start, end - start)) << "]\n";

		start = next;
	}
	if (temp_folder.size()) {
		char c = temp_folder[temp_folder.size() - 1];
		if (c != '/' && c != '\\') temp_folder += '/';
#ifdef CPP17
		if (!std::experimental::filesystem::exists(std::experimental::filesystem::path(temp_folder))) {
			dsout << "Cannot find the temp folder " << temp_folder << ". Specify an existing folder\n";
			exit(-1);
		}
#endif
	}

	if (prefix.length()) for (auto it = files.begin(); it != files.end(); it++) *it = prefix + *it;
	if (ext.length()) for (auto it = files.begin(); it != files.end(); it++) *it += ext;
	for (int i = 0; i < files.size(); i++) {
		auto extension = get_extension(files[i]);
#ifdef WIFFREADER
		if (extension == std::string(".wiff")) {
			if (wiff_dll == NULL && !skip_wiff) {
				wiff_dll = LoadLibrary(LPSTR("DIA-NN.Wiff.dll"));
				if (wiff_dll == NULL) {
					dsout << "Cannot load DIA-NN.Wiff.dll. Wiff files will be skipped.\n";
					skip_wiff = true;
				}
			}
			if (skip_wiff) continue;
		}
#endif
		auto ext = std::lower_bound(MSFormatsExt.begin(), MSFormatsExt.end(), extension);
		if (*ext != extension) {
			dsout << "WARNING: skipping " << files[i] << " - invalid raw MS data format\n";
			continue;
		}
		if (extension == std::string(".dia")) {
			if (!Convert) ms_files.push_back(files[i]);
			continue;
		}
		ms_files.push_back(files[i]);
	}
	if (!ms_files.size() && GenSpecLib) ExportLibrary = true, GenSpecLib = false;
	if (ms_files.size() < 2) RTProfiling = false;
	if (Predictor || ExportProsit) ExportLibrary = false;

	if (DisableRT) RTWindowedSearch = false;
	if (!ExportLibrary && !GenSpecLib && !Predictor && !ExportProsit) out_lib_file.clear();
	if (out_file.find_first_not_of(' ') == std::string::npos || out_file.find_first_of('\r') != std::string::npos || out_file.find_first_of('\n') != std::string::npos) out_file.clear();
	if (out_gene_file.find_first_not_of(' ') == std::string::npos || out_gene_file.find_first_of('\r') != std::string::npos || out_gene_file.find_first_of('\n') != std::string::npos) out_gene_file.clear();
	if (out_lib_file.find_first_not_of(' ') == std::string::npos || out_lib_file.find_first_of('\r') != std::string::npos || out_lib_file.find_first_of('\n') != std::string::npos) out_lib_file.clear();
	if (!out_file.size() && !out_gene_file.size() && !out_lib_file.size() && !IndividualReports) {
		IndividualReports = true;
		dsout << "No output files specified: reports will be generated separately for different runs (in the respective folders)\n";
	}
	if ((ExportLibrary || GenSpecLib || Predictor || ExportProsit) && !out_lib_file.size()) {
		ExportLibrary = GenSpecLib = Predictor = ExportProsit = false;
		dsout << "No output library file specified, library generation/export will not be performed\n";
	}
	
	if (UseQuant || QuantOnly) UseRTInfo = true;
	if (VarMods.size()) {
		MaxVarMods = Max(1, MaxVarMods);
		if (VarMods.size() >= 32) {
			VarMods.resize(31);
			dsout << "WARNING: only the first 31 variable modifications specified will be searched for\n";
		}
	}

	nnIter = Min(Max(Max(RTWindowIter, CalibrationIter), nnIter), iN);
	if (nnCrossVal) {
		int bagging = Max(2, nnBagging / 4) * 4;
		if (bagging != nnBagging) nnBagging = bagging, dsout << "Bagging factor set to " << bagging << " (required for NN cross validation)\n";
	}
	if (BatchMode) nnIter = Max(iN - 1, nnIter);
	MinBatch = Max(MinBatch, Min(MinClassifier, MinCal));
	if (ForceMassAcc) {
		dsout << "Mass accuracy will be fixed to " << GlobalMassAccuracy << " (MS2) and "
			<< GlobalMassAccuracyMs1 << " (MS1)\n";
	}

	if (!fasta_files.size()) {
		PGLevel = 0;
		if (FastaSearch) {
			FastaSearch = false;
			dsout << "WARNING: no FASTA provided, nothing to digest\n";
		}
	} else if (!FastaSearch && !lib_file.size() && !Convert) dsout << "WARNING: no spectral library provided, FASTA digest not enabled: nothing to do\n";
	if (fasta_files.size() && FastaSearch) {
		if (ExcludeSharedFragments) {
			dsout << "Exclusion of fragments shared between heavy and light peptides from quantification is not supported in library-free mode - disabled\n";
			ExcludeSharedFragments = false;
		}
		nnIter = Max(nnIter, iN - 1);
	}
	if (Reannotate && !Convert) {
		if (!lib_file.size()) Reannotate = false;
		if (!fasta_files.size() && lib_file.size()) Reannotate = false, dsout << "WARNING: no FASTA, cannot reannotate the library\n";
	}

	QFilter = Max(Max(Max(ReportQValue, MaxQuantQvalue), Max(ProteinQuantQvalue, MaxProfilingQvalue)), Max(NormalisationQvalue, ProteinIDQvalue));

	init_unimod();
	init_prediction();
	dsout << "\n";
}

template <class F, class T> void write_vector(F &out, std::vector<T> &v) {
	int size = v.size();
	out.write((char*)&size, sizeof(int));
	if (size) out.write((char*)&(v[0]), size * sizeof(T));
}

template <class F, class T> void read_vector(F &in, std::vector<T> &v) {
	int size = 0; in.read((char*)&size, sizeof(int));
	if (size) {
		v.resize(size);
		in.read((char*)&(v[0]), size * sizeof(T));
	}
}

template<class F> void write_string(F &out, std::string &s) {
	int size = s.size();
	out.write((char*)&size, sizeof(int));
	out.write((char*)&(s[0]), size);
}

template<class F> void read_string(F &in, std::string &s) {
	int size = 0; in.read((char*)&size, sizeof(int));
	if (size) {
		s.resize(size);
		in.read((char*)&(s[0]), size);
	}
}

template <class F, class T> void write_array(F &out, std::vector<T> &a) {
	int size = a.size();
	out.write((char*)&size, sizeof(int));
	for (int i = 0; i < size; i++) a[i].write(out);
}

template <class F, class T> void read_array(F &in, std::vector<T> &a) {
	int size = 0; in.read((char*)&size, sizeof(int)); 
	if (size) {
		a.resize(size);
		for (int i = 0; i < size; i++) a[i].read(in);
	}
}

template<class F> void write_strings(F &out, std::vector<std::string> &strs) {
	int size = strs.size();
	out.write((char*)&size, sizeof(int));
	for (int i = 0; i < size; i++) write_string(out, strs[i]);
}

template<class F> void read_strings(F &in, std::vector<std::string> &strs) {
	int size = 0; in.read((char*)&size, sizeof(int));
	if (size) {
		strs.resize(size);
		for (int i = 0; i < size; i++) read_string(in, strs[i]);
	}
}

char char_from_type[3] = { '?', 'b', 'y' };
std::string name_from_loss[loss_other + 1] = { std::string("noloss"), std::string("H2O"), std::string("NH3"), std::string("CO"), std::string("unknown"), std::string("unknown") };

class Ion {
public:
	char type, index, loss, charge; // index - number of AAs to the left from the cleavage site
	float mz;

	Ion() {};
	Ion(int _type, int _index, int _loss, int _charge, double _mz) {
		type = _type;
		index = _index;
		loss = _loss;
		charge = _charge;
		mz = _mz;
	}

	template<class T> void init(T &other) {
		type = other.type & 3;
		index = other.index;
		loss = other.loss;
		charge = other.charge;
		mz = other.mz;
	}

	bool operator == (const Ion &other) { return index == other.index && type == other.type && charge == other.charge && loss == other.loss; }
};

std::vector<Ion> generate_fragments(const std::vector<double> &sequence, int charge, int loss, int *cnt, int min_aas = 0) {
	int i;
	double c = (double)charge, curr, s = sum(&(sequence[0]), sequence.size());
	std::vector<Ion> result;

	for (i = 0, curr = 0.0; i < sequence.size() - 1; i++) {
		curr += sequence[i];
		double b = curr;
		double y = s - curr + proton + OH;

		if (i + 1 >= min_aas) result.push_back(Ion(fTypeB, i + 1, loss, charge, (b + c * proton - Loss[loss]) / c)), (*cnt)++;
		if (sequence.size() - i - 1 >= min_aas) result.push_back(Ion(fTypeY, i + 1, loss, charge, (y + c * proton - Loss[loss]) / c)), (*cnt)++;
	}
	return result;
}

std::vector<Ion> generate_all_fragments(const std::vector<double> &sequence, int charge, int max_loss, int *cnt) {
	int i;
	double c = (double)charge, curr, s = sum(&(sequence[0]), sequence.size());
	std::vector<Ion> result;

	for (int loss = loss_none; loss <= max_loss; loss++) {
		for (i = 0, curr = 0.0; i < sequence.size() - 1; i++) {
			curr += sequence[i];
			double b = curr;
			double y = s - curr + proton + OH;

			result.push_back(Ion(fTypeB, i + 1, loss, charge, (b + c * proton - Loss[loss]) / c)), (*cnt)++;
			result.push_back(Ion(fTypeY, i + 1, loss, charge, (y + c * proton - Loss[loss]) / c)), (*cnt)++;
		}
	}
	return result;
}

double max_gen_acc = 0.0;
Lock gen_acc_lock;
std::vector<Ion> recognise_fragments(const std::vector<double> &sequence, const std::vector<Product> &fragments, float &gen_acc, int pr_charge, bool full_spectrum = false, int max_charge = 19, int loss_cap = loss_N) {
	int i, j, cnt, tot, charge, loss, index;
	std::vector<Ion> result(fragments.size());
	std::vector<Ion> v;
	double delta, min, s = sum(&(sequence[0]), sequence.size());
	gen_acc = GeneratorAccuracy;

anew:
	cnt = tot = index = 0, charge = 1;
	loss = loss_none;
	for (i = 0; i < result.size(); i++) result[i].charge = 0;

start:
	v = generate_fragments(sequence, charge, loss, &cnt);
	if (charge == 1) v.push_back(Ion(fTypeY, sequence.size(), loss, pr_charge, (s + OH + (1.0 + double(pr_charge)) * proton - Loss[loss]) / double(pr_charge))); // non-fragmented

	for (i = 0; i < fragments.size(); i++) {
		if (result[i].charge) continue;
		double margin = fragments[i].mz * gen_acc;
		for (j = 0, min = margin; j < v.size(); j++) {
			delta = Abs(v[j].mz - fragments[i].mz);
			if (delta < min) min = delta, index = j;
		}
		if (min < margin) {
			tot++;
			result[i].init(v[index]);
			if (tot >= fragments.size()) goto stop;
		}
	}

	loss++;
	if (loss == loss_cap) loss = loss_none, charge++;
	if (charge <= max_charge && tot < fragments.size()) goto start;

stop:
	if (tot < fragments.size()) {
		if (full_spectrum) return result;
		gen_acc *= 2.0;
		if (Verbose >= 1 && gen_acc > max_gen_acc)
			dsout << "WARNING: not all fragments recognised; generator accuracy threshold increased to " << gen_acc << "\n";
		if (gen_acc > max_gen_acc) {
			while (!gen_acc_lock.set()) {}
			max_gen_acc = gen_acc;
			gen_acc_lock.free();
		}
		v.clear();
		goto anew;
	}

	return result;
}

std::vector<Ion> recognise_fragments(const std::string &name, const std::vector<Product> &fragments, float &gen_acc, int pr_charge, bool full_spectrum = false, int max_charge = 19, int loss_cap = loss_N) {
	int i;
	if (fragments.size()) if (fragments[0].type & 3) {
		std::vector<Ion> result(fragments.size());
		for (i = 0; i < result.size(); i++) 
			if (fragments[i].type & 3) result[i].init(fragments[i]); else break;
		if (i == result.size()) return result;
	}
	auto sequence = get_sequence(name);
	return recognise_fragments(sequence, fragments, gen_acc, pr_charge, full_spectrum, max_charge, loss_cap);
}


std::vector<Ion> generate_fragments(const std::vector<double> &sequence, const std::vector<Ion> &pattern, int pr_charge) {
	int i, j, charge = 1, loss = loss_none, cnt = 0, tot = 0;
	double s = sum(&(sequence[0]), sequence.size());
	std::vector<Ion> result(pattern.size());
	std::vector<Ion> v;

start:
	v = generate_fragments(sequence, charge, loss, &cnt); i = 0;
	if (charge == 1) v.push_back(Ion(fTypeY, sequence.size(), loss, pr_charge, (s + OH + (1.0 + double(pr_charge)) * proton - Loss[loss]) / double(pr_charge))); // non-fragmented

	for (auto pos = pattern.begin(); pos != pattern.end(); pos++, i++) {
		for (j = 0; j < v.size(); j++) if (v[j] == *pos) {
			result[i].init(v[j]); tot++;
			break;
		}
	}
	loss++;
	if (loss == loss_N) loss = loss_none, charge++;
	if (charge < 20 && tot < pattern.size()) goto start;

	return result;
}

bool ProfileSpectrumWarning = false;

struct FI {
	float value;
	int index;

	friend bool inline operator < (const FI& left, const FI& right) { return left.value < right.value; }
	friend bool inline operator > (const FI& left, const FI& right) { return left.value > right.value; }
};

struct Scan {
	double RT = 0.0, window_low = 0.0, window_high = 0.0;
	int n = 0, type = 0;
	Peak * peaks = NULL;

	Scan() { }

	inline int size() { return n; }

	inline friend bool operator < (const Scan& left, const Scan& right) { return left.RT < right.RT; }

	inline bool has(float mz) { return (window_low <= mz && window_high > mz); }

	inline bool has(float mz, float margin) { return (window_low <= mz + margin && window_high > mz - margin); }

	template <bool get_mz> inline double level(int &low, int &high, float mz, float accuracy, float * peak_mz = NULL) {
		int i;
		float v, s, margin, min, max;
		margin = mz * accuracy;
		min = mz - margin, max = mz + margin;
		if (!size()) {
			if (get_mz) *peak_mz = 0.0;
			return 0.0;
		}

		while (high > low) {
			int middle = (high + low) >> 1;
			v = peaks[middle].mz;
			if (v < max) {
				if (v > min) {
					s = peaks[middle].height;
					if (get_mz) *peak_mz = peaks[middle].mz;
					for (i = middle + 1; i < high && peaks[i].mz < max; i++) if (peaks[i].height > s) {
						s = peaks[i].height;
						if (get_mz) *peak_mz = peaks[i].mz;
					}
					if (i < high) high = i;
					for (i = middle - 1; i >= low && peaks[i].mz > min; i--) if (peaks[i].height > s) {
						s = peaks[i].height;
						if (get_mz) *peak_mz = peaks[i].mz;
					}
					if (i >= low) low = i + 1;
					return s;
				}
				low = middle + 1;
			} else high = middle;
		}
		return 0.0;
	}

	template <bool get_mz> inline double level(float mz, float accuracy, float * peak_mz = NULL) {
		int low = 0, high = size();
		return level<get_mz>(low, high, mz, accuracy, peak_mz);
	}

	inline double level(float mz, float accuracy, float ref_height) { // finds peak with height closest to the ref_height
		int i, low = 0, high = size();
		float v, s, delta, u, margin, min, max;
		margin = mz * accuracy;
		min = mz - margin, max = mz + margin;
		if (!size()) return 0.0;

		while (high > low) {
			int middle = (high + low) >> 1;
			v = peaks[middle].mz;
			if (v < max) {
				if (v > min) {
					s = peaks[middle].height;
					delta = Abs(s - ref_height);
					for (i = middle + 1; i < high && peaks[i].mz < max; i++) if ((u = Abs(peaks[i].height - ref_height)) < delta) {
						s = peaks[i].height;
						delta = u;
					}
					if (i < high) high = i;
					for (i = middle - 1; i >= low && peaks[i].mz > min; i--) if ((u = Abs(peaks[i].height - ref_height)) < delta) {
						s = peaks[i].height;
						delta = u;
					}
					if (i >= low) low = i + 1;
					return s;
				}
				low = middle + 1;
			} else high = middle;
		}
		return 0.0;
	}

#if (HASH > 0)
	unsigned int hash() {
		unsigned int res = hashS((float)RT) ^ hashS((float)window_low) ^ hashS((float)window_high);
		for (int i = 0; i < n; i++) res ^= peaks[i].hash();
		return res;
	}
#endif

	void bin_peaks(FI * binned, int bins_per_Da) {
		int i;
		FI * fi = NULL;
		float bpd = (float)bins_per_Da, bpdi = 1.0 / bpd, max = -INF;

		if (!n) return;
		for (i = 0; i < n; i++) {
			auto &p = peaks[i];
			if (p.mz >= max) {
				int bin = floor(p.mz * bpd);
				max = bpdi * float(bin + 1);
				fi = binned + bin;
			}
			if (p.height > fi->value) fi->value = p.height, fi->index = i;
		}
	}

	void bin_peaks(std::vector<FI>& binned, int bins_per_Da) {
		if (!n) return;
		binned.resize(floor(peaks[n - 1].mz * (float)bins_per_Da) + 1);
		bin_peaks(&(binned[0]), bins_per_Da);
	}

	void top_peaks(Peak * top, int N, std::vector<Peak>& temp) {
		if (N >= n) {
			N = n;
			memcpy(top, peaks, N * sizeof(Peak));
			return;
		}

		temp.resize(n);
		memcpy(&(temp[0]), peaks, n * sizeof(Peak));
		std::sort(temp.begin(), temp.end(), [&](const Peak &left, const Peak &right) { return left.height > right.height; });
		temp.resize(N);
		std::sort(temp.begin(), temp.end());
		memcpy(top, &(temp[0]), N * sizeof(Peak));
	}

	void top_peaks(std::vector<Peak> &top, int N, std::vector<Peak>& temp) {
		top.resize(N);
		top_peaks(&(top[0]), N, temp);
	}
};

inline int mass_bin(float mz, double acc) { return floor(log(mz) / acc); }
inline float bin_to_mass(int bin, double acc) { return exp(acc * (double)bin); }
inline float binned_float(float mz, double acc) { return bin_to_mass(mass_bin(mz, acc), acc); }

class Isoform {
public:
	std::string id;
	mutable std::string name, gene, description;
	mutable std::set<int> precursors; // precursor indices in the library
	int name_index = 0, gene_index = 0;
	bool swissprot = true;

	Isoform() {}
	Isoform(const std::string &_id) { id = _id; }
	Isoform(const std::string &_id, const std::string &_name, const std::string &_gene, const std::string &_description, bool _swissprot) {
		id = _id;
		name = _name;
		gene = _gene;
		description = _description;
		swissprot = _swissprot;
	}

	friend inline bool operator < (const Isoform &left, const Isoform &right) { return left.id < right.id; }

	template <class F> void write(F &out) {
		int sp = swissprot, size = precursors.size();
		out.write((char*)&sp, sizeof(int));
		out.write((char*)&size, sizeof(int));
		write_string(out, id);
		write_string(out, name);
		write_string(out, gene);
		out.write((char*)&name_index, sizeof(int));
		out.write((char*)&gene_index, sizeof(int));
		for (auto it = precursors.begin(); it != precursors.end(); it++) out.write((char*)&(*it), sizeof(int));
	}

	template <class F> void read(F &in) {
		int sp = 0, size = 0;
		in.read((char*)&sp, sizeof(int));
		in.read((char*)&size, sizeof(int));
		swissprot = sp;

		read_string(in, id);
		read_string(in, name);
		read_string(in, gene);
		in.read((char*)&name_index, sizeof(int));
		in.read((char*)&gene_index, sizeof(int));
		precursors.clear();
		for (int i = 0; i < size; i++) {
			int pr = -1;
			in.read((char*)&pr, sizeof(int));
			if (pr >= 0) precursors.insert(pr);
		}
	}
};

class PG {
public:
	std::string ids;
	mutable std::string names, genes;
	mutable std::vector<int> precursors; // precursor indices in the library
	mutable std::set<int> proteins;
	std::vector<int> name_indices, gene_indices;

	PG() {}

	PG(std::string _ids) {
		ids = _ids;
	}

	friend inline bool operator < (const PG &left, const PG &right) { return left.ids < right.ids; }

	void annotate(std::vector<Isoform> &_proteins, std::vector<std::string> &_names, std::vector<std::string> &_genes) {
		std::set<int> name, gene;
		std::string word;
		if (_names.size()) {
			for (auto &p : proteins) if (_names[_proteins[p].name_index].size()) name.insert(_proteins[p].name_index);
			if (names.size() && !proteins.size()) {
				std::stringstream list(names);
				while (std::getline(list, word, ';')) {
					word = fast_trim(word);
					auto pos = std::lower_bound(_names.begin(), _names.end(), word);
					if (pos != _names.end()) if (*pos == word) name.insert(std::distance(_names.begin(), pos));
				}
			}
			names.clear(); name_indices.clear(); name_indices.insert(name_indices.begin(), name.begin(), name.end());
			if (name_indices.size()) for (int i = 0; i < name_indices.size(); i++) names += (names.size() ? std::string(";") : std::string("")) + _names[name_indices[i]];
		}
		if (_genes.size()) {
			for (auto &p : proteins) if (_genes[_proteins[p].gene_index].size()) gene.insert(_proteins[p].gene_index);
			if (genes.size() && !proteins.size()) {
				std::stringstream list(genes);
				while (std::getline(list, word, ';')) {
					word = fast_trim(word);
					auto pos = std::lower_bound(_genes.begin(), _genes.end(), word);
					if (pos != _genes.end()) if (*pos == word) gene.insert(std::distance(_genes.begin(), pos));
				}
			}
			genes.clear(); gene_indices.clear(); gene_indices.insert(gene_indices.begin(), gene.begin(), gene.end());
			if (gene_indices.size()) for (int i = 0; i < gene_indices.size(); i++) genes += (genes.size() ? std::string(";") : std::string("")) + _genes[gene_indices[i]];
		}
	}

	template <class F> void write(F &out) {
		int size_p = proteins.size();
		out.write((char*)&size_p, sizeof(int));
		write_string(out, ids);
		write_string(out, names);
		write_string(out, genes);
		write_vector(out, name_indices);
		write_vector(out, gene_indices);
		write_vector(out, precursors);
		for (auto it = proteins.begin(); it != proteins.end(); it++) out.write((char*)&(*it), sizeof(int));
	}

	template <class F> void read(F &in) {
		int sp = 0, size_p = 0;
		in.read((char*)&size_p, sizeof(int));

		read_string(in, ids);
		read_string(in, names);
		read_string(in, genes);
		read_vector(in, name_indices);
		read_vector(in, gene_indices);
		read_vector(in, precursors);
		for (int i = 0; i < size_p; i++) {
			int p = -1;
			in.read((char*)&p, sizeof(int));
			if (p >= 0) proteins.insert(p);
		}
	}
};

enum {
	qTotal, qFiltered,
	qN
};

class Fragment {
public:
	mutable int index; // index of the fragment in the library (i.e. its number among the fragments of the precursor)
	mutable float quantity[qN];
	mutable float corr;
};

class PrecursorEntry {
public:
	int index; // index of the precursor in the spectral library
	int run_index;
	mutable bool decoy_found;
	mutable int apex, peak, best_fragment, peak_width;
	mutable float RT, iRT, predicted_RT, predicted_iRT, best_fr_mz, profile_qvalue, qvalue, decoy_qvalue, protein_qvalue, quantity, ratio, level, norm;
	mutable float pg_quantity, pg_norm, gene_quantity, gene_norm, gene_quantity_u, gene_norm_u;
	mutable float evidence, decoy_evidence, ms1_corr, cscore, decoy_cscore;
#if REPORT_SCORES
	float scores[pN], decoy_scores[pN];
#endif
	friend inline bool operator < (const PrecursorEntry &left, const PrecursorEntry &right) { return left.index < right.index; }
};

class QuantEntry {
public:
	int index = -1, window, fr_n = 0;
    mutable PrecursorEntry pr;
    mutable Fragment fr[TopF];
#if ELUTION_PROFILE
	mutable std::pair<float, float> ms1_elution_profile[2 * ElutionProfileRadius + 1]; // (retention time, intensity) pairs
#endif

	QuantEntry() {  }

	template <class F> void write(F &out) { out.write((char*)this, sizeof(QuantEntry)); }

	template <class F> void read(F &in) { in.read((char*)this, sizeof(QuantEntry)); }

	friend inline bool operator < (const QuantEntry &left, const QuantEntry &right) { return left.index < right.index; }
};

class DecoyEntry {
public:
	int index = -1;
	float qvalue;

	DecoyEntry() {  }

	template <class F> void write(F &out) {
		out.write((char*)this, sizeof(DecoyEntry));
	}

	template <class F> void read(F &in) {
		in.read((char*)this, sizeof(DecoyEntry));
	}

	friend inline bool operator < (const DecoyEntry &left, const DecoyEntry &right) { return left.index < right.index; }
};

class Quant {
public:
	RunStats RS;
	std::vector<QuantEntry> entries;
	std::vector<DecoyEntry> decoys; // used only for cross-run q-value calculation - for spectral library generation from DIA data; 
	std::vector<int> proteins; // protein ids at <= ProteinIDQvalue
	double weights[pN], guide_weights[pN];
	double tandem_min, tandem_max;

	int run_index, lib_size;
	double MassAccuracy = GlobalMassAccuracy, MassAccuracyMs1 = GlobalMassAccuracyMs1;
	std::vector<double> MassCorrection, MassCorrectionMs1, MassCalSplit, MassCalSplitMs1, MassCalCenter, MassCalCenterMs1;
	double Q1Correction[2];

	template <class F> void write(F &out) {
		out.write((char*)&run_index, sizeof(int));
		out.write((char*)&lib_size, sizeof(int));
		out.write((char*)(&RS), sizeof(RunStats));
		out.write((char*)weights, pN * sizeof(double));
		out.write((char*)guide_weights, pN * sizeof(double));
		out.write((char*)&tandem_min, sizeof(double));
		out.write((char*)&tandem_max, sizeof(double));
		out.write((char*)&MassAccuracy, sizeof(double));
		out.write((char*)&MassAccuracyMs1, sizeof(double));
		out.write((char*)Q1Correction, 2 * sizeof(double));
		write_vector(out, MassCorrection);
		write_vector(out, MassCorrectionMs1);
		write_vector(out, MassCalSplit);
		write_vector(out, MassCalSplitMs1);
		write_vector(out, MassCalCenter);
		write_vector(out, MassCalCenterMs1);

		write_vector(out, proteins);

        int size = entries.size();
        out.write((char*)&size, sizeof(int));
        for (int i = 0; i < size; i++) entries[i].write(out);

		size = decoys.size();
		out.write((char*)&size, sizeof(int));
		for (int i = 0; i < size; i++) decoys[i].write(out);
    }

	template <class F> void read_meta(F &in, int _lib_size = 0) {
		if (in.fail()) {
			dsout << "ERROR: cannot read the quant file\n";
			exit(0);
		}
		in.read((char*)&run_index, sizeof(int));
		in.read((char*)&lib_size, sizeof(int));
		if (_lib_size > 0 && lib_size != _lib_size) {
			dsout << "ERROR: a .quant file was obtained using a different spectral library / different library-free search settings\n";
			exit(0);
		}
		in.read((char*)(&RS), sizeof(RunStats));
		in.read((char*)weights, pN * sizeof(double));
		in.read((char*)guide_weights, pN * sizeof(double));
		in.read((char*)&tandem_min, sizeof(double));
		in.read((char*)&tandem_max, sizeof(double));
		in.read((char*)&MassAccuracy, sizeof(double));
		in.read((char*)&MassAccuracyMs1, sizeof(double));
		in.read((char*)Q1Correction, 2 * sizeof(double));
		read_vector(in, MassCorrection);
		read_vector(in, MassCorrectionMs1);
		read_vector(in, MassCalSplit);
		read_vector(in, MassCalSplitMs1);
		read_vector(in, MassCalCenter);
		read_vector(in, MassCalCenterMs1);
	}

	void read_meta(const std::string &file, int _lib_size = 0) {
		std::ifstream in(file, std::ifstream::binary);
		if (in.fail() || temp_folder.size()) {
			auto pos = file.find_last_of('.');
			if (pos >= 1) in = std::ifstream(location_to_file_name(file.substr(0, pos - 1)) + std::string(".quant"), std::ifstream::binary);
		}
		read_meta(in, _lib_size);
		in.close();
	}

	template <class F> void read(F &in, int _lib_size = 0) {
		read_meta(in, _lib_size);
		read_vector(in, proteins);
		int size; in.read((char*)&size, sizeof(int));
		entries.resize(size);
		for (int i = 0; i < size; i++) entries[i].read(in);
		in.read((char*)&size, sizeof(int));
		decoys.resize(size);
		for (int i = 0; i < size; i++) decoys[i].read(in);
	}
};

std::vector<Quant> quants;

class Profile {
public:
	std::vector<QuantEntry> entries;

	Profile(std::vector<std::string> &files, int lib_size = 0) {
		std::vector<float> tq, dq;
		entries.resize(MaxLibSize);

		bool decoy_list = FastaSearch || GenSpecLib;
		if (decoy_list) dq.resize(MaxLibSize, 1.0), tq.resize(MaxLibSize, 1.0);

		for (int i = 0; i < files.size(); i++) {
			Quant Qt;
			auto Q = &Qt;
			if (!QuantInMem) {
				std::ifstream in(files[i] + std::string(".quant"), std::ifstream::binary);
				if (in.fail() || temp_folder.size()) in = std::ifstream(location_to_file_name(files[i]) + std::string(".quant"), std::ifstream::binary);
				Q->read(in, lib_size);
				in.close();
			} else Q = &(quants[i]);

			for (auto it = Q->entries.begin(); it != Q->entries.end(); it++) {
				auto pos = it->pr.index;
				if (pos >= entries.size()) entries.resize(entries.size() + entries.size() / 2);
				if (entries[pos].index < 0) entries[pos] = *it;
				else if (it->pr.qvalue < entries[pos].pr.qvalue) entries[pos].pr = it->pr;

				if (decoy_list) {
					if (pos >= tq.size()) tq.resize(pos + 1 + pos / 2, 1.0);
					if (it->pr.qvalue <= ReportQValue && it->pr.qvalue < tq[pos]) tq[pos] = it->pr.qvalue;
				}
			}
			if (decoy_list) {
				for (auto it = Q->decoys.begin(); it != Q->decoys.end(); it++) {
					auto pos = it->index;

					if (pos >= dq.size()) dq.resize(pos + 1 + pos / 2, 1.0);
					if (it->qvalue <= ReportQValue && it->qvalue < dq[pos]) dq[pos] = it->qvalue;
				}
			}
			if (!QuantInMem) std::vector<QuantEntry>().swap(Q->entries), std::vector<DecoyEntry>().swap(Q->decoys);
		}
		if (decoy_list) {
			std::sort(dq.begin(), dq.end());
			auto pos = std::lower_bound(dq.begin(), dq.end(), ReportQValue);
			dq.resize(std::distance(dq.begin(), pos));
			std::sort(tq.begin(), tq.end());
			pos = std::lower_bound(tq.begin(), tq.end(), ReportQValue);
			tq.resize(std::distance(tq.begin(), pos));

			std::map<float, float> cs_qv;

			for (int i = 0; i < tq.size(); i++) {
				int n_targets = i + 1;
				int n_decoys = std::distance(dq.begin(), std::lower_bound(dq.begin(), dq.end(), tq[i]));
				float q = Min(1.0, ((double)Max(1, n_decoys)) / (double)Max(1, n_targets));

				auto pair = std::pair<float, float>(-tq[i], q);
				auto pos = cs_qv.insert(pair);
				if (pos.second) {
					if (pos.first != cs_qv.begin() && std::prev(pos.first)->second < pair.second) pos.first->second = std::prev(pos.first)->second;
					else for (auto jt = std::next(pos.first); jt != cs_qv.end() && jt->second > pair.second; jt++) jt->second = pair.second;
				}
			}

			for (auto &e : entries) {
				if (e.index < 0) continue;
				auto pos = cs_qv.lower_bound(-e.pr.qvalue);
				if (pos == cs_qv.end()) e.pr.profile_qvalue = 1.0;
				else e.pr.profile_qvalue = pos->second;
			}
		}
	}
};

struct Elution {
	int peak, apex, best_fragment;
	float score;
};

class Peptide {
public:
	int index = 0, charge = 0, length = 0, no_cal = 0;
	float mz = 0.0, iRT = 0.0, sRT = 0.0, lib_qvalue = 0.0;
	std::vector<Product> fragments;

	void init(float _mz, float _iRT, int _charge, int _index) {
		mz = _mz;
		iRT = _iRT;
		sRT = 0.0;
		charge = _charge;
		index = _index;
	}

	inline void free() {
		std::vector<Product>().swap(fragments);
	}

	template <class F> void write(F &out) {
		out.write((char*)&index, sizeof(int));
		out.write((char*)&charge, sizeof(int));
		out.write((char*)&length, sizeof(int));

		out.write((char*)&mz, sizeof(float));
		out.write((char*)&iRT, sizeof(float));
		out.write((char*)&sRT, sizeof(float));

		write_vector(out, fragments);
	}

	template <class F> void read(F &in) {
		in.read((char*)&index, sizeof(int));
		in.read((char*)&charge, sizeof(int));
		in.read((char*)&length, sizeof(int));

		in.read((char*)&mz, sizeof(float));
		in.read((char*)&iRT, sizeof(float));
		in.read((char*)&sRT, sizeof(float));

		read_vector(in, fragments);
	}

#if (HASH > 0)
	unsigned int hash() {
		unsigned int res = hashS(mz) ^ hashS(iRT) ^ hashS(sRT) ^ hashS(lib_qvalue);
		for (int i = 0; i < fragments.size(); i++) res ^= fragments[i].hash();
		return res;
	}
#endif
};

class FastaEntry {
public:
	std::string stripped;
	mutable std::string ids, names, genes;
	FastaEntry(const std::string &_stripped, const std::string &_ids, const std::string &_names, const std::string &_genes) { stripped = _stripped, ids = _ids, names = _names, genes = _genes; }
	friend inline bool operator < (const FastaEntry &left, const FastaEntry &right) { return left.stripped < right.stripped; }
};

bool name_included(const std::string &names, const std::string &name) {
	int start = 0, m = names.size(), n = name.size();
	if (!m || !n) return false;
	const char *p = &(name[0]);
	while (start + n <= m) {
		if (!memcmp(&(names[start]), p, n)) {
			if (start + n >= m) return true;
			if (names[start + n] == ';') return true;
		}
		while (start < m - n) {
			start++;
			if (names[start] == ';') break;
		}
		start++;
	}
	return false;
}

class Fasta {
public:
	std::string name;
	std::vector<std::string> peptides, sequences, ids;
	std::vector<std::pair<std::string, std::vector<int> > > dict;

	std::vector<FastaEntry> entries;
	std::vector<Isoform> proteins;

	void annotation_sgd(std::string &header, std::string &id, std::string &name, std::string &gene, std::string &description, bool * swissprot) {
		*swissprot = (header.find("Verified ORF") != std::string::npos);
		int start = 0, end2;
		auto end = header.find_first_of(' ');
		if (end == std::string::npos) id = fast_trim(header.substr(start + 1)), name = description = "";
		else {
			id = fast_trim(header.substr(start + 1, end - start - 1));
			for (end2 = end + 1; end2 < header.size() && header[end2] != ' '; end2++);
			end2 = header.find_first_of(' ', end2);
			name = fast_trim(header.substr(end + 1, end2 - end - 1));
			end = header.find_first_of('"', end2);
			if (end == std::string::npos) description = "";
			else {
				end2 = header.find_first_of('"', end + 1);
				if (end2 == std::string::npos) description = "";
				else description = header.substr(end + 1, end2 - end - 1);
			}
		}
		gene = name;
	}

	void annotation(std::string &header, std::string &id, std::string &name, std::string &gene, std::string &description, bool * swissprot) {
		if (header.size() < 3) {
			id = "", name = "", gene = "", description = "", *swissprot = false;
			return;
		}
		auto start = header.find("|");
		if (start == std::string::npos || (start != 3 && !(header[1] == 's' && header[2] == 'p') && !(header[1] == 't' && header[2] == 'r'))) {
			annotation_sgd(header, id, name, gene, description, swissprot);
			return;
		}
		auto end = header.find("|", start + 1);
		if (end == std::string::npos) end = header.find(" ", start + 1), name = description = "";
		else {
			auto end2 = header.find(" ", end + 1);
			if (end2 != std::string::npos) {
				name = fast_trim(header.substr(end + 1, end2 - end - 1));
				auto end3 = header.find("=", end2 + 1);
				if (end3 != std::string::npos) {
					while (header[end3] != ' ' && end3 > end2 + 2) end3--;
					if (end3 > end2 + 2) description = fast_trim(header.substr(end2 + 1, end3 - end2 - 1));
					else description = "";
				} else description = fast_trim(header.substr(end2 + 1));
			} else name = fast_trim(header.substr(end + 1)), description = "";
		}
		if (end != std::string::npos) id = fast_trim(header.substr(start + 1, end - start - 1));
		else id = fast_trim(header.substr(start + 1));
		start = header.find("GN=");
		if (start == std::string::npos) gene = "";
		else {
			end = header.find(" ", start + 3);
			if (end != std::string::npos) gene = fast_trim(header.substr(start + 3, end - start - 3));
			else gene = fast_trim(header.substr(start + 3));
		}
		*swissprot = (header[1] == 's' && header[2] == 'p');
		if (SpeciesGenes && gene.size() && name.size()) {
			int pos = name.find_last_of('_');
			if (pos != std::string::npos) gene += name.substr(pos);
		}
	}

	bool load_proteins(std::vector<std::string> &files) {
		std::set<std::string> protein_ids;
		for (auto file : files) {
			if (Verbose >= 1) Time(), dsout << "Loading protein annotations from FASTA " << file << "\n";
			std::ifstream in(file);
			if (in.fail()) {
				dsout << "cannot read the file\n";
				return false;
			}
			bool swissprot;
			std::string line, id, name, gene, description;
			while (getline(in, line)) if (line[0] == '>') {
				annotation(line, id, name, gene, description, &swissprot);
				auto ins = protein_ids.insert(id);
				if (ins.second) proteins.push_back(Isoform(id, name, gene, description, swissprot));
			}
			in.close();
		}
		std::sort(proteins.begin(), proteins.end());
		return true;
	}

	bool load(std::vector<std::string> &files, std::vector<std::string> &filter_files) {
		std::set<FastaEntry> unique;
		std::set<std::string> protein_ids, filter_peptides;
		std::vector<std::string> filters;
		std::vector<int> cuts, vec(1);
		std::map<std::string, std::vector<int> > dic;
		std::string line, sequence, peptide, id, name, gene, description, next_id, next_name, next_gene;
		bool swissprot;

		for (auto &filter : filter_files) {
			std::ifstream in(filter);
			if (!in.fail()) {
				while (getline(in, line)) if (line[0] != '>') {
					auto trimmed = fast_trim(line);
					trimmed.erase(std::remove(trimmed.begin(), trimmed.end(), '\t'), trimmed.end());
					filter_peptides.insert(trimmed);
				}
				in.close();
				if (Verbose >= 1) Time(), dsout << "Loaded FASTA filter " << filter << "\n";
			} else dsout << "WARNING: couldn't open the FASTA filter file " << filter << "\n";
			filters.insert(filters.begin(), filter_peptides.begin(), filter_peptides.end());
			filter_peptides.clear();
		}

		for (auto file : files) {
			if (Verbose >= 1) Time(), dsout << "Loading FASTA " << file << "\n";
			std::ifstream in(file);
			if (in.fail()) {
				dsout << "cannot read the file\n";
				return false;
			}

			sequence.clear(), id.clear(), name.clear(), gene.clear(), peptide.clear();
			bool skip = false, eof = !getline(in, line);
			if (!eof) while (true) {
				bool new_prot = (line.size() != 0);
				if (new_prot) new_prot &= (line[0] == '>');
				if (new_prot || eof) {
					if (new_prot) {
						annotation(line, next_id, next_name, next_gene, description, &swissprot);
						auto ins = protein_ids.insert(next_id);
						skip = !ins.second; // dublicate protein?
						if (!skip) proteins.push_back(Isoform(next_id, next_name, next_gene, description, swissprot));
						if (FastaProtDuplicates) skip = false;
						if (ForceSwissProt) skip = skip || !swissprot;
					}

					cuts.clear();
					cuts.push_back(-1);
					vec[0] = sequences.size();
					for (int i = 0; i < sequence.length(); i++) {
						char aa = sequence[i];
						if (aa >= 'A' && aa <= 'Z') {
							bool cut = (i == sequence.length() - 1);
							if (!cut) {
								for (auto &s : CutAfter) if (s == aa) { cut = true; break; }
								if (cut && i + 1 < sequence.length()) {
									char next_aa = sequence[i + 1];
									for (auto &s : NoCutBefore) if (s == next_aa) { cut = false; break; }
								}
							}
							if (cut) cuts.push_back(i);
						}
						if (i < sequence.length() - 4) {
							auto ins = dic.insert(std::pair<std::string, std::vector<int> >(sequence.substr(i, 5), vec));
							if (!ins.second) ins.first->second.push_back(vec[0]);
						}
					}
					sequences.push_back(sequence);
					ids.push_back(id);

					int max_miss = MissedCleavages;
					if (max_miss >= cuts.size()) max_miss = cuts.size() - 1;
					for (int miss = 0; miss <= max_miss; miss++) {
						for (int start = 0; start < cuts.size() - miss - 1; start++) {
							bool first = (start == 0);
							int len = cuts[start + miss + 1] - cuts[start];
							if (len < MinPeptideLength || len > MaxPeptideLength) continue;
							auto stripped = sequence.substr(cuts[start] + 1, cuts[start + miss + 1] - cuts[start]);
							bool skip_peptide = false;
							for (auto s : stripped) if (s < 'A' || s > 'Z') skip_peptide = true;
							if (skip_peptide) continue; // handle undefined amino acids in SGD FASTA
							if (filters.size()) {
								auto pos = std::lower_bound(filters.begin(), filters.end(), stripped);
								if (pos == filters.end()) continue;
								if (*pos != stripped) continue;
							}
						insert:
							auto ins = unique.insert(FastaEntry(stripped, id, name, gene));
							if (!ins.second) {
								if (!name_included(ins.first->ids, id)) ins.first->ids += std::string(";") + id;
								if (!name_included(ins.first->names, name)) ins.first->names += std::string(";") + name;
								if (!name_included(ins.first->genes, gene)) ins.first->genes += std::string(";") + gene;
							}

							if (first && NMetExcision && (stripped[0] == 'M' || stripped[0] == 'X') && stripped.length() >= MinPeptideLength + 1) {
								stripped = stripped.substr(1);
								first = false;
								goto insert;
							}
						}
					}
					sequence.clear();
					if (eof) break;
					id = next_id, name = next_name, gene = next_gene;
				} else if (!skip) sequence += line;
				if (!getline(in, line)) eof = true, line.clear();
			}
			in.close();
		}
		dict.insert(dict.begin(), dic.begin(), dic.end()); dic.clear();
		entries.insert(entries.begin(), unique.begin(), unique.end());

		std::vector<int> ind(MaxVarMods), type(MaxVarMods), mod, mod_cnt;
		for (auto &seq : unique) {
			peptide.clear();
			auto &s = seq.stripped;
			if (!FixedMods.size()) peptide = s;
			else {
				char lc = to_lower(s[0]);
				for (auto &fmod : FixedMods) for (int x = 0; x < fmod.aas.size(); x++) if (fmod.aas[x] == lc) { // N-terminal modification: encoded with a lower-case letter for the amino acid
					peptide += '(' + fmod.name + ')';
					break;
				}
				for (auto &aa : s) {
					peptide.push_back(aa);
					for (auto &fmod : FixedMods) for (int x = 0; x < fmod.aas.size(); x++) if (fmod.aas[x] == aa) {
						peptide += '(' + fmod.name + ')';
						break;
					}
				}
			}
			peptides.push_back(peptide);

			if (VarMods.size()) {
				int i, j, k, cnt, l = s.size(), m = Min(MaxVarMods, l);
				mod.resize(l), mod_cnt.resize(l);
				char lc = to_lower(s[0]);
				for (int y = 0; y < VarMods.size(); y++) { // N-terminal modification: encoded with a lower-case letter for the amino acid
					auto &vmod = VarMods[y];
					for (int x = 0; x < vmod.aas.size(); x++) if (vmod.aas[x] == lc) {
						mod[0] |= (1 << y), mod_cnt[0]++;
						break;
					}
				}
				for (i = 0; i < l; i++) {
					mod[i] = mod_cnt[i] = 0;
					for (int y = 0; y < VarMods.size(); y++) {
						auto &vmod = VarMods[y];
						for (int x = 0; x < vmod.aas.size(); x++) if (vmod.aas[x] == s[i]) {
							mod[i] |= (1 << y), mod_cnt[i]++;
							break;
						}
					}
				}
				for (cnt = 1; cnt <= m; cnt++) {
					for (i = 0; i < cnt - 1; i++) ind[i] = i;
					ind[cnt - 1] = cnt - 2;
					while (ind[0] <= l - cnt) {
						for (i = cnt - 1; i >= 0; i--) {
							if (ind[i] + cnt - i < l) {
								ind[i]++;
								if (i < cnt - 1) for (j = i + 1; j < cnt; j++) ind[j] = ind[i] + j - i;
								break;
							}
						}
						if (i < 0) break;
						for (i = 0; i < cnt; i++) if (!mod[ind[i]]) break;
						if (i < cnt) continue;
						for (i = 0; i < cnt - 1; i++) type[i] = 0;
						type[cnt - 1] = -1;
						while (type[0] < mod_cnt[ind[0]]) {
							for (i = cnt - 1; i >= 0; i--) {
								if (type[i] + 1 < mod_cnt[ind[i]]) {
									type[i]++;
									for (j = i + 1; j < cnt; j++) type[j] = 0;
									break;
								}
							}
							if (i < 0) break;
							peptide.clear();
							char lc = to_lower(s[0]);
							for (auto &fmod : FixedMods) for (int x = 0; x < fmod.aas.size(); x++) if (fmod.aas[x] == lc) { // N-terminal modification: encoded with a lower-case letter for the amino acid
								peptide += '(' + fmod.name + ')';
								break;
							}
							for (i = k = 0; i < l; i++) {
								char aa = s[i];
								peptide.push_back(aa);
								if (FixedMods.size()) for (auto &fmod : FixedMods) for (int x = 0; x < fmod.aas.size(); x++) if (fmod.aas[x] == aa) {
									peptide += '(' + fmod.name + ')';
									break;
								}
								if (i > ind[cnt - 1]) continue;
								while (k < cnt - 1 && i > ind[k]) k++;
								if (i < ind[k]) continue;
								int code = mod[i];
								for (j = 0; j < type[k]; j++) code &= code - 1;
								for (j = 0; j < 2; j++) if (code & (1 << j)) break;
								assert(j < 2);
								peptide += '(' + VarMods[j].name + ')';
							}
							peptides.push_back(peptide);
						}
					}
				}
			}
		}
		std::sort(proteins.begin(), proteins.end());
		return true;
	}

	void free() {
		std::vector<std::string>().swap(peptides);
		std::vector<FastaEntry>().swap(entries);
		std::vector<Isoform>().swap(proteins);
	}
};

std::vector<float> norm_totals, norm_shares, norm_ratios, norm_saved_ratios;
struct NormInfo {
	int index = 0;
	float signal = 0.0, par = 0.0;

	NormInfo() {}
	NormInfo(float _par, int _index, float _signal) { index = _index; signal = _signal; par = _par; }
	inline friend bool operator < (const NormInfo &left, const NormInfo &right) { return left.par < right.par; }
};
std::vector<NormInfo> norm_ind;

struct XIC { // for visualisation
	float RT = 0.0; // apex RT
	float qvalue = 0.0;
	float pr_mz = 0.0, fr_mz = 0.0; // masses used for extraction (i.e. mass corrected), not reference masses
	int run = -1, pr = -1, level = 0; // run index, precursor index in the library, MS level
	Product fr; // fragment info, if MS2 XIC
	std::vector<std::pair<float, float> > peaks; // RT, intensity pairs

	XIC() {}
	inline friend bool operator < (const XIC& left, const XIC& right) { 
		if (left.pr < right.pr) return true;
		if (left.pr == right.pr) {
			if (left.run < right.run) return true;
			if (left.run == right.run && left.fr < right.fr) return true;
		}
		return false;
	}
	inline friend bool operator == (const XIC& left, const XIC& right) { 
		return (left.pr == right.pr && left.run == right.run && !(left.fr < right.fr) && !(right.fr < left.fr)); 
	}
};

std::set<XIC> XICs;

const int fFromFasta = 1 << 0;
const int fPredictedSpectrum = 1 << 1;
const int fPredictedRT = 1 << 2;

class Library {
public:
	std::string name, fasta_names;
	std::vector<Isoform> proteins;
	std::vector<PG> protein_ids;
	std::vector<PG> protein_groups;
	std::vector<std::string> gene_groups;
	std::vector<int> gg_index;
	std::vector<std::string> precursors; // precursor IDs in canonical format; library-based entries only
	std::vector<std::string> names;
	std::vector<std::string> genes;
	int skipped = 0;
	double iRT_min = 0.0, iRT_max = 0.0;
	bool gen_decoys = true, gen_charges = true, infer_proteotypicity = true, from_speclib = false;

	std::map<std::string, int> eg;
	std::vector<int> elution_groups; // same length as entries, indices of the elution groups
	std::vector<int> co_elution;
	std::vector<std::pair<int, int> > co_elution_index;

	class Entry {
	public:
		Lock lock;
		Library * lib;
		Peptide target, decoy;
		int entry_flags = 0, proteotypic = 0;
		std::string name; // precursor id
		std::set<PG>::iterator prot;
		int pid_index = 0, pg_index = 0, best_run = -1, peak = 0, apex = 0, window = 0;
		float qvalue = 0.0, protein_qvalue = 0.0, best_fr_mz = 0.0;

		void init() {
			std::sort(target.fragments.begin(), target.fragments.end(), [](const Product &left, const Product &right) { return left.height > right.height; });
			std::sort(decoy.fragments.begin(), decoy.fragments.end(), [](const Product &left, const Product &right) { return left.height > right.height; });
			if (target.fragments.size() > MaxF) target.fragments.resize(MaxF);
			if (decoy.fragments.size() > MaxF) decoy.fragments.resize(MaxF);

			float max = 0.0;
			for (int i = 0; i < target.fragments.size(); i++) if (target.fragments[i].height > max) max = target.fragments[i].height;
			for (int i = 0; i < target.fragments.size(); i++) target.fragments[i].height /= max;

			max = 0.0;
			for (int i = 0; i < decoy.fragments.size(); i++) if (decoy.fragments[i].height > max) max = decoy.fragments[i].height;
			for (int i = 0; i < decoy.fragments.size(); i++) decoy.fragments[i].height /= max;

			target.length = decoy.length = peptide_length(name);
		}

		friend bool operator < (const Entry &left, const Entry &right) { return left.name < right.name; }

		inline void generate_decoy() {
			int i;

			decoy = target;
			auto aas = get_aas(name);
			int m = aas.size() - 2;
			if (m < 0) return;

			bool recognise = false;
			float gen_acc = 0.0;
			std::vector<Ion> pattern;
			for (auto &f : decoy.fragments) if (!(f.type & 3) || f.charge <= 0) {
				recognise = true;
				break;
			}
			if (recognise || ReverseDecoys) {
				pattern = recognise_fragments(name, target.fragments, gen_acc, target.charge, false, MaxRecCharge, MaxRecLoss);
				for (i = 0; i < pattern.size(); i++) {
					if (lib->gen_charges) target.fragments[i].charge = pattern[i].charge;
					target.fragments[i].type = (target.fragments[i].type & ~3) | pattern[i].type;
					target.fragments[i].loss = pattern[i].loss;
					target.fragments[i].index = pattern[i].index;
				}
				decoy.fragments = target.fragments;
				if (gen_acc > GeneratorAccuracy + E && ForceFragRec) {
					decoy.fragments.clear();
					lib->skipped++;
					return;
				}
			}

			if (!ReverseDecoys) {
				int pos;
				float N_shift, C_shift;
				auto mod_state = get_mod_state(name, aas.length());

				if (!mod_state[1]) pos = 1;
				else if (!mod_state[Min(2, m + 1)]) pos = Min(2, m + 1);
				else if (!mod_state[0]) pos = 0;
				else pos = 1;
				N_shift = AA[MutateAAto[AA_index[aas[pos]]]] - AA[aas[pos]];

				if (!mod_state[m]) pos = m;
				else if (!mod_state[Max(0, m - 1)]) pos = Max(0, m - 1);
				else if (!mod_state[m + 1]) pos = m + 1;
				else pos = m;
				C_shift = AA[MutateAAto[AA_index[aas[pos]]]] - AA[aas[pos]];
				
				for (i = 0; i < decoy.fragments.size(); i++) {
					double c = (double)Max(decoy.fragments[i].charge, 1);
					if (decoy.fragments[i].type & fTypeY) decoy.fragments[i].mz += C_shift / c;
					else decoy.fragments[i].mz += N_shift / c;
				}
			} else {
				auto seq = get_sequence(name);
				auto inv_seq = seq;

				for (i = 1; i < seq.size() - 1; i++) inv_seq[i] = seq[seq.size() - i - 1];
				for (i = 0; i < seq.size(); i++) if (Abs(inv_seq[i] - seq[i]) > 1200.0 * GeneratorAccuracy) break;
				if (i == seq.size()) inv_seq[1] += 12.0, inv_seq[seq.size() - 2] -= 12.0;

				auto dfr = generate_fragments(inv_seq, pattern, target.charge);
				for (i = 0; i < dfr.size(); i++) decoy.fragments[i].mz = dfr[i].mz;
			}
		}

		inline void annotate_charges() {
			int i;

			auto seq = get_sequence(name, &target.no_cal);
			if (!seq.size()) return;

			float gen_acc = 0.0;
			auto pattern = recognise_fragments(seq, target.fragments, gen_acc, target.charge, false, MaxRecCharge, MaxRecLoss);
			for (i = 0; i < pattern.size(); i++) {
				target.fragments[i].charge = pattern[i].charge;
				if (i < decoy.fragments.size()) decoy.fragments[i].charge = pattern[i].charge;
			}
		}

		std::vector<Ion> fragment() {
			assert(MaxF > 1000);

			int cnt = 0;
			auto seq = get_sequence(name);

			std::vector<Ion> gen;
			for (int charge = 1; charge <= 2; charge++) {
				for (int loss = loss_none; loss <= loss_NH3; loss++) {
					auto frs = generate_fragments(seq, charge, loss, &cnt, MinFrAAs);
					gen.insert(gen.end(), frs.begin(), frs.end());
				}
			}
			return gen;
		}

		void generate() {
			int i, cnt = 0, cnt_noloss = 0;
			auto seq = get_sequence(name);
			auto aas = get_aas(name);
			auto gen = generate_all_fragments(seq, 1, AddLosses ? loss_NH3 : loss_none, &cnt);
			std::vector<double> scores, H2O_scores, NH3_scores;
			y_scores(scores, target.charge, aas), to_exp(scores);
			if (AddLosses) y_loss_scores(H2O_scores, aas, false), y_loss_scores(NH3_scores, aas, true), to_exp(H2O_scores), to_exp(NH3_scores);
			for (i = cnt = 0; i < gen.size(); i++)
				if (gen[i].mz >= MinFrMz && gen[i].mz <= MaxFrMz) if ((gen[i].type & fTypeY) && seq.size() - gen[i].index >= MinFrAAs) cnt++;
			target.fragments.resize(cnt);
			for (i = cnt = 0; i < gen.size(); i++)
				if (gen[i].mz >= MinFrMz && gen[i].mz <= MaxFrMz) if ((gen[i].type & fTypeY) && seq.size() - gen[i].index >= MinFrAAs) {
					target.fragments[cnt].mz = gen[i].mz;
					if (gen[i].loss == loss_none) target.fragments[cnt].height = scores[gen[i].index];
					else if (gen[i].loss == loss_H2O) target.fragments[cnt].height = scores[gen[i].index] * H2O_scores[gen[i].index];
					else target.fragments[cnt].height = scores[gen[i].index] * NH3_scores[gen[i].index];
					target.fragments[cnt].charge = gen[i].charge;
					target.fragments[cnt].type = gen[i].type;
					target.fragments[cnt].index = gen[i].index;
					target.fragments[cnt].loss = gen[i].loss;
					if (gen[i].loss == loss_none) cnt_noloss++;
					cnt++;
				}
			generate_decoy();
			init();
		}

		inline void free() {
			target.free();
			decoy.free();
		}

		template <class F> void write(F &out) {
			target.write(out);
			int dc = !lib->gen_decoys;
			out.write((char*)&dc, sizeof(int));
			if (dc) decoy.write(out);

			int ff = entry_flags, prt = proteotypic;
			out.write((char*)&ff, sizeof(int));
			out.write((char*)&prt, sizeof(int));
			out.write((char*)&pid_index, sizeof(int));
			write_string(out, name);
		}

		template <class F> void read(F &in) {
			target.read(in);
			int dc = 0; in.read((char*)&dc, sizeof(int));
			if (dc) decoy.read(in);

			int ff = 0, prt = 0;
			in.read((char*)&ff, sizeof(int));
			in.read((char*)&prt, sizeof(int));
			entry_flags = ff, proteotypic = prt;
			in.read((char*)&pid_index, sizeof(int));
			read_string(in, name);
		}

#if (HASH > 0)
		unsigned int hash() { return target.hash() ^ decoy.hash(); }
#endif
	};

	std::vector<Entry> entries;

	template<class F> void write(F &out) {
		int gd = gen_decoys, gc = gen_charges, ip = infer_proteotypicity;
		out.write((char*)&gd, sizeof(int));
		out.write((char*)&gc, sizeof(int));
		out.write((char*)&ip, sizeof(int));

		write_string(out, name);
		write_string(out, fasta_names);
		write_array(out, proteins);
		write_array(out, protein_ids);
		write_strings(out, precursors);
		write_strings(out, names);
		write_strings(out, genes);
		out.write((char*)&iRT_min, sizeof(double));
		out.write((char*)&iRT_max, sizeof(double));
		write_array(out, entries);
		write_vector(out, elution_groups);
	}

	template<class F> void read(F &in) {
		int gd = 0, gc = 0, ip = 0;
		in.read((char*)&gd, sizeof(int));
		in.read((char*)&gc, sizeof(int));
		in.read((char*)&ip, sizeof(int));
		gen_decoys = gd, gen_charges = gc, infer_proteotypicity = ip;

		read_string(in, name);
		read_string(in, fasta_names);
		read_array(in, proteins);
		read_array(in, protein_ids);
		read_strings(in, precursors);
		read_strings(in, names);
		read_strings(in, genes);
		in.read((char*)&iRT_min, sizeof(double));
		in.read((char*)&iRT_max, sizeof(double));
		read_array(in, entries);
		for (auto &e : entries) e.lib = this;
		if (in.peek() != std::char_traits<char>::eof()) read_vector(in, elution_groups);
	}

	void generate_spectra() {
		for (int i = 0; i < entries.size(); i++)
			if (entries[i].lock.set()) entries[i].generate();
	}
	
	void generate_all() {
		int i;
		if (Threads > 1) {
			std::vector<std::thread> threads;
			for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Library::generate_spectra, this));
			for (i = 0; i < Threads; i++) threads[i].join();
		} else generate_spectra();
		for (i = 0; i < entries.size(); i++) entries[i].lock.free();
	}

#ifdef PREDICTOR
	void smp_generate_codes(std::vector<std::pair<std::vector<long long>, int> > *_codes, std::vector<std::pair<std::string, int> > *_dict, std::vector<Lock> *_locks) {
		std::string temp_s;
		std::vector<std::string> temp_sv;
		auto &codes = *_codes;
		auto &dict = *_dict;
		auto &locks = *_locks;

		int skipped = 0;
		for (int i = 0; i < entries.size(); i++) if (locks[i].set()) {
			auto &e = entries[i];
			predictor::code_from_precursor(codes[i].first, to_charged_eg(e.name, e.target.charge), temp_s, temp_sv, dict);
			codes[i].second = i;
		}
	}

	void smp_decode_spectra(std::vector<std::vector<float> > *_spectra, std::vector<int> *_decoder, int max_frc, std::vector<Lock> *_locks) {
		auto &spectra = *_spectra;
		auto &decoder = *_decoder;
		auto &locks = *_locks;
		int i, j, cnt;

		std::vector<Ion> ions, filtered;
		for (i = 0; i < entries.size(); i++) if (locks[i].set()) {
			int ind = decoder[i];
			if (ind < 0) continue;

			auto &e = entries[i];
			auto &sp = spectra[ind];

			int tot = sp.size(), L = tot / (2 * max_frc);
			if (tot < MinGenFrNum) continue;

			ions.resize(tot); memset(&(ions[0]), 0, tot * sizeof(Ion));
			for (int c = 1; c <= max_frc; c++) for (j = 0; j < L; j++) {
				int index = L * (c - 1) + j;
				ions[index] = Ion(fTypeY, j + 1, loss_none, c, sp[index]); // trick with writing intensity in the m/z field of class Ion for later sorting
				index += L * max_frc;
				ions[index] = Ion(fTypeB, j + 3, loss_none, c, sp[index]);
			}

			std::sort(ions.begin(), ions.end(), [&](const Ion &x, const Ion &y) { return x.mz > y.mz; });
			if (ions[MinGenFrNum - 1].mz - E <= E) continue;
			double margin = Min(ions[0].mz * 0.001, ions[MinGenFrNum - 1].mz - E);
			filtered.clear(); for (auto &ion : ions) if (ion.mz > margin && (ion.type == fTypeY ? (L + 3 - ion.index) : ion.index) >= MinFrAAs) filtered.push_back(ion);
			auto frs = generate_fragments(get_sequence(e.name), filtered, e.target.charge);
			for (j = cnt = 0; j < frs.size(); j++) if (frs[j].mz >= MinFrMz && frs[j].mz <= MaxFrMz) cnt++;
			e.target.fragments.resize(Min(cnt, auxF));
			for (j = cnt = 0; j < frs.size(); j++) if (frs[j].mz >= MinFrMz && frs[j].mz <= MaxFrMz) {
				e.target.fragments[cnt++] = Product(frs[j].mz, filtered[j].mz, filtered[j].charge, filtered[j].type, filtered[j].index, loss_none);
				if (cnt >= auxF) break;
			}
			e.entry_flags |= fPredictedSpectrum;

			e.generate_decoy();
			e.init();
		}
	}

	void generate_predictions() {
		if (!entries.size()) return;
		predictor::init_predictor(1, Threads);
		P.set_instance(0);
		int i, j;

		if (Verbose >= 1) Time(), std::cout << "Encoding peptides for spectra and RTs prediction\n";
		auto dict_map = P.get_aa_indices();
		std::vector<std::pair<std::string, int> > dict(dict_map.begin(), dict_map.end());
		std::vector<std::pair<std::vector<long long>, int> > codes(entries.size());
		std::vector<Lock> locks(entries.size());

		{
			std::vector<std::thread> threads;
			for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Library::smp_generate_codes, this, &codes, &dict, &locks));
			for (i = 0; i < Threads; i++) threads[i].join();
			for (auto &L : locks) L.free();
		}

		int skipped = 0;
		for (auto &c : codes) if (!c.first.size()) skipped++;
		if (skipped) cout << "WARNING: skipping " << skipped << " precursors; unrecognised modifications?\n";

		auto rt_codes = codes;
		std::sort(codes.begin(), codes.end(), [](auto &x, auto &y) { return x.first < y.first; });
		for (auto &c : rt_codes) if (c.first.size()) c.first[0] = 0;
		std::sort(rt_codes.begin(), rt_codes.end(), [](auto &x, auto &y) { return x.first < y.first; });
		std::vector<int> fr_decoder(entries.size(), -1), rt_decoder(entries.size(), -1);
		int last_fr = 0, last_rt = 0;
		auto *code = &codes[0], *rt_code = &rt_codes[0];
		for (i = 0; i < codes.size(); i++) {
			auto *curr = &codes[i];
			bool flag = (i == 0);
			if (!curr->first.size()) goto rt;
			if (!flag) if (curr->first != code->first) flag = true;
			if (!flag) curr->first.clear();
			else last_fr = i, code = curr;
			fr_decoder[curr->second] = last_fr;

		rt:
			curr = &rt_codes[i];
			if (!curr->first.size()) continue;
			flag = (i == 0);
			if (!flag) if (curr->first != rt_code->first) flag = true;
			if (!flag) curr->first.clear();
			else last_rt = i, rt_code = curr;
			rt_decoder[curr->second] = last_rt;
		}

		std::vector<std::vector<float> > spectra, rts;
		try {
			if (Verbose >= 1) Time(), std::cout << "Predicting spectra\n";
			P.predict(spectra, codes, 1, Verbose >= 3);
			if (Verbose >= 1) Time(), std::cout << "Predicting RTs\n";
			P.predict(rts, rt_codes, 2, Verbose >= 3);
		} catch (std::exception &e) { std::cout << "ERROR: " << e.what() << '\n'; }

		if (Verbose >= 1) Time(), std::cout << "Decoding predicted spectra\n";
		{
			std::vector<std::thread> threads;
			for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Library::smp_decode_spectra, this, &spectra, &fr_decoder, 3, &locks));
			for (i = 0; i < Threads; i++) threads[i].join();
		}

		if (Verbose >= 1) Time(), std::cout << "Decoding RTs\n";
		for (i = 0; i < entries.size(); i++) {
			int ind = rt_decoder[i];
			if (ind >= 0) if (rts[ind].size()) {
				entries[i].target.iRT = entries[i].target.sRT = rts[ind][0] * 243.19 - 58.4;
				entries[i].entry_flags |= fPredictedRT;
			}
		}
	}
#endif

	class Info {
	public:
		Library * lib;
		int n_s; // number of samples (runs)
		int n_entries; // total number of entries
		std::map<int, std::vector<std::pair<int, QuantEntry> > > map;
		std::vector<std::vector<int> > proteins;
		std::vector<RunStats> RS;

		Info() { }

		void clear() { map.clear(); std::vector<std::vector<int> >().swap(proteins); }

		void load(Library * parent, std::vector<std::string> &files, std::vector<Quant> * _quants = NULL) {
			if (Verbose >= 1) Time(), dsout << "Reading quantification information: " << (n_s = files.size()) << " files\n";
			clear();

			lib = parent;
			n_entries = 0;
			proteins.clear();
			proteins.resize(n_s);
			RS.resize(n_s);
			for (int i = 0; i < n_s; i++) {
				Quant Qt;
				auto Q = &Qt;
				if (_quants == NULL || !QuantInMem) {
					std::ifstream in(std::string(files[i]) + std::string(".quant"), std::ifstream::binary);
					if (in.fail() || temp_folder.size()) in = std::ifstream(location_to_file_name(files[i]) + std::string(".quant"), std::ifstream::binary);
					Qt.read(in, lib->entries.size());
					in.close();
				} else Q = &((*_quants)[i]);
				memcpy(&(RS[i]), &(Q->RS), sizeof(RunStats));
				proteins[i] = Q->proteins;
				for (auto it = Q->entries.begin(); it != Q->entries.end(); it++) {
					n_entries++;
					int index = it->pr.index;
					auto pos = map.find(index);
					if (pos != map.end())
						pos->second.push_back(std::pair<int, QuantEntry>(i, *it));
					else {
						std::vector<std::pair<int, QuantEntry> > v;
						v.push_back(std::pair<int, QuantEntry>(i, *it));
						map.insert(std::pair<int, std::vector<std::pair<int, QuantEntry> > >(index, v));
					}
				}
				if (_quants == NULL || !QuantInMem) std::vector<QuantEntry>().swap(Q->entries);
			}
		}

		void load(Library * parent, Quant &quant) {
			clear();

			lib = parent;
			n_entries = 0;
			proteins.clear();
			proteins.resize(n_s = 1);

			Quant &Q = quant;
			proteins[0] = Q.proteins;
			for (auto it = Q.entries.begin(); it != Q.entries.end(); it++) {
				n_entries++;
				int index = it->pr.index;
				auto pos = map.find(index);
				if (pos != map.end())
					pos->second.push_back(std::pair<int, QuantEntry>(0, *it));
				else {
					std::vector<std::pair<int, QuantEntry> > v;
					v.push_back(std::pair<int, QuantEntry>(0, *it));
					map.insert(std::pair<int, std::vector<std::pair<int, QuantEntry> > >(index, v));
				}
			}
		}

		void quantify() {
			int i, j, k, m, pos;

			if (Verbose >= 1) Time(), dsout << "Quantifying peptides\n";
			std::vector<std::pair<float, int> > fr_score(1000), fr_ordered(1000);
			std::vector<bool> exclude(1024);
			for (auto it = map.begin(); it != map.end(); it++) {
				auto v = &(it->second);
				int pr_index = it->first;
				auto &lib_e = lib->entries[pr_index].target;

				fr_score.clear();
				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) if (jt->second.pr.qvalue <= MaxQuantQvalue) {
					for (i = 0; i < jt->second.fr_n; i++) {
						int ind = jt->second.fr[i].index;
						if (ind >= fr_score.size()) fr_score.resize(ind + 1, std::pair<float, int>(0.0, 0));
						fr_score[ind].first += jt->second.fr[i].corr;
						fr_score[ind].second++;
					}
				}
				for (auto &s : fr_score) if (s.second) s.first /= (double)s.second;

				if (ExcludeSharedFragments) { // exclude (from quantification) fragments shared by peptides in the same elution group
					int eg = lib->elution_groups[pr_index];
					auto &ce = lib->co_elution_index[eg];
					for (int pos = ce.first; pos < ce.first + ce.second; pos++) {
						int next = lib->co_elution[pos];
						if (next != pr_index && lib->entries[next].target.charge == lib->entries[pr_index].target.charge) { // difference only in the label
							auto &lf = lib_e;
							auto &ls = lib->entries[next].target;
							for (i = 0; i < lf.fragments.size(); i++) {
								double margin = GeneratorAccuracy * lf.fragments[i].mz;
								bool flag = false;
								for (j = 0; j < ls.fragments.size(); j++)
									if (Abs(lf.fragments[i].mz - ls.fragments[j].mz) < margin) {
										flag = true;
										break;
									}
								if (flag) fr_score[i].first -= 2.0;
							}
						}
					}
				}

				if (RestrictFragments) { // exclude fragments from quantification based on their library annotation
					auto &frs = lib_e.fragments;
					assert(fr_score.size() <= frs.size());
					for (j = 0; j < fr_score.size(); j++) if (frs[j].type & fExclude) fr_score[j].first -= INF;
				}

				fr_ordered.resize(fr_score.size());
				fr_ordered.assign(fr_score.begin(), fr_score.end());
				std::sort(fr_ordered.begin(), fr_ordered.end());
				for (pos = fr_ordered.size() - 3; pos < fr_ordered.size() - 1; pos++) if (fr_ordered[pos].first > -1.0 + E) break;
				double margin = Max(-INF / 2.0, fr_ordered[pos].first - E);
				if (NoFragmentSelectionForQuant) margin = -INF / 2.0;
				if (GenFrExclusionInfo) for (j = 0; j < fr_score.size(); j++) if (fr_score[j].second) {
					auto &fr = lib_e.fragments[j];
					if (fr_score[j].first >= margin) fr.type &= ~fExclude;
					else fr.type |= fExclude;
				}

				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					for (i = 0, jt->second.pr.quantity = jt->second.pr.ratio = 0.0; i < jt->second.fr_n; i++) {
						int ind = jt->second.fr[i].index;
						if (fr_score[ind].first >= margin) jt->second.pr.quantity += jt->second.fr[i].quantity[qFiltered];
					}
					jt->second.pr.norm = jt->second.pr.level = jt->second.pr.quantity;
				}
			}

			if (TranslatePeaks) { // calculate ratios between labelled and unlabelled peptides; requires TranslatePeaks to ensure co-elution
				for (auto it = map.begin(); it != map.end(); it++) {
					auto v = &(it->second);
					int pr_index = it->first;
					auto &lib_e = lib->entries[pr_index].target;

					if (lib->entries[pr_index].entry_flags & fFromFasta) continue; // only spectral library entries

					int eg = lib->elution_groups[pr_index];
					auto &ce = lib->co_elution_index[eg];

					exclude.clear();
					if (ExcludeSharedFragments) { // exclude (from quantification) fragments shared by peptides in the same elution group
						int eg = lib->elution_groups[pr_index];
						auto &ce = lib->co_elution_index[eg];
						for (int pos = ce.first; pos < ce.first + ce.second; pos++) {
							int next = lib->co_elution[pos];
							if (next != pr_index && lib->entries[next].target.charge == lib->entries[pr_index].target.charge) { // difference only in the label
								auto &lf = lib_e;
								auto &ls = lib->entries[next].target;
								for (i = 0; i < lf.fragments.size(); i++) {
									double margin = GeneratorAccuracy * lf.fragments[i].mz;
									bool flag = false;
									for (j = 0; j < ls.fragments.size(); j++)
										if (Abs(lf.fragments[i].mz - ls.fragments[j].mz) < margin) {
											flag = true;
											break;
										}
									if (flag) {
										if (i >= exclude.size()) exclude.resize(i + 1, false);
										exclude[i] = true;
									}
								}
							}
						}
					}

					const int fr_n = 3;
					std::vector<float> fratio(fr_n);
					for (int pos = ce.first; pos < ce.first + ce.second; pos++) {
						int next = lib->co_elution[pos];
						if (next != pr_index) if (lib->entries[next].target.charge == lib->entries[pr_index].target.charge) {
							auto p = map.find(next);
							if (p != map.end()) if (p->first == next) {
								auto nv = &(p->second);
								auto kt = (*nv).begin();
								for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
									for (; kt != (*nv).end(); kt++) if (kt->first >= jt->first) break;
									if (kt->first == jt->first) {
										fr_score.clear();
										for (i = 0; i < jt->second.fr_n; i++) {
											int ind = jt->second.fr[i].index;
											if (ind >= fr_score.size()) fr_score.resize(ind + 1, std::pair<float, int>(0.0, 0));
											fr_score[ind].first = jt->second.fr[i].corr, fr_score[ind].second = lib->entries[pr_index].target.fragments[ind].ion_code();
										}
										for (i = 0; i < kt->second.fr_n; i++) {
											int ind = kt->second.fr[i].index;
											if (ind >= fr_score.size()) fr_score.resize(ind + 1, std::pair<float, int>(0.0, 0));
											if (fr_score[ind].second != lib->entries[next].target.fragments[ind].ion_code()) fr_score[i].first = -INF;
											else if (kt->second.fr[i].corr < fr_score[ind].first) fr_score[ind].first = kt->second.fr[i].corr;
										}
										for (i = 0; i < fr_score.size(); i++) {
											fr_score[i].second = i;
											if (i < exclude.size()) if (exclude[i]) fr_score[i].first -= INF * 0.5;
										}
										std::sort(fr_score.begin(), fr_score.end());
										if (fr_score.size()) if (fr_score[fr_score.size() - 1].first > E) {
											int fcnt = 0;
											float f1 = 0.0, f2 = 0.0, f1f[fr_n], f2f[fr_n];
											for (i = 0; i < fr_n; i++) f1f[i] = f2f[i] = 0.0;
											for (int ifr = fr_score.size() - 1; ifr >= 0 && ifr >= fr_score.size() - fr_n; ifr--) {
												if (ifr < fr_score.size() - 1) if (fr_score[ifr].first < 0.8) break;
												int fr = fr_score[ifr].second;
												for (i = 0; i < jt->second.fr_n; i++) if (jt->second.fr[i].index == fr) {
													f1f[fcnt] += jt->second.fr[i].quantity[qFiltered];
													break;
												}
												for (i = 0; i < kt->second.fr_n; i++) if (kt->second.fr[i].index == fr) {
													f2f[fcnt] += kt->second.fr[i].quantity[qFiltered];
													break;
												}
												fcnt++;
												if (fcnt >= fr_n) break;
											}
											fratio.resize(fr_n);
											for (i = fcnt = 0; i < fr_n; i++) if (f1f[i] > E) fratio[fcnt++] = f2f[i] / f1f[i];
											if (fcnt <= 2) {
												for (i = 0; i < fr_n; i++) f1 += f1f[i], f2 += f2f[i];
												if (f1 > E) jt->second.pr.ratio += (jt->second.pr.quantity * f2) / f1;
											} else {
												fratio.resize(fcnt); std::sort(fratio.begin(), fratio.end());
												jt->second.pr.ratio += jt->second.pr.quantity * fratio[fcnt / 2];
											}
											continue;
										}
										jt->second.pr.ratio += kt->second.pr.quantity;
									}
								}
							}
						}
					}
					for (auto jt = (*v).begin(); jt != (*v).end(); jt++) jt->second.pr.ratio = (jt->second.pr.ratio > E ? (jt->second.pr.quantity / jt->second.pr.ratio) : -1.0);
				}
			}

			if (ms_files.size() <= 1 || IndividualReports) return;

			// simple normalisation by the total signal first
			double av = 0.0;
			std::vector<float> sums(ms_files.size());
			for (auto it = map.begin(); it != map.end(); it++, i++) {
				auto v = &(it->second);
				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					int index = jt->first;
					if (jt->second.pr.qvalue <= NormalisationQvalue) sums[index] += jt->second.pr.quantity;
				}
			}
			for (i = k = 0; i < sums.size(); i++) if (sums[i] > E) av += sums[i], k++;
			if (k) av /= (double)k;
			if (av <= E) {
				Warning("not enough peptides for normalisation");
				return;
			}

			for (i = 0; i < sums.size(); i++) if (sums[i] <= E) sums[i] = av;
			for (auto it = map.begin(); it != map.end(); it++) {
				auto v = &(it->second);
				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					int index = jt->first;
					jt->second.pr.level = (jt->second.pr.quantity * av) / sums[index];
				}
			}

			// advanced normalisation
			std::vector<float> score(map.size());
			std::vector<float> x(ms_files.size());
			i = m = 0;
			for (auto it = map.begin(); it != map.end(); it++, i++) {
				auto v = &(it->second);

				k = 0;
				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					int index = jt->first;
					if (jt->second.pr.qvalue <= NormalisationQvalue) x[index] = jt->second.pr.level, k++, m++;
					else x[index] = 0.0;
				}
				if (k >= 2) {
					double u = mean(&(x[0]), x.size());
					if (u > E) score[i] = sqrt(var(&(x[0]), x.size())) / u;
					else score[i] = INF;
				}
				else score[i] = INF;
			}
			if (!m) {
				Warning("not enough peptides for normalisation");
				return;
			}

			auto ordered = score;
			std::sort(ordered.begin(), ordered.end());
			int average_number = m / ms_files.size();
			int used = Max(1, int(NormalisationPeptidesFraction * (double)average_number));
			double margin = ordered[used] + E;

			for (i = 0; i < sums.size(); i++) sums[i] = 0.0;
			i = 0;
			for (auto it = map.begin(); it != map.end(); it++, i++) {
				auto v = &(it->second);
				if (score[i] > margin) continue;

				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					int index = jt->first;
					if (jt->second.pr.qvalue <= NormalisationQvalue) sums[index] += jt->second.pr.level;
				}
			}

			for (i = k = 0, av = 0.0; i < sums.size(); i++) if (sums[i] > E) av += sums[i], k++;
			if (k) av /= (double)k;
			if (av < E) {
				Warning("cannot perform normalisation");
				return;
			}
			for (i = 0; i < sums.size(); i++) if (sums[i] <= E) sums[i] = av;
			for (auto it = map.begin(); it != map.end(); it++) {
				auto v = &(it->second);

				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					int index = jt->first;
					jt->second.pr.norm = (jt->second.pr.level * av) / sums[index];
				}
			}

			if (LocalNormalisation) for (int iter = 0; iter < 2; iter++) { // iter = 0: RT-dependent, iter = 1: signal-dependent normalisation
				if (iter == 0 && NoRTDepNorm) continue;
				if (iter == 1 && NoSigDepNorm) continue;

				norm_totals.clear(), norm_shares.clear();
				norm_totals.resize(map.size(), 0.0), norm_shares.resize(map.size(), 0.0);
				i = 0;
				for (auto it = map.begin(); it != map.end(); it++, i++) {
					auto v = &(it->second);
					int cnt = 0;
					for (auto jt = (*v).begin(); jt != (*v).end(); jt++)
						if (jt->second.pr.qvalue <= NormalisationQvalue) norm_totals[i] += jt->second.pr.norm, cnt++;
					if (cnt) norm_shares[i] = 1.0 / (double)cnt;
				}

				NormInfo nie;

				norm_ind.resize(map.size()), norm_saved_ratios.resize(map.size());
				norm_ratios.resize(2 * LocNormRadius + 1);
				for (k = 0; k < ms_files.size(); k++) {
					i = 0;
					norm_ind.clear(), norm_saved_ratios.clear();
					for (auto it = map.begin(); it != map.end(); it++, i++) {
						auto v = &(it->second);
						for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
							int index = jt->first;
							if (index == k && jt->second.pr.qvalue <= NormalisationQvalue && norm_shares[i] > E) if (norm_shares[i] * (double)ms_files.size() <= 2.0)
								norm_ind.push_back(NormInfo(iter == 0 ? jt->second.pr.RT : jt->second.pr.norm, i, jt->second.pr.norm));
						}
					}
					std::sort(norm_ind.begin(), norm_ind.end());
					norm_saved_ratios.resize(norm_ind.size(), 0.0);

					for (auto it = map.begin(); it != map.end(); it++, i++) {
						auto v = &(it->second);
						for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
							int index = jt->first;
							if (index == k) {
								auto pos = std::lower_bound(norm_ind.begin(), norm_ind.end(), NormInfo(iter == 0 ? jt->second.pr.RT : jt->second.pr.norm, 0, 0.0));
								float run = 0.0, all = 0.0, ratio = 0.0, w = 0.0;
								int ind = std::distance(norm_ind.begin(), pos), low = Max(0, ind - LocNormRadius), high = low + 2 * LocNormRadius;
								if (high > norm_ind.size()) high = norm_ind.size(), low = Max(0, high - 2 * LocNormRadius);
								if (norm_saved_ratios[Min(ind, norm_ind.size() - 1)] > E) ratio = norm_saved_ratios[Min(ind, norm_ind.size() - 1)];
								else {
									norm_ratios.clear();
									for (j = low; j < high; j++) {
										float tot = norm_totals[norm_ind[j].index] * norm_shares[norm_ind[j].index];
										float r = tot / Max(E, norm_ind[j].signal);
										run += norm_ind[j].signal, all += tot, ratio += r, w += 1.0;
										norm_ratios.push_back(r);
									}
									std::sort(norm_ratios.begin(), norm_ratios.end());
									ratio = 0.5 * (norm_ratios[norm_ratios.size() / 2] + norm_ratios[(norm_ratios.size() / 2) + 1]);
									norm_saved_ratios[Min(ind, norm_ind.size() - 1)] = ratio = Min(LocNormMax, Max(1.0 / LocNormMax, ratio));
								}
								jt->second.pr.norm *= ratio;
							}
						}
					}
				}
			}
		}
	};

	Info info;

	void extract_proteins() {
		std::set<Isoform> prot;
		std::string word;

		for (auto &pg : protein_ids) {
			std::stringstream list(pg.ids);
			while (std::getline(list, word, ';')) {
				auto ins = prot.insert(Isoform(fast_trim(word)));
				for (auto &pr : pg.precursors) ins.first->precursors.insert(pr);
			}
		}

		proteins.clear();
		proteins.insert(proteins.begin(), prot.begin(), prot.end());
		prot.clear();

		for (int i = 0; i < proteins.size(); i++) {
			auto &p = proteins[i];
			for (auto &pr : p.precursors) {
				auto &pid = protein_ids[entries[pr].pid_index];
				pid.proteins.insert(i);
			}
		}

		if (infer_proteotypicity) {
			if (Verbose >= 1 && (!FastaSearch || GuideLibrary)) Time(), dsout << "Finding proteotypic peptides (assuming that the list of UniProt ids provided for each peptide is complete)\n";
			for (auto &e : entries) if (protein_ids[e.pid_index].proteins.size() <= 1) e.proteotypic = true;
		}
	}

	void elution_group_index() {
		int i, egt = 0, peg;
		for (i = 0; i < elution_groups.size(); i++) if (elution_groups[i] > egt) egt = elution_groups[i]; egt++;
		co_elution_index.resize(egt, std::pair<int, int>(0, 1));
		std::set<std::pair<int, int> > ce;
		for (i = 0; i < elution_groups.size(); i++) ce.insert(std::pair<int, int>(elution_groups[i], i));
		co_elution.resize(ce.size()); i = peg = 0;
		if (co_elution_index.size()) co_elution_index[0].second = 0;
		for (auto it = ce.begin(); it != ce.end(); it++, i++) {
			co_elution[i] = it->second;
			if (it->first != peg) co_elution_index[it->first].first = i;
			else co_elution_index[it->first].second++;
			peg = it->first;
		}
		ce.clear();
	}

	bool load(const char * file_name) {
		int colInd[libCols];

		name = std::string(file_name);
		if (Verbose >= 1) Time(), dsout << "Loading spectral library " << name << "\n";

		elution_groups.clear();
		if (get_extension(name) == std::string(".speclib")) {
			std::ifstream speclib(file_name, std::ifstream::binary);
			if (speclib.fail()) {
				dsout << "cannot read the file\n";
				return false;
			}
			read(speclib);
			speclib.close();

			if (!elution_groups.size()) { // no elution groups have been added: old version of the .speclib format
				for (auto &e : entries) {
					auto name = to_eg(e.name);
					auto egp = eg.insert(std::pair<std::string, int>(name, eg.size()));
					elution_groups.push_back(egp.first->second);
				}
			}

			from_speclib = true;
			if (PGLevel != PGLevelSet) {
				if (PGLevelSet == 2 && genes.size() >= 2) PGLevel = 2;
				else if (PGLevelSet == 1 && names.size() >= 2) PGLevel = 1;
			}
			if (fasta_names.size()) if (Verbose >= 1) Time(), dsout << "Library annotated with sequence database(s): " << fasta_names << "\n";
			if (!fasta_files.size() && (genes.size() >= 2 || names.size() >= 2)) library_protein_stats();
		} else {
			std::ifstream csv(file_name, std::ifstream::in);
			if (csv.fail()) {
				dsout << "cannot read the file\n";
				return false;
			}

			int i, cnt = 0;
			bool ftw = false, ftmw = false, flw = false, fragment_num_info = false, fragment_loss_info = false, fragment_type = false;
			std::map<std::string, Entry> map;
			Entry e; e.lib = this;

			std::string line, prev_id = "", id;
			std::vector<std::string> words(1000);
			char delim = '\t';
			if (std::string(file_name).find(".csv") != std::string::npos) delim = ','; // not for OpenSWATH libraries
			auto ins = map.insert(std::pair<std::string, Entry>("", e)).first; // insert a dummy entry to create a variable (ins) of the desired type
			map.clear(); // remove the dummy entry

			std::set<PG> prot;
			while (std::getline(csv, line)) {
				if (!cnt) {
					if (delim == '\t' && line.find("\t") == std::string::npos) delim = ',';
					if (delim == ',' && line.find(",") == std::string::npos) delim = '\t';
				}
				int in_quote = 0, nw = 0, start = 0, end = 0;
				for (; end < line.size(); end++) {
					if (line[end] == delim && !in_quote) {
						int en = end, st = start;
						if (line[st] == '\"') st++, en--;
						if (en > st) {
							if (nw > words.size()) words.resize(words.size() * 2 + 1);
							words[nw].resize(en - st);
							for (i = st; i < en; i++) words[nw][i - st] = line[i];
						} else words[nw].clear();
						start = end + 1, nw++;
					} else if (line[end] == '\"') in_quote ^= 1;
				}
				if (line[start] == '\"') start++, end--;
				if (end > start) {
					if (nw > words.size()) words.resize(words.size() * 2 + 1);
					words[nw].resize(end - start);
					for (i = start; i < end; i++) words[nw][i - start] = line[i];
					nw++;
				}
				cnt++;

				if (cnt == 1) { // header
					std::vector<std::string>::iterator it, loc;
					for (i = 0; i < libCols; i++) colInd[i] = -1;
					for (i = 0; i < nw; i++) {
						for (auto jt = library_headers.begin(); jt != library_headers.end(); jt++)
							if (jt->find(words[i]) != std::string::npos) {
								int ind = std::distance(library_headers.begin(), jt);
								if (colInd[ind] < 0) colInd[ind] = i;
								else if (words[colInd[ind]].size() < words[i].size()) colInd[ind] = i;
								break;
							}
					}
					for (i = 0; i < libCols; i++) if (colInd[i] < 0 && i < libPID) {
						if (Verbose >= 1) dsout << "WARNING: cannot find column " + library_headers[i] << "\n";
						csv.close();
						return false;
					}
					if (colInd[libPT] >= 0) infer_proteotypicity = false;
					if (colInd[libFrCharge] >= 0) {
						gen_charges = false;
						if (colInd[libFrType] >= 0) {
							fragment_type = true;
							if (colInd[libFrNumber] >= 0) {
								fragment_num_info = true;
								if (colInd[libFrLoss] >= 0) fragment_loss_info = true;
								else dsout << "WARNING: no neutral loss information found in the library - assuming fragments without losses\n";
							} else dsout << "WARNING: no fragment number information found in the library\n";
						} else dsout << "WARNING: no fragment type information found in the library\n";
					} else dsout << "WARNING: no fragment charge information found in the library - assuming fragments with charge 1\n";
					continue;
				}

				Product p(std::stof(words[colInd[libFrMz]]), std::stof(words[colInd[libFrI]]), (!gen_charges ? std::stof(words[colInd[libFrCharge]]) : 1));
				if (fragment_type) {
					auto &wt = words[colInd[libFrType]];
					if (wt.size()) {
						if (wt[0] == 'y') p.type = fTypeY;
						else if (wt[0] == 'b') p.type = fTypeB;
						else if (!ftw) ftw = true, p.type = 0, dsout << "WARNING: unknown fragment type " << wt << "\n";
					} else if (!ftmw) ftmw = true, p.type = 0, dsout << "WARNING: fragment type missing for row number " << cnt << "\n";
					if (p.type && fragment_num_info) {
						p.index = std::stoi(words[colInd[libFrNumber]]);
						if (p.type != fTypeB) p.index = peptide_length(words[colInd[libPr]]) - p.index;
						if (fragment_loss_info) {
							auto &wl = words[colInd[libFrLoss]];
							if (wl == "noloss") p.loss = loss_none;
							else if (wl == "H2O") p.loss = loss_H2O;
							else if (wl == "NH3") p.loss = loss_NH3;
							else if (wl == "CO") p.loss = loss_CO;
							else p.loss = loss_other;
						} else p.loss = loss_none;
					}
				}
				if (colInd[libFrExc] >= 0) {
					auto &we = words[colInd[libFrExc]];
					if (we.size()) if (we[0] == 'T' || we[0] == '1' || we[0] == 't') p.type |= fExclude;
				}

				bool decoy_fragment = false;
				if (colInd[libIsDecoy] >= 0) {
					auto &w = words[colInd[libIsDecoy]];
					if (w.length()) if (w[0] == '1' || w[0] == 't' || w[0] == 'T') decoy_fragment = true;
				}
				if (decoy_fragment && FastaSearch) continue; // do not use supplied decoys for the guide library
				if (decoy_fragment) gen_decoys = false;
				int charge = std::stoi(words[colInd[libCharge]]);
				id = to_canonical(words[colInd[libPr]], charge);
				if (id != prev_id || !prev_id.length()) {
					ins = map.insert(std::pair<std::string, Entry>(id, e)).first;
					prev_id = id;
					ins->second.name = ins->first;

					auto eg_name = (colInd[libEG] >= 0 ? words[colInd[libEG]] : to_eg(words[colInd[libPr]]));
					auto egp = eg.insert(std::pair<std::string, int>(eg_name, eg.size()));
					elution_groups.push_back(egp.first->second);

					auto prot_id = (colInd[libPID] >= 0 ? words[colInd[libPID]] : "");
					ins->second.prot = prot.insert(PG(prot_id)).first;
					if (colInd[libPN] >= 0) ins->second.prot->names = words[colInd[libPN]];
					if (colInd[libGenes] >= 0) ins->second.prot->genes = words[colInd[libGenes]];
					if (colInd[libPT] >= 0 && words[colInd[libPT]].size()) {
						char pt = words[colInd[libPT]][0];
						if (pt == 'T' || pt == 't' || pt == '1' || pt == 'Y' || pt == 'y') ins->second.proteotypic = true;
						else ins->second.proteotypic = false;
					}
				}

				auto pep = !decoy_fragment ? (&(ins->second.target)) : (&(ins->second.decoy));
				pep->mz = std::stof(words[colInd[libPrMz]]);
				pep->iRT = std::stof(words[colInd[libiRT]]);
				if (colInd[libQ] >= 0) pep->lib_qvalue = std::stof(words[colInd[libQ]]);
				pep->charge = charge;
				pep->fragments.push_back(p);
			}
			csv.close();

			eg.clear();
			entries.resize(map.size());
			precursors.resize(map.size());
			i = 0;
			for (auto it = map.begin(); it != map.end(); it++, i++) {
				precursors[i] = it->first;
				entries[i] = it->second;
				entries[i].target.index = entries[i].decoy.index = i;
				if (entries[i].target.iRT < iRT_min) iRT_min = entries[i].target.iRT;
				if (entries[i].target.iRT > iRT_max) iRT_max = entries[i].target.iRT;
				it->second.prot->precursors.push_back(i);
			}
			protein_ids.insert(protein_ids.begin(), prot.begin(), prot.end());

			for (i = 0; i < protein_ids.size(); i++)
				for (auto it = protein_ids[i].precursors.begin(); it != protein_ids[i].precursors.end(); it++)
					entries[*it].pid_index = i;

			extract_proteins();
		}

		int ps = proteins.size(), pis = protein_ids.size();
		if (ps) if (!proteins[0].id.size()) ps--;
		if (pis) if (!protein_ids[0].ids.size()) pis--;

		if (!FastaSearch) eg.clear();
		elution_group_index();

		if (Verbose >= 1) Time(), dsout << "Spectral library loaded: "
			<< ps << " protein isoforms, "
			<< pis << " protein groups and "
			<< entries.size() << " precursors in " << co_elution_index.size() << " elution groups.\n";

		return true;
	}

	void export_prosit() {
		prosit_file = remove_extension(out_lib_file) + std::string(".prosit.csv");
		std::ofstream out(prosit_file);
		if (out.fail()) {
			dsout << "ERROR: cannot save precursors in the Prosit format to " << prosit_file << ". Check if the destination folder is write-protected or the file is in use";
			return;
		}
		out << "modified_sequence,collision_energy,precursor_charge\n";
		std::set<std::pair<std::pair<std::string, int>, float> > prs;
		for (auto &e : entries) prs.insert(std::pair<std::pair<std::string, int>, float>(std::pair<std::string, int>(to_prosit(e.name), e.target.charge), e.target.mz));
		for (auto &s : prs) out << s.first.first << ",30," << s.first.second << '\n';
		out.close();
		if (Verbose >= 1) Time(), dsout << "Prosit input saved to " << prosit_file << "\n";
	}

	void load_protein_annotations(Fasta &fasta) {
		std::set<PG> pid;
		std::string s;
		std::pair<std::string, std::vector<int> > key;
		for (int i = 0; i < entries.size(); i++) {
			bool prot_from_lib = false;
			if (OverwriteLibraryPGs || (entries[i].entry_flags & fFromFasta)) {
				auto stripped = get_aas(entries[i].name);
				if (stripped.size() < 5) prot_from_lib = true;
				else {
					key.first = stripped.substr(0, 5);
					auto pos = std::lower_bound(fasta.dict.begin(), fasta.dict.end(), key);
					if (pos != fasta.dict.end()) if (pos->first == key.first) { // checks needed when !(entries[i].entry_flags & fFromFasta)
						std::string ids;
						for (auto &seq : pos->second) {
							auto &sequence = fasta.sequences[seq];
							for (int shift = 0; shift <= sequence.size() - stripped.size(); shift++) {
								int found = sequence.find(stripped, shift);
								if (found != std::string::npos) {
									bool l_cut = false, r_cut = false;
									shift = found;
									if (found == 0) l_cut = true;
									else if (found == 1 && NMetExcision && (sequence[0] == 'M' || sequence[0] == 'X')) l_cut = true;
									else if (CutAfter.find(sequence[found - 1]) != std::string::npos) l_cut = true;
									if (l_cut) {
										if (found >= sequence.size() - stripped.size()) r_cut = true;
										else if (NoCutBefore.find(sequence[found + stripped.size()]) == std::string::npos) r_cut = true;
									}
									if (l_cut && r_cut) {
										if (!ids.size()) ids = fasta.ids[seq];
										else if (!name_included(ids, fasta.ids[seq])) ids += std::string(";") + fasta.ids[seq];
										break;
									}
								} else break;
							}
						}
						if (ids.size()) {
							auto ins = pid.insert(PG(ids));
							ins.first->precursors.push_back(i);
						} else prot_from_lib = true;
					} else prot_from_lib = true;
				}
			}
			if (prot_from_lib || (!(entries[i].entry_flags & fFromFasta) && !OverwriteLibraryPGs)) {
				auto ins = pid.insert(PG(protein_ids[entries[i].pid_index].ids));
				ins.first->precursors.push_back(i);
			}
		}

		protein_ids.clear();
		protein_ids.insert(protein_ids.begin(), pid.begin(), pid.end());
		pid.clear();

		for (int i = 0; i < protein_ids.size(); i++)
			for (auto it = protein_ids[i].precursors.begin(); it != protein_ids[i].precursors.end(); it++)
				entries[*it].pid_index = i;

		extract_proteins();
		if (Verbose >= 1) Time(), dsout << entries.size() << " precursors generated\n";

		if (ExportProsit) export_prosit();
	}

	void load(Fasta &fasta) {
		if (Verbose >= 1) Time(), dsout << "Processing FASTA\n";

		int i = entries.size(), cnt = 0;
		if (!i && !InSilicoRTPrediction) iRT_min = -1.0, iRT_max = 1.0;
		double default_iRT = (iRT_min + iRT_max) * 0.5;
		Entry e;
		e.lib = this; e.entry_flags = fFromFasta;
		for (auto &pep : fasta.peptides) {
			auto seq = get_sequence(pep);
			auto fragments = generate_fragments(seq, 1, loss_none, &cnt, MinFrAAs);
			auto aas = get_aas(pep);
			double mass = sum(seq) + proton + OH;
			int min_charge = Max(1, (int)(mass / MaxPrMz));
			int max_charge = Max(min_charge, Min(MaxGenCharge, (int)(mass / MinPrMz)));
			double iRT = InSilicoRTPrediction ? predict_irt(pep) : default_iRT;
			if (iRT < iRT_min) iRT_min = iRT;
			if (iRT > iRT_max) iRT_max = iRT;
			for (int charge = min_charge; charge <= max_charge; charge++) {
				double mz = proton + mass / (double)charge;
				if (mz < MinPrMz || mz > MaxPrMz) continue;
				e.target.init(mz, iRT, charge, i);
				e.target.fragments.clear();
				cnt = 0;
				for (auto &fr : fragments)
					if (fr.mz >= MinFrMz && fr.mz <= MaxFrMz)
						if (fr.type & fTypeY) cnt++;
				if (cnt >= MinGenFrNum) {
					e.name = to_canonical(pep, charge);
					auto pos = std::lower_bound(precursors.begin(), precursors.end(), e.name);
					if (pos == precursors.end() || *pos != e.name) {
						entries.push_back(e), i++;
						auto eg_name = to_eg(e.name);
						auto egp = eg.insert(std::pair<std::string, int>(eg_name, eg.size()));
						elution_groups.push_back(egp.first->second);
					}
				}
			}
		}
		eg.clear();
		elution_group_index();

		load_protein_annotations(fasta);
	}

	void library_protein_stats() {
		int ns = names.size(), gs = genes.size();
		if (ns) if (!names[0].size()) {
			ns--;
			if (Verbose >= 1) Time(), dsout << "Protein names missing for some isoforms\n";
		}
		if (gs) {
			if (!genes[0].size()) gs--;
			if (Verbose >= 1) Time(), dsout << "Gene names missing for some isoforms\n";
		}
		if (Verbose >= 1) Time(), dsout << "Library contains " << ns << " proteins, and " << gs << " genes\n";
	}

	void annotate() {
		std::set<std::string> name(names.begin(), names.end()), gene(genes.begin(), genes.end());
		for (auto &p : proteins) name.insert(p.name), gene.insert(p.gene);
		names.clear(), genes.clear();
		names.insert(names.begin(), name.begin(), name.end());
		genes.insert(genes.begin(), gene.begin(), gene.end());
		for (auto &p : proteins) {
			p.name_index = std::distance(names.begin(), std::lower_bound(names.begin(), names.end(), p.name));
			p.gene_index = std::distance(genes.begin(), std::lower_bound(genes.begin(), genes.end(), p.gene));
		}
		library_protein_stats();
	}

	void annotate_pgs(std::vector<PG> &pgs) {
		for (auto &pid : pgs) pid.annotate(proteins, names, genes);
	}

	void save(const std::string &file_name, std::vector<int> * ref, bool searched, bool decoys = false) {
		if (Verbose >= 1) Time(), dsout << "Saving spectral library to " << file_name << "\n";
		std::ofstream out(file_name, std::ofstream::out);
		if (out.fail()) { dsout << "ERROR: cannot write to " << file_name << ". Check if the destination folder is write-protected or the file is in use\n"; return; }
		out << "FileName\tPrecursorMz\tProductMz\tTr_recalibrated\ttransition_name\tLibraryIntensity\ttransition_group_id\tdecoy\tPeptideSequence\tProteotypic\tQValue\t";
		out << "ProteinGroup\tProteinName\tGenes\tFullUniModPeptideName\tModifiedPeptide\t";
		out << "PrecursorCharge\tPeptideGroupLabel\tUniprotID\tFragmentType\tFragmentCharge\tFragmentSeriesNumber\tFragmentLossType\tExcludeFromAssay\n";
		out.precision(8);

		auto &prot = (InferPGs && searched) ? protein_groups : protein_ids;

		bool rec_info = false;
		int cnt = -1, skipped = 0, empty = 0, targets_written = 0, decoys_written = 0;
		for (auto &it : entries) {
			cnt++;
			if (searched && it.best_run < 0 && ((FastaSearch && (it.entry_flags & fFromFasta)) || !SaveOriginalLib)) continue;
			if (!searched && it.target.fragments.size() < MinF) continue;
			if (ref != NULL) {
				auto pos = std::lower_bound(ref->begin(), ref->end(), cnt);
				if (pos == ref->end()) continue;
				if (*pos != cnt) continue;
			}
			auto seq = get_aas(it.name);
			float gen_acc = 0;

			bool recognise = false;
			std::vector<Ion> ions;
			for (auto &f : it.target.fragments) if (!(f.type & 3) || f.charge <= 0) {
				recognise = true;
				if (!rec_info) {
					rec_info = true;
					if (Verbose >= 1) Time(), dsout << "Some precursors lack fragment annotation; annotating...\n";
				}
				break;
			}
			if (recognise) {
				ions = recognise_fragments(it.name, it.target.fragments, gen_acc, it.target.charge, ExportRecFragOnly, MaxRecCharge, MaxRecLoss);
				if (ExportRecFragOnly && !searched) {
					int rec = 0;
					for (int fr = 0; fr < ions.size(); fr++) {
						ions[fr].type = it.target.fragments[fr].type;
						if (ions[fr].charge) rec++;
					}
					if (rec < MinF) continue;
				}
			} else {
				auto &frs = it.target.fragments;
				ions.resize(frs.size());
				for (int i = 0; i < ions.size(); i++) ions[i].init(frs[i]), ions[i].type = frs[i].type;
			}
			auto &pep = it.target;
			if (!pep.fragments.size()) {
				empty++;
				continue;
			}
			auto pep_name = to_canonical(it.name);
			auto name = pep_name + std::to_string(pep.charge);
			for (int dc = 0; dc <= 1; dc++) {
				if (!decoys && dc) continue;
				std::string prefix = dc ? std::string("DECOY") : "";
				auto &pep = dc ? it.decoy : it.target;
				int fr_cnt = 0;
				for (int fr = 0; fr < pep.fragments.size(); fr++) {
					int fr_type, fr_loss, fr_num, fr_charge = ions[fr].charge;
					float fr_mz = pep.fragments[fr].mz;
					if (fr_charge) {
						fr_type = ions[fr].type;
						fr_loss = ions[fr].loss;
						fr_num = ions[fr].index;
						if (ExportRecFragOnly && !dc) fr_mz = ions[fr].mz;
					} else { // unrecognised
						fr_type = fr_loss = fr_num = 0;
						if (ExportRecFragOnly) continue;
					}
					if (!(it.entry_flags & fFromFasta)) {
						bool skip = false;
						for (auto pos = ions.begin(); pos < ions.begin() + fr; pos++)
							if (pos->index == fr_num && pos->type == fr_type && pos->charge == fr_charge && pos->loss == fr_loss) {
								skip = true;
								break;
							}
						if (skip) continue;
					}
					if (pep.fragments[fr].height < MinGenFrInt) continue;
					int pg = (InferPGs && searched) ? it.pg_index : it.pid_index;
					auto fr_name = name + std::string("_") + std::to_string(char_from_type[fr_type & 3])
						+ std::string("_") + std::to_string(fr_charge) + std::string("_")
						+ std::to_string(fr_loss) + std::string("_") + std::to_string(fr_num);
					out << ((it.best_run >= 0) ? ms_files[it.best_run].c_str() : lib_file.c_str()) << '\t'
						<< pep.mz << '\t'
						<< fr_mz << '\t'
						<< pep.iRT << '\t'
						<< (prefix + fr_name).c_str() << '\t'
						<< pep.fragments[fr].height << '\t'
						<< (prefix + name).c_str() << '\t'
						<< dc << '\t'
						<< seq << '\t'
						<< (int)it.proteotypic << '\t'
						<< it.qvalue << '\t'
						<< (prefix + prot[pg].ids).c_str() << '\t'
						<< prot[pg].names.c_str() << '\t'
						<< prot[pg].genes.c_str() << '\t'
						<< pep_name.c_str() << '\t'
						<< pep_name.c_str() << '\t'
						<< pep.charge << '\t'
						<< (prefix + pep_name).c_str() << '\t'
						<< protein_ids[it.pid_index].ids.c_str() << '\t'
						<< char_from_type[fr_type & 3] << '\t'
						<< fr_charge << '\t'
						<< ((fr_type & fTypeB) ? fr_num : seq.size() - fr_num) << '\t'
						<< name_from_loss[fr_loss].c_str() << "\t"
						<< ((fr_type & fExclude) ? "True" : "False") << "\n";
					fr_cnt++;
				}
				if (!fr_cnt) empty++;
				else {
					if (!dc) targets_written++;
					else decoys_written++;
				}
			}
		}

		out.close();
		if (Verbose >= 1) {
			if (!decoys_written) Time(), dsout << targets_written << " precursors saved\n";
			else Time(), dsout << targets_written << " target and " << decoys_written << " decoy precursors saved\n";
		}
		if (skipped) dsout << "WARNING: " << skipped << " precursors with unrecognised library fragments were skipped\n";
		if (empty) dsout << "WARNING: " << empty << " precursors without any fragments annotated were skipped\n";
	}

	bool save(const std::string &file_name) {
		if (Verbose >= 1) Time(), dsout << "Saving the library to " << file_name << "\n";
		std::ofstream speclib(file_name, std::ofstream::binary);
		if (speclib.fail()) {
			dsout << "Could not save " << file_name << "\n";
			return false;
		} else {
			write(speclib);
			speclib.close();
			return true;
		}
	}

	void init(bool decoys) {
		bool corrected_charge = false;
		for (auto it = entries.begin(); it != entries.end(); it++) if (it->lock.set()) {
			if (it->target.charge <= 0) it->target.charge = 2, corrected_charge = true;
			if (gen_decoys && decoys) it->generate_decoy();
			if (it->decoy.charge <= 0 && decoys) it->decoy.charge = 2, corrected_charge = true;
			else if (UseIsotopes && gen_charges && !(it->entry_flags & fFromFasta)) it->annotate_charges();
			it->init();
		}
		if (corrected_charge) dsout << "WARNING: at least one library precursor had non-positive charge; corrected to charge = 2\n";
	}

	void initialise(bool decoys) {
		int i;

		if (Predictor) {
#ifdef PREDICTOR
			generate_predictions();
			predictor_file = remove_extension(out_lib_file) + std::string(".predicted.speclib");
			PredictorSaved = save(predictor_file);
#endif
		}
		skipped = 0;
		if (Verbose >= 1) Time(), dsout << "Initialising library\n";
		if (Threads > 1) {
			std::vector<std::thread> threads;
			for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Library::init, this, decoys));
			for (i = 0; i < Threads; i++) threads[i].join();
		} else init(decoys);

		for (i = 0; i < entries.size(); i++) entries[i].lock.free();
		if (skipped) dsout << "WARNING: " << skipped
			<< " library precursors were skipped due to unrecognised fragments when generating decoys.\n";
		if (ExportLibrary) save(out_lib_file, NULL, false, ExportDecoys);
	}

	void infer_proteins() { // filters precursors and proteins at q_cutoff level and infers protein groups
		int i;
		if (PGsInferred || !protein_ids.size() || !proteins.size()) return;
		protein_groups.clear();

		if (Verbose >= 1) Time(), dsout << "Assembling protein groups\n";

		// get IDs
		std::vector<std::pair<int, int> > IDs;
		std::vector<std::vector<int> > sets(proteins.size());
		for (auto it = info.map.begin(); it != info.map.end(); it++) {
			auto v = &(it->second);

			for (auto jt = v->begin(); jt != v->end(); jt++) {
				int s = jt->first;
				auto pr = &(jt->second.pr);
				if (pr->qvalue <= ProteinIDQvalue) {
					auto &pid = protein_ids[entries[pr->index].pid_index];
					bool first = true, swissprot = false;
					if (SwissProtPriority) for (auto &p : pid.proteins) if (proteins[p].swissprot) swissprot = true;
					for (auto &p : pid.proteins) {
						if (PGLevel == 0) if (!proteins[p].id.size()) continue;
						if (PGLevel == 1) if (!proteins[p].name.size()) continue;
						if (PGLevel == 2) if (!proteins[p].gene.size()) continue;
						if (swissprot && !proteins[p].swissprot) continue;
						auto pos = std::lower_bound(info.proteins[s].begin(), info.proteins[s].end(), p);
						if (pos != info.proteins[s].end()) if (*pos == p) {
							if (first) {
								sets[p].push_back(IDs.size());
								IDs.push_back(std::pair<int, int>(pr->index, s));
							} else sets[p].push_back(IDs.size() - 1);
							first = false;
						}
					}
				}
			}
		}

		auto groups = bipartite_set_cover(sets, IDs.size());

		// for each precursor merge all groups associated with it
		std::vector<std::set<int> > pgs(entries.size());
		for (auto &g : groups) for (auto &e : g.e) {
			auto &pr = pgs[IDs[e].first];
			for (auto &s : g.s)
				pr.insert(s);
		}

		std::set<PG> pg;
		std::string name;
		for (i = 0; i < entries.size(); i++) {
			auto &pid = protein_ids[entries[i].pid_index];
			if (!pgs[i].size()) name = pid.ids;
			else {
				auto it = pgs[i].begin();
				name = proteins[*it].id;
				for (it++; it != pgs[i].end(); it++)
					name += ';' + proteins[*it].id;
			}
			auto ins = pg.insert(PG(name));
			if (!pgs[i].size() && ins.second) ins.first->proteins = pid.proteins;
			for (auto &p : pgs[i]) ins.first->proteins.insert(p);
			ins.first->precursors.push_back(i);
		}
		protein_groups.clear();
		protein_groups.insert(protein_groups.begin(), pg.begin(), pg.end());

		for (i = 0; i < protein_groups.size(); i++)
			for (auto it = protein_groups[i].precursors.begin(); it != protein_groups[i].precursors.end(); it++)
				entries[*it].pg_index = i;

		bool annotate = false;
		if (names.size()) if (names[names.size() - 1] != "" && names[names.size() - 1] != " ") annotate = true;
		if (genes.size()) if (genes[genes.size() - 1] != "" && genes[genes.size() - 1] != " ") annotate = true;
		if (annotate) annotate_pgs(protein_groups);
	}

	void quantify_proteins(int N, double q_cutoff, int stage = 0) { // top N method for protein quantification
		if (InferPGs && !stage) infer_proteins();
		auto &prot = (InferPGs ? protein_groups : protein_ids);
		int i, j, pass, n = stage ? genes.size() : prot.size();

		if (stage == 2) {
			gg_index.clear(); gg_index.resize(entries.size(), 0);
			std::set<std::string> gg;
			for (auto &pg : prot) gg.insert(pg.genes);
			gene_groups.clear();  gene_groups.insert(gene_groups.begin(), gg.begin(), gg.end()); gg.clear();
			if (!gene_groups.size()) gene_groups.resize(1);

			for (auto &pg : prot) {
				int ind = std::distance(gene_groups.begin(), std::lower_bound(gene_groups.begin(), gene_groups.end(), pg.genes));
				for (auto &pr : pg.precursors) gg_index[pr] = ind;
			}
			n = gene_groups.size();
		}
		if (!n) {
			if (stage < 2) quantify_proteins(N, q_cutoff, stage + 1);
			return;
		}

		if (Verbose >= 1 && stage == 0) Time(), dsout << "Quantifying proteins\n";
		std::vector<float> max_q(n * info.n_s, 1.0), top_l(n * info.n_s * N, 0.0), quant(n * info.n_s, 0.0), norm(n * info.n_s, 0.0);

		for (pass = 0; pass <= 3; pass++) {
			for (auto it = info.map.begin(); it != info.map.end(); it++) {
				auto v = &(it->second);

				for (auto jt = v->begin(); jt != v->end(); jt++) {
					int s = jt->first;
					auto pr = &(jt->second.pr);

					int ggi = (stage == 2) ? gg_index[pr->index] : 0;
					if (stage == 2) if (!gene_groups[gg_index[pr->index]].size()) {
						if (pass == 0) pr->gene_quantity = pr->gene_norm = 0.0;
						continue;
					}
					if (stage == 1 && protein_ids[entries[pr->index].pid_index].gene_indices.size() != 1) {
						if (pass == 0) pr->gene_quantity_u = pr->gene_norm_u = 0.0;
						continue;
					}
					int pg = InferPGs ? entries[pr->index].pg_index : entries[pr->index].pid_index;
					int gi = (stage == 1) ? protein_ids[entries[pr->index].pid_index].gene_indices[0] : 0;
					if (stage == 1 && !genes[gi].size()) {
						if (pass == 0) pr->gene_quantity_u = pr->gene_norm_u = 0.0;
						continue;
					}

					float quantity = pr->norm;
					float q = pr->qvalue;

					int index = (stage == 2 ? ggi : (stage == 1 ? gi : pg)) * info.n_s + s;
					auto pos_q = &(max_q[index]);
					auto pos_l = &(top_l[index * N]);

					if (pass == 0) {
						if (q < *pos_q && *pos_q > q_cutoff) *pos_q = Max(q_cutoff, q);
					} else if (pass == 1) {
						if (q <= *pos_q) for (i = 0; i < N; i++) if (quantity > pos_l[i]) {
							for (j = N - 1; j > i; j--) pos_l[j] = pos_l[j - 1];
							pos_l[i] = quantity;
							break;
						}
					} else if (pass == 2) {
						if (q <= *pos_q && quantity >= pos_l[N - 1]) quant[index] += pr->quantity, norm[index] += pr->norm;
					} else if (pass == 3) {
						if (stage == 1) pr->gene_quantity_u = quant[index], pr->gene_norm_u = norm[index];
						else if (stage == 2) pr->gene_quantity = quant[index], pr->gene_norm = norm[index];
						else pr->pg_quantity = quant[index], pr->pg_norm = norm[index];
					}
				}
			}
		}
		if (stage < 2) quantify_proteins(N, q_cutoff, stage + 1);
	}

	void report(const std::string &file_name) {
		if (Verbose >= 1) Time(), dsout << "Writing report\n";
		auto &prot = InferPGs ? protein_groups : protein_ids;

		std::ofstream out(file_name, std::ofstream::out);
		if (out.fail()) { dsout << "ERROR: cannot write to " << file_name << ". Check if the destination folder is write-protected or the file is in use\n"; return; }
		out << oh[outFile] << '\t' << oh[outPG] << '\t' << oh[outPID] << '\t' << oh[outPNames] << '\t' << oh[outGenes] << '\t'
			<< oh[outPGQ] << '\t' << oh[outPGN] << '\t' << oh[outGQ] << '\t' << oh[outGN] << '\t' << oh[outGQP] << '\t' << oh[outGNP] << '\t' << oh[outModSeq] << '\t' << oh[outStrSeq] << '\t'
			<< oh[outPrId] << '\t' << oh[outCharge] << '\t' << oh[outQv] << '\t' << oh[outPQv] << '\t' << oh[outPPt] << '\t'
			<< oh[outPrQ] << '\t' << oh[outPrN] << '\t' << oh[outPrLR] << '\t'
			<< oh[outRT] << '\t' << oh[outiRT] << '\t' << oh[outpRT] << '\t' << oh[outpiRT];
		if (ExtendedReport) out << (fasta_files.size() ? "\tFirst.Protein.Description\t" : "\t") << "Lib.Q.Value\tMs1.Profile.Corr\tEvidence\tCScore\tDecoy.Evidence\tDecoy.CScore\tFragment.Quant.Raw\tFragment.Quant.Corrected\tFragment.Correlations";
#if REPORT_SCORES
		for (int i = 0; i < pN; i++) out << "\tScore." << i;
		for (int i = 0; i < pN; i++) out << "\tDecoy.Score." << i;
#endif
		out << "\n";

		for (auto it = info.map.begin(); it != info.map.end(); it++) {
			auto v = &(it->second);
			auto entry = &(entries[it->first]);
			int pg = InferPGs ? entry->pg_index : entry->pid_index;

			for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
				if (jt->second.pr.qvalue > ReportQValue) continue;
				if (jt->second.pr.protein_qvalue > ReportProteinQValue) continue;

				out << ms_files[jt->first].c_str() << '\t'
					<< prot[pg].ids.c_str() << '\t'
					<< protein_ids[entry->pid_index].ids.c_str() << '\t'
					<< prot[pg].names.c_str() << '\t'
					<< prot[pg].genes.c_str() << '\t'
					<< jt->second.pr.pg_quantity << '\t'
					<< jt->second.pr.pg_norm << '\t';
				if (jt->second.pr.gene_quantity > E) {
					out << jt->second.pr.gene_quantity << '\t'
						<< jt->second.pr.gene_norm << '\t';
					if (jt->second.pr.gene_quantity_u > E) {
						out << jt->second.pr.gene_quantity_u << '\t'
							<< jt->second.pr.gene_norm_u << '\t';
					} else out << "\t\t";
				} else out << "\t\t\t\t";
				out << pep_name(entry->name).c_str() << '\t'
					<< get_aas(entry->name).c_str() << '\t'
					<< entry->name.c_str() << '\t'
					<< entry->target.charge << '\t'
					<< jt->second.pr.qvalue << '\t'
					<< jt->second.pr.protein_qvalue << '\t'
					<< (int)entry->proteotypic << '\t'
					<< jt->second.pr.quantity << '\t'
					<< jt->second.pr.norm << '\t';
				
				if (jt->second.pr.ratio >= 0.0) out << jt->second.pr.ratio << '\t';
				else out << "NA\t";

				out << jt->second.pr.RT << '\t'
					<< jt->second.pr.iRT << '\t'
					<< jt->second.pr.predicted_RT << '\t'
					<< jt->second.pr.predicted_iRT;
				if (ExtendedReport) {
					out << '\t';
					if (fasta_files.size()) {
						if (prot[pg].proteins.size()) out << proteins[*(prot[pg].proteins.begin())].description.c_str() << '\t';
						else out << '\t';
					}
					out << entry->target.lib_qvalue << '\t'
						<< jt->second.pr.ms1_corr << '\t'
						<< jt->second.pr.evidence << '\t'
						<< jt->second.pr.cscore << '\t'
						<< jt->second.pr.decoy_evidence << '\t'
						<< jt->second.pr.decoy_cscore;
					out << '\t'; for (int fr = 0; fr < jt->second.fr_n; fr++) out << jt->second.fr[fr].quantity[qTotal] << ";";
					out << '\t'; for (int fr = 0; fr < jt->second.fr_n; fr++) out << jt->second.fr[fr].quantity[qFiltered] << ";";
					out << '\t'; for (int fr = 0; fr < jt->second.fr_n; fr++) out << jt->second.fr[fr].corr << ";";
				}
#if REPORT_SCORES
				for (int i = 0; i < pN; i++) out << '\t' << jt->second.pr.scores[i];
				for (int i = 0; i < pN; i++) out << '\t' << jt->second.pr.decoy_scores[i];
#endif
				out << "\n";
			}
		}

		out.close();
		if (Verbose >= 1) Time(), dsout << "Report saved to " << file_name << ".\n";
	}

	void gene_report(const std::string &file_name) {
		int s = ms_files.size(), size = gene_groups.size() * s;
		std::vector<float> gq(size), gn(size), gqu(size), gnu(size);
		if (!size || !gg_index.size()) return;

		if (Verbose >= 1) Time(), dsout << "Writing gene report\n";
		for (auto it = info.map.begin(); it != info.map.end(); it++) {
			auto v = &(it->second);
			auto entry = &(entries[it->first]);
			int pg = InferPGs ? entry->pg_index : entry->pid_index;

			for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
				if (jt->second.pr.qvalue > ReportQValue) continue;
				if (jt->second.pr.protein_qvalue > ReportProteinQValue) continue;
				int index = gg_index[it->first] * s + jt->first;
				gq[index] = jt->second.pr.gene_quantity;
				gn[index] = jt->second.pr.gene_norm;
				gqu[index] = jt->second.pr.gene_quantity_u;
				gnu[index] = jt->second.pr.gene_norm_u;
			}
		}

		std::ofstream out(file_name, std::ofstream::out);
		if (out.fail()) { dsout << "ERROR: cannot write to " << file_name << ". Check if the destination folder is write-protected or the file is in use\n"; return; }
		out << "File.Name\tGenes\tGene.Group.Quantity\tGene.Group.Normalised\tGene.Quantity.Unique\tGene.Normalised.Unique\n";

		for (int i = 0; i < gene_groups.size(); i++) for (int j = 0; j < s; j++) {
			int index = i * s + j;
			if (gq[index] > E) {
				out << ms_files[j].c_str() << '\t'
					<< gene_groups[i] << '\t'
					<< gq[index] << '\t'
					<< gn[index] << '\t';
				if (gqu[index] > E) {
					out << gqu[index] << '\t'
						<< gnu[index] << '\n';
				} else out << "\t\n";
			}
		}

		out.close();
		if (Verbose >= 1) Time(), dsout << "Gene report saved to " << file_name << ".\n";
	}

	void stats_report(const std::string &file_name) {
		if (!SaveRunStats) return;
		for (auto &rs : info.RS) {
			double mul = 1.0 / Max(1.0, (double)rs.cnt_valid);
			rs.fwhm_scans *= mul;
			rs.fwhm_RT *= mul;
			rs.pep_length *= mul;
			rs.pep_charge *= mul;
			rs.pep_mc *= mul;
			rs.MS1_acc *= 1000000.0;
			rs.MS1_acc_corrected *= 1000000.0;
			rs.MS2_acc *= 1000000.0;
			rs.MS2_acc_corrected *= 1000000.0;
		}

		std::ofstream out(file_name, std::ofstream::out);
		if (out.fail()) { dsout << "ERROR: cannot write to " << file_name << ". Check if the destination folder is write-protected or the file is in use\n"; return; }
		out << "File.Name\tPrecursors.Identified\tProteins.Identified\tTotal.Quantity\tMS1.Signal\tMS2.Signal\tFWHM.Scans\tFWHM.RT"
			<< "\tMedian.Mass.Acc.MS1\tMedian.Mass.Acc.MS1.Corrected\tMedian.Mass.Acc.MS2\tMedian.Mass.Acc.MS2.Corrected\tMedian.RT.Prediction.Acc"
			<< "\tAverage.Peptide.Length\tAverage.Peptide.Charge\tAverage.Missed.Tryptic.Cleavages\n";

		for (int i = 0; i < info.RS.size(); i++) {
			auto &rs = info.RS[i];
			out << ms_files[i] << "\t" << rs.precursors << "\t" << rs.proteins << "\t" << rs.tot_quantity << "\t" << rs.MS1_tot << "\t" << rs.MS2_tot
				<< "\t" << rs.fwhm_scans << "\t" << rs.fwhm_RT << "\t" << rs.MS1_acc << "\t" << rs.MS1_acc_corrected
				<< "\t" << rs.MS2_acc << "\t" << rs.MS2_acc_corrected << "\t" << rs.RT_acc << "\t" << rs.pep_length << "\t" << rs.pep_charge << "\t" << rs.pep_mc << "\n";
		}
		out.close();

		if (Verbose >= 1) Time(), dsout << "Stats report saved to " << file_name << "\n";
	}

	void tic_report(const std::string &file_name) {
		if (!SaveRunStats) return;

		std::ofstream out(file_name, std::ofstream::out);
		if (out.fail()) { dsout << "ERROR: cannot write to " << file_name << ". Check if the destination folder is write-protected or the file is in use\n"; return; }
		for (int i = 0; i < TIC_n; i++) out << (((float)i) / (float)TIC_n) << "\t"; out << "\n";
		for (auto &rs : info.RS) { for (int i = 0; i < TIC_n; i++) out << rs.TIC[i] << "\t"; out << "\n"; }
		out.close();
		if (Verbose >= 1) Time(), dsout << "TICs saved to " << file_name << "\n";
	}

	void XIC_report(const std::string& file_name) {
		if (!XICs.size()) {
			dsout << "No XICs recorded\n";
			return;
		}
		std::ofstream out(file_name, std::ofstream::out);
		if (out.fail()) { dsout << "ERROR: cannot write to " << file_name << ". Check if the destination folder is write-protected or the file is in use\n"; return; }

		int i, cnt = 0, in, W = VisWindowRadius * 2 + 1;
		out << "File.Name\tPrecursor.Id\tModified.Sequence\tStripped.Sequence\tQ.Value\tMS.Level\tIntensities\tRetention.Times\tTheoretical.Mz\tFragmentType\tFragmentCharge\tFragmentSeriesNumber\tFragmentLossType";
		for (i = 0; i < W; i++) out << '\t' << i; out << '\n';

		for (auto& xic : XICs) {
			if (!xic.peaks.size()) continue;
			auto& le = entries[xic.pr];
			auto aas = get_aas(le.name);
			for (in = 0; in <= 1; in++) {
				out << ms_files[xic.run] << '\t'
					<< le.name << '\t'
					<< pep_name(le.name) << '\t'
					<< aas << '\t'
					<< xic.qvalue << '\t'
					<< xic.level << '\t'
					<< in << '\t'
					<< (in ^ 1) << '\t';
				if (xic.level == 1) out << xic.pr_mz << "\tNA\tNA\tNA\tNA";
				else out << xic.fr_mz << '\t' << char_from_type[xic.fr.type & 3] << '\t' << ((int)xic.fr.charge) << '\t'
					<< ((xic.fr.type & fTypeB) ? ((int)xic.fr.index) : aas.size() - ((int)xic.fr.index)) << '\t' << name_from_loss[xic.fr.loss].c_str();

				if (!in) { for (i = 0; i < W; i++) out << '\t' << xic.peaks[i].first; out << '\n'; }
				if (in) { for (i = 0; i < W; i++) out << '\t' << xic.peaks[i].second; out << '\n'; }
			}
			cnt++;
		}

		out.close();
		if (Verbose >= 1) Time(), dsout << cnt << " XICs saved to " << file_name << "\n";
	}

#if (HASH > 0)
	unsigned int hash() {
		unsigned int res = 0;
		for (int i = 0; i < entries.size(); i++) res ^= entries[i].hash();
		return res;
	}
#endif
};

void annotate_library(Library &lib, Fasta &fasta) {
	for (auto &e : lib.proteins) {
		auto pos = std::lower_bound(fasta.proteins.begin(), fasta.proteins.end(), e);
		if (pos != fasta.proteins.end()) if (pos->id == e.id) e.name = pos->name, e.gene = pos->gene, e.description = pos->description, e.swissprot = pos->swissprot;
	}
	lib.annotate();
	lib.annotate_pgs(lib.protein_ids);
	for (auto &le : lib.entries) {
		if (lib.protein_ids[le.pid_index].proteins.size() == 1) le.proteotypic = true;
		else {
			if (PGLevel == 0) le.proteotypic = false;
			else if (PGLevel == 1) le.proteotypic = (lib.protein_ids[le.pid_index].name_indices.size() == 1);
			else if (PGLevel == 2) le.proteotypic = (lib.protein_ids[le.pid_index].gene_indices.size() == 1);
		}
	}
}

void train(int thread_id, std::vector<NN> * net) {
	for (int i = 0; i * Threads + thread_id < nnBagging; i++)
		(*net)[i * Threads + thread_id].optimise();
}

void bagging_score(int thread_id, std::vector<NN> * net, int in, std::vector<float*> * input, std::vector<float> * nn_sc) {
	int i, j, cnt, batch_size = 2000, bagging = net->size();
	for (int n = 0; n * Threads + thread_id < bagging; n++) {
		int ind = n * Threads + thread_id;

		for (i = 0; i < input->size(); i += cnt) {
			cnt = Min(batch_size, input->size() - i);
			auto dataset = createDataSet(cnt, in, &((*input)[i]));
			forwardPassDataSet((*net)[ind].network, dataset);
			auto result = getOuput((*net)[ind].network);
			for (j = 0; j < cnt; j++) (*nn_sc)[(i + j) * 2 * bagging + ind] = result->data[2 * j], (*nn_sc)[(i + j) * 2 * bagging + ind + 1] = result->data[2 * j + 1];
			free(dataset);
		}
	}
}

class NNClassifier {
public:
	int in, hidden;
	std::vector<NN> net;
	std::vector<int> size;

	NNClassifier(int bagging, int _in, int _hidden, std::vector<int> &_size, std::vector<std::vector<float*> > &data, std::vector<std::vector<float*> > &classes) {
		int i;
		in = _in, hidden = _hidden, size = _size;

		int bg = Min(Min(size.size(), data.size()), classes.size());
		if (!bg) {
			dsout << "ERROR: no data provided to the NNClassifier() constructor, aborting\n";
			return;
		}
		if (bagging > bg) {
			dsout << "ERROR: bagging factor reduced to " << bg << " to match the data provided\n";
			bagging = bg;
		}

		size_t* hiddenSize = (size_t*)alloca(hidden * sizeof(size_t));
		Activation* hiddenActivation = (Activation*)alloca(hidden * sizeof(Activation));
		for (i = 0; i < hidden; i++) {
			hiddenSize[i] = (hidden - i) * 5;
			hiddenActivation[i] = tanH;
		}

		net.resize(bagging);
		for (i = 0; i < bagging; i++) {
			int si = Min(i, size.size() - 1), di = Min(i, data.size() - 1), ci = Min(i, classes.size() - 1);
			DataSet* trainingData = createDataSet(size[si], in, &(data[di][0]));
			DataSet* trainingClasses = createDataSet(size[si], 2, &(classes[ci][0]));

			net[i].network = createNetwork(in, hidden, hiddenSize, hiddenActivation, 2, softmax, net[i].random);
			net[i].seed(i);
			net[i].lossFunction = CROSS_ENTROPY_LOSS;
			net[i].batchSize = Min(50, Max(1, size[si] / 100));
			net[i].learningRate = nnLearning;
			net[i].searchTime = 0;
			net[i].regularizationStrength = Regularisation;
			net[i].B1 = 0.9, net[i].B2 = 0.999;
			net[i].shuffle = 1;
			net[i].verbose = (Verbose >= 5);
			net[i].data = trainingData;
			net[i].classes = trainingClasses;
		}
	}

	void run(int epochs) {
		int i;
		for (i = 0; i < net.size(); i++) net[i].maxIters = (size[Min(i, size.size() - 1)] * epochs) / net[i].batchSize;

		if (net.size() > 1 && Threads > 1) {
			std::vector<std::thread> thr;
			for (i = 0; i < Threads; i++) thr.push_back(std::thread(train, i, &net));
			for (i = 0; i < Threads; i++) thr[i].join();
		} else train(0, &net);
	}

	void predict(std::vector<float*> &data, std::vector<float*> &prediction) {
		int i, j, bagging = net.size();

		std::vector<float> nn_sc(data.size() * bagging * 2, -1.0);
		if (bagging > 1 && Threads > 1) {
			std::vector<std::thread> thr;
			for (i = 0; i < Threads; i++) thr.push_back(std::thread(&bagging_score, i, &net, in, &data, &nn_sc));
			for (i = 0; i < Threads; i++) thr[i].join();
		} else bagging_score(0, &net, in, &data, &nn_sc);

		for (i = 0; i < data.size(); i++) {
			float sc_one = 0.0, sc_two = 0.0;
			int cnt = 0;
			for (j = 0; j < bagging; j++) if (nn_sc[i * bagging * 2 + j] >= 0.0) sc_one += nn_sc[i * bagging * 2 + j], sc_two += nn_sc[i * bagging * 2 + j + 1], cnt++;
			if (cnt) prediction[i][0] = sc_one / (double)cnt, prediction[i][1] = sc_two / (double)cnt;
			else prediction[i][0] = prediction[i][1] = -1.0;
		}
	}

	~NNClassifier() {
		for (int i = 0; i < net.size(); i++) {
			free(net[i].data);
			free(net[i].classes);
			destroyNetwork(net[i].network);
		}
	}
};

void learn_from_library(const std::string &file_name) {
	int i, j, n, P_y = y_delta_size(), P_l = y_loss_size(), M = in_silico_rt_size(), pos_y = 0, pos_l = 0, pos_RT = 0;
	Library lib;
	if (!lib.load(&(file_name[0]))) Throw("Cannot load the library");
	if (Verbose >= 1) Time(), dsout << "Learning peptide characteristics\n";

	yDelta.resize(P_y);
	InSilicoRT.resize(M);
	std::vector<int> y_index;

	std::vector<double> b_y, b_l, b_RT;
	typedef Eigen::Triplet<double> T;
	std::vector<T> TL_y, TL_l, TL_RT;
	std::vector<double> count(M);
	for (auto &e : lib.entries) {
		auto aas = get_aas(e.name);
		float gen_acc = 0;
		auto frs = recognise_fragments(e.name, e.target.fragments, gen_acc, e.target.charge, true, 1, AddLosses ? (loss_NH3 + 1) : 1);
		if (gen_acc > GeneratorAccuracy + E) continue;

		if (InSilicoRTPrediction) {
			count[0] = 1;
			for (i = 1; i < M; i++) count[i] = 0;
			count_rt_aas(count, e.name);
			for (i = 0; i < M; i++) if (Abs(count[i]) > E) TL_RT.push_back(T(pos_RT, i, (double)count[i]));
			b_RT.push_back(e.target.iRT);
			pos_RT++;
		}

		y_index.resize(n = aas.size());
		for (i = 0; i < y_index.size(); i++) y_index[i] = -1;
		for (i = 0; i < frs.size(); i++) if (frs[i].charge) {
			auto &fr = frs[i];
			if (fr.charge == 1 && (fr.type & fTypeY)) if (fr.loss == loss_none) y_index[fr.index] = i;
		}

		if (AddLosses) for (i = 0; i < frs.size(); i++) if (frs[i].charge) {
			auto &fr = frs[i];
			if (fr.charge == 1 && (fr.type & fTypeY)) if (y_index[fr.index] >= 0) if (fr.loss == loss_H2O || fr.loss == loss_NH3) {
				double ratio = log(e.target.fragments[i].height / e.target.fragments[y_index[fr.index]].height);
				b_l.push_back(ratio);

				TL_l.push_back(T(pos_l, y_loss(fr.loss == loss_NH3), 1.0));
				for (j = fr.index; j < aas.length(); j++) TL_l.push_back(T(pos_l, y_loss_aa(aas[j], fr.loss == loss_NH3), 1.0));
				pos_l++;
			}
		}
		if (AddLosses && !pos_l) {
			AddLosses = false;
			if (Verbose >= 1) Time(), dsout << "No H2O/NH3 neutral losses in the library, losses will not be used\n";
		}

		for (i = 2; i < y_index.size(); i++) if (y_index[i] >= 0 && y_index[i - 1] >= 0) {
			double ratio = log(e.target.fragments[y_index[i]].height / e.target.fragments[y_index[i - 1]].height);
			if (Abs(ratio) >= 5.0) continue;
			b_y.push_back(ratio);

			for (j = 0; j < aas.length(); j++) TL_y.push_back(T(pos_y, y_delta_composition_index(aas[j]), 1.0));
			TL_y.push_back(T(pos_y, y_delta_charge_index(), (double)e.target.charge));
			for (j = Max(i - yDeltaS, 0); j <= Min(i + yDeltaS, n - 1); j++)
				TL_y.push_back(T(pos_y, y_delta_index(aas[j], j - i), 1.0));
			TL_y.push_back(T(pos_y, y_nterm_index(i - 2), 1.0));
			TL_y.push_back(T(pos_y, y_cterm_index(aas[n - 1], n - 1 - i), 1.0));

			pos_y++;
		}
	}

	if (pos_y) {
		auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b_y.data(), b_y.size());
		Eigen::SparseMatrix<double, Eigen::ColMajor> A(pos_y, P_y);
		A.setFromTriplets(TL_y.begin(), TL_y.end());
		Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
		SQR.compute(A);
		auto X = SQR.solve(B);
		for (i = 0; i < P_y; i++) yDelta[i] = X[i];
	}

	if (AddLosses) {
		auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b_l.data(), b_l.size());
		Eigen::SparseMatrix<double, Eigen::ColMajor> A(pos_l, P_l);
		A.setFromTriplets(TL_l.begin(), TL_l.end());
		Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
		SQR.compute(A);
		auto X = SQR.solve(B);
		for (i = 0; i < P_l; i++) yLoss[i] = X[i];
	}

	if (InSilicoRTPrediction && pos_RT) {
		auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b_RT.data(), b_RT.size());
		Eigen::SparseMatrix<double, Eigen::ColMajor> A(pos_RT, M);
		A.setFromTriplets(TL_RT.begin(), TL_RT.end());
		Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
		SQR.compute(A);
		auto X = SQR.solve(B);
		for (i = 0; i < M; i++) InSilicoRT[i] = X[i];
	}

	// Pearson correlation between predicted and actual spectra
	double r_y = 0.0, s2_y = 0.0;
	int cnt = 0, s2_cnt = 0;
	std::vector<double> actual_y, predicted_y, rt_d, rm_y;
#if (HASH > 0)
	unsigned int learn_hash = 0;
#endif
	for (auto &e : lib.entries) {
		auto aas = get_aas(e.name);
		float gen_acc = 0;
		auto frs = recognise_fragments(e.name, e.target.fragments, gen_acc, e.target.charge, true, 1, 1);
		if (gen_acc > GeneratorAccuracy + E) continue;

		if (InSilicoRTPrediction) {
			double nrt = predict_irt(e.name);
			rt_d.push_back(Abs(nrt - e.target.iRT));
		}

		predicted_y.resize(n = aas.size());
		actual_y.resize(n);
		for (int i = 0; i < n; i++) actual_y[i] = 0.0;

		for (int i = 0; i < frs.size(); i++) if (frs[i].charge) {
			auto &fr = frs[i];
			if (fr.charge == 1 && fr.loss == loss_none)
				if (fr.type & fTypeY) actual_y[fr.index] = e.target.fragments[i].height;
		}
		y_scores(predicted_y, e.target.charge, aas);
		for (i = 2; i < aas.size(); i++) if (actual_y[i] > E && actual_y[i - 1] > E) s2_y += Sqr(log(actual_y[i] / actual_y[i - 1]) - (predicted_y[i] - predicted_y[i - 1])), s2_cnt++;
		to_exp(predicted_y);
#if (HASH > 0)
		learn_hash ^= hashA((int*)&(predicted_y[0]), 2 * (aas.size() - 1));
#endif
		double c = corr(&(predicted_y[1]), &(actual_y[1]), aas.size() - 1);
		r_y += c;
		rm_y.push_back(c);
		cnt++;
	}
#if (HASH > 0)
	dsout << "Fragmentation prediction hash: " << learn_hash << "\n";
#endif
	std::sort(rm_y.begin(), rm_y.end());
	Time(), dsout << "y-series fragmentation prediction: ratio SD = " << sqrt(s2_y / (double)Max(1, s2_cnt - 1))
		<< ", Pearson correlation = " << r_y / (double)cnt << " average, " << rm_y[cnt / 2] << " median\n";
	if (InSilicoRTPrediction) {
		std::sort(rt_d.begin(), rt_d.end());
		Time(), dsout << "iRT prediction: median error = " << rt_d[rt_d.size() / 2] << "\n";
#if (HASH > 0)
		dsout << "LRT prediction hash: " << hashA((int*)&(rt_d[0]), 2 * (rt_d.size())) << "\n";
#endif
		RTLearnLib = true;
	}
}

const int fFound = 1;
const int fDecoy = 1 << 1;
const int fInit = 1 << 2;
const int fTranslated = 1 << 3;

class Run {
public:
	Lock lock;
	volatile int ms1_cnt = 0, ms2_cnt = 0, total_spectra = 0;
	volatile long long peaks_cnt;
	std::atomic<int> sp_alloc;
	bool scanning = false, no_report = false;

	int run_index = 0, LDA = 0;
	RunStats RS;
    std::string name; // run name
	float weights[pN], guide_weights[pN], best_weights[pN], best_guide_weights[pN], selection_weights[pN], selection_guide_weights[pN];
	bool default_weights = true;
	std::vector<Scan> ms1h, ms2h; // headers
	std::vector<Peak> peaks;
	std::list<std::vector<Peak> > peak_lists;

	std::vector<std::vector<int> > IDs;
	std::vector<std::pair<int, float> > best_egs;

	Library * lib;

    int n_scans = 0, cycle_length = 0;
	std::vector<double> RT_coeff, RT_points, iRT_coeff, iRT_points;
    std::vector<float> ms1_RT, scan_RT; // list of RTs of scans
	std::vector<int> scan_cycle; // list of SWATH cycle numbers
	double RT_window, RT_min = 0.0, RT_max = 0.0;
	double tandem_min = 0.0, tandem_max = INF;
	bool full_spectrum = false;

	std::vector<std::vector<float> > NF, MS1, MS2, MS2_H, MS2_min, MS2_tight_one, MS2_tight_two, Shadow;
	std::vector<std::vector<bool> > Covered;
	std::vector<std::vector<int> > PeakList, BestFrList, ScanIndex, ChromIndex;
	std::vector<std::vector<float> > CorrSumList;
	std::vector<std::vector<Product> > FR;
	std::vector<std::vector<int> > FRI;

	std::vector<int> rt_stats, rt_ref;
	std::vector<float> rt_coo, q1_diff, q1_diff_mz;
	std::vector<std::pair<float, float> > rt_data;
	std::vector<float> rt_delta, mass_acc;
	std::vector<double> b, b_r;

	int Ids01, Ids1, Ids10, Ids50, IdsCal;

	bool nn_trained = false, all_iters = false, nn_validated = false;
    int curr_iter, curr_batch, linear_classifier_ids1 = 0, linear_classifier_ids10 = 0;
    double iRT_cscore, iRT_ref_score, PeakWidth = 0.0, MassAccuracy = GlobalMassAccuracy, MassAccuracyMs1 = GlobalMassAccuracyMs1;
	double MassAccOutlier = INF, MassAccMs1Outlier = INF;
	bool RemoveMassAccOutliers = false, recalibrate = false;
	std::vector<double> MassCorrection, MassCorrectionMs1, MassCalSplit, MassCalSplitMs1, MassCalCenter, MassCalCenterMs1;
	double Q1Correction[2] = { 0.0, 0.0 };

    std::vector<Parameter> pars;
    bool par_seek[pN], par_learn[pN], par_limit = false;

	bool in_ref_run = false, RT_windowed_search = false, mz_calibrated = false, acc_calibrated = false, acc_ms1_calibrated = false, window_calculated = false;

	void reset_weights() {
		for (int i = 1; i < pN; i++) weights[i] = guide_weights[i] = 0.0;
		weights[0] = guide_weights[0] = 1.0;
		default_weights = true;
	}

	void copy_weights(float *dst, float *src) { for (int i = 0; i < pN; i++) dst[i] = src[i]; }

	void check_weights() {
		if (LCAllScores && curr_iter == nnIter - 1) return;
		if (weights[pdRT] > 0.0) weights[pdRT] = 0.0;
		if (guide_weights[pdRT] > 0.0) guide_weights[pdRT] = 0.0;
		for (int i = 0; i < TopF; i++) {
			if (weights[pAcc + i] > 0.0) weights[pAcc + i] = 0.0;
			if (guide_weights[pAcc + i] > 0.0) guide_weights[pAcc + i] = 0.0;
		}
	}

	Run(int _run_index) {
		run_index = _run_index;

		auto p_base = Parameter(0, 0);
		auto p_all = Parameter(1, 0);
		auto p_end = Parameter(CalibrationIter + 2, CalibrationIter + 2);
		auto p_fit = Parameter(CalibrationIter + 2, 0);
		auto p_learn = Parameter(iN, CalibrationIter + 2);
		auto p_none = Parameter(iN, iN);

		pars.push_back(p_base); // pTimeCorr
		pars.push_back(p_all); // pLocCorr
		pars.push_back(p_end); // pMinCorr
		pars.push_back(p_end); // pTotal
		pars.push_back(p_end); // pCos
		pars.push_back(p_all); // pCosCube
		pars.push_back(p_all); // pMs1TimeCorr
		pars.push_back(p_end); // pNFCorr
		pars.push_back(p_end); // pdRT
		pars.push_back(p_none); // pResCorr
		pars.push_back(p_none); // pResCorrNorm
		pars.push_back(TightMassAauxForCal ? p_fit : p_end); // pTightCorrOne
		pars.push_back(TightMassAauxForCal ? p_fit : p_end); // pTightCorrTwo
		pars.push_back(p_end); // pShadow
		pars.push_back(p_end); // pHeavy
		pars.push_back(p_end); // pMs1TightOne
		pars.push_back(p_end); // pMs1TightTwo
		pars.push_back(p_end); // pMs1Iso
		pars.push_back(p_end); // pMs1IsoOne
		pars.push_back(p_end); // pMs1IsoTwo
		pars.push_back(p_none); // pMs1Ratio
		pars.push_back(p_none); // pBestCorrDelta
		pars.push_back(p_none); // pTotCorrSum
		pars.push_back(p_none); // pMz
		pars.push_back(p_none); // pCharge
		pars.push_back(p_none); // pLength
		pars.push_back(p_none); // pFrNum
#if AAS
		pars.push_back(p_none); // pMods
		for (int i = 0; i < 20; i++) pars.push_back(p_none); // pAAs
#endif
		pars.push_back(p_none); // pRT
		for (int i = 0; i < TopF; i++) pars.push_back(p_end); // pAcc
		for (int i = 1; i < auxF; i++) pars.push_back(p_none); // pRef
		for (int i = 0; i < TopF; i++) pars.push_back(i < TopF - 1 ? p_end : p_none); // pSig
		for (int i = 0; i < auxF; i++) pars.push_back(p_none); // pCorr
		for (int i = 0; i < TopF; i++) pars.push_back(p_none); // pShadowCorr
		for (int i = 0; i < nnW; i++) pars.push_back(p_none); // pShape
#if Q1
		for (int i = 0; i < qL; i++) pars.push_back(p_learn); // pQLeft
		for (int i = 0; i < qL; i++) pars.push_back(p_learn); // pQRight
		for (int i = 0; i < qL; i++) pars.push_back(p_learn); // pQPos
		for (int i = 0; i < qL; i++) pars.push_back(p_learn); // PQNFCorr
		for (int i = 0; i < qL; i++) pars.push_back(p_learn); // pQCorr
#endif
		assert(pars.size() == pN);

		reset_weights();
		copy_weights(selection_weights, weights), copy_weights(selection_guide_weights, weights), copy_weights(best_weights, weights), copy_weights(best_guide_weights, weights);
		training.resize(nnBagging);
		training_classes.resize(nnBagging);
		nn_mod.resize(nnBagging), nn_shift.resize(nnBagging), nn_scale.resize(nnBagging);

		NF.resize(Threads), MS1.resize(Threads), MS2.resize(Threads), MS2_H.resize(Threads), MS2_min.resize(Threads);
		MS2_tight_one.resize(Threads), MS2_tight_two.resize(Threads), Shadow.resize(Threads);
		PeakList.resize(Threads), BestFrList.resize(Threads), CorrSumList.resize(Threads);
		ScanIndex.resize(Threads), ChromIndex.resize(Threads), FR.resize(Threads), FRI.resize(Threads), Covered.resize(Threads);

		if (Calibrate && !QuantOnly) MassAccuracy = MassAccuracyMs1 = CalibrationMassAccuracy;
		MassCorrection.resize(1 + MassCalBinsMax * 2, 0.0); MassCorrectionMs1.resize(1 + MassCalBinsMax * 2, 0.0);
		MassCalSplit.resize(MassCalBinsMax + 1, 0.0); MassCalSplitMs1.resize(MassCalBinsMax + 1, 0.0);
		MassCalCenter.resize(MassCalBinsMax + 1, 0.0); MassCalCenterMs1.resize(MassCalBinsMax + 1, 0.0);

		memset(&(RS), 0, sizeof(RunStats));
	}

	inline double predicted_mz(double * t, double mz, double rt) {
		double s = t[0] * Sqr(mz);
		if (rt <= MassCalCenter[0]) s += t[1] + t[2] * mz;
		else if (rt >= MassCalCenter[MassCalBins - 1]) s += t[1 + (MassCalBins - 1) * 2] + t[2 + (MassCalBins - 1) * 2] * mz;
		else for (int i = 1; i < MassCalBins; i++) if (rt < MassCalCenter[i]) {
			double u = rt - MassCalCenter[i - 1], v = MassCalCenter[i] - rt, w = u + v;
			if (w > E) s += ((t[1 + (i - 1) * 2] + t[2 + (i - 1) * 2] * mz) * v + (t[1 + i * 2] + t[2 + i * 2] * mz) * u) / w;
			break;
		}
		return s + mz;
	}

	inline double predicted_mz_ms1(double * t, double mz, double rt) {
		double s = t[0] * Sqr(mz);
		if (rt <= MassCalCenterMs1[0]) s += t[1] + t[2] * mz;
		else if (rt >= MassCalCenterMs1[MassCalBinsMs1 - 1]) s += t[1 + (MassCalBinsMs1 - 1) * 2] + t[2 + (MassCalBinsMs1 - 1) * 2] * mz;
		else for (int i = 1; i < MassCalBinsMs1; i++) if (rt < MassCalCenterMs1[i]) {
			double u = rt - MassCalCenterMs1[i - 1], v = MassCalCenterMs1[i] - rt, w = u + v;
			if (w > E) s += ((t[1 + (i - 1) * 2] + t[2 + (i - 1) * 2] * mz) * v + (t[1 + i * 2] + t[2 + i * 2] * mz) * u) / w;
			break;
		}
		return s + mz;
	}

    template <bool seek> inline double cscore(float * scores, bool guide = false) {
		if (seek) return scalar(scores, guide ? guide_weights : weights, par_seek, pN);
		else return scalar(scores, guide ? guide_weights : weights, par_learn, pN);
    }

#if Q1
#define WinCenter(x) (0.5 * (ms2h[x].window_low + ms2h[x].window_high))

	void q1_profiles(float * profiles, float * coo, int apex, float * fmz, int QS, int m) {
		int i, j, k, pos, QW = 2 * QS + 1;
		float rt = scan_RT[apex];
		float wc = WinCenter(apex);
		for (i = 0; i < QW * m; i++) profiles[i] = 0.0;
		for (i = QS + 1; i < QW; i++) coo[i] = INF;
		for (i = 0; i < QS; i++) coo[i] = -INF;
		for (j = 0; j < m; j++) profiles[j * QW + QS] = ms2h[apex].level<false>(fmz[j], MassAccuracy);
		coo[QS] = wc;

		for (int inc = -1; inc <= 1; inc += 2) {
			int lsgn = Sgn(WinCenter(apex + inc) - wc);
			for (pos = apex + inc; pos >= 0 && pos < ms2h.size(); pos += inc) {
				float w = WinCenter(pos);
				if (Sgn(w - wc) != lsgn) break;
				if (w - wc >= 0.0) for (i = QS + 1; i < QW; i++) if (w < coo[i]) {
					for (int j = QW - 1; j > i; j--) {
						coo[j] = coo[j - 1];
						for (k = 0; k < m; k++) profiles[k * QW + j] = profiles[k * QW + j - 1];
					}
					coo[i] = w;
					for (k = 0; k < m; k++) profiles[k * QW + i] = ms2h[pos].level<false>(fmz[k], MassAccuracy);
					break;
				}
				if (w - wc < 0.0) for (i = QS - 1; i >= 0; i--) if (w > coo[i]) {
					for (int j = 0; j < i; j++) {
						coo[j] = coo[j + 1];
						for (k = 0; k < m; k++) profiles[k * QW + j] = profiles[k * QW + j + 1];
					}
					coo[i] = w;
					for (k = 0; k < m; k++) profiles[k * QW + i] = ms2h[pos].level<false>(fmz[k], MassAccuracy);
					break;
				}
			}
		}
	}

	double q1_delta(int apex, float pr_mz, float mz, int QS) {
		float _mz = mz;
		int QW = 2 * QS + 1;
		float *profiles = (float*)alloca(QW * sizeof(float)), *coo = (float*)alloca(QW * sizeof(float));
		q1_profiles(profiles, coo, apex, &_mz, QS, 1);
		return centroid_coo(profiles, coo, QW) - pr_mz;
	}
#endif

	void ms1_features(std::vector<std::vector<Peak> > &features, int N = 500, float accuracy = 12.0 / 1000000.0) {
		std::vector<Peak> top(N * ms1h.size()), temp;
		memset(&(top[0]), 0, top.size() * sizeof(Peak));
		int S;
		if (WindowRadius) S = WindowRadius;
		else if (PeakWidth > E) S = Max(1, int(ScanScale * PeakWidth));
		else S = Max(1, ms1h.size() / ScanFactor); 
		int W = 2 * S + 1;
		std::vector<float> ch(W * N), elution(W);
		std::vector<int> ind(W, 0);

		features.resize(ms1h.size());
		int i, j, pos, jj, jjN, k, l, SN = S * N, iN;
		for (i = 0; i < ms1h.size(); i++) ms1h[i].top_peaks(&(top[i * N]), N, temp);
		for (i = S; i < ms1h.size() - S; i++) {
			iN = i * N;
			for (j = 0; j < W; j++) ind[j] = 0;
			memset(&(ch[0]), 0, ch.size() * sizeof(float));
			for (k = 0; k < N; k++) {
				float mz = top[iN + k].mz, margin = mz * accuracy, low = mz - margin, high = mz + margin;
				ch[SN + k] = top[iN + k].height;
				for (j = 0; j < W; j++) {
					jj = i - S + j, jjN = jj * N;
					for (l = ind[j]; l < N; l++) {
						if (top[jjN + l].mz <= low) continue;
						if (top[jjN + l].mz >= high) break;
						if (top[jjN + l].height > ch[j * N + k]) ch[j * N + k] = top[jjN + l].height;
					}
					ind[j] = l;
				}
			}

			for (k = 0; k < N; k++) {
				elution[0] = (2.0 / 3.0) * ch[0 + k] + (1.0 / 3.0) * ch[N + k];
				elution[W - 1] = (2.0 / 3.0) * ch[(W - 1) * N + k] + (1.0 / 3.0) * ch[(W - 2) * N + k];
				for (j = 1; j < W - 1; j++) elution[j] = 0.5 * ch[j * N + k] + 0.25 * (ch[(j  -1) * N + k] + ch[(j + 1) * N + k]);

				bool skip_flag = false, inc_flag = false;
				float curr_evidence = elution[S];
				for (pos = S - Max(S / 3, 1); pos <= S + Max(S / 3, 1); pos++) if (elution[pos] > curr_evidence) {
					skip_flag = true;
					break;
				}
				if (skip_flag) continue;

				float max_evidence = curr_evidence, margin = curr_evidence / 3.0;
				for (pos = k - S + 1; pos <= k + S - 1; pos++) {
					int ind = pos - k + S;
					float evidence = elution[ind];
					if (evidence > max_evidence) max_evidence = evidence;
					if (evidence < margin) inc_flag = true;
				}
				if (!inc_flag) continue;
				if (curr_evidence < max_evidence * PeakApexEvidence) continue;

				features[i].push_back(top[iN + k]);
			}
		}
	}

	void feature_spectra(std::vector<Scan> &spectra, std::vector<std::vector<Peak> > &features, int N = 100, 
		float accuracy = 20.0 / 1000000.0, int min_fragments = 4, float corr_threshold = 0.8) {

	}

	XIC ms1_XIC(float mz, float RT, int S) {
		XIC xic;
		xic.pr_mz = mz, xic.RT = RT, xic.level = 1;
		if (!ms1h.size()) return xic;

		Scan query; query.RT = RT;
		int j, W = 2 * S + 1, n = ms1h.size(), index = std::distance(ms1h.begin(), std::lower_bound(ms1h.begin(), ms1h.end(), query));
		if (index >= n) index--;
		if (index > 0) if (RT - ms1h[index - 1].RT < ms1h[index].RT - RT) index--;

		xic.peaks.resize(W);
		float query_mz = predicted_mz_ms1(&(MassCorrectionMs1[0]), mz, RT);

		int start = Max(index - S, 0);
		for (int j = start; j <= Min(index + S, n - 1); j++) xic.peaks[j - start].first = ms1h[j].RT, xic.peaks[j - start].second = ms1h[j].level<false>(mz, MassAccuracyMs1);
		return xic;
	}

	XIC ms2_XIC(float pr_mz, float fr_mz, int apex, int S) {
		XIC xic;
		int cnt, j, W = 2 * S + 1, n = ms2h.size();

		xic.pr_mz = pr_mz, xic.fr_mz = fr_mz, xic.RT = scan_RT[apex], xic.level = 2;
		xic.peaks.resize(W);

		float query_mz = predicted_mz(&(MassCorrection[0]), fr_mz, scan_RT[apex]);
		for (j = apex, cnt = 0; j >= 0; j--) {
			if (cnt >= S + 1) break;
			if (ms2h[j].has(pr_mz)) xic.peaks[S - cnt].first = scan_RT[j], xic.peaks[S - cnt].second = ms2h[j].level<false>(query_mz, MassAccuracy), cnt++;
		}
		for (j = apex + 1, cnt = 0; j < n; j++) {
			if (cnt >= S) break;
			if (ms2h[j].has(pr_mz)) xic.peaks[S + cnt + 1].first = scan_RT[j], xic.peaks[S + cnt + 1].second = ms2h[j].level<false>(query_mz, MassAccuracy), cnt++;
		}
		return xic;
	}

	class Search {
	public:
		std::vector<Elution> peaks;
		std::vector<float> scoring;
		float scores[pN];
		float nn_sc[2];

		float best_corr_sum = 0.0, total_corr_sum = 0.0;

		std::vector<Fragment> quant;
		float RT, RT_start, RT_stop, evidence, ms1_corr, cscore, quantity, qvalue = 1.0, protein_qvalue = 0.0, best_fr_mz, q1_shift = 0.0;
		int apex, peak_width, best_peak, nn_index, nn_inc;
		int peak_pos, best_fragment;
		float mass_delta = 0.0, mass_delta_mz = 0.0, mass_delta_ms1 = 0.0;

#if (HASH > 0)
		unsigned int hash() {
			unsigned int res = 0;
			res ^= hashA((int*)scores, pN);
			return res;
		}
#endif
	};

    class Precursor {
    public:
		Lock lock;
        Run * run;
		Peptide * pep;

		int flags = 0, S, thread_id = 0;

		std::vector<Search> info; // used only for precursors detected in the sample
#if ELUTION_PROFILE
		std::vector<std::pair<float, float> > ms1_elution_profile;
#endif

		void init(Run * _run, Peptide * p) {
			run = _run;
			pep = p;
		}

		void build_index() {
			int i, k, n = run->ms2h.size();

			float mz = pep->mz, qmz = pep->mz, dmz = run->Q1Correction[0] + run->Q1Correction[1] * qmz;
			if (Abs(dmz) < Min(Abs(run->tandem_max - mz), Abs(mz - run->tandem_min)) * 0.5) qmz += dmz;

			auto scan_index = &(run->ScanIndex[thread_id]);
			auto chrom_index = &(run->ChromIndex[thread_id]);
			for (i = 0, k = 0; i < n; i++) if (run->ms2h[i].has(qmz)) k++;
			scan_index->resize(k);
			for (i = 0, k = 0; i < n; i++) if (run->ms2h[i].has(qmz)) (*scan_index)[k++] = i;

			if (!QuantOnly) {
				if (WindowRadius) S = WindowRadius;
				else if (run->PeakWidth > E) S = Max(1, int(ScanScale * run->PeakWidth));
				else S = Max(1, k / ScanFactor);
			}
		}

		class Searcher {
		public:
			Run * run;
			Precursor * pr;
			Peptide * pep;

			int n, m, S, W;

			float *nf, *ms1, *ms2, *ms2_H, *ms2_e, *ms2_min, *ms2_tight_one, *ms2_tight_two, *shadow, mz, qmz;
			std::vector<bool> * covered;
			std::vector<Product> * fragments;
			std::vector<int> *scan_index, *chrom_index, *fri;
			std::vector<Search> * info;

			Searcher(Precursor * precursor) {
				pr = precursor;
				run = pr->run;
				pep = pr->pep;

				auto &le = run->lib->entries[pep->index];
				if ((le.entry_flags & fFromFasta) && !le.target.fragments.size()) le.generate();

				mz = pep->mz;
				n = run->ms2h.size();
				S = pr->S, W = 2 * S + 1;
				scan_index = &(run->ScanIndex[pr->thread_id]);
				chrom_index = &(run->ChromIndex[pr->thread_id]);
 				info = &(pr->info);

				fragments = &(run->FR[pr->thread_id]);
				fragments->clear();
				fri = &(run->FRI[pr->thread_id]);
				fri->clear();
				int i = 0;
				for (auto &fr : pr->pep->fragments) {
					if (UseIsotopes) assert(fr.charge > 0);
					if (fr.mz > run->tandem_min && fr.mz < run->tandem_max) fragments->push_back(fr), fri->push_back(i);
					i++;
				}
				if (fragments->size() < MinSearchFrNum) fragments->clear(), fri->clear();
				m = fragments->size();
			}

			void chromatogram(int from = 0, int to = 0) { // extract full chromatograms for the top 6 fragments
				if (!m) return;
				if (!to) to = (*scan_index).size();
				int i, k, fr, l = to - from, ind, ms1_i = 1;
				int get_ms1 = (MS1PeakSelection && (run->ms1h.size() != 0)), force_ms1 = (ForceMs1 && (run->ms1h.size() != 0)), use_ms1 = get_ms1 || force_ms1;

				run->MS2[pr->thread_id].resize(l * m);
				ms2 = &(run->MS2[pr->thread_id][0]);
				for (i = 0; i < l * m; i++) ms2[i] = 0.0;
				if (get_ms1) {
					run->MS1[pr->thread_id].resize(l);
					ms1 = &(run->MS1[pr->thread_id][0]);
					for (i = 0; i < l; i++) ms1[i] = 0.0;
				}

				double pred_RT, RT_min, RT_max;
				if (run->RT_windowed_search) 
					pred_RT = calc_spline(run->RT_coeff, run->RT_points, pep->iRT),
					RT_min = pred_RT - run->RT_window, RT_max = pred_RT + run->RT_window;

				for (k = from; k < to; k++) {
					i = (*scan_index)[k];
					ind = k - from;
					float rt = run->scan_RT[i];
					if (use_ms1) {
						while (ms1_i < run->ms1_RT.size()) {
							if (run->ms1_RT[ms1_i] >= rt) break;
							ms1_i++;
						}
						if (ms1_i < run->ms1h.size()) {
							double query_mz = run->predicted_mz_ms1(&(run->MassCorrectionMs1[0]), mz, run->scan_RT[i]);
							double ms1_left = run->ms1h[ms1_i - 1].level<false>(query_mz, run->MassAccuracyMs1);
							double ms1_right = run->ms1h[ms1_i].level<false>(query_mz, run->MassAccuracyMs1);
							if (force_ms1) if (ms1_left < MinPeakHeight && ms1_right < MinPeakHeight) continue;
							if (get_ms1) ms1[ind] = ((rt - run->ms1h[ms1_i - 1].RT) * ms1_right + (run->ms1h[ms1_i].RT - rt) * ms1_left) / Max(E, run->ms1h[ms1_i].RT - run->ms1h[ms1_i - 1].RT);
						}
					}

					if (run->RT_windowed_search) { // main search phase only
						if (run->scan_RT[(*scan_index)[Min(to - 1, k + W)]] < RT_min) continue;
						if (run->scan_RT[(*scan_index)[Max(from, k - W)]] > RT_max) break;
					}

					for (fr = 0; fr < Min(m, TopF); fr++) {
						float query_mz = run->predicted_mz(&(run->MassCorrection[0]), (*fragments)[fr].mz, run->scan_RT[i]);
						int low = 0, high = run->ms2h[i].size();
						ms2[l * fr + ind] = run->ms2h[i].level<false>(query_mz, run->MassAccuracy);
						if (ms2[l * fr + ind] < MinPeakHeight) ms2[l * fr + ind] = 0.0;
					}
				}
			}

			void quant_chromatogram(int from, int to) { // extract chromatogram in the vicinity of the peak - for the quantify() function
				if (!m) return;
				int i, k, fr, l = to - from, ind, pos, center = (from + to) / 2;
				float ms1_left = 0.0, ms1_right = 0.0, ms1_left_mz = 0.0, ms1_right_mz = 0.0, rt, ms1_left_RT = 0.0, ms1_right_RT = 0.0, peak_mz = 0.0;

				run->MS2[pr->thread_id].resize(l * m);
				ms2 = &(run->MS2[pr->thread_id][0]);
				for (i = 0; i < l * m; i++) ms2[i] = 0.0;

				k = center, i = (*scan_index)[k];
				rt = run->scan_RT[i];

				auto ms1_ptr = std::lower_bound(run->ms1_RT.begin(), run->ms1_RT.end(), rt);
				for (pos = std::distance(run->ms1_RT.begin(), ms1_ptr); pos < run->ms1_RT.size()
					&& (run->ms1h[pos].window_low >= mz || run->ms1h[pos].window_high <= mz); pos++);
				if (pos == run->ms1_RT.size()) ms1_right_RT = INF;
				else {
					ms1_right_RT = run->ms1_RT[pos];
					double query_mz = run->predicted_mz_ms1(&(run->MassCorrectionMs1[0]), mz, rt);
					run->ms1h[pos].level<true>(query_mz, run->MassAccuracyMs1, &ms1_right_mz);
				}
				for (pos = std::distance(run->ms1_RT.begin(), ms1_ptr); pos >= 0
					&& (run->ms1h[pos].window_low >= mz || run->ms1h[pos].window_high <= mz || run->ms1_RT[pos] > rt); pos--);
				if (pos == -1) ms1_left_RT = 0.0;
				else {
					ms1_left_RT = run->ms1_RT[pos];
					double query_mz = run->predicted_mz_ms1(&(run->MassCorrectionMs1[0]), mz, rt);
					run->ms1h[pos].level<true>(query_mz, run->MassAccuracyMs1, &ms1_left_mz);
				}

				for (k = from; k < to; k++) {
					i = (*scan_index)[k];
					ind = k - from;

					for (fr = 0; fr < m; fr++) {
						float query_mz = run->predicted_mz(&(run->MassCorrection[0]), (*fragments)[fr].mz, run->scan_RT[i]);
						ms2[l * fr + ind] = run->ms2h[i].level<true>(query_mz, run->MassAccuracy, &peak_mz);
						if (ms2[l * fr + ind] < MinPeakHeight) ms2[l * fr + ind] = 0.0;

						if (k == center && info->size() && fr == (*info)[0].best_fragment && (!pr->pep->no_cal || !MassCalFilter)) {
							(*info)[0].mass_delta = ms2[l * fr + ind] >= MinPeakHeight ? (peak_mz - (*fragments)[fr].mz) : 0.0;
							(*info)[0].mass_delta_mz = (*fragments)[fr].mz;
							rt = run->scan_RT[i];
							if ((*info)[0].scores[pMs1TimeCorr] >= MassCalMs1Corr && ms1_left_mz > E && ms1_right_mz > E) {
								if (ms1_right_RT - ms1_left_RT > E) (*info)[0].mass_delta_ms1 = ((rt - ms1_left_RT) * ms1_right_mz + (ms1_right_RT - rt) * ms1_left_mz) / (ms1_right_RT - ms1_left_RT) - mz;
								else (*info)[0].mass_delta_ms1 = 0.5 * (ms1_right_mz + ms1_left_mz) - mz;
							} else (*info)[0].mass_delta_ms1 = 0.0;
						}
					}
				}
			}

			void extra_chromatogram(bool full_spectrum, int from = 0, int to = 0) { // obtain extra information in the vicinity of the putative elution peaks
				if (!m || !pr->info.size()) return;
				if (!to) to = scan_index->size();
				int i, k, fr, l = to - from, pos, ind;
				float ms1_left[9], ms1_right[9], rt, ms1_left_RT = 0.0, ms1_right_RT = 0.0;
				for (i = 0; i < 9; i++) ms1_left[i] = ms1_right[i] = 0.0;
				bool iso = full_spectrum, get_ms1 = true;

				run->Covered[pr->thread_id].resize(l);
				covered = &(run->Covered[pr->thread_id]);
				for (auto &p : (*info)[0].peaks)
					for (k = p.peak - S; k <= p.peak + S; k++) (*covered)[k] = false;

				run->MS2_min[pr->thread_id].resize(l * m);
				if (iso) run->MS2_H[pr->thread_id].resize(l * m);
				run->Shadow[pr->thread_id].resize(l * m);
				run->NF[pr->thread_id].resize(l);
				if (TightMassAcc) run->MS2_tight_one[pr->thread_id].resize(l * m), run->MS2_tight_two[pr->thread_id].resize(l * m);
				if (!run->ms1h.size()) get_ms1 = false;
				if (get_ms1) run->MS1[pr->thread_id].resize(l * 9);
				nf = &(run->NF[pr->thread_id][0]);
				if (get_ms1) ms1 = &(run->MS1[pr->thread_id][0]);
				ms2 = &(run->MS2[pr->thread_id][0]);
				ms2_min = &(run->MS2_min[pr->thread_id][0]);
				if (iso) ms2_H = &(run->MS2_H[pr->thread_id][0]);
				shadow = &(run->Shadow[pr->thread_id][0]);
				if (TightMassAcc) ms2_tight_one = &(run->MS2_tight_one[pr->thread_id][0]), ms2_tight_two = &(run->MS2_tight_two[pr->thread_id][0]);
				for (i = 0; i < l * m; i++) ms2_min[i] = shadow[i] = 0.0;
				if (TightMassAcc) for (i = 0; i < l * m; i++) ms2_tight_one[i] = ms2_tight_two[i] = 0.0;
				if (iso) for (i = 0; i < l * m; i++) ms2_H[i] = 0.0;
				if (get_ms1) for (i = 0; i < l * 9; i++) ms1[i] = 0.0;
				for (i = 0; i < l; i++) nf[i] = 0.0;

				for (auto &p : (*info)[0].peaks) {
					for (k = p.peak - S; k <= p.peak + S; k++) {
						if ((*covered)[k]) continue;

						i = (*scan_index)[k];
						ind = k - from;

						if (run->RT_windowed_search) { // main search phase only
							double margin = run->RT_window;
							double pred_RT = calc_spline(run->RT_coeff, run->RT_points, pep->iRT);
							if (run->scan_RT[(*scan_index)[Min(to - 1, k + W)]] < pred_RT - margin) continue;
							if (run->scan_RT[(*scan_index)[Max(from, k - W)]] > pred_RT + margin) break;
						}

						if (get_ms1) {
							rt = run->scan_RT[i];
							if (rt > ms1_right_RT) {
								auto ms1_ptr = std::lower_bound(run->ms1_RT.begin(), run->ms1_RT.end(), rt);
								for (pos = std::distance(run->ms1_RT.begin(), ms1_ptr); pos < run->ms1_RT.size()
									&& (run->ms1h[pos].window_low >= mz || run->ms1h[pos].window_high <= mz); pos++);
								if (pos == run->ms1_RT.size()) ms1_right_RT = INF, ms1[ind] = ms1_right[0], ms1[ind + l] = ms1_right[1], ms1[ind + 2 * l] = ms1_right[2];
								else {
									for (int j = 0; j < 9; j++) ms1_left[j] = ms1_right[j];
									ms1_left_RT = ms1_right_RT;
									ms1_right_RT = *ms1_ptr;

									double query_mz = run->predicted_mz_ms1(&(run->MassCorrectionMs1[0]), mz, run->scan_RT[i]);
									int low = 0, high = run->ms1h[pos].size();
									ms1_right[0] = run->ms1h[pos].level<false>(low, high, query_mz, run->MassAccuracyMs1);
									if (ms1_right[0] >= MinPeakHeight) {
										ms1_right[1] = run->ms1h[pos].level<false>(low, high, query_mz, run->MassAccuracyMs1 * TightMassAccRatioOne);
										ms1_right[2] = run->ms1h[pos].level<false>(low, high, query_mz, run->MassAccuracyMs1 * TightMassAccRatioTwo);
									} else ms1_right[1] = ms1_right[2] = 0.0;

									if (UseIsotopes) for (int isotope = 1; isotope < 3; isotope++) {
										int ii = 3 * isotope;
										query_mz = run->predicted_mz_ms1(&(run->MassCorrectionMs1[0]), mz + (C13delta * (double)isotope) / (double)pr->pep->charge, run->scan_RT[i]);
										low = 0, high = run->ms1h[pos].size();
										ms1_right[ii] = run->ms1h[pos].level<false>(low, high, query_mz, run->MassAccuracyMs1);
										if (ms1_right[ii] >= MinPeakHeight) {
											ms1_right[ii + 1] = run->ms1h[pos].level<false>(low, high, query_mz, run->MassAccuracyMs1 * TightMassAccRatioOne);
											ms1_right[ii + 2] = run->ms1h[pos].level<false>(low, high, query_mz, run->MassAccuracyMs1 * TightMassAccRatioTwo);
										} else ms1_right[ii + 1] = ms1_right[ii + 2] = 0.0;
									}
								}
							}
							for (int p = 0; p < 9; p++) ms1[ind + p * l] = ((rt - ms1_left_RT) * ms1_right[p] + (ms1_right_RT - rt) * ms1_left[p]) / Max(E, ms1_right_RT - ms1_left_RT);
						}

						for (fr = 0; fr < (full_spectrum ? m : Min(m, TopF)); fr++) {
							float query_mz = run->predicted_mz(&(run->MassCorrection[0]), (*fragments)[fr].mz, run->scan_RT[i]);
							int low = 0, high = run->ms2h[i].size();
							if (fr >= TopF) ms2[l * fr + ind] = run->ms2h[i].level<false>(low, high, query_mz, run->MassAccuracy);
							if (ms2[l * fr + ind] >= MinPeakHeight && TightMassAcc) {
								ms2_tight_one[l * fr + ind] = run->ms2h[i].level<false>(low, high, query_mz, run->MassAccuracy * TightMassAccRatioOne);
								ms2_tight_two[l * fr + ind] = run->ms2h[i].level<false>(low, high, query_mz, run->MassAccuracy * TightMassAccRatioTwo);
							}

							// shadow
							if (iso && ms2[l * fr + ind] >= MinPeakHeight) {
								query_mz = run->predicted_mz(&(run->MassCorrection[0]), (*fragments)[fr].mz - C13delta, run->scan_RT[i]); // assume charge 1 for the interfering fragment
								shadow[l * fr + ind] = run->ms2h[i].level<false>(query_mz, run->MassAccuracy);

								float hmz = mz + (C13delta / (double)pr->pep->charge) + run->Q1Correction[0] + run->Q1Correction[1] * mz;
								query_mz = run->predicted_mz(&(run->MassCorrection[0]), (*fragments)[fr].mz + C13delta / (double)(*fragments)[fr].charge, run->scan_RT[i]);
								int hi = -1;
								if (run->ms2h[i].has(hmz)) hi = i;
								if (i < run->ms2h.size() - 1 && hi < 0) if (run->ms2h[i + 1].has(hmz)) hi = i + 1;
								if (i > 0 && hi < 0) if (run->ms2h[i - 1].has(hmz)) hi = i - 1;
								if (hi >= 0) ms2_H[l * fr + ind] = run->ms2h[hi].level<false>(query_mz, run->MassAccuracy);
							}
						}

						// no fragmentation
						if (full_spectrum) {
							float query_mz = run->predicted_mz(&(run->MassCorrection[0]), mz, run->scan_RT[i]);
							int low = 0, high = run->ms2h[i].size();
							nf[ind] = run->ms2h[i].level<false>(low, high, query_mz, run->MassAccuracy);
						}
					}

					for (k = p.peak - S; k <= p.peak + S; k++) {
						if ((*covered)[k]) continue;
						else (*covered)[k] = true;

						for (fr = 0; fr < m; fr++)
							ms2_min[l * fr + k] = Min(Min(ms2[l * fr + k - 1], ms2[l * fr + k + 1]), ms2[l * fr + k]);
					}
				}
			}

#if ELUTION_PROFILE
			void get_ms1_elution_profile(std::vector<std::pair<float, float> > &ms1_elution_profile, float rt) {
				ms1_elution_profile.clear(); ms1_elution_profile.resize(2 * ElutionProfileRadius + 1, std::pair<float, float>(0.0, 0.0));
				if (!run->ms1h.size()) return;

				int pos, low = 0, high = run->ms1h.size();
				while (true) {
					pos = (low + high) >> 1;
					if (high <= low + 1) break;
					if (run->ms1h[pos].RT > rt) high = pos;
					else low = pos;
				}
				if (pos < run->ms1h.size() - 1) if (rt - run->ms1h[pos].RT > run->ms1h[pos + 1].RT - rt) pos++;
				int step = Max(1, S / ElutionProfileRadius);

				for (int i = 0; i <= ElutionProfileRadius; i++) {
					int ind = pos + i * step;
					if (ind >= run->ms1h.size()) break;
					double query_mz = run->predicted_mz_ms1(&(run->MassCorrectionMs1[0]), mz, run->ms1h[ind].RT);
					float height = run->ms1h[ind].level<false>(query_mz, run->MassAccuracyMs1);
					ms1_elution_profile[i + ElutionProfileRadius] = std::pair<float, float>(run->ms1h[ind].RT, height);
				}
				for (int i = 1; i <= ElutionProfileRadius; i++) {
					int ind = pos - i * step;
					if (ind < 0) break;
					double query_mz = run->predicted_mz_ms1(&(run->MassCorrectionMs1[0]), mz, run->ms1h[ind].RT);
					float height = run->ms1h[ind].level<false>(query_mz, run->MassAccuracyMs1);
					ms1_elution_profile[ElutionProfileRadius - i] = std::pair<float, float>(run->ms1h[ind].RT, height);
				}
			}
#endif

			noninline void peaks(int from = 0, int to = 0, float bcs = -INF, float tcs = -INF) {
				pr->info.clear();
				if (!m) return;
				if (!to) to = scan_index->size();
				int i, j, k, fr, added, pos, next, best_fr, next_fr, l = to - from, best_n = 0, n_peaks = 0, p = Min(m, TopF);
				float max, best = -INF, s;
				double total_corr_sum;
				double * corr_matrix = (double*)alloca(p * p * sizeof(double));
				float * elution = (float*)alloca(W * sizeof(float)), *score = (float*)alloca(p * sizeof(float));
				int * index = (int*)alloca(p * sizeof(int));

				run->PeakList[pr->thread_id].resize(l);
				run->BestFrList[pr->thread_id].resize(l);
				run->CorrSumList[pr->thread_id].resize(l);
				int * peak_list = &(run->PeakList[pr->thread_id][0]);
				int * best_fr_list = &(run->BestFrList[pr->thread_id][0]);
				float * corr_sum_list = &(run->CorrSumList[pr->thread_id][0]);

				if (p >= 2) for (k = from + S + 1, total_corr_sum = 0.0; k < to - S - 1; k++) {
					int ind = k - from;
					best_fr = next_fr = -1;
					for (fr = i = j = 0; fr < p; fr++) {
						bool flag = false;
						if (ms2[fr * l + ind] >= MinPeakHeight) i++, flag = true;
						else if (ms2[fr * l + ind - 1] >= MinPeakHeight || ms2[fr * l + ind + 1] >= MinPeakHeight) flag = true;
						if (flag) {
							if (!j) best_fr = fr;
							else if (j == 1) next_fr = fr;
							j++;
						}
					}
					if (i == 0 || j < 2) continue;

					for (fr = 0; fr < p; fr++) score[fr] = -1.0, index[fr] = -1;
					if (j == 2) {
						max = score[0] = score[1] = corr(&(ms2[best_fr * l + ind - S]), &(ms2[next_fr * l + ind - S]), W);
						if (max <= E) continue;
						double s1 = sum(&(ms2[best_fr * l + ind - S]), W);
						double s2 = sum(&(ms2[next_fr * l + ind - S]), W);
						if (s2 > s1) index[1] = best_fr, index[0] = best_fr = next_fr;
						else index[0] = best_fr, index[1] = next_fr;
					} else {
						best_fr = -1;
						get_corr_matrix(corr_matrix, &(ms2[ind - S]), p, l, W);
						for (fr = added = 0, max = E; fr < p; fr++) {
							s = 0.0;
							for (next = 0; next < p; next++) if (next != fr) s += corr_matrix[fr * p + next];
							if (s > E) {
								added++;
								if (s > max) max = s, best_fr = fr;
								for (int f = 0; f < p; f++) {
									if (index[f] < 0) { index[f] = fr; score[f] = s; break; }
									if (s > score[f]) {
										for (int g = added - 1; g > f; g--) index[g] = index[g - 1], score[g] = score[g - 1];
										index[f] = fr; score[f] = s; break;
									}
								}
							}
						}
					}

					bool found = false;
					for (next = 0; next < p; next++) {
						best_fr = index[next], max = score[next];
						if (best_fr < 0 || max <= E) break;

						if (MS1PeakSelection && run->ms1h.size()) {
							double ms1_corr = corr(&(ms2[best_fr * l + ind - S]), &(ms1[ind - S]), W);
							if (ms1_corr < MinMs1Corr) continue;
							max += ms1_corr;
						}
						if (max < MinCorrScore) continue;

						// smooth the curve
						smooth(elution, &(ms2[best_fr * l + ind - S]), W);
						if (elution[S] < MinPeakHeight) continue;

						// check if k is at the apex
						bool skip_flag = false;
						double curr_evidence = elution[S];
						for (pos = S - Max(S / 3, 1); pos <= S + Max(S / 3, 1); pos++) if (elution[pos] > curr_evidence) {
							skip_flag = true;
							break;
						}
						if (skip_flag) continue;

						double max_evidence = curr_evidence;
						int R = S - 1;
						for (pos = k - R; pos <= k + R; pos++) {
							double evidence = elution[pos - k + S];
							if (evidence > max_evidence) max_evidence = evidence;
						}
						if (curr_evidence < max_evidence * PeakApexEvidence) continue;

						if (run->full_spectrum && ForceQ1Apex) {
							float query_mz = run->predicted_mz(&(run->MassCorrection[0]), (*fragments)[best_fr].mz, run->scan_RT[(*scan_index)[k]]);
							if (run->ms2h[(*scan_index)[k] - 1].level<false>(query_mz, run->MassAccuracy) <= ms2[best_fr * l + ind]
								&& run->ms2h[(*scan_index)[k] + 1].level<false>(query_mz, run->MassAccuracy) <= ms2[best_fr * l + ind]) {
							} else continue;
						}

						total_corr_sum += max;
						peak_list[n_peaks] = k, best_fr_list[n_peaks] = best_fr, corr_sum_list[n_peaks] = max;
						if (max > best) best = max, best_n = n_peaks;
						n_peaks++, found = true;
						break;
					}
				}
				if (n_peaks == 0) return;

				pr->info.resize(1);
				pr->info[0].total_corr_sum = total_corr_sum;
				pr->info[0].best_corr_sum = best;
				if (bcs > -INF / 2.0) best = pr->info[0].best_corr_sum = bcs, pr->info[0].total_corr_sum = tcs;

				double min_score = best - MaxCorrDiff;
				for (i = pos = 0; i < n_peaks; i++) if (corr_sum_list[i] >= min_score) pos++;
				pr->info[0].peaks.resize(pos);
				pr->info[0].scoring.resize(pos * pN, 0.0);
				for (i = pos = 0; i < n_peaks; i++) if (corr_sum_list[i] >= min_score) {
					pr->info[0].peaks[pos].peak = peak_list[i];
					pr->info[0].peaks[pos].apex = (*scan_index)[peak_list[i]];
					pr->info[0].peaks[pos].best_fragment = best_fr_list[i];
					pos++;
				}
			}

			noninline void score(int peak, int from = 0, int to = 0) {
				assert(peak >= 0);
				assert(info->size());
				Search * inf = &((*info)[0]);
				assert(peak < inf->peaks.size());
				if (!m) return;
				if (!to) to = scan_index->size();
					
				int k = inf->peaks[peak].peak, pos, fr, l = to - from, ind = k - from, i, best_fr = inf->peaks[peak].best_fragment, apex = inf->peaks[peak].apex;
				float *x = (float*)alloca(m * sizeof(float)), *ref = (float*)alloca(m * sizeof(float));
				float *elution = (float*)alloca(W * sizeof(float)), *sc = &(inf->scoring[peak * pN]);
				double w, weight;	

				assert(ind >= pr->S);
				assert(ind < l - pr->S);

				// elution curve
				smooth(elution, &(ms2[best_fr * l + ind - S]), W);
				for (i = 0; i < pN; i++) sc[i] = 0.0;

				// pRef
				for (i = 0; i < m; i++) ref[i] = (*fragments)[i].height;
				for (fr = 1; fr < Min(m, auxF); fr++)
					sc[pRef + fr - 1] = ref[fr];

				// peptide characteristics
				sc[pMz] = (mz - 600.0) * 0.002;
				sc[pCharge] = -2.0 + (float)pep->charge;
				sc[pLength] = (-10.0 + (float)pep->length) / 10.0;
				sc[pFrNum] = m;

#if AAS
				// amino acid composition & number of modifications
				auto &name = run->lib->entries[pr->pep->index].name;
				for (i = 0; i < name.size(); i++) {
					char symbol = name[i];
					if (symbol < 'A' || symbol > 'Z') {
						if (symbol != '(' && symbol != '[') continue;
						i = closing_bracket(name, symbol, i);
						sc[pMods] += 1.0;
						continue;
					}
					sc[pAAs + AA_index[symbol]] += 1.0;
				}
#endif

				// MS2 to MS1 signal ratio
				float heights[TopF];
				for (i = 0; i < TopF; i++) heights[i] = 0.0;
				for (fr = 0; fr < Min(m, TopF); fr++) {
					float h = ms2[fr * l + ind];
					for (pos = 0; pos < TopF; pos++) {
						if (h > heights[pos]) {
							for (int j = TopF - 1; j > pos; j--) heights[j] = heights[j - 1];
							heights[pos] = h;
							break;
						}
					}
				}
				for (i = 1; i < TopF; i++) sc[pTotal] += heights[i];
				if (sc[pTotal] < MinPeakHeight) sc[pTotal] += heights[0];
				sc[pTotal] = log(Max(E, sc[pTotal]));
				sc[pMs1Ratio] = sc[pTotal] - log(Max(E, ms1[ind]));

				// pCos
				for (pos = ind - S, w = 0.0; pos <= ind + S; pos++) {
					weight = Sqr(elution[pos - ind + S]);
					if (weight < E) continue;
					w += weight;

					for (fr = 0; fr < Min(m, TopF); fr++) x[fr] = ms2[fr * l + pos];
					double cv = cos(&(x[0]), &(ref[0]), Min(m, TopF));
					sc[pCos] += cv * weight;
					sc[pCosCube] += Cube(cv) * weight;
				}
				assert(w > E);
				sc[pCos] /= w;
				sc[pCosCube] /= w;

				// time corr
				int T = Max(S / 2, 1);
				for (fr = 0; fr < Min(m, auxF); fr++) {
					double tu, u = corr(elution, &(ms2[fr * l + ind - S]), W);
					if (fr < TopF) {
						if (TightMassAcc) {
							sc[pTightCorrOne] += (tu = corr(elution, &(ms2_tight_one[fr * l + ind - S]), W)); if (tu > u) u = tu;
							sc[pTightCorrTwo] += (tu = corr(elution, &(ms2_tight_two[fr * l + ind - S]), W)); if (tu > u) u = tu;
						}
						sc[pTimeCorr] += u;

						sc[pLocCorr] += corr(&(elution[S - T]), &(ms2[fr * l + ind - T]), 2 * T + 1);
						sc[pMinCorr] += corr(elution, &(ms2_min[fr * l + ind - S]), W);
						if (run->full_spectrum && UseIsotopes) {
							double v = corr(elution, &(shadow[fr * l + ind - S]), W);
							sc[pShadow] += v;
							sc[pShadowCorr + fr] += v;
							sc[pHeavy] += corr(elution, &(ms2_H[fr * l + ind - S]), W);
						}

						float peak_mz = 0.0, query_mz = run->predicted_mz(&(run->MassCorrection[0]), (*fragments)[fr].mz, run->scan_RT[i]);
						run->ms2h[inf->peaks[peak].apex].level<true>(query_mz, run->MassAccuracy, &peak_mz);
						if (Abs(peak_mz) > E) sc[pAcc + fr] = (Abs(peak_mz - query_mz) / ((*fragments)[fr].mz * run->MassAccuracy)) * u;
					}
					sc[pCorr + fr] += u;
				}

#if Q1
				if (UseQ1 && run->full_spectrum && (run->scanning || ForceQ1)) {
					int QS = QSL[qL - 1], QW = 2 * QS + 1;
					float *qprofile = (float*)alloca(QW * (Min(m, TopF) + 1) * sizeof(float)), *ws = (float*)alloca(QW * sizeof(float)), *fmz = (float*)alloca((Min(m, TopF) + 1) * sizeof(float));
					float rt = run->scan_RT[apex];
					for (i = 0; i < Min(m, TopF); i++) fmz[i] = run->predicted_mz(&(run->MassCorrection[0]), (*fragments)[i].mz, rt);
					fmz[Min(m, TopF)] = run->predicted_mz(&(run->MassCorrection[0]), mz, rt);
					run->q1_profiles(qprofile, ws, apex, fmz, QS, Min(m, TopF) + 1);
					for (i = 0; i < qL; i++) sc[pQPos + i] = centroid_coo(qprofile + best_fr * QW + QS - QSL[i], ws + QS - QSL[i], 2 * QSL[i] + 1) - mz - run->Q1Correction[0] - run->Q1Correction[1] * mz;
					for (fr = 0; fr < Min(m, TopF); fr++) {
						if (fr != best_fr)
							for (i = 0; i < qL; i++) sc[pQCorr + i] += corr(qprofile + fr * QW + QS - QSL[i], qprofile + best_fr * QW + QS - QSL[i], 2 * QSL[i] + 1);
						for (i = 0; i < qL; i++) {
							sc[pQLeft + i] += signal_level(qprofile + fr * QW + QS - QSL[i], QSL[i], QSL[i] + 1) * sc[pCorr + fr];
							sc[pQRight + i] += signal_level(qprofile + fr * QW + QS, 0, QSL[i] + 1) * sc[pCorr + fr];
						}
					}
					for (i = 0; i < qL; i++) sc[pQNFCorr + i] += corr(qprofile + Min(m, TopF) * QW + QS - QSL[i], qprofile + best_fr * QW + QS - QSL[i], 2 * QSL[i] + 1);
				}
#endif

				if (m > TopF) {
					for (fr = TopF; fr < m; fr++) {
						double u = corr(elution, &(ms2[fr * l + ind - S]), W);
						sc[pResCorr] += u;
					}
					sc[pResCorrNorm] = sc[pResCorr] / (double)(m - TopF);
				}

				sc[pNFCorr] = corr(elution, &(nf[ind - S]), W);
				sc[pMs1TimeCorr] += corr(elution, &(ms1[ind - S]), W);
				sc[pMs1TightOne] += corr(elution, &(ms1[ind - S + l]), W);
				sc[pMs1TightTwo] += corr(elution, &(ms1[ind - S + 2 * l]), W);
				if (UseIsotopes) for (int isotope = 1; isotope < 3; isotope++) {
					int ii = isotope * 3;
					sc[pMs1Iso] += corr(elution, &(ms1[ind - S + ii * l]), W);
					sc[pMs1IsoOne] += corr(elution, &(ms1[ind - S + (ii + 1) * l]), W);
					sc[pMs1IsoTwo] += corr(elution, &(ms1[ind - S + (ii + 2) * l]), W);
				}

				// signals
				for (fr = 0, w = 0.0; fr < Min(m, TopF); fr++) {
					for (pos = ind - S; pos <= ind + S; pos++) {
						w += ms2[fr * l + pos];
						sc[pSig + fr] += ms2[fr * l + pos];
					}
				}
				if (w > E) for (fr = 0; fr < TopF; fr++) sc[pSig + fr] /= w;

				// shape
				int nnD = 2 * (S / (nnS * 2)) + 1;
				for (pos = k - S, w = 0.0; pos <= k + S; pos++) {
					int dist = Min(nnS, (Abs(pos - k) + nnD / 2) / nnD);
					int bin = pos >= k ? (nnS + dist) : (nnS - dist);
					w += elution[pos - k + S];
					sc[pShape + bin] += elution[pos - k + S];
				}
				for (i = 0; i < nnW; i++) sc[pShape + i] /= w;

				// corr sums
				sc[pBestCorrDelta] = sc[pTimeCorr] - inf->best_corr_sum;
				sc[pTotCorrSum] = log(Max(1.0, sc[pTimeCorr]) / (inf->total_corr_sum + 1.0));

				// extra transform to improve the performance of the linear classifier
#if LOGSC
#define LogEnhance(x,cap) x -= log(Max(E, -x + ((double)cap)));
				LogEnhance(sc[pTimeCorr], TopF);
				LogEnhance(sc[pLocCorr], TopF);
				LogEnhance(sc[pMinCorr], TopF);
				LogEnhance(sc[pMs1TimeCorr], 1.0);
				LogEnhance(sc[pNFCorr], 1.0);
				LogEnhance(sc[pResCorrNorm], 1.0);
				LogEnhance(sc[pTightCorrOne], TopF);
				LogEnhance(sc[pTightCorrTwo], TopF);
				LogEnhance(sc[pHeavy], TopF);
				LogEnhance(sc[pMs1TightOne], 1.0);
				LogEnhance(sc[pMs1TightTwo], 1.0);
				LogEnhance(sc[pMs1Iso], 2.0);
				LogEnhance(sc[pMs1IsoOne], 2.0);
				LogEnhance(sc[pMs1IsoTwo], 2.0);
				LogEnhance(sc[pCos], 1.0);
				LogEnhance(sc[pCosCube], 1.0);
#endif
			}
		};

		void score_RT(int peak) {
			assert(peak >= 0);
			assert(info.size());
			assert(peak < info[0].peaks.size());
			assert(info[0].peaks[peak].apex >= 0);
			assert(info[0].peaks[peak].apex < run->scan_RT.size());

			if (DisableRT) return;
			double span = Max(E, run->RT_max - run->RT_min);
			double mes_RT = run->scan_RT[info[0].peaks[peak].apex];
			double pred_RT = calc_spline(run->RT_coeff, run->RT_points, pep->iRT);
			double pos = Min(pred_RT - run->RT_min, span);
			double delta = mes_RT - pred_RT;
			float *sc = &(info[0].scoring[peak * pN]);
			sc[pRT] = pos / span;
			sc[pdRT] = sqrt(Min(Abs(delta), span) / span);
		}

		void score_RT() {
			assert(info.size());
			if (!(flags & fFound)) return;
			if (DisableRT) return;
			double mes_RT = run->scan_RT[info[0].apex];
			double pred_RT = calc_spline(run->RT_coeff, run->RT_points, pep->iRT);
			double span = Max(E, run->RT_max - run->RT_min);
			double pos = Min(pred_RT - run->RT_min, span);
			double delta = mes_RT - pred_RT;
			float *sc = &(info[0].scores[0]);
			sc[pRT] = pos / span;
			sc[pdRT] = sqrt(Min(Abs(delta), span) / span);
		}

		void search() {
			if (!(flags & fInit)) {
				build_index();
				Searcher searcher(this);
				searcher.chromatogram();
				searcher.peaks();
				searcher.extra_chromatogram(run->full_spectrum);
				flags |= fInit;
				if (info.size()) for (int peak = 0; peak < info[0].peaks.size(); peak++) searcher.score(peak);
			}
			if (info.size()) for (int peak = 0; peak < info[0].peaks.size(); peak++) score_RT(peak);
		}

		void free() {
			flags &= ~fInit;
			if (info.size()) {
				std::vector<Elution>().swap(info[0].peaks);
				std::vector<float>().swap(info[0].scoring);
			}
		}

		noninline void quantify(int peak, int stage, bool filter = false) {
			assert(peak >= 0);

			if (!(flags & fFound)) return;
			if (flags & fTranslated) {
				if (stage < 3) return;
			} else if (stage >= 3) return;
			if (stage) {
				if (filter && info[0].qvalue > QFilter) return;
				auto &le = run->lib->entries[pep->index];
				if (!pep->fragments.size()) le.generate();
			}
			build_index();

			int fr, best_fr, pos, k = peak, start, stop, vpos, width;
            double w, r, e;

			Searcher searcher(this);
			int p = Min(TopF, searcher.m), R, low, high, len;
			float scaled[TopF], *elution;

			R = ((QuantMode == 1 || stage < 2) ? Max(1, S / 2) : Max(1, S));
			low = Max(0, k - R);
			high = Min(searcher.scan_index->size() - 1, k + R);
			len = high - low + 1;
			elution = (float*)alloca(len * sizeof(float));
			searcher.quant_chromatogram(low, high + 1);

			if (stage) info[0].quant.resize(p);
			best_fr = info[0].best_fragment;
			smooth(&(elution[0]), &(searcher.ms2[best_fr * len]), len);

			e = elution[k - low] * 0.5;
			for (pos = info[0].peak_width = 0; pos < len; pos++) if (elution[pos] >= e)
				info[0].peak_width++;
			if (flags & fTranslated) { // translate peak boundaries
				int eg = run->lib->elution_groups[pep->index];
				if (run->best_egs[eg].second <= -INF) return;
				int tr_index = run->best_egs[eg].first;
				auto &tr_pr = run->entries[tr_index].target;
				assert(tr_pr.pep->index != pep->index);

				for (start = k - low - 1; start > 0; start--) if (run->ms2h[(*searcher.scan_index)[start]].RT < tr_pr.info[0].RT_start - E) {
					start = Min(k - low - 1, start + 1);
					break;
				}
				for (stop = k - low + 1; stop < high - low; stop++) if (run->ms2h[(*searcher.scan_index)[stop]].RT > tr_pr.info[0].RT_stop + E) {
					stop = Max(k - low + 1, stop - 1);
					break;
				}
			} else {
				start = 0, stop = len - 1;
				if (QuantMode == 0 && stage >= 2) {
					float valley, boundary = elution[k - low] / PeakBoundary;
					for (start = vpos = k - low - 1, valley = elution[k - low]; start > 0; start--) {
						if (elution[start] < valley) valley = elution[start], vpos = start;
						else if (valley < elution[k - low] / 3.0 && valley < elution[start] / 2.0) { start = vpos;  break; }
						if (Max(searcher.ms2[best_fr * len + start], searcher.ms2[best_fr * len + start + 1]) < boundary) break;
					}
					for (stop = vpos = k - low + 1, valley = elution[k - low]; stop < high - low; stop++) {
						if (elution[stop] < valley) valley = elution[stop], vpos = stop;
						else if (valley < elution[k - low] / 3.0 && valley < elution[stop] / 2.0) { stop = vpos;  break; }
						if (Max(searcher.ms2[best_fr * len + stop], searcher.ms2[best_fr * len + stop - 1]) < boundary) break;
					}
				}
			}
			width = stop - start + 1;
			info[0].RT_start = run->ms2h[(*searcher.scan_index)[start]].RT;
			info[0].RT_stop = run->ms2h[(*searcher.scan_index)[stop]].RT;
#if Q1
			float bmz = run->predicted_mz(&(run->MassCorrection[0]), info[0].best_fr_mz, info[0].RT);
			info[0].q1_shift = run->q1_delta(info[0].apex, pep->mz, bmz, QSL[qL - 1]);
#endif
			if (!stage) return;

			for (fr = 0; fr < p; fr++) {
				info[0].quant[fr].index = (*searcher.fri)[fr];
				info[0].quant[fr].corr = 0.0;
				for (int q = 0; q < qN; q++) info[0].quant[fr].quantity[q] = 0.0;
				scaled[fr] = 0.0;
			}
			for (pos = k - low - 1, w = 0.0; pos <= k - low + 1; pos++) w += Sqr(elution[pos]);
			assert(w > E);

			for (fr = 0; fr < p; fr++) {
				info[0].quant[fr].corr = corr(&(elution[start]), &(searcher.ms2[fr * len + start]), width);
				info[0].quant[fr].quantity[qTotal] = sum(&(searcher.ms2[fr * len + start]), width);
				if (NoIfsRemoval) info[0].quant[fr].quantity[qFiltered] = info[0].quant[fr].quantity[qTotal];
				else {
					for (pos = k - low - 1; pos <= k - low + 1; pos++)
						scaled[fr] += searcher.ms2[fr * len + pos] * Sqr(elution[pos]);
					scaled[fr] /= w;
				}
			}

			if (!NoIfsRemoval) {
				assert(best_fr < p);
				assert(scaled[best_fr] > E);

				info[0].quant[best_fr].quantity[qFiltered] = info[0].quant[best_fr].quantity[qTotal];
				for (fr = 0; fr < p; fr++) {
					if (fr == best_fr) continue;
					r = scaled[fr] / scaled[best_fr];
					if (info[0].quant[fr].corr < 0.8) {
						r = Min(r, searcher.ms2[fr * len + k - low] / Max(E, searcher.ms2[best_fr * len + k - low]));

						double radj = r;
						e = elution[k - low] * 0.5;
						for (pos = k - low - 1; pos >= 0; pos--) {
							if (elution[pos] < e) break;
							else radj = Min(radj, searcher.ms2[fr * len + pos] / Max(E, searcher.ms2[best_fr * len + pos]));
						}
						for (pos = k - low + 1; pos < len; pos++) {
							if (elution[pos] < e) break;
							else radj = Min(radj, searcher.ms2[fr * len + pos] / Max(E, searcher.ms2[best_fr * len + pos]));
						}

						r = Max(r * 0.5, Min(r, radj));
					}

					float query_mz = run->predicted_mz(&(run->MassCorrection[0]), searcher.fragments->at(fr).mz, run->ms2h[(*searcher.scan_index)[k]].RT);
					for (pos = start; pos <= stop; pos++) {
						float value;
						if (!QuantFitProfiles) value = searcher.ms2[fr * len + pos];
						else value = run->ms2h[(*searcher.scan_index)[low + pos]].level(query_mz, run->MassAccuracy, r * searcher.ms2[best_fr * len + pos]);
						info[0].quant[fr].quantity[qFiltered] += Min(value, FilterFactor * r * elution[pos]);
					}
				}
			}
			for (fr = 0, info[0].quantity = 0.0; fr < p; fr++) info[0].quantity += info[0].quant[fr].quantity[qFiltered];

#if ELUTION_PROFILE
			searcher.get_ms1_elution_profile(ms1_elution_profile, info[0].RT);
#endif

			if (info[0].qvalue <= 0.01 && (e = searcher.ms2[best_fr * len + k - low] * 0.5) >= MinPeakHeight) { // for run stats
				double fwhm_sc = 0.0, fwhm_rt = 0.0;
				int flag;
				for (pos = k, flag = 0; pos < high; pos++) {
					int i = pos - low;
					double mul, rtd = run->scan_RT[(*searcher.scan_index)[pos + 1]] - run->scan_RT[(*searcher.scan_index)[pos]];
					if (searcher.ms2[best_fr * len + i + 1] >= e) mul = 1.0;
					else {
						mul = (searcher.ms2[best_fr * len + i] - e) / Max(E, searcher.ms2[best_fr * len + i] - searcher.ms2[best_fr * len + i + 1]);
						flag = 1;
					}
					fwhm_sc += mul, fwhm_rt += rtd * mul;
					if (flag) break;
				}
				for (pos = k, flag = 0; pos > low; pos--) {
					int i = pos - low;
					double mul, rtd = run->scan_RT[(*searcher.scan_index)[pos]] - run->scan_RT[(*searcher.scan_index)[pos - 1]];
					if (searcher.ms2[best_fr * len + i - 1] >= e) mul = 1.0;
					else {
						mul = (searcher.ms2[best_fr * len + i] - e) / Max(E, searcher.ms2[best_fr * len + i] - searcher.ms2[best_fr * len + i - 1]);
						flag = 1;
					}
					fwhm_sc += mul, fwhm_rt += rtd * mul;
					if (flag) break;
				}
				run->RS.fwhm_scans += fwhm_sc;
				run->RS.fwhm_RT += fwhm_rt;

				run->RS.tot_quantity += info[0].quantity;

				run->RS.pep_length += pep->length;
				run->RS.pep_charge += pep->charge;
				run->RS.pep_mc += missed_KR_cleavages(run->lib->entries[pep->index].name);

				run->RS.cnt_valid++;
			}

			if (stage == 2 && Visualise.size()) {
				auto &le = run->lib->entries[pep->index];
				auto aas = get_aas(le.name);
				if (std::find(Visualise.begin(), Visualise.end(), aas) != Visualise.end()) {
					if (run->ms1h.size()) {
						auto xic = run->ms1_XIC(pep->mz, info[0].RT, VisWindowRadius); 
						xic.qvalue = info[0].qvalue, xic.run = run->run_index, xic.pr = pep->index; XICs.insert(xic);
					}
					for (auto& fr : pep->fragments) {
						auto xic = run->ms2_XIC(pep->mz, fr.mz, (*searcher.scan_index)[peak], VisWindowRadius); 
						xic.qvalue = info[0].qvalue, xic.run = run->run_index, xic.pr = pep->index;
						xic.fr = fr; XICs.insert(xic);
					}
				}
			}
        }

		noninline void intensities(int apex, std::vector<Ion> &gen) { // get fragment intensities
			int i, j, fr, cnt, low, high;
			double u, max, mz = pep->mz;

			std::vector<Product>().swap(pep->fragments);
			auto &le = run->lib->entries[pep->index];
			auto scan_index = &(run->ScanIndex[thread_id]);
			int S = Max(1, le.window / 2);
			scan_index->resize(2 * S + 1);
			(*scan_index)[S] = apex;

			float qmz = mz, dmz = run->Q1Correction[0] + run->Q1Correction[1] * mz;
			if (Abs(dmz) < Min(Abs(run->tandem_max - mz), Abs(mz - run->tandem_min)) * 0.5) qmz += dmz;
			for (high = S, i = apex + 1; i < run->ms2h.size(); i++) {
				if (run->ms2h[i].has(qmz)) (*scan_index)[++high] = i;
				if (high >= 2 * S) break;
			}
			for (low = S, i = apex - 1; i >= 0; i--) {
				if (run->ms2h[i].has(qmz)) (*scan_index)[--low] = i;
				if (low <= 0) break;
			}

			int l = high - low + 1;
			float bmz = le.best_fr_mz;
			float *x = (float*)alloca(l * sizeof(float)), *elution = (float*)alloca(l * sizeof(float));

			volatile float query_mz = run->predicted_mz(&(run->MassCorrection[0]), bmz, run->scan_RT[apex]);
			for (i = 0; i < l; i++) x[i] = run->ms2h[(*scan_index)[low + i]].level<false>(query_mz, run->MassAccuracy);
			smooth(elution, x, l);
			float best_height = elution[S - low];

			Peak * new_fr = (Peak*)alloca(gen.size() * sizeof(Peak));
			float * actual_mz = (float*)alloca(gen.size() * sizeof(float)), *err_mz = (float*)alloca(gen.size() * sizeof(float));;
			for (fr = cnt = 0, max = 0.0; fr < gen.size(); fr++) {
				float fmz = gen[fr].mz, amz = 0.0;
				if (fmz < MinFrMz || fmz > MaxFrMz) continue;
				if (fmz <= run->tandem_min || fmz >= run->tandem_max) continue;
				float query_mz = run->predicted_mz(&(run->MassCorrection[0]), fmz, run->scan_RT[apex]);
				x[S - low] = run->ms2h[apex].level<true>(query_mz, run->MassAccuracy, &amz);
				for (i = 0; i < l; i++) if (i != S - low) x[i] = run->ms2h[(*scan_index)[low + i]].level<false>(query_mz, run->MassAccuracy);
				u = x[S - low];
				if (u >= best_height * MinRelFrHeight && u >= MinGenFrInt) {
					float c = corr(x, elution, l);
					if (gen[fr].type & fTypeB) c -= FrCorrBSeriesPenalty;
					if (c >= MinGenFrCorr && (c >= MinRareFrCorr || ((gen[fr].charge == 1 || pep->charge >= 3) && gen[fr].loss == loss_none))) {
						bool skip = false;
						for (j = 0; j < cnt; j++) if (Abs(amz - actual_mz[j]) < E) {
							skip = true;
							if (Abs(amz - query_mz) < err_mz[j]) new_fr[j].mz = fmz, err_mz[j] = Abs(amz - query_mz);
						}
						if (skip) continue;
						new_fr[cnt].mz = fmz, new_fr[cnt].height = u;
						actual_mz[cnt] = amz, err_mz[cnt] = Abs(amz - query_mz);
						if (u > max) max = u;
						cnt++;
					}
				}
			}
			if (cnt < MinOutFrNum) run->lib->entries[pep->index].best_run = -1;
			else {
				pep->fragments.resize(cnt);
				for (fr = 0; fr < cnt; fr++) pep->fragments[fr].mz = new_fr[fr].mz, pep->fragments[fr].height = new_fr[fr].height / max;
			}
		}

        void seek() {
			if (!info.size()) return;

            int k, peak;
            double score;
            info[0].cscore = -INF, info[0].best_peak = -1;

			auto fragments = &(run->FR[thread_id]);
			fragments->clear();
			for (auto &fr : pep->fragments) if (fr.mz > run->tandem_min && fr.mz < run->tandem_max) (*fragments).push_back(fr);

			bool guide = GuideLibrary && GuideClassifier && !(run->lib->entries[pep->index].entry_flags & fFromFasta);
            for (peak = 0; peak < info[0].peaks.size(); peak++) {
                score = run->cscore<true>(&(info[0].scoring[peak * pN]), guide);
                if (score > info[0].cscore) {
					if (run->RT_windowed_search && Abs(run->scan_RT[info[0].peaks[peak].apex] - calc_spline(run->RT_coeff, run->RT_points, pep->iRT)) > run->RT_window) continue;
					info[0].best_peak = peak;
					info[0].cscore = score;
                }
            }
			bool found = (info[0].best_peak >= 0);
			if (found) {
				info[0].best_fragment = info[0].peaks[info[0].best_peak].best_fragment;
				info[0].best_fr_mz = (*fragments)[info[0].best_fragment].mz;
				info[0].apex = info[0].peaks[info[0].best_peak].apex;
				info[0].peak_pos = info[0].peaks[info[0].best_peak].peak;
				info[0].RT = run->scan_RT[info[0].apex];
				for (k = 0; k < pN; k++) info[0].scores[k] = info[0].scoring[info[0].best_peak * pN + k];
				info[0].evidence = info[0].scores[pTimeCorr];
				info[0].ms1_corr = info[0].scores[pMs1TimeCorr];
				flags |= fFound;
				if (!(flags & fDecoy) && run->curr_iter == CalibrationIter && (InferWindow || Calibrate))
					quantify(info[0].peaks[info[0].best_peak].peak, 0);
			} else flags &= ~fFound;
        }

        void find(int _thread_id, bool _free, bool skip) {
			if (lock.set()) {
				if (skip) {
					if (info.size()) score_RT();
					return;
				}
				thread_id = _thread_id;
				if (!QuantOnly && run->curr_iter < iN) {
					search();
					seek();
				} else if ((flags & fFound) && !(flags & fDecoy)) quantify(info[0].peak_pos, 2);
				if (_free) free();
			}
        }

		noninline void translate() {
			if (!(flags & fFound)) return;
			int eg = run->lib->elution_groups[pep->index];
			if (run->best_egs[eg].second <= -INF) return;
			int tr_index = run->best_egs[eg].first;
			auto &tr_pr = run->entries[tr_index].target;
			if (tr_pr.pep->index == pep->index) return;
			if (tr_pr.info[0].qvalue > RubbishFilter) return;

			build_index();
			Searcher searcher(this);
			int p = Min(TopF, searcher.m), S = searcher.S, W = searcher.W;
			if (searcher.scan_index->size() <= W + 3) return;

			int k, low = S + 1, high = searcher.scan_index->size() - S - 2;
			float target_RT = run->entries[tr_index].target.info[0].RT;
			while (high - low > 1) {
				k = (low + high) >> 1;
				if (run->scan_RT[searcher.scan_index->at(k)] <= target_RT) low = k;
				else high = k;
			}
			k = low;
			if (run->scan_RT[searcher.scan_index->at(k + 1)] - target_RT < target_RT - run->scan_RT[searcher.scan_index->at(k)]) k++;

			if (Abs(k - info[0].peak_pos) <= Max(S / 2, 1)) return;
			float bcs = info[0].best_corr_sum, tcs = info[0].total_corr_sum;

			int min = Max(0, k - S - 1), max = Min(k + S + 2, searcher.scan_index->size());
			searcher.chromatogram(min, max);
			info.clear();
			searcher.peaks(min, max, bcs, tcs);
			if (info.size()) {
				searcher.extra_chromatogram(run->full_spectrum, min, max);
				for (int peak = 0; peak < info[0].peaks.size(); peak++) searcher.score(peak, min, max), score_RT(peak);
				seek();
				flags |= fTranslated;
			}
		}

#if (HASH > 0)
		unsigned int hash() {
			unsigned int res = 0;
			if (flags & fFound) if (info.size()) res ^= info[0].hash();
			return res;
		}
#endif
    };

	class Target {
	public:
		Lock lock;
		Precursor target;
		Precursor decoy;
		int batch = 0;

		Target() { }

#if (HASH > 0)
		unsigned int hash() { return target.hash() ^ decoy.hash(); }
#endif
	};

	std::vector<Target> entries;

	void write(char * file_name) {
		std::ofstream out(file_name, std::ofstream::binary);
		if (out.fail()) { dsout << "ERROR: cannot write to " << file_name << ". Check if the destination folder is write-protected or the file is in use\n"; return; }

		int size = -2; // second version of the .dia file format
		out.write((char*)&size, sizeof(int));
		int size_ms1 = ms1h.size();
		out.write((char*)&size_ms1, sizeof(int));
		int size_ms2 = ms2h.size();
		out.write((char*)&size_ms2, sizeof(int));
		long long size_peaks = peaks.size();
		if (!size_peaks) {
			for (int i = 0; i < size_ms1; i++) size_peaks += ms1h[i].n;
			for (int i = 0; i < size_ms2; i++) size_peaks += ms2h[i].n;
		}
		out.write((char*)&size_peaks, sizeof(long long));

		std::vector<Scan> ms1hc(ms1h.size()), ms2hc(ms2h.size());
		memcpy((char*)&(ms1hc[0]), (char*)&(ms1h[0]), size_ms1 * sizeof(Scan));
		memcpy((char*)&(ms2hc[0]), (char*)&(ms2h[0]), size_ms2 * sizeof(Scan));
		long long pos = 0;
		for (int i = 0; i < ms1hc.size(); i++) ms1hc[i].peaks = (Peak*)pos, pos += ms1hc[i].n;
		for (int i = 0; i < ms2hc.size(); i++) ms2hc[i].peaks = (Peak*)pos, pos += ms2hc[i].n;
		
		out.write((char*)&(ms1hc[0]), size_ms1 * sizeof(Scan));
		out.write((char*)&(ms2hc[0]), size_ms2 * sizeof(Scan));
		for (int i = 0; i < size_ms1; i++) out.write((char*)ms1h[i].peaks, ms1h[i].n * sizeof(Peak));
		for (int i = 0; i < size_ms2; i++) out.write((char*)ms2h[i].peaks, ms2h[i].n * sizeof(Peak));

		out.close();
	}

	void read_scan(std::ifstream &in, Scan &h, std::vector<Peak> &s) {
		double RT = 0.0, window_low = 0.0, window_high = 0.0;

		in.read((char*)&RT, sizeof(double));
		in.read((char*)&window_low, sizeof(double));
		in.read((char*)&window_high, sizeof(double));
		h.RT = RT, h.window_low = window_low, h.window_high = window_high;

		int size = 0; in.read((char*)&size, sizeof(int));
		if (size) {
			s.resize(size);
			in.read((char*)&(s[0]), size * sizeof(Peak));
		}
	}

	bool read(char * file_name) {
		std::ifstream in(file_name, std::ifstream::binary);
		if (in.fail()) { dsout << "ERROR: cannot open " << file_name << "\n"; return false; }

		int size; in.read((char*)&size, sizeof(int));
		if (size == 0) { dsout << "ERROR: no MS2 spectra in the file\n"; return false; }
		if (size > 0) { // first version of the .dia format
			const int Block = 5;
			std::vector<std::vector<Peak> > S(Block);
			ms2h.resize(size);
			int i, j, tot;
			for (i = tot = 0; i < size; i++) {
				int ind = i % Block; auto &s = S[ind];
				read_scan(in, ms2h[i], s), ms2h[i].type = 2, tot += s.size();

				if (ind == Block - 1 || i == size - 1) {
					peak_lists.push_back(std::vector<Peak>(tot));
					auto &pl = peak_lists.back();
					for (j = tot = 0; j <= ind; j++) memcpy(&(pl[tot]), &(S[j][0]), S[j].size() * sizeof(Peak)), ms2h[i - ind + j].peaks = &(pl[tot]), ms2h[i - ind + j].n = S[j].size(), tot += S[j].size();
				}
			}

			in.read((char*)&size, sizeof(int));
			ms1h.resize(size);
			for (i = tot = 0; i < size; i++) {
				int ind = i % Block; auto &s = S[ind];
				read_scan(in, ms1h[i], s), ms1h[i].type = 1, tot += s.size();
				ms1h[i].window_low = ms1h[i].window_high = 0.0;

				if (ind == Block - 1 || i == size - 1) {
					peak_lists.push_back(std::vector<Peak>(tot));
					auto &pl = peak_lists.back();
					for (j = tot = 0; j <= ind; j++) memcpy(&(pl[tot]), &(S[j][0]), S[j].size() * sizeof(Peak)), ms1h[i - ind + j].peaks = &(pl[tot]), ms1h[i - ind + j].n = S[j].size(), tot += S[j].size();
				}
			}
		} else if (size == -2) { // second version of the .dia format
			int size_ms1 = 0, size_ms2 = 0;
			long long size_peaks = 0;
			in.read((char*)&size_ms1, sizeof(int));
			in.read((char*)&size_ms2, sizeof(int));
			in.read((char*)&size_peaks, sizeof(long long));

			ms1h.resize(size_ms1);
			ms2h.resize(size_ms2);
			peaks.resize(size_peaks);
			if (size_ms1) in.read((char*)&(ms1h[0]), size_ms1 * sizeof(Scan));
			if (size_ms2) in.read((char*)&(ms2h[0]), size_ms2 * sizeof(Scan));
			if (size_peaks) in.read((char*)&(peaks[0]), size_peaks * sizeof(Peak));

			for (int i = 0; i < ms1h.size(); i++) ms1h[i].peaks = &(peaks[(long long)ms1h[i].peaks]);
			for (int i = 0; i < ms2h.size(); i++) ms2h[i].peaks = &(peaks[(long long)ms2h[i].peaks]);
		} else { dsout << "ERROR: unknown version of the .dia file format\n"; return false; }

		in.close();
		return true;
	}

	void init() { // called after read()
		int i, cycle = 0;
		long long tot_peaks = 0;
		n_scans = ms2h.size();

		if (ExportWindows) {
			if (Verbose >= 1) Time(), dsout << "Exporting window acquisition scheme\n";
			std::ofstream out(name + std::string(".txt"), std::ofstream::out);
			out << "lower_offset\tupper_offset\n";
			for (i = 0; i < Min(n_scans, MaxCycleLength); i++) out << ms2h[i].window_low << "\t" << ms2h[i].window_high << "\n";
			out.close();
		}

		scan_cycle.resize(n_scans); if (n_scans) scan_cycle[0] = 0;
		for (i = 1; i < n_scans; i++) { // handle overlapping windows
			if (!ForceNormalSWATH && !ForceScanningSWATH) if (!scanning && ms2h[i].window_high > ms2h[i - 1].window_high + E && Abs(ms2h[i].window_low - ms2h[i - 1].window_high) < E && ms2h[i].RT - ms2h[i - 1].RT < (1.0 / (60.0 * 70.0)))
				scanning = true;
			if (ms2h[i - 1].window_low < ms2h[i].window_low - E && ms2h[i - 1].window_high < ms2h[i].window_high - E && ms2h[i - 1].window_high > ms2h[i].window_low + E)
				ms2h[i - 1].window_high = ms2h[i].window_low = (ms2h[i - 1].window_high + ms2h[i].window_low) * 0.5;
			if (ms2h[i - 1].window_low > ms2h[i].window_low + E && ms2h[i - 1].window_high > ms2h[i].window_high + E && ms2h[i - 1].window_low < ms2h[i].window_high - E)
				ms2h[i - 1].window_low = ms2h[i].window_high = (ms2h[i - 1].window_low + ms2h[i].window_high) * 0.5;
			if (ms2h[i - 1].window_high < ms2h[i].window_high && ms2h[i - 1].window_low < ms2h[i].window_low) scan_cycle[i] = cycle;
			else scan_cycle[i] = ++cycle;
		}
		for (cycle_length = 0; cycle_length < scan_cycle.size(); cycle_length++) if (scan_cycle[cycle_length]) break;

		if (ForceNormalSWATH) scanning = false;
		if (ForceScanningSWATH) scanning = true;
		if (Verbose >= 1 && scanning && !Convert) Time(), dsout << "Analysing as Scanning SWATH run\n";

		scan_RT.resize(n_scans);
		double t_min = INF, t_max = -INF;
		for (i = 0; i < n_scans; i++) {
			tot_peaks += ms2h[i].size();
			scan_RT[i] = ms2h[i].RT;
			if (ms2h[i].size()) {
				if (ms2h[i].peaks[0].mz < t_min) 
					t_min = ms2h[i].peaks[0].mz;
				if (ms2h[i].peaks[ms2h[i].size() - 1].mz > t_max) 
					t_max = ms2h[i].peaks[ms2h[i].size() - 1].mz;
			}
		}

		if (MS2Range) {
			if (t_min < t_max) tandem_min = t_min, tandem_max = t_max;
			if (Verbose >= 3) Time(), dsout << "Detected MS/MS range: " << tandem_min << " - " << tandem_max << "\n";
		}

		ms1_RT.resize(ms1h.size());
		for (i = 0; i < ms1_RT.size(); i++) ms1_RT[i] = ms1h[i].RT;

		for (int i = 0; i < ms1h.size(); i++) { // calculate ms1 windows
			if (!ms1h[i].size()) ms1h[i].window_low = ms1h[i].window_high = 0.0;
			else {
				ms1h[i].window_low = Min(ms1h[i].peaks[0].mz, ms1h[i].peaks[ms1h[i].size() - 1].mz) + MinMs1RangeOverlap;
				ms1h[i].window_high = Max(ms1h[i].peaks[0].mz, ms1h[i].peaks[ms1h[i].size() - 1].mz) - MinMs1RangeOverlap;
			}
		}
	}

#ifdef MSTOOLKIT
	void load_raw(int thread_id, char * file_name, std::vector<Scan> * spectra) {
		const int Block = 5;
		MSToolkit::MSReader r;
		std::vector<MSToolkit::Spectrum> s(Block);

		r.setFilter(MSToolkit::MS1);
		r.addFilter(MSToolkit::MS2);

		bool first = true, finish = false;
		int start = 0, stop = 0, next = 0;
		if (first) while (!lock.set()) {}
		while (true) {
			int curr = sp_alloc.fetch_add(Block);
			start = curr + 1, stop = curr + Block;
			for (next = start; next <= stop; next++) {
				if (first) {
					first = false;
					if (!r.readFile(file_name, s[next - start], next)) {
						sp_alloc.fetch_sub(Block);
						lock.free();
						return;
					}
					total_spectra = r.getLastScan();
					if (total_spectra > 1) spectra->reserve(total_spectra + 2);
					lock.free();
				} else { if (!r.readFile(NULL, s[next - start], next)) { finish = true; break; } }
				if (s[next - start].getScanNumber() == 0) { finish = true; break; }
			}

			if (next > start) {
				int tot = 0, pos = 0;
				for (int i = start; i < next; i++) tot += s[i - start].size();

				while (!lock.set()) {}
				peak_lists.push_back(std::vector<Peak>(tot));
				auto &pl = peak_lists.back();
				if (next - 1 > spectra->size()) spectra->resize(next - 1);
				for (int i = start; i < next; i++) {
					auto &sc = spectra->at(i - 1);
					auto &sp = s[i - start];
					int level = sp.getMsLevel();
					if (level != 1 && level != 2) {
						sc.type = 0;
						continue;
					}

					sc.n = sp.size();
					sc.type = level;
					sc.RT = sp.getRTime();
					sc.window_low = sp.getSelWindowLower(), sc.window_high = sp.getSelWindowUpper();
					sc.peaks = &(pl[pos]);
					for (int k = 0; k < sc.n; k++) pl[pos + k].mz = sp.at(k).mz, pl[pos + k].height = sp.at(k).intensity;
					pos += sc.n;

					if (level == 1) ms1_cnt++;
					else ms2_cnt++;
				}
				lock.free();
			}

			if (finish) break;
		}
	}
#endif

	bool load(char * file_name) { // Spectra in the file should be ordered by the acquisition time
		name = std::string(file_name);
		if (Verbose >= 1) Time(), dsout << "Loading run " << name << "\n";

		auto ext = name.substr(name.find_last_of('.'));
		if (ext == std::string(".dia")) {
			if (!read(file_name)) return false;
			goto finalise;
		}

#ifdef WIFFREADER
		typedef HANDLE(*LPLOAD)(char*, bool, int);
		if (ext == std::string(".wiff")) {
			try {
				auto load_func = (LPLOAD)GetProcAddress(wiff_dll, wiff_load_func);
				if (load_func == NULL) {
					dsout << "Cannot load the reader from DIA-NN.Wiff.dll\n";
					return false;
				}
				auto data = load_func(file_name, !fast_wiff, Max(1, Threads - 1));
				if (data == NULL) { dsout << "Unhandled exception when reading the .wiff file\n"; return false; }
				int * ptr = (int*)data;
				int i, j, pos = 0, n = *(ptr++), nms1 = 0, nms2 = 0;
				long long *lptr = (long long*)ptr, tot_peaks = *(lptr++);
				ptr = (int*)lptr;
				ms1h.reserve(n), ms2h.reserve(n), peaks.resize(tot_peaks);
				Scan S;
				for (i = 0; i < n; i++) {
					S.n = *(ptr++);
					S.type = *(ptr++);
					double *dptr = (double*)ptr;
					S.RT = *(dptr++), S.window_low = *(dptr++), S.window_high = *(dptr++);
					S.peaks = &(peaks[pos]);
					float *fptr = (float*)dptr;
					for (j = 0; j < S.n; j++) S.peaks[j].mz = *(fptr++), S.peaks[j].height = *(fptr++);
					ptr = (int*)fptr;
					pos += S.n;

					if (S.type == 1) ms1h.push_back(S);
					else ms2h.push_back(S);
				}
				GlobalFree(data);
			} catch (std::exception &e) {
				dsout << "Unhandled exception when reading the .wiff file\n";
				return false;
			};
			goto finalise;
		}
#endif

#ifdef MSTOOLKIT
		{
			std::vector<Scan> spectra;
			peaks.clear();

			sp_alloc.store(0);
			ms1_cnt = ms2_cnt = total_spectra = 0;
			std::vector<std::thread> threads;
			int max_threads = Max(1, Threads - 1);
			for (int i = 0; i < max_threads; i++) threads.push_back(std::thread(&Run::load_raw, this, i, file_name, &spectra));
			for (int i = 0; i < max_threads; i++) threads[i].join();

			if (!ms2_cnt) {
				dsout << "No MS2 spectra: aborting\n";
				return false;
			}
			if (ms1_cnt + ms2_cnt != total_spectra) dsout << "WARNING: spectra other than MS1 and MS2 detected\n";

			ms1h.resize(ms1_cnt);
			ms2h.resize(ms2_cnt);

			int i, pos = 0, n_ms1 = 0, n_ms2 = 0;
			for (auto &s : spectra) {
				if (s.type == 1) ms1h[n_ms1++] = s;
				else if (s.type == 2) ms2h[n_ms2++] = s;
			}
			std::vector<Scan>().swap(spectra);
		}
#endif

		finalise:
		if (!Convert) init();

#if (HASH > 0)
		dsout << "Raw hash: " << raw_hash() << "\n";
#endif

		if (Verbose >= 2) Time(), dsout << "Run loaded\n";
		return true;
	}

    void load_library(Library * _lib) {
		lib = _lib;
		int i, fr, frn, pos = 0, cnt = 0, gcnt = 0, guide = 0, guide_batches = 0, fasta_batches = 0;
		Target entry;
		std::vector<bool> has(lib->entries.size());
		for (auto it = lib->entries.begin(); it != lib->entries.end(); it++, pos++) {
			double mz = it->target.mz;
			int max = Min(ms2h.size(), MaxCycleLength);
			for (i = 0; i < max; i++)
				if (ms2h[i].has(mz)) {
					if (!(it->entry_flags & fFromFasta)) break;
					gcnt = 0;
					auto gen = generate_all_fragments(get_sequence(it->name), 1, AddLosses ? loss_NH3 : loss_none, &gcnt);
					for (fr = frn = 0; fr < gen.size(); fr++) {
						if (gen[fr].mz > tandem_min && gen[fr].mz < tandem_max) frn++;
						if (frn >= MinGenFrNum) break;
					}
					if (frn >= MinGenFrNum) break;
				}
			if (i < max) has[pos] = true, cnt++;
		}
		entries.resize(cnt);
		if (Verbose >= 1) Time(), dsout << cnt << " library precursors are potentially detectable\n";
		cnt = pos = 0;
		for (auto it = lib->entries.begin(); it != lib->entries.end(); it++, pos++) {
			if (has[pos]) {
				if (!(it->entry_flags & fFromFasta)) guide++;
				entries[cnt].target.init(this, &(it->target));
				entries[cnt].decoy.init(this, &(it->decoy));
				entries[cnt].decoy.flags |= fDecoy;
				cnt++;
			}
		}

		std::mt19937_64 gen(1);
		if (BatchMode) Batches = Min(MaxBatches, Max(1, entries.size() / MinBatch));
		if (Batches == 1) for (pos = 0; pos < entries.size(); pos++) entries[pos].batch = 0;
		else {
			if (GuideLibrary) guide_batches = Min(MaxBatches / 2 + 1, Max(1, guide / MinBatch)), fasta_batches = Max(1, Batches - guide_batches), Batches = guide_batches + fasta_batches;
			std::vector<int> index(entries.size());
			for (pos = 0; pos < index.size(); pos++) index[pos] = pos;
			std::shuffle(index.begin(), index.end(), gen);
			for (pos = 0; pos < entries.size(); pos++) {
				if (!GuideLibrary)
					entries[pos].batch = index[pos] % Batches;
				else {
					if (!(lib->entries[entries[pos].target.pep->index].entry_flags & fFromFasta))
						entries[pos].batch = index[pos] % guide_batches;
					else entries[pos].batch = guide_batches + (index[pos] % fasta_batches);
				}
			}
		}
	}

	Target * find_entry(int index) {
		int low = 0, high = entries.size();
		while (high > low) {
			int mid = (low + high) / 2;
			int i = entries[mid].target.pep->index;
			if (i < index) low = mid + 1;
			else if (i > index) high = mid;
			else return &(entries[mid]);
		}
		return NULL;
	}

    void process_precursors(int thread_id, bool free, bool free_fragments, int min_batch, int max_batch) { // skip before min_batch (only update RT scores), ignore after max_batch
        for (int i = 0; i < entries.size(); i++) {
			if (entries[i].batch > max_batch) continue;
			if (entries[i].lock.set()) {
				if ((!QuantOnly && curr_iter < iN) || (entries[i].target.flags & fFound)) entries[i].target.find(thread_id, free, entries[i].batch < min_batch);
				if (!QuantOnly && curr_iter < iN) entries[i].decoy.find(thread_id, free, entries[i].batch < min_batch);

				auto &le = lib->entries[entries[i].target.pep->index];
				if (free_fragments && (le.entry_flags & fFromFasta) && le.target.fragments.size() && !(le.entry_flags & fPredictedSpectrum)) le.free();
				entries[i].lock.free();
			}
        }
    }

    void seek_precursors(bool clean, bool free, bool free_fragments, int min_batch, int max_batch) {
		if (Verbose >= 2) Time(), dsout << "Precursor search\n";

        int i;
        if (Threads > 1) {
            std::vector<std::thread> threads;
            for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Run::process_precursors, this, i, free, free_fragments, min_batch, max_batch));
            for (i = 0; i < Threads; i++) threads[i].join();
        } else process_precursors(0, free, free_fragments, min_batch, max_batch);
		for (i = 0; i < entries.size(); i++) {
			entries[i].target.lock.free(), entries[i].decoy.lock.free();
			if (clean) entries[i].target.free(), entries[i].decoy.free();
		}
    }

	void quantify_precursors(int thread_id, int stage, bool filter) {
		for (int i = 0; i < entries.size(); i++) {
			if (entries[i].target.flags & fFound) if (entries[i].target.lock.set()) {
				entries[i].target.thread_id = thread_id;
				entries[i].target.quantify(entries[i].target.info[0].peak_pos, stage, filter);
			}
		}
	}

	void quantify_all(int stage, bool filter) {
		int i;
		if (Threads > 1) {
			std::vector<std::thread> threads;
			for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Run::quantify_precursors, this, i, stage, filter));
			for (i = 0; i < Threads; i++) threads[i].join();
		}
		else quantify_precursors(0, stage, filter);
		for (i = 0; i < entries.size(); i++) entries[i].target.lock.free();
	}

	void free_precursors() {
		for (int i = 0; i < entries.size(); i++)
			entries[i].target.free(), entries[i].decoy.free();
	}

	void reset_precursors() {
		for (int i = 0; i < entries.size(); i++) {
			entries[i].target.free(), entries[i].decoy.free();
			entries[i].target.flags &= ~fFound, entries[i].decoy.flags &= ~fFound;
			if (entries[i].target.info.size()) std::vector<Search>().swap(entries[i].target.info);
			if (entries[i].decoy.info.size()) std::vector<Search>().swap(entries[i].decoy.info);
		}
	}

	void remove_rubbish() {
		double smin = INF, smin_g = INF;
		int cnt = 0, tot = 0;
		for (auto &e : entries)
			if (e.target.flags & fFound)
				if (e.target.info[0].qvalue > RubbishFilter)
					e.target.flags &= ~fFound, tot++;
				else {
					if (GuideLibrary && !(lib->entries[e.target.pep->index].entry_flags & fFromFasta)) {
						if (e.target.info[0].cscore < smin_g) smin_g = e.target.info[0].cscore;
					} else if (e.target.info[0].cscore < smin) smin = e.target.info[0].cscore;
					cnt++;
				}
		
		tot += cnt;
		if (cnt >= MinNonRubbish && tot > MinRubbishFactor * cnt) { // this is likely indicative of library-free or Prosit library search
			if (Verbose >= 1) Time(), dsout << "Removing low confidence identifications\n";
			for (auto &e : entries)
				if (e.decoy.flags & fFound) {
					if (GuideLibrary && !(lib->entries[e.decoy.pep->index].entry_flags & fFromFasta)) {
						if (e.decoy.info[0].cscore < smin_g) e.decoy.flags &= ~fFound;
					} else if (e.decoy.info[0].cscore < smin) e.decoy.flags &= ~fFound;
				}
		}
	}

	void map_ids() {
		int i;

		IDs.clear();
		IDs.resize(ms2h.size());
		for (i = 0; i < entries.size(); i++)
			if (entries[i].target.flags & fFound)  
				IDs[entries[i].target.info[0].apex].push_back(i);

		for (i = 0; i < ms2h.size(); i++) {
			if (IDs[i].size() < 2) continue;
			std::sort(IDs[i].begin(), IDs[i].end(), [&](const int &left, const int &right) { return entries[left].target.info[0].cscore > entries[right].target.info[0].cscore; });
		}
	}

	inline void remove_ifs(int i, int start, Precursor &es) {
		int fr;
		float smz, nmz, error, tot_corr, ifs_corr, smax = es.info[0].cscore;

		for (int next = start + 1; next < IDs[i].size(); next++) {
			if (IDs[i][next] < 0) continue;
			auto &en = entries[IDs[i][next]].target;
			if (en.info[0].cscore > smax) continue; 
			if (lib->elution_groups[en.pep->index] == lib->elution_groups[es.pep->index]) continue; // same elution group => might actually share fragments 

			float masses[TopF];
			for (fr = 0; fr < TopF; fr++) masses[fr] = 0.0;
			for (fr = 0, tot_corr = ifs_corr = 0.0; fr < Min(TopF, en.info[0].quant.size()); fr++) tot_corr += en.info[0].quant[fr].corr;
			for (int l = 0; l < Min(TopF, en.info[0].quant.size()); l++) {
				int li = en.info[0].quant[l].index;
				if (masses[l] > E) nmz = masses[l];
				else nmz = en.pep->fragments[li].mz;
				for (int k = 0; k < Min(TopF, es.info[0].quant.size()); k++) {
					int ki = es.info[0].quant[k].index;
					smz = es.pep->fragments[ki].mz;
					error = smz * 2.0 * MassAccuracy;
					if (nmz < smz - error) continue;
					if (nmz > C13delta + smz + error) continue;
					smz = 0.0;
					if (Abs(nmz - es.pep->fragments[ki].mz) < error)
						ms2h[i].level<true>(predicted_mz(&(MassCorrection[0]), es.pep->fragments[ki].mz, es.info[0].RT), MassAccuracy, &smz);
					else {
						if (!UseIsotopes) continue;
						if (Abs(nmz - es.pep->fragments[ki].mz - C13delta / (double)es.pep->fragments[ki].charge) >= error) continue;
						ms2h[i].level<true>(predicted_mz(&(MassCorrection[0]), es.pep->fragments[ki].mz + C13delta / (double)es.pep->fragments[ki].charge, es.info[0].RT), MassAccuracy, &smz);
					}
					if (masses[l] <= E) {
						ms2h[i].level<true>(predicted_mz(&(MassCorrection[0]), en.pep->fragments[li].mz, en.info[0].RT), MassAccuracy, &nmz);
						masses[l] = nmz;
					}
					if (Abs(smz - nmz) < E && smz > E && nmz > E) {
						ifs_corr += en.info[0].quant[l].corr;
						break;
					}
				}
				if (ifs_corr > E && ifs_corr >= tot_corr - InterferenceCorrMargin) break;
			}
			if (ifs_corr > E && ifs_corr >= tot_corr - InterferenceCorrMargin) {
				IDs[i][next] = -1;
				en.free();
				en.flags &= ~fFound;
			}
		}
	}

	void refine_ids() {
		int i, j, cnt;

		std::vector<Elution> peaks;
		for (i = 0; i < ms2h.size(); i++) {
			if (IDs[i].size() < 2) continue;
			for (int start = 0; start < IDs[i].size() - 1; start++) {
				if (IDs[i][start] < 0) continue; // already removed
				auto &es = entries[IDs[i][start]].target;
				remove_ifs(i, start, es);
				
				float qmz = es.pep->mz + Q1Correction[0] + Q1Correction[1] * es.pep->mz, delta = C13delta / (double)es.pep->charge;
				for (j = Max(0, i - MaxIfsSpan); j < ms2h.size() && j <= i + MaxIfsSpan; j++) {
					if (j == i || Sgn(ms2h[j].window_low - ms2h[i].window_low) != Sgn(j - i)) continue;
					if ((scanning && Abs(j - i) <= MaxQ1Bins + 1) || ms2h[j].has(qmz, QuadrupoleError) || (UseIsotopes &&
						(ms2h[j].has(qmz + delta, QuadrupoleError) || ms2h[j].has(qmz + 2.0 * delta, QuadrupoleError)
							|| ms2h[j].has(qmz + 3.0 * delta, QuadrupoleError) || ms2h[j].has(qmz + 4.0 * delta, QuadrupoleError)))) remove_ifs(j, -1, es);
				}

				if (StrictIntRemoval) {
					for (j = i + 1, cnt = 0; j < ms2h.size(); j++) if (ms2h[j].has(es.pep->mz)) {
						remove_ifs(j, -1, es);
						if ((++cnt) >= es.S / 4) break;
					}
					for (j = i - 1, cnt = 0; j >= 0; j--) if (ms2h[j].has(es.pep->mz)) {
						remove_ifs(j, -1, es);
						if ((++cnt) >= es.S / 4) break;
					}
				}
			}
		}
	}

	void map_egs() {
		int i = 0;
		best_egs.clear(); best_egs.resize(lib->co_elution_index.size(), std::pair<int, float>(0, -INF));
		for (auto &e : entries) {
			if (e.target.flags & fFound) {
				int index = lib->elution_groups[e.target.pep->index];
				if (e.target.info[0].cscore > best_egs[index].second) best_egs[index].first = i, best_egs[index].second = e.target.info[0].cscore;
			}
			i++;
		}
	}

	void translate_peaks(int thread_id) {
		for (int i = 0; i < entries.size(); i++) {
			if (entries[i].target.lock.set()) {
				entries[i].target.thread_id = thread_id;
				entries[i].target.translate();
			}
		}
	}

	void translate_all() {
		int i;
		if (Threads > 1) {
			std::vector<std::thread> threads;
			for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Run::translate_peaks, this, i));
			for (i = 0; i < Threads; i++) threads[i].join();
		} else translate_peaks(0);
		for (i = 0; i < entries.size(); i++) entries[i].target.lock.free();
	}

	void standardise() {
		int i, pos;
		double e[pN], s2[pN], mul;

		if (Verbose >= 2) Time(), dsout << "Standardising scores\n";
		pos = 0;
		for (i = 0; i < pN; i++) e[i] = s2[i] = 0.0;
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->target.flags & fFound) { for (i = 0; i < pN; i++) e[i] += it->target.info[0].scores[i]; pos++; }
			if (it->decoy.flags & fFound) { for (i = 0; i < pN; i++) e[i] += it->decoy.info[0].scores[i]; pos++; }
		}

		if (pos < 2) return;
		mul = 1.0 / (double)pos;
		for (i = 0; i < pN; i++) e[i] *= mul;

		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->target.flags & fFound) { for (i = 0; i < pN; i++) s2[i] += Sqr(it->target.info[0].scores[i] - e[i]); }
			if (it->decoy.flags & fFound) { for (i = 0; i < pN; i++) s2[i] += Sqr(it->decoy.info[0].scores[i] - e[i]);  }
		}
		mul = 1.0 / (double)(pos - 1);
		double scale = 1.0 / StandardisationScale;
		for (i = 0; i < pN; i++) s2[i] = sqrt(Max(E, s2[i] * mul)) * scale;

		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->target.flags & fFound)
				for (i = 0; i < pN; i++) it->target.info[0].scores[i] = (it->target.info[0].scores[i] - e[i]) / s2[i];
			if (it->decoy.flags & fFound)
				for (i = 0; i < pN; i++) it->decoy.info[0].scores[i] = (it->decoy.info[0].scores[i] - e[i]) / s2[i];
		}
	}

	noninline void fit_weights(int max_batch, bool all = false) {
		int i, j, pos, t10 = 0, t50 = 0, dall = 0;
		double av[pN], mt[pN], md[pN];
		Eigen::MatrixXd A((int)pN, (int)pN);
		Eigen::VectorXd b((int)pN);
		bool non_guide_present = false, limit_ids10 = false;

#if (HASH > 1)
		dsout << "Search hash: " << search_hash() << "\n";
#endif

		if (!GuideLibrary || !GuideClassifier || curr_iter >= nnIter) all = true;

		int maxTraining = 0, tr_t = 0, tr_d = 0;
		if (curr_iter >= nnIter) {
			nn_cnt.clear();
			nn_cnt.resize(nnBagging, 0);

			for (auto it = entries.begin(); it != entries.end(); it++) {
				if (it->batch > max_batch) continue;
				if (it->target.flags & fFound) tr_t++;
				if (it->decoy.flags & fFound) tr_d++;
			}
			maxTraining = tr_t + tr_d;

			std::mt19937_64 gen(1);
			std::vector<int> index_t(tr_t), index_d(tr_d);
			for (i = 0; i < index_t.size(); i++) index_t[i] = i;
			for (i = 0; i < index_d.size(); i++) index_d[i] = i;
			std::shuffle(index_t.begin(), index_t.end(), gen);
			std::shuffle(index_d.begin(), index_d.end(), gen);

			i = 0, j = 0;
			for (auto it = entries.begin(); it != entries.end(); it++) {
				if (it->batch > max_batch) continue;
				if (it->target.flags & fFound) it->target.info[0].nn_index = index_t[i++];
				if (it->decoy.flags & fFound) it->decoy.info[0].nn_index = index_d[j++];
			}

			for (i = 0; i < nnBagging; i++) training[i].resize(maxTraining), training_classes[i].resize(maxTraining);

			for (i = j = 0; i < nnBagging; i++) {
				nn_mod[i] = i % nn_mod_n;
				if (nn_mod[i] == 0) nn_scale[i] = ((++j) % (nn_scale_max + 1 - nn_scale_min)) + nn_scale_min;
				else nn_scale[i] = nn_scale[i - 1];
				if (nn_mod[i] == 0) nn_shift[i] = gen() % nn_scale[i];
				else nn_shift[i] = nn_shift[i - 1];
			}
		}

		if (Verbose >= 2) Time(), dsout << "Optimising weights\n";

		pos = 0;
		for (i = 0; i < pN; i++) {
			av[i] = mt[i] = md[i] = 0.0;
			for (j = 0; j < pN; j++) A(i, j) = 0.0;
		}
		int nnCnt = 0, nnTest = 0;
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->batch > max_batch) continue;
			if (GuideLibrary && !all && (lib->entries[it->target.pep->index].entry_flags & fFromFasta)) {
				non_guide_present = true;
				continue;
			}

			if (curr_iter >= nnIter && (!nnFilter || !it->target.pep->no_cal)) {
				if (it->target.flags & fFound) {
					if (it->target.info[0].qvalue <= 0.5) {
						t50++;
						if (it->target.info[0].qvalue <= 0.1) t10++;
					}
					for (i = 0; i < nnBagging; i++) if (train_on_index(it->target.info[0].nn_index, i))
						training[i][nn_cnt[i]] = it->target.info[0].scores, training_classes[i][nn_cnt[i]] = target_nn, nn_cnt[i]++;
					nnCnt++;
				}
				if (it->decoy.flags & fFound) {
					dall++;
					for (i = 0; i < nnBagging; i++) if (train_on_index(it->decoy.info[0].nn_index, i))
						training[i][nn_cnt[i]] = it->decoy.info[0].scores, training_classes[i][nn_cnt[i]] = decoy_nn, nn_cnt[i]++;
					nnCnt++;
				}
			}

			if (!(it->target.flags & fFound) || !(it->decoy.flags & fFound)) continue;

			for (i = 0; i < pN; i++) {
				double x = par_learn[i] ? it->target.info[0].scores[i] : 0.0;
				double y = par_learn[i] ? it->decoy.info[0].scores[i] : 0.0;
				av[i] += x - y;
			}
			if (LDA) for (i = 0; i < pN; i++) mt[i] += it->target.info[0].scores[i], md[i] += it->decoy.info[0].scores[i];

			pos++;
		}
		if (pos < 20) {
			reset_weights();
			if (Verbose >= 3) Warning("too few training precursors, classifier will not be used");
			goto try_nn;
		}

		for (i = 0; i < pN; i++) {
			av[i] /= (double)pos;
			mt[i] = par_learn[i] ? (mt[i] / (double)pos) : 0.0;
			md[i] = par_learn[i] ? (md[i] / (double)pos) : 0.0;
		}

		pos = 0;
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->batch > max_batch) continue;
			if (!(it->target.flags & fFound) || !(it->decoy.flags & fFound)) continue;

			if (!LDA) for (i = 0; i < pN; i++) if (par_learn[i])
				for (j = i; j < pN; j++) if (par_learn[j]) A(j, i) = (A(i, j) +=
					(it->target.info[0].scores[i] - it->decoy.info[0].scores[i] - av[i]) * (it->target.info[0].scores[j] - it->decoy.info[0].scores[j] - av[j]));
			if (LDA) for (i = 0; i < pN; i++) if (par_learn[i])
					for (j = i; j < pN; j++) if (par_learn[j]) A(j, i) = (A(i, j) += 0.5 *
						((it->decoy.info[0].scores[i] - md[i]) * (it->decoy.info[0].scores[j] - md[j]) +
						(it->target.info[0].scores[i] - mt[i]) * (it->target.info[0].scores[j] - mt[j])));
			pos++;
		}

		if (pos > 1) for (i = 0; i < pN; i++) {
			for (j = i; j < pN; j++) A(j, i) = (A(i, j) /= (double)(pos - 1));
			A(i, i) += E;
			if (!par_learn[i]) for (j = 0; j < pN; j++) A(i, j) = A(j, i) = 0.0;
		}

		if (Verbose >= 3) {
			dsout << "Averages: \n";
			for (i = 0; i < pN; i++) dsout << av[i] << " ";
			dsout << "\n";
		}

		if (curr_iter < nnIter || !default_weights) {
			for (i = 0; i < pN; i++) b(i) = av[i];
			if (all) {
				solve(weights, A, b, pN);
				for (i = 0; i < pN; i++) guide_weights[i] = weights[i];
			}
			else solve(guide_weights, A, b, pN);
			default_weights = false;
			check_weights();
		}

		if (Verbose >= 3) {
			if (all) {
				dsout << "Weights: \n";
				for (i = 0; i < pN; i++) dsout << weights[i] << " ";
			} else {
				dsout << "Guide library classifier weights: \n";
				for (i = 0; i < pN; i++) dsout << guide_weights[i] << " ";
			}
			dsout << "\n";
		}

		try_nn:
		if (curr_iter != nnIter) goto finish;
		if ((t50 < 200 && t10 < 100) || dall < 100) {
			if (Verbose >= 1) Time(), dsout << "Too few confident identifications, neural network will not be used\n";
			goto finish;
		}
		if (Verbose >= 1) Time(), dsout << "Training the neural network: " << nnCnt - dall << " targets, " << dall << " decoys\n";

		{
			NNClassifier NNC(nnBagging, pN, nnHidden, nn_cnt, training, training_classes);
			NNC.run(nnEpochs);
			nn_trained = true;
			nn_scores(all, &NNC.net);

			if (TranslatePeaks) {
				calculate_qvalues();
				nn_validated = true;
				if (Verbose >= 1) Time(), dsout << "Translating peaks within elution groups\n";
				map_egs();
				translate_all();
				nn_scores(all, &NNC.net, true);
			}
		}

		finish:
		if (!all && non_guide_present) fit_weights(max_batch, true);
	}

	void nn_score(int thread_id, bool all, std::vector<NN> * net, std::vector<float> * nn_sc, bool translated) {
		int i, j, cnt, batch_size = 1000;
		std::vector<float*> data(2 * batch_size);
		std::vector<int> index(2 * batch_size);
		for (int n = 0; n * Threads + thread_id < nnBagging; n++) {
			int ind = n * Threads + thread_id;

			i = cnt = 0;
			for (auto it = entries.begin(); it != entries.end(); it++, i++) {
				if (it->batch > curr_batch) continue;
				if (GuideLibrary && !all && (lib->entries[it->target.pep->index].entry_flags & fFromFasta)) continue;
				if (it->target.flags & fFound) if (!translated || (it->target.flags & fTranslated))
					if (test_on_index(it->target.info[0].nn_index, ind))
						index[cnt] = 2 * i, data[cnt++] = it->target.info[0].scores;
				if (it->decoy.flags & fFound) if (!translated)
					if (test_on_index(it->decoy.info[0].nn_index, ind))
						index[cnt] = 2 * i + 1, data[cnt++] = it->decoy.info[0].scores;

				if (cnt >= batch_size || (cnt > 0 && it + 1 == entries.end())) {
					auto dataset = createDataSet(cnt, pN, &(data[0]));
					forwardPassDataSet((*net)[ind].network, dataset);
					Matrix* result = getOuput((*net)[ind].network);
					for (j = 0; j < cnt; j++) {
						bool decoy = index[j] & 1;
						int pos = index[j] >> 1;
						if (!decoy) (*nn_sc)[pos * (2 * nnBagging) + ind] = result->data[2 * j];
						else (*nn_sc)[pos * (2 * nnBagging) + nnBagging + ind] = result->data[2 * j];
					}
					free(dataset);
					cnt = 0;
				}
			}
		}
	}

	void nn_scores(bool all, std::vector<NN> * net, bool translated = false) {
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->batch > curr_batch) continue;
			if (it->target.flags & fFound) if (!translated || (it->target.flags & fTranslated)) it->target.info[0].cscore = 0.0, it->target.info[0].nn_inc = 0;
			if (!translated) if (it->decoy.flags & fFound) it->decoy.info[0].cscore = 0.0, it->decoy.info[0].nn_inc = 0;
		}

		std::vector<float> nn_sc(entries.size() * nnBagging * 2, -1.0);
		std::vector<std::thread> thr;
		for (int i = 0; i < Threads; i++) thr.push_back(std::thread(&Run::nn_score, this, i, all, net, &nn_sc, translated));
		for (int i = 0; i < Threads; i++) thr[i].join();

		for (int i = 0; i < entries.size(); i++) {
			auto &e = entries[i];
			if (e.batch > curr_batch) continue;
			if (e.target.flags & fFound) if (!translated || (e.target.flags & fTranslated)) {
				e.target.info[0].cscore = 0.0;
				int cnt = 0;
				for (int j = 0; j < nnBagging; j++) if (nn_sc[i * (nnBagging * 2) + j] >= 0.0) e.target.info[0].cscore += nn_sc[i * (nnBagging * 2) + j], cnt++;
				assert(cnt > 0);
				e.target.info[0].cscore /= (double)cnt;
			}
			if (e.decoy.flags & fFound) if (!translated) {
				e.decoy.info[0].cscore = 0.0;
				int cnt = 0;
				for (int j = 0; j < nnBagging; j++) if (nn_sc[i * (nnBagging * 2) + nnBagging + j] >= 0.0) e.decoy.info[0].cscore += nn_sc[i * (nnBagging * 2) + nnBagging + j], cnt++;
				assert(cnt > 0);
				e.decoy.info[0].cscore /= (double)cnt;
			}
		}
	}

	std::vector<float> target_scores, decoy_scores, guide_target_scores, guide_decoy_scores;
    int calculate_qvalues(bool use_evidence = true, bool recalc = true) {
		start:
		if (Verbose >= 2) Time(), dsout << "Calculating q-values\n";
		iRTTopPrecursors = Max(2, iRTTargetTopPrecursors);

        int i = 0, ids = 0;
		Ids01 = Ids1 = Ids10 = Ids50 = IdsCal = 0;
        iRT_cscore = iRT_ref_score = INF;
		if (curr_iter >= nnIter && nn_trained) use_evidence = false;
		if (curr_iter < nnIter - 2) use_evidence = false;

		target_scores.clear(), decoy_scores.clear(), guide_target_scores.clear(), guide_decoy_scores.clear();
        for (auto it = entries.begin(); it != entries.end(); it++, i++) {
			if (it->batch > curr_batch) continue;

			if (curr_iter < nnIter || !nn_trained) {
				bool guide = GuideLibrary && GuideClassifier && !(lib->entries[it->target.pep->index].entry_flags & fFromFasta);
				if (it->target.flags & fFound) it->target.info[0].cscore = !use_evidence ? cscore<false>(it->target.info[0].scores, guide) : it->target.info[0].scores[0];
				if (it->decoy.flags & fFound) it->decoy.info[0].cscore = !use_evidence ? cscore<false>(it->decoy.info[0].scores, guide) : it->decoy.info[0].scores[0];
			}

			if (it->target.flags & fFound) if (!(it->target.flags & fTranslated)) {
				if (!GuideLibrary || (lib->entries[it->target.pep->index].entry_flags & fFromFasta))
					target_scores.push_back(it->target.info[0].cscore);
				else
					guide_target_scores.push_back(it->target.info[0].cscore);
			}
			if (it->decoy.flags & fFound) {
				if (!GuideLibrary || (lib->entries[it->target.pep->index].entry_flags & fFromFasta))
					decoy_scores.push_back(it->decoy.info[0].cscore);
				else
					guide_decoy_scores.push_back(it->decoy.info[0].cscore);
			}
        }

		std::sort(target_scores.begin(), target_scores.end());
		std::sort(decoy_scores.begin(), decoy_scores.end());

		if (GuideLibrary) {
			std::sort(guide_target_scores.begin(), guide_target_scores.end());
			std::sort(guide_decoy_scores.begin(), guide_decoy_scores.end());
		}

		std::map<float, float> cs_qv, guide_cs_qv;

        for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->batch > curr_batch) continue;
            if (it->target.flags & fFound) {
				auto &vd = (!GuideLibrary || (lib->entries[it->target.pep->index].entry_flags & fFromFasta)) ? decoy_scores : guide_decoy_scores;
				auto &vt = (!GuideLibrary || (lib->entries[it->target.pep->index].entry_flags & fFromFasta)) ? target_scores : guide_target_scores;

				int n_targets = 0, n_better = 0, n_decoys = 0;
				double sc = it->target.info[0].cscore;
				auto pT = std::lower_bound(vt.begin(), vt.end(), sc);
				n_better = std::distance(pT, vt.end());

				auto pD = std::lower_bound(vd.begin(), vd.end(), sc);
				if (pD > vd.begin()) sc = *(--pD);
				pT = std::lower_bound(vt.begin(), vt.end(), sc);

				n_targets = std::distance(pT, vt.end());
				n_decoys = std::distance(pD, vd.end());
				it->target.info[0].qvalue = Min(1.0, ((double)Max(1, n_decoys)) / (double)Max(1, n_targets));

				if (n_better <= iRTTopPrecursors && (!GuideLibrary || !(lib->entries[it->target.pep->index].entry_flags & fFromFasta)))
					if (it->target.info[0].cscore < iRT_cscore) iRT_cscore = it->target.info[0].cscore;
				if (n_better <= iRTRefTopPrecursors && it->target.info[0].cscore < iRT_ref_score) iRT_ref_score = it->target.info[0].cscore;

				auto pair = std::pair<float, float>(it->target.info[0].cscore, it->target.info[0].qvalue);
				if (!GuideLibrary || (lib->entries[it->target.pep->index].entry_flags & fFromFasta)) {
					auto pos = cs_qv.insert(pair);
					if (pos.second) {
						if (pos.first != cs_qv.begin() && std::prev(pos.first)->second < pair.second) pos.first->second = std::prev(pos.first)->second;
						else for (auto jt = std::next(pos.first); jt != cs_qv.end() && jt->second > pair.second; jt++) jt->second = pair.second;
					}
				} else {
					auto pos = guide_cs_qv.insert(pair);
					if (pos.second) {
						if (pos.first != guide_cs_qv.begin() && std::prev(pos.first)->second < pair.second) pos.first->second = std::prev(pos.first)->second;
						else for (auto jt = std::next(pos.first); jt != guide_cs_qv.end() && jt->second > pair.second; jt++) jt->second = pair.second;
					}
				}
            }
        }

		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->batch > curr_batch) continue;
			if (curr_iter >= iN - 1) if (it->decoy.flags & fFound) {
				auto pos = (!GuideLibrary || (lib->entries[it->target.pep->index].entry_flags & fFromFasta)) ?
					cs_qv.lower_bound(it->decoy.info[0].cscore) : guide_cs_qv.lower_bound(it->decoy.info[0].cscore);
				it->decoy.info[0].qvalue = pos->second;
				if (it->decoy.info[0].qvalue > 1.0) it->decoy.info[0].qvalue = 1.0;
			}
			if (!(it->target.flags & fFound)) continue;
			auto pos = (!GuideLibrary || (lib->entries[it->target.pep->index].entry_flags & fFromFasta)) ?
				cs_qv.lower_bound(it->target.info[0].cscore) : guide_cs_qv.lower_bound(it->target.info[0].cscore);
			it->target.info[0].qvalue = pos->second;
			if (it->target.info[0].qvalue > 1.0) it->target.info[0].qvalue = 1.0;
			if (it->target.info[0].qvalue <= 0.01) ids++;
			if (it->target.info[0].qvalue <= MassCalQvalue && Abs(it->target.info[0].mass_delta) > E) IdsCal++;
			if (it->target.info[0].qvalue <= 0.5) Ids50++; else continue;
			if (it->target.info[0].qvalue <= 0.1) Ids10++; else continue;
			if (it->target.info[0].qvalue <= 0.01) Ids1++; else continue;
			if (it->target.info[0].qvalue <= 0.001) Ids01++;
		}
		RS.precursors = Ids1;
		if (Verbose >= 2 || (Verbose >= 1 && curr_iter >= iN - 1)) {
			if (Verbose < 3) Time(), dsout << "Number of IDs at 0.01 FDR: " << ids << "\n";
			if (Verbose >= 3) Time(), dsout << "Number of IDs at 50%, 10%, 1%, 0.1% FDR: "
				<< Ids50 << ", " << Ids10 << ", " << Ids1 << ", " << Ids01 << "\n";
		}
		if (!nn_validated) if (curr_iter >= nnIter && nn_trained && (Ids10 < linear_classifier_ids10 ||
			(ids < linear_classifier_ids1 && linear_classifier_ids1 >= 50 && ((double)ids/(double)linear_classifier_ids1) < ((double)linear_classifier_ids10 / (double)Max(1, Ids10)))
			|| Ids10 < 100)) {
			if (Verbose >= 1) Time(), dsout << "Too low number of IDs with NNs: reverting to the linear classifier\n";
			nn_trained = false;
			goto start;
		}

		if (Ids10 < 500 && Ids1 < 100 && !all_iters && !use_evidence) reset_weights();
		else if (use_evidence && recalc && !default_weights) {
			int Ids50o = Ids50, Ids10o = Ids10, Ids1o = Ids1, Ids01o = Ids01, idso = ids;
			ids = calculate_qvalues(false, false);
			if (Ids10 < Ids10o) {
				reset_weights();
				Ids50 = Ids50o, Ids10 = Ids10o, Ids1 = Ids1o, Ids01 = Ids01o, ids = idso;
				calculate_qvalues(true, false);
			}
		}
		if (curr_iter == nnIter - 1) linear_classifier_ids1 = Ids1, linear_classifier_ids10 = Ids10;
		return ids;
    }

	int update_classifier(bool fit, bool clean, bool free, bool free_fragments, int min_batch, int max_batch) {
		for (int i = 0; i < pN; i++) {
			par_seek[i] = (curr_iter >= pars[i].min_iter_seek);
			par_learn[i] = (curr_iter >= pars[i].min_iter_learn);
		}
		if (par_limit) {
			par_seek[pNFCorr] = par_learn[pNFCorr] = false;
			par_seek[pShadow] = par_learn[pShadow] = false;
			par_seek[pHeavy] = par_learn[pHeavy] = false;
			par_seek[pResCorr] = par_learn[pResCorr] = false;
			par_seek[pResCorrNorm] = par_learn[pResCorrNorm] = false;
		}
		seek_precursors(clean, free, free_fragments, min_batch, max_batch);
		copy_weights(selection_weights, weights);
		copy_weights(selection_guide_weights, guide_weights);
		if (LCAllScores && curr_iter == nnIter - 1) {
			for (int i = 0; i < pN; i++) par_learn[i] = true;
			par_learn[pTimeCorr] = par_learn[pShadow] = par_learn[pSig + TopF - 1] = false;
			for (int i = pShape; i < pN; i++) par_learn[i] = false; // remove pShape scores to avoid potential linear dependency
		}
		if (fit) fit_weights(max_batch);
		return calculate_qvalues();
	}

    void update(bool set_RT_window, bool set_scan_window, bool calibrate, bool calibrate_q1, bool gen_ref) {
		if (Verbose >= 2) Time(), dsout << "Calibrating retention times\n";
		int peak_width = 0, peak_cnt = 0, mass_cnt = 0, mass_cnt_ms1 = 0, min_mass_cnt = in_ref_run ? MinMassDeltaCalRef : MinMassDeltaCal;

		if (!RemoveMassAccOutliers) {
			rt_stats.clear();
			rt_ref.clear();
			rt_coo.clear();
#if Q1
			if (calibrate_q1) q1_diff.clear(), q1_diff_mz.clear();
#endif

			for (auto it = entries.begin(); it != entries.end(); it++) {
				if (!(it->target.flags & fFound)) continue;
				if (it->target.info[0].qvalue <= MassCalQvalue && Abs(it->target.info[0].mass_delta) > E) {
					mass_cnt++;
					if (calibrate) rt_coo.push_back(it->target.info[0].RT);
					if (Abs(it->target.info[0].mass_delta_ms1) > E) mass_cnt_ms1++;
				}
				if (GuideLibrary && (lib->entries[it->target.pep->index].entry_flags & fFromFasta)) continue;
				if (it->target.info[0].qvalue <= iRTMaxQvalue || it->target.info[0].cscore >= iRT_cscore) {
					int index = std::distance(entries.begin(), it);
					rt_stats.push_back(index);
					if (it->target.info[0].qvalue <= iRTMaxQvalue && it->target.info[0].cscore >= iRT_ref_score) {
						peak_width += it->target.info[0].peak_width, peak_cnt++;
#if Q1
						if (calibrate_q1) q1_diff.push_back(it->target.info[0].q1_shift), q1_diff_mz.push_back(it->target.pep->mz);
#endif
					}
					if (gen_ref && it->target.info[0].cscore >= iRT_ref_score) rt_ref.push_back(it->target.pep->index);
				}
			}
			std::sort(rt_stats.begin(), rt_stats.end(), [&](const auto& lhs, const auto& rhs) { return entries[lhs].target.pep->iRT < entries[rhs].target.pep->iRT; });
			rt_data.resize(rt_stats.size());
			for (int i = 0; i < rt_stats.size(); i++) rt_data[i] = std::pair<float, float>(entries[rt_stats[i]].target.pep->iRT, entries[rt_stats[i]].target.info[0].RT);
			if (Verbose >= 2) Time(), dsout << rt_data.size() << " precursors used for iRT estimation.\n";

			if (rt_data.size() >= 2) {
				int segments = Min(RTSegments, Max(1, 2.0 * sqrt(rt_data.size() / MinRTPredBin)));
				map_RT(RT_coeff, RT_points, rt_data, segments);
				if (curr_iter >= iN - 1) {
					flip_map(rt_data);
					std::sort(rt_data.begin(), rt_data.end());
					map_RT(iRT_coeff, iRT_points, rt_data, segments);
				}
			} else return;

			if (set_RT_window) {
				rt_delta.resize(rt_stats.size());
				for (int i = 0; i < rt_stats.size(); i++) {
					if (!GuideLibrary || (lib->entries[entries[rt_stats[i]].target.pep->index].entry_flags & fFromFasta)) rt_delta[i] = -Abs(entries[rt_stats[i]].target.info[0].RT - calc_spline(RT_coeff, RT_points, entries[rt_stats[i]].target.pep->iRT));
					else rt_delta[i] = -Abs(entries[rt_stats[i]].target.info[0].RT - calc_spline(RT_coeff, RT_points, entries[rt_stats[i]].target.pep->sRT));
				}
				std::sort(rt_delta.begin(), rt_delta.end());
				std::sort(rt_stats.begin(), rt_stats.end(), [&](const auto& lhs, const auto& rhs) { return entries[lhs].target.info[0].RT < entries[rhs].target.info[0].RT; });
				RT_min = entries[rt_stats[rt_stats.size() / 100]].target.info[0].RT;
				RT_max = entries[rt_stats[rt_stats.size() - 1 - rt_stats.size() / 100]].target.info[0].RT;
				RT_window = Max(-RTWindowMargin * rt_delta[(int)(RTWindowLoss * (double)rt_delta.size())], MinRTWinFactor * (RT_max - RT_min) / double(in_ref_run ? RTRefWindowFactor : RTWindowFactor));
				RT_windowed_search = (in_ref_run ? rt_delta.size() >= RTWinSearchMinRef : rt_delta.size() >= RTWinSearchMinCal);
				RS.RT_acc = -rt_delta[(int)(0.5 * (double)rt_delta.size())];
				if (!RT_windowed_search && Verbose >= 3) { Warning("not enough confidently identified precursors for RT-windowed search"); }
				else if (Verbose >= 1) Time(), dsout << "RT window set to " << RT_window << "\n";
			}

			if (gen_ref) {
				std::sort(rt_ref.begin(), rt_ref.end());
				lib->save(gen_ref_file, &rt_ref, false, false);
			}

			if (set_scan_window) {
				PeakWidth = (double(peak_width)) / (double)Max(1, peak_cnt);
				if (Verbose >= 1) Time(), dsout << "Peak width: " << PeakWidth << "\n";

				int window_size = Max(1, int(ScanScale * PeakWidth));
				if (Verbose >= 1) Time(), dsout << "Scan window radius set to " << window_size << "\n";
				if (!IndividualWindows) WindowRadius = window_size, InferWindow = false;

				window_calculated = true;
			}

#if Q1
			if (calibrate_q1) {
				int k;
				float diff = 0.0;

				if (q1_diff.size() >= 500 && Q1CalLinear) {
					typedef Eigen::Triplet<double> T;
					std::vector<T> TL;

					auto qd = q1_diff;
					std::sort(qd.begin(), qd.end());
					double dmin = qd[(q1_diff.size() / 5)], dmax = q1_diff.size() - (q1_diff.size() / 5) - 1;

					b.clear(), k = 0;
					for (int i = 0; i < q1_diff.size(); i++) if (q1_diff[i] >= dmin && q1_diff[i] <= dmax) {
						b.push_back(q1_diff[i]);
						TL.push_back(T(k, 0, 1)), TL.push_back(T(k++, 1, q1_diff_mz[i]));
					}

					auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), k);
					Eigen::SparseMatrix<double, Eigen::RowMajor> A(k, 2);
					A.setFromTriplets(TL.begin(), TL.end());
					Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
					SQR.compute(A);
					auto X = SQR.solve(B);
					for (int i = 0; i < 2; i++) Q1Correction[i] = X(i);
				}

				std::sort(q1_diff.begin(), q1_diff.end()), k = 0;
				for (int i = (q1_diff.size() / 5); i < q1_diff.size() - (q1_diff.size() / 5); i++) diff += q1_diff[i], k++;
				Q1Correction[0] = diff / (float)Max(1, k);

				if (Verbose >= 1) Time(), dsout << "Q1 correction: " << Q1Correction[0] << "Th + " << Q1Correction[1] << "(m/z)" << "\n";
			}
#endif
		}

		if (calibrate) {
			if (!recalibrate) {
				MassAccuracy = GlobalMassAccuracy;
				MassAccuracyMs1 = GlobalMassAccuracyMs1;
			}
			if (mass_cnt >= min_mass_cnt || RemoveMassAccOutliers) {
				if (!RemoveMassAccOutliers) b.resize(mass_cnt);
				typedef Eigen::Triplet<double> T;
				std::vector<T> TL, TL_ms1;

				if (!RemoveMassAccOutliers) {
					std::sort(rt_coo.begin(), rt_coo.end());
					MassCalBins = Max(1, Min(MassCalBinsMax, mass_cnt / MinMassCalBin));
					MassCorrection.resize(1 + 2 * MassCalBins);
					double split = ((double)(mass_cnt - 1)) / (double)MassCalBins;
					for (int i = 0; i < MassCalBins; i++) MassCalSplit[i] = rt_coo[Min(split * (double)i, mass_cnt - 1)]; MassCalSplit[MassCalBins] = rt_coo[mass_cnt - 1];
					for (int i = 0; i < MassCalBins; i++) MassCalCenter[i] = 0.5 * (MassCalSplit[i] + MassCalSplit[i + 1]);
					for (int i = 1; i < MassCalBins; i++) if (MassCalCenter[i] - MassCalCenter[i - 1] < E) {
						for (int k = i - 1; k < MassCalBins - 1; k++) MassCalCenter[k] = MassCalCenter[k + 1];
						MassCalBins--;
					}
				}

				mass_cnt = 0;
				for (auto it = entries.begin(); it != entries.end(); it++) {
					if (!(it->target.flags & fFound)) continue;
					if (it->target.info[0].qvalue <= MassCalQvalue && Abs(it->target.info[0].mass_delta) > E
						&& (!RemoveMassAccOutliers
							|| Abs(it->target.info[0].mass_delta_mz + it->target.info[0].mass_delta -
								predicted_mz(&(MassCorrection[0]), it->target.info[0].mass_delta_mz, it->target.info[0].RT))
							< MassAccOutlier * it->target.info[0].mass_delta_mz)) {
						double mz = it->target.info[0].mass_delta_mz, mzw = 1.0 / mz;
						double rt = it->target.info[0].RT;
						b[mass_cnt] = mzw * it->target.info[0].mass_delta / mz;
						TL.push_back(T(mass_cnt, 0, mzw * mz));

						if (rt <= MassCalCenter[0])
							TL.push_back(T(mass_cnt, 1, mzw * 1.0 / mz)), TL.push_back(T(mass_cnt, 2, mzw * 1.0));
						else if (rt >= MassCalCenter[MassCalBins - 1])
							TL.push_back(T(mass_cnt, 1 + (MassCalBins - 1) * 2, mzw * 1.0 / mz)), TL.push_back(T(mass_cnt, 2 + (MassCalBins - 1) * 2, mzw * 1.0));
						else for (int i = 1; i < MassCalBins; i++) if (rt < MassCalCenter[i]) {
							double u = rt - MassCalCenter[i - 1], v = MassCalCenter[i] - rt, w = u + v;
							if (w > E) {
								double cl = v / w, cr = u / w;
								TL.push_back(T(mass_cnt, 1 + (i - 1) * 2, mzw * cl / mz)), TL.push_back(T(mass_cnt, 2 + (i - 1) * 2, mzw * cl));
								TL.push_back(T(mass_cnt, 1 + i * 2, mzw * cr / mz)), TL.push_back(T(mass_cnt, 2 + i * 2, mzw * cr));
							}
							break;
						}

						mass_cnt++;
					}
				}
				if (mass_cnt < min_mass_cnt / 2) return;

				auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), mass_cnt);
				Eigen::SparseMatrix<double, Eigen::RowMajor> A(mass_cnt, MassCorrection.size());
				A.setFromTriplets(TL.begin(), TL.end());
				Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
				SQR.compute(A);
				auto X = SQR.solve(B);
				for (int i = 0; i < MassCorrection.size(); i++) MassCorrection[i] = X(i);

				if (Verbose >= 3) {
					Time(), dsout << "Mass correction transform (" << mass_cnt << " precursors): ";
					for (int i = 0; i < MassCorrection.size(); i++) dsout << MassCorrection[i] << " ";
					dsout << "\n";
				}

				b_r.clear(), b_r.reserve(2 * mass_cnt);
				for (auto it = entries.begin(); it != entries.end(); it++) {
					if (!(it->target.flags & fFound)) continue;
					if (it->target.info[0].qvalue <= MassCalQvalue && Abs(it->target.info[0].mass_delta) > E)
						b_r.push_back(Abs(it->target.info[0].mass_delta_mz + it->target.info[0].mass_delta -
							predicted_mz(&(MassCorrection[0]), it->target.info[0].mass_delta_mz, it->target.info[0].RT)) / it->target.info[0].mass_delta_mz);
				}

				if (Verbose >= 3) Time(), dsout << "M/z SD: " << sqrt(var(b_r)) * 1000000.0 << " ppm\n";
				std::sort(b_r.begin(), b_r.end());
				MassAccuracy = b_r[0.7 * (double)b_r.size()];
				RS.MS2_acc_corrected = b_r[0.5 * (double)b_r.size()];
				if (Verbose >= 2) Time(), dsout << "Top 70% mass accuracy: " << MassAccuracy * 1000000.0 << " ppm\n";
				acc_calibrated = true;

				if (SaveCalInfo) {
					std::ofstream out("cal.tsv");
					out.precision(10);
					out << "Lib\tDelta\tCorr\tRT\tQ\n";
					for (auto it = entries.begin(); it != entries.end(); it++) {
						if (!(it->target.flags & fFound)) continue;
						if (it->target.info[0].qvalue <= MassCalQvalue && Abs(it->target.info[0].mass_delta) > E)
							out << it->target.info[0].mass_delta_mz << "\t" << it->target.info[0].mass_delta << "\t"
								<< predicted_mz(&(MassCorrection[0]), it->target.info[0].mass_delta_mz, it->target.info[0].RT) << "\t"
								<< it->target.info[0].RT << "\t" << it->target.info[0].qvalue << "\n";
					}
					out.close();
				}

				// try without calibration to check if it is better
				b_r.clear();
				for (auto it = entries.begin(); it != entries.end(); it++) {
					if (!(it->target.flags & fFound)) continue;
					if (it->target.info[0].qvalue <= MassCalQvalue && Abs(it->target.info[0].mass_delta) > E)
						b_r.push_back(Abs(it->target.info[0].mass_delta) / it->target.info[0].mass_delta_mz);
				}
				std::sort(b_r.begin(), b_r.end());
				double acc_no_cal = b_r[0.7 * (double)b_r.size()];
				RS.MS2_acc = b_r[0.5 * (double)b_r.size()];
				if (Verbose >= 2) Time(), dsout << "Top 70% mass accuracy without correction: " << acc_no_cal * 1000000.0 << "ppm\n";
				if (acc_no_cal < MassAccuracy) {
					if (Verbose >= 2) Time(), dsout << "No mass correction required\n";
					MassAccuracy = acc_no_cal;
					for (int i = 0; i < MassCorrection.size(); i++) MassCorrection[i] = 0.0;
				}
			} else if (Verbose >= 1) Time(), dsout << "Cannot perform mass calibration, too few confidently identified precursors\n";

			if (mass_cnt_ms1 >= min_mass_cnt || RemoveMassAccOutliers) {
				if (!RemoveMassAccOutliers) b.resize(mass_cnt_ms1);
				typedef Eigen::Triplet<double> T;
				std::vector<T> TL, TL_ms1;

				if (!RemoveMassAccOutliers) {
					std::sort(rt_coo.begin(), rt_coo.end());
					MassCalBinsMs1 = Max(1, Min(MassCalBinsMax, mass_cnt_ms1 / MinMassCalBin));
					MassCorrectionMs1.resize(1 + 2 * MassCalBinsMs1);
					double split = ((double)(mass_cnt_ms1 - 1)) / (double)MassCalBinsMs1;
					for (int i = 0; i < MassCalBinsMs1; i++) MassCalSplitMs1[i] = rt_coo[Min(split * (double)i, mass_cnt_ms1 - 1)]; MassCalSplitMs1[MassCalBinsMs1] = rt_coo[mass_cnt_ms1 - 1];
					for (int i = 0; i < MassCalBinsMs1; i++) MassCalCenterMs1[i] = 0.5 * (MassCalSplitMs1[i] + MassCalSplitMs1[i + 1]);
					for (int i = 1; i < MassCalBinsMs1; i++) if (MassCalCenterMs1[i] - MassCalCenterMs1[i - 1] < E) {
						for (int k = i - 1; k < MassCalBinsMs1 - 1; k++) MassCalCenterMs1[k] = MassCalCenterMs1[k + 1];
						MassCalBinsMs1--;
					}
				}

				mass_cnt_ms1 = 0;
				for (auto it = entries.begin(); it != entries.end(); it++) {
					if (!(it->target.flags & fFound)) continue;
					if (it->target.info[0].qvalue <= MassCalQvalue && Abs(it->target.info[0].mass_delta) > E && Abs(it->target.info[0].mass_delta_ms1) > E
						&& (!RemoveMassAccOutliers
							|| Abs(it->target.pep->mz + it->target.info[0].mass_delta_ms1 -
								predicted_mz_ms1(&(MassCorrectionMs1[0]), it->target.pep->mz, it->target.info[0].RT))
							< MassAccMs1Outlier * it->target.pep->mz)) {
						double mz = it->target.pep->mz, mzw = 1.0 / mz;
						double rt = it->target.info[0].RT;
						b[mass_cnt_ms1] = mzw * it->target.info[0].mass_delta_ms1 / mz;
						TL.push_back(T(mass_cnt_ms1, 0, mzw * mz));

						if (rt <= MassCalCenterMs1[0])
							TL.push_back(T(mass_cnt_ms1, 1, mzw * 1.0 / mz)), TL.push_back(T(mass_cnt_ms1, 2, mzw * 1.0));
						else if (rt >= MassCalCenterMs1[MassCalBinsMs1 - 1])
							TL.push_back(T(mass_cnt_ms1, 1 + (MassCalBinsMs1 - 1) * 2, mzw * 1.0 / mz)), TL.push_back(T(mass_cnt_ms1, 2 + (MassCalBinsMs1 - 1) * 2, mzw * 1.0));
						else for (int i = 1; i < MassCalBinsMs1; i++) if (rt < MassCalCenterMs1[i]) {
							double u = rt - MassCalCenterMs1[i - 1], v = MassCalCenterMs1[i] - rt, w = u + v;
							if (w > E) {
								double cl = v / w, cr = u / w;
								TL.push_back(T(mass_cnt_ms1, 1 + (i - 1) * 2, mzw * cl / mz)), TL.push_back(T(mass_cnt_ms1, 2 + (i - 1) * 2, mzw * cl));
								TL.push_back(T(mass_cnt_ms1, 1 + i * 2, mzw * cr / mz)), TL.push_back(T(mass_cnt_ms1, 2 + i * 2, mzw * cr));
							}
							break;
						}

						mass_cnt_ms1++;
					}
				}
				if (mass_cnt_ms1 < min_mass_cnt / 2) return;

				auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), mass_cnt_ms1);
				Eigen::SparseMatrix<double, Eigen::RowMajor> A(mass_cnt_ms1, MassCorrectionMs1.size());
				A.setFromTriplets(TL.begin(), TL.end());
				Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
				SQR.compute(A);
				auto X = SQR.solve(B);
				for (int i = 0; i < MassCorrectionMs1.size(); i++) MassCorrectionMs1[i] = X(i);

				if (Verbose >= 3) {
					Time(), dsout << "MS1 mass correction transform (" << mass_cnt_ms1 << " precursors): ";
					for (int i = 0; i < MassCorrectionMs1.size(); i++) dsout << MassCorrectionMs1[i] << " ";
					dsout << "\n";
				}

				b_r.clear(), b_r.reserve(2 * mass_cnt_ms1);
				for (auto it = entries.begin(); it != entries.end(); it++) {
					if (!(it->target.flags & fFound)) continue;
					if (it->target.info[0].qvalue <= MassCalQvalue
						&& Abs(it->target.info[0].mass_delta) > E && Abs(it->target.info[0].mass_delta_ms1) > E)
						b_r.push_back(Abs(it->target.pep->mz + it->target.info[0].mass_delta_ms1 -
							predicted_mz_ms1(&(MassCorrectionMs1[0]), it->target.pep->mz, it->target.info[0].RT)) / it->target.pep->mz);
				}

				std::sort(b_r.begin(), b_r.end());
				MassAccuracyMs1 = b_r[0.7 * (double)b_r.size()];
				RS.MS1_acc_corrected = b_r[0.5 * (double)b_r.size()];
				double acc = MassAccuracyMs1;
				if (!acc_calibrated || RemoveMassAccOutliers) MassAccuracyMs1 *= 5.0;
				if (Verbose >= 2) Time(), dsout << "Top 70% MS1 mass accuracy: " << acc * 1000000.0 << " ppm\n";
				acc_ms1_calibrated = true;

				if (SaveCalInfo) {
					std::ofstream out("cal_ms1h.tsv");
					out.precision(10);
					out << "Lib\tDelta\tCorr\tRT\tQ\n";
					for (auto it = entries.begin(); it != entries.end(); it++) {
						if (!(it->target.flags & fFound)) continue;
						if (it->target.info[0].qvalue <= MassCalQvalue && Abs(it->target.info[0].mass_delta) > E && Abs(it->target.info[0].mass_delta_ms1) > E)
							out << it->target.pep->mz << "\t" << it->target.info[0].mass_delta_ms1 << "\t"
							<< predicted_mz_ms1(&(MassCorrectionMs1[0]), it->target.pep->mz, it->target.info[0].RT) << "\t"
							<< it->target.info[0].RT << "\t" << it->target.info[0].qvalue << "\n";
					}
					out.close();
				}

				b_r.clear();
				for (auto it = entries.begin(); it != entries.end(); it++) {
					if (!(it->target.flags & fFound)) continue;
					if (it->target.info[0].qvalue <= MassCalQvalue && Abs(it->target.info[0].mass_delta) > E && Abs(it->target.info[0].mass_delta_ms1) > E)
						b_r.push_back(Abs(it->target.info[0].mass_delta_ms1) / it->target.pep->mz);
				}
				std::sort(b_r.begin(), b_r.end());
				double acc_no_cal = b_r[0.7 * (double)b_r.size()];
				RS.MS1_acc = b_r[0.5 * (double)b_r.size()];
				if (Verbose >= 2) Time(), dsout << "Top 70% MS1 mass accuracy without correction: " << acc_no_cal * 1000000.0 << "ppm\n";
				if (acc_no_cal < acc) {
					if (Verbose >= 2) Time(), dsout << "No MS1 mass correction required\n";
					MassAccuracyMs1 = acc_no_cal;
					if (!acc_calibrated || RemoveMassAccOutliers) MassAccuracyMs1 *= 5.0;
					for (int i = 0; i < MassCorrectionMs1.size(); i++) MassCorrectionMs1[i] = 0.0;
				}
				if (Verbose >= 2 && MassAccuracyMs1 > Min(acc, acc_no_cal) + E) Time(), dsout << "Recommended extraction MS1 mass accuracy: "
					<< MassAccuracyMs1 * 1000000.0 << " ppm\n";
			} else if (Verbose >= 1) Time(), dsout << "Cannot perform MS1 mass calibration, too few confidently identified precursors\n";
			mz_calibrated = true;
		}
    }

	std::vector<int> calculate_protein_qvalues() {
		int i, pi, np = (PGLevel == 2 ? lib->genes.size() : (PGLevel == 1 ? lib->names.size() : lib->proteins.size()));
		std::vector<float> tsc(np, -INF), dsc(np, -INF), q(np);
		std::vector<bool> identified(np);
		std::vector<int> ids;
		if (!np) {
			if (Verbose >= 1) Time(), dsout << "No protein annotation, skipping protein q-value calculation\n";
			return ids;
		}

		if (Verbose >= 1) Time(), dsout << "Calculating protein q-values\n";

		for (auto &e : entries) {
			auto &le = lib->entries[e.target.pep->index];
			if (!le.proteotypic) continue;
			auto &pid = lib->protein_ids[le.pid_index];
			if (e.target.flags & fFound) {
				float sc = e.target.info[0].cscore;
				if (PGLevel == 0) for (auto &p : pid.proteins) if (lib->proteins[p].id.size()) {
					if (sc > tsc[p]) tsc[p] = sc;
					if (e.target.info[0].qvalue <= 0.01) identified[p] = true;
				}
				if (PGLevel) for (auto &p : pid.proteins) {
					if (PGLevel == 1) {
						pi = lib->proteins[p].name_index;
						if (!lib->names[pi].size()) continue;
					} else {
						pi = lib->proteins[p].gene_index;
						if (!lib->genes[pi].size()) continue;
					}
					if (sc > tsc[pi]) tsc[pi] = sc;
					if (e.target.info[0].qvalue <= 0.01) identified[pi] = true;
				}
			}
			if (e.decoy.flags & fFound) {
				float sc = e.decoy.info[0].cscore;
				if (PGLevel == 0) for (auto &p : pid.proteins) if (lib->proteins[p].id.size()) if (sc > dsc[p]) dsc[p] = sc;
				if (PGLevel) for (auto &p : pid.proteins) {
					if (PGLevel == 1) {
						pi = lib->proteins[p].name_index;
						if (!lib->names[pi].size()) continue;
					} else if (PGLevel == 2) {
						pi = lib->proteins[p].gene_index;
						if (!lib->genes[pi].size()) continue;
					}
					if (sc > dsc[pi]) dsc[pi] = sc;
				}
			}
		}

		auto targets = tsc;
		auto decoys = dsc;
		std::sort(targets.begin(), targets.end());
		std::sort(decoys.begin(), decoys.end());

		std::map<float, float> cs_qv;
		for (i = 0; i < np; i++) {
			if (tsc[i] < -INF + 1.0) {
				q[i] = 1.0;
				continue;
			}
			float sc = tsc[i];
			auto pD = std::lower_bound(decoys.begin(), decoys.end(), sc);
			if (pD > decoys.begin()) sc = *(--pD);
			auto pT = std::lower_bound(targets.begin(), targets.end(), sc);

			q[i] = Min(1.0, ((double)std::distance(pD, decoys.end())) / Max(1.0, (double)std::distance(pT, targets.end())));

			auto pair = std::pair<float, float>(tsc[i], q[i]);
			auto pos = cs_qv.insert(pair);
			if (pos.second) {
				if (pos.first != cs_qv.begin() && std::prev(pos.first)->second < pair.second) pos.first->second = std::prev(pos.first)->second;
				else for (auto jt = std::next(pos.first); jt != cs_qv.end() && jt->second > pair.second; jt++) jt->second = pair.second;
			}
		}

		int ids01 = 0, ids01p = 0;
		for (i = 0; i < np; i++) {
			if (q[i] < 1.0) {
				auto pos = cs_qv.lower_bound(tsc[i]);
				q[i] = pos->second;
				if (q[i] <= 0.01) ids01p++;
			}
			if (identified[i]) ids01++;
		}
		if (Verbose >= 1) Time(), dsout << "Number of " << (PGLevel == 0 ? "protein isoforms" : (PGLevel == 1 ? "proteins" : "genes")) << " identified at 1% FDR: "
			<< ids01 << " (precursor-level), " << ids01p << " (protein-level) (inference performed using proteotypic peptides only)\n";
		RS.proteins = ids01p;

		for (auto &e : entries) {
			auto &le = lib->entries[e.target.pep->index];
			auto &pid = lib->protein_ids[le.pid_index];
			if (e.target.flags & fFound) {
				e.target.info[0].protein_qvalue = 1.0;
				if (PGLevel == 0) for (auto &p : pid.proteins) if (q[p] < e.target.info[0].protein_qvalue) e.target.info[0].protein_qvalue = q[p];
				if (PGLevel) for (auto &p : pid.proteins) {
					pi = (PGLevel == 1 ? lib->proteins[p].name_index : lib->proteins[p].gene_index);
					if (q[pi] < e.target.info[0].protein_qvalue) e.target.info[0].protein_qvalue = q[pi];
				}
			}
		}

		for (i = 0; i < lib->proteins.size(); i++) {
			auto &p = lib->proteins[i];
			if (PGLevel == 0 && q[i] <= ProteinIDQvalue) ids.push_back(i);
			else if (PGLevel == 1 && q[p.name_index] <= ProteinIDQvalue) ids.push_back(i);
			else if (PGLevel == 2 && q[p.gene_index] <= ProteinIDQvalue) ids.push_back(i);
		}
		std::sort(ids.begin(), ids.end());
		return ids; // identified protein isoforms at <= ProteinIDQvalue
	}

	bool reference_run(Library * ref) {
		if (!RTWindowedSearch && (!RefCal || (!InferWindow && !Calibrate))) return false;
		if (Verbose >= 1) Time(), dsout << "Calibrating retention times using a set of reference precursors.\n";

		in_ref_run = true, curr_batch = Batches;
		load_library(ref);
		for (curr_iter = (entries.size() >= 2 * pN ? 0 : RTWindowIter); curr_iter <= RTWindowIter; curr_iter++) {
			for (int i = 0; i < pN; i++) {
				par_seek[i] = (curr_iter >= pars[i].min_iter_seek);
				par_learn[i] = (curr_iter >= pars[i].min_iter_learn);
			}
			seek_precursors(false, false, false, 0, Batches);
			fit_weights(Batches);
			calculate_qvalues();
		}
		update(RTWindowedSearch, RefCal && InferWindow, RefCal && Calibrate, Q1Cal && scanning, false);
		free_precursors();
		in_ref_run = false;

		return true;
	}

	void process(std::vector<QuantEntry> * result = NULL, float q_filtering = 0.01) {
		int i, ids = 0, best_ids = 0, best1, best50, cal_batch = 0;

		if (QuantOnly) goto report;

		{
			recalibrate = (Calibrate && !mz_calibrated);
			bool do_reset = (InferWindow && !window_calculated) || recalibrate;
			bool do_calibrate = do_reset || RTWindowedSearch;
			int start_batch = 0;
			full_spectrum = true;
			if (Verbose == 1) Time(), dsout << "Processing...\n";
			if (do_calibrate) {
			calibrate:
				full_spectrum = false;
				for (curr_batch = 0; curr_batch < Batches; curr_batch++) {
					if (curr_batch >= 5 && curr_batch < Batches - 1) {
						int first_batch = curr_batch;
						curr_batch = Min(Batches - 1, curr_batch + curr_batch / 5);
						if (Verbose >= 2) Time(), dsout << "Processing batches #" << first_batch + 1 << "-" << curr_batch + 1 << " out of " << Batches << " \n";
					} else if (Verbose >= 2) {
						if (Batches == 1) Time(), dsout << "Processing\n";
						else Time(), dsout << "Processing batch #" << curr_batch + 1 << " out of " << Batches << " \n";
					}

					curr_iter = curr_batch ? 1 : 0;
					update_classifier(true, false, false, false, 0, curr_batch); ids = Ids10;
					update(false, false, false, false, false);
					cal_batch = curr_batch;
					if (ids >= MinCal && do_calibrate) {
						recalibrate = false;
						curr_iter = CalibrationIter;
						update_classifier(true, false, false, false, 0, curr_batch);
						curr_iter = 2;
						if (IdsCal >= MinCal / 2) break;
					}
					if (recalibrate && ids >= MinCalRec && do_calibrate) break;
				}
				if (recalibrate) {
					curr_iter = CalibrationIter;
					update_classifier(true, false, false, false, 0, curr_batch);
					RemoveMassAccOutliers = false, recalibrate = true;
					update(false, false, true, false, false), recalibrate = false;
					MassAccuracy = Min(CalibrationMassAccuracy, MassAccuracy * 5.0);
					MassAccuracyMs1 = Min(CalibrationMassAccuracy, MassAccuracyMs1 * 5.0);
					if (Verbose >= 2) Time(), dsout << "Recalibrating with mass accuracy " << MassAccuracy << ", " << MassAccuracyMs1 << " (MS2, MS1)\n";
					mz_calibrated = acc_calibrated = false;
					reset_precursors();
					reset_weights();
					goto calibrate;
				}
				for (curr_iter = curr_iter + 1; curr_iter <= CalibrationIter; curr_iter++)
					update_classifier(true, false, false, false, 0, curr_batch);
				RemoveMassAccOutliers = false;
				update(RTWindowedSearch, InferWindow && !window_calculated, Calibrate && !mz_calibrated, Q1Cal && scanning, false);
				if (acc_calibrated) {
					if (Verbose >= 3) Time(), dsout << "Refining mass correction\n";
					RemoveMassAccOutliers = true;
					MassAccOutlier = MassAccuracy;
					MassAccMs1Outlier = MassAccuracyMs1;
					update(false, false, true, false, false);
					RemoveMassAccOutliers = false;
				}
				if (ForceMassAcc) {
					MassAccuracy = GlobalMassAccuracy, MassAccuracyMs1 = GlobalMassAccuracyMs1;
					if (Verbose >= 2) Time(), dsout << "Using mass accuracy " << MassAccuracy << ", " << MassAccuracyMs1 << " (MS2, MS1)\n";
				}

				if (do_reset) reset_precursors();
				else start_batch = curr_batch;
			}

			full_spectrum = !(acc_calibrated && !ForceMassAcc);
			for (curr_batch = start_batch; curr_batch < Batches; curr_batch++) {
				if (curr_batch >= start_batch + 5 && curr_batch < Batches - 1) {
					int first_batch = curr_batch;
					curr_batch = Min(Batches - 1, curr_batch + (curr_batch - start_batch) / 5);
					if (Verbose >= 2) Time(), dsout << "Processing batches #" << first_batch + 1 << "-" << curr_batch + 1 << " out of " << Batches << " \n";
				} else if (Verbose >= 2) {
					if (Batches == 1) Time(), dsout << "Processing\n";
					else Time(), dsout << "Processing batch #" << curr_batch + 1 << " out of " << Batches << " \n";
				}

				curr_iter = CalibrationIter + 1;
				update_classifier(true, false, false, false, 0, curr_batch); ids = Ids10;
				copy_weights(best_weights, selection_weights), copy_weights(best_guide_weights, selection_guide_weights), best_ids = ids;
				update(false, false, false, false, false);
				if (ids >= MinClassifier) break;
				if (acc_calibrated && curr_batch >= cal_batch) break;
			}
			curr_iter = CalibrationIter + 2;
			if (acc_calibrated && !ForceMassAcc) {
				double best_acc = MassAccuracy, start_acc = MassAccuracy;
				reset_weights();
				update_classifier(true, false, false, false, 0, curr_batch), best_ids = ids = Ids10, best1 = Ids1, best50 = Ids50;
				int fail = 0;
				while (true) {
					MassAccuracy *= 1.2;
					if (MassAccuracy > start_acc * 10.0) break;
					if (Verbose >= 3) Time(), dsout << "Trying mass accuracy " << MassAccuracy * 1000000.0 << " ppm\n";
					reset_precursors();
					reset_weights();
					update_classifier(true, false, true, false, 0, curr_batch); ids = Ids10;
					if (ids > best_ids) best_acc = MassAccuracy, best_ids = ids, best1 = Ids1, best50 = Ids50, fail = 0;
					else fail++;
					if (fail >= 3) break;
				}
				MassAccuracy = best_acc, fail = 0;
				if (Verbose >= 1) Time(), dsout << "Optimised mass accuracy: " << MassAccuracy * 1000000.0 << " ppm\n";
				reset_precursors();

				reset_weights();
				full_spectrum = true;
				update_classifier(true, false, false, false, 0, curr_batch), ids = Ids10;
				if (ids < best_ids) {
					reset_weights(), par_limit = true;
					if (Verbose >= 3) Time(), dsout << "Resetting weights and preventing the linear classifier from using the full set of weights in the future\n";
				} else update(false, false, false, false, false);

				GlobalMassAccuracy = MassAccuracy;
				GlobalMassAccuracyMs1 = MassAccuracyMs1;
			}
			if (!IndividualMassAcc) ForceMassAcc = true;

			int processed_batches = curr_batch, fail = 0, switched = 0;
			copy_weights(best_weights, selection_weights), copy_weights(best_guide_weights, selection_guide_weights), best_ids = Ids10;
			for (; curr_iter < nnIter - 1; curr_iter++) {
				update_classifier(true, false, false, false, 0, curr_batch);
				if (Ids10 < best_ids && switched <= 2) {
					if (Verbose >= 3) Time(), dsout << "Trying the other linear classifier\n";
					LDA ^= 1, switched++;
					int last10 = Ids10;
					copy_weights(weights, best_weights), copy_weights(guide_weights, best_guide_weights);
					update_classifier(true, false, false, false, 0, curr_batch);
					if (Ids10 <= last10) {
						if (Verbose >= 3) Time(), dsout << "Switching back\n";
						LDA ^= 1;
					}
				}
				if (Ids10 >= best_ids) {
					copy_weights(best_weights, selection_weights), copy_weights(best_guide_weights, selection_guide_weights), best_ids = Ids10;
					update(false, false, false, false, false);
				} else {
					copy_weights(weights, best_weights), copy_weights(guide_weights, best_guide_weights);
					if (Verbose >= 3) Time(), dsout << "Reverting weights\n";
					if (curr_iter == CalibrationIter + 2 && !par_limit) {
						par_limit = true;
						if (Verbose >= 3) Time(), dsout << "Preventing the linear classifier from using the full set of weights in the future\n";
					} else fail++;
					update_classifier(true, false, false, false, 0, curr_batch);
					if (Ids10 < best_ids) {
						if (Verbose >= 3) Time(), dsout << "Reverting weights\n";
						copy_weights(weights, best_weights), copy_weights(guide_weights, best_guide_weights);
					} else best_ids = Ids10, update(false, false, false, false, false);
					if (fail >= 1) {
						curr_iter = nnIter - 1;
						break;
					}
				}
			}
			if (Batches > 1) {
				assert(nnIter >= iN - 1);
				free_precursors();
			}
			curr_batch = Batches;
			all_iters = true;
			copy_weights(weights, best_weights), copy_weights(guide_weights, best_guide_weights);
			update_classifier(true, false, true, true, processed_batches + 1, Batches);
			update(false, false, false, false, false);

			remove_rubbish();
			if (IDsInterference) {
				if (Verbose >= 1) Time(), dsout << "Removing interfering precursors\n";
				quantify_all(1, false);
				for (int ir = 0; ir < IDsInterference; ir++) {
					map_ids();
					refine_ids();
					calculate_qvalues();
				}
			}

			curr_iter = nnIter;
			if (curr_iter >= iN - 1) {
				free_precursors();
				if (curr_iter >= iN) {
					if (TranslatePeaks) {
						if (Verbose >= 1) Time(), dsout << "Translating peaks within elution groups\n";
						map_egs();
						translate_all();
						calculate_qvalues();
					}
					goto report;
				}
			}

			if (Standardise) standardise();
			fit_weights(Batches); // train the neural network
			calculate_qvalues();
			update(false, false, false, false, GenRef && curr_iter == iN - 1);
		}

	report:
		curr_iter = iN;
		Quant quant;
		quant.proteins = calculate_protein_qvalues();
		if (Verbose >= 1) Time(), dsout << "Quantification\n";
		quantify_all(2, true);
		if (TranslatePeaks) quantify_all(3, true);
		if (no_report) return;

		if (args.find("--!list-peptides") != std::string::npos) {
			std::map<std::string, float> pep;
			std::pair<std::string, float> p;
			for (auto &e : entries) {
				float abundance = 0.0;
				if (e.target.flags & fFound) if (e.target.info[0].qvalue <= 0.01) abundance = e.target.info[0].quantity;
				p.first = get_aas(lib->entries[e.target.pep->index].name);
				p.second = abundance;
				auto ins = pep.insert(p);
				if (!ins.second) if (ins.first->second < abundance) ins.first->second = abundance;
			}
			std::ofstream out("peptides_list.tsv");
			out << "Peptide\tQuantity\n";
			for (auto &it : pep) out << it.first << "\t" << it.second << "\n";
			out.close();
		}

		for (i = 0; i < ms1h.size(); i++) {
			float sum = 0; for (auto j = 0; j < ms1h[i].n; j++) sum += ms1h[i].peaks[j].height;
			RS.MS1_tot += sum;
		}
		if (ms2h.size()) {
			float RT_min = ms2h[0].RT, RT_max = ms2h[ms2h.size() - 1].RT;
			for (i = 0; i < ms2h.size(); i++) {
				int bin = Max(0, Min(TIC_n - 1, floor(((ms2h[i].RT - RT_min) / Max(E, RT_max - RT_min)) * (float)TIC_n)));
				float sum = 0; for (auto j = 0; j < ms2h[i].n; j++) sum += ms2h[i].peaks[j].height;
				RS.TIC[bin] += sum;
				RS.MS2_tot += sum;
			}
		}

		QuantEntry qe;
		DecoyEntry de;
		for (i = 0; i < pN; i++) quant.weights[i] = weights[i], quant.guide_weights[i] = guide_weights[i];
		quant.run_index = run_index;
		quant.lib_size = lib->entries.size();
		quant.RS = RS;
		quant.tandem_min = tandem_min, quant.tandem_max = tandem_max;
		quant.Q1Correction[0] = Q1Correction[0], quant.Q1Correction[1] = Q1Correction[1];
		quant.MassAccuracy = MassAccuracy, quant.MassAccuracyMs1 = MassAccuracyMs1;
		quant.MassCorrection = MassCorrection, quant.MassCorrectionMs1 = MassCorrectionMs1;
		quant.MassCalSplit = MassCalSplit, quant.MassCalSplitMs1 = MassCalSplitMs1, quant.MassCalCenter = MassCalCenter, quant.MassCalCenterMs1 = MassCalCenterMs1;

		double q_threshold = (result == NULL ? QFilter : q_filtering);
		if (result != NULL) {
			int tot_entries = 0;
			for (auto it = entries.begin(); it != entries.end(); it++) if (it->target.flags & fFound) if (it->target.info[0].qvalue <= q_threshold) tot_entries++;
			result->clear(); result->reserve(tot_entries);
		}
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (result == NULL) if (it->decoy.flags & fFound) if (it->decoy.info[0].qvalue <= QFilter) {
				de.index = it->target.pep->index;
				de.qvalue = it->decoy.info[0].qvalue;
				quant.decoys.push_back(de);
			}

			if (!(it->target.flags & fFound)) continue;
			if (it->target.info[0].qvalue > q_threshold) continue;

			qe.index = std::distance(entries.begin(), it);
			qe.window = it->target.S;
			qe.pr.decoy_found = (it->decoy.flags & fFound) != 0;
			qe.pr.apex = it->target.info[0].apex;
			qe.pr.peak = it->target.info[0].peak_pos;
			qe.pr.best_fragment = it->target.info[0].best_fragment;
			qe.pr.best_fr_mz = it->target.info[0].best_fr_mz;
			qe.pr.peak_width = it->target.info[0].peak_width;
			qe.pr.run_index = run_index;
			qe.pr.index = it->target.pep->index;
			qe.pr.iRT = it->target.pep->iRT;
			qe.pr.predicted_RT = calc_spline(RT_coeff, RT_points, qe.pr.iRT);
			qe.pr.quantity = it->target.info[0].quantity;
			qe.pr.qvalue = it->target.info[0].qvalue;
			qe.pr.protein_qvalue = it->target.info[0].protein_qvalue;
			qe.pr.RT = scan_RT[it->target.info[0].apex];
			qe.pr.predicted_iRT = calc_spline(iRT_coeff, iRT_points, qe.pr.RT);
			qe.pr.evidence = it->target.info[0].evidence;
			qe.pr.ms1_corr = it->target.info[0].ms1_corr;
			qe.pr.cscore = it->target.info[0].cscore;
			qe.pr.decoy_evidence = (it->decoy.flags & fFound) ? it->decoy.info[0].evidence : 0.0;
			qe.pr.decoy_cscore = (it->decoy.flags & fFound) ? it->decoy.info[0].cscore : -INF;
			if (qe.pr.decoy_found) qe.pr.decoy_qvalue = it->decoy.info[0].qvalue;
			else qe.pr.decoy_qvalue = 1.0;
			qe.pr.profile_qvalue = 1.0;

			qe.fr_n = it->target.info[0].quant.size();
			for (i = 0; i < qe.fr_n; i++) qe.fr[i] = it->target.info[0].quant[i];
			if (qe.fr_n < TopF) memset(&(qe.fr[qe.fr_n]), 0, (TopF - qe.fr_n) * sizeof(Fragment));
#if REPORT_SCORES
			for (i = 0; i < pN; i++) qe.pr.scores[i] = it->target.info[0].scores[i];
			if (qe.pr.decoy_found) { for (i = 0; i < pN; i++) qe.pr.decoy_scores[i] = it->decoy.info[0].scores[i]; } else for (i = 0; i < pN; i++) qe.pr.decoy_scores[i] = 0.0;
#endif
#if ELUTION_PROFILE
			assert(it->target.ms1_elution_profile.size() == 2 * ElutionProfileRadius + 1);
			for (i = 0; i <= 2 * ElutionProfileRadius; i++) qe.ms1_elution_profile[i] = it->target.ms1_elution_profile[i];
#endif
			if (result == NULL) quant.entries.push_back(qe);
			else result->push_back(qe);
		}
		if (result != NULL) return;

		if (QuantInMem) { quants.push_back(quant); dsout << "\n"; }
		else {
			std::ofstream out;
			std::string out_quant;
			if (temp_folder.size()) out = std::ofstream(out_quant = location_to_file_name(name) + std::string(".quant"), std::ofstream::binary);
			else {
				out = std::ofstream(out_quant = name + std::string(".quant"), std::ofstream::binary);
				if (out.fail()) {
					dsout << "WARNING: cannot save the .quant file to the raw data folder; using the current working folder instead\n";
					out = std::ofstream(out_quant = location_to_file_name(name) + std::string(".quant"), std::ofstream::binary);
				}
			}
			if (out.fail()) dsout << "ERROR: cannot save the .quant file to " << out_quant << ". Check if the destination folder is write-protected or the file is in use\n";
			else {
				quant.write(out);
				long long fsize = out.tellp();
				out.close();
				if (Verbose >= 1) Time(), dsout << "Quantification information saved to " << out_quant << ".\n\n";
#ifdef _MSC_VER
				if (args.find("--!mmap-quant") != std::string::npos) {
					std::string mapping_name(1000, 0);
					sprintf(&(mapping_name[0]), "DIA-NN_%d_%d", GetCurrentProcessId(),run_index);
					for (i = 0; i < 1000; i++) if (mapping_name[i] == 0) break; mapping_name.resize(i);
					if (Verbose >= 1) Time(), dsout << "Creating a file mapping " << mapping_name << " with the quantification information\n";
					auto map = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, fsize, mapping_name.c_str());
					if (map == NULL) dsout << "ERROR: Cannot create file mapping\n";
					else {
						auto buf = (LPTSTR)MapViewOfFile(map, FILE_MAP_ALL_ACCESS, 0, 0, fsize);
						if (buf == NULL) dsout << "ERROR: Cannot access file mapping\n";
						else {
							std::stringstream ss(buf);
							quant.write(ss);
							quant.read(ss);
						}
					}
				}
#endif
			}
		}

		if (IndividualReports) {
			lib->info.load(lib, quant);
			lib->info.quantify();
			lib->quantify_proteins(TopN, ProteinQuantQvalue);
			lib->report(ms_files[run_index] + std::string(".tsv"));
			if (lib->genes.size()) lib->gene_report(ms_files[run_index] + std::string(".genes.tsv"));
		}
    }

	noninline void update_library() { // update fragment intensities in the spectral library
		int cnt = 0;
		for (auto &e : entries) {
			auto &le = lib->entries[e.target.pep->index];
			if (le.best_run != run_index || le.qvalue > ReportQValue || le.protein_qvalue > ReportProteinQValue) continue;
			e.target.thread_id = 0;
			auto gen = le.fragment();
			e.target.intensities(le.apex, gen);
			if (le.best_run >= 0) cnt++;
		}
		if (Verbose >= 1) Time(), dsout << cnt << " precursors added to the library\n";
	}

#if (HASH > 0)
	unsigned int raw_hash() {
		unsigned int res = 0;
		for (int i = 0; i < ms2h.size(); i++) res ^= ms2h[i].hash();
		for (int i = 0; i < ms1h.size(); i++) res ^= ms1h[i].hash();
		return res;
	}

	unsigned int search_hash() {
		unsigned int res = 0;
		for (int i = 0; i < entries.size(); i++) res ^= entries[i].hash();
		return res;
	}
#endif
};

#if (EXTERNAL > 0)
std::vector<QuantEntry> analyse(char * lib_bytes, long long lib_bytes_size, int threads, double acc_ms1, double acc_ms2, int win_size, float q_value_threshold,
	int n_ms1, float *ms1_rt, int *ms1_peak_n, float *ms1_peaks,
	int n_ms2, float *ms2_rt, float *ms2_low, float *ms2_high, int *ms2_peak_n, float *ms2_peaks) {

	if (acc_ms2 > E) {
		ForceMassAcc = true;
		GlobalMassAccuracy = acc_ms2;
		if (acc_ms1 > E) GlobalMassAccuracyMs1 = acc_ms1;
	}
	if (win_size) WindowRadius = win_size, InferWindow = false;

	auto str = std::string(lib_bytes, lib_bytes_size);
	std::stringstream ls(str);
	Library lib;
	lib.read(ls);
	lib.initialise(true);

	Threads = threads;
	Run run(0);

	run.ms1h.resize(n_ms1);
	run.ms2h.resize(n_ms2);

	long long pos;
	for (int i = pos = 0; i < n_ms1; i++) run.ms1h[i].peaks = pos + ((Peak*)ms1_peaks), pos += ms1_peak_n[i], run.ms1h[i].n = ms1_peak_n[i], run.ms1h[i].RT = ms1_rt[i], 
		run.ms1h[i].window_low = run.ms1h[i].window_high = 0.0;
	for (int i = pos = 0; i < n_ms2; i++) run.ms2h[i].peaks = pos + ((Peak*)ms2_peaks), pos += ms2_peak_n[i], run.ms2h[i].n = ms2_peak_n[i], run.ms2h[i].RT = ms2_rt[i],
		run.ms2h[i].window_low = ms2_low[i], run.ms2h[i].window_high = ms2_high[i];

	run.init();
	run.load_library(&lib);

	std::vector<QuantEntry> result;
	run.process(&result, q_value_threshold);
	return result;
}

int analyse_unit_test(int threads) {
	std::vector<char> lib_bytes;
	std::vector<float> ms1_rt, ms2_rt, ms2_low, ms2_high;
	std::vector<int> ms1_peak_n, ms2_peak_n;

	std::ifstream lib("yeast_ref.tsv.speclib", std::ifstream::binary);
	if (lib.fail()) { dsout << "Cannot open the spectral library file\n"; return -1; }
	lib.seekg(0, lib.end);
	long long size = lib.tellg();
	lib.seekg(0, lib.beg);
	lib_bytes.resize(size);
	lib.read(&(lib_bytes[0]), size);
	lib.close();

	Run run(0);
	run.load("D:/Raw/USP/20181113_normalswath_4.wiff.dia");
	run.init();

	ms1_rt.resize(run.ms1h.size()), ms1_peak_n.resize(run.ms1h.size());
	ms2_rt.resize(run.ms2h.size()), ms2_low.resize(run.ms2h.size()), ms2_high.resize(run.ms2h.size()), ms2_peak_n.resize(run.ms2h.size());

	for (int i = 0; i < run.ms1h.size(); i++) ms1_rt[i] = run.ms1h[i].RT, ms1_peak_n[i] = run.ms1h[i].n;
	for (int i = 0; i < run.ms2h.size(); i++) ms2_rt[i] = run.ms2h[i].RT, ms2_peak_n[i] = run.ms2h[i].n,
		ms2_low[i] = run.ms2h[i].window_low, ms2_high[i] = run.ms2h[i].window_high;

	auto result = analyse(&(lib_bytes[0]), lib_bytes.size(), threads, 20.0 / 1000000.0, 20.0 / 1000000.0, 8, 0.01, run.ms1h.size(), &(ms1_rt[0]), &(ms1_peak_n[0]), (float*)(run.ms1h[0].peaks),
		run.ms2h.size(), &(ms2_rt[0]), &(ms2_low[0]), &(ms2_high[0]), &(ms2_peak_n[0]), (float*)(run.ms2h[0].peaks));

	//dsout << result.size() << " precursors at 1% FDR:\n";
	//dsout << "Precursor.Index.In.Library\tPrecursor.RT\tPrecursor.iRT\tPrecursor.Quantity\tPrecursor.QValue\n";
	//for (auto &qe : result) dsout << qe.pr.index << '\t' << qe.pr.RT << '\t' << qe.pr.iRT << '\t' << qe.pr.quantity << '\t' << qe.pr.qvalue << '\n';
	return result.size();
}
#endif

int main(int argc, char *argv[]) {
#ifdef _MSC_VER
	DWORD mode;
	HANDLE console = GetStdHandle(STD_INPUT_HANDLE);
	GetConsoleMode(console, &mode);
	SetConsoleMode(console, (mode & ~(ENABLE_QUICK_EDIT_MODE | ENABLE_MOUSE_INPUT | ENABLE_WINDOW_INPUT)) | ENABLE_PROCESSED_INPUT);
	_set_FMA3_enable(0); // essential to make calculations reproducible between CPUs with (new ones) and without (old ones) FMA support
#endif
	std::cout.setf(std::ios::unitbuf);
	auto curr_time = time(0);
	dsout << std::ctime(&curr_time) << "DIA-NN (Data Independent Acquisition by Neural Networks)\nCompiled on " << __DATE__ << " " << __TIME__ << "\n";
#ifdef _MSC_VER
	cpu_info(dsout);
#endif
	dsout << "Logical CPU cores: " << std::thread::hardware_concurrency() << "\n";
#if (HASH > 0)
	init_hash();
#endif
	init_aas();
	arguments(argc, argv);
	Eigen::setNbThreads(Threads);

#if (EXTERNAL > 0)
	if (args.find("--!test") != std::string::npos) {
		dsout << "Analyse() function unit test:\n"; analyse_unit_test(Threads);
		dsout << "\n";
	}
#endif

	if (Verbose >= 1) dsout << ms_files.size() << " files will be processed\n";
	if (Convert) {
		for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
			Run run(std::distance(ms_files.begin(), it));
			if (!run.load(&((*it)[0]))) continue;
			std::string(dia_name) = run.name + std::string(".dia");
			if (out_dir.size()) dia_name = out_dir + '\\' + get_file_name(dia_name);
			run.write(&(dia_name[0]));
		}
		dsout << "Finished\n\n";
		return 0;
	}

	Fasta fasta;
	if (FastaSearch && fasta_files.size()) {
		if (learn_lib_file.size() && !ReportOnly) learn_from_library(learn_lib_file);
		if (!fasta.load(fasta_files, fasta_filter_files)) Throw("Cannot load FASTA");
	}

	Library lib;
	if (fasta_files.size()) {
		all_fastas = fasta_files[0];
		for (int i = 1; i < fasta_files.size(); i++) all_fastas += std::string("; ") + fasta_files[i];
		lib.fasta_names = all_fastas;
	}
	if (FastaSearch && fasta_files.size()) {
		if (lib_file.size() && lib.load(&(lib_file[0]))) {
			GuideLibrary = true;
			double default_iRT = (lib.iRT_min + lib.iRT_max) * 0.5;
			if (learn_lib_file.size()) for (auto &e : lib.entries) e.target.sRT = (InSilicoRTPrediction ? predict_irt(e.name) : default_iRT);
		}
		lib.load(fasta);
	} else {
		if (!lib_file.size()) return 0;
		if (!lib.load(&(lib_file[0]))) Throw("Cannot load spectral library");
		if (ExportProsit && !fasta_files.size()) {
			if (Verbose >= 1) Time(), dsout << "Preparing Prosit input from the spectral library provided\n";
			lib.export_prosit();
		}
	}
	if (!FastaSearch && fasta_files.size()) {
		if (lib.protein_ids.size() > 1 && !Reannotate) fasta.load_proteins(fasta_files);
		else if (!fasta.load(fasta_files, fasta_filter_files)) Throw("Cannot load FASTA");
	}
	if (fasta_files.size()) {
		if ((lib.protein_ids.size() > 1 || FastaSearch) && !Reannotate) {
			if (Verbose >= 1 && !FastaSearch) Time(), dsout << "Annotating library proteins with information from the FASTA database\n";
			annotate_library(lib, fasta);
		} else {
			if (Verbose >= 1) {
				if (!Reannotate) 
					Time(), dsout << "No protein information found in the library. Annotating library precursors with information from the FASTA database\n";
				else {
					Time(), dsout << "Reannotating library precursors with information from the FASTA database\n";
				}
				if (!lib.protein_ids.size() || Reannotate) {
					PG pg;
					lib.protein_ids.clear(); lib.protein_ids.push_back(pg);
					lib.names.clear(); lib.genes.clear(); lib.proteins.clear();
					for (auto &e : lib.entries) e.pid_index = 0;
				}
			}
			lib.load_protein_annotations(fasta);
			annotate_library(lib, fasta);
		}
		fasta.free();
	}
	if (GuideLibrary && FastaSearch && nnIter >= iN)
		ReportProteinQValue = 1.0, Time(),
		dsout << "WARNING: protein q-values cannot be calculated without NNs when a guide library is used for FASTA search, protein q-value filter reset to 1.0\n";
	if (!ReportOnly) lib.initialise(!QuantOnly);
	if (!FastaSearch && !lib.from_speclib && lib.gen_decoys && !PredictorSaved) lib.save(lib_file + std::string(".speclib"));

	Library ref;
	if (ref_file.size()) {
		if (!ref.load(&(ref_file[0]))) ref_file.clear();
		else ref.initialise(true);
	}

	if (ReportOnly) goto cross_run;
	if (FastaSearch && QuantOnly) goto gen_spec_lib;
	if (QuantOnly) goto quant_only;

#if (HASH > 0)
	dsout << "Library hash: " << lib.hash() << "\n";
#endif

	// first loop
	dsout << "\n";
	if (RTProfiling) if (Verbose >= 1) Time(), dsout << "First pass: collecting retention time information\n";
	for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
		if ((UseRTInfo && RTProfiling) || UseQuant) {
			std::ifstream in((*it) + std::string(".quant"), std::ifstream::binary);
			if (in.fail() || temp_folder.size()) in = std::ifstream(location_to_file_name(*it) + std::string(".quant"), std::ifstream::binary);
			if (in.is_open()) {
				if (IndividualReports) {
					Quant Q;
					Q.read(in, lib.entries.size());
					lib.info.load(&(lib), Q);
					lib.info.quantify();
					lib.quantify_proteins(TopN, ProteinQuantQvalue);
					lib.report((*it) + std::string(".tsv"));
					if (lib.genes.size()) lib.gene_report((*it) + std::string(".genes.tsv"));
				}
				in.close();
				continue;
			}
		}
		if (Verbose >= 1) Time(), dsout << "File #" << std::distance(ms_files.begin(), it) + 1 << "/" << ms_files.size() << "\n";

		Run run(std::distance(ms_files.begin(), it));
		if (!run.load(&((*it)[0]))) dsout << "ERROR: cannot load the file, skipping\n", failed_files.insert(*it);
		if (ref_file.size()) run.reference_run(&ref);
		run.load_library(&lib);
		run.process();

		if (GenRef) {
			GenRef = false, ref_file = gen_ref_file;
			if (!ref.load(&(ref_file[0]))) ref_file.clear();
			else ref.initialise(true);
		}
	}

	if (!RTProfiling) goto gen_spec_lib;

quant_only:
	if (QuantOnly) {
		if (Verbose >= 1) Time(), dsout << "Quantification\n";
		for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
			Quant Q; std::ifstream in(*it + std::string(".quant"), std::ifstream::binary);
			if (in.fail() || temp_folder.size()) in = std::ifstream(location_to_file_name(*it) + std::string(".quant"), std::ifstream::binary);
			Q.read(in, lib.entries.size()); in.close();

			int runN = Q.run_index;
			Run run(runN);
			if (!run.load(&((*it)[0]))) dsout << "ERROR: cannot load the file, skipping\n", failed_files.insert(*it);
			run.load_library(&lib);
			for (auto pos = Q.entries.begin(); pos != Q.entries.end(); pos++) {
				auto e = &(run.entries[pos->index]);
				e->target.info.resize(1);
				e->target.flags |= fFound;
				e->target.S = pos->window;
				e->target.info[0].apex = pos->pr.apex;
				e->target.info[0].RT = run.scan_RT[pos->pr.apex];
				e->target.info[0].peak_pos = pos->pr.peak;
				e->target.info[0].best_fragment = pos->pr.best_fragment;
				e->target.info[0].peak_width = pos->pr.peak_width;
				e->target.info[0].qvalue = pos->pr.qvalue;
				e->target.info[0].protein_qvalue = pos->pr.protein_qvalue;
			}
			run.Q1Correction[0] = Q.Q1Correction[0], run.Q1Correction[1] = Q.Q1Correction[1];
			run.MassAccuracy = Q.MassAccuracy, run.MassAccuracyMs1 = Q.MassAccuracyMs1;
			run.MassCorrection = Q.MassCorrection, run.MassCorrectionMs1 = Q.MassCorrectionMs1;
			run.MassCalSplit = Q.MassCalSplit, run.MassCalSplitMs1 = Q.MassCalSplitMs1;
			run.MassCalCenter = Q.MassCalCenter, run.MassCalCenterMs1 = Q.MassCalCenterMs1;
			run.process();
		}
		lib.info.clear();
		goto cross_run;
	}

	{
		// second loop: refine iRT predictions
		if (Verbose >= 1) Time(), dsout << "Second pass: indentification and quantification\n";
		Profile profile(ms_files, lib.entries.size());
		if (RTProfiling) for (auto jt = profile.entries.begin(); jt != profile.entries.end(); jt++)
			if (jt->index >= 0) if (jt->pr.qvalue < MaxProfilingQvalue) {
				lib.entries[jt->pr.index].target.iRT = lib.entries[jt->pr.index].decoy.iRT = (!FastaSearch || GuideLibrary || RTLearnLib) ? jt->pr.predicted_iRT : jt->pr.RT;
				lib.entries[jt->pr.index].best_run = jt->pr.run_index;
				lib.entries[jt->pr.index].qvalue = jt->pr.qvalue;
				lib.entries[jt->pr.index].protein_qvalue = jt->pr.protein_qvalue;
			}

		for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
			if (UseQuant) {
				std::ifstream in((*it) + std::string(".quant"), std::ifstream::binary);
				if (in.fail() || temp_folder.size()) in = std::ifstream(location_to_file_name(*it) + std::string(".quant"), std::ifstream::binary);
				if (in.is_open()) {
					if (IndividualReports) {
						Quant Q;
						Q.read(in, lib.entries.size());
						lib.info.load(&(lib), Q);
						lib.info.quantify();
						lib.quantify_proteins(TopN, ProteinQuantQvalue);
						lib.report((*it) + std::string(".tsv"));
						if (lib.genes.size()) lib.gene_report((*it) + std::string(".genes.tsv"));
					}
					in.close();
					continue;
				}
			}
			if (Verbose >= 1) Time(), dsout << "Second pass: file #" << std::distance(ms_files.begin(), it) + 1 << "/" << ms_files.size() << "\n";

			Run run(std::distance(ms_files.begin(), it));
			if (!run.load(&((*it)[0]))) dsout << "ERROR: cannot load the file, skipping\n", failed_files.insert(*it);
			if (ref_file.size()) run.reference_run(&ref);
			run.load_library(&lib);
			run.process();

			if (GenRef) {
				GenRef = false, ref_file = gen_ref_file;
				if (!ref.load(&(ref_file[0]))) ref_file.clear();
				else ref.initialise(true);
			}
		}
	}

cross_run:
gen_spec_lib:
	if (failed_files.size()) {
		dsout << "ERROR: DIA-NN tried but failed to load the following files: ";
		int cnt = 0;
		for (auto &f : failed_files) {
			dsout << f;
			cnt++;
			if (cnt != failed_files.size()) dsout << ", ";
		}
		dsout << "\n";
	}

	if (IndividualReports || !(out_file.size() || out_gene_file.size()) || !ms_files.size()) goto gen_lib;
	if (Verbose >= 1) Time(), dsout << "Cross-run analysis\n";
	if (QuantInMem) lib.info.load(&(lib), ms_files, &quants);
	else lib.info.load(&(lib), ms_files);
	lib.info.quantify();
	if (GenFrExclusionInfo) {
		auto ann_lib = remove_extension(lib_file) + std::string("_fr_exclusion.tsv");
		lib.save(ann_lib, NULL, false, false);
	}
	lib.quantify_proteins(TopN, ProteinQuantQvalue);
	if (out_file.size()) {
		lib.report(out_file), lib.stats_report(remove_extension(out_file) + std::string(".stats.tsv"));
		if (Visualise.size()) lib.XIC_report(remove_extension(out_file) + std::string(".XIC.tsv"));
	}
	if (out_gene_file.size()) lib.gene_report(out_gene_file);

gen_lib:
	if (out_lib_file.size() && ms_files.size()) if (FastaSearch || GenSpecLib) { // generating spectral library
		MaxF = INF;
		if (Verbose >= 1) Time(), dsout << "Generating spectral library:\n";

		for (auto &le : lib.entries) le.best_run = -1;
		// determine the best run for each precursor and collect information from this run
		{
			Profile profile(ms_files);
			for (auto jt = profile.entries.begin(); jt != profile.entries.end(); jt++) if (jt->index >= 0) {
				if (jt->pr.qvalue <= ReportQValue) {
					lib.entries[jt->pr.index].qvalue = jt->pr.qvalue;
					if (ProfileQValue && jt->pr.profile_qvalue > jt->pr.qvalue) lib.entries[jt->pr.index].qvalue = jt->pr.profile_qvalue;
					if (lib.entries[jt->pr.index].qvalue > ReportQValue) lib.entries[jt->pr.index].best_run = -1;
					else {
						lib.entries[jt->pr.index].best_run = jt->pr.run_index;
						lib.entries[jt->pr.index].target.iRT = (GuideLibrary || (GenSpecLib && !FastaSearch) || RTLearnLib) ? jt->pr.predicted_iRT : jt->pr.RT;
						lib.entries[jt->pr.index].peak = jt->pr.peak;
						lib.entries[jt->pr.index].apex = jt->pr.apex;
						lib.entries[jt->pr.index].protein_qvalue = jt->pr.protein_qvalue;
						lib.entries[jt->pr.index].best_fr_mz = jt->pr.best_fr_mz;
						lib.entries[jt->pr.index].window = jt->window;
					}
				} else lib.entries[jt->pr.index].best_run = -1;
			}
		}

		if (InferPGs) {
			if (QuantInMem) lib.info.load(&(lib), ms_files, &quants);
			else lib.info.load(&(lib), ms_files);
			lib.infer_proteins();
			lib.info.clear();
			PGsInferred = true;
		}

		// extract spectra from runs
		for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
			Quant Qt;
			auto Q = &Qt;
			if (!QuantInMem) {
				std::ifstream in((*it) + std::string(".quant"), std::ifstream::binary);
				if (in.fail() || temp_folder.size()) in = std::ifstream(location_to_file_name(*it) + std::string(".quant"), std::ifstream::binary);
				if (in.is_open()) {
					Q->read_meta(in, lib.entries.size());
					in.close();
				}
			} else Q = &(quants[std::distance(ms_files.begin(), it)]);

			Run run(Q->run_index);
			if (!run.load(&((*it)[0]))) dsout << "ERROR: cannot load the file, skipping\n", failed_files.insert(*it);
			run.load_library(&lib);
			run.full_spectrum = true;

			run.Q1Correction[0] = Q->Q1Correction[0], run.Q1Correction[1] = Q->Q1Correction[1];
			run.MassAccuracy = Q->MassAccuracy, run.MassAccuracyMs1 = Q->MassAccuracyMs1;
			run.MassCorrection = Q->MassCorrection, run.MassCorrectionMs1 = Q->MassCorrectionMs1;
			run.MassCalSplit = Q->MassCalSplit, run.MassCalSplitMs1 = Q->MassCalSplitMs1;
			run.MassCalCenter = Q->MassCalCenter, run.MassCalCenterMs1 = Q->MassCalCenterMs1;

			run.update_library();
		}

		lib.save(out_lib_file, NULL, true, false);

		if (Verbose >= 1) Time(), dsout << "Loading the generated library and saving it in the .speclib format\n";
		Library new_lib;
		new_lib.load(&(out_lib_file[0]));
		Fasta fst; fst.load_proteins(fasta_files);
		annotate_library(new_lib, fst);
		new_lib.fasta_names = all_fastas;
		new_lib.save(out_lib_file + std::string(".speclib"));
	}

	remove:
#ifdef CPP17
	if (RemoveQuant & ms_files.size()) {
		if (Verbose >= 1) Time(), dsout << "Removing .quant files\n";
		for (auto &file : ms_files) try { std::experimental::filesystem::remove(file + std::string(".quant")); } catch (std::exception &e) { }
	}
#endif

	if ((out_file.size() && ms_files.size()) || out_lib_file.size()) {
		auto log_name = (out_file.size() && ms_files.size()) ? (remove_extension(out_file) + std::string(".log.txt")) : (remove_extension(out_lib_file) + std::string(".log.txt"));
		std::ofstream out(log_name, std::ofstream::out);
		if (out.fail()) dsout << "ERROR: cannot write log to " << log_name << "\n";
		else {
			out << dsout.s.rdbuf();
			out.close();
			Time(), dsout << "Log saved to " << log_name << "\n";
		}
	}

	dsout << "Finished\n\n";
	return 0;
}

