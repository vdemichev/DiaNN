/*
Copyright 2018, Vadim Demichev

This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license,
visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/

#define _HAS_ITERATOR_DEBUGGING 0
#define _ITERATOR_DEBUG_LEVEL 0  
#define _CRT_SECURE_NO_WARNINGS

#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4305)

#ifdef _MSC_VER
#define INPUT INPUT_
#include <windows.h>
#undef INPUT
#endif

// comment if no MKL installation available
#define CRANIUM_USE_MKL

#include "../cranium/src/cranium.h"
#include "../eigen/Eigen/Dense"
#include "../eigen/Eigen/Sparse"
#include <filesystem> // requires C++17
#include <iostream>
#include <fstream>
#include <sstream>
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

#define MZML

#ifdef MZML
#include "MSReader.h"
#endif

#define Throw(x) { std::cerr << __FILE__ << ": " << __LINE__ << ": " << x << "\n"; throw std::invalid_argument(x); }
#define Warning(x) { if (Verbose >= 1) std::cout << "WARNING: " << x << "\n"; }

#define Min(x, y) ((x) < (y) ? (x) : (y))
#define Max(x, y) ((x) > (y) ? (x) : (y))
#define Abs(x) ((x) >= 0 ? (x) : (-(x)))
#define Sqr(x) ((x) * (x))

const double E = 0.000000001;
const double INF = 10000000.0;

int Threads = 1;
int iN = 12;

bool BatchMode = false;
int MinBatch = 5000;
int MaxBatches = 1000;
int Batches = 1;

const int MaxLibSize = 1000000;

bool RTProfiling = false;

bool TestDataset = true;
double TestDatasetSize = 0.25;

int ScanFactor = 150;
int WindowRadius = 0;
double ScanScale = 2.2;
bool InferWindow = true;
bool IndividualWindows = false;
bool Calibrate = true;

int MaxDpPP = INF;

const int MaxCycleLength = 1000;
const double MinPeakHeight = 0.01;

int nRTTargetTopPrecursors = 50;
int nRTTopPrecursors;
int nRTRefTopPrecursors = 250;
double nRTMaxQvalue = 0.001; // only precursors with lower q-value are used for nRT estimation
double MassCalQvalue = 0.001;
double MassCalMs1Corr = 0.95;
int MassCalBins = 1; 
int MassCalBinsMs1 = 1;
int MinMassCalBin = 500;

int RTSegments = 20;
int MinRTPredBin = 20;

bool RefCal = false;
int nRTWindowIter = 4;
const int nRTRefWindowFactor = 10;
const int nRTWindowFactor = 20;
const double nRTWindowLoss = 0.05;
const double nRTWindowMargin = 1.5;
int nRTWinSearchMinRef = 11;
int nRTWinSearchMinCal = 100;

int nnIter = 10;
int nnBagging = 12;
int nnEpochs = 50;
double nnLearning = 0.0005;
bool GlobalNN = false;
bool GlobalTraining = false;
int nnGlobalEpochs = 50;
int nnHidden = 2;
int nnGlobalHidden = 5;
float Regularisation = 0.0;

bool nnStandardise = false;
const double nnStScale = 0.5;

bool nRTWindowedSearch = true;

const int CalibrationIter = 4;
bool ForceMassAcc = false;
bool IndividualMassAcc = false;
double CalibrationMassAccuracy = 100.0 / 1000000.0;
double GlobalMassAccuracy = 20.0 / 1000000.0; 
double GlobalMassAccuracyMs1 = 20.0 / 1000000.0;
double MinMassAcc = 1.0;
double MinMassAccMs1 = 1.0;
double GeneratorAccuracy = 5.0 / 1000000.0;
double MaxProfilingQvalue = 0.001;
double MaxQuantQvalue = 0.01;
double ProteinQuantQvalue = 0.01;
double PeakApexEvidence = 0.9; // force peak apexes to be close to the detection evidence maxima
double MinCorrScore = -INF;
double MaxCorrDiff = INF;
int MinMassDeltaCal = 50;
int MinMassAccCal = 100;
int MinFreeMassAcc = 500;
int MinCal = 5000;
int MinClassifier = 5000;

double FilterFactor = 1.5;

bool Convert = false; // .mzML to .dia conversion mode
bool UseRTInfo = false; // use .quant files created previously for RT profiling
bool UseQuant = false; // use .quant files created previously; implies UseRTInfo
bool QuantOnly = false; // quantification will be performed anew using identification info from .quant files
bool ReportOnly = false; // generate report from quant files
bool GenRef = false; // update the .ref file with newly computed data
std::string args, lib_file, learn_lib_file, fasta_file, out_file = "quant.tsv", ref_file, gen_ref_file, exclude_from_training, include_for_training;
std::vector<std::string> ms_files;
bool FastaSearch = false;

int Verbose = 1;
bool ExportWindows = false;
bool ExportLibrary = false;
bool GuideLibrary = false;
bool UseLibnRT = false;
bool InSilicoRTPrediction = true;
bool GenDecoys = true;

bool MissedCleavage = false, UniMod4 = true, UniMod35 = false, UniMod21 = false;
double MinFrIntensity = 0.01;
double MinFrMz = 300.0, MaxFrMz = 1800.0, MinOutFrMz = 300.0, MaxOutFrMz = 1800.0;
double MinPrMz = 400.0, MaxPrMz = 1200.0;
double SpectralLibraryGenQvalue = 0.01;

bool ExtraReportInfo = true;

enum {
	libPG, libPr, libCharge, libPrMz,
	libnRT, libFrMz, libFrI, libIsDecoy,
	libCols
};

std::vector<std::string> library_headers = {
	" UniProtIds UniprotID \"UniprotID\" ",
	" IntModifiedPeptide FullUniModPeptideName \"FullUniModPeptideName\" ",
	" PrecursorCharge Charge \"PrecursorCharge\" ",
	" PrecursorMz \"PrecursorMz\" ",
	" nRT iRT RetentionTime NormalizedRetentionTime Tr_recalibrated \"Tr_recalibrated\" ",
	" FragmentMz ProductMz \"ProductMz\" ",
	" RelativeIntensity RelativeFragmentIntesity LibraryIntensity \"LibraryIntensity\" ",
	" Decoy decoy \"Decoy\" \"decoy\" "
};

const double proton = 1.007825035;
const double OH = 17.003288;

std::vector<std::pair<std::string, float> > Modifications = {
	std::pair<std::string, float>("UniMod:4", (float)57.021464),
	std::pair<std::string, float>("Carbamidomethyl (C)", (float)57.021464),
	std::pair<std::string, float>("UniMod:5", (float)43.005814),
	std::pair<std::string, float>("Carbamylation (KR)", (float)43.005814),
	std::pair<std::string, float>("UniMod:7", (float)0.984016),
	std::pair<std::string, float>("Deamidation (NQ)", (float)0.984016),
	std::pair<std::string, float>("UniMod:35", (float)15.994915),
	std::pair<std::string, float>("Oxidation (M)", (float)15.994915),
	std::pair<std::string, float>("UniMod:1", (float)42.010565),
	std::pair<std::string, float>("Acetyl (Protein N-term)", (float)42.010565),
	std::pair<std::string, float>("Phosphorylation (ST)", (float)79.966331),
	std::pair<std::string, float>("UniMod:21", (float)79.966331),
	std::pair<std::string, float>("UniMod:259", (float)8.014199),
	std::pair<std::string, float>("UniMod:267", (float)10.008269),
	std::pair<std::string, float>("UniMod:268", (float)6.013809),
	std::pair<std::string, float>("UniMod:269", (float)10.027228)
};
std::vector<std::pair<std::string, std::string> > UniMod;
std::vector<std::pair<std::string, int> > UniModIndex;
std::vector<int> UniModIndices;

void init_unimod() {
	std::set<int> indices;
	auto mods = Modifications;
	std::sort(mods.begin(), mods.end());
	UniMod.resize(mods.size());
	UniModIndex.resize(mods.size());
	for (int i = 0; i < mods.size(); i++) {
		UniMod[i].first = UniModIndex[i].first = mods[i].first;
		if (!std::memcmp(&(mods[i].first[0]), "UniMod", 6)) UniMod[i].second = mods[i].first;
		else {
			double min = INF;
			int min_index;
			for (int j = 0; j < mods.size(); j++) if (j != i) {
				double delta = Abs(((double)mods[i].second) - (double)mods[j].second);
				if (delta < min) min = delta, min_index = j;
			}
			if (min > 0.00001 && Verbose >= 1) std::cout << "WARNING: potentially incorrect UniMod modification match for " << mods[i].first << ": " << min << "\n";
			UniMod[i].second = mods[min_index].first;
		}
		indices.insert(UniModIndex[i].second = std::stoi(&(UniMod[i].second[7])));
	}
	UniModIndices.insert(UniModIndices.begin(), indices.begin(), indices.end());
}

std::string to_unimod(std::string &mod) {
	for (int i = 0; i < UniMod.size(); i++) if (mod == UniMod[i].first) return UniMod[i].second;
	Warning("cannot convert to UniMod: unknown modification");
	return "";
}

int unimod_index(std::string &mod) {
	for (int i = 0; i < UniModIndex.size(); i++) if (mod == UniModIndex[i].first) return UniModIndex[i].second;
	Warning("cannot convert to UniMod: unknown modification");
	return 0;
}

int unimod_index_number(int index) {
	return std::distance(UniModIndices.begin(), std::lower_bound(UniModIndices.begin(), UniModIndices.end(), index));
}

enum {
	outFile, outPG, outPGQ, outModSeq,
	outPrId, outCharge, outQv, outPrQ, 
	outPrQRaw, outRT, outnRT,
	outCols
};

std::vector<std::string> oh = { // output headers
	"File.Name",
	"Protein.Group",
	"PG.Quantity",
	"Modified.Sequence",
	"Precursor.Id",
	"Precursor.Charge",
	"Q.Value",
	"Precursor.Quantity",
	"Precursor.Quantity.Raw",
	"RT",
	"nRT"
};

bool Normalisation = true;
double NormalisationQvalue = 0.001;
double NormalisationPeptidesFraction = 0.4;

class Lock {
public:
	std::atomic_bool lock;
	Lock() :lock() {}

	Lock(const std::atomic_bool &flag) :lock(flag.load()) {}

	Lock(const Lock &another) :lock(another.lock.load()) {}

	Lock &operator=(const Lock &another) { lock.store(another.lock.load()); return *this; }

	bool set() { return !lock.exchange(true); }
	void free() { lock.store(false); }
};

std::vector<NN> net;
std::vector<Lock> netLock;
std::vector<std::vector<float*> > training, training_classes;
std::vector<float*> test, test_classes;
float target_nn[2] = { 1.0, 0.0 };
float decoy_nn[2] = { 0.0, 1.0 };

class Parameter {
public:
	int min_iter_seek, min_iter_learn;

	Parameter(int _min_iter_seek, int _min_iter_learn) { min_iter_seek = _min_iter_seek, min_iter_learn = _min_iter_learn; }
};

const int nnS = 2, nnW = (2 * nnS) + 1, nnF = 6;
enum {
	pTimeCorr, pMinCorr, pCos, pMs1TimeCorr, pRT, pnRT, pCharge,
	pRef,
	pSig = pRef + nnF,
	pCorr = pSig + nnF,
	pArray = pCorr + nnF,
	pN = pArray + nnF * nnW
};
int MaxF = nnF;

std::string trim(std::string& input) { return std::regex_replace(input, std::regex("^ +| +$"), ""); }

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

	double ex = 0.0, ey = 0.0, s2x = 0.0, s2y = 0.0, r = 0.0;
	for (int i = 0; i < n; i++) ex += x[i], ey += y[i];
	if (n) ex /= (double)n, ey /= (double)n;
	for (int i = 0; i < n; i++) s2x += Sqr(x[i] - ex), s2y += Sqr(y[i] - ey), r += (x[i] - ex) * (y[i] - ey);
	if (s2x < E || s2y < E) return 0.0;
	return r / sqrt(s2x * s2y);
}

std::vector<int> rx, ry;
inline double spearman_corr(double * x, double * y, int n) {
	rx.resize(n);
	ry.resize(n);
	for (int i = 0; i < n; i++) rx[i] = ry[i] = i;
	std::sort(rx.begin(), rx.end(), [&](const int &l, const int &r) { return x[l] < x[r]; });
	std::sort(ry.begin(), ry.end(), [&](const int &l, const int &r) { return y[l] < y[r]; });
	return corr(&(rx[0]), &(ry[0]), n);
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

template<class T> inline void smooth(T * dst, T * src, int n) {
	assert(n >= 2);

	dst[0] = (2.0 / 3.0) * src[0] + (1.0 / 3.0) * src[1];
	dst[n - 1] = (2.0 / 3.0) * src[n - 1] + (1.0 / 3.0) * src[n - 2];
	for (int i = 1; i < n - 1; i++)
		dst[i] = 0.5 * src[i] + 0.25 * (src[i - 1] + src[i + 1]);
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

std::vector<double> rt_p, rt_e, rt_s, rt_c;
typedef Eigen::Triplet<double> Triplet;
std::vector<Triplet> TrL;
template <class T, class Tx, class Tfunc> void spline(std::vector<T> * result, std::vector<Tx> * coo, std::vector<double> &coeff, std::vector<std::pair<Tx, Tfunc> > &data, int segments) {
	int i, k, m, n = data.size();
	auto pava = PAVA(data);
	double w;
	
	rt_p.resize(n);
	rt_e.resize(n);
	m = segments;
	coeff.resize(2 * m);
	rt_s.resize(m + 1), rt_c.resize(m + 1);
	for (i = 0; i < n; i++) rt_p[i] = data[i].second, rt_e[i] = data[i].first;

	double split = ((double)(n - 1)) / (double)m;
	for (i = 0; i < m; i++) rt_s[i] = rt_e[Min(split * (double)i, n - 1)]; rt_s[m] = rt_e[n - 1];
	for (i = 0; i < m; i++) rt_c[i] = 0.5 * (rt_s[i] + rt_s[i + 1]);
	for (i = 1; i < m; i++) if (rt_c[i] - rt_c[i - 1] < E) {
		for (k = i - 1; k < m - 1; k++) rt_c[k] = rt_c[k + 1];
		m--;
	}

	TrL.clear();
	for (i = k = 0; i < n; i++) {
		double x = rt_e[i];
		if (!k && x <= rt_c[0]) {
			TrL.push_back(Triplet(i, 0, 1.0));
			TrL.push_back(Triplet(i, 1, x - rt_c[0]));
		}
		else {
			for (; k < m && rt_c[k] < x; k++);
			if (k == m) {
				TrL.push_back(Triplet(i, (k - 1) * 2, 1.0));
				TrL.push_back(Triplet(i, (k - 1) * 2 + 1, x - rt_c[k - 1]));
			}
			else {
				w = rt_c[k] - rt_c[k - 1];
				double u = (x - rt_c[k - 1]) / w, v = u - 1.0;
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

	if (result == NULL || coo == NULL) return;

	result->clear(); result->resize(coo->size(), 0.0);
	for (i = k = 0; i < coo->size(); i++) {
		double x = (*coo)[i];
		if (!k && x <= rt_c[0]) {
			(*result)[i] += coeff[0];
			(*result)[i] += coeff[1] * (x - rt_c[0]);
		}
		else {
			for (; k < m && rt_c[k] < x; k++);
			if (k == m) {
				(*result)[i] += coeff[(k - 1) * 2];
				(*result)[i] += coeff[(k - 1) * 2 + 1] * (x - rt_c[k - 1]);
			}
			else {
				w = rt_c[k] - rt_c[k - 1];
				double u = (x - rt_c[k - 1]) / w, v = u - 1.0;
				(*result)[i] += coeff[(k - 1) * 2] * (1.0 + 2.0 * u) * Sqr(v);
				(*result)[i] += coeff[(k - 1) * 2 + 1] * w * u * Sqr(v);

				(*result)[i] += coeff[k * 2] * Sqr(u) * (1.0 - 2.0 * v);
				(*result)[i] += coeff[k * 2 + 1] * w * Sqr(u) * v;
			}
		}
	}
}

bool st_calculated = false;
double st_e[pN], st_sd[pN]; // for standardisation

class Score {
public:
    float x[pN];
    
    Score() { for (int i = 0; i < pN; i++) x[i] = 0.0; }
};

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
};

double AA[256];
int AA_index[256];
const char * AAs = "GAVLIFMPWSCTYHKRQEND";

void init_aas() {
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

	for (int i = 0; i < 20; i++)
		AA_index[AAs[i]] = i;
}

std::vector<double> yDelta;
int yDeltaS = 2, yDeltaW = (yDeltaS * 2 + 1);
inline int y_delta_composition_index(char aa) { return AA_index[aa]; }
inline int y_delta_charge_index() { return y_delta_composition_index(AAs[19]) + 1; }
inline int y_delta_index(char aa, int shift) {
	assert(Abs(shift) <= yDeltaS);
	return y_delta_charge_index() + AA_index[aa] * yDeltaW + shift + yDeltaS + 1;
}
int yCTermD = 8, yNTermD = 5;
inline int y_cterm_index(char aa, int shift) {
	assert(shift <= yCTermD - 1);
	int row = (aa == 'K' ? 0 : 1);
	return row * yCTermD + shift + y_delta_index(AAs[19], yDeltaS) + 1;
}
inline int y_nterm_index(int shift) {
	assert(shift <= yNTermD - 1);
	return shift + y_cterm_index('G', yCTermD - 1) + 1;
}
inline int y_delta_size() { return y_nterm_index(yNTermD - 1) + 1; }

inline double y_ratio(int i, int charge, std::string &aas) {
	int n = aas.length(), j;
	double v = yDelta[y_delta_charge_index()] * (double)charge;
	for (j = 0; j < n; j++) v += yDelta[y_delta_composition_index(aas[j])];
	for (j = Max(i - yDeltaS, 0); j <= Min(i + yDeltaS, n - 1); j++)
		v += yDelta[y_delta_index(aas[j], j - i)];
	if (i <= yNTermD) v += yDelta[y_nterm_index(i - 1)];
	if (i >= n - yCTermD) v += yDelta[y_cterm_index(aas[n - 1], n - 1 - i)];

	return v;
}

void y_scores(std::vector<double> &s, int charge, std::string &aas) {
	int n = aas.length(), i;
	s.resize(n);
	s[0] = s[1] = 0.0;
	for (i = 2; i < n; i++)
		s[i] = s[i - 1] + y_ratio(i, charge, aas);
}

void to_exp(std::vector<double> &s) { for (int i = 0; i < s.size(); i++) s[i] = exp(s[i]); }

std::vector<double> InSilicoRT;
int RTTermD = 1, RTTermDScaled = 5;
inline int aa_rt_nterm(char aa, int shift) { return 1 + Min(shift, RTTermD - 1)  * 20 + AA_index[aa]; }
inline int aa_rt_cterm(char aa, int shift) { return aa_rt_nterm(AAs[19], RTTermD - 1) + 1 + Min(shift, RTTermD - 1) * 20 + AA_index[aa]; }
inline int aa_rt_nterm_scaled(char aa, int shift) { return aa_rt_cterm(AAs[19], RTTermD - 1) + 1 + Min(shift, RTTermDScaled - 1) * 20 + AA_index[aa]; }
inline int aa_rt_cterm_scaled(char aa, int shift) { return aa_rt_nterm_scaled(AAs[19], RTTermDScaled - 1) + 1 + Min(shift, RTTermDScaled - 1) * 20 + AA_index[aa]; }
inline int mod_rt_index(int mod) { return aa_rt_cterm_scaled(AAs[19], RTTermDScaled - 1) + 1 + unimod_index_number(mod); }
inline int mod_rt_index(std::string &mod) { return mod_rt_index(unimod_index(mod)); }
inline int in_silico_rt_size() { return mod_rt_index(UniModIndices[UniModIndices.size() - 1]) + 1; }

void init_prediction() {
	yDelta.resize(y_delta_size(), 0.0);
	InSilicoRT.resize(in_silico_rt_size(), 0.0);
}

inline int peptide_length(std::string &name) {
	int i, n = 0;
	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol >= 'A' && symbol <= 'Y') n++;
	}
	return n;
}

std::string get_aas(std::string &name) {
	int i;
	std::string result;

	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			int end = name.find(symbol == '(' ? ")" : "]", i + 1);
			i = end;
			continue;
		}
		result.push_back(symbol);
	}
	return result;
}

std::vector<double> get_sequence(std::string &name) {
	int i, j;
	double add = 0.0;
	std::vector<double> result;

	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i++;

			int end = name.find(symbol == '(' ? ")" : "]", i);
			if (end == std::string::npos) Throw(std::string("incorrect peptide name format: ") + name);

			std::string mod = name.substr(i, end - i);
			for (j = 0; j < Modifications.size(); j++) if (Modifications[j].first == mod) {
				if (!result.size()) add += Modifications[j].second;
				else result[result.size() - 1] += Modifications[j].second;
				break;
			}
			if (j == Modifications.size())
				Throw(std::string("unknown modification: ") + name);
			i = end;
			continue;
		}

		result.push_back(AA[symbol] + add);
		add = 0.0;
	}
	return result;
}

std::string to_canonical(std::string &name) {
	int i, j;
	std::string result;

	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i++;

			int end = name.find(symbol == '(' ? ")" : "]", i);
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
inline std::string to_canonical(std::string &name, int charge) { return to_canonical(name) + std::to_string(charge); }

void count_rt_aas(std::vector<double> &count, std::string &name) {
	int i, pos, l = peptide_length(name);
	double scale = 1.0 / (double)l;

	for (i = pos = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i++;

			int end = name.find(symbol == '(' ? ")" : "]", i);
			if (end == std::string::npos) Throw(std::string("incorrect peptide name format: ") + name);
			count[mod_rt_index(name.substr(i, end - i))]++;
			i = end;
			continue;
		}
		count[aa_rt_nterm(symbol, pos)] += 1.0;
		count[aa_rt_cterm(symbol, l - 1 - pos)] += 1.0;
		count[aa_rt_nterm_scaled(symbol, pos)] += scale;
		count[aa_rt_cterm_scaled(symbol, l - 1 - pos)] += scale;
		pos++;
	}
}

double predict_nrt(std::string &name) {
	int i, pos, l = peptide_length(name);
	double nRT = InSilicoRT[0], scale = 1.0 / (double)l;

	for (i = pos = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol < 'A' || symbol > 'Z') {
			if (symbol != '(' && symbol != '[') continue;
			i++;

			int end = name.find(symbol == '(' ? ")" : "]", i);
			if (end == std::string::npos) Throw(std::string("incorrect peptide name format: ") + name);
			nRT += InSilicoRT[mod_rt_index(name.substr(i, end - i))];
			i = end;
			continue;
		}
		nRT += InSilicoRT[aa_rt_nterm(symbol, pos)];
		nRT += InSilicoRT[aa_rt_cterm(symbol, l - 1 - pos)];
		nRT += InSilicoRT[aa_rt_nterm_scaled(symbol, pos)] * scale;
		nRT += InSilicoRT[aa_rt_cterm_scaled(symbol, l - 1 - pos)] * scale;
		pos++;
	}
	return nRT;
}

void arguments(int argc, char *argv[]) {
	int i, start, next, end;
	std::string prefix, ext;
	std::vector<std::string> files;

	for (i = 0; i < argc; i++) args += argv[i] + std::string(" ");
	args = std::regex_replace(args, std::regex("\n"), " ");
	start = args.find("--", 0);
	while (start != std::string::npos) {
		start += 2;
		next = args.find("--", start);
		end = (next == std::string::npos) ? args.length() : next;

		if (!memcmp(&(args[start]), "cfg", 3)) { // parse the contents of the file
			std::string name = args.substr(start + 4, end - start - 4);
			std::ifstream in(name);
			if (in.is_open()) {
				std::stringstream stream;
				stream << in.rdbuf();
				args = stream.str() + std::string(" ")
					+ args.substr(end, std::string::npos);
				in.close();

				args = std::regex_replace(args, std::regex("\n"), " ");
				start = args.find("--", 0);
				continue;
			}
			else std::cerr << "failed to open cfg file\n";
		}
		else if (!memcmp(&(args[start]), "threads ", 8)) Threads = std::stoi(args.substr(start + 8, std::string::npos)), std::cout << "Thread number set to " << Threads << "\n";
		else if (!memcmp(&(args[start]), "batch-mode ", 11)) BatchMode = true, std::cout << "Batch mode will be used to reduce memory usage. NN classifier can be used for the last iteration only.\n";
		else if (!memcmp(&(args[start]), "verbose ", 8)) Verbose = std::stoi(args.substr(start + 8, std::string::npos));
		else if (!memcmp(&(args[start]), "export-windows ", 15)) ExportWindows = true;
		else if (!memcmp(&(args[start]), "export-library ", 15)) ExportLibrary = true;
		else if (!memcmp(&(args[start]), "f ", 2)) files.push_back(trim(args.substr(start + 2, end - start - 2)));
		else if (!memcmp(&(args[start]), "dir ", 4)) {
			string dir = trim(args.substr(start + 4, end - start - 4));
			for (auto &file : std::experimental::filesystem::directory_iterator(dir)) files.push_back(file.path().string());
		} 
		else if (!memcmp(&(args[start]), "lib ", 4)) lib_file = trim(args.substr(start + 4, end - start - 4));
		else if (!memcmp(&(args[start]), "fasta ", 6)) fasta_file = trim(args.substr(start + 6, end - start - 6));
		else if (!memcmp(&(args[start]), "ref ", 4)) ref_file = trim(args.substr(start + 4, end - start - 4));
		else if (!memcmp(&(args[start]), "out ", 4)) out_file = trim(args.substr(start + 4, end - start - 4));
		else if (!memcmp(&(args[start]), "learn-lib ", 10)) learn_lib_file = trim(args.substr(start + 10, end - start - 10));
		else if (!memcmp(&(args[start]), "use-lib-nrt ", 12)) UseLibnRT = true, std::cout << "Spectral library will be saved with normalised RTs used instead of experimental RTs\n";
		else if (!memcmp(&(args[start]), "library-headers ", 16)) {
			std::string word;
			std::stringstream list(trim(args.substr(start + 16, end - start - 16)));
			for (i = 0; std::getline(list, word, ','); i++) {
				if (i >= libCols) {
					std::cout << "WARNING: " << word << ": extra headers will be ignored\n";
					break;
				}
				library_headers[i] = std::string(" ") + word + std::string(" ");
			}
		}
		else if (!memcmp(&(args[start]), "output-headers ", 15)) {
			std::string word;
			std::stringstream list(trim(args.substr(start + 15, end - start - 15)));
			for (i = 0; std::getline(list, word, ','); i++) {
				if (i >= outCols) {
					std::cout << "WARNING: " << word << ": extra headers will be ignored\n";
					break;
				}
				oh[i] = std::string(" ") + word + std::string(" ");
			}
		}
		else if (!memcmp(&(args[start]), "mod ", 4)) {
			std::string name, mass;
			std::stringstream list(trim(args.substr(start + 4, end - start - 4)));
			if (!std::getline(list, name, ',')) std::cout << "WARNING: no modification name, modification ignored\n";
			else if (!std::getline(list, mass, ',')) std::cout << "WARNING: no modification mass, modification ignored\n";
			else Modifications.push_back(std::pair<std::string, float>(name, (float)std::stod(mass)));
		}
		else if (!memcmp(&(args[start]), "exclude ", 8)) exclude_from_training = trim(args.substr(start + 8, end - start - 8)),
			std::cout << "Precursors corresponding to proteins with [" << exclude_from_training << "] in the ID will be ignored when training the classifier\n";
		else if (!memcmp(&(args[start]), "include ", 8)) include_for_training = trim(args.substr(start + 8, end - start - 8)),
			std::cout << "Precursors corresponding to proteins without [" << include_for_training << "] in the ID will be ignored when training the classifier\n";
		else if (!memcmp(&(args[start]), "ref-cal ", 8)) RefCal = true, std::cout << "Reference peptides will be used for calibration\n";
		else if (!memcmp(&(args[start]), "gen-ref ", 8)) GenRef = true, gen_ref_file = trim(args.substr(start + 8, end - start - 8)),
			std::cout << "A library of reference peptides will be generated\n";
		else if (!memcmp(&(args[start]), "window ", 7)) {
			WindowRadius = std::stoi(args.substr(start + 7, std::string::npos));
			if (WindowRadius <= 0) std::cout << "WARNING: scan window radius should be a positive integer\n";
			else InferWindow = false, std::cout << "Scan window radius set to " << WindowRadius << "\n";
		}
		else if (!memcmp(&(args[start]), "no-rt-window ", 13)) nRTWindowedSearch = false, std::cout << "Full range of retention times will be considered\n";
		else if (!memcmp(&(args[start]), "no-window-inference ", 20)) InferWindow = false, std::cout << "Scan window inference turned off\n";
		else if (!memcmp(&(args[start]), "individual-windows ", 19)) IndividualWindows = true, std::cout << "Scan windows will be inferred separately for different runs\n";
		else if (!memcmp(&(args[start]), "individual-mass-acc ", 20)) IndividualMassAcc = true, std::cout << "Mass accuracy will be determined separately for different runs\n";
		else if (!memcmp(&(args[start]), "convert ", 8)) Convert = true, std::cout << ".mzML to .dia conversion\n";
		else if (!memcmp(&(args[start]), "use-rt ", 7)) UseRTInfo = true, std::cout << "Existing .quant files will be used for RT profiling\n";
		else if (!memcmp(&(args[start]), "use-quant ", 10)) UseQuant = true, std::cout << "Existing .quant files will be used\n";
		else if (!memcmp(&(args[start]), "quant-only ", 11)) QuantOnly = true, std::cout << "Quantification will be performed anew using existing identification info\n";
		else if (!memcmp(&(args[start]), "report-only ", 12)) ReportOnly = true, std::cout << "Report will be generated using .quant files\n";
		else if (!memcmp(&(args[start]), "iter ", 5)) iN = Max(CalibrationIter + 4, std::stoi(args.substr(start + 5, std::string::npos))),
			std::cout << "Number of iterations set to " << iN << "\n";
		else if (!memcmp(&(args[start]), "profiling-qvalue ", 17)) MaxProfilingQvalue = std::stod(args.substr(start + 17, std::string::npos)),
			std::cout << "RT profiling q-value threshold set to " << MaxProfilingQvalue << "\n";
		else if (!memcmp(&(args[start]), "quant-qvalue ", 13)) MaxQuantQvalue = std::stod(args.substr(start + 13, std::string::npos)),
			std::cout << "Q-value threshold for cross-run quantification set to " << MaxQuantQvalue << "\n";
		else if (!memcmp(&(args[start]), "lib-qvalue ", 11)) SpectralLibraryGenQvalue = std::stod(args.substr(start + 11, std::string::npos)),
			std::cout << "Q-value threshold spectral library generation set to " << SpectralLibraryGenQvalue << "\n";
		else if (!memcmp(&(args[start]), "rt-profiling ", 13)) RTProfiling = true, std::cout << "RT profiling turned on\n";
		else if (!memcmp(&(args[start]), "prefix ", 7)) prefix = trim(args.substr(start + 7, end - start - 7)); // prefix added to input file names
		else if (!memcmp(&(args[start]), "ext ", 4)) ext = trim(args.substr(start + 4, end - start - 4)); // extension added to input file names
		else if (!memcmp(&(args[start]), "no-test-dataset ", 16)) TestDataset = false, std::cout << "Library will not be split into training and test datasets\n";
		else if (!memcmp(&(args[start]), "test-proportion ", 16)) TestDatasetSize = std::stod(args.substr(start + 16, std::string::npos)),
			std::cout << "The " << TestDatasetSize << " fraction of the precursors will be used as the test dataset\n";
		else if (!memcmp(&(args[start]), "no-nn ", 6)) nnIter = INF, std::cout << "Neural network classifier turned off\n";
		else if (!memcmp(&(args[start]), "global-nn ", 10)) GlobalNN = GlobalTraining = true, std::cout << "Global NN training will be used\n";
		else if (!memcmp(&(args[start]), "nn-iter ", 8)) nnIter = Max(CalibrationIter + 3, std::stoi(args.substr(start + 8, std::string::npos))),
			std::cout << "Neural network classifier will be used starting from the interation number " << nnIter << "\n";
		else if (!memcmp(&(args[start]), "nn-bagging ", 11)) nnBagging = std::stoi(args.substr(start + 11, std::string::npos)),
			std::cout << "Neural network bagging set to " << nnBagging << "\n";
		else if (!memcmp(&(args[start]), "nn-epochs ", 10)) nnEpochs = std::stoi(args.substr(start + 10, std::string::npos)),
			std::cout << "Neural network epochs number set to " << nnEpochs << "\n";
		else if (!memcmp(&(args[start]), "nn-global-epochs ", 17)) nnGlobalEpochs = std::stoi(args.substr(start + 17, std::string::npos)),
			std::cout << "Neural network global epochs number set to " << nnGlobalEpochs << "\n";
		else if (!memcmp(&(args[start]), "nn-learning-rate ", 17)) nnLearning = std::stod(args.substr(start + 17, std::string::npos)),
			std::cout << "Neural network learning rate set to " << nnLearning << "\n";
		else if (!memcmp(&(args[start]), "nn-reg ", 7)) Regularisation = std::stod(args.substr(start + 7, std::string::npos)),
			std::cout << "Neural network regularisation set to " << Regularisation << "\n";
		else if (!memcmp(&(args[start]), "nn-hidden ", 10)) nnHidden = Max(1, std::stoi(args.substr(start + 10, std::string::npos))),
			std::cout << "Number of hidden layers set to " << nnHidden << "\n";
		else if (!memcmp(&(args[start]), "nn-global-hidden ", 17)) nnGlobalHidden = Max(1, std::stoi(args.substr(start + 17, std::string::npos))),
			std::cout << "Number of hidden layers in the global NN set to " << nnGlobalHidden << "\n";
		else if (!memcmp(&(args[start]), "standardise ", 12)) nnStandardise = true, std::cout << "Standardisation of scores will be performed\n";
		else if (!memcmp(&(args[start]), "mass-acc-cal ", 13)) CalibrationMassAccuracy = std::stod(args.substr(start + 13, std::string::npos)) / 1000000.0,
			std::cout << "Calibration mass accuracy set to " << CalibrationMassAccuracy << "\n";
		else if (!memcmp(&(args[start]), "mass-acc ", 9)) GlobalMassAccuracy = std::stod(args.substr(start + 9, std::string::npos)) / 1000000.0,
			std::cout << "Mass accuracy set to " << GlobalMassAccuracy << "\n";
		else if (!memcmp(&(args[start]), "mass-acc-ms1 ", 13)) GlobalMassAccuracyMs1 = std::stod(args.substr(start + 13, std::string::npos)) / 1000000.0,
			std::cout << "MS1 mass accuracy set to " << GlobalMassAccuracyMs1 << "\n";
		else if (!memcmp(&(args[start]), "force-mass-acc ", 15)) ForceMassAcc = true, std::cout << "Mass accuracy will be fixed to the specified values\n";
		else if (!memcmp(&(args[start]), "gen-acc ", 8)) GeneratorAccuracy = std::stod(args.substr(start + 8, std::string::npos)) / 1000000.0,
			std::cout << "Fragmentation spectrum generator accuracy set to " << GeneratorAccuracy << "\n";
		else if (!memcmp(&(args[start]), "min-corr ", 9)) MinCorrScore = Min(std::stod(args.substr(start + 9, std::string::npos)), -1.1 + (double)nnF),
			std::cout << "Only peaks with correlation sum exceeding " << MinCorrScore << " will be considered\n";
		else if (!memcmp(&(args[start]), "corr-diff ", 10)) MaxCorrDiff = Max(std::stod(args.substr(start + 10, std::string::npos)), E),
			std::cout << "Peaks with correlation sum below " << MaxCorrDiff << " from maximum will not be considered\n";
		else if (!memcmp(&(args[start]), "peak-apex ", 10)) PeakApexEvidence = Min(std::stod(args.substr(start + 10, std::string::npos)), 0.99),
			std::cout << "Peaks must have apex height which is at least " << MinCorrScore << " of the maximum in the m/z scan window\n";
		else if (!memcmp(&(args[start]), "max-dppp ", 9)) MaxDpPP = Max(1, std::stoi(args.substr(start + 9, std::string::npos))),
			std::cout << "Maximum number of data points per peak set to " << MaxDpPP << "\n";
		else if (!memcmp(&(args[start]), "no-norm ", 8)) Normalisation = false, std::cout << "Cross-run normalisation turned off\n";
		else if (!memcmp(&(args[start]), "norm-qvalue ", 12)) NormalisationQvalue = std::stod(args.substr(start + 11, std::string::npos)),
			std::cout << "Q-value threshold for cross-run normalisation set to " << NormalisationQvalue << "\n";
		else if (!memcmp(&(args[start]), "norm-fraction ", 14)) NormalisationPeptidesFraction = std::stod(args.substr(start + 14, std::string::npos)),
			std::cout << "Normalisation peptides fraction set to " << NormalisationPeptidesFraction << "\n";
		else if (!memcmp(&(args[start]), "no-calibration ", 15)) Calibrate = false, std::cout << "Mass calibration turned off\n";
		else if (!memcmp(&(args[start]), "mass-cal-bins ", 14)) MassCalBins = MassCalBinsMs1 = std::stoi(args.substr(start + 14, std::string::npos)),
			std::cout << "Maximum number of mass calibration bins set to " << MassCalBins << "\n";
		else std::cerr << "WARNING: unrecognised option [--" << trim(args.substr(start, end - start)) << "]\n";

		start = next;
	}
	if (prefix.length()) for (auto it = files.begin(); it != files.end(); it++) *it = prefix + *it;
	if (ext.length()) for (auto it = files.begin(); it != files.end(); it++) *it += ext;
	for (int i = 0; i < files.size(); i++) {
		int pos = files[i].find_last_of('.'), j;
		auto extension = files[i].substr(pos);
		std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
		if (extension == std::string(".dia")) {
			if (!Convert) ms_files.push_back(files[i]);
			continue;
		}
		auto converted = files[i] + std::string(".dia");
		for (j = 0; j < files.size(); j++) if (j != i) if (files[j] == converted) break;
		if (j == files.size()) ms_files.push_back(files[i]);
	}
	if (ms_files.size() < 2) RTProfiling = Normalisation = false;
	if (UseQuant || QuantOnly) UseRTInfo = true;
	if (GlobalNN) TestDataset = false;
	nnIter = Min(Max(Max(nRTWindowIter, CalibrationIter), nnIter), iN);
	if (BatchMode) nnIter = Max(iN - 1, nnIter);
	MinBatch = Max(MinBatch, Min(MinClassifier, MinCal));
		
	net.resize(nnBagging);
	netLock.resize(nnBagging);
	for (i = 0; i < nnBagging; i++) net[i].network = NULL, net[i].seed(i);

	if (fasta_file.size()) {
		FastaSearch = true;
		MinCorrScore = Max(MinCorrScore, 1.0); // essential to reduce memory usage
		PeakApexEvidence = 0.99;
		nnIter = Max(nnIter, iN - 1);
	}
	if (!learn_lib_file.size() && !GuideLibrary) UseLibnRT = false;

	init_unimod();
	init_prediction();
}

enum {
	type_b, type_y,
	type_N
};

enum {
	loss_none, loss_H2O, loss_NH3, loss_CO,
	loss_N
};

const double Loss[loss_N] = { 0.0, 18.011113035, 17.026548, 27.994915 };

char char_from_type[2] = { 'b', 'y' };
std::string name_from_loss[4] = { std::string("noloss"), std::string("H2O"), std::string("NH3"), std::string("CO") };

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

	void init(Ion &other) {
		type = other.type;
		index = other.index;
		loss = other.loss;
		charge = other.charge;
		mz = other.mz;
	}

	bool operator == (const Ion &other) { return index == other.index && type == other.type && charge == other.charge && loss == other.loss; }
};

std::vector<Ion> generate_fragments(std::vector<double> &sequence, int charge, int loss, int *cnt) {
	int i;
	double c = (double)charge, curr, s = sum(&(sequence[0]), sequence.size());
	std::vector<Ion> result;

	for (i = 0, curr = 0.0; i < sequence.size() - 1; i++) {
		curr += sequence[i];
		double b = curr;
		double y = s - curr + proton + OH;

		result.push_back(Ion(type_b, i + 1, loss, charge, (b + c * proton - Loss[loss]) / c)), (*cnt)++;
		result.push_back(Ion(type_y, i + 1, loss, charge, (y + c * proton - Loss[loss]) / c)), (*cnt)++;
	}
	return result;
}

std::vector<Ion> recognise_fragments(std::vector<double> &sequence, std::vector<Peak> &fragments, bool full_spectrum = false, int max_charge = 19, int loss_cap = loss_N) {
	int i, j, cnt, tot, charge, loss, index;
	std::vector<Ion> result(fragments.size());
	std::vector<Ion> v;
	double delta, min;

anew:
	cnt = tot = index = charge = 0;
	loss = loss_none;
	for (i = 0; i < result.size(); i++) result[i].charge = 0;

start:
	v = generate_fragments(sequence, charge, loss, &cnt);
	for (i = 0; i < fragments.size(); i++) {
		if (result[i].charge) continue;
		double margin = fragments[i].mz * GeneratorAccuracy;
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
		GeneratorAccuracy *= 2.0;
		if (Verbose >= 1) std::cerr << "WARNING: not all fragments recognised; generator accuracy threshold increased to " << GeneratorAccuracy << "\n";
		v.clear();
		goto anew;
	}

	return result;
}

std::vector<Ion> generate_fragments(std::vector<double> &sequence, std::vector<Ion> pattern) {
	int i, j, charge = 1, loss = loss_none, cnt = 0, tot = 0;
	std::vector<Ion> result(pattern.size());
	std::vector<Ion> v;

start:
	v = generate_fragments(sequence, charge, loss, &cnt); i = 0;
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

class Tandem {
public:
    double RT, window_low, window_high;
    std::vector<Peak> peaks;
    inline int size() { return peaks.size(); }
    inline void resize(int size) { peaks.resize(size); }
    
	Tandem() {

	}

#ifdef MZML
    Tandem(MSToolkit::Spectrum &s) {
        RT = s.getRTime();
        window_low = s.getSelWindowLower();
        window_high = s.getSelWindowUpper();
        
        resize(s.size());
        for (int i = 0; i < s.size(); i++) peaks[i].init(s.at(i).mz, s.at(i).intensity);
    }
#endif
    
    inline bool has(float mz) { return (window_low <= mz && window_high > mz); }
    
	template <bool get_mz> inline double level(float mz, float accuracy, float * peak_mz = NULL) {
		int i, low = 0, high = size();
		float v, s, margin, min, max;
		margin = mz * accuracy;
		min = mz - margin, max = mz + margin;

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
					for (i = middle - 1; i >= low && peaks[i].mz > min; i--) if (peaks[i].height > s) {
						s = peaks[i].height;
						if (get_mz) *peak_mz = peaks[i].mz;
					}
					return s;
				}
				low = middle + 1;
			}
			else high = middle;
		}
		return 0.0;
	}

	void write(std::ofstream &out) {
		out.write((char*)&RT, sizeof(double));
		out.write((char*)&window_low, sizeof(double));
		out.write((char*)&window_high, sizeof(double));

		int size = peaks.size();
		out.write((char*)&size, sizeof(int));
		out.write((char*)&(peaks[0]), peaks.size() * sizeof(Peak));
	}

	void read(std::ifstream &in) {
		in.read((char*)&RT, sizeof(double));
		in.read((char*)&window_low, sizeof(double));
		in.read((char*)&window_high, sizeof(double));

		int size; in.read((char*)&size, sizeof(int));
		peaks.resize(size);
		in.read((char*)&(peaks[0]), peaks.size() * sizeof(Peak));
	}
};

class Protein {
public:
    std::string name;
    mutable std::list<int> precursors; // precursor indices in the library
    
    Protein(std::string _name) {
        name = _name;
    }
    
    friend inline bool operator < (const Protein &left, const Protein &right) { return left.name < right.name; }
};

enum {
	qScaled, qTotal, qFiltered,
	qN
};

class Fragment {
public:
	mutable float quantity[qN];
	mutable float corr;
};

class PrecursorEntry {
public:
	int index, run_index;
	mutable bool decoy_found;
	mutable int apex, peak, best_fragment, peak_width;
	mutable float RT, nRT, predicted_nRT, qvalue, quantity, normalised, pg_quantity;

	friend inline bool operator < (const PrecursorEntry &left, const PrecursorEntry &right) { return left.index < right.index; }
};

struct ScoreEntry {
public:
	mutable float target_scores[pN];
	mutable float decoy_scores[pN];
};

class QuantEntry {
public:
	int index, stored_scores, window;
	mutable ScoreEntry sc;
    mutable PrecursorEntry pr;
    mutable std::vector<Fragment> fr;

	QuantEntry() { index = -1; }
    
    void write(std::ofstream &out) {
        int size = fr.size();
		out.write((char*)&index, sizeof(int));
		out.write((char*)&stored_scores, sizeof(int));
		out.write((char*)&window, sizeof(int));
        out.write((char*)&pr, sizeof(PrecursorEntry));
		if (stored_scores) out.write((char*)&sc, sizeof(ScoreEntry));
        out.write((char*)&size, sizeof(int));
        for (int i = 0; i < size; i++) out.write((char*)&(fr[i]), sizeof(Fragment));
    }
    
    void read(std::ifstream &in) {
		in.read((char*)&index, sizeof(int));
		in.read((char*)&stored_scores, sizeof(int));
		in.read((char*)&window, sizeof(int));
        in.read((char*)&pr, sizeof(PrecursorEntry));
		if (stored_scores) in.read((char*)&sc, sizeof(ScoreEntry));
        int size; in.read((char*)&size, sizeof(int));
        fr.resize(size);
        for (int i = 0; i < size; i++) in.read((char*)&(fr[i]), sizeof(Fragment));
    }

	friend inline bool operator < (const QuantEntry &left, const QuantEntry &right) { return left.index < right.index; }
};

class Quant {
public:
	std::vector<QuantEntry> entries;
	double weights[pN];
    
    void write(std::ofstream &out) {
        int size = entries.size();
        out.write((char*)&size, sizeof(int));
        for (int i = 0; i < size; i++) entries[i].write(out);
		out.write((char*)weights, pN * sizeof(double));
    }
    
    void read(std::ifstream &in) {
        int size; in.read((char*)&size, sizeof(int));
        entries.resize(size);
        for (int i = 0; i < size; i++) entries[i].read(in);
		in.read((char*)weights, pN * sizeof(double));
    }
};

class Profile {
public:
	std::vector<QuantEntry> entries;

	Profile(std::vector<std::string> &files) {
		entries.resize(MaxLibSize);
		for (int i = 0; i < files.size(); i++) {
			std::ifstream in(files[i] + std::string(".quant"), std::ifstream::binary);
			Quant Q;
			Q.read(in);
			in.close();

			for (auto it = Q.entries.begin(); it != Q.entries.end(); it++) {
				auto pos = it->pr.index;
				if (pos >= entries.size()) entries.resize(entries.size() + entries.size() / 2);
				if (entries[pos].index < 0) entries[pos] = *it;
				else if (it->pr.qvalue < entries[pos].pr.qvalue) entries[pos].pr = it->pr;
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
	int index, charge;
	double mz, nRT;
	std::vector<Peak> fragments;

	void init(double _mz, double _nRT, int _charge, int _index) {
		mz = _mz;
		nRT = _nRT;
		charge = _charge;
		index = _index;
	}
};

class FastaEntry {
public:
	std::string name;
	mutable std::string proteins;
	FastaEntry(std::string _name, std::string _proteins) { name = _name; proteins = _proteins; }
	friend inline bool operator < (const FastaEntry &left, const FastaEntry &right) { return left.name < right.name; }
};

class Fasta {
public:
	int min_peptide_length = 7, max_peptide_length = 30; // max <= 63
	std::string name;
	std::vector<std::string> peptides;

	std::vector<FastaEntry> proteins;

	bool load(const char * file_name) {
		name = std::string(file_name);
		if (Verbose >= 1) std::cout << "Loading FASTA " << name << "\n";
		std::ifstream in(file_name);
		if (in.fail()) {
			std::cerr << "cannot read the file\n";
			return false;
		}

		std::set<FastaEntry> unique;
		std::string line, sequence, merged, previous, stripped, peptide, modified, protein, next_protein;
		while (getline(in, line)) {
			if (line[0] == '>') {
				auto start = line.find("|");
				if (start == std::string::npos) start = 1;
				auto end = line.find("|", start + 1);
				if (end == std::string::npos) end = line.find(" ", start + 1);
				if (end != std::string::npos) next_protein = line.substr(start + 1, end - start - 1);
				else next_protein = line.substr(start + 1);

				stripped.clear();
				peptide.clear();
				previous.clear();
				bool missed = false;
				for (int i = 0; i < sequence.length(); i++) {
					if (sequence[i] >= 'A' && sequence[i] <= 'Z') {
						stripped.push_back(sequence[i]);
						peptide.push_back(sequence[i]);
						if (UniMod4 && sequence[i] == 'C') peptide += "(UniMod:4)";
						if (sequence[i] == 'R' || sequence[i] == 'K' || i == sequence.length() - 1) {
							add_peptide:
							if (stripped.length() >= min_peptide_length && stripped.length() <= max_peptide_length) {
								auto ins = unique.insert(FastaEntry(stripped, protein));
								if (!ins.second) ins.first->proteins += std::string(";") + protein;
								else {
									peptides.push_back(peptide);
									int in_mod = 0;
									for (int j = 0; j < peptide.size(); j++) {
										if (peptide[j] == '(') in_mod++;
										else if (peptide[j] == ')') in_mod--;
										else if (!in_mod) {
											if (UniMod35 && peptide[j] == 'M') {
												modified = peptide.substr(0, j + 1) + std::string("(UniMod:35)") + peptide.substr(j + 1);
												peptides.push_back(modified);
											}
											if (UniMod21 && (peptide[j] == 'S' || peptide[j] == 'T')) {
												modified = peptide.substr(0, j + 1) + std::string("(UniMod:21)") + peptide.substr(j + 1);
												peptides.push_back(modified);
											}
										}
									}
								}
							}
							if (!missed) {
								merged = previous + peptide;
								previous = peptide;
								if (merged.length() > peptide.length() && MissedCleavage) {
									peptide = merged;
									stripped = get_aas(peptide);
									missed = true;
									goto add_peptide;
								}
							}
							missed = false;
							stripped.clear();
							peptide.clear();
						}
					}
				}
				sequence.clear();
				protein = next_protein;
			} else sequence += line;
		}
		proteins.insert(proteins.begin(), unique.begin(), unique.end());

		in.close();
		return true;
	}
};

class Library {
public:
	std::string name;
	std::vector<Protein> protein_list; // proteins
	std::vector<std::string> precursors; // precursor IDs in canonical format
	float nRT_min = INF, nRT_max = -INF;

	class Entry {
	public:
		Library * lib;
		Peptide target, decoy;
		bool exclude = false, from_fasta = false;
		std::string name; // precursor id
		std::string pep_name; // modified peptide id
		std::string pg_name; // protein group id
		std::vector<Ion> gen; // generated fragments
		int pg_index = 0, best_run = -1, peak;
		float qvalue;

		void init() {
			std::sort(target.fragments.begin(), target.fragments.end(), [](const Peak &left, const Peak &right) { return left.height > right.height; });
			std::sort(decoy.fragments.begin(), decoy.fragments.end(), [](const Peak &left, const Peak &right) { return left.height > right.height; });
			if (target.fragments.size() > MaxF) target.fragments.resize(MaxF);
			if (decoy.fragments.size() > MaxF) decoy.fragments.resize(MaxF);

			// L1 normalisation
			float sum = 0.0;
			for (int i = 0; i < target.fragments.size(); i++) sum += target.fragments[i].height;
			for (int i = 0; i < target.fragments.size(); i++) target.fragments[i].height /= sum;

			sum = 0.0;
			for (int i = 0; i < decoy.fragments.size(); i++) sum += decoy.fragments[i].height;
			for (int i = 0; i < decoy.fragments.size(); i++) decoy.fragments[i].height /= sum;
		}

		inline void generate_decoy() {
			int i;
			decoy = target;

			auto seq = get_sequence(pep_name);

			auto pattern = recognise_fragments(seq, target.fragments);
			auto inv_seq = seq;

			for (i = 2; i < seq.size() - 1; i++) inv_seq[i] = seq[seq.size() - i - 1];
			for (i = 0; i < seq.size(); i++) if (Abs(inv_seq[i] - seq[i]) > 1200.0 * GeneratorAccuracy) break;
			if (i == seq.size()) inv_seq[seq.size() / 2] += 12.0;

			auto dfr = generate_fragments(inv_seq, pattern);
			for (i = 0; i < dfr.size(); i++) decoy.fragments[i].mz = dfr[i].mz;
		}
	};

	std::vector<Entry> entries;

	class Info {
	public:
		Library * lib;
		int n_s; // number of samples (runs)
		int n_entries; // total number of entries
		std::vector<Quant> raw;
		std::map<int, std::vector<std::pair<int, QuantEntry> > > map;

		Info() { }

		void clear() {
			raw.clear();
			map.clear();
		}

		void load(Library * parent, std::vector<std::string> &files) {
			if (Verbose >= 1) std::cout << "Reading quantification information: " << (n_s = files.size()) << " files...\n";

			lib = parent;
			raw.resize(n_s);
			n_entries = 0;
			for (int i = 0; i < n_s; i++) {
				std::ifstream in(std::string(files[i]) + std::string(".quant"), std::ifstream::binary);
				raw[i].read(in);
				for (auto it = raw[i].entries.begin(); it != raw[i].entries.end(); it++) {
					n_entries++;
					int index = it->pr.index;
					auto pos = map.find(index);
					if (pos != map.end()) pos->second.push_back(std::pair<int, QuantEntry>(i, *it));
					else {
						std::vector<std::pair<int, QuantEntry> > v;
						v.push_back(std::pair<int, QuantEntry>(i, *it));
						map.insert(std::pair<int, std::vector<std::pair<int, QuantEntry> > >(index, v));
					}
				}
				in.close();
			}
		}

		void quantify() {
			int i, m, k;

			if (Verbose >= 1) std::cout << "Quantifying peptides...\n";
			for (auto it = map.begin(); it != map.end(); it++) {
				auto v = &(it->second);
				m = (*v)[0].second.fr.size(), k = 0;

				std::vector<double> score(m);
				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) if (jt->second.pr.qvalue <= MaxQuantQvalue) {
					for (i = 0; i < m; i++) score[i] += jt->second.fr[i].corr;
					k++;
				}
				if (k) for (i = 0; i < m; i++) score[i] /= (double)k;

				std::vector<double> ordered(m);
				ordered = score;
				std::sort(ordered.begin(), ordered.end());
				double margin = ordered[Max(0, m - 3)] - E;

				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					for (i = 0, jt->second.pr.quantity = 0.0; i < m; i++)
						if (score[i] >= margin) jt->second.pr.quantity += jt->second.fr[i].quantity[qFiltered];
					jt->second.pr.normalised = jt->second.pr.quantity;
				}
			}

			if (!Normalisation) return;

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
			for (i = 0; i < sums.size(); i++) if (sums[i] <= E) sums[i] = av;
			for (auto it = map.begin(); it != map.end(); it++) {
				auto v = &(it->second);
				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					int index = jt->first;
					jt->second.pr.normalised = (jt->second.pr.quantity * av) / sums[index];
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
					if (jt->second.pr.qvalue < NormalisationQvalue) x[index] = jt->second.pr.normalised, k++, m++;
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
					if (jt->second.pr.qvalue <= NormalisationQvalue) sums[index] += jt->second.pr.normalised;
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
					jt->second.pr.normalised  *= av / sums[index];
				}
			}
		}
	};

	Info info;

	bool load_spectra(const char * file_name, std::map<std::string, Entry> &map) {
		if (Verbose >= 1) std::cout << "Reading library as a collection of spectra\n";
		name = std::string(file_name);
		if (Verbose >= 1) std::cout << "Loading file " << name << "...\n";

		std::ifstream in(file_name, std::ifstream::in);
		if (in.fail()) {
			std::cerr << "cannot read the file\n";
			return false;
		}

		std::string line, pep_name, dat;
		Entry e; e.lib = this;
		std::pair<std::map<std::string, Entry>::iterator,bool> ins;

		while (std::getline(in, line)) {
			if (line[0] >= 'A' && line[0] <= 'Z') {
				int first_space = line.find(" ");
				dat = line.substr(first_space + 1);
				double mz = std::stod(line.substr(line.find(" ")));
				int charge = std::stoi(dat.substr(dat.find(" ")));
				pep_name = line.substr(0, first_space) + to_string(charge);
				ins = map.insert(std::pair<std::string, Entry>(to_canonical(pep_name, charge), e));
				ins.first->second.name = ins.first->first;
				ins.first->second.pg_name = "PG";
				ins.first->second.pep_name = pep_name;
				auto &pep = ins.first->second.target;
				pep.mz = mz;
				pep.nRT = 0.0;
				pep.charge = charge;
			} else if (line[0] >= '1' && line[0] <= '9' && pep_name.length()) {
				double mz = std::stod(line);
				double height = std::stod(line.substr(line.find(" ")));
				auto &pep = ins.first->second.target;
				pep.fragments.push_back(Peak(mz, height));
			}
		}

		in.close();
		InSilicoRTPrediction = false;
		return true;
	}

	bool load(char * file_name) {
		int colInd[libCols];

		name = std::string(file_name);
		if (Verbose >= 1) std::cout << "Loading spectral library " << name << "...\n";

		std::ifstream csv(file_name, std::ifstream::in);
		if (csv.fail()) {
			std::cerr << "cannot read the file\n";
			return false;
		}

		int i, cnt = 0;
		std::map<std::string, Entry> map;
		Entry e; e.lib = this;

		std::string line, word;
		char delim = '\t';
		if (std::string(file_name).find(".csv") != std::string::npos) delim = ','; // not for OpenSWATH libraries
		while (std::getline(csv, line)) {
			std::vector<std::string> words;
			std::stringstream list(line);
			while (std::getline(list, word, delim)) words.push_back(word);
			if (!list && word.empty()) words.push_back("");
			cnt++;

			if (cnt == 1) { // header
				std::vector<std::string>::iterator it, loc;
				for (i = 0; i < libCols; i++) colInd[i] = -1;
				for (it = words.begin(); it != words.end(); it++) {
					for (auto jt = library_headers.begin(); jt != library_headers.end(); jt++)
						if (jt->find(*it) != std::string::npos) {
							colInd[std::distance(library_headers.begin(), jt)] = std::distance(words.begin(), it);
							break;
						}
				}
				for (i = 0; i < libCols; i++) if (colInd[i] < 0 && i < libIsDecoy) {
					if (Verbose >= 1) std::cout << "WARNING: cannot find column " + library_headers[i] << "\n";
					csv.close();
					load_spectra(file_name, map);
					goto finish;
				}
				continue;
			}

			Peak p(std::stof(words[colInd[libFrMz]]), std::stof(words[colInd[libFrI]]));

			bool decoy_fragment = colInd[libIsDecoy] < 0 ? false : std::stoi(words[colInd[libIsDecoy]]);
			if (decoy_fragment) GenDecoys = false;
			int charge = std::stoi(words[colInd[libCharge]]);
			auto ins = map.insert(std::pair<std::string, Entry>(to_canonical(words[colInd[libPr]], charge), e));
			ins.first->second.name = ins.first->first;
			ins.first->second.pg_name = words[colInd[libPG]];
			ins.first->second.pep_name = words[colInd[libPr]];

			auto pep = !decoy_fragment ? (&(ins.first->second.target)) : (&(ins.first->second.decoy));
			pep->mz = std::stof(words[colInd[libPrMz]]);
			pep->nRT = std::stof(words[colInd[libnRT]]);
			pep->charge = charge;
			pep->fragments.push_back(p);
		}

		csv.close();

		finish:
		entries.resize(map.size());
		precursors.resize(map.size());
		i = 0;
		for (auto it = map.begin(); it != map.end(); it++, i++) {
			precursors[i] = it->first;
			entries[i] = it->second;
			entries[i].target.index = entries[i].decoy.index = i;
			if (entries[i].target.nRT < nRT_min) nRT_min = entries[i].target.nRT;
			if (entries[i].target.nRT > nRT_max) nRT_max = entries[i].target.nRT;
		}

		std::set<Protein> pg;
		for (i = 0; i < entries.size(); i++) {
			auto ins = pg.insert(Protein(entries[i].pg_name));
			ins.first->precursors.push_back(i);
			if (exclude_from_training.size()) if (entries[i].pg_name.find(exclude_from_training, 0) != std::string::npos) entries[i].exclude = true;
			if (include_for_training.size()) if (entries[i].pg_name.find(include_for_training, 0) == std::string::npos) entries[i].exclude = true;
		}
		protein_list.insert(protein_list.begin(), pg.begin(), pg.end());

		for (i = 0; i < protein_list.size(); i++)
			for (auto it = protein_list[i].precursors.begin(); it != protein_list[i].precursors.end(); it++)
				entries[*it].pg_index = i;

		if (Verbose >= 1) std::cout << "Spectral library loaded: "
			<< protein_list.size() << " protein groups and "
			<< entries.size() << " precursors.\n";

		return true;
	}

	void load(Fasta &fasta) {
		int i = entries.size(), cnt;
		if (!i && !InSilicoRTPrediction) nRT_min = -1.0, nRT_max = 1.0;
		double default_nRT = (nRT_min + nRT_max) * 0.5;
		Entry e;
		e.lib = this; e.from_fasta = true;
		std::vector<double> scores_y;
		for (auto &pep : fasta.peptides) {
			auto seq = get_sequence(pep);
			auto fragments = generate_fragments(seq, 1, loss_none, &cnt);
			auto aas = get_aas(pep);
			double mass = sum(seq) + proton + OH;
			int min_charge = Max(1, (int)(mass / MaxPrMz));
			int max_charge = Max(min_charge, Min(4, (int)(mass / MinPrMz)));
			double nRT = InSilicoRTPrediction ? predict_nrt(pep) : default_nRT;
			if (nRT < nRT_min) nRT_min = nRT;
			if (nRT > nRT_max) nRT_max = nRT;
			for (int charge = min_charge; charge <= max_charge; charge++) {
				y_scores(scores_y, charge, aas);
				to_exp(scores_y);
				e.target.init(proton + mass / (double)charge, nRT, charge, i);
				e.target.fragments.clear();
				for (auto &fr : fragments)
					if (fr.mz >= MinFrMz && fr.mz <= MaxFrMz)
						if (fr.type == type_y) e.target.fragments.push_back(Peak(fr.mz, scores_y[fr.index]));
				if (e.target.fragments.size() >= 4) {
					e.pep_name = pep;
					e.name = to_canonical(pep, charge);
					auto pos = std::lower_bound(precursors.begin(), precursors.end(), e.name);
					if (pos == precursors.end() || *pos != e.name)
						entries.push_back(e), i++;
				}
			}
		}

		std::set<Protein> pg;
		std::string empty_string;
		for (i = 0; i < entries.size(); i++) {
			if (entries[i].from_fasta) {
				auto stripped = get_aas(entries[i].pep_name);
				auto pos = std::lower_bound(fasta.proteins.begin(), fasta.proteins.end(), FastaEntry(stripped, empty_string));
				entries[i].pg_name = pos->proteins;
			}
			auto ins = pg.insert(Protein(entries[i].pg_name));
			ins.first->precursors.push_back(i);
		}
		protein_list.clear();
		protein_list.insert(protein_list.begin(), pg.begin(), pg.end());

		for (i = 0; i < protein_list.size(); i++)
			for (auto it = protein_list[i].precursors.begin(); it != protein_list[i].precursors.end(); it++)
				entries[*it].pg_index = i;

		if (Verbose >= 1) std::cout << entries.size() << " precursors generated\n";
	}

	void save(std::string &file_name, std::vector<int> * ref = NULL, bool decoys = false) {
		if (Verbose >= 1) std::cout << "Saving spectral library to " << file_name << "...\n";
		std::ofstream out(file_name, std::ofstream::out);
		out << "FileName\tPrecursorMz\tProductMz\tTr_recalibrated\ttransition_name\tLibraryIntensity\ttransition_group_id\tdecoy\tPeptideSequence\tProteinName\tFullUniModPeptideName\tModifiedPeptide\t";
		out << "PrecursorCharge\tPeptideGroupLabel\tUniprotID\tFragmentType\tFragmentCharge\tFragmentSeriesNumber\tFragmentLossType\n";
		out.precision(7);

		int cnt = -1;
		for (auto &it : entries) {
			cnt++;
			if (it.best_run < 0 && it.from_fasta) continue;
			if (ref != NULL) {
				auto pos = std::lower_bound(ref->begin(), ref->end(), cnt);
				if (pos == ref->end()) continue;
				if (*pos != cnt) continue;
			}
			auto seq = get_aas(it.pep_name);
			auto masses = get_sequence(it.pep_name);
			auto &ions = it.from_fasta ? it.gen : recognise_fragments(masses, it.target.fragments);
			auto &pep = it.target;
			auto pep_name = to_canonical(it.pep_name);
			auto name = pep_name + std::to_string(pep.charge);
			for (int dc = 0; dc <= 1; dc++) {
				std::string prefix = dc ? std::string("DECOY") : "";
				if (!decoys && dc) continue;
				auto &pep = dc ? it.decoy : it.target;
				for (int fr = 0; fr < pep.fragments.size(); fr++) {
					int fr_type = ions[fr].type;
					int fr_loss = ions[fr].loss;
					int fr_num = ions[fr].index;
					int fr_charge = ions[fr].charge;
					if (!it.from_fasta) {
						bool skip = false;
						for (auto pos = ions.begin(); pos < ions.begin() + fr; pos++)
							if (pos->index == fr_num && pos->type == fr_type && pos->charge == fr_charge && pos->loss == fr_loss) {
								skip = true;
								break;
							}
						if (skip) continue;
					}
					if (pep.fragments[fr].height < E) continue;
					auto fr_name = name + std::string("_") + std::to_string(char_from_type[fr_type])
						+ std::string("_") + std::to_string(fr_charge) + std::string("_")
						+ std::to_string(fr_loss) + std::string("_") + std::to_string(fr_num);
					out << (it.from_fasta ? ms_files[it.best_run] : lib_file) << "\t"
						<< pep.mz << "\t"
						<< pep.fragments[fr].mz << "\t"
						<< pep.nRT << "\t"
						<< prefix + fr_name << "\t"
						<< pep.fragments[fr].height << "\t"
						<< prefix + name << "\t"
						<< dc << "\t"
						<< seq << "\t"
						<< prefix + it.pg_name << "\t"
						<< pep_name << "\t"
						<< pep_name << "\t"
						<< pep.charge << "\t"
						<< prefix + pep_name << "\t"
						<< it.pg_name << "\t"
						<< char_from_type[fr_type] << "\t"
						<< fr_charge << "\t"
						<< (fr_type == type_b ? fr_num : masses.size() - fr_num) << "\t"
						<< name_from_loss[fr_loss] << "\n";
				}
			}
		}

		out.close();
	}

	void generate_decoys() {
		if (Verbose >= 1) std::cout << "Generating decoys...\n";
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (GenDecoys) it->generate_decoy();
			it->init();
		}
		if (ExportLibrary) save(out_file.substr(0, out_file.find_last_of(".")) + std::string(".lib.tsv"), NULL, true);
	}

	void fragment() {
		assert(FastaSearch);
		assert(MaxF > 1000);
		int i, cnt;
		for (auto &e : entries) {
			if (!e.from_fasta) continue; // do not overwrite existing fragmentation patterns
			auto pep = &(e.target);

			auto seq = get_sequence(e.pep_name);
			std::vector<Ion> gen;
			for (int charge = 1; charge < Max(2, pep->charge); charge++) {
				for (int loss = loss_none; loss < loss_CO; loss++) {
					auto frs = generate_fragments(seq, charge, loss, &cnt);
					gen.insert(gen.end(), frs.begin(), frs.end());
				}
			}
			for (i = cnt = 0; i < gen.size(); i++) if (gen[i].mz >= MinOutFrMz && gen[i].mz <= MaxOutFrMz) cnt++;
			pep->fragments.resize(cnt);
			e.gen.resize(cnt);
			for (i = cnt = 0; i < gen.size(); i++) if (gen[i].mz >= MinOutFrMz && gen[i].mz <= MaxOutFrMz) {
				pep->fragments[cnt].mz = gen[i].mz;
				pep->fragments[cnt].height = 1.0;
				e.gen[cnt] = gen[i];
				cnt++;
			}
		}
	}

	void quantify_proteins(int N, double q_cutoff) { // top N method for protein quantification
		int i, j, pass, n = protein_list.size();

		if (Verbose >= 1) std::cout << "Quantifying proteins...\n";
		std::vector<float> max_q(n * info.n_s, 1.0), level(n * info.n_s * N), quant(n * info.n_s);

		for (pass = 0; pass <= 3; pass++) {
			for (auto it = info.map.begin(); it != info.map.end(); it++) {
				auto v = &(it->second);

				for (auto jt = v->begin(); jt != v->end(); jt++) {
					int s = jt->first;
					auto pr = &(jt->second.pr);
					int pg = entries[pr->index].pg_index;

					float quantity = pr->quantity;
					float q = pr->qvalue;

					int index = pg * info.n_s + s;
					auto pos_q = &(max_q[index]);
					auto pos_l = &(level[index * N]);

					if (pass == 0) {
						if (q < *pos_q && *pos_q > q_cutoff) *pos_q = Max(q_cutoff, q);
					}
					else if (pass == 1 && q <= *pos_q) {
						for (i = 0; i < N; i++) if (quantity > pos_l[i]) {
							for (j = N - 1; j > i; j--) pos_l[j] = pos_l[j - 1];
							pos_l[i] = quantity;
							break;
						}
					}
					else if (pass == 2 && q <= *pos_q && quantity >= pos_l[N - 1]) quant[index] += quantity;
					else if (pass == 3) pr->pg_quantity = quant[index];
				}
			}
		}
	}

	void report(std::string &file_name) {
		if (Verbose >= 1) std::cout << "Writing report... ";

		std::ofstream out(file_name, std::ofstream::out);
		out << oh[outFile] << "\t" << oh[outPG] << "\t" << oh[outPGQ] << "\t" << oh[outModSeq] << "\t"
			<< oh[outPrId] << "\t" << oh[outCharge] << "\t" << oh[outQv] << "\t" << oh[outPrQ] << "\t"
			<< oh[outPrQRaw] << "\t" << oh[outRT] << "\t" << oh[outnRT];
		if (ExtraReportInfo) out << "\tPredicted.nRT";
		out << "\n";

		for (auto it = info.map.begin(); it != info.map.end(); it++) {
			auto v = &(it->second);
			auto entry = &(entries[it->first]);

			for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
				out << ms_files[jt->first] << "\t"
					<< entry->pg_name << "\t"
					<< jt->second.pr.pg_quantity << "\t"
					<< entry->pep_name << "\t"
					<< entry->name << "\t"
					<< entry->target.charge << "\t"
					<< jt->second.pr.qvalue << "\t"
					<< jt->second.pr.normalised << "\t"
					<< jt->second.pr.quantity << "\t"
					<< jt->second.pr.RT << "\t"
					<< jt->second.pr.nRT;
				if (ExtraReportInfo) out << "\t" << jt->second.pr.predicted_nRT;
				out << "\n";
			}
		}

		out.close();
		if (Verbose >= 1) std::cout << "Report saved to " << file_name << ".\n";
	}
};

void learn_from_library(std::string &file_name) {
	int i, j, n, P_y = y_delta_size(), M = in_silico_rt_size(), pos_y = 0, pos_RT = 0;
	Library lib;
	if (!lib.load(&(file_name[0]))) Throw("Cannot load the library");
	if (Verbose >= 1) std::cout << "Learning peptide characteristics...\n";

	yDelta.resize(P_y);
	InSilicoRT.resize(M);
	std::vector<int> y_index;

	std::vector<double> b_y, b_b, b_RT;
	typedef Eigen::Triplet<double> T;
	std::vector<T> TL_y, TL_b, TL_RT;
	std::vector<double> count(M);
	for (auto &e : lib.entries) {
		auto aas = get_aas(e.pep_name);
		auto seq = get_sequence(e.pep_name);
		auto frs = recognise_fragments(seq, e.target.fragments, true, 1, 1);

		if (InSilicoRTPrediction) {
			count[0] = 1;
			for (i = 1; i < M; i++) count[i] = 0;
			count_rt_aas(count, e.pep_name);
			for (i = 0; i < M; i++) if (Abs(count[i]) > E) TL_RT.push_back(T(pos_RT, i, (double)count[i]));
			b_RT.push_back(e.target.nRT);
			pos_RT++;
		}

		y_index.resize(n = seq.size());
		for (i = 0; i < y_index.size(); i++) y_index[i] = -1;
		for (i = 0; i < frs.size(); i++) if (frs[i].charge) {
			auto &fr = frs[i];
			if (fr.charge == 1 && fr.type == type_y && fr.loss == loss_none) y_index[fr.index] = i;
		}
		for (i = 2; i < y_index.size(); i++) if (y_index[i] >= 0 && y_index[i - 1] >= 0) {
			double ratio = log(e.target.fragments[y_index[i]].height / e.target.fragments[y_index[i - 1]].height);
			if (Abs(ratio) >= 5.0) continue;
			b_y.push_back(ratio);

			for (j = 0; j < aas.length(); j++) TL_y.push_back(T(pos_y, y_delta_composition_index(aas[j]), 1.0));
			TL_y.push_back(T(pos_y, y_delta_charge_index(), (double)e.target.charge));
			for (j = Max(i - yDeltaS, 0); j <= Min(i + yDeltaS, n - 1); j++)
				TL_y.push_back(T(pos_y, y_delta_index(aas[j], j - i), 1.0));
			if (i <= yNTermD) TL_y.push_back(T(pos_y, y_nterm_index(i - 1), 1.0));
			if (i >= n - yCTermD) TL_y.push_back(T(pos_y, y_cterm_index(aas[n - 1], n - 1 - i), 1.0));
			TL_y.push_back(T(pos_y, 0, 1.0));

			pos_y++;
		}
	}

	{
		auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b_y.data(), b_y.size());
		Eigen::SparseMatrix<double, Eigen::RowMajor> A(pos_y, P_y);
		A.setFromTriplets(TL_y.begin(), TL_y.end());
		Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
		SQR.compute(A);
		auto X = SQR.solve(B);
		for (i = 0; i < P_y; i++) yDelta[i] = X[i];
	}

	if (InSilicoRTPrediction) {
		auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b_RT.data(), b_RT.size());
		Eigen::SparseMatrix<double, Eigen::RowMajor> A(pos_RT, M);
		A.setFromTriplets(TL_RT.begin(), TL_RT.end());
		Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
		SQR.compute(A);
		auto X = SQR.solve(B);
		for (i = 0; i < M; i++) InSilicoRT[i] = X[i];
	}

	// Pearson correlation between predicted and actual spectra
	double r_y = 0.0;
	int cnt = 0;
	std::vector<double> actual_y, predicted_y, rt_d;
	for (auto &e : lib.entries) {
		auto aas = get_aas(e.pep_name);
		auto seq = get_sequence(e.pep_name);
		auto frs = recognise_fragments(seq, e.target.fragments, true, 1, 1);

		if (InSilicoRTPrediction) {
			double nrt = predict_nrt(e.pep_name);
			rt_d.push_back(Abs(nrt - e.target.nRT));
		}

		predicted_y.resize(n = seq.size());
		actual_y.resize(n);
		for (int i = 0; i < n; i++) actual_y[i] = 0.0;

		for (int i = 0; i < frs.size(); i++) if (frs[i].charge) {
			auto &fr = frs[i];
			if (fr.charge == 1 && fr.loss == loss_none)
				if (fr.type == type_y) actual_y[fr.index] = e.target.fragments[i].height;
		}
		y_scores(predicted_y, e.target.charge, aas);
		to_exp(predicted_y);
		r_y += corr(&(predicted_y[1]), &(actual_y[1]), seq.size() - 1);
		cnt++;
	}
	std::cout << "y-series fragmentation prediction: Pearson correlation = " << r_y / (double)cnt << "\n";
	if (InSilicoRTPrediction) {
		std::sort(rt_d.begin(), rt_d.end());
		std::cout << "nRT prediction: median error = " << rt_d[rt_d.size() / 2] << "\n";
	}
}

void train(int thread_id) {
	for (int i = 0; i * Threads + thread_id < nnBagging; i++)
		net[i * Threads + thread_id].optimise();
}

void global_training(Library * lib) {
	int i, nnCnt = 0;

	lib->info.load(lib, ms_files);
	if (Verbose >= 1) std::cout << "Training the neural network...\n";

	net.resize(nnBagging);
	netLock.resize(nnBagging);
	training.resize(nnBagging);
	training_classes.resize(nnBagging);
	for (i = 0; i < nnBagging; i++) {
		if (net[i].network != NULL) destroyNetwork(net[i].network);
		training[i].resize(2 * lib->info.n_entries);
		training_classes[i].resize(2 * lib->info.n_entries);
	}

	for (auto it = lib->info.map.begin(); it != lib->info.map.end(); it++) {
		auto v = &(it->second);
		for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
			for (i = 0; i < nnBagging; i++) training[i][nnCnt] = jt->second.sc.target_scores, training_classes[i][nnCnt] = target_nn;
			nnCnt++;
			if (jt->second.pr.decoy_found) {
				for (i = 0; i < nnBagging; i++) training[i][nnCnt] = jt->second.sc.decoy_scores, training_classes[i][nnCnt] = decoy_nn;
				nnCnt++;
			}
		}
	}

	size_t* hiddenSize = (size_t*)alloca(nnGlobalHidden * sizeof(size_t));
	Activation* hiddenActivation = (Activation*)alloca(nnGlobalHidden * sizeof(Activation));
	for (i = 0; i < nnGlobalHidden; i++) {
		hiddenSize[i] = (nnGlobalHidden - i) * 5;
		hiddenActivation[i] = tanH;
	}

	for (i = 0; i < nnBagging; i++) {
		DataSet* trainingData = createDataSet(nnCnt, pN, &(training[i][0]));
		DataSet* trainingClasses = createDataSet(nnCnt, 2, &(training_classes[i][0]));

		net[i].network = createNetwork(pN, nnGlobalHidden, hiddenSize, hiddenActivation, 2, softmax, net[i].random);
		net[i].lossFunction = CROSS_ENTROPY_LOSS;
		net[i].batchSize = Min(50, Max(1, nnCnt / 100));
		net[i].learningRate = nnLearning;
		net[i].searchTime = 0;
		net[i].regularizationStrength = Regularisation;
		net[i].momentumFactor = 0.9;
		net[i].shuffle = 1;
		net[i].verbose = (Verbose >= 5); // large amount of RAM required
		net[i].data = trainingData;
		net[i].classes = trainingClasses;
		net[i].maxIters = (nnCnt * nnGlobalEpochs) / net[i].batchSize;
	}

	std::vector<std::thread> thr;
	for (i = 0; i < Threads; i++) thr.push_back(std::thread(train, i));
	for (i = 0; i < Threads; i++) thr[i].join();

	lib->info.clear();
}

std::vector<int> rt_stats, rt_ref;
std::vector<float> rt_coo;
std::vector<std::pair<float, float> > rt_data;
std::vector<float> rt_delta, mass_acc;

class Run {
public:
	int run_index;
    std::string name; // run name
    float weights[pN], best_weights[pN];
    std::vector<Tandem> scans;
    std::vector<Tandem> ms1;
    std::vector<float> ms1_RT;

	Library * lib;
    
    int n_scans;
    std::vector<float> predicted_nRT; // nRT from RT prediction 
    std::vector<float> scan_RT; // list of RTs of scans
	double nRT_window;

	std::vector<std::vector<float> > MS1, MS2, MS2_min;
	std::vector<std::vector<int> > PeakList, BestFrList;
	std::vector<std::vector<float> > CorrSumList;

    int curr_iter, curr_batch;
    double min_target_decoy_score = 0.0, nRT_cscore, nRT_ref_score, PeakWidth = 0.0, MassAccuracy = GlobalMassAccuracy, MassAccuracyMs1 = GlobalMassAccuracyMs1;
	std::vector<double> MassCorrection, MassCorrectionMs1, MassCalSplit, MassCalSplitMs1, MassCalCenter, MassCalCenterMs1;
    
    std::vector<Parameter> pars;
    bool par_seek[pN], par_learn[pN];
	int max_ids = 0;

	bool standardised = false;
	bool use_nn = true;
	bool in_ref_run = false, nRT_windowed_search = false, mz_calibrated = false, acc_calibrated = false, acc_ms1_calibrated = false, window_calculated = false;

    Run(int _run_index) {
		run_index = _run_index;
		pars.push_back(Parameter(0, 0)); // pTimeCorr
		pars.push_back(Parameter(1, 0)); // pMinCorr
        pars.push_back(Parameter(1, 0)); // pCos
        pars.push_back(Parameter(1, 0)); // pMs1TimeCorr
		pars.push_back(Parameter(1, 0)); // pRT
		pars.push_back(Parameter(1, 0)); // pnRT
		pars.push_back(Parameter(iN, iN)); // pCharge
		for (int i = 0; i < nnF; i++) pars.push_back(Parameter(iN, iN)); // pRef
		for (int i = 0; i < nnF; i++) pars.push_back(i < nnF - 1 ? Parameter(1, 0) : Parameter(iN, iN)); // pSig
		for (int i = 0; i < nnF; i++) pars.push_back(Parameter(iN, iN)); // pCorr
		for (int i = 0; i < nnF * nnW; i++) pars.push_back(Parameter(iN, iN)); // pArray
        
        for (int i = 0; i < pN; i++) weights[i] = 0.0;
        weights[0] = 1.0;

		if (!GlobalNN) for (int i = 0; i < nnBagging; i++) {
			if (net[i].network != NULL) destroyNetwork(net[i].network);
			net[i].network = NULL;
		}
		training.resize(nnBagging);
		training_classes.resize(nnBagging);

		MS1.resize(Threads), MS2.resize(Threads), MS2_min.resize(Threads);
		PeakList.resize(Threads), BestFrList.resize(Threads), CorrSumList.resize(Threads);

		if (Calibrate) MassAccuracy = MassAccuracyMs1 = CalibrationMassAccuracy;
		MassCorrection.resize(1 + MassCalBins * 2, 0.0); MassCorrectionMs1.resize(1 + MassCalBinsMs1 * 2, 0.0);
		MassCalSplit.resize(MassCalBins + 1, 0.0); MassCalSplitMs1.resize(MassCalBinsMs1 + 1, 0.0);
		MassCalCenter.resize(MassCalBins + 1, 0.0); MassCalCenterMs1.resize(MassCalBinsMs1 + 1, 0.0);
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
    
    template <bool seek> inline double combined_score(float * scores) {
		if (curr_iter < nnIter + (int)seek || GlobalTraining || !use_nn) {
			if (seek) return scalar(scores, weights, par_seek, pN);
			else return scalar(scores, weights, par_learn, pN);
		} else {
			Matrix* mat = createMatrix(1, pN, scores);
			double r = 0.0;
			for (int i = 0; i < nnBagging; i++) {
				if (seek) while (!netLock[i].set()) {}
				forwardPass(net[i].network, mat);
				Matrix* result = getOuput(net[i].network);
				if (seek) netLock[i].free();
				r += result->data[0];
			}
			return r;
		}
    }
    
    class Precursor {
    public:
		Lock lock;

        Run * run;
        bool found = false, decoy = false, initialised = false, training = false, test = false, exclude = false;
        double mz, nRT;
        int index, charge, length;
		int S, W; // scan window "radius" and size
        std::vector<Peak> fragments;
		std::vector<float> Mz, Ref;

		int thread_id = 0;
        int scan_number, peak_number;
        std::vector<int> scan_index;
        std::vector<Elution> peaks; 
		std::vector<float> scoring;
        
        // information reported
        std::vector<Fragment> quant; // quantification information for individual fragments
        double RT, predicted_nRT, combined_score, quantity, qvalue;
		float scores[pN];
        int apex, peak_width, best_peak;
		int peak_pos, best_fragment;
		float mass_delta = 0.0, mass_delta_mz = 0.0, mass_delta_ms1 = 0.0;
        
		void init(Run * _run, Peptide &p) {
			run = _run;
			index = p.index;
			charge = p.charge;
			mz = p.mz;
			nRT = p.nRT;
			fragments = p.fragments;
		}

		void build_index() {
			int i, k, n = run->scans.size();

			for (i = k = 0; i < n; i++) if (run->scans[i].has(mz)) k++;
			scan_index.resize(scan_number = k);
			for (i = k = 0; i < n; i++) if (run->scans[i].has(mz)) scan_index[k++] = i;

			if (!QuantOnly) {
				if (WindowRadius) S = WindowRadius;
				else if (run->PeakWidth > E) S = Max(1, int(ScanScale * run->PeakWidth));
				else S = Max(1, k / ScanFactor);
			}
			W = 2 * S + 1;
		}

		class Searcher {
		public:
			Run * run;
			Precursor * pr;

			int n, m;

			float Mz, *ms1, *ms2, *ms2_e, *ms2_min;

			float *mz, *ref;

			Searcher(Precursor * precursor) {
				pr = precursor;
				run = pr->run;

				n = run->scans.size();
				m = pr->fragments.size();

				Mz = pr->mz;

				pr->Mz.resize(m);
				pr->Ref.resize(m);
				mz = &(pr->Mz[0]);
				ref = &(pr->Ref[0]);
				for (int fr = 0; fr < m; fr++) mz[fr] = pr->fragments[fr].mz, ref[fr] = pr->fragments[fr].height;
			}

			template <bool get_mz> void chromatogram(int from, int to, bool get_ms1, bool quant) { // [from, to)
				int i, k, fr, l = to - from, pos, ind, center = (from + to) / 2;
				float ms1_left = 0.0, ms1_right = 0.0, ms1_left_mz = 0.0, ms1_right_mz = 0.0, rt, ms1_left_RT = 0.0, ms1_right_RT = 0.0, peak_mz = 0.0;

				run->MS2[pr->thread_id].resize(l * m);
				run->MS2_min[pr->thread_id].resize(l * m);
				if (get_ms1) run->MS1[pr->thread_id].resize(l);
				if (!run->MS1[pr->thread_id].size()) get_ms1 = false;
				ms1 = &(run->MS1[pr->thread_id][0]);
				ms2 = &(run->MS2[pr->thread_id][0]);
				ms2_min = &(run->MS2_min[pr->thread_id][0]);
				for (i = 0; i < l * m; i++) ms2[i] = ms2_min[i] = 0.0;
				if (get_ms1) for (i = 0; i < l; i++) ms1[i] = 0.0;

				for (k = from; k < to; k++) {
					i = pr->scan_index[k];
					ind = k - from;

					if (!quant && run->nRT_windowed_search) { // main search phase only
						double margin = run->nRT_window;
						if (run->predicted_nRT[pr->scan_index[Min(to - 1, k + pr->W)]] < pr->nRT - margin) continue;
						if (run->predicted_nRT[pr->scan_index[Max(from, k - pr->W)]] > pr->nRT + margin) break;
					}

					if (get_ms1) {
						rt = run->scan_RT[i];
						if (rt > ms1_right_RT) {
							auto ms1_ptr = std::lower_bound(run->ms1_RT.begin(), run->ms1_RT.end(), rt);
							for (pos = std::distance(run->ms1_RT.begin(), ms1_ptr); pos < run->ms1_RT.size()
								&& (run->ms1[pos].window_low >= Mz || run->ms1[pos].window_high <= Mz); pos++);
							if (pos == run->ms1_RT.size()) ms1_right_RT = INF, ms1[k] = ms1_right;
							else {
								ms1_left = ms1_right, ms1_left_RT = ms1_right_RT, ms1_left_mz = ms1_right_mz;
								ms1_right_RT = *ms1_ptr;
								double query_mz = run->predicted_mz_ms1(&(run->MassCorrectionMs1[0]), Mz, run->scan_RT[i]);
								ms1_right = run->ms1[pos].level<get_mz>(query_mz, run->MassAccuracyMs1, &ms1_right_mz);
							}
						}
						ms1[ind] = ((rt - ms1_left_RT) * ms1_right + (ms1_right_RT - rt) * ms1_left) / Max(E, ms1_right_RT - ms1_left_RT);
					}
					for (fr = 0; fr < m; fr++) {
						double query_mz = run->predicted_mz(&(run->MassCorrection[0]), mz[fr], run->scan_RT[i]);
						ms2[l * fr + ind] = run->scans[i].level<get_mz>(query_mz, run->MassAccuracy, &(peak_mz));
						if (ms2[l * fr + ind] < MinPeakHeight) ms2[l * fr + ind] = 0.0;
						if (get_mz && k == center && fr == pr->best_fragment) {
							pr->mass_delta = ms2[l * fr + ind] > E ? (peak_mz - mz[fr]) : 0.0;
							pr->mass_delta_mz = mz[fr];
							rt = run->scan_RT[i];
							pr->mass_delta_ms1 = (pr->scores[pMs1TimeCorr] >= MassCalMs1Corr && ms1_left_mz > E && ms1_right_mz > E) ?
								(((rt - ms1_left_RT) * ms1_right_mz + (ms1_right_RT - rt) * ms1_left_mz) / Max(E, ms1_right_RT - ms1_left_RT) - Mz) : 0.0;
						}
					}
				}

				for (fr = 0; fr < m; fr++)
					for (k = 1; k < l - 1; k++)
						ms2_min[l * fr + k] = Min(Min(ms2[l * fr + k - 1], ms2[l * fr + k + 1]), ms2[l * fr + k]);
			}

			void peaks() {
				int i, k, fr, pos, next, best_fr = -1, l = pr->scan_number, n_peaks = 0;
				double max, best = -INF, s;
				double * corr_matrix = (double*)alloca(m * m * sizeof(double));
				float * elution = (float*)alloca(pr->W * sizeof(float));

				run->PeakList[pr->thread_id].resize(l);
				run->BestFrList[pr->thread_id].resize(l);
				run->CorrSumList[pr->thread_id].resize(l);
				int * peak_list = &(run->PeakList[pr->thread_id][0]);
				int * best_fr_list = &(run->BestFrList[pr->thread_id][0]);
				float * corr_sum_list = &(run->CorrSumList[pr->thread_id][0]);

				for (k = pr->S, pr->peak_number = 0; k < l - pr->S; k++) {
					for (fr = 0; fr < m; fr++) if (ms2[fr * l + k]) break;
					if (fr >= m) continue; // no signal

					// find the best fragment in the window
					for (fr = 0, max = 0.0; fr < m; fr++) {
						s = sum(&(ms2[fr * l + k - pr->S]), pr->W);
						if (s > max) max = s, best_fr = fr;
						for (next = fr + 1; next < m; next++)
							corr_matrix[fr * m + next] = corr_matrix[next * m + fr]
								= corr(&(ms2[fr * l + k - pr->S]), &(ms2[next * l + k - pr->S]), pr->W);
					}
					for (fr = 0, max = 0.0; fr < m; fr++) {
						s = 0.0;
						for (next = 0; next < m; next++) if (next != fr) s += corr_matrix[fr * m + next];
						if (s > max) max = s, best_fr = fr;
					}
					if (max < MinCorrScore) continue;

					// smooth the curve
					smooth(elution, &(ms2[best_fr * l + k - pr->S]), pr->W);

					// check if k is at the apex
					double curr_evidence = Min(elution[pr->S], Max(elution[pr->S - 1], elution[pr->S + 1])); // suppress random noise this way
					double max_evidence = curr_evidence;
					for (pos = k - pr->S + 1; pos <= k + pr->S - 1; pos++) {
						int ind = pos - k + pr->S;
						double evidence = Min(elution[ind], Max(elution[ind - 1], elution[ind + 1]));
						if (evidence > max_evidence) max_evidence = evidence;
					}
					if (curr_evidence < max_evidence * PeakApexEvidence) continue;

					peak_list[n_peaks] = k, best_fr_list[n_peaks] = best_fr, corr_sum_list[n_peaks++] = max;
					if (max > best) best = max;
				}

				for (i = pos = 0; i < n_peaks; i++) if (corr_sum_list[i] >= best - MaxCorrDiff) pos++;
				pr->peaks.resize(pr->peak_number = pos);
				pr->scoring.resize(pos * pN, 0.0);
				for (i = pos = 0; i < n_peaks; i++) if (corr_sum_list[i] >= best - MaxCorrDiff) {
					pr->peaks[pos].peak = peak_list[i];
					pr->peaks[pos].apex = pr->scan_index[peak_list[i]];
					pr->peaks[pos].best_fragment = best_fr_list[i];
					pos++;
				}
			}

			void score(int peak) {
				assert(peak >= 0);
				assert(peak < pr->peaks.size());

				int k = pr->peaks[peak].peak, pos, fr, l = pr->scan_number, ind, i;
				float *x = (float*)alloca(m * sizeof(float)), *elution = (float*)alloca(pr->W * sizeof(float)), *sc = &(pr->scoring[peak * pN]);
				double w, weight;

				assert(k >= 0);
				assert(k < l - pr->S);

				// elution curve
				smooth(elution, &(ms2[pr->peaks[peak].best_fragment * l + k - pr->S]), pr->W);
				for (i = 0; i < pN; i++) sc[i] = 0.0;

				// pRef
				for (fr = 0; fr < Min(m, nnF); fr++)
					sc[pRef + fr] = ref[fr];

				// pCharge
				sc[pCharge] = ((double)pr->charge - 2.0);

				// pCos
				for (pos = k - pr->S, w = 0.0; pos <= k + pr->S; pos++) {
					ind = pos - k + pr->S;
					weight = Sqr(elution[ind]);
					if (weight < E) continue;
					w += weight;

					for (fr = 0; fr < m; fr++) x[fr] = ms2[fr * l + pos];
					sc[pCos] += cos(&(x[0]), &(ref[0]), m) * weight;
				}
				assert(w > E);
				sc[pCos] /= w;

				// time corr
				for (fr = 0; fr < m; fr++) {
					double u = corr(elution, &(ms2[fr * l + k - pr->S]), pr->W);
					if (fr < nnF) sc[pCorr + fr] += u;
					sc[pTimeCorr] += u;
					sc[pMinCorr] += corr(elution, &(ms2_min[fr * l + k - pr->S]), pr->W);
				}

				sc[pMs1TimeCorr] = corr(elution, &(ms1[k - pr->S]), pr->W);

				// raw data
				for (fr = 0; fr < Min(m, nnF); fr++) {
					int nnD = 2 * (pr->S / (nnS * 2)) + 1;
					for (pos = k - pr->S; pos <= k + pr->S; pos++) {
						int dist = Min(nnS, (Abs(pos - k) + nnD / 2) / nnD);
						ind = pos >= k ? (nnS + dist) : (nnS - dist);
						sc[pArray + nnW * fr + ind] += ms2[fr * l + pos];
						sc[pSig + fr] += ms2[fr * l + pos];
					}
				}
				w = sum(sc + pArray, nnF * nnW);
				if (w > E) {
					for (fr = 0; fr < nnF; fr++) {
						if (sc[pSig + fr] <= E) continue;
						for (i = 0; i < nnW; i++) sc[pArray + nnW * fr + i] /= sc[pSig + fr];
						sc[pSig + fr] /= w;
					}
				}
			}
		};
        
		void score_RT(int peak) {
			assert(peak >= 0);
			assert(peak < peaks.size());
			assert(peaks[peak].apex >= 0);
			assert(peaks[peak].apex < run->predicted_nRT.size());

			double span = Max(E, run->lib->nRT_max - run->lib->nRT_min);
			double pos = Min(run->predicted_nRT[peaks[peak].apex] - run->lib->nRT_min, span);
			double delta = run->predicted_nRT[peaks[peak].apex] - nRT;
			float *sc = &(scoring[peak * pN]);
			sc[pRT] = pos / span;
			sc[pnRT] = sqrt(Min(Abs(delta), span) / span);
		}

		void score_RT() {
			double span = Max(E, run->lib->nRT_max - run->lib->nRT_min);
			double pos = Min(run->predicted_nRT[apex] - run->lib->nRT_min, span);
			double delta = run->predicted_nRT[apex] - nRT;
			float *sc = &(scores[0]);
			sc[pRT] = pos / span;
			sc[pnRT] = sqrt(Min(Abs(delta), span) / span);
		}

		void search() {
			if (!initialised) {
				Searcher searcher(this);
				build_index();
				searcher.chromatogram<true>(0, scan_number, true, false);
				searcher.peaks();
				for (int peak = 0; peak < peak_number; peak++) searcher.score(peak);
				initialised = true;
			}
			for (int peak = 0; peak < peak_number; peak++) score_RT(peak);
		}

		void free() {
			initialised = false;
			std::vector<int>().swap(scan_index);
			std::vector<Elution>().swap(peaks);
			std::vector<float>().swap(scoring);
		}
        
        void quantify(int peak) {
			assert(peak >= 0);
			if (!initialised) { // if we are doing quantification only (identification information read from .quant files)
				if (!found) return;
				build_index();
			}

			int fr, pos, k = peak, best_fr = best_fragment;
            double w, r, e;
			Searcher searcher(this);

			int low = Max(0, k - Min(MaxDpPP / 2 - (1 ^ (MaxDpPP & 1)), Max(1, S / 2)));
			int high = Min(scan_number - 1, k + Min(MaxDpPP / 2, Max(1, S / 2)));
			int len = high - low + 1;
			float *elution = (float*)alloca(len * sizeof(float));
			float *signal = (float*)alloca(len * sizeof(float));

			searcher.chromatogram<true>(low, high + 1, true, true);
			for (pos = 0; pos < len; pos++) signal[pos] = searcher.ms2[best_fr * len + pos];
			smooth(&(elution[0]), &(signal[0]), len);

			for (pos = 0, w = 0.0; pos < len; pos++) w += Sqr(elution[pos]);
			assert(w > E);

			e = elution[k - low] * 0.5;
			for (pos = peak_width = 0; pos < len; pos++) if (elution[pos] >= e) 
				peak_width++;

			for (fr = 0; fr < searcher.m; fr++) {
				quant[fr].corr = corr(&(elution[0]), &(searcher.ms2[fr * len]), len);
				quant[fr].quantity[qTotal] = sum(&(searcher.ms2[fr * len]), len);
				for (pos = 0, quant[fr].quantity[qScaled] = 0.0; pos < len; pos++)
					quant[fr].quantity[qScaled] += searcher.ms2[fr * len + pos] * Sqr(elution[pos]);
				quant[fr].quantity[qScaled] /= w;
			}

			for (fr = 0, quantity = 0.0; fr < searcher.m; fr++) {
				assert(quant[best_fr].quantity[qScaled] > E);
				r = quant[fr].quantity[qScaled] / quant[best_fr].quantity[qScaled];
				for (pos = 0, quant[fr].quantity[qFiltered] = 0.0; pos < len; pos++)
					quant[fr].quantity[qFiltered] += Min(searcher.ms2[fr * len + pos], FilterFactor * r * elution[pos]);
				quantity += quant[fr].quantity[qFiltered];
			}
			free();
        }

		void intensities(int peak) { // get fragment intensities
			assert(peak >= 0);
			if (!initialised) build_index();
			int fr, k = peak, len = 1;
			double w;
			Searcher searcher(this);
			searcher.chromatogram<false>(k, k + 1, false, true);
			for (fr = 0, w = 0.0; fr < fragments.size(); fr++) w += (fragments[fr].height = searcher.ms2[fr * len]);
			assert(w > E);
			for (fr = 0; fr < fragments.size(); fr++) fragments[fr].height /= w;
		}

		void combine_scores(int _thread_id) {
			int i, j;
			if (run->curr_iter <= nnIter || GlobalTraining)
				for (i = 0; i < peaks.size(); i++) 
					peaks[i].score = scalar(scores, run->weights, run->par_learn, pN);
			else {
				Matrix* mat = createMatrix(peaks.size(), pN, &(scoring[0]));
				for (i = 0; i * Threads + _thread_id < nnBagging; i++) {
					forwardPassOnData(net[i * Threads + _thread_id].network, mat);
					Matrix* result = getOuput(net[i * Threads + _thread_id].network);
					while (!lock.set()) {}
					for (j = 0; j < peaks.size(); j++) peaks[j].score += result->data[2 * j];
					lock.free();
				}
			}
		}

        void seek() {
            int k, peak;
            double score;
            combined_score = -INF, best_peak = -1;

            for (peak = 0; peak < peak_number; peak++) {
                if (run->curr_iter <= nnIter || GlobalTraining) score = run->combined_score<true>(&(scoring[peak * pN]));
				else score = peaks[peak].score;
                if (score > combined_score) {
					if (run->nRT_windowed_search && Abs(run->predicted_nRT[peaks[peak].apex] - nRT) > run->nRT_window) continue;
                    best_peak = peak;
                    combined_score = score;
                }
            }
            found = (best_peak >= 0);
			if (found) {
				best_fragment = peaks[best_peak].best_fragment;
				apex = peaks[best_peak].apex;
				peak_pos = peaks[best_peak].peak;
				predicted_nRT = run->predicted_nRT[apex];
				RT = run->scan_RT[apex];
				for (k = 0; k < pN; k++) scores[k] = scoring[best_peak * pN + k];
				if (!decoy) if (run->curr_iter >= iN - 1 || (run->curr_iter == CalibrationIter && (InferWindow || Calibrate))) quantify(peaks[best_peak].peak);
			}
			
        }
        
        void find(int _thread_id, bool _free, bool skip) {
			if (lock.set()) {
				if (skip) {
					score_RT();
					return;
				}
				thread_id = _thread_id;
				if (!initialised) quant.resize(fragments.size());
				if (!QuantOnly) {
					search();
					seek();
				} else quantify(peak_pos);
				if (_free) free();
			}
        }
    };

	class Target {
	public:
		Precursor target;
		Precursor decoy;
		bool training, test, exclude;
		int batch;
	};

	std::vector<Target> entries;

	void write(char * file_name) {
		std::ofstream out(file_name, std::ofstream::binary);

		int size = scans.size();
		out.write((char*)&size, sizeof(int));
		for (int i = 0; i < scans.size(); i++) scans[i].write(out);
		size = ms1.size();
		out.write((char*)&size, sizeof(int));
		for (int i = 0; i < ms1.size(); i++) ms1[i].write(out);

		out.close();
	}

	void read(char * file_name) {
		std::ifstream in(file_name, std::ifstream::binary);

		int size; in.read((char*)&size, sizeof(int));;
		scans.resize(size);
		for (int i = 0; i < scans.size(); i++) scans[i].read(in);
		in.read((char*)&size, sizeof(int));
		ms1.resize(size);
		for (int i = 0; i < ms1.size(); i++) 
			ms1[i].read(in);

		in.close();
	}

	void init() { // called after read()
		int i;
		n_scans = scans.size();
		if (ExportWindows) {
			if (Verbose >= 1) std::cout << "Exporting window acquisition scheme...\n";
			std::ofstream out(name + std::string(".txt"), std::ofstream::out);
			out << "lower_offset\tupper_offset\n";
			for (i = 0; i < Min(n_scans, MaxCycleLength); i++) out << scans[i].window_low << "\t" << scans[i].window_high << "\n";
			out.close();
		}
		for (i = 1; i < n_scans; i++) { // handle overlapping windows
			if (scans[i - 1].window_low < scans[i].window_low - E && scans[i - 1].window_high < scans[i].window_high - E && scans[i - 1].window_high > scans[i].window_low + E)
				scans[i - 1].window_high = scans[i].window_low = (scans[i - 1].window_high + scans[i].window_low) * 0.5;
			if (scans[i - 1].window_low > scans[i].window_low + E && scans[i - 1].window_high > scans[i].window_high + E && scans[i - 1].window_low < scans[i].window_high - E)
				scans[i - 1].window_low = scans[i].window_high = (scans[i - 1].window_low + scans[i].window_high) * 0.5;
		}
		predicted_nRT.resize(n_scans), scan_RT.resize(n_scans);
		for (i = 0; i < n_scans; i++) predicted_nRT[i] = 0.0, scan_RT[i] = scans[i].RT;

		ms1_RT.resize(ms1.size());
		for (i = 0; i < ms1_RT.size(); i++) ms1_RT[i] = ms1[i].RT;

		for (int i = 0; i < ms1.size(); i++) { // calculate ms1 windows
			ms1[i].window_low = Min(ms1[i].peaks[0].mz, ms1[i].peaks[ms1[i].peaks.size() - 1].mz) - E * 0.25;
			ms1[i].window_high = Max(ms1[i].peaks[0].mz, ms1[i].peaks[ms1[i].peaks.size() - 1].mz) + E * 0.25;
			if (i > 0) {
				if (ms1[i - 1].window_low < ms1[i].window_low - E && ms1[i - 1].window_high < ms1[i].window_high - E && ms1[i - 1].window_high > ms1[i].window_low + E)
					ms1[i - 1].window_high = ms1[i].window_low = (ms1[i - 1].window_high + ms1[i].window_low) * 0.5;
				if (ms1[i - 1].window_low > ms1[i].window_low + E && ms1[i - 1].window_high > ms1[i].window_high + E && ms1[i - 1].window_low < ms1[i].window_high - E)
					ms1[i - 1].window_low = ms1[i].window_high = (ms1[i - 1].window_low + ms1[i].window_high) * 0.5;
			}
		}
	}

	bool load(char * file_name) { // Spectra in the file should be ordered by the acquisition time
		name = std::string(file_name);
		if (Verbose >= 1) std::cout << "Loading run " << name << "...\n";

		if (name.find(std::string(".dia")) != std::string::npos) {
			read(file_name);
			goto finalise;
		}

#ifdef MZML
		{
			MSToolkit::MSReader r;
			MSToolkit::Spectrum s;

			r.setFilter(MSToolkit::MS1);
			r.addFilter(MSToolkit::MS2);
			if (!r.readFile(file_name, s, 1)) {
				std::cerr << "cannot read the file\n";
				return false;
			}

			while (true) {
				if (!r.readFile(NULL, s)) break;
				if (s.getScanNumber() == 0) break;

				if (s.getMsLevel() == 2) scans.push_back(Tandem(s));
				else if (s.getMsLevel() == 1) ms1.push_back(Tandem(s));
			}
		}
#endif

		finalise:
		init();

		if (Verbose >= 2) std::cout << "Run loaded\n";
		return true;
	}
    
    void load_library(Library * _lib) { 
		lib = _lib;
		int i, ls, pos = 0, cnt = 0;
		Target entry;
		std::vector<bool> has(lib->entries.size());
		for (auto it = lib->entries.begin(); it != lib->entries.end(); it++, pos++) {
			double mz = it->target.mz;
			int max = Min(scans.size(), MaxCycleLength);
			for (i = 0; i < max; i++) 
				if (scans[i].has(mz)) break;
			if (i < max) has[pos] = true, cnt++;
		}
		entries.resize(cnt);
		if (Verbose >= 1) std::cout << cnt << " library precursors are potentially detectable\n";
		cnt = pos = 0;
		for (auto it = lib->entries.begin(); it != lib->entries.end(); it++, pos++) {
			if (has[pos]) {
				entries[cnt].target.init(this, it->target);
				entries[cnt].decoy.init(this, it->decoy);
				entries[cnt].decoy.decoy = true;
				entries[cnt].target.exclude = entries[cnt].decoy.exclude = it->exclude;
				entries[cnt].target.length = entries[cnt].decoy.length = peptide_length(it->name);
				cnt++;
			}
		}

		// select precursors that will not be used for classifier training and assign batches
		cnt = 0;
		std::mt19937_64 gen(1);
		if (TestDataset && !in_ref_run) {
			for (pos = ls = 0; pos < entries.size(); pos++) if (!entries[pos].exclude) ls++;
			int max = Max(1, Min(ls / 2, (int)(TestDatasetSize * (double)ls)));
			while (cnt < max) {
				pos = gen() % (unsigned long long)ls;
				if (entries[pos].exclude) continue;
				if (!entries[pos].test) entries[pos].test = true, cnt++;
			}
		}

		if (BatchMode) Batches = Min(MaxBatches, Max(1, entries.size() / MinBatch));
		std::vector<int> index(entries.size());
		for (pos = 0; pos < index.size(); pos++) index[pos] = pos;
		std::shuffle(index.begin(), index.end(), gen);
		for (pos = 0; pos < entries.size(); pos++) {
			if (!entries[pos].exclude && !entries[pos].test) entries[pos].training = true;
			entries[pos].batch = index[pos] % Batches;
		}
	}

	Target * find_entry(int index) {
		int low = 0, high = entries.size();
		while (high > low) {
			int mid = (low + high) / 2;
			int i = entries[mid].target.index;
			if (i < index) low = mid + 1;
			else if (i > index) high = mid;
			else return &(entries[mid]);
		}
		return NULL;
	}
    
    void process_precursors(int thread_id, int min_batch) {
        for (int i = 0; i < entries.size(); i++) {  
			if (entries[i].batch > curr_batch) continue;
			if (!QuantOnly || entries[i].target.found) entries[i].target.find(thread_id, min_batch, entries[i].batch < min_batch);
			if (!QuantOnly) entries[i].decoy.find(thread_id, min_batch, entries[i].batch < min_batch);
        }
    }

	void combine_scores(int thread_id) {
		for (int i = 0; i < entries.size(); i++) {
			if (entries[i].target.found) entries[i].target.combine_scores(thread_id);
			if (entries[i].decoy.found) entries[i].decoy.combine_scores(thread_id);
		}
	}
    
    void seek_precursors(bool free, int min_batch) {
		if (Verbose >= 2) std::cout << "Precursor search...\n";
        
        int i;
        if (Threads > 1) {
            std::vector<std::thread> threads;
            for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Run::process_precursors, this, i, min_batch));
            for (i = 0; i < Threads; i++) threads[i].join();
        } else process_precursors(0, min_batch);
		for (i = 0; i < entries.size(); i++) {
			entries[i].target.lock.free(), entries[i].decoy.lock.free();
			if (free) entries[i].target.free(), entries[i].decoy.free();
		}
    }

	void free_precursors() {
		for (int i = 0; i < entries.size(); i++)
			entries[i].target.free(), entries[i].decoy.free();
	}

	void reset_precursors() {
		for (int i = 0; i < entries.size(); i++)
			entries[i].target.free(), entries[i].decoy.free(), entries[i].target.found = entries[i].decoy.found = false;
	}

	void standardise_scores() {
		if (st_calculated || standardised) goto scale;
		if (Verbose >= 1) std::cout << "Calculating standardisation coefficients...\n";

		int i, cnt = 0;
		for (i = 0; i < pN; i++) st_e[i] = st_sd[i] = 0.0;
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->target.found) {
				cnt++;
				for (i = 0; i < pN; i++) st_e[i] += it->target.scores[i];
			}
			if (it->decoy.found) {
				cnt++;
				for (i = 0; i < pN; i++) st_e[i] += it->decoy.scores[i];
			}
		}
		if (cnt <= 2) return;
		for (i = 0; i < pN; i++) st_e[i] /= (double)cnt;
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->target.found) for (i = 0; i < pN; i++) st_sd[i] += Sqr(it->target.scores[i] - st_e[i]);
			if (it->decoy.found) for (i = 0; i < pN; i++) st_sd[i] += Sqr(it->decoy.scores[i] - st_e[i]);
		}
		for (i = 0; i < pN; i++) st_sd[i] = sqrt(Max(0.0, st_sd[i] / (double)(cnt - 1)));

	scale:
		if (Verbose >= 1) std::cout << "Applying standardisation...\n";
		for (i = 0; i < pN; i++) {
			if (standardised && i != pRT && i != pnRT) continue; // only RT scores are updated each iteration
			if (st_sd[i] < E) {
				for (auto it = entries.begin(); it != entries.end(); it++) {
					it->target.scores[i] = nnStScale * (it->target.scores[i] - st_e[i]);
					it->decoy.scores[i] = nnStScale * (it->decoy.scores[i] - st_e[i]);
					for (int j = 0; j < it->target.peaks.size(); j++) it->target.scoring[j * pN + i] = nnStScale * (it->target.scoring[j * pN + i] - st_e[i]);
					for (int j = 0; j < it->decoy.peaks.size(); j++) it->decoy.scoring[j * pN + i] = nnStScale * (it->decoy.scoring[j * pN + i] - st_e[i]);
				}
			} else {
				for (auto it = entries.begin(); it != entries.end(); it++) {
					it->target.scores[i] = nnStScale * (it->target.scores[i] - st_e[i]) / st_sd[i];
					it->decoy.scores[i] = nnStScale * (it->decoy.scores[i] - st_e[i]) / st_sd[i];
					for (int j = 0; j < it->target.peaks.size(); j++) it->target.scoring[j * pN + i] = nnStScale * (it->target.scoring[j * pN + i] - st_e[i]) / st_sd[i];
					for (int j = 0; j < it->decoy.peaks.size(); j++) it->decoy.scoring[j * pN + i] = nnStScale * (it->decoy.scoring[j * pN + i] - st_e[i]) / st_sd[i];
				}
			}
		}
		standardised = true;
		if (GlobalNN) st_calculated = true; // calculate standardisation coefficients using the first run only
	}

    void fit_weights() { // LDA
        int i, j, pos;
        double av[pN];
        std::vector<float> delta(entries.size() * pN);
        Eigen::MatrixXd A((int) pN, (int) pN);
        Eigen::VectorXd b((int) pN);

		int maxTraining = entries.size() * 2;
		if (curr_iter >= nnIter && !GlobalNN) if (training[0].size() < maxTraining) {
			for (i = 0; i < nnBagging; i++) training[i].resize(maxTraining), training_classes[i].resize(maxTraining);
			test.resize(maxTraining), test_classes.resize(maxTraining);
		}

		if (Verbose >= 2) std::cout << "Optimising weights...\n";
        
		pos = 0;
        for (i = 0; i < pN; i++) {
            av[i] = 0.0;
            for (j = 0; j < pN; j++) A(i, j) = 0.0;
        }
		int nnCnt = 0, nnTest = 0;
        for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->batch > curr_batch) continue;
			if (!GlobalNN) {
				if (curr_iter >= nnIter && it->training) {
					if (it->target.found) {
						for (i = 0; i < nnBagging; i++) training[i][nnCnt] = it->target.scores, training_classes[i][nnCnt] = target_nn;
						nnCnt++;
					}
					if (it->decoy.found) {
						for (i = 0; i < nnBagging; i++) training[i][nnCnt] = it->decoy.scores, training_classes[i][nnCnt] = decoy_nn;
						nnCnt++;
					}
				}

				if (curr_iter >= nnIter && Verbose >= 3 && TestDataset && it->test) {
					if (it->target.found) test[nnTest] = it->target.scores, test_classes[nnTest++] = target_nn;
					if (it->decoy.found) test[nnTest] = it->decoy.scores, test_classes[nnTest++] = decoy_nn;
				}
			}

            if (!it->training) continue;
            if (!it->target.found || !it->decoy.found) continue;

            for (i = 0; i < pN; i++) {
                double x = par_learn[i] ? it->target.scores[i] : 0.0;
				double y = par_learn[i] ? it->decoy.scores[i] : 0.0;
                av[i] += (delta[pos * pN + i] = x - y);
            }

            pos++;
        }
		if (!pos) {
			Warning("no training precursors");
			if (curr_iter == nnIter) use_nn = false;
			return;
		}
		if (in_ref_run && entries.size() < 2 * pN) return;
	
        for (i = 0; i < pN; i++) av[i] /= (double)pos;
		pos = 0;
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (!it->training) continue;
			if (it->batch > curr_batch) continue;
			if (!it->target.found || !it->decoy.found) continue;

			for (i = 0; i < pN; i++)
				for (j = i; j < pN; j++) A(j, i) = (A(i, j) += (delta[pos * pN + i] - av[i]) * (delta[pos * pN + j] - av[j]));
			pos++;
		}
		if (pos > 1) for (i = 0; i < pN; i++)
			for (j = i; j < pN; j++) A(j, i) = (A(i, j) /= (double)(pos - 1));
        
		if (Verbose >= 3) {
			std::cout << "Averages: \n";
			for (i = 0; i < pN; i++) std::cout << av[i] << " ";
			std::cout << "\n";
		}

		for (i = 0; i < pN; i++) b(i) = av[i];
		solve(weights, A, b, pN);
      
		if (Verbose >= 3) {
			std::cout << "Weights: \n";
			for (i = 0; i < pN; i++) std::cout << weights[i] << " ";
			std::cout << "\n";
		}

		if (curr_iter != nnIter || GlobalNN) return;
		if (Verbose >= 1) std::cout << "Training the neural network...\n";

		size_t* hiddenSize = (size_t*)alloca(nnHidden * sizeof(size_t));
		Activation* hiddenActivation = (Activation*)alloca(nnHidden * sizeof(Activation));
		for (i = 0; i < nnHidden; i++) {
			hiddenSize[i] = (nnHidden - i) * 5;
			hiddenActivation[i] = tanH;
		}

		for (i = 0; i < nnBagging; i++) {
			DataSet* trainingData = createDataSet(nnCnt, pN, &(training[i][0]));
			DataSet* trainingClasses = createDataSet(nnCnt, 2, &(training_classes[i][0]));

			if (net[i].network == NULL)
				net[i].network = createNetwork(pN, nnHidden, hiddenSize, hiddenActivation, 2, softmax, net[i].random);
			net[i].seed(i);
			net[i].lossFunction = CROSS_ENTROPY_LOSS;
			net[i].batchSize = Min(50, Max(1, nnCnt / 100));
			net[i].learningRate = Min(10.0, Max(0.1, nnLearning * (double)nnCnt));
			net[i].searchTime = 0;
			net[i].regularizationStrength = Regularisation;
			net[i].momentumFactor = 0.9;
			net[i].shuffle = 1;
			net[i].verbose = (Verbose >= 4);
			net[i].data = trainingData;
			net[i].classes = trainingClasses;
			net[i].maxIters = (nnCnt * nnEpochs) / net[i].batchSize;
		}

		std::vector<std::thread> thr;
		for (i = 0; i < Threads; i++) thr.push_back(std::thread(train, i));
		for (i = 0; i < Threads; i++) thr[i].join();

		if (Verbose >= 3 && TestDataset) for (i = 0; i < nnBagging; i++) {
			DataSet* testData = createDataSet(nnTest, pN, &(test[0]));
			DataSet* testClasses = createDataSet(nnTest, 2, &(test_classes[0]));
			float loss = crossEntropyLoss(net[i].network, testData, testClasses);
			std::cout << "Test loss: " << loss << "\n";
		}
	}
    
    int calculate_qvalues() {
		if (Verbose >= 2) std::cout << "Calculating q-values...\n";
		nRTTopPrecursors = (!TestDataset || in_ref_run) ? nRTTargetTopPrecursors : (TestDatasetSize * (double)nRTTargetTopPrecursors);
        
        int i = 0, ids = 0;
        std::set<std::pair<float, int> > decoy_scores;
        std::set<std::pair<float, int> > target_scores;
        min_target_decoy_score = INF;
        nRT_cscore = nRT_ref_score = INF;
        
        for (auto it = entries.begin(); it != entries.end(); it++, i++) {
			if (it->batch > curr_batch) continue;
			it->target.combined_score = combined_score<false>(it->target.scores);
			it->decoy.combined_score = combined_score<false>(it->decoy.scores);
            
            if (it->target.found && it->target.combined_score < min_target_decoy_score) 
                min_target_decoy_score = it->target.combined_score;
            if (it->decoy.found && it->decoy.combined_score < min_target_decoy_score) 
                min_target_decoy_score = it->decoy.combined_score;
            
            if (!TestDataset || in_ref_run || it->test) {
                if (it->target.found) target_scores.insert(std::pair<float, int>(-it->target.combined_score, i));
                if (it->decoy.found) decoy_scores.insert(std::pair<float, int>(-it->decoy.combined_score, i));
            }
        }
        
        for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->batch > curr_batch) continue;
            if (!it->target.found) it->target.combined_score = min_target_decoy_score;
            if (!it->decoy.found) it->decoy.combined_score = min_target_decoy_score;
        }

        std::vector<std::pair<float, int> > decoys(decoy_scores.begin(), decoy_scores.end());
        std::vector<std::pair<float, int> > targets(target_scores.begin(), target_scores.end());
        
		std::map<float, float> cs_qv;

        for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->batch > curr_batch) continue;
            if (!it->target.found) it->target.qvalue = 1.0;
            else {
                auto pD = std::upper_bound(decoys.begin(), decoys.end(), 
                        std::pair<float, int>(-it->target.combined_score, 0));
                auto pT = std::upper_bound(targets.begin(), targets.end(), 
                        std::pair<float, int>(-it->target.combined_score, 0));
				int n_targets = std::distance(targets.begin(), pT);
                it->target.qvalue = Min(1.0, (double(std::distance(decoys.begin(), pD))) / (double) Max(1, n_targets));
				if (n_targets <= nRTTopPrecursors && (!GuideLibrary || !lib->entries[it->target.index].from_fasta))
					if (it->target.combined_score < nRT_cscore) nRT_cscore = it->target.combined_score;
				if (n_targets <= nRTRefTopPrecursors && it->target.combined_score < nRT_ref_score) nRT_ref_score = it->target.combined_score;

				auto pair = std::pair<float, float>(it->target.combined_score, it->target.qvalue);
				auto pos = cs_qv.insert(pair);
				if (pos.second) {
					if (pos.first != cs_qv.begin() && std::prev(pos.first)->second < pair.second) pair.second = std::prev(pos.first)->second;
					else for (auto jt = std::next(pos.first); jt != cs_qv.end() && jt->second > pair.second; jt++) jt->second = pair.second;
				}
            }
        }

		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (it->batch > curr_batch) continue;
			if (!it->target.found) continue;
			if (it->target.qvalue < 1.0) {
				auto pos = cs_qv.lower_bound(it->target.combined_score);
				it->target.qvalue = pos->second;
			}
			if (!TestDataset || in_ref_run || it->test) if (it->target.qvalue <= 0.01) ids++;
		}
		if (Verbose >= 2) std::cout << "Number of test set IDs at 0.01 FDR: " << ids << "\n";
		if (TestDataset) {
			ids = (int)(((double)ids) / TestDatasetSize);
			if (Verbose >= 1 && curr_iter == iN - 1) std::cout << "Estimated number of precursors at 0.01 FDR: " << ids << "\n";
		}
		return ids;
    }

	int update_classifier(bool fit, bool free, int min_batch) {
		seek_precursors(free, min_batch);
		if (fit) fit_weights();
		int i, ids = calculate_qvalues();
		if (ids > max_ids) {
			max_ids = ids;
			for (i = 0; i < pN; i++) best_weights[i] = weights[i];
		} else if (ids < max_ids && fit) {
			if (Verbose >= 3) std::cout << "Trying the previous set of weights\n";
			for (i = 0; i < pN; i++) weights[i] = best_weights[i];
			ids = calculate_qvalues();
			if (ids > max_ids) {
				max_ids = ids;
				for (i = 0; i < pN; i++) best_weights[i] = weights[i];
			}
		}
		return ids;
	}
    
    void update(bool set_nRT_window, bool set_scan_window, bool calibrate, bool gen_ref) {
		if (Verbose >= 2) std::cout << "Calibrating retention times...\n";
		int peak_width = 0, peak_cnt = 0, mass_cnt = 0, mass_cnt_ms1 = 0;

		rt_stats.clear();  
		rt_stats.reserve(entries.size());
		rt_ref.clear();
		rt_ref.reserve(entries.size());
		rt_coo.clear();
		rt_coo.reserve(entries.size());
        
		for (auto it = entries.begin(); it != entries.end(); it++) {
            if (!it->target.found) continue;
			if ((it->target.qvalue <= MassCalQvalue || it->target.combined_score >= nRT_cscore) && Abs(it->target.mass_delta) > E) {
				mass_cnt++;
				if (calibrate) rt_coo.push_back(it->target.RT);
				if (Abs(it->target.mass_delta_ms1) > E) mass_cnt_ms1++;
			}
			if (GuideLibrary && lib->entries[it->target.index].from_fasta) continue;
			if (it->target.qvalue <= nRTMaxQvalue || it->target.combined_score >= nRT_cscore) {
				int index = std::distance(entries.begin(), it);
				rt_stats.push_back(index);
				if (it->target.qvalue <= nRTMaxQvalue && it->target.combined_score >= nRT_cscore)
					peak_width += it->target.peak_width, peak_cnt++;
				if (gen_ref && it->target.combined_score >= nRT_ref_score) rt_ref.push_back(it->target.index);
 			}
        }
		std::sort(rt_stats.begin(), rt_stats.end(), [&](const auto& lhs, const auto& rhs) { return entries[lhs].target.RT < entries[rhs].target.RT; });
        rt_data.resize(rt_stats.size());
		for (int i = 0; i < rt_stats.size(); i++) rt_data[i] = std::pair<float, float>(entries[rt_stats[i]].target.RT, entries[rt_stats[i]].target.nRT);
		if (Verbose >= 2) std::cout << rt_data.size() << " precursors used for nRT estimation.\n";

		if (rt_data.size() >= 2) {
			std::vector<double> coeff;
			spline(&predicted_nRT, &scan_RT, coeff, rt_data, Min(RTSegments, Max(1, 2.0 * sqrt(rt_data.size() / MinRTPredBin))));
		} else return;

		if (set_nRT_window) {
			rt_delta.resize(rt_stats.size());
			for (int i = 0; i < rt_stats.size(); i++)
				rt_delta[i] = -Abs(entries[rt_stats[i]].target.nRT - predicted_nRT[entries[rt_stats[i]].target.apex]);
			std::sort(rt_delta.begin(), rt_delta.end());
			nRT_window = Max(-nRTWindowMargin * rt_delta[(int)(nRTWindowLoss * (double)rt_delta.size())], (lib->nRT_max - lib->nRT_min) / double(in_ref_run ? nRTRefWindowFactor : nRTWindowFactor));
			nRT_windowed_search = (in_ref_run ? rt_delta.size() >= nRTWinSearchMinRef : rt_delta.size() >= nRTWinSearchMinCal);
			if (!nRT_windowed_search) { Warning("not enough confidently identified precursors for RT-windowed search"); }
			else if (Verbose >= 1) std::cout << "nRT_window set to " << nRT_window << "\n";
		}

		if (gen_ref) {
			std::sort(rt_ref.begin(), rt_ref.end());
			lib->save(gen_ref_file, &rt_ref, false);
		}

		if (set_scan_window) {
			PeakWidth = (double(peak_width)) / (double)Max(1, peak_cnt);
			if (Verbose >= 1) std::cout << "Peak width: " << PeakWidth << "\n";

			if (!IndividualWindows) WindowRadius = Max(1, int(ScanScale * PeakWidth)), InferWindow = false;
			if (Verbose >= 1) std::cout << "Scan window radius set to " << WindowRadius << "\n";
			window_calculated = true;
		}

		if (calibrate) {
			MassAccuracy = GlobalMassAccuracy;
			MassAccuracyMs1 = GlobalMassAccuracyMs1;
			if (mass_cnt >= MinMassDeltaCal) {
				std::vector<double> b(mass_cnt);
				typedef Eigen::Triplet<double> T;
				std::vector<T> TL, TL_ms1;

				std::sort(rt_coo.begin(), rt_coo.end());
				MassCalBins = Max(1, Min(MassCalBins, mass_cnt / MinMassCalBin));
				MassCorrection.resize(1 + 2 * MassCalBins);
				double split = ((double)(mass_cnt - 1)) / (double)MassCalBins;
				for (int i = 0; i < MassCalBins; i++) MassCalSplit[i] = rt_coo[Min(split * (double)i, mass_cnt - 1)]; MassCalSplit[MassCalBins] = rt_coo[mass_cnt - 1];
				for (int i = 0; i < MassCalBins; i++) MassCalCenter[i] = 0.5 * (MassCalSplit[i] + MassCalSplit[i + 1]);
				for (int i = 1; i < MassCalBins; i++) if (MassCalCenter[i] - MassCalCenter[i - 1] < E) {
					for (int k = i - 1; k < MassCalBins - 1; k++) MassCalCenter[k] = MassCalCenter[k + 1];
					MassCalBins--;
				}

				mass_cnt = 0;
				for (auto it = entries.begin(); it != entries.end(); it++) {
					if (!it->target.found) continue;
					if ((it->target.qvalue <= MassCalQvalue || it->target.combined_score >= nRT_cscore) && Abs(it->target.mass_delta) > E) {
						b[mass_cnt] = it->target.mass_delta;
						double mz = it->target.mass_delta_mz;
						double rt = it->target.RT;
						TL.push_back(T(mass_cnt, 0, Sqr(mz)));

						if (rt <= MassCalCenter[0])
							TL.push_back(T(mass_cnt, 1, 1.0)), TL.push_back(T(mass_cnt, 2, mz));
						else if (rt >= MassCalCenter[MassCalBins - 1])
							TL.push_back(T(mass_cnt, 1 + (MassCalBins - 1) * 2, 1.0)), TL.push_back(T(mass_cnt, 2 + (MassCalBins - 1) * 2, mz));
						else for (int i = 1; i < MassCalBins; i++) if (rt < MassCalCenter[i]) {
							double u = rt - MassCalCenter[i - 1], v = MassCalCenter[i] - rt, w = u + v;
							if (w > E) {
								double cl = v / w, cr = u / w;
								TL.push_back(T(mass_cnt, 1 + (i - 1) * 2, cl)), TL.push_back(T(mass_cnt, 2 + (i - 1) * 2, mz * cl));
								TL.push_back(T(mass_cnt, 1 + i * 2, cr)), TL.push_back(T(mass_cnt, 2 + i * 2, mz * cr));
							}
							break;
						}

						mass_cnt++;
					}
				}

				auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), b.size());
				Eigen::SparseMatrix<double, Eigen::RowMajor> A(mass_cnt, MassCorrection.size());
				A.setFromTriplets(TL.begin(), TL.end());
				Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
				SQR.compute(A);
				auto X = SQR.solve(B);
				for (int i = 0; i < MassCorrection.size(); i++) MassCorrection[i] = X(i);

				if (Verbose >= 4) {
					std::cout << "Mass correction transform (" << mass_cnt << " precursors): \n";
					for (int i = 0; i < MassCorrection.size(); i++) std::cout << MassCorrection[i] << " ";
					std::cout << "\n";
				}

				if (mass_cnt < MinMassAccCal) { Warning("not enough confidently identified precursors for mass accuracy calibration"); }
				else if (!ForceMassAcc) { // mass accuracy calibration
					std::vector<double>b_r(mass_cnt);
					mass_cnt = 0;
					for (auto it = entries.begin(); it != entries.end(); it++) {
						if (!it->target.found) continue;
						if ((it->target.qvalue <= MassCalQvalue || it->target.combined_score >= nRT_cscore) && Abs(it->target.mass_delta) > E)
							b_r[mass_cnt++] = Abs(it->target.mass_delta_mz + it->target.mass_delta - predicted_mz(&(MassCorrection[0]), it->target.mass_delta_mz, it->target.RT)) / it->target.mass_delta_mz;
					}
					if (Verbose >= 4) std::cout << "M/z SD: " << sqrt(var(b_r)) * 1000000.0 << " ppm\n";
					std::sort(b_r.begin(), b_r.end());
					MassAccuracy = MinMassAcc = b_r[0.8 * (double)b_r.size()];
					if (mass_cnt <= MinFreeMassAcc) MassAccuracy = Min(GlobalMassAccuracy, MassAccuracy);
					if (Verbose >= 1) std::cout << "Mass accuracy: " << MassAccuracy * 1000000.0 << " ppm\n";
					acc_calibrated = true;
				}
			} else Warning("cannot perform mass calibration, too few confidently identified precursors");

			if (mass_cnt_ms1 >= MinMassDeltaCal) {
				std::vector<double> b(mass_cnt_ms1);
				typedef Eigen::Triplet<double> T;
				std::vector<T> TL, TL_ms1;

				std::sort(rt_coo.begin(), rt_coo.end());
				MassCalBinsMs1 = Max(1, Min(MassCalBinsMs1, mass_cnt_ms1 / MinMassCalBin));
				MassCorrectionMs1.resize(1 + 2 * MassCalBinsMs1);
				double split = ((double)(mass_cnt_ms1 - 1)) / (double)MassCalBinsMs1;
				for (int i = 0; i < MassCalBinsMs1; i++) MassCalSplitMs1[i] = rt_coo[Min(split * (double)i, mass_cnt_ms1 - 1)]; MassCalSplitMs1[MassCalBinsMs1] = rt_coo[mass_cnt_ms1 - 1];
				for (int i = 0; i < MassCalBinsMs1; i++) MassCalCenterMs1[i] = 0.5 * (MassCalSplitMs1[i] + MassCalSplitMs1[i + 1]);
				for (int i = 1; i < MassCalBinsMs1; i++) if (MassCalCenterMs1[i] - MassCalCenterMs1[i - 1] < E) {
					for (int k = i - 1; k < MassCalBinsMs1 - 1; k++) MassCalCenterMs1[k] = MassCalCenterMs1[k + 1];
					MassCalBinsMs1--;
				}

				mass_cnt_ms1 = 0;
				for (auto it = entries.begin(); it != entries.end(); it++) {
					if (!it->target.found) continue;
					if ((it->target.qvalue <= MassCalQvalue || it->target.combined_score >= nRT_cscore) && Abs(it->target.mass_delta) > E && Abs(it->target.mass_delta_ms1) > E) {
						b[mass_cnt_ms1] = it->target.mass_delta_ms1;
						double mz = it->target.mz;
						double rt = it->target.RT;
						TL.push_back(T(mass_cnt_ms1, 0, Sqr(mz)));

						if (rt <= MassCalCenterMs1[0])
							TL.push_back(T(mass_cnt_ms1, 1, 1.0)), TL.push_back(T(mass_cnt_ms1, 2, mz));
						else if (rt >= MassCalCenterMs1[MassCalBinsMs1 - 1])
							TL.push_back(T(mass_cnt_ms1, 1 + (MassCalBinsMs1 - 1) * 2, 1.0)), TL.push_back(T(mass_cnt_ms1, 2 + (MassCalBinsMs1 - 1) * 2, mz));
						else for (int i = 1; i < MassCalBinsMs1; i++) if (rt < MassCalCenterMs1[i]) {
							double u = rt - MassCalCenterMs1[i - 1], v = MassCalCenterMs1[i] - rt, w = u + v;
							if (w > E) {
								double cl = v / w, cr = u / w;
								TL.push_back(T(mass_cnt_ms1, 1 + (i - 1) * 2, cl)), TL.push_back(T(mass_cnt_ms1, 2 + (i - 1) * 2, mz * cl));
								TL.push_back(T(mass_cnt_ms1, 1 + i * 2, cr)), TL.push_back(T(mass_cnt_ms1, 2 + i * 2, mz * cr));
							}
							break;
						}

						mass_cnt_ms1++;
					}
				}

				auto B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), b.size());
				Eigen::SparseMatrix<double, Eigen::RowMajor> A(mass_cnt_ms1, MassCorrectionMs1.size());
				A.setFromTriplets(TL.begin(), TL.end());
				Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SQR;
				SQR.compute(A);
				auto X = SQR.solve(B);
				for (int i = 0; i < MassCorrectionMs1.size(); i++) MassCorrectionMs1[i] = X(i);

				if (Verbose >= 4) {
					std::cout << "MS1 mass correction transform (" << mass_cnt_ms1 << " precursors): \n";
					for (int i = 0; i < MassCorrectionMs1.size(); i++) std::cout << MassCorrectionMs1[i] << " ";
					std::cout << "\n";
				}

				if (mass_cnt_ms1 < MinMassAccCal) { Warning("not enough confidently identified precursors for MS1 mass accuracy calibration"); }
				else if (!ForceMassAcc) { // mass accuracy calibration
					std::vector<double>b_r(mass_cnt_ms1);
					mass_cnt_ms1 = 0;
					for (auto it = entries.begin(); it != entries.end(); it++) {
						if (!it->target.found) continue;
						if ((it->target.qvalue <= MassCalQvalue || it->target.combined_score >= nRT_cscore) && Abs(it->target.mass_delta) > E && Abs(it->target.mass_delta_ms1) > E)
							b_r[mass_cnt_ms1++] = Abs(it->target.mz + it->target.mass_delta_ms1 - predicted_mz_ms1(&(MassCorrectionMs1[0]), it->target.mz, it->target.RT)) / it->target.mz;
					}

					std::sort(b_r.begin(), b_r.end());
					MassAccuracyMs1 = MinMassAccMs1 = b_r[0.8 * (double)b_r.size()];
					if (mass_cnt_ms1 <= MinFreeMassAcc) MassAccuracyMs1 = Min(GlobalMassAccuracyMs1, MassAccuracyMs1);
					if (Verbose >= 1) std::cout << "MS1 mass accuracy: " << MassAccuracyMs1 * 1000000.0 << " ppm\n";
					acc_ms1_calibrated = true;
				}
			} else Warning("cannot perform MS1 mass calibration, too few confidently identified precursors");
			if (!IndividualMassAcc) ForceMassAcc = true;
			mz_calibrated = true;
		}
    }

	bool reference_run(Library * ref) {
		if (!nRTWindowedSearch && (!RefCal || (!InferWindow && !Calibrate))) return false;
		if (Verbose >= 1) std::cout << "Calibrating retention times using a set of reference precursors.\n";

		in_ref_run = true, curr_batch = Batches;
		load_library(ref);
		for (curr_iter = (entries.size() >= 2 * pN ? 0 : nRTWindowIter); curr_iter <= nRTWindowIter; curr_iter++) {
			for (int i = 0; i < pN; i++) {
				par_seek[i] = (curr_iter >= pars[i].min_iter_seek);
				par_learn[i] = (curr_iter >= pars[i].min_iter_learn);
			}
			seek_precursors(false, 0);
			fit_weights();
			calculate_qvalues();
		}
		update(nRTWindowedSearch, RefCal && InferWindow, RefCal && Calibrate, false);
		free_precursors();
		in_ref_run = false;

		return true;
	}
    
    void process() {
		int i, ids = 0, cal_batch = 0;

		if (QuantOnly) goto report;

		bool do_calibrate = (InferWindow && !window_calculated) || (Calibrate && !mz_calibrated);
		for (curr_batch = 0; curr_batch < Batches; curr_batch++) {
			if (Verbose >= 1) {
				if (Batches == 1) std::cout << "Processing...\n";
				else std::cout << "Processing batch #" << curr_batch + 1 << " out of " << Batches << " ...\n";
			}

			curr_iter = curr_batch ? 1 : 0;
			for (i = 0; i < pN; i++) {
				par_seek[i] = (curr_iter >= pars[i].min_iter_seek);
				par_learn[i] = (curr_iter >= pars[i].min_iter_learn);
			}

			ids = update_classifier(true, false, 0);
			update(false, false, false, false);
			cal_batch = curr_batch;
			if (ids >= MinCal && do_calibrate) break;
		}
		for (curr_iter = curr_iter + 1; curr_iter <= CalibrationIter; curr_iter++) update_classifier(true, false, 0);
		if (do_calibrate) {
			update(nRTWindowedSearch, InferWindow && !window_calculated, Calibrate && !mz_calibrated, false);
			reset_precursors();
		}
		max_ids = 0;
		for (curr_batch = 0; curr_batch < Batches; curr_batch++) {
			if (Verbose >= 1) {
				if (Batches == 1) std::cout << "Processing...\n";
				else std::cout << "Processing batch #" << curr_batch + 1 << " out of " << Batches << " ...\n";
			}

			curr_iter = CalibrationIter + 1;
			for (i = 0; i < pN; i++) {
				par_seek[i] = (curr_iter >= pars[i].min_iter_seek);
				par_learn[i] = (curr_iter >= pars[i].min_iter_learn);
			}

			ids = update_classifier(true, false, 0);
			update(false, false, false, false);
			if (ids >= MinClassifier) break;
			if (acc_calibrated && curr_batch >= cal_batch) break;
		}
		if (acc_calibrated) {
			double best_acc = MassAccuracy, best_acc_ms1 = MassAccuracyMs1, start_acc = MassAccuracy, start_acc_ms1 = MassAccuracyMs1;
			int best_ids = update_classifier(false, false, 0), fail = 0;
			while (true) {
				MassAccuracy *= 1.2;
				if (MassAccuracy > start_acc * 10.0) break;
				if (Verbose >= 3) std::cout << "Trying mass accuracy " << MassAccuracy * 1000000.0 << " ppm\n";
				reset_precursors();
				ids = update_classifier(false, true, 0);
				if (ids > best_ids) best_acc = MassAccuracy, best_ids = ids, fail = 0;
				else fail++;
				if (fail >= 3) break;
			}
			MassAccuracy = best_acc, fail = 0;
			if (Verbose >= 1) std::cout << "Optimised mass accuracy: " << MassAccuracy * 1000000.0 << " ppm\n";

			if (acc_ms1_calibrated) while (true) {
				MassAccuracyMs1 *= 1.2;
				if (MassAccuracyMs1 > start_acc_ms1 * 10.0) break;
				if (Verbose >= 3) std::cout << "Trying MS1 mass accuracy " << MassAccuracyMs1 * 1000000.0 << " ppm\n";
				reset_precursors();
				ids = update_classifier(false, true, 0);
				if (ids > best_ids) best_acc_ms1 = MassAccuracyMs1, best_ids = ids, fail = 0;
				else fail++;
				if (fail >= 3) break;
			}
			MassAccuracyMs1 = best_acc_ms1;
			if (Verbose >= 1) std::cout << "Optimised MS1 mass accuracy: " << MassAccuracyMs1 * 1000000.0 << " ppm\n";
			reset_precursors();

			max_ids = update_classifier(true, false, 0);
			update(false, false, false, false);

			GlobalMassAccuracy = MassAccuracy;
			GlobalMassAccuracyMs1 = MassAccuracyMs1;
		}

		int processed_batches = curr_batch;
		for (curr_iter = curr_iter + 1; curr_iter < nnIter - 1; curr_iter++) update_classifier(true, false, 0), update(false, false, false, false);
		if (Batches > 1) {
			assert(nnIter >= iN - 1);
			free_precursors();
		}
		curr_batch = Batches;
		update_classifier(true, false, processed_batches + 1);
		update(false, false, false, false);
		curr_iter = nnIter;
		if (curr_iter >= iN - 1) {
			free_precursors();
			if (curr_iter >= iN) goto report;
		}
		if (!GlobalNN || GlobalTraining) {
			fit_weights(); // train the neural network
			calculate_qvalues();
			update(false, false, false, GenRef && curr_iter == iN - 1);
		}
		for (curr_iter = curr_iter + 1; curr_iter < iN; curr_iter++) {
			if (!GlobalTraining && use_nn) {
				for (int i = 0; i < nnBagging; i++) destroyMatrix(net[i].network->layers[0]->input);
				std::vector<std::thread> thr;
				for (i = 0; i < Threads; i++) thr.push_back(std::thread(&Run::combine_scores, this, i));
				for (i = 0; i < Threads; i++) thr[i].join();
				for (int i = 0; i < nnBagging; i++) net[i].network->layers[0]->input = NULL;
			}
			update_classifier(true, false, 0), update(false, false, false, GenRef && curr_iter == iN - 1);
		}

	report:
		curr_iter = iN - 1;
		seek_precursors(false, 0);
		Quant quant;
		QuantEntry qe;
		for (i = 0; i < pN; i++) quant.weights[i] = weights[i];

        for (auto it = entries.begin(); it != entries.end(); it++) {
            if (!it->target.found) continue;

			qe.index = std::distance(entries.begin(), it);
			qe.stored_scores = GlobalNN ? 1 : 0;
			qe.window = it->target.S;
			qe.pr.decoy_found = it->decoy.found;
			qe.pr.apex = it->target.apex;
			qe.pr.peak = it->target.peak_pos;
			qe.pr.best_fragment = it->target.best_fragment;
			qe.pr.peak_width = it->target.peak_width;
			qe.pr.run_index = run_index;
			qe.pr.index = it->target.index;
			qe.pr.nRT = it->target.nRT;
			qe.pr.predicted_nRT = !QuantOnly ? it->target.predicted_nRT : predicted_nRT[it->target.apex];
			qe.pr.quantity = it->target.quantity;
			qe.pr.qvalue = it->target.qvalue;
			qe.pr.RT = scan_RT[it->target.apex];
			qe.fr.resize(it->target.fragments.size());
			for (i = 0; i < qe.fr.size(); i++) qe.fr[i] = it->target.quant[i];
			if (!QuantOnly) for (i = 0; i < pN; i++) {
				qe.sc.target_scores[i] = it->target.scores[i];
				qe.sc.decoy_scores[i] = it->decoy.scores[i];
			}
			quant.entries.push_back(qe);
        }

		std::ofstream out(name + std::string(".quant"), std::ofstream::binary);
		quant.write(out);
		out.close();
		if (Verbose >= 1) std::cout << "Quantification information saved to " << name + std::string(".quant") << ".\n\n";
    }

	void update_library() { // update fragment intensities in the spectral library
		for (auto &e : entries) {
			auto &le = lib->entries[e.target.index];
			if (!le.from_fasta) continue; // do not overwrite existing learn_from_library patterns
			if (le.best_run != run_index || le.qvalue > SpectralLibraryGenQvalue) continue;
			e.target.intensities(le.peak);
			for (int fr = 0; fr < e.target.fragments.size(); fr++) le.target.fragments[fr].height = e.target.fragments[fr].height;
		}
	}
};

int main(int argc, char *argv[]) {  
#ifdef _MSC_VER
	DWORD mode;
	HANDLE console = GetStdHandle(STD_INPUT_HANDLE);
	GetConsoleMode(console, &mode);
	SetConsoleMode(console, (mode & ~(ENABLE_QUICK_EDIT_MODE | ENABLE_MOUSE_INPUT | ENABLE_WINDOW_INPUT)) | ENABLE_PROCESSED_INPUT);
#endif
	std::cout.setf(std::ios::unitbuf);

	init_aas();
    arguments(argc, argv);

	if (Convert) {
        for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
			Run run(std::distance(ms_files.begin(), it));
			if (!run.load(&((*it)[0]))) continue;
			std::string(dia_name) = run.name + std::string(".dia");
			run.write(&(dia_name[0]));
		}
		return 0;
	}

	Fasta fasta;
	if (FastaSearch) {
		if (learn_lib_file.size()) learn_from_library(learn_lib_file);
		if (!fasta.load(&(fasta_file[0]))) Throw("Cannot load FASTA");
	}

	Library lib;
	if (FastaSearch) {
		if (lib_file.size() && lib.load(&(lib_file[0]))) GuideLibrary = true;
		lib.load(fasta);
	} else if (!lib.load(&(lib_file[0]))) Throw("Cannot load spectral library");
	if (!QuantOnly) lib.generate_decoys();

	Library ref;
	if (ref_file.size()) {
		if (!ref.load(&(ref_file[0]))) ref_file.clear();
		else ref.generate_decoys();
	}

	if (ReportOnly) goto cross_run;
	if (FastaSearch && QuantOnly) goto gen_spec_lib;
	if (QuantOnly) goto quant_only;

	// first loop
	first_loop:
	if (RTProfiling) if (Verbose >= 1) std::cout << "First pass: collecting retention time information\n";
    for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
        if ((UseRTInfo && RTProfiling) || UseQuant) {
            std::ifstream in((*it) + std::string(".quant"));
            if (in.is_open()) {
                in.close();
                continue;
            }
        }
		if (Verbose >= 1) std::cout << "File #" << std::distance(ms_files.begin(), it) << "\n";
        
		Run run(std::distance(ms_files.begin(), it));
        if (!run.load(&((*it)[0]))) Throw("Cannot load the file");
		if (ref_file.size()) run.reference_run(&ref);
		run.load_library(&lib);
        run.process();

		if (GenRef) {
			GenRef = false, ref_file = gen_ref_file;
			if (!ref.load(&(ref_file[0]))) ref_file.clear();
			else ref.generate_decoys();
		}
    }

	if (GlobalTraining) {
		global_training(&lib);
		GlobalTraining = false;
		goto first_loop;
	}
	if (!RTProfiling) goto gen_spec_lib;

	quant_only:
	if (QuantOnly) {
		if (Verbose >= 1) std::cout << "Quantification...\n";
		lib.info.load(&(lib), ms_files);
		for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
			int runN = std::distance(ms_files.begin(), it);
			Run run(runN);
			if (!run.load(&((*it)[0]))) Throw("Cannot load the file");
			run.load_library(&lib);
			for (auto pos = lib.info.raw[runN].entries.begin(); pos != lib.info.raw[runN].entries.end(); pos++) {
				auto e = &(run.entries[pos->index]);
				e->target.found = true;
				e->target.S = pos->window;
				e->target.apex = pos->pr.apex;
				e->target.RT = run.scan_RT[pos->pr.apex];
				e->target.peak_pos = pos->pr.peak;
				e->target.best_fragment = pos->pr.best_fragment;
				e->target.peak_width = pos->pr.peak_width;
				e->target.qvalue = pos->pr.qvalue;
			}
			run.process();
		}
		lib.info.clear();
		goto cross_run;
	}

	{
		// second loop: refine nRT predictions
		if (Verbose >= 1) std::cout << "Second pass: indentification and quantification\n";
		Profile profile(ms_files);
		if (RTProfiling) for (auto jt = profile.entries.begin(); jt != profile.entries.end(); jt++)
			if (jt->index >= 0) if (jt->pr.qvalue < MaxProfilingQvalue) {
				lib.entries[jt->pr.index].target.nRT = lib.entries[jt->pr.index].decoy.nRT = jt->pr.predicted_nRT;
				lib.entries[jt->pr.index].best_run = jt->pr.run_index;
				lib.entries[jt->pr.index].peak = jt->pr.peak;
				lib.entries[jt->pr.index].qvalue = jt->pr.qvalue;
			}

        for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
			if (UseQuant) {
				std::ifstream in((*it) + std::string(".quant"));
				if (in.is_open()) {
					in.close();
					continue;
				}
			}
			if (Verbose >= 1) std::cout << "Second pass: file #" << std::distance(ms_files.begin(), it) << "\n";

			Run run(std::distance(ms_files.begin(), it));
			if (!run.load(&((*it)[0]))) Throw("Cannot load the file");
			if (ref_file.size()) run.reference_run(&ref);
			run.load_library(&lib);
			run.process();

			if (GenRef) {
				GenRef = false, ref_file = gen_ref_file;
				if (!ref.load(&(ref_file[0]))) ref_file.clear();
				else ref.generate_decoys();
			}
		}
	}

gen_spec_lib:
	if (FastaSearch) { // generating spectral library
		MaxF = INF;
		if (Verbose >= 1) std::cout << "Generating spectral library...\n";
		lib.fragment(); // generate all likely detectable fragments

		// determine the best run for each precursor and collect information from this run
		Profile profile(ms_files);
		for (auto jt = profile.entries.begin(); jt != profile.entries.end(); jt++)
			if (jt->index >= 0) if (jt->pr.qvalue <= SpectralLibraryGenQvalue) {
				lib.entries[jt->pr.index].target.nRT = (GuideLibrary || UseLibnRT) ? jt->pr.predicted_nRT : jt->pr.RT;
				lib.entries[jt->pr.index].best_run = jt->pr.run_index;
				lib.entries[jt->pr.index].peak = jt->pr.peak;
				lib.entries[jt->pr.index].qvalue = jt->pr.qvalue;
			}

		// extract spectra from runs
		for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
			Run run(std::distance(ms_files.begin(), it));
			if (!run.load(&((*it)[0]))) Throw("Cannot load the file");
			run.load_library(&lib);
			run.update_library();
		}

		lib.save(out_file.substr(0, out_file.find_last_of(".")) + std::string(".lib.tsv"), NULL, false);
	}
    
cross_run:
	if (Verbose >= 1) std::cout << "Cross-run analysis...\n";
	lib.info.load(&(lib), ms_files);
	lib.info.quantify();
	lib.quantify_proteins(3, ProteinQuantQvalue);
	lib.report(out_file);

    return 0;
}


