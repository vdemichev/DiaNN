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

#include "../cranium/src/cranium.h"
#include "../eigen/Eigen/Dense"
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

#define Min(x, y) ((x) < (y) ? (x) : (y))
#define Max(x, y) ((x) > (y) ? (x) : (y))
#define Abs(x) ((x) >= 0 ? (x) : (-(x)))
#define Sqr(x) ((x) * (x))

const double E = 0.000000001;
const double INF = 10000000.0;

int Threads = 1;
int iN = 8;

const int MaxLibSize = 1000000;

bool RTProfiling = false;

bool TestDataset = true;
double TestDatasetSize = 0.25;

int ScanFactor = 150;
int WindowRadius = 0;
double ScanScale = 2.2;
bool InferWindow = true;
bool IndividualWindows = false;

int MaxDpPP = INF;

const int MaxCycleLength = 1000;
const double MaxWindowOverlap = 2.0;

const int nRTSmoothingFactor = 50;
const int nRTTopPrecursors = 10; // use at least this number of training precursors for nRT estimation
const int nRTRefPrecursors = 500; // min number of precursors for the reference precursor set generation
const double nRTMaxQvalue = 0.01; // only precursors with lower q value are used for nRT estimation

const int nRTWindowIter = 3;
const int nRTWindowFactor = 20;
const double nRTWindowLoss = 0.01;
const double nRTWindowMargin = 1.5;
double nRTWindow = 0.0;

int nnIter = 6;
int nnBagging = 12;
int nnEpochs = 50;
double nnLearning = 10.0;
bool GlobalNN = false;
bool GlobalTraining = false;
int nnGlobalEpochs = 50;
int nnHidden = 2;
int nnGlobalHidden = 5;
float Regularisation = 0.0;

bool nnStandardise = false;
const double nnStScale = 0.5;

double MassAccuracy = 20.0 / 1000000.0; 
double GeneratorAccuracy = 20.0 / 1000000.0;
double MaxProfilingQvalue = 0.01;
double MaxQuantQvalue = 0.01;
const double PeakApexEvidence = 0.9; // force peak apexes to be close to the detection evidence maxima

double FilterFactor = 1.5;

bool Convert = false; // .mzML to .dia conversion mode
bool UseRTInfo = false; // use .quant files created previously for RT profiling
bool UseQuant = false; // use .quant files created previously; implies UseRTInfo
bool QuantOnly = false; // quantification will be perfromed anew using identification info from .quant files
bool ReportOnly = false; // generate report from quant files
bool UseRefFile = false; // use a .ref library file created previously
bool UpdateRefFile = false; // update the .ref file with newly computed data
std::string args, lib_file, out_file = "quant.csv", ref_file, exclude_from_training, include_for_training;
std::vector<std::string> ms_files;

int Verbose = 1;

enum {
	libPG, libPr, libCharge, libPrMz,
	libnRT, libFrMz, libFrI, libIsDecoy,
	libCols
};

std::vector<std::string> library_headers = {
	" UniProtIds UniprotID \"UniprotID\" ",
	" IntModifiedPeptide FullUniModPeptideName \"FullUniModPeptideName\" ",
	" PrecursorCharge \"PrecursorCharge\" ",
	" PrecursorMz \"PrecursorMz\" ",
	" nRT iRT Tr_recalibrated \"Tr_recalibrated\" ",
	" FragmentMz ProductMz \"ProductMz\" ",
	" RelativeIntensity LibraryIntensity \"LibraryIntensity\" ",
	" Decoy decoy \"Decoy\" \"decoy\" "
};

std::vector<std::pair<std::string, float> > Modifications = {
	std::pair<std::string, float>("UniMod:4", (float)57.021464),
	std::pair<std::string, float>("Carbamidomethyl (C)", (float)57.021464),
	std::pair<std::string, float>("UniMod:35", (float)15.994915),
	std::pair<std::string, float>("Oxidation (M)", (float)15.994915),
	std::pair<std::string, float>("UniMod:1", (float)42.010565),
	std::pair<std::string, float>("Acetyl (Protein N-term)", (float)42.010565)
};

enum {
	outFile, outPG, outPGQ, outModSeq,
	outPrId, outCharge, outQv, outPrQ,
	outRT,
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
	"RT"
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

std::string trim(std::string& input) { return std::regex_replace(input, std::regex("^ +| +$"), ""); }

void arguments(int argc, char *argv[]) {
    int i, start, next, end;
	std::string prefix, ext;
    
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
            } else std::cerr << "failed to open cfg file\n";
        }
        else if (!memcmp(&(args[start]), "threads ", 8)) Threads = std::stoi(args.substr(start + 8, std::string::npos)), std::cout << "Thread number set to " << Threads << "\n";
		else if (!memcmp(&(args[start]), "verbose ", 8)) Verbose = std::stoi(args.substr(start + 8, std::string::npos));
        else if (!memcmp(&(args[start]), "f ", 2)) ms_files.push_back(trim(args.substr(start + 2, end - start - 2)));
        else if (!memcmp(&(args[start]), "lib ", 4)) lib_file = trim(args.substr(start + 4, end - start - 4));
		else if (!memcmp(&(args[start]), "ref ", 4)) ref_file = trim(args.substr(start + 4, end - start - 4)), UseRefFile = true;
		else if (!memcmp(&(args[start]), "out ", 4)) out_file = trim(args.substr(start + 4, end - start - 4));
		else if (!memcmp(&(args[start]), "library-headers ", 16)) {
			std::string word;
			std::stringstream list(trim(args.substr(start + 16, end - start - 16)));
			for (i = 0; std::getline(list, word, ','); i++) {
				if (i >= libCols) {
					std::cout << "\nWARNING: " << word << ": extra headers will be ignored\n";
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
					std::cout << "\nWARNING: " << word << ": extra headers will be ignored\n";
					break;
				}
				oh[i] = std::string(" ") + word + std::string(" ");
			}
		}
		else if (!memcmp(&(args[start]), "mod ", 4)) {
			std::string name, mz;
			std::stringstream list(trim(args.substr(start + 4, end - start - 4)));
			if (!std::getline(list, name, ',')) std::cout << "\nWARNING: no modification name, modification ignored\n";
			else if (!std::getline(list, mz, ',')) std::cout << "\nWARNING: no modification m/z, modification ignored\n";
			else Modifications.push_back(std::pair<std::string, float>(name, (float)std::stod(mz)));
		}
		else if (!memcmp(&(args[start]), "exclude ", 8)) exclude_from_training = trim(args.substr(start + 8, end - start - 8)),
			std::cout << "Precursors corresponding to proteins with [" << exclude_from_training << "] in the ID will be ignored when training the classifier\n";
		else if (!memcmp(&(args[start]), "include ", 8)) include_for_training = trim(args.substr(start + 8, end - start - 8)),
			std::cout << "Precursors corresponding to proteins without [" << include_for_training << "] in the ID will be ignored when training the classifier\n";
        else if (!memcmp(&(args[start]), "use-ref ", 8)) UseRefFile = true, std::cout << "Use of reference peptides turned on\n";
        else if (!memcmp(&(args[start]), "update-ref ", 11)) UpdateRefFile = true, std::cout << "The set of reference peptides will be updated\n";
		else if (!memcmp(&(args[start]), "window ", 7)) {
			WindowRadius = std::stoi(args.substr(start + 7, std::string::npos));
			if (WindowRadius <= 0) std::cout << "\nWARNING: scan window radius should be a positive integer\n";
			else InferWindow = false, std::cout << "Scan window radius set to " << WindowRadius << "\n";
		}
		else if (!memcmp(&(args[start]), "nRT-window ", 11)) {
			nRTWindow = std::stod(args.substr(start + 11, std::string::npos));
			std::cout << "RT window radius set to " << nRTWindow << "\n";
		}
		else if (!memcmp(&(args[start]), "no-window-inference ", 20)) InferWindow = false, std::cout << "Scan window inference turned off\n";
		else if (!memcmp(&(args[start]), "individual-windows ", 19)) IndividualWindows = true, std::cout << "Scan windows will be inferred separately for different runs\n";
        else if (!memcmp(&(args[start]), "convert ", 8)) Convert = true, std::cout << ".mzML to .dia conversion\n";
        else if (!memcmp(&(args[start]), "use-rt ", 7)) UseRTInfo = true, std::cout << "Existing .quant files will be used for RT profiling\n";
        else if (!memcmp(&(args[start]), "use-quant ", 10)) UseQuant = true, std::cout << "Existing .quant files will be used\n";
		else if (!memcmp(&(args[start]), "quant-only ", 11)) QuantOnly = true, std::cout << "Quantification will be performed anew using existing identification info\n";
		else if (!memcmp(&(args[start]), "report-only ", 12)) ReportOnly = true, std::cout << "Report will be generated using .quant files\n";
		else if (!memcmp(&(args[start]), "iter ", 5)) iN = Max(nRTWindowIter + 2, std::stoi(args.substr(start + 5, std::string::npos))),
			std::cout << "Number of iterations set to " << iN << "\n";
		else if (!memcmp(&(args[start]), "profiling-qvalue ", 17)) MaxProfilingQvalue = std::stod(args.substr(start + 17, std::string::npos)),
			std::cout << "RT profiling q-value threshold set to " << MaxProfilingQvalue << "\n";
		else if (!memcmp(&(args[start]), "quant-qvalue ", 13)) MaxQuantQvalue = std::stod(args.substr(start + 13, std::string::npos)),
			std::cout << "Q-value threshold for cross-run quantification set to " << MaxQuantQvalue << "\n";
		else if (!memcmp(&(args[start]), "rt-profiling ", 13)) RTProfiling = true, std::cout << "RT profiling turned on\n";
		else if (!memcmp(&(args[start]), "prefix ", 7)) prefix = trim(args.substr(start + 7, end - start - 7)); // prefix added to input file names
		else if (!memcmp(&(args[start]), "ext ", 4)) ext = trim(args.substr(start + 4, end - start - 4)); // extension added to input file names
		else if (!memcmp(&(args[start]), "no-test-dataset ", 16)) TestDataset = false, std::cout << "Library will not be split into training and test datasets\n";
		else if (!memcmp(&(args[start]), "test-proportion ", 16)) TestDatasetSize = std::stod(args.substr(start + 16, std::string::npos)),
			std::cout << "The " << TestDatasetSize << " fraction of the precursors will be used as the test dataset\n";
		else if (!memcmp(&(args[start]), "no-nn ", 6)) nnIter = iN, std::cout << "Neural network classifier turned off\n";
		else if (!memcmp(&(args[start]), "global-nn ", 10)) GlobalNN = GlobalTraining = true, std::cout << "Global NN training will be used\n";
		else if (!memcmp(&(args[start]), "nn-iter ", 8)) nnIter = Max(nRTWindowIter + 3, std::stoi(args.substr(start + 8, std::string::npos))),
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
		else if (!memcmp(&(args[start]), "mass-acc ", 9)) MassAccuracy = std::stod(args.substr(start + 9, std::string::npos)) / 1000000.0,
			std::cout << "Mass accuracy set to " << MassAccuracy << "\n";
		else if (!memcmp(&(args[start]), "gen-acc ", 8)) GeneratorAccuracy = std::stod(args.substr(start + 8, std::string::npos)) / 1000000.0,
			std::cout << "Fragmentation spectrum generator accuracy set to " << GeneratorAccuracy << "\n";
		else if (!memcmp(&(args[start]), "max-dppp ", 9)) MaxDpPP = Max(1, std::stoi(args.substr(start + 9, std::string::npos))),
			std::cout << "Maximum number of data points per peak set to " << MaxDpPP << "\n";
		else if (!memcmp(&(args[start]), "no-norm ", 8)) Normalisation = false, std::cout << "Cross-run normalisation turned off\n";
		else if (!memcmp(&(args[start]), "norm-qvalue ", 12)) NormalisationQvalue = std::stod(args.substr(start + 11, std::string::npos)), 
			std::cout << "Q-value threshold for cross-run normalisation set to " << NormalisationQvalue << "\n";
		else if (!memcmp(&(args[start]), "norm-fraction ", 14)) NormalisationPeptidesFraction = std::stod(args.substr(start + 14, std::string::npos)),
			std::cout << "Normalisation peptides fraction set to " << NormalisationPeptidesFraction << "\n";
		else std::cerr << "\nWARNING: unrecognised option [--" << trim(args.substr(start, end - start)) << "]\n\n";
        
        start = next;
    }
	if (prefix.length()) for (auto it = ms_files.begin(); it != ms_files.end(); it++) *it = prefix + *it;
	if (ext.length()) for (auto it = ms_files.begin(); it != ms_files.end(); it++) *it += ext;
	if (ms_files.size() < 2) RTProfiling = Normalisation = false;
    if (UseQuant || QuantOnly) UseRTInfo = true;
	if (GlobalNN) TestDataset = false;
	nnIter = Min(nnIter, iN);

	net.resize(nnBagging);
	netLock.resize(nnBagging);
	for (i = 0; i < nnBagging; i++) net[i].network = NULL, net[i].seed(i);
}

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

template <class Tx, class Ty> inline double corr(Tx * x, Ty * y, int n) {
	assert(n >= 1);

	double ex = 0.0, ey = 0.0, s2x = 0.0, s2y = 0.0, r = 0.0;
	for (int i = 0; i < n; i++) ex += x[i], ey += y[i];
	if (n) ex /= (double)n, ey /= (double)n;
	for (int i = 0; i < n; i++) s2x += Sqr(x[i] - ex), s2y += Sqr(y[i] - ey), r += (x[i] - ex) * (y[i] - ey);
	if (s2x < E || s2y < E) return 0.0;
	return r / sqrt(s2x * s2y);
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
template <class Tx, class Tfunc> std::vector<std::pair<Tx, Tfunc> > PAVA(std::vector<std::pair<Tx, Tfunc> > &data) {
	int i, j, k, l, n = data.size();
	std::vector<int> index(n + 1);
	std::vector<double> w(n), a(n);

	std::vector<std::pair<Tx, Tfunc> > result(n);
	for (i = 0; i < n; i++) result[i].first = data[i].first;

	a[0] = data[0].second;
	w[0] = 1.0;
	j = 0;
	index[0] = 0;
	index[1] = 1;

	for (i = 1; i < n; i++) {
		j++;
		a[j] = data[i].second;
		w[j] = 1.0;
		while (j > 0 && a[j] < a[j - 1]) {
			a[j - 1] = (w[j] * a[j] + w[j - 1] * a[j - 1])
				/ (w[j] + w[j - 1]);
			w[j - 1] = w[j] + w[j - 1];
			j--;
		}
		index[j + 1] = i + 1;
	}

	for (k = 0; k <= j; k++)
		for (l = index[k]; l < index[k + 1]; l++)
			result[l].second = a[k];

	return result;
}

template <class T, class Tx, class Tfunc> void
piecewise_linear(T * result, std::vector<Tx> &x, std::vector<std::pair<Tx, Tfunc> > &data, T min, T max) {
	int k, start, stop = 1, n = x.size(), m = data.size();

	assert(n >= 2);
	assert(m >= 2);

	result[0] = min;
	if (data[0].first > x[0]) for (stop = 1; stop < n && x[stop] <= data[0].first; stop++)
		result[stop] = result[0] + ((data[0].second - result[0]) * (x[stop] - x[0])) / (data[0].first - x[0]);
	for (start = stop - 1, k = 1; k < m; k++) {
		for (stop = start + 1; stop < n && x[stop] <= data[k].first; stop++)
			result[stop] = result[start] + ((data[k].second - result[start]) * (x[stop] - x[start])) / (data[k].first - x[start]);
		start = stop - 1;
	}
	if (start < n - 1) for (stop = start + 1; stop < n; stop++)
		result[stop] = result[start] + ((max - result[start]) * (x[stop] - x[start])) / (x[n - 1] - x[start]);

	return;
}

template <class T, class Tx, class Tfunc> void
smooth_isotonic_regression(T * result, std::vector<Tx> &x, std::vector<std::pair<Tx, Tfunc> > &data, T min, T max, int factor) {
	int i, j, start, stop, inc, n = x.size(), step = Max(8, n / factor);
	while (step % 8) step++;
	double s, w, u, mul = 1.0 / double(step);

	assert(n >= factor);
	assert(data.size() >= 2);

	std::vector<T> linear(n);
	std::vector<std::pair<Tx, Tfunc> > pava = PAVA(data);
	piecewise_linear(&(linear[0]), x, pava, min, max);

	for (i = 0; i < n; i++) {
		s = w = 0.0;
		start = Max(0, i - step);
		stop = Min(n - 1, i + step);
		inc = step / 8;
		for (j = start; j <= stop; j += inc) {
			u = exp(-2.0 * Sqr((double(j - i)) * mul));
			s += u * linear[j];
			w += u;
		}
		result[i] = s / w;
	}
}

class Parameter {
public:
    int min_iter_seek, min_iter_learn;

    Parameter(int _min_iter_seek, int _min_iter_learn) { min_iter_seek = _min_iter_seek, min_iter_learn = _min_iter_learn; }
};

const int nnS = 2, nnW = (2 * nnS) + 1, nnF = 6;
enum {
	pTimeCorr, pMinCorr, pCos, pMs1TimeCorr, pRT, pnRT,
	pRef,
	pSig = pRef + nnF,
	pCorr = pSig + nnF,
	pArray = pCorr + nnF,
	pN = pArray + nnF * nnW
};

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

const double proton = 1.007825035;
const double OH = 17.003288;

double AA[256];

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
}

inline int peptide_length(std::string &name) {
	int i, n = 0;
	for (i = 0; i < name.size(); i++) {
		char symbol = name[i];
		if (symbol >= 'A' && symbol <= 'Y') n++;
	}
	return n;
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

enum {
	type_b, type_y,
	type_N
};

enum {
	loss_none, loss_H2O, loss_NH3, loss_CO,
	loss_N
};

const double Loss[loss_N] = { 0.0, 18.011113035, 17.026548, 27.994915 };

class Ion {
public:
	char type, index, loss, charge; // index - number of AAs to the left from the cleavage site
	float mz;

	Ion(int _type, int _index, int _loss, int _charge, double _mz) {
		type = _type;
		index = _index;
		loss = _loss;
		charge = _charge;
		mz = _mz;
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

std::vector<Ion> recognise_fragments(std::vector<double> &sequence, std::vector<Peak> &fragments) {
	int i, j, cnt = 0, charge = 1, loss = loss_none;
	std::vector<Ion> result;
	std::vector<Ion> v;
	double delta, min, index = 0;

start:
	v = generate_fragments(sequence, charge, loss, &cnt);
	for (i = 0; i < fragments.size(); i++) {
		double margin = fragments[i].mz * GeneratorAccuracy;
		for (j = 0, min = margin; j < v.size(); j++) {
			delta = Abs(v[j].mz - fragments[i].mz);
			if (delta < min) min = delta, index = j;
		}
		if (min < margin) {
			result.push_back(v[index]);
			if (result.size() >= fragments.size()) goto stop;
		}
	}

	loss++;
	if (loss == loss_N) loss = loss_none, charge++;
	if (charge < 5 && result.size() < fragments.size()) goto start;

stop:
	if (result.size() < fragments.size())
		std::cerr << "Not all fragments recognised\n";
	return result;
}

std::vector<Ion> generate_fragments(std::vector<double> &sequence, std::vector<Ion> pattern) {
	int j, charge = 1, loss = loss_none, cnt = 0;
	std::vector<Ion> result;
	std::vector<Ion> v;
	auto pos = pattern.begin();

start:
	v = generate_fragments(sequence, charge, loss, &cnt);
	for (; pos != pattern.end(); pos++) {
		for (j = 0; j < v.size(); j++) if (v[j] == *pos) {
			result.push_back(v[j]);
			break;
		}
		if (j == v.size()) break;
	}
	loss++;
	if (loss == loss_N) loss = loss_none, charge++;
	if (charge < 4 && pos != pattern.end()) goto start;

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
    
	template <bool get_mz> inline double level(float mz, float * peak_mz = NULL) {
		int i, low = 0, high = size();
		float v, s, margin, min, max;
		margin = mz * MassAccuracy;
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
	int index;
	mutable bool decoy_found;
	mutable int apex, peak, best_fragment, peak_width;
	mutable float RT, predicted_nRT, qvalue, quantity, pg_quantity;

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

class Library {
public:
	std::string name;
	std::vector<Protein> protein_list; // proteins
	std::vector<int> reference; // indices of reference peptides used for initial nRT window determination
	bool gen_decoys = true;
	float nRT_min = INF, nRT_max = -INF;

	class Entry {
	public:
		Library * lib;
		Peptide target, decoy;
		bool exclude = false;
		std::string name; // precursor id
		std::string pep_name; // modified peptide id
		std::string pg_name; // protein group id
		int pg_index;

		void init() {
			std::sort(target.fragments.begin(), target.fragments.end(), [](const Peak &left, const Peak &right) { return left.height > right.height; });
			std::sort(decoy.fragments.begin(), decoy.fragments.end(), [](const Peak &left, const Peak &right) { return left.height > right.height; });
			if (target.fragments.size() > nnF) target.fragments.resize(nnF);
			if (decoy.fragments.size() > nnF) decoy.fragments.resize(nnF);

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

			auto seq = get_sequence(name);

			auto pattern = recognise_fragments(seq, target.fragments);
			auto inv_seq = seq;

			for (i = 1; i < seq.size() - 1; i++) inv_seq[i] = seq[seq.size() - i - 1];
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

				for (auto jt = (*v).begin(); jt != (*v).end(); jt++)
					for (i = 0, jt->second.pr.quantity = 0.0; i < m; i++)
						if (score[i] >= margin) jt->second.pr.quantity += jt->second.fr[i].quantity[qFiltered];
			}

			if (!Normalisation) return;

			std::vector<float> score(map.size());
			std::vector<float> x(ms_files.size());
			i = m = 0;
			for (auto it = map.begin(); it != map.end(); it++, i++) {
				auto v = &(it->second);

				k = 0;
				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					int index = jt->first;
					if (jt->second.pr.qvalue < NormalisationQvalue) x[index] = jt->second.pr.quantity, k++, m++;
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
				std::cout << "\nWARNING: not enough peptides for normalisation\n";
				return;
			}

			auto ordered = score;
			std::sort(ordered.begin(), ordered.end());
			int average_number = m / ms_files.size();
			int used = Max(1, int(NormalisationPeptidesFraction * (double)average_number));
			double margin = ordered[used] + E;

			std::vector<float> sums(ms_files.size());
			i = 0;
			for (auto it = map.begin(); it != map.end(); it++, i++) {
				auto v = &(it->second);
				if (score[i] > margin) continue;

				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					int index = jt->first;
					if (jt->second.pr.qvalue <= NormalisationQvalue) sums[index] += jt->second.pr.quantity;
				}
			}

			double av = 0.0;
			for (i = k = 0; i < sums.size(); i++) if (sums[i] > E) av += sums[i], k++;
			if (k) av /= (double)k;
			if (av < E) {
				if (Verbose >= 1) std::cout << "\nWARNING: Cannot perform normalisation\n\n";
				return;
			}
			for (i = 0; i < sums.size(); i++) if (sums[i] <= E) sums[i] = av;
			for (auto it = map.begin(); it != map.end(); it++) {
				auto v = &(it->second);

				for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
					int index = jt->first;
					jt->second.pr.quantity *= av / sums[index];
				}
			}
		}
	};

	Info info;

	bool load(const char * file_name) {
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
				for (it = words.begin(); it != words.end(); it++)
					for (auto jt = library_headers.begin(); jt != library_headers.end(); jt++)
						if (jt->find(*it) != std::string::npos) {
							colInd[std::distance(library_headers.begin(), jt)] = std::distance(words.begin(), it);
							break;
						}
				for (i = 0; i < libCols; i++) if (colInd[i] < 0 && i < libIsDecoy)
					Throw("Cannot find column " + library_headers[i]);
				continue;
			}

			Peak p(std::stof(words[colInd[libFrMz]]), std::stof(words[colInd[libFrI]]));

			auto ins = map.insert(std::pair<std::string, Entry>(words[colInd[libPr]] + words[colInd[libCharge]], e));
			ins.first->second.name = ins.first->first;
			ins.first->second.pg_name = words[colInd[libPG]];
			ins.first->second.pep_name = words[colInd[libPr]];

			bool decoy_fragment = colInd[libIsDecoy] < 0 ? false : std::stoi(words[colInd[libIsDecoy]]);
			if (decoy_fragment) gen_decoys = false;
			auto pep = !decoy_fragment ? (&(ins.first->second.target)) : (&(ins.first->second.decoy));
			pep->mz = std::stof(words[colInd[libPrMz]]);
			pep->nRT = std::stof(words[colInd[libnRT]]);
			pep->charge = std::stoi(words[colInd[libCharge]]);
			pep->fragments.push_back(p);
		}

		csv.close();

		entries.resize(map.size());
		i = 0;
		for (auto it = map.begin(); it != map.end(); it++, i++) {
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

		if (!ref_file.size()) {
			std::ifstream in(name + std::string(".ref"), std::ifstream::binary);
			if (in.is_open()) {
				int size; in.read((char*)&size, sizeof(int));
				reference.resize(size);
				in.read((char*)&(reference[0]), size * sizeof(int));
				in.close();
			}
		}

		return true;
	}

	void generate_decoys() {
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (gen_decoys) it->generate_decoy();
			it->init();
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

		std::ofstream out(file_name, std::ifstream::out);
		out << oh[outFile] << "," << oh[outPG] << "," << oh[outPGQ] << "," << oh[outModSeq] << "," 
			<< oh[outPrId] << "," << oh[outCharge] << "," << oh[outQv] << "," << oh[outPrQ] << "," << oh[outRT] << "\n";

		for (auto it = info.map.begin(); it != info.map.end(); it++) {
			auto v = &(it->second);
			auto entry = &(entries[it->first]);

			for (auto jt = (*v).begin(); jt != (*v).end(); jt++) {
				out << ms_files[jt->first] << ","
					<< entry->pg_name << ","
					<< jt->second.pr.pg_quantity << ","
					<< entry->pep_name << ","
					<< entry->name << ","
					<< entry->target.charge << ","
					<< jt->second.pr.qvalue << ","
					<< jt->second.pr.quantity << ","
					<< jt->second.pr.RT << "\n";
			}
		}

		out.close();
		if (Verbose >= 1) std::cout << "Report saved to " << file_name << ".\n";
	}
};

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
std::vector<std::pair<float, float> > rt_data;
std::vector<float> rt_delta;

class Run {
public:
    std::string name; // run name
    float weights[pN];
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

    int curr_iter;
    double min_target_decoy_score = 0.0, nRT_cscore, nRT_ref_score, PeakWidth = 0.0;
    
    std::vector<Parameter> pars;
    bool par_seek[pN], par_learn[pN];

	bool standardised = false;

    Run() {
		pars.push_back(Parameter(0, 0)); // pTimeCorr
		pars.push_back(Parameter(1, 0)); // pMinCorr
        pars.push_back(Parameter(1, 0)); // pCos
        pars.push_back(Parameter(1, 0)); // pMs1TimeCorr
		pars.push_back(Parameter(nRTWindowIter + 1, nRTWindowIter + 1)); // pRT
		pars.push_back(Parameter(nRTWindowIter + 1, nRTWindowIter + 1)); // pnRT
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
		PeakList.resize(Threads), BestFrList.resize(Threads);
    }
    
    template <bool seek> inline double combined_score(float * scores) {
		if (curr_iter < nnIter + (int)seek || GlobalTraining) {
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

			float Mz, *ms1, *ms2, *ms2_min;

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

			void chromatogram(int from, int to, bool get_ms1 = false) { // [from, to)
				int i, k, fr, l = to - from, pos, ind;
				float ms1_left = 0.0, ms1_right = 0.0, rt, ms1_left_RT = 0.0, ms1_right_RT = 0.0;

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

					if (get_ms1 && run->curr_iter > nRTWindowIter && UseRefFile) { // main search phase only
						double margin = nRTWindowMargin * run->nRT_window;
						if (run->predicted_nRT[pr->scan_index[Min(to - 1, k + pr->W)]] < pr->nRT - margin) continue;
						if (run->predicted_nRT[pr->scan_index[Max(from, k - pr->W)]] > pr->nRT + margin) break;
					}

					if (get_ms1) {
						rt = run->scan_RT[i];
						if (rt > ms1_right_RT) {
							auto ms1_ptr = std::lower_bound(run->ms1_RT.begin(), run->ms1_RT.end(), rt);
							pos = std::distance(run->ms1_RT.begin(), ms1_ptr);
							if (ms1_ptr == run->ms1_RT.end()) ms1_right_RT = INF, ms1[k] = ms1_right;
							else {
								ms1_left = ms1_right, ms1_left_RT = ms1_right_RT;
								ms1_right_RT = *ms1_ptr;
								ms1_right = run->ms1[pos].level<false>(Mz);
							}
						}
						ms1[ind] = ((rt - ms1_left_RT) * ms1_right + (ms1_right_RT - rt) * ms1_left) / Max(E, ms1_right_RT - ms1_left_RT);
					}

					for (fr = 0; fr < m; fr++)
						ms2[l * fr + ind] = run->scans[i].level<false>(mz[fr]);
				}

				for (fr = 0; fr < m; fr++)
					for (k = 1; k < l - 1; k++)
						ms2_min[l * fr + k] = Min(Min(ms2[l * fr + k - 1], ms2[l * fr + k + 1]), ms2[l * fr + k]);
			}

			void peaks() {
				int i, k, fr, pos, next, best_fr = -1, l = pr->scan_number, n_peaks = 0;
				double max, s;
				double * corr_matrix = (double*)alloca(m * m * sizeof(double));
				float * elution = (float*)alloca(pr->W * sizeof(float));

				run->PeakList[pr->thread_id].resize(l);
				run->BestFrList[pr->thread_id].resize(l);
				int * peak_list = &(run->PeakList[pr->thread_id][0]);
				int * best_fr_list = &(run->BestFrList[pr->thread_id][0]);

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

					peak_list[n_peaks] = k, best_fr_list[n_peaks++] = best_fr;
				}

				pr->peaks.resize(pr->peak_number = n_peaks);
				pr->scoring.resize(n_peaks * pN, 0.0);
				for (i = 0; i < pr->peak_number; i++) {
					pr->peaks[i].peak = peak_list[i];
					pr->peaks[i].apex = pr->scan_index[peak_list[i]];
					pr->peaks[i].best_fragment = best_fr_list[i];
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

		void search() {
			if (!initialised) {
				Searcher searcher(this);
				build_index();
				searcher.chromatogram(0, scan_number, true);
				searcher.peaks();
				for (int peak = 0; peak < peak_number; peak++) searcher.score(peak);
				initialised = true;
			}
			for (int peak = 0; peak < peak_number; peak++) score_RT(peak);
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

			int low = Max(0, k - Min(MaxDpPP / 2 - (1 ^ (MaxDpPP & 1)), S / 2));
			int high = Min(scan_number - 1, k + Min(MaxDpPP / 2, S / 2));
			int len = high - low + 1;
			float *elution = (float*)alloca(len * sizeof(float));
			float *signal = (float*)alloca(len * sizeof(float));

			searcher.chromatogram(low, high + 1, false);
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
        }

		void combine_scores(int _thread_id) {
			int i, j;
			if (run->curr_iter <= nnIter || GlobalTraining)
				for (i = 0; i < peaks.size(); i++) 
					peaks[i].score = scalar(scores, run->weights, run->par_learn, pN);
			else {
				Matrix* mat = createMatrix(peaks.size(), pN, &(scoring[0]));
				for (i = 0; i * Threads + _thread_id < nnBagging; i++) {
					forwardPass(net[i * Threads + _thread_id].network, mat);
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
					if (UseRefFile && run->curr_iter > nRTWindowIter 
						&& Abs(run->predicted_nRT[peaks[peak].apex] - nRT) > nRTWindowMargin * run->nRT_window) continue;
                    best_peak = peak;
                    combined_score = score;
                }
            }
            found = (best_peak >= 0);
			if (found) {
				best_fragment = peaks[best_peak].best_fragment;
				apex = peaks[best_peak].apex;
				predicted_nRT = run->predicted_nRT[apex];
				RT = run->scan_RT[apex];
				for (k = 0; k < pN; k++) scores[k] = scoring[best_peak * pN + k];
				if (!decoy) if (run->curr_iter == iN - 1 || (run->curr_iter == nRTWindowIter && InferWindow)) quantify(peaks[best_peak].peak);
			}
			
        }

		void free() {
			initialised = false;
			std::vector<int>().swap(scan_index);
			std::vector<Elution>().swap(peaks);
		}
        
        void find(int _thread_id) {
			if (lock.set()) {
				thread_id = _thread_id;
				if (!initialised) quant.resize(fragments.size());
				if (!QuantOnly) {
					search();
					seek();
				} else quantify(peak_pos);
			}
        }
    };

	class Target {
	public:
		Precursor target;
		Precursor decoy;
	};

	std::vector<Target> entries;
	std::vector<int> reference;

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
		for (i = 1; i < n_scans; i++) { // handle overlapping windows
			if (scans[i - 1].window_low < scans[i].window_low && scans[i - 1].window_high > scans[i].window_low) {
				float middle = (scans[i - 1].window_high + scans[i].window_low) * 0.5;
				scans[i].window_low = scans[i - 1].window_high = middle;
			}
		}
		predicted_nRT.resize(n_scans), scan_RT.resize(n_scans);
		for (i = 0; i < n_scans; i++) predicted_nRT[i] = 0.0, scan_RT[i] = scans[i].RT;

		ms1_RT.resize(ms1.size());
		for (i = 0; i < ms1_RT.size(); i++) ms1_RT[i] = ms1[i].RT;
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
			for (i = 0; i < max; i++) if (scans[i].has(mz)) break;
			if (i < max) has[pos] = true, cnt++;
		}
		entries.resize(cnt);
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

		reference.clear();
		for (i = 0; i < lib->reference.size(); i++) 
			if (has[lib->reference[i]]) reference.push_back(lib->reference[i]);

		// select precursors that will not be used for classifier training
		cnt = 0;
		if (TestDataset) {
			std::mt19937_64 gen(1);
			for (pos = ls = 0; pos < entries.size(); pos++) if (!entries[pos].target.exclude) ls++;
			int max = Max(1, Min(ls / 2, (int)(TestDatasetSize * (double)ls)));
			while (cnt < max) {
				pos = gen() % (unsigned long long)ls;
				if (entries[pos].target.exclude) continue;
				if (!entries[pos].target.test) entries[pos].target.test = true, cnt++;
			}
		}
		for (pos = 0; pos < entries.size(); pos++) 
			if (!entries[pos].target.exclude && !entries[pos].target.test) entries[pos].target.training = true;
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
    
    void process_precursors(int thread_id) {
		if (curr_iter <= nRTWindowIter && reference.size() && UseRefFile && !QuantOnly) {
			for (int i = 0; i < reference.size(); i++) {
				auto it = find_entry(reference[i]);
				if (it == NULL) continue;

				it->target.find(thread_id);
				it->decoy.find(thread_id);
			}
			return;
		}

        for (int i = 0; i < entries.size(); i++) {  
            if (!QuantOnly || entries[i].target.found) entries[i].target.find(thread_id);
			if (!QuantOnly) entries[i].decoy.find(thread_id);
        }
    }

	void combine_scores(int thread_id) {
		for (int i = 0; i < entries.size(); i++) {
			if (entries[i].target.found) entries[i].target.combine_scores(thread_id);
			if (entries[i].decoy.found) entries[i].decoy.combine_scores(thread_id);
		}
	}
    
    void seek_precursors() {
		if (Verbose >= 1) std::cout << "Precursor search...\n";
        
        int i;
        if (Threads > 1) {
            std::vector<std::thread> threads;
            for (i = 0; i < Threads; i++) threads.push_back(std::thread(&Run::process_precursors, this, i));
            for (i = 0; i < Threads; i++) threads[i].join();
        } else process_precursors(0);
		for (i = 0; i < entries.size(); i++) {
			entries[i].target.lock.free(), entries[i].decoy.lock.free();
			if (curr_iter == nRTWindowIter && InferWindow) entries[i].target.free(), entries[i].decoy.free();
		}
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

		if (Verbose >= 1) std::cout << "Optimising weights...\n";
        
		pos = 0;
        for (i = 0; i < pN; i++) {
            av[i] = 0.0;
            for (j = 0; j < pN; j++) A(i, j) = 0.0;
        }
		int nnCnt = 0, nnTest = 0;
        for (auto it = entries.begin(); it != entries.end(); it++) {
			if (!GlobalNN) {
				if (curr_iter >= nnIter && it->target.training) {
					if (it->target.found) {
						for (i = 0; i < nnBagging; i++) training[i][nnCnt] = it->target.scores, training_classes[i][nnCnt] = target_nn;
						nnCnt++;
					}
					if (it->decoy.found) {
						for (i = 0; i < nnBagging; i++) training[i][nnCnt] = it->decoy.scores, training_classes[i][nnCnt] = decoy_nn;
						nnCnt++;
					}
				}

				if (curr_iter >= nnIter && Verbose >= 3 && TestDataset && it->target.test) {
					if (it->target.found) test[nnTest] = it->target.scores, test_classes[nnTest++] = target_nn;
					if (it->decoy.found) test[nnTest] = it->decoy.scores, test_classes[nnTest++] = decoy_nn;
				}
			}

            if (!it->target.training) continue;
            if (!it->target.found || !it->decoy.found) continue;

            for (i = 0; i < pN; i++) {
                double x = par_learn[i] ? it->target.scores[i] : 0.0;
				double y = par_learn[i] ? it->decoy.scores[i] : 0.0;
                av[i] += (delta[pos * pN + i] = x - y);
            }

            pos++;
        }
		if (!pos) {
			std::cout << "\nWARNING: no training precursors\n";
			return;
		}

        for (i = 0; i < pN; i++) av[i] /= (double)pos;
		pos = 0;
		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (!it->target.training) continue;
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
			net[i].learningRate = nnLearning;
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
    
    void calculate_qvalues() {
		if (Verbose >= 2) std::cout << "Calculating q-values...\n";
        
        int i = 0, ids = 0;
        std::set<std::pair<float, int> > decoy_scores;
        std::set<std::pair<float, int> > target_scores;
        min_target_decoy_score = INF;
        nRT_cscore = nRT_ref_score = INF;
        
        for (auto it = entries.begin(); it != entries.end(); it++, i++) {
			it->target.combined_score = combined_score<false>(it->target.scores);
			it->decoy.combined_score = combined_score<false>(it->decoy.scores);
            
            if (it->target.found && it->target.combined_score < min_target_decoy_score) 
                min_target_decoy_score = it->target.combined_score;
            if (it->decoy.found && it->decoy.combined_score < min_target_decoy_score) 
                min_target_decoy_score = it->decoy.combined_score;
            
            if (!TestDataset || it->target.test) {
                if (it->target.found) target_scores.insert(std::pair<float, int>(-it->target.combined_score, i));
                if (it->decoy.found) decoy_scores.insert(std::pair<float, int>(-it->decoy.combined_score, i));
            }
        }
        
        for (auto it = entries.begin(); it != entries.end(); it++) {
            if (!it->target.found) it->target.combined_score = min_target_decoy_score;
            if (!it->decoy.found) it->decoy.combined_score = min_target_decoy_score;
        }

        std::vector<std::pair<float, int> > decoys(decoy_scores.begin(), decoy_scores.end());
        std::vector<std::pair<float, int> > targets(target_scores.begin(), target_scores.end());
        
		std::map<float, float> cs_qv;

        for (auto it = entries.begin(); it != entries.end(); it++) {
            if (!it->target.found) it->target.qvalue = 1.0;
            else {
                auto pD = std::upper_bound(decoys.begin(), decoys.end(), 
                        std::pair<float, int>(-it->target.combined_score, 0));
                auto pT = std::upper_bound(targets.begin(), targets.end(), 
                        std::pair<float, int>(-it->target.combined_score, 0));
				int n_targets = std::distance(targets.begin(), pT);
                it->target.qvalue = Min(1.0, (double(std::distance(decoys.begin(), pD))) / (double) Max(1, n_targets));
				if (n_targets <= nRTTopPrecursors && it->target.combined_score < nRT_cscore) nRT_cscore = it->target.combined_score;
				if (n_targets <= nRTRefPrecursors && it->target.combined_score < nRT_ref_score) nRT_ref_score = it->target.combined_score;

				auto pair = std::pair<float, float>(it->target.combined_score, it->target.qvalue);
				auto pos = cs_qv.insert(pair);
				if (pos.second) {
					if (pos.first != cs_qv.begin() && std::prev(pos.first)->second < pair.second) pair.second = std::prev(pos.first)->second;
					else for (auto jt = std::next(pos.first); jt != cs_qv.end() && jt->second > pair.second; jt++) jt->second = pair.second;
				}
            }
        }

		for (auto it = entries.begin(); it != entries.end(); it++) {
			if (!it->target.found) continue;
			if (it->target.qvalue < 1.0) {
				auto pos = cs_qv.lower_bound(it->target.combined_score);
				it->target.qvalue = pos->second;
			}
			if (!TestDataset || it->target.test) if (it->target.qvalue <= 0.01) ids++;
		}
		if (Verbose >= 2) std::cout << "Number of test set IDs at 0.01 FDR: " << ids << "\n";
    }
    
    void estimate_nRT() {
		if (Verbose >= 1) std::cout << "Estimating nRT, score threshold = " << nRT_cscore << "...\n";

		int peak_width = 0, peak_cnt = 0;
		bool gen_ref = false;
		if (curr_iter == iN - 1 && ((!reference.size() && UseRefFile) || UpdateRefFile)) gen_ref = true;

		rt_stats.clear();  
		rt_stats.reserve(entries.size());
		rt_ref.clear();
		rt_ref.reserve(entries.size());
        
        for (auto it = entries.begin(); it != entries.end(); it++) {
            if (!it->target.found) continue;
			if (it->target.qvalue <= nRTMaxQvalue || it->target.combined_score >= nRT_cscore) {
				int index = std::distance(entries.begin(), it);
				rt_stats.push_back(index);
				if (curr_iter == nRTWindowIter && it->target.qvalue < nRTMaxQvalue && it->target.combined_score >= nRT_cscore)
					peak_width += it->target.peak_width, peak_cnt++;
				if (gen_ref && it->target.combined_score >= nRT_ref_score) rt_ref.push_back(it->target.index);
 			}
        }
		std::sort(rt_stats.begin(), rt_stats.end(), [&](const auto& lhs, const auto& rhs) { return entries[lhs].target.RT < entries[rhs].target.RT; });
        rt_data.resize(rt_stats.size());
		for (int i = 0; i < rt_stats.size(); i++) rt_data[i] = std::pair<float, float>(entries[rt_stats[i]].target.RT, entries[rt_stats[i]].target.nRT);
		if (Verbose >= 2) std::cout << rt_data.size() << " precursors used for nRT estimation.\n";

		if (rt_data.size() >= 2)
			smooth_isotonic_regression(&(predicted_nRT[0]), scan_RT, rt_data, lib->nRT_min, lib->nRT_max, Min(scan_RT.size(), nRTSmoothingFactor));
		else return;

		if (curr_iter == nRTWindowIter) {
			if (nRTWindow > E) nRT_window = nRTWindow;
			else {
				if (Verbose >= 1) std::cout << "Calculating nRT windows...\n";
				rt_delta.resize(rt_stats.size());
				for (int i = 0; i < rt_stats.size(); i++) rt_delta[i] = Abs(entries[rt_stats[i]].target.nRT - entries[rt_stats[i]].target.predicted_nRT);
				std::sort(rt_delta.begin(), rt_delta.end());
				nRT_window = Max(-rt_delta[(int)(nRTWindowLoss * (double)rt_delta.size())], (lib->nRT_max - lib->nRT_min) / (double)nRTWindowFactor);
			}
		}

		if (gen_ref) {
			std::sort(rt_ref.begin(), rt_ref.end());
			std::ofstream out(lib->name + std::string(".ref"), std::ofstream::binary);
			int size = rt_ref.size();
			out.write((char*)&size, sizeof(int));
			out.write((char*)&(rt_ref[0]), size * sizeof(int));
			out.close();
		}

		if (curr_iter == nRTWindowIter && InferWindow) {
			PeakWidth = (double(peak_width)) / (double)Max(1, peak_cnt);
			if (Verbose >= 1) std::cout << "Peak width: " << PeakWidth << "\n";

			if (!IndividualWindows) WindowRadius = Max(1, int(ScanScale * PeakWidth)), InferWindow = false;
			if (Verbose >= 1) std::cout << "Scan window radius set to " << WindowRadius << "\n";
		}
    }
    
    void process() {
        int i;

        for (curr_iter = 0; curr_iter < iN; curr_iter++) {
			if (Verbose >= 2) std::cout << "\n[ iteration " << curr_iter << " ]\n";
            for (i = 0; i < pN; i++) {
                par_seek[i] = (curr_iter >= pars[i].min_iter_seek);
                par_learn[i] = (curr_iter >= pars[i].min_iter_learn);
            }
            
			if (curr_iter > nnIter && !GlobalTraining && !QuantOnly) {
				std::vector<std::thread> thr;
				for (i = 0; i < Threads; i++) thr.push_back(std::thread(&Run::combine_scores, this, i));
				for (i = 0; i < Threads; i++) thr[i].join();
			}
            seek_precursors();
			if (!QuantOnly) {
				if (curr_iter >= nnIter && nnIter < iN && nnStandardise) standardise_scores();
				fit_weights();
				calculate_qvalues();
			}
            estimate_nRT();
			if (QuantOnly) break;
        }

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
			qe.pr.peak = !QuantOnly ? it->target.peaks[it->target.best_peak].peak : it->target.peak_pos;
			qe.pr.best_fragment = !QuantOnly ? it->target.peaks[it->target.best_peak].best_fragment : it->target.best_fragment;
			qe.pr.peak_width = it->target.peak_width;
			qe.pr.index = it->target.index;
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
			Run run;
			if (!run.load(&((*it)[0]))) exit(-1);
			std::string(dia_name) = run.name + std::string(".dia");
			run.write(&(dia_name[0]));
		}
		return 0;
	}

	Library lib;
	if (!lib.load(&(lib_file[0]))) Throw("Cannot load spectral library");
	if (ref_file.size()) {
		Library ref;
		if (!ref.load(&(ref_file[0]))) std::cerr << "\nWARNING: Cannot load reference peptides\n\n", UseRefFile = false;
		else {
			int target_size = lib.entries.size();
			lib.entries.insert(lib.entries.end(), ref.entries.begin(), ref.entries.end());
			lib.reference.resize(ref.entries.size());
			std::string ref_prot_name = std::string("ReferencePG"), ref_pr_name = std::string("ref_pr_");
			lib.protein_list.push_back(Protein(ref_prot_name));
			for (int i = 0; i < lib.reference.size(); i++) {
				lib.reference[i] = lib.entries[target_size + i].target.index = target_size + i;
				lib.entries[target_size + i].name = ref_pr_name + lib.entries[target_size + i].name;
				lib.entries[target_size + i].pep_name = ref_pr_name + lib.entries[target_size + i].pep_name;
				lib.entries[target_size + i].pg_name = ref_prot_name;
				lib.entries[target_size + i].pg_index = lib.protein_list.size() - 1;
				lib.protein_list[lib.protein_list.size() - 1].precursors.push_back(target_size + i);
			}
		}
	}
	lib.generate_decoys();

	if (ReportOnly) goto cross_run;
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
        
        Run run;
        if (!run.load(&((*it)[0]))) exit(-1);
		run.load_library(&lib);
        run.process();

		if (UpdateRefFile) UpdateRefFile = false, UseRefFile = true;
    }
	if (GlobalTraining) {
		global_training(&lib);
		GlobalTraining = false;
		goto first_loop;
	}
	if (!RTProfiling) goto cross_run;

	quant_only:
	if (QuantOnly) {
		if (Verbose >= 1) std::cout << "Quantification...\n";
		lib.info.load(&(lib), ms_files);
		int runN = 0;
		for (auto it = ms_files.begin(); it != ms_files.end(); it++, runN++) {
			Run run;
			if (!run.load(&((*it)[0]))) exit(-1);
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
			if (jt->index >= 0) if (jt->pr.qvalue < MaxProfilingQvalue)
				lib.entries[jt->index].target.nRT = lib.entries[jt->index].decoy.nRT = jt->pr.predicted_nRT;

        for (auto it = ms_files.begin(); it != ms_files.end(); it++) {
			if (UseQuant) {
				std::ifstream in((*it) + std::string(".quant"));
				if (in.is_open()) {
					in.close();
					continue;
				}
			}
			if (Verbose >= 1) std::cout << "Second pass: file #" << std::distance(ms_files.begin(), it) << "\n";

			Run run;
			if (!run.load(&((*it)[0]))) exit(-1);
			run.load_library(&lib);
			run.process();

			if (UpdateRefFile) UpdateRefFile = false, UseRefFile = true;
		}
	}
    
cross_run:
	if (Verbose >= 1) std::cout << "Cross-run analysis...\n";
	lib.info.load(&(lib), ms_files);
	lib.info.quantify();
	lib.quantify_proteins(3, 0.01);
	lib.report(out_file);

    return 0;
}


