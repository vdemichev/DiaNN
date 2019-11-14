/*
Copyright 2019, Vadim Demichev

This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license,
visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
*/

#define _HAS_ITERATOR_DEBUGGING 0
#define _ITERATOR_DEBUG_LEVEL 0  
#define _CRT_SECURE_NO_WARNINGS

#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4305)

const double E = 0.000000001;
const double MinPeakHeight = 0.01;
#define Min(x, y) ((x) <= (y) ? (x) : (y))
#define Max(x, y) ((x) >= (y) ? (x) : (y))

#include <iostream>
#include <windows.h>
#include <vector>
#include <algorithm>
#include <atomic>

#using "../Clearcore2.Data.dll"
#using "../Clearcore2.Data.AnalystDataProvider.dll"
#using "../Sciex.Data.XYData.dll"
#using "../Clearcore2.RawXYProcessing.dll" // for vendor centroiding

using namespace Clearcore2::Data;
using namespace Clearcore2::Data::AnalystDataProvider;
using namespace Clearcore2::Data::DataAccess;
using namespace Clearcore2::Data::DataAccess::SampleData;
using namespace System::Threading;

class Lock {
public:
	std::atomic_bool lock = ATOMIC_VAR_INIT(false);
	Lock() :lock() {}

	Lock(const std::atomic_bool &flag) :lock(flag.load()) {}

	Lock(const Lock &another) :lock(another.lock.load()) {}

	Lock &operator=(const Lock &another) { lock.store(another.lock.load()); return *this; }

	bool set() { return !lock.exchange(true); }
	void free() { lock.store(false); }
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

const double qspread = 5.0;
__forceinline double quant_func(double t) {
	double t2 = t * t;
	return (1.0 - t2) * exp(-0.5 * t2);
}

class Spectrum {
public:
	int MS_level = 0;
	double RT = 0.0, window_low = 0.0, window_high = 0.0;
	std::vector<Peak> peaks;

	Spectrum() {}
	void init(Spectrum &S) { RT = S.RT, window_low = S.window_low, window_high = S.window_high, peaks = std::move(S.peaks); }
	void init(double * x, double * y, int n, bool centroid, double scan_width, double quant_width, int level, double _RT, double _low, double _high) {
		peaks.clear();
		MS_level = level, RT = _RT, window_low = _low, window_high = _high;

		if (!n) return;
		if (!centroid) {
			int i, m;
			for (i = m = 0; i < n; i++) if (y[i] >= MinPeakHeight) m++;
			peaks.resize(m);
			for (i = m = 0; i < n; i++) if (y[i] >= MinPeakHeight) {
				peaks[m].mz = x[i];
				peaks[m++].height = y[i];
			}
			return;
		}

		int i, back, forward;
		bool up = true;
		double sum = 0.0, win = x[0] * scan_width, min = x[0] - win, max = x[0] + win, low, high;
		double qwin = x[0] * quant_width, qmin = x[0] - qwin, qmax = x[0] + qwin, qlow, qhigh, qscale = 1.0 / qwin;
		for (i = 0; i < n; i++) {
			if (x[i] >= max) break;
			sum += y[i];
		}
		back = 0, forward = i;

		for (i = 1; i < n; i++) {
			low = min, high = max, qlow = qmin, qhigh = qmax;
			win = x[i] * scan_width, qwin = x[i] * quant_width * qspread;
			min = x[i] - win, max = x[i] + win, qmin = x[i] - qwin, qmax = x[i] + qwin;

			double delta = 0.0;
			while (true) {
				if (x[back] > min) break;
				delta -= y[back++];
			}
			while (forward < n) {
				if (x[forward] >= max) break;
				delta += y[forward++];
			}

			if (delta <= 0.0) {
				if (up) {
					sum = 0.0;
					for (int j = i; j < n; j++) {
						if (x[j] >= qhigh) break;
						sum += y[j] * quant_func((x[j] - x[i - 1]) * qscale);
					}
					for (int j = i - 1; j >= 0; j--) {
						if (x[j] <= qlow) break;
						sum += y[j] * quant_func((x[j] - x[i - 1]) * qscale);
					}
					if (sum >= MinPeakHeight) peaks.push_back(Peak(x[i - 1], sum));
				}
				up = false;
			} else {
				up = true;
				sum = Max(0.0, sum + delta);
			}
			qscale = 1.0 / (x[i] * quant_width);
		}
	}

	int * put(int * to) {
		int * iptr = to;
		(*iptr++) = peaks.size();
		(*iptr++) = MS_level;
		double * dptr = (double*)iptr;
		(*dptr++) = RT; (*dptr++) = window_low; (*dptr++) = window_high;
		float * fptr = (float*)dptr;
		for (int i = 0; i < peaks.size(); i++) (*fptr++) = peaks[i].mz, (*fptr++) = peaks[i].height;
		return (int*)fptr;
	}

	inline long long mem_size() { return 2 * sizeof(int) + 3 * sizeof(double) + 2 * sizeof(float) * (long long)peaks.size(); }
};

public ref class MTLoad {
public:
	std::vector<Spectrum> * spectra;
	std::vector<Lock> * locks;
	char * file;
	bool vendor = true;

	MTLoad(char * _file, bool _vendor, std::vector<Spectrum> * _spectra, std::vector<Lock> * _locks) {
		file = _file;
		spectra = _spectra;
		locks = _locks;
		vendor = _vendor;
	}

	void load_spectra() {
		try {
			int pos = 0, e = 0, ns = 0, s = 0, m;
			auto provider = gcnew AnalystWiffDataProvider();
			while (!(*locks)[0].set()) {}
			auto batch = AnalystDataProviderFactory::CreateBatch(gcnew System::String(file), provider);
			auto sample = batch->GetSample(0);
			auto mss = sample->MassSpectrometerSample;
			int en = mss->ExperimentCount;
			(*locks)[0].free();
			
			try {
				std::vector<double> X, Y;
				for (e = 0; e < en; e++) {
					auto exp = mss->GetMSExperiment(e);
					int n = exp->Details->NumberOfScans;
					double res = Max(10000.0, exp->Details->DefaultResolution);
					double wl = 0.0, wh = 0.0, win, centre;
					if (exp->Details->ExperimentType == (ExperimentType)1 && exp->Details->MassRangeInfo->Length > 0) try {
						win = (double)((FragmentBasedScanMassRange^)(exp->Details->MassRangeInfo[0]))->IsolationWindow;
						centre = (double)((FragmentBasedScanMassRange^)(exp->Details->MassRangeInfo[0]))->FixedMasses[0];
					} catch (std::exception &exc) {}
					for (s = 0; s < n; s++) {
						if ((*locks)[pos].set()) {
							auto spi = exp->GetMassSpectrumInfo(s);
							wl = centre - 0.5 * win, wh = centre + 0.5 * win;
							if (!vendor) {
								try {
									auto sp = exp->GetMassSpectrum(s);
									auto x = sp->GetActualXValues();
									auto y = sp->GetActualYValues();
									m = x->Length;
									X.resize(m), Y.resize(m);
									for (int i = 0; i < m; i++) X[i] = x[i], Y[i] = y[i];
									if (spi->MSLevel == 1) (*spectra)[pos].init(&(X[0]), &(Y[0]), m, !spi->CentroidMode, 0.6 / res, 2.0 / res, 1, 0.5 * (spi->StartRT + spi->EndRT), 0.0, 0.0);
									else if (spi->MSLevel == 2) (*spectra)[pos].init(&(X[0]), &(Y[0]), m, !spi->CentroidMode, 0.6 / res, 2.0 / res, 2, 0.5 * (spi->StartRT + spi->EndRT), wl, wh);
								} catch (System::Exception^ e) { std::cout << "ERROR: cannot read the .wiff file. Perhaps the respective .wiff.scan file is absent or corrupted?\n"; std::flush(std::cout);  provider->Close(); return; }
							} else {
								try {
									auto pl = exp->GetPeakArray(s);
									m = pl->Length;
									int i, k;
									auto &se = (*spectra)[pos];
									se.peaks.clear();
									se.MS_level = spi->MSLevel;
									se.RT = 0.5 * (spi->StartRT + spi->EndRT);
									if (spi->MSLevel == 1) se.window_low = se.window_high = 0.0;
									else se.window_low = wl, se.window_high = wh;
									for (i = k = 0; i < m; i++) if (pl[i]->area >= MinPeakHeight) k++;
									se.peaks.resize(k);
									for (i = k = 0; i < m; i++) if (pl[i]->area >= MinPeakHeight) {
										se.peaks[k].mz = pl[i]->xValue;
										se.peaks[k++].height = pl[i]->area;
									}
								} catch (System::Exception^ e) { std::cout << "ERROR: cannot read the .wiff file. Perhaps the respective .wiff.scan file is absent or corrupted?\n"; std::flush(std::cout);  provider->Close(); return; }
							}
						}
						pos++;
					}
				}
			} catch (System::Exception^ e) { std::cout << "ERROR: cannot read the .wiff file. Perhaps the respective .wiff.scan file is absent or corrupted?\n"; std::flush(std::cout); }
			provider->Close();
		} catch (System::Exception^ e) { std::cout << "ERROR: cannot read the .wiff file. Perhaps the respective .wiff.scan file is absent or corrupted?\n"; std::flush(std::cout);  return; }
	}
};

__declspec(dllexport) HANDLE diann_wiff_load(char * file, bool vendor, int Threads) {
	int i, e, ns = 0, verbose = false;
	long long mem, tot_peaks;
	if (verbose) {
		std::cout.setf(std::ios::unitbuf);
		std::cout << "Loading wiff file\n";
	}
	auto provider = gcnew AnalystWiffDataProvider();
	try {
		auto batch = AnalystDataProviderFactory::CreateBatch(gcnew System::String(file), provider);
		auto sample = batch->GetSample(0);
		auto mss = sample->MassSpectrometerSample;
		int en = mss->ExperimentCount;

		if (verbose) std::cout << en << " experiments\n";

		for (e = 0; e < en; e++) {
			auto exp = mss->GetMSExperiment(e);
			ns += exp->Details->NumberOfScans;
			if (verbose) std::cout << "Experiment " << e << ": " << exp->Details->NumberOfScans << " scans\n";
		}
		if (verbose) std::cout << ns << " total scans\n";
	} catch (System::Exception^ e) { std::cout << "ERROR: cannot read the .wiff file. Perhaps the respective .wiff.scan file is absent or corrupted?\n"; std::flush(std::cout); provider->Close(); return false; }

	provider->Close();
	std::vector<Spectrum> spectra(Max(1, ns));
	std::vector<Lock> locks(Max(1, ns));

	auto threads = gcnew System::Collections::Generic::List<Thread^>;
	for (i = 0; i < Threads; i++)
		threads->Add(gcnew Thread(gcnew ThreadStart(gcnew MTLoad(file, vendor, &spectra, &locks), &MTLoad::load_spectra)));
	for each(Thread ^th in threads) th->Start();
	for each(Thread ^th in threads) th->Join();

	if (verbose) std::cout << "File read successfully\n";

	for (i = mem = tot_peaks = 0; i < ns; i++) mem += spectra[i].mem_size(), tot_peaks += spectra[i].peaks.size();
	if (verbose) std::cout << "Total memory required: " << mem << " bytes\n";
	HANDLE data = GlobalAlloc(0, mem + 16);
	if (verbose) std::cout << "Memory allocated\n";
	int * ptr = (int*)data;
	*(ptr++) = ns;
	long long* lptr = (long long*)ptr; 
	*(lptr++) = tot_peaks; ptr = (int*)lptr;

	std::vector<int> index(ns);
	for (i = 0; i < ns; i++) index[i] = i;
	std::sort(index.begin(), index.end(), [&](const int &left, const int &right) { return spectra[left].RT < spectra[right].RT; });
	for (i = 0; i < ns; i++) ptr = spectra[index[i]].put(ptr);

	return data;
}