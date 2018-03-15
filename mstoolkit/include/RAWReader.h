/*
Copyright 2005-2016, Michael R. Hoopmann

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
#ifndef _RAWREADER_H
#define _RAWREADER_H

#ifdef _MSC_VER
#ifndef _NO_THERMORAW

#include "MSToolkitTypes.h"
#include "Spectrum.h"
#include <algorithm>
#include <iostream>
#include <vector>

#include <objbase.h>

#import "libid:F0C5F3E3-4F2A-443E-A74D-0AABE3237494" rename_namespace("XRawfile")
//#import "libid:5FE970A2-29C3-11D3-811D-00104B304896" rename_namespace("XRawfile")
using namespace XRawfile;
using namespace std;

namespace MSToolkit {

typedef struct rawPrecursorInfo{
  double  dIsoMZ;
  double  dMonoMZ;
  long    charge;
  long    parScanNum;
  rawPrecursorInfo(){
    dIsoMZ=0;
    dMonoMZ=0;
    charge=0;
    parScanNum=0;
  }
  void clear(){
    dIsoMZ = 0;
    dMonoMZ = 0;
    charge = 0;
    parScanNum = 0;
  }
} rawPrecursorInfo;

class RAWReader {
public:
	//Constructors & Destructors
  RAWReader();
  ~RAWReader();

	//Public Functions
  void getInstrument(char* str);
	long getLastScanNumber();
  void getManufacturer(char* str);
	long getScanCount();
	bool getStatus();

	bool lookupRT(char* c, int scanNum, float& rt);
	bool readRawFile(const char* c, Spectrum& s, int scNum=0);
  void setAverageRaw(bool b, int width=1, long cutoff=1000);
  void setLabel(bool b);  //label data contains all centroids (including noise and excluded peaks)
	void setMSLevelFilter(vector<MSSpectrumType>* v);
  void setRawFilter(char* c);
  void setRawFilterExact(bool b);
	

private:

	//Data Members
  bool bRaw;
	bool rawAvg;
	bool rawFileOpen;
	bool rawLabel;
	bool rawUserFilterExact;

  char rawCurrentFile[256];
  char rawInstrument[256];
  char rawManufacturer[256];
  char rawUserFilter[256];

	int rawAvgWidth;

	long rawAvgCutoff;
	long rawCurSpec;
	long rawTotSpec;
  
  IXRawfilePtr m_Raw;

	vector<MSSpectrumType>* msLevelFilter;

	//Private Functions
  int							calcChargeState(double precursormz, double highmass, VARIANT* varMassList, long nArraySize);
  double					calcPepMass(int chargestate, double precursormz);
  MSSpectrumType	evaluateFilter(long scan, char* chFilter, vector<double>& MZs, bool& bCentroid, double& cv, MSActivation& act);
  double          evaluateTrailerDouble(const char* id);
  int             evaluateTrailerInt(const char* id);
	bool						initRaw();
  

};

}
#endif
#endif
#endif

