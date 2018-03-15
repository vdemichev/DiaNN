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
#ifndef _MSTOOLKITTYPES_H
#define _MSTOOLKITTYPES_H

#include <stddef.h>

namespace MSToolkit {

enum MSSpectrumType {
  MS1,
  MS2,
  MS3,
  ZS,
  UZS,
  IonSpec,
  SRM,
  REFERENCE,
  Unspecified,
	MSX
};

enum MSFileFormat {
  bms1,
  bms2,
  cms1,
  cms2,
  mgf,
  ms1,
  ms2,
  msmat_ff,
  mzXML,
  mz5,
	mzML,
  raw,
  sqlite,
  psm,
  uzs,
  zs,
	mzXMLgz,
	mzMLgz,
  dunno
};

enum MSTag {
  no,
  D,
  H,
  I,
  S,
  Z
};

enum MSActivation {
  mstCID,
  mstECD,
  mstETD,
  mstETDSA,
  mstPQD,
  mstHCD,
	mstIRMPD,
  mstNA
};

struct MSHeader {
	char header[16][128];
};

struct MSScanInfo {
	int scanNumber[2];
	int numDataPoints;
  int numEZStates;
	int numZStates;
  float rTime;
  float IIT;
  float BPI;
  double* mz;
	int mzCount;
  double convA;
  double convB;
	double convC;
	double convD;
	double convE;
	double convI;
  double TIC;
  double BPM;
	MSScanInfo(){
    mz=NULL;
		scanNumber[0]=scanNumber[1]=0;
		numDataPoints=numEZStates=numZStates=0;
		rTime=IIT=BPI=0.0f;
		TIC=BPM=0.0;
		convA=convB=convC=convD=convE=convI=0.0;
    mzCount=0;
	}
	~MSScanInfo(){
		if(mz!=NULL) delete [] mz;
	}
};

struct DataPeak {
	double dMass;
	double dIntensity;
}; //For RAW files

struct Peak_T {
  double mz;
  float intensity;
};

struct ZState {
  int z;
  double mh;   //M+H, not mz
};

struct EZState {
  int z;
  double mh;      //M+H
  float pRTime;   //precursor area
  float pArea;    //precursor retention time
};

}

#endif


