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
#ifndef _MZMLWRITER_H
#define _MZMLWRITER_H

#include <string>
#include <iostream>
#include "MSObject.h"
#include "MSReader.h"
#include "Spectrum.h"
#include "mzParser.h"

using namespace std;

namespace MSToolkit {

typedef struct sMzMLIndex{
  string id;
  f_off offset;
} sMzMLIndex;

class MzMLWriter {
public:

  MzMLWriter();
  ~MzMLWriter();

  bool  closeList(); //false is chromatogram list
  bool  closeMzML();
  bool  createList(bool specList=true); //false is chromatogram list
  bool  createMzML(char* fn);
  int   checkState();
  void  setNumpress(bool b);
  void  setTabs(bool b);
  void  setZlib(bool b);
  bool  writeRunInformation();
  bool  writeSpectra(MSObject& o);
  bool  writeSpectra(Spectrum& s);
  bool  writeChromatogram(BasicChromatogram& c);
  bool  writeIndex();

private:
  bool exportActivation(Spectrum& s, int tabs=0);
  bool exportAnalyzer();
  bool exportBinary(char* str, int len, int tabs=0);
  bool exportBinaryDataArray(BasicChromatogram& c, bool bRT, int tabs = 0);
  bool exportBinaryDataArray(Spectrum& s, bool bMZ, int tabs=0);
  bool exportBinaryDataArrayList(BasicChromatogram& c, int tabs = 0);
  bool exportBinaryDataArrayList(Spectrum& s, int tabs=0);
  bool exportChromatogram(BasicChromatogram& c, int tabs);
  bool exportChromatogramList();
  bool exportComponentList();
  bool exportContact();
  bool exportCv();
  bool exportCvList();
  bool exportCvParam(string ac, string ref, string name, string unitAc="", string unitRef="", string unitName="", string value="", int tabs=0);
  bool exportDataProcessing();
  bool exportDataProcessingList();
  bool exportFileContent();
  bool exportFileDescription();
  bool exportInstrumentConfiguration();
  bool exportInstrumentConfigurationList();
  bool exportIsolationWindow(BasicChromatogram& c, bool bPre, int tabs = 0);
  bool exportIsolationWindow(Spectrum& s, int tabs=0);
  bool exportMzML();
  bool exportOffset(string idRef, f_off offset, int tabs=0);
  bool exportPrecursor(BasicChromatogram& c, int tabs = 0);
  bool exportPrecursor(Spectrum& s, int tabs=0);
  bool exportPrecursorList(Spectrum& s, int tabs=0);
  bool exportProcessingMethod();
  bool exportProduct(BasicChromatogram&c, int tabs = 0);
  bool exportProductList();
  bool exportReferencableParamGroup();
  bool exportReferenceableParamGroupList();
  bool exportReferenceableParamGroupRef();
  bool exportRun();
  bool exportSample();
  bool exportSampleList();
  bool exportScan(Spectrum& s, int tabs=0);
  bool exportScanList(Spectrum& s, int tabs=0);
  bool exportScanSettings();
  bool exportScanSettingsList();
  bool exportScanWindow(Spectrum& s, int tabs = 0);
  bool exportScanWindowList(Spectrum& s, int tabs = 0);
  bool exportSelectedIon(BasicChromatogram& c, int tabs = 0);
  bool exportSelectedIon(Spectrum& s, int tabs=0);
  bool exportSelectedIonList(BasicChromatogram& c, int tabs = 0);
  bool exportSelectedIonList(Spectrum& s, int tabs=0);
  bool exportSoftware();
  bool exportSoftwareList();
  bool exportSoftwareRef();
  bool exportSource();
  bool exportSouceFile();
  bool exportSourceFileList();
  bool exportSourceFileRef();
  bool exportSourceFileRefList();
  bool exportSpectrum(Spectrum& s, int tabs=0);
  bool exportSpectrumList();
  void exportTabs(int tabs);
  bool exportTarget();
  bool exportTargetList();
  bool exportDetector();
  bool exportUserParam();

  int index;
  int chromIndex;
  FILE* fptr;
  f_off fSpecList;
  f_off fChromList;

  int iSpecList;
  int iChromList;
  bool bTabs;
  bool bFileOpen;
  bool bZlib;
  bool bNumpress;

  vector<sMzMLIndex> vIndex;
  vector<sMzMLIndex> vChromIndex;

};
}

#endif
