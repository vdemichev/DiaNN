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
#ifndef _NO_THERMORAW
#include "RAWReader.h"

using namespace MSToolkit;

// ==========================
// Constructors & Destructors
// ==========================
RAWReader::RAWReader(){

  CoInitialize( NULL );
	bRaw = initRaw();
  rawCurSpec=0;
  rawTotSpec=0;
  rawAvg=false;
  rawAvgWidth=1;
  rawAvgCutoff=1000;
	rawFileOpen=false;
  rawLabel=false;
  rawUserFilterExact=true;
  strcpy(rawCurrentFile,".");
  strcpy(rawInstrument,"unknown");
  strcpy(rawManufacturer,"Thermo Scientific");
  strcpy(rawUserFilter,"");
	msLevelFilter=NULL;

}

RAWReader::~RAWReader(){

  if(bRaw){
    if(rawFileOpen) m_Raw->Close();
    m_Raw.Release();
    m_Raw=NULL;
  }
	msLevelFilter=NULL;

}

int RAWReader::calcChargeState(double precursormz, double highmass, VARIANT* varMassList, long nArraySize) {
// Assumes spectrum is +1 or +2.  Figures out charge by
// seeing if signal is present above the parent mass
// indicating +2 (by taking ratio above/below precursor)

	bool bFound;
	long i, iStart;
	double dLeftSum,dRightSum;
	double FractionWindow;
	double CorrectionFactor;

	dLeftSum = 0.00001;
	dRightSum = 0.00001;

	DataPeak* pDataPeaks = NULL;
	SAFEARRAY FAR* psa = varMassList->parray;
	SafeArrayAccessData( psa, (void**)(&pDataPeaks) );

//-------------
// calc charge
//-------------
	bFound=false;
	i=0;
	while(i<nArraySize && !bFound){
    if(pDataPeaks[i].dMass < precursormz - 20){
			//do nothing
		} else {
			bFound = true;
      iStart = i;
    }
    i++;
	}
	if(!bFound) iStart = nArraySize;

	for(i=0;i<iStart;i++)	dLeftSum = dLeftSum + pDataPeaks[i].dIntensity;

	bFound=false;
	i=0;
	while(i<nArraySize && !bFound){
    if(pDataPeaks[i].dMass < precursormz + 20){
			//do nothing
		} else {
      bFound = true;
      iStart = i;
    }
    i++;
	}

	if(!bFound) {
		SafeArrayUnaccessData( psa );
		psa = NULL;
		pDataPeaks = NULL;
		return 1;
	}
	if(iStart = 0) iStart++;

	for(i=iStart;i<nArraySize;i++) dRightSum = dRightSum + pDataPeaks[i].dIntensity;

	if(precursormz * 2 < highmass){
    CorrectionFactor = 1;
	} else {
    FractionWindow = (precursormz * 2) - highmass;
    CorrectionFactor = (precursormz - FractionWindow) / precursormz;
	}

	if(dLeftSum > 0 && (dRightSum / dLeftSum) < (0.2 * CorrectionFactor)){
		SafeArrayUnaccessData( psa );
		psa=NULL;
		pDataPeaks=NULL;
		return 1;
	} else {
		SafeArrayUnaccessData( psa );
		psa=NULL;
		pDataPeaks=NULL;
    return 0;  //Set charge to 0 to indicate that both +2 and +3 spectra should be created
	}

  //When all else fails, return 0
  return 0;
}

MSSpectrumType RAWReader::evaluateFilter(long scan, char* chFilter, vector<double>& MZs, bool& bCentroid, double& cv, MSActivation& act) {

  BSTR Filter = NULL;
	char cStr[256];
	string tStr;
	string mzVal;
	int stop;
  bool bSA=false;

	//For non-ATL and non-MFC conversions
	int sl;

  //Initialize raw values to default
	MZs.clear();
  cv=0;

  m_Raw->GetFilterForScanNum(scan, &Filter);
	sl = SysStringLen(Filter)+1;
	WideCharToMultiByte(CP_ACP,0,Filter,-1,chFilter,sl,NULL,NULL);
	SysFreeString(Filter);

	strcpy(cStr,chFilter);
	MSSpectrumType mst=Unspecified;
	char* tok;
	tok=strtok(cStr," \n");
	while(tok!=NULL){

		if(strcmp(tok,"c")==0){
      bCentroid=true;
		} else if(strlen(tok)>2 && tok[0]=='c' && tok[1]=='v'){
      cv=atof(tok+3);
		} else if(strcmp(tok,"d")==0){
    } else if(strcmp(tok,"E")==0){ //enhanced 
		} else if(strcmp(tok,"ESI")==0){
		} else if(strcmp(tok,"FTMS")==0){
		} else if(strcmp(tok,"Full")==0){
		} else if(strcmp(tok,"ITMS")==0){
		} else if(strcmp(tok,"lock")==0){
		} else if(strcmp(tok,"ms")==0){
			mst=MS1;
		} else if(strcmp(tok,"msx")==0){
			mst=MSX;
		} else if(strcmp(tok,"ms2")==0){
			if(mst!=MSX) mst=MS2;
		} else if(strcmp(tok,"ms3")==0){
			if(mst!=MSX) mst=MS3;
		} else if(strcmp(tok,"NSI")==0){
		} else if(strcmp(tok,"p")==0){
      bCentroid=false;
    } else if(strcmp(tok,"r")==0){ //appears in fusion data, no documentation
    } else if(strcmp(tok,"sa")==0){
      bSA=true;
		} else if(strncmp(tok,"sid",3)==0){
		} else if(strcmp(tok,"SRM")==0){
			mst=SRM;
    } else if(strcmp(tok,"t")==0){ //turbo scan 
		} else if(strcmp(tok,"u")==0){
			mst=UZS;
		} else if(strcmp(tok,"w")==0){ //wideband activation?
		} else if(strcmp(tok,"Z")==0){
			if(mst!=UZS) mst=ZS;
		} else if(strcmp(tok,"+")==0){
		} else if(strcmp(tok,"-")==0){
		} else if(strchr(tok,'@')!=NULL){
			tStr=tok;
			stop=(int)tStr.find("@");
			mzVal=tStr.substr(0,stop);
			MZs.push_back(atof(&mzVal[0]));
      mzVal=tStr.substr(stop+1,3);
      if(mzVal.compare("cid")==0){
        act=mstCID;
      } else if(mzVal.compare("etd")==0){
        if(bSA) act=mstETDSA;
        else act=mstETD;
      } else if(mzVal.compare("hcd")==0){
        act=mstHCD;
      } else {
        cout << "Unknown activation method: " << &mzVal[0] << endl;
        act=mstNA;
      }

		} else if(strchr(tok,'[')!=NULL){
		} else {
			cout << "Unknown token: " << tok << endl;
		}

		tok=strtok(NULL," \n");
	}

	return mst;

}

double RAWReader::evaluateTrailerDouble(const char* str){
  VARIANT v;
  BSTR bs;
  int sl;
  double ret;

  VariantInit(&v);
  sl = lstrlenA(str);
  bs = SysAllocStringLen(NULL, sl);
  MultiByteToWideChar(CP_ACP, 0, str, sl, bs, sl);
  m_Raw->GetTrailerExtraValueForScanNum(rawCurSpec, bs, &v);
  SysFreeString(bs);
  if (v.vt == VT_R4) ret = (double)v.fltVal;
  else if (v.vt == VT_R8) ret = (double)v.dblVal;
  else ret=0;

  VariantClear(&v);
  return ret;
}

int RAWReader::evaluateTrailerInt(const char* str){
  VARIANT v;
  BSTR bs;
  int sl;
  int ret;

  VariantInit(&v);
  sl = lstrlenA(str);
  bs = SysAllocStringLen(NULL, sl);
  MultiByteToWideChar(CP_ACP, 0, str, sl, bs, sl);
  m_Raw->GetTrailerExtraValueForScanNum(rawCurSpec, bs, &v);
  SysFreeString(bs);
  if (v.vt == VT_I2) ret = (int)v.iVal;
  else if (v.vt == VT_I4) ret = (int)v.lVal;
  else if (v.vt == VT_I1) ret = (int)v.cVal;
  else if (v.vt == VT_UI4) ret = (int)v.ulVal;
  else if (v.vt == VT_UI2) ret = (int)v.uiVal;
  else if (v.vt == VT_UI1) ret = (int)v.bVal;
  else if (v.vt == VT_INT) ret = (int)v.intVal;
  else if (v.vt == VT_UINT) ret = (int)v.uintVal;
  else ret = 0;

  VariantClear(&v);
  return ret;
}

void RAWReader::getInstrument(char* str){
  strcpy(str,rawInstrument);
}

long RAWReader::getLastScanNumber(){
	return rawCurSpec;
}

void RAWReader::getManufacturer(char* str){
  strcpy(str,rawManufacturer);
}

long RAWReader::getScanCount(){
	return rawTotSpec;
}

bool RAWReader::getStatus(){
	return bRaw;
}

bool RAWReader::initRaw(){

	int raw=0;

	IXRawfile2Ptr m_Raw2;
	IXRawfile3Ptr m_Raw3;
	IXRawfile4Ptr m_Raw4;
	IXRawfile5Ptr m_Raw5;

	//Example of Xcalibur/Foundation first
	//if(FAILED(m_Raw5.CreateInstance("XRawfile.XRawfile.1"))){'

	//Try MSFileReader - using ProteoWizard strategy
  if(FAILED(m_Raw5.CreateInstance("MSFileReader.XRawfile.1"))){
		if(FAILED(m_Raw4.CreateInstance("MSFileReader.XRawfile.1"))){
			if(FAILED(m_Raw3.CreateInstance("MSFileReader.XRawfile.1"))){
				if(FAILED(m_Raw2.CreateInstance("MSFileReader.XRawfile.1"))){
					if(FAILED(m_Raw.CreateInstance("MSFileReader.XRawfile.1"))){
            raw=0;
						//cout << "Cannot load Thermo MSFileReader. Cannot read .RAW files." << endl;
					} else {
						raw=1;
					}
				} else {
					m_Raw=m_Raw2;
					raw=2;
				}
			} else {
				m_Raw=m_Raw3;
				raw=3;
			}
		} else {
			m_Raw=m_Raw4;
			raw=4;
		}
	} else {
		m_Raw=m_Raw5;
		raw=5;
	}
	
	if(raw>0) return true;
	return false;
}

bool RAWReader::readRawFile(const char *c, Spectrum &s, int scNum){

	//General purpose function members
	bool bCheckNext;
  bool bNewFile;

	char chFilter[256];
  char curFilter[256];
	
	double dRTime;
	double highmass=0.0;
	double pm1;
	double pw;

	long i;
  long j;
	long lArraySize=0;
  long ret;

	vector<double> MZs;

	DataPeak* pDataPeaks = NULL;
  HRESULT lRet;
	MSSpectrumType MSn;
	SAFEARRAY FAR* psa;
  TCHAR pth[MAX_PATH];
  VARIANT varMassList;
	VARIANT varPeakFlags;

  //Members for gathering averaged scans
	int charge;
  int sl;
	int widthCount;
	
	long FirstBkg1=0;
	long FirstBkg2=0;
  long LastBkg1=0;
  long LastBkg2=0;
	long lowerBound;
  long upperBound;
    
  BSTR rawFilter=NULL;
  BSTR testStr;

	//Additional members for Scan Information
  bool bCentroid;

  double cv;    //Compensation Voltage
  double BPI;   //Base peak intensity
	double BPM;   //Base peak mass
	double td;    //temp double value
	double TIC;
  float IIT;    //ion injection time

	long tl;      //temp long value
  MSActivation act;


  if(!bRaw) return false;

	//Clear spectrum object
  s.clear();

  if(c==NULL){
		//continue reading current file
    if(!rawFileOpen) return false;
    if(scNum>0) rawCurSpec=scNum;
    else rawCurSpec++;
    if(rawCurSpec>rawTotSpec) return false;
    bNewFile=false;
  } else {
	
    //check if requested file is already open
    if(rawFileOpen) {
      if(strcmp(c,rawCurrentFile)==0){
        if(scNum>0) rawCurSpec=scNum;
        else rawCurSpec++;
        if(rawCurSpec>rawTotSpec) return false;
        bNewFile=false;
      } else {
        //new file requested, so close the existing one
        lRet = m_Raw->Close();
        rawFileOpen=false;
        bNewFile=true;
      }
    } else {
      bNewFile=true;
    }

    if(bNewFile){
      MultiByteToWideChar(CP_ACP,0,c,-1,(LPWSTR)pth,MAX_PATH);
      lRet = m_Raw->Open((LPWSTR)pth);
		  if(lRet != ERROR_SUCCESS) {
			  cerr << "Cannot open " << c << endl;
			  return false;
		  }
	    else lRet = m_Raw->SetCurrentController(0,1);
      rawFileOpen=true;
      m_Raw->GetNumSpectra(&rawTotSpec);
      testStr=NULL;
      m_Raw->GetInstModel(&testStr);
      sl = SysStringLen(testStr)+1;
	    WideCharToMultiByte(CP_ACP,0,testStr,-1,rawInstrument,sl,NULL,NULL);
      SysFreeString(testStr);
      strcpy(rawCurrentFile,c);

		  //if scan number is requested, grab it
      if(scNum>0) rawCurSpec=scNum;
      else rawCurSpec=1;
      if(rawCurSpec>rawTotSpec) return false;
    }
  }

	//Initialize members
	strcpy(chFilter,"");
  strcpy(curFilter,"");
	VariantInit(&varMassList);
	VariantInit(&varPeakFlags);

  rawPrecursorInfo preInfo;

	//Rather than grab the next scan number, get the next scan based on a user-filter (if supplied).
  //if the filter was set, make sure we pass the filter
  while(true){

	  MSn = evaluateFilter(rawCurSpec, curFilter, MZs, bCentroid,cv,act);

    //check for spectrum filter (string)
    if(strlen(rawUserFilter)>0){
      bCheckNext=false;
      if(rawUserFilterExact) {
        if(strcmp(curFilter,rawUserFilter)!=0) bCheckNext=true;
      } else {
        if(strstr(curFilter,rawUserFilter)==NULL) bCheckNext=true;
      }

      //if string doesn't match, get next scan until it does match or EOF
      if(bCheckNext){
        if(scNum>0) return false;
        rawCurSpec++;
        if(rawCurSpec>rawTotSpec) return false;
        continue;
      }
    }

    //check for msLevel filter
    if(msLevelFilter->size()>0 && find(msLevelFilter->begin(), msLevelFilter->end(), MSn) == msLevelFilter->end()) {
      if(scNum>0) return false;
      rawCurSpec++;
      if(rawCurSpec>rawTotSpec) return false;
    } else {
      break;
    }
  }

  //Get basic spectrum metadata. It will be replaced/supplemented later, if available
  preInfo.clear();
  m_Raw->GetScanHeaderInfoForScanNum(rawCurSpec, &tl, &td, &td, &td, &TIC, &BPM, &BPI, &tl, &tl, &td);
  m_Raw->RTFromScanNum(rawCurSpec, &dRTime);
  preInfo.charge = evaluateTrailerInt("Charge State:");
  preInfo.dMonoMZ = evaluateTrailerDouble("Monoisotopic M/Z:");
  IIT = (float)evaluateTrailerDouble("Ion Injection Time (ms):");

  //Get more sig digits for isolation mass
  if (raw > 3 && (MSn==MS2 || MSn==MS3)){
    IXRawfile4Ptr raw4 = (IXRawfile4Ptr)m_Raw;
    if (MSn==MS2) tl=2;
    else tl=3;
    ret = raw4->GetPrecursorMassForScanNum(rawCurSpec,tl,&td);
    if (ret == 0) preInfo.dIsoMZ = td;
    raw4=NULL;

    //Correct precursor mono mass if it is more than 5 13C atoms away from isolation mass
    if (preInfo.dMonoMZ > 0 && preInfo.charge>0){
      td = preInfo.dIsoMZ-preInfo.dMonoMZ;
      if (td>5.01675/preInfo.charge) preInfo.dMonoMZ=preInfo.dIsoMZ;
    }

  }

  //Get the peaks
	//Average raw files if requested by user
  if(rawAvg){
    widthCount=0;
    lowerBound=0;
    upperBound=0;
    for(i=rawCurSpec-1;i>0;i--){
      evaluateFilter(i, chFilter, MZs, bCentroid,cv,act);
      if(strcmp(curFilter,chFilter)==0){
        widthCount++;
        if(widthCount==rawAvgWidth) {
          lowerBound=i;
          break;
        }
      }
    }
    if(lowerBound==0) lowerBound=rawCurSpec; //this will have "edge" effects

    widthCount=0;
    for(i=rawCurSpec+1;i<rawTotSpec;i++){
      evaluateFilter(i, chFilter, MZs, bCentroid,cv,act);
      if(strcmp(curFilter,chFilter)==0){
        widthCount++;
        if(widthCount==rawAvgWidth) {
          upperBound=i;
          break;
        }
      }
    }
    if(upperBound==0) upperBound=rawCurSpec; //this will have "edge" effects

    m_Raw->GetFilterForScanNum(i, &rawFilter);
    j=m_Raw->GetAverageMassList(&lowerBound, &upperBound, &FirstBkg1, &LastBkg1, &FirstBkg2, &LastBkg2,
      rawFilter, 1, rawAvgCutoff, 0, FALSE, &pw, &varMassList, &varPeakFlags, &lArraySize );
    SysFreeString(rawFilter);
    rawFilter=NULL;

  } else {

		//Get regular spectrum data
		sl=lstrlenA("");
		testStr = SysAllocStringLen(NULL,sl);
		MultiByteToWideChar(CP_ACP,0,"",sl,testStr,sl);
		j=m_Raw->GetMassListFromScanNum(&rawCurSpec,testStr,0,0,0,FALSE,&pw,&varMassList,&varPeakFlags,&lArraySize);
		SysFreeString(testStr);
  }

	//Handle MS2 and MS3 files differently to create Z-lines
	if(MSn==MS2 || MSn==MS3){

		//if charge state is assigned to spectrum, add Z-lines.
    if (preInfo.dIsoMZ>0.01) MZs[0]=preInfo.dIsoMZ;  //overwrite isolation mass if we have more sig figs.
		if(preInfo.charge>0){ //if(Charge.iVal>0){
			if(preInfo.dMonoMZ>0.01) { //if(MonoMZ.dblVal>0.01) {
				//pm1 = MonoMZ.dblVal * Charge.iVal - ((Charge.iVal-1)*1.007276466);
        pm1 = preInfo.dMonoMZ * preInfo.charge - ((preInfo.charge - 1)*1.007276466);
        s.setMZ(MZs[0], preInfo.dMonoMZ);
			}	else {
        pm1 = MZs[0] * preInfo.charge - ((preInfo.charge - 1)*1.007276466);
        s.setMZ(MZs[0]);
			}
      s.addZState(preInfo.charge, pm1);
      s.setCharge(preInfo.charge);
    } else {
      s.setMZ(preInfo.dIsoMZ); s.setMZ(MZs[0]);
      charge = calcChargeState(MZs[0], highmass, &varMassList, lArraySize);

      //Charge greater than 0 means the charge state is known
      if(charge>0){
        pm1 = MZs[0]*charge - ((charge-1)*1.007276466);
  	    s.addZState(charge,pm1);

      //Charge of 0 means unknown charge state, therefore, compute +2 and +3 states.
      } else {
        pm1 = MZs[0]*2 - 1.007276466;
        s.addZState(2,pm1);
        pm1 = MZs[0]*3 - 2*1.007276466;
        s.addZState(3,pm1);
      }

    }

  } //endif MS2 and MS3

	if(MSn==MSX){
		for(i=0;i<(int)MZs.size();i++){
			if(i==0) s.setMZ(MZs[i],0);
			else s.addMZ(MZs[i],0);
		}
		s.setCharge(0);
	}

	//Set basic scan info
  if(bCentroid) s.setCentroidStatus(1);
  else s.setCentroidStatus(0);
  s.setRawFilter(curFilter);
  s.setActivationMethod(act);
	s.setScanNumber((int)rawCurSpec);
  s.setScanNumber((int)rawCurSpec,true);
	s.setRTime((float)dRTime);
	s.setFileType(MSn);
  s.setBPI((float)BPI);
  s.setBPM(BPM);
  s.setCompensationVoltage(cv);
//  s.setConversionA(ConversionA.dblVal);  //Deprecating Conversion values
//  s.setConversionB(ConversionB.dblVal);
//	s.setConversionC(ConversionC.dblVal);
//  s.setConversionD(ConversionD.dblVal);
//	s.setConversionE(ConversionE.dblVal);
//  s.setConversionI(ConversionI.dblVal);
  s.setTIC(TIC);
  s.setIonInjectionTime(IIT);
	if(MSn==SRM) s.setMZ(MZs[0]);
  switch(MSn){
    case MS1: s.setMsLevel(1); break;
    case MS2: s.setMsLevel(2); break;
    case MS3: s.setMsLevel(3); break;
		case MSX: s.setMsLevel(2); break;
    default: s.setMsLevel(0); break;
  }
	psa = varMassList.parray;
  SafeArrayAccessData( psa, (void**)(&pDataPeaks) );
  for(j=0;j<lArraySize;j++)	s.add(pDataPeaks[j].dMass,(float)pDataPeaks[j].dIntensity);
  SafeArrayUnaccessData( psa );

  //Clean up memory
	VariantClear(&varMassList);
	VariantClear(&varPeakFlags);

  return true;

}

void RAWReader::setMSLevelFilter(vector<MSSpectrumType>* v){
	msLevelFilter=v;
}

void RAWReader::setRawFilter(char *c){
  strcpy(rawUserFilter,c);
}

#endif
