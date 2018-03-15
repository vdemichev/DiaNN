/*  
  BasicSpectrum.cpp
  Copyright (C) 2010-2012, Mike Hoopmann
  Institute for Systems Biology
  Version 1.1, Mar. 28, 2012
*/

#include "mzParser.h"


//------------------------------------------
//  Constructors & Destructors
//------------------------------------------
BasicSpectrum::BasicSpectrum() {
  activation=none;
  basePeakIntensity=0.0;
  basePeakMZ=0.0;
  centroid=false;
  filterLine[0]='\0';
  highMZ=0.0;
  ionInjectionTime=0.0;
  lowMZ=0.0;
  msLevel=1;
  peaksCount=0;
  positiveScan=true;
  //precursorCharge[0]=0;
  //precursorCharge[1] = 0;
  //precursorCharge[2] = 0;
  //precursorCharge[3] = 0;
  //precursorCharge[4] = 0;
  //precursorChargeCount=0;
  //precursorIntensity=0.0;
  compensationVoltage=0.0;
  //precursorMonoMZ=0.0;
  //precursorMZ=0.0;
  precursorScanNum=-1;
  rTime=0.0f;
  scanIndex=0;
  scanNum=-1;
  totalIonCurrent=0.0;
  idString[0]='\0';
  vData=new vector<specDP>;
  vPrecursor=new vector<sPrecursorIon>;
}
BasicSpectrum::BasicSpectrum(const BasicSpectrum& s){
  vData=new vector<specDP>;
  vPrecursor=new vector<sPrecursorIon>;
  unsigned int i;
  for(i=0;i<s.vData->size();i++) vData->push_back(s.vData->at(i));
  for(i=0;i<s.vPrecursor->size();i++) vPrecursor->push_back(s.vPrecursor->at(i));
  //for(i=0;i<5;i++) precursorCharge[i] = s.precursorCharge[i];
  activation=s.activation;
  basePeakIntensity=s.basePeakIntensity;
  basePeakMZ=s.basePeakMZ;
  centroid=s.centroid;
  highMZ=s.highMZ;
  ionInjectionTime=s.ionInjectionTime;
  lowMZ=s.lowMZ;
  msLevel=s.msLevel;
  peaksCount=s.peaksCount;
  positiveScan=s.positiveScan;
  //precursorChargeCount=s.precursorChargeCount;
  //precursorIntensity=s.precursorIntensity;
  compensationVoltage=s.compensationVoltage;
  //precursorMonoMZ=s.precursorMonoMZ;
  //precursorMZ=s.precursorMZ;
  precursorScanNum=s.precursorScanNum;
  rTime=s.rTime;
  scanIndex=s.scanIndex;
  scanNum=s.scanNum;
  totalIonCurrent=s.totalIonCurrent;
  strcpy(idString,s.idString);
  strcpy(filterLine,s.filterLine);
}
BasicSpectrum::~BasicSpectrum() {
  delete vData;
  delete vPrecursor;
}

//------------------------------------------
//  Operator overloads
//------------------------------------------
BasicSpectrum& BasicSpectrum::operator=(const BasicSpectrum& s){
  if (this != &s) {
    delete vData;
    delete vPrecursor;
    vData=new vector<specDP>;
    vPrecursor=new vector<sPrecursorIon>;
    unsigned int i;
    for(i=0;i<s.vData->size();i++) vData->push_back(s.vData->at(i));
    for(i=0;i<s.vPrecursor->size();i++) vPrecursor->push_back(s.vPrecursor->at(i));
    //for(i=0;i<5;i++) precursorCharge[i] = s.precursorCharge[i];
    activation=s.activation;
    basePeakIntensity=s.basePeakIntensity;
    basePeakMZ=s.basePeakMZ;
    centroid=s.centroid;
    highMZ=s.highMZ;
    ionInjectionTime=s.ionInjectionTime;
    lowMZ=s.lowMZ;
    msLevel=s.msLevel;
    peaksCount=s.peaksCount;
    positiveScan=s.positiveScan;
    //precursorChargeCount=s.precursorChargeCount;
    //precursorIntensity=s.precursorIntensity;
    compensationVoltage=s.compensationVoltage;
    //precursorMonoMZ=s.precursorMonoMZ;
    //precursorMZ=s.precursorMZ;
    precursorScanNum=s.precursorScanNum;
    rTime=s.rTime;
    scanIndex=s.scanIndex;
    scanNum=s.scanNum;
    totalIonCurrent=s.totalIonCurrent;
    strcpy(filterLine,s.filterLine);
    strcpy(idString,s.idString);
  }
  return *this;
}
specDP& BasicSpectrum::operator[ ](const unsigned int index) {
  return vData->at(index);
}

//------------------------------------------
//  Modifiers
//------------------------------------------
void BasicSpectrum::addDP(specDP dp) { vData->push_back(dp);}
void BasicSpectrum::clear(){
  activation=none;
  basePeakIntensity=0.0;
  basePeakMZ=0.0;
  centroid=false;
  filterLine[0]='\0';
  highMZ=0.0;
  idString[0]='\0';
  ionInjectionTime=0.0;
  lowMZ=0.0;
  msLevel=1;
  peaksCount=0;
  positiveScan=true;
  //precursorCharge[0]=0;
  //precursorChargeCount=0;
  //precursorIntensity=0.0;
  compensationVoltage=0.0;
  //precursorMZ=0;
  precursorScanNum=-1;
  rTime=0.0f;
  scanIndex=0;
  scanNum=-1;
  totalIonCurrent=0.0;
  vData->clear();
  vPrecursor->clear();
}
void BasicSpectrum::setActivation(int a){ activation=a;}
void BasicSpectrum::setBasePeakIntensity(double d){ basePeakIntensity=d;}
void BasicSpectrum::setBasePeakMZ(double d){ basePeakMZ=d;}
void BasicSpectrum::setCentroid(bool b){ centroid=b;}
void BasicSpectrum::setCollisionEnergy(double d){ collisionEnergy=d;}
void BasicSpectrum::setCompensationVoltage(double d){ compensationVoltage=d; }
void BasicSpectrum::setFilterLine(char* str) { 
  strncpy(filterLine,str,127);
  filterLine[127]='\0';
}
void BasicSpectrum::setHighMZ(double d){ highMZ=d;}
void BasicSpectrum::setIDString(char* str) { 
  strncpy(idString,str,127); 
  idString[127]='\0';
}
void BasicSpectrum::setIonInjectionTime(double d){ ionInjectionTime=d;}
void BasicSpectrum::setLowMZ(double d){ lowMZ=d;}
void BasicSpectrum::setMSLevel(int level){ msLevel=level;}
void BasicSpectrum::setPeaksCount(int i){ peaksCount=i;}
void BasicSpectrum::setPositiveScan(bool b){ positiveScan=b;}
void BasicSpectrum::setPrecursorIon(sPrecursorIon& p){ vPrecursor->push_back(p);}
void BasicSpectrum::setPrecursorScanNum(int i){ precursorScanNum=i;}
void BasicSpectrum::setRTime(float f){ rTime=f;}
void BasicSpectrum::setScanIndex(int num) { scanIndex=num;}
void BasicSpectrum::setScanNum(int num){ scanNum=num;}
void BasicSpectrum::setTotalIonCurrent(double d){ totalIonCurrent=d;}

//------------------------------------------
//  Accessors
//------------------------------------------
int BasicSpectrum::getActivation(){ return activation;}
double BasicSpectrum::getBasePeakIntensity(){ return basePeakIntensity;}
double BasicSpectrum::getBasePeakMZ(){ return basePeakMZ;}
bool BasicSpectrum::getCentroid(){ return centroid;}
double BasicSpectrum::getCollisionEnergy(){ return collisionEnergy;}
double BasicSpectrum::getCompensationVoltage(){ return compensationVoltage;}
int BasicSpectrum::getFilterLine(char* str) {
  strcpy(str,filterLine);
  return (int)strlen(str);
}
double BasicSpectrum::getHighMZ(){ return highMZ;}
int BasicSpectrum::getIDString(char* str) { 
  strcpy(str,idString);
  return (int)strlen(str);
}
double BasicSpectrum::getIonInjectionTime(){return ionInjectionTime;}
double BasicSpectrum::getLowMZ(){ return lowMZ;}
int BasicSpectrum::getMSLevel(){ return msLevel;}
int BasicSpectrum::getPeaksCount(){ return peaksCount;}
bool BasicSpectrum::getPositiveScan(){ return positiveScan;}
int BasicSpectrum::getPrecursorCharge(int i){ //legacy function. Always returns first charge of requested precursor, or first charge of first precursor
  if (vPrecursor->size()==0) return 0;
  if (i>=vPrecursor->size()) return 0;
  return vPrecursor->at(i).charge;
}
//int BasicSpectrum::getPrecursorChargeCount(){ return precursorChargeCount;}
//double BasicSpectrum::getPrecursorIntensity(){ return precursorIntensity;}
//double BasicSpectrum::getPrecursorMonoMZ(){ return precursorMonoMZ;}
double BasicSpectrum::getPrecursorMZ(int i){ //legacy function. Always returns first precursor mz
  if (vPrecursor->size()==0) return 0;
  if (i >= vPrecursor->size()) return 0;
  return vPrecursor->at(i).mz;
}
sPrecursorIon BasicSpectrum::getPrecursorIon(int i){ return vPrecursor->at(i); }
int BasicSpectrum::getPrecursorIonCount() { return (int)vPrecursor->size(); }
int BasicSpectrum::getPrecursorScanNum(){ return precursorScanNum;}
float BasicSpectrum::getRTime(bool min){
  if(min) return rTime;
  else return rTime*60;
}
int BasicSpectrum::getScanIndex(){ return scanIndex;}
int BasicSpectrum::getScanNum(){ return scanNum;}
double BasicSpectrum::getTotalIonCurrent(){ return totalIonCurrent;}
size_t BasicSpectrum::size(){ return vData->size();}

