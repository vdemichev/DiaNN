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
#include "Spectrum.h"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace MSToolkit;

Spectrum::Spectrum(){
  //cout<<"in Spectrum constructor!"<<endl;

  rTime=0;
  charge=0;
  scanNumber=0;
  scanNumber2=0;
  msLevel = 2;
  monoMZ=new vector<double>;
  mz=new vector<double>;
  TIC=0;
  IIT=0;
  compensationVoltage=0;
  convA=0;
  convB=0;
	convC=0;
	convD=0;
	convE=0;
	convI=0;
  BPI=0;
  BPM=0;
  selectionWinLower = 0;
  selectionWinUpper = 0;
  centroidStatus=2;

  fileType=Unspecified;
  vPeaks = new vector<Peak_T>;
  vEZ = new vector<EZState>;
  vZ = new vector<ZState>;
  actMethod=mstNA;

  strcpy(rawFilter,"");
  strcpy(nativeID,"");
}



Spectrum::~Spectrum(){
  if(vPeaks) delete vPeaks;
  if(vEZ) delete vEZ;
  if(vZ) delete vZ;
	if(mz) delete mz;
  if(monoMZ) delete monoMZ;
}

Spectrum::Spectrum(const Spectrum& s){

  unsigned int i;

  rTime = s.rTime;
  charge = s.charge;
  scanNumber = s.scanNumber;
  scanNumber2 = s.scanNumber2;
  msLevel = s.msLevel;
  monoMZ = new vector<double>;
  for(i=0;i<s.monoMZ->size();i++){
		monoMZ->push_back(s.monoMZ->at(i));
	}
  mz = new vector<double>;
	for(i=0;i<s.mz->size();i++){
		mz->push_back(s.mz->at(i));
	}
  fileType = s.fileType;
  IIT = s.IIT;
  TIC = s.TIC;
  compensationVoltage = s.compensationVoltage;
  convA = s.convA;
  convB = s.convB;
  convC = s.convC;
  convD = s.convD;
  convE = s.convE;
  convI = s.convI;
  BPI = s.BPI;
  BPM = s.BPM;
  selectionWinLower=s.selectionWinLower;
  selectionWinUpper=s.selectionWinUpper;
  centroidStatus = s.centroidStatus;
  vPeaks = new vector<Peak_T>;
  for(i=0;i<s.vPeaks->size();i++){
    vPeaks->push_back(s.vPeaks->at(i));
  }
  vEZ = new vector<EZState>;
  for(i=0;i<s.vEZ->size();i++){
    vEZ->push_back(s.vEZ->at(i));
  }
  vZ = new vector<ZState>;
  for(i=0;i<s.vZ->size();i++){
    vZ->push_back(s.vZ->at(i));
  }
  strcpy(rawFilter,s.rawFilter);
  strcpy(nativeID,s.nativeID);
}

Spectrum& Spectrum::operator=(const Spectrum& s){
	//cout<<"in Spectrum ="<<endl;
  unsigned int i;
  if (this != &s) {
    delete vPeaks;
    delete vEZ;
    delete vZ;
    delete monoMZ;
		delete mz;
    monoMZ = new vector<double>;
    for(i=0;i<s.monoMZ->size();i++){
		  monoMZ->push_back(s.monoMZ->at(i));
	  }
		mz = new vector<double>;
		for(i=0;i<s.mz->size();i++){
			mz->push_back(s.mz->at(i));
		}
    vPeaks = new vector<Peak_T>;
    for(i=0;i<s.vPeaks->size();i++){
      vPeaks->push_back(s.vPeaks->at(i));
    }
    vEZ = new vector<EZState>;
    for(i=0;i<s.vEZ->size();i++){
      vEZ->push_back(s.vEZ->at(i));
    }
    vZ = new vector<ZState>;
    for(i=0;i<s.vZ->size();i++){
      vZ->push_back(s.vZ->at(i));
    }
    rTime = s.rTime;
    charge = s.charge;
    scanNumber = s.scanNumber;
    scanNumber2 = s.scanNumber2;
    msLevel = s.msLevel;
    BPI = s.BPI;
    BPM = s.BPM;
    compensationVoltage = s.compensationVoltage;
    convA = s.convA;
    convB = s.convB;
    convC = s.convC;
    convD = s.convD;
    convE = s.convE;
    convI = s.convI;
    TIC = s.TIC;
    IIT = s.IIT;
    fileType = s.fileType;
    centroidStatus = s.centroidStatus;
    selectionWinLower = s.selectionWinLower;
    selectionWinUpper = s.selectionWinUpper;
    strcpy(rawFilter,s.rawFilter);
    strcpy(nativeID,s.nativeID);
  }
  return *this;
}

Peak_T& Spectrum::operator[](const int& i) {
	return vPeaks->operator[](i);
}


/* ----- Begin Functions ----- */


/* Adds Result struct to end of spectrum. */
void Spectrum::add(Peak_T& p){
  vPeaks->push_back(p);
}

void Spectrum::add(double d1, float d2){
  Peak_T p;
  p.mz=d1;
  p.intensity=d2;
  vPeaks->push_back(p);
}

void Spectrum::addEZState(EZState& z){
	vEZ->push_back(z);
}

void Spectrum::addEZState(int i, double d, float f1, float f2){
	EZState z;
	z.z=i;
	z.mh=d;
  z.pRTime=f1;
  z.pArea=f2;
	vEZ->push_back(z);
}

void Spectrum::addMZ(double d, double mono){
	mz->push_back(d);
  monoMZ->push_back(mono);
}

void Spectrum::addZState(ZState& z){
	vZ->push_back(z);
}

void Spectrum::addZState(int i, double d){
	ZState z;
	z.z=i;
	z.mh=d;
	vZ->push_back(z);
}

/* Returns Result struct of single element in the spectrum. */
Peak_T& Spectrum::at(const int& i){
  return vPeaks->operator [](i);
}

Peak_T& Spectrum::at(const unsigned int& i){
  return vPeaks->operator [](i);
}

EZState& Spectrum::atEZ(const int& i){
	return vEZ->operator [](i);
}

EZState& Spectrum::atEZ(const unsigned int& i){
	return vEZ->operator [](i);
}

ZState& Spectrum::atZ(const int& i){
	return vZ->operator [](i);
}

ZState& Spectrum::atZ(const unsigned int& i){
	return vZ->operator [](i);
}

/* Clears the spectrum */
void Spectrum::clear(){
	delete vPeaks;
	vPeaks = new vector<Peak_T>;
  delete vEZ;
  vEZ = new vector<EZState>;
	delete vZ;
	vZ = new vector<ZState>;
	delete mz;
	mz = new vector<double>;
  delete monoMZ;
  monoMZ = new vector<double>;
	scanNumber = 0;
  scanNumber2 = 0;
	rTime = 0;
	charge = 0;
	msLevel = 0;
  convA = 0;
  convB = 0;
  TIC = 0;
  IIT = 0;
  BPI = 0;
  BPM = 0;
  selectionWinLower=0;
  selectionWinUpper=0;
	fileType = Unspecified;
  actMethod=mstNA;
}

void Spectrum::clearMZ(){
	delete mz;
	mz = new vector<double>;
  delete monoMZ;
	monoMZ = new vector<double>;
}

void Spectrum::clearPeaks(){
	delete vPeaks;
	vPeaks = new vector<Peak_T>;
}

/* Erases element i in the spectrum. */
void Spectrum::erase(unsigned int i){
  vector<Peak_T>::iterator vi;
  vi=vPeaks->begin()+i;
  vPeaks->erase(vi);
}

/* Erases element i to element j, inclusive, in the spectrum. */
void Spectrum::erase(unsigned int i, unsigned int j){
  vector<Peak_T>::iterator vi1;
  vector<Peak_T>::iterator vi2;
  vi1=vPeaks->begin()+i;
  vi2=vPeaks->begin()+j+1;
  vPeaks->erase(vi1,vi2);
}

void Spectrum::eraseEZ(unsigned int i){
  vector<EZState>::iterator vi;
  vi=vEZ->begin()+i;
  vEZ->erase(vi);
}

/* Erases element i to element j, inclusive, in the spectrum. */
void Spectrum::eraseEZ(unsigned int i, unsigned int j){
  vector<EZState>::iterator vi1;
  vector<EZState>::iterator vi2;
  vi1=vEZ->begin()+i;
  vi2=vEZ->begin()+j+1;
  vEZ->erase(vi1,vi2);
}

void Spectrum::eraseZ(unsigned int i){
  vector<ZState>::iterator vi;
  vi=vZ->begin()+i;
  vZ->erase(vi);
}

/* Erases element i to element j, inclusive, in the spectrum. */
void Spectrum::eraseZ(unsigned int i, unsigned int j){
  vector<ZState>::iterator vi1;
  vector<ZState>::iterator vi2;
  vi1=vZ->begin()+i;
  vi2=vZ->begin()+j+1;
  vZ->erase(vi1,vi2);
}

MSActivation Spectrum::getActivationMethod(){
  return actMethod;
}

float Spectrum::getBPI(){
  return BPI;
}

double Spectrum::getBPM(){
  return BPM;
}

int Spectrum::getCentroidStatus(){
  return centroidStatus;
}

int Spectrum::getCharge(){
  return charge;
}

double Spectrum::getCompensationVoltage(){
  return compensationVoltage;
}

double Spectrum::getConversionA(){
  return convA;
}

double Spectrum::getConversionB(){
  return convB;
}

double Spectrum::getConversionC(){
  return convC;
}

double Spectrum::getConversionD(){
  return convD;
}

double Spectrum::getConversionE(){
  return convE;
}

double Spectrum::getConversionI(){
  return convI;
}

MSSpectrumType Spectrum::getFileType(){
	return fileType;
}

float Spectrum::getIonInjectionTime(){
  return IIT;
}

double Spectrum::getMonoMZ(int index){
	if(index>=(int)monoMZ->size()) return -1.0;
  return monoMZ->at(index);
}

double Spectrum::getMZ(int index){
	if(index>=(int)mz->size()) return -1.0;
  return mz->at(index);
}

bool Spectrum::getNativeID(char* c, int sz){
  if(sz<(int)strlen(nativeID)) {
    cout << "Buffer too small to retrieve spectrumNativeID. " << sizeof(c) << " " << strlen(nativeID) << endl;
    return false;
  } else {
    strcpy(c,nativeID);
    return true;
  }
}

bool Spectrum::getRawFilter(char* c, int sz, bool bLock){
  if(sz<(int)strlen(rawFilter)) {
    cout << "Buffer too small to retrieve RAW filter. " << sizeof(c) << " " << strlen(rawFilter) << endl;
    return false;
  } else {
		strcpy(c,rawFilter);
		char* chp=strstr(c,"lock");
		if(!bLock && chp!=NULL) strcpy(chp,chp+5);
    return true;
  }
}

float Spectrum::getRTime(){
  return rTime;
}

int Spectrum::getScanNumber(bool second){
  if(second) return scanNumber2;
  else return scanNumber;
}

double Spectrum::getScanWindowLower(){
  return scanWinLower;
}

double Spectrum::getScanWindowUpper(){
  return scanWinUpper;
}

double Spectrum::getSelWindowLower(){
  return selectionWinLower;
}

double Spectrum::getSelWindowUpper(){
  return selectionWinUpper;
}

double Spectrum::getTIC(){
  return TIC;
}

void Spectrum::setBPI(float f){
  BPI=f;
}

void Spectrum::setBPM(double d){
  BPM=d;
}

void Spectrum::setCentroidStatus(int i){
  if(i>2) centroidStatus=2;
  else centroidStatus=i;
}

void Spectrum::setCharge(int i){
  charge=i;
}

void Spectrum::setCompensationVoltage(double d){
  compensationVoltage=d;
}

void Spectrum::setConversionA(double d){
  convA=d;
}

void Spectrum::setConversionB(double d){
  convB=d;
}

void Spectrum::setConversionC(double d){
  convC=d;
}

void Spectrum::setConversionD(double d){
  convD=d;
}

void Spectrum::setConversionE(double d){
  convE=d;
}

void Spectrum::setConversionI(double d){
  convI=d;
}

void Spectrum::setFileType(MSSpectrumType f){
	fileType=f;
}

void Spectrum::setIonInjectionTime(float f){
  IIT=f;
}

void Spectrum::setMZ(double d, double mono){
	clearMZ();
	mz->push_back(d);
  monoMZ->push_back(mono);
}

void Spectrum::setNativeID(char* c){
  if(strlen(c)>256) cout << "Error - spectrumNativeID filter larger than 256 characters." << endl;
  else strcpy(nativeID,c);
}

void Spectrum::setRawFilter(char* c){
  if(strlen(c)>256) cout << "Error - RAW filter larger than 256 characters." << endl;
  else strcpy(rawFilter,c);
}

void Spectrum::setRTime(float d){
  rTime=d;
}

void Spectrum::setScanNumber(int i, bool second){
  if(second)scanNumber2=i;
  else scanNumber=i;
}

void Spectrum::setScanWindow(double lower, double upper){
  scanWinLower = lower;
  scanWinUpper = upper;
}

void Spectrum::setSelWindow(double lower, double upper){
  selectionWinLower=lower;
  selectionWinUpper=upper;
}

void Spectrum::setTIC(double d){
  TIC=d;
}

void Spectrum::setMsLevel(int level) {
  msLevel = level;
}

int Spectrum::getMsLevel() {
  return msLevel;
}

int Spectrum::getScanID(){
  return scanID;
}

void Spectrum::setScanID(int scanid){
  scanID = scanid;
}

/* Returns the number of elements in the spectrum. */
int Spectrum::size(){
  return (int)vPeaks->size();
}

int Spectrum::sizeEZ(){
	return (int)vEZ->size();
}

int Spectrum::sizeMZ(){
	return (int)mz->size();
}

int Spectrum::sizeZ(){
	return (int)vZ->size();
}

float Spectrum::getTotalIntensity(){
  float totalIntensity = 0;
  for(unsigned int i=0; i<vPeaks->size(); i++) totalIntensity += (vPeaks->at(i)).intensity;
  return totalIntensity;
}

/* Sorts the spectrum by Data. */
void Spectrum::sortIntensity(){
  qsort(&vPeaks->at(0),vPeaks->size(),sizeof(Peak_T),compareIntensity);
}

/* Sorts the spectrum by Mass. */
void Spectrum::sortMZ(){
  qsort(&vPeaks->at(0),vPeaks->size(),sizeof(Peak_T),compareMZ);
}

/* Sorts the spectrum by Data. */
void Spectrum::sortIntensityRev(){
  qsort(&vPeaks->at(0),vPeaks->size(),sizeof(Peak_T),compareIntensityRev);
}

/* Sorts the spectrum by Mass. */
void Spectrum::sortMZRev(){
  qsort(&vPeaks->at(0),vPeaks->size(),sizeof(Peak_T),compareMZRev);
}

//const vector<Peak_T>* Spectrum::getPeaks(){
//	return vPeaks;
//};

vector<Peak_T>* Spectrum::getPeaks(){
	return vPeaks;
}

void Spectrum::setPeaks(vector<Peak_T> peaks) {
  if(!vPeaks->empty()) vPeaks->clear();
  for(unsigned int i=0; i<peaks.size(); i++){
    vPeaks->push_back(peaks.at(i));
  }
}

void Spectrum::setActivationMethod(MSActivation m){
  actMethod=m;
}

void Spectrum::printMe() {
  cout << "Scan Number: " << getScanNumber() << endl
       << "Mass to charge: " << getMZ() << endl
       << "S Charge: " << getCharge() << endl 
       << "RT: " << getRTime() << endl;
    
  cout << fixed;
 
  for(unsigned int i=0; i<vPeaks->size(); i++) {
    cout << setprecision(10) << vPeaks->at(i).mz << "  " << vPeaks->at(i).intensity << endl;
  }
}

//Private Functions

/* For the qsort */
int Spectrum::compareIntensity(const void *p1, const void *p2){
  const Peak_T d1 = *(Peak_T *)p1;
  const Peak_T d2 = *(Peak_T *)p2;
  if(d1.intensity<d2.intensity) return -1;
  else if(d1.intensity>d2.intensity) return 1;
  else return 0;
}

/* For the qsort */
int Spectrum::compareMZ(const void *p1, const void *p2){
  const Peak_T d1 = *(Peak_T *)p1;
  const Peak_T d2 = *(Peak_T *)p2;
  if(d1.mz<d2.mz) return -1;
  else if(d1.mz>d2.mz) return 1;
  else return 0;
}

/* For the qsort */
int Spectrum::compareIntensityRev(const void *p1, const void *p2){
  const Peak_T d1 = *(Peak_T *)p1;
  const Peak_T d2 = *(Peak_T *)p2;
  if(d1.intensity>d2.intensity) return -1;
  else if(d1.intensity<d2.intensity) return 1;
  else return 0;
}

/* For the qsort */
int Spectrum::compareMZRev(const void *p1, const void *p2){
  const Peak_T d1 = *(Peak_T *)p1;
  const Peak_T d2 = *(Peak_T *)p2;
  if(d1.mz>d2.mz) return -1;
  else if(d1.mz<d2.mz) return 1;
  else return 0;
}

