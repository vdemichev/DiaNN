/*
BasicChromatogram - The code is
open source under the FreeBSD License, please see LICENSE file
for detailed information.

Copyright (C) 2011, Mike Hoopmann, Institute for Systems Biology
Version 1.0, January 4, 2011.
Version 1.1, March 14, 2012.
*/

#include "mzParser.h"

BasicChromatogram::BasicChromatogram(){}

BasicChromatogram::BasicChromatogram(const BasicChromatogram& c){
  vData.clear();
  charge=c.charge;
  precursorMZ=c.precursorMZ;
  precursorOffsetLower=c.precursorOffsetLower;
  precursorOffsetUpper=c.precursorOffsetUpper;
  productMZ=c.productMZ;
  productOffsetLower=c.productOffsetLower;
  productOffsetUpper=c.productOffsetUpper;
  for(unsigned int i=0;i<c.vData.size();i++) vData.push_back(c.vData[i]);
  strcpy(idString,c.idString);
}

BasicChromatogram::~BasicChromatogram(){}

BasicChromatogram& BasicChromatogram::operator=(const BasicChromatogram& c){
  if(this != &c){
    vData.clear();
    charge = c.charge;
    precursorMZ = c.precursorMZ;
    precursorOffsetLower = c.precursorOffsetLower;
    precursorOffsetUpper = c.precursorOffsetUpper;
    productMZ = c.productMZ;
    productOffsetLower = c.productOffsetLower;
    productOffsetUpper = c.productOffsetUpper;
    for(unsigned int i=0;i<c.vData.size();i++) vData.push_back(c.vData[i]);
    strcpy(idString,c.idString);
  }
  return *this;
}
  
TimeIntensityPair& BasicChromatogram::operator[ ](const unsigned int index){ return vData[index];  }

void BasicChromatogram::addTIP(TimeIntensityPair tip){ vData.push_back(tip); }
void BasicChromatogram::clear(){
  vData.clear();
  strcpy(idString,"");
  charge=0;
  precursorMZ=0;
  precursorOffsetLower=0;
  precursorOffsetUpper=0;
  productMZ=0;
  productOffsetLower=0;
  productOffsetUpper=0;
}
int BasicChromatogram::getCharge(){return charge;}
vector<TimeIntensityPair>&  BasicChromatogram::getData() { return vData; }
int BasicChromatogram::getIDString(char* str){
  strcpy(str,idString);
  return (int)strlen(str);
}
double BasicChromatogram::getPreMZ(){ return precursorMZ;}
double BasicChromatogram::getPreOffsetLower(){ return precursorOffsetLower; }
double BasicChromatogram::getPreOffsetUpper(){ return precursorOffsetUpper; }
double BasicChromatogram::getProdMZ(){ return productMZ; }
double BasicChromatogram::getProdOffsetLower(){return productOffsetLower;}
double BasicChromatogram::getProdOffsetUpper(){ return productOffsetUpper; }
void BasicChromatogram::setIDString(char* str) { strcpy(idString,str); }
void BasicChromatogram::setPrecursor(double mz, int z, double offLow, double offHigh){
  precursorMZ=mz;
  charge=z;
  precursorOffsetLower=offLow;
  precursorOffsetUpper=offHigh;
}
void BasicChromatogram::setProduct(double mz, double offLow, double offHigh){
  productMZ = mz;
  productOffsetLower = offLow;
  productOffsetUpper = offHigh;
}
size_t BasicChromatogram::size(){  return vData.size(); }

