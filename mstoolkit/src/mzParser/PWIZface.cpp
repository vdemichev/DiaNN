/*
PWIZface - The code is
open source under the FreeBSD License, please see LICENSE file
for detailed information.

Copyright (C) 2011, Mike Hoopmann, Institute for Systems Biology
Version 1.0, January 4, 2011.
Version 1.1, March 14, 2012.
*/
#include "mzParser.h"

Chromatogram::Chromatogram(){
  bc = NULL;
}

Chromatogram::~Chromatogram(){
  bc = NULL;
}

void Chromatogram::getTimeIntensityPairs(vector<TimeIntensityPair>& v){
  if(bc==NULL) cerr << "Null chromatogram" << endl;
  else v=bc->getData();
}

ChromatogramList::ChromatogramList(){
}

ChromatogramList::ChromatogramList(mzpSAXMzmlHandler* ml, void* m5, BasicChromatogram* bc){
  mzML=ml;
  #ifdef MZP_MZ5
    mz5=(mzpMz5Handler*)m5;
  #endif
  chromat=new Chromatogram();
  chromat->bc=bc;
}

ChromatogramList::~ChromatogramList(){
  mzML=NULL;
  vChromatIndex=NULL;
  #ifdef MZP_MZ5
  mz5=NULL;
  vMz5Index=NULL;
  #endif
  delete chromat;
}

ChromatogramPtr ChromatogramList::chromatogram(int index, bool binaryData) {
  char str[128];
  if(mzML!=NULL) {
    mzML->readChromatogram(index);
    chromat->bc->getIDString(str);
    chromat->id=str;
    return chromat;
  #ifdef MZP_MZ5
  }  else if(mz5!=NULL) {
    mz5->readChromatogram(index);
    chromat->bc->getIDString(str);
    chromat->id=str;
    return chromat;
  #endif
  }
  return NULL;
}
bool ChromatogramList::get() {
  if(mzML!=NULL) vChromatIndex = mzML->getChromatIndex();
  #ifdef MZP_MZ5
  else if(mz5!=NULL) vMz5Index = mz5->getChromatIndex();
  #endif
  else return false;
  return true;
}

size_t ChromatogramList::size() {
  if(vChromatIndex==NULL) {
    cerr << "Get chromatogram list first." << endl;
    return 0;
  }
  if(mzML!=NULL) return vChromatIndex->size();
  #ifdef MZP_MZ5
  else if(mz5!=NULL) return vMz5Index->size();
  #endif
  else return 0;
}

PwizRun::PwizRun(){
  chromatogramListPtr = new ChromatogramList();
}

PwizRun::PwizRun(mzpSAXMzmlHandler* ml, void* m5, BasicChromatogram* b){
  mzML=ml;
#ifdef MZP_MZ5
  mz5=(mzpMz5Handler*)m5;
#endif
  bc=b;
  chromatogramListPtr = new ChromatogramList(ml, m5, b);
}

PwizRun::~PwizRun(){
  mzML=NULL;
#ifdef MZP_MZ5
  mz5=NULL;
#endif
  bc=NULL;
  delete chromatogramListPtr;
}

void PwizRun::set(mzpSAXMzmlHandler* ml, void* m5, BasicChromatogram* b){
  mzML=ml;
#ifdef MZP_MZ5
  mz5=(mzpMz5Handler*)m5;
#endif
  bc=b;
  delete chromatogramListPtr;
  chromatogramListPtr = new ChromatogramList(ml, m5, b);
}

MSDataFile::MSDataFile(string s){
  int i=checkFileType(&s[0]);
  if(i==0){
    cerr << "Cannot identify file type." << endl;
  } else {
    bs = new BasicSpectrum();
    bc = new BasicChromatogram();
    switch(i){
      case 1: //mzML
      case 3:
        mzML=new mzpSAXMzmlHandler(bs,bc);
        if(i==3) mzML->setGZCompression(true);
        else mzML->setGZCompression(false);
        if(!mzML->load(&s[0])){
          cerr << "Failed to load file." << endl;
          delete mzML;
          delete bs;
          delete bc;
        }
        run.chromatogramListPtr->vChromatIndex=mzML->getChromatIndex();
        break;
      case 2: //mzXML
      case 4:
        cerr << "mzXML not supported in this interface." << endl;
        delete bs;
        delete bc;
        break;
#ifdef MZP_MZ5
      case 5: //mz5
        mz5Config = new mzpMz5Config();
        mz5=new mzpMz5Handler(mz5Config, bs, bc);
        if(!mz5->readFile(&s[0])){
          cerr << "Failed to load file." << endl;
          delete mz5;
          delete mz5Config;
          delete bs;
          delete bc;
        }
        break;
#endif
      default:
        break;
    }
#ifdef MZP_MZ5
    run.set(mzML,mz5,bc);
#else
    run.set(mzML,NULL,bc);
#endif
  }
}

MSDataFile::~MSDataFile(){
  if(mzML!=NULL) delete mzML;
#ifdef MZP_MZ5
  if(mz5!=NULL){
    delete mz5;
    delete mz5Config;
  }
#endif
  delete bs;
  delete bc;
}
