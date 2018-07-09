/*
RAMPface - The code is
open source under the FreeBSD License, please see LICENSE file
for detailed information.

Copyright (C) 2011, Mike Hoopmann, Institute for Systems Biology
Version 1.0, January 4, 2011.
Version 1.1, March 14, 2012.
*/

// THIS FILE HAS BEEN MODIFIED BY VADIM DEMICHEV

#include "mzParser.h"

int checkFileType(const char* fname){
  char file[256];
  char ext[256];
  char *tok;
  char preExt[256];
  unsigned int i;

  if (strlen(fname) < 4) {
    cerr << "Incomplete file name. No file loaded: " << fname << endl;
    return 0;
  }

  strcpy(ext,"");

  strcpy(file,fname);
  tok=strtok(file,".\n");
  while(tok!=NULL){
    strcpy(preExt,ext);
    strcpy(ext,tok);
    tok=strtok(NULL,".\n");
  }

  for(i=0;i<strlen(ext);i++) ext[i]=toupper(ext[i]);
  for(i=0;i<strlen(preExt);i++) preExt[i]=toupper(preExt[i]);

  if (!strcmp(ext, "MZML")) return 1;
  if(!strcmp(ext,"MSTOOLKIT")) return 1;
  if(!strcmp(ext,"MZXML")) return 2;
  if(!strcmp(ext,"MZ5")) return 5;
  if(!strcmp(ext,"GZ")) {
    if(!strcmp(preExt,"MSTOOLKIT")) return 3;
    if(!strcmp(preExt,"MZXML")) return 4;
    cerr << "Unknown .gz file. Only .mzML.gz and .mzXML.gz allowed. No file loaded." << endl;
    return 0;
  }
  cerr << "Unknown file type. No file loaded." << fname << endl;
  return 0;
}

void getPrecursor(const struct ScanHeaderStruct *scanHeader,int index,double &mz,double &monoMZ,double &intensity,int &charge,int &possibleCharges,int *&possibleChargeArray){
  int i,j,k;
  double d;
  
  if(index==0){
    mz=scanHeader->precursorMZ;
    d=scanHeader->selectionWindowLower; //these are place holders until the values are actually returned;
    d=scanHeader->selectionWindowUpper;
    monoMZ=scanHeader->precursorMonoMZ;
    intensity=scanHeader->precursorIntensity;
    charge=scanHeader->precursorCharge;
    possibleCharges=scanHeader->numPossibleCharges;
    if(possibleChargeArray!=NULL) delete[] possibleChargeArray;
    if(possibleCharges>0){
      possibleChargeArray = new int[possibleCharges];
      for(i=0;i<possibleCharges;i++) memcpy(&possibleChargeArray[i],&scanHeader->possibleCharges[i*4],sizeof(int));
    } else {
      possibleChargeArray = NULL;
    }
  } else {
    j=0;
    for(i=1;i<scanHeader->precursorCount;i++){
      memcpy(&mz,&scanHeader->additionalPrecursors[j+=8],sizeof(double));
      memcpy(&d, &scanHeader->additionalPrecursors[j += 8], sizeof(double)); //selectionWindowLower
      memcpy(&d, &scanHeader->additionalPrecursors[j += 8], sizeof(double)); //selectionWindowUpper
      memcpy(&monoMZ,&scanHeader->additionalPrecursors[j+=8],sizeof(double));
      memcpy(&intensity,&scanHeader->additionalPrecursors[j+=8],sizeof(double));
      memcpy(&charge,&scanHeader->additionalPrecursors[j+=4],sizeof(int));
      memcpy(&possibleCharges,&scanHeader->additionalPrecursors[j+=4],sizeof(int));
      if(possibleChargeArray!=NULL) delete[] possibleChargeArray;
      if(possibleCharges>0){
        possibleChargeArray = new int[possibleCharges];
        for(k=0;k<possibleCharges;k++) memcpy(&possibleChargeArray[k],&scanHeader->additionalPrecursors[j+=4],sizeof(int));
      } else {
        possibleChargeArray = NULL;
      }
      if(i==index) break;
    }
  }
}

ramp_fileoffset_t getIndexOffset(RAMPFILE *pFI){
  switch(pFI->fileType){
    case 1:
    case 3:
      return (ramp_fileoffset_t)pFI->mzML->getIndexOffset();
      break;
    case 2:
    case 4:
      return (ramp_fileoffset_t)pFI->mzXML->getIndexOffset();
      break;
    default:
      return -1;
  }
}

InstrumentStruct* getInstrumentStruct(RAMPFILE *pFI){
  InstrumentStruct* r=(InstrumentStruct *) calloc(1,sizeof(InstrumentStruct));
  if(r==NULL) {
    printf("Cannot allocate memory\n");
    return NULL;
  } else {
    strcpy(r->analyzer,"UNKNOWN");
    strcpy(r->detector,"UNKNOWN");
    strcpy(r->ionisation,"UNKNOWN");
    strcpy(r->manufacturer,"UNKNOWN");
    strcpy(r->model,"UNKNOWN");
  }

  switch(pFI->fileType){
    case 1:
    case 3:
      if(pFI->mzML->getInstrument()->size()>0){
        if(pFI->mzML->getInstrument()->at(0).analyzer.size()>1) strcpy(r->analyzer,&pFI->mzML->getInstrument()->at(0).analyzer[0]);
        if(pFI->mzML->getInstrument()->at(0).detector.size()>1) strcpy(r->detector,&pFI->mzML->getInstrument()->at(0).detector[0]);
        if(pFI->mzML->getInstrument()->at(0).ionization.size()>1) strcpy(r->ionisation,&pFI->mzML->getInstrument()->at(0).ionization[0]);
        if(pFI->mzML->getInstrument()->at(0).manufacturer.size()>1) strcpy(r->manufacturer,&pFI->mzML->getInstrument()->at(0).manufacturer[0]);
        if(pFI->mzML->getInstrument()->at(0).model.size()>1) strcpy(r->model,&pFI->mzML->getInstrument()->at(0).model[0]);
      }
      break;

    case 2:
    case 4:
      if(pFI->mzXML->getInstrument().analyzer.size()>1) strcpy(r->analyzer,&pFI->mzXML->getInstrument().analyzer[0]);
      if(pFI->mzXML->getInstrument().detector.size()>1) strcpy(r->detector,&pFI->mzXML->getInstrument().detector[0]);
      if(pFI->mzXML->getInstrument().ionization.size()>1) strcpy(r->ionisation,&pFI->mzXML->getInstrument().ionization[0]);
      if(pFI->mzXML->getInstrument().manufacturer.size()>1) strcpy(r->manufacturer,&pFI->mzXML->getInstrument().manufacturer[0]);
      if(pFI->mzXML->getInstrument().model.size()>1) strcpy(r->model,&pFI->mzXML->getInstrument().model[0]);
      break;

    case 5:
    default:
      break;
  }

  return r;
}

void getScanSpanRange(const struct ScanHeaderStruct *scanHeader, int *startScanNum, int *endScanNum) {
   if (0 == scanHeader->mergedResultStartScanNum || 0 == scanHeader->mergedResultEndScanNum) {
      *startScanNum = scanHeader->acquisitionNum;
      *endScanNum = scanHeader->acquisitionNum;
   } else {
      *startScanNum = scanHeader->mergedResultStartScanNum;
      *endScanNum = scanHeader->mergedResultEndScanNum;
   }
}

void rampCloseFile(RAMPFILE *pFI){
  if(pFI!=NULL) {
    delete pFI;
    pFI=NULL;
  }
}

string rampConstructInputFileName(const string &basename){
  int len;
  char *buf = new char[len = (int)(basename.length()+100)]; 
  rampConstructInputPath(buf, len, "", basename.c_str());
  string result(buf);
  delete[] buf;
  return result;
}

char* rampConstructInputFileName(char *buf,int buflen,const char *basename){
  return rampConstructInputPath(buf, buflen, "", basename);
}

char* rampConstructInputPath(char *buf, int inbuflen, const char *dir_in, const char *basename){

  if( (int)(strlen(dir_in)+strlen(basename)+1) > inbuflen ){
    //Can't output error messages in TPP software that use this function
    //printf("rampConstructInputPath error: buffer too small for file\n");
    return NULL;
  }

  FILE* f;
  char* result = NULL;
  char base[512];
  strcpy(base,basename);

  //Try opening the base name first, then with directory:
  for(int j=0;j<2;j++){
    for(int i=0;i<5;i++){

      if(j==1){
        strcpy(buf,dir_in);
        strcat(buf,"/");
        strcat(buf,base);
      } else {
        strcpy(buf,base);
      }

      switch(i){
        case 0:  strcat(buf,".mzML");  break;
        case 1:  strcat(buf,".mzXML"); break;
        case 2:  strcat(buf,".mzML.gz"); break;
        case 3:  strcat(buf,".mzXML.gz");  break;
        case 4: strcat(buf,".mz5"); break;
        default: break;
      }
      
      f=fopen(buf,"r");
      if(f==NULL) continue;
        
      fclose(f);
      result=buf;
      return result;
    }
  }

  //Can't output error messages in TPP software that use this function
  //printf("rampConstructInputPath: file not found\n");
  strcpy(buf,"");
  result=NULL;
  return result;

}


const char** rampListSupportedFileTypes(){
  if (!data_Ext.size()) { // needs init
    data_Ext.push_back(".mzXML");
    data_Ext.push_back(".mzML");
    data_Ext.push_back(".mzXML.gz");
    data_Ext.push_back(".mzML.gz");
    data_Ext.push_back(".mz5");
    data_Ext.push_back(NULL); // end of list
  }
  return &(data_Ext[0]);
}

RAMPFILE* rampOpenFile(const char* filename){
  int i=checkFileType(filename);
  if(i==0){
    return NULL;
  } else {
    RAMPFILE* r=new RAMPFILE();
    r->bs = new BasicSpectrum();
    r->fileType=i;
    switch(i){
      case 1: //mzML
      case 3:
        r->mzML=new mzpSAXMzmlHandler(r->bs);
        if(i==3)r->mzML->setGZCompression(true);
        else r->mzML->setGZCompression(false);
        if(!r->mzML->load(filename)){
          delete r;
          return NULL;
        } else {
          return r;
        }
      case 2: //mzXML
      case 4:
        r->mzXML=new mzpSAXMzxmlHandler(r->bs);
        if(i==4) r->mzXML->setGZCompression(true);
        else r->mzXML->setGZCompression(false);
        if(!r->mzXML->load(filename)){
          delete r;
          return NULL;
        } else {
          return r;
        }
#ifdef MZP_MZ5
      case 5: //mz5
        r->mz5Config = new mzpMz5Config();
        r->mz5=new mzpMz5Handler(r->mz5Config, r->bs);
        if(!r->mz5->readFile(filename)){
          delete r;
          return NULL;
        } else {
          return r;
        }
#endif
      default:
        delete r;
        return NULL;
    }
  }

}

char* rampValidFileType(const char *buf){
  char ext[256];
  char preExt[256];

  const char* result=NULL;
  const char* result2=NULL;
  const char* p;

  unsigned int i;

  p=strchr(buf,'.');
  while(p!=NULL){
    result2=result;
    result=p;
    p=strchr(p+1,'.');
  }

  if(result==NULL) return (char*) result;

  strcpy(ext,result);
  for(i=0;i<strlen(ext);i++) ext[i]=toupper(ext[i]);

  if(result2){
    strcpy(preExt,result2);
    for(i=0;i<strlen(preExt);i++) preExt[i]=toupper(preExt[i]);
  }

  if(!strcmp(ext,".MSTOOLKIT")) return (char*) result;
  if(!strcmp(ext,".MZXML")) return (char*) result;
  if(!strcmp(ext,".MZ5")) return (char*) result;
  if(!strcmp(ext,".GZ")) {
    if(!strcmp(preExt,".MSTOOLKIT.GZ")) return (char*) result2;
    if(!strcmp(preExt,".MZXML.GZ")) return (char*) result2;
    cout << "Unknown .gz file. Only .mzML.gz and .mzXML.gz allowed. No file loaded." << endl;
  }
  if(!strcmp(ext,".MZDATA")) {
    cout << ".mzData is not supported. Please convert to mz5, mzXML, or mzML." << endl;
  }
  result=NULL;
  return (char*) result;
}

//MH: Read header is redundant with readPeaks, which automatically reads the header.
//But due to legacy issues, this function must exist.
void readHeader(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex, struct ScanHeaderStruct *scanHeader){
    return readHeader(pFI, lScanIndex, scanHeader, 0);
}

void readHeader(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex, struct ScanHeaderStruct *scanHeader, unsigned long scanI){

  vector<cindex>* v;
  sPrecursorIon p;
  unsigned int i;

#ifdef MZP_MZ5
  vector<cMz5Index>* v2;
#endif

  //memset(scanHeader,0,sizeof(struct ScanHeaderStruct));
  scanHeader->acquisitionNum=-1;
  scanHeader->activationMethod[0]='\0';
  scanHeader->basePeakIntensity=0.0;
  scanHeader->basePeakMZ=0.0;
  scanHeader->centroid=false;
  scanHeader->collisionEnergy=0.0;
  scanHeader->compensationVoltage=0.0;
  scanHeader->filePosition=0;
  scanHeader->filterLine[0]='\0';
  scanHeader->highMZ=0.0;
  scanHeader->idString[0]='\0';
  scanHeader->ionisationEnergy=0.0;
  scanHeader->lowMZ=0.0;
  scanHeader->mergedScan=0;
  scanHeader->mergedResultScanNum=0;
  scanHeader->mergedResultStartScanNum=0;
  scanHeader->mergedResultEndScanNum=0;
  scanHeader->msLevel=0;
  scanHeader->numPossibleCharges=0;
  scanHeader->precursorCharge=0;
  scanHeader->precursorIntensity=0.0;
  scanHeader->precursorMonoMZ=0.0;
  scanHeader->precursorMZ=0.0;
  scanHeader->precursorScanNum=-1;
  scanHeader->retentionTime=0.0;
  scanHeader->scanType[0]='\0';
  scanHeader->selectionWindowLower=0;
  scanHeader->selectionWindowUpper=0;
  scanHeader->totIonCurrent=0.0;
  scanHeader->scanIndex=0;
  scanHeader->seqNum=-1;

  if(lScanIndex<0) return;
  
  switch(pFI->fileType){
    case 1:
    case 3:
      v=pFI->mzML->getSpecIndex();
      if(!scanI){
        for(i=0;i<v->size();i++) {
          if(v->at(i).offset==(f_off)lScanIndex) {
            if(!pFI->mzML->readHeader(v->at(i).scanNum)){
              v=NULL;
              return;
            }
            break;
          }
        }
      } else {
        if(v->at(scanI).offset==(f_off)lScanIndex){
          if(!pFI->mzML->readHeader(v->at(scanI).scanNum)){
            v=NULL;
            return;
          }
        }
      }
      break;
    case 2:
    case 4:
      v=pFI->mzXML->getIndex();
      if(!scanI){
        for(i=0;i<v->size();i++) {
          if(v->at(i).offset==(f_off)lScanIndex) {
            if(!pFI->mzXML->readHeader(v->at(i).scanNum)){
              v=NULL;
              return;
            }
            break;
          }
        }
      } else {
        if(v->at(scanI).offset==(f_off)lScanIndex){
          if (!pFI->mzXML->readHeader(v->at(scanI).scanNum)){
            v = NULL;
            return;
          }
        }
      }
      break;
#ifdef MZP_MZ5
    case 5:
      v2=pFI->mz5->getSpecIndex();
      if(!scanI){
        for(i=0;i<v2->size();i++) {
          if(v2->at(i).offset==(f_off)lScanIndex) {
            if(!pFI->mz5->readHeader(v2->at(i).scanNum)){
              v2=NULL;
              return;
            }
            break;
          }
        }
      } else {
        if (v2->at(scanI).offset == (f_off)lScanIndex) {
          if (!pFI->mz5->readHeader(v2->at(scanI).scanNum)){
            v2 = NULL;
            return;
          }
          break;
        }
      }
      break;
#endif
    default:
      pFI->bs->clear();
      v=NULL;
#ifdef MZP_MZ5
      v2=NULL;
#endif
      return;
  }
  v=NULL;
#ifdef MZP_MZ5
  v2=NULL;
#endif

  scanHeader->acquisitionNum=pFI->bs->getScanNum();
  scanHeader->basePeakIntensity=pFI->bs->getBasePeakIntensity();
  scanHeader->basePeakMZ=pFI->bs->getBasePeakMZ();
  scanHeader->centroid=pFI->bs->getCentroid();
  scanHeader->collisionEnergy=pFI->bs->getCollisionEnergy();
  scanHeader->compensationVoltage=pFI->bs->getCompensationVoltage();
  scanHeader->highMZ=pFI->bs->getHighMZ();
  scanHeader->ionInjectionTime = pFI->bs->getIonInjectionTime();
  scanHeader->lowMZ=pFI->bs->getLowMZ();
  scanHeader->msLevel=pFI->bs->getMSLevel();
  scanHeader->peaksCount=pFI->bs->getPeaksCount();
  scanHeader->precursorScanNum=pFI->bs->getPrecursorScanNum();
  scanHeader->retentionTime=(double)pFI->bs->getRTime(false);
  scanHeader->totIonCurrent=pFI->bs->getTotalIonCurrent();
  scanHeader->scanIndex=pFI->bs->getScanIndex();
  scanHeader->seqNum=pFI->bs->getScanIndex();

  int j=0;
  int k;
  double d;
  scanHeader->precursorCount=pFI->bs->getPrecursorIonCount();
  for(i=0;i<(unsigned int)scanHeader->precursorCount;i++){
    p=pFI->bs->getPrecursorIon(i);
    if(i==0){
      scanHeader->precursorCharge=p.charge;
      scanHeader->precursorIntensity=p.intensity;
      scanHeader->precursorMonoMZ=p.monoMZ;
      if(p.isoMZ==0) scanHeader->precursorMZ=p.mz;
      else scanHeader->precursorMZ=p.isoMZ;
      scanHeader->selectionWindowLower=p.isoMZ-p.isoLowerMZ;
      scanHeader->selectionWindowUpper=p.isoMZ+p.isoUpperMZ;
      scanHeader->numPossibleCharges=(int)p.possibleCharges->size();
      for(k=0;k<scanHeader->numPossibleCharges;k++){
        memcpy(&scanHeader->possibleCharges[k*4],&p.possibleCharges->at(k),sizeof(int));
        if(k==7){
          cout << "Warning: too many possible charges for precursor in scan " << scanHeader->acquisitionNum << endl;
          break;
        }
      }
    } else {
      if( (j+32+p.possibleCharges->size()*4) > PRECURSORARRAY_LENGTH-1) {
        cout << "Warning: too many precursors. Must improve RAMP interface." << endl;
        break;
      }
      if (p.isoMZ == 0) memcpy(&scanHeader->additionalPrecursors[j += 8], &p.mz, sizeof(double));
      else memcpy(&scanHeader->additionalPrecursors[j+=8],&p.isoMZ,sizeof(double));
      d = p.isoMZ - p.isoLowerMZ;
      memcpy(&scanHeader->additionalPrecursors[j+=8],&d, sizeof(double));
      d = p.isoMZ + p.isoUpperMZ;
      memcpy(&scanHeader->additionalPrecursors[j+=8],&d, sizeof(double));
      memcpy(&scanHeader->additionalPrecursors[j+=8],&p.monoMZ,sizeof(double));
      memcpy(&scanHeader->additionalPrecursors[j+=8],&p.intensity,sizeof(double));
      memcpy(&scanHeader->additionalPrecursors[j+=4],&p.charge,sizeof(int));
      k=(int)p.possibleCharges->size();
      memcpy(&scanHeader->additionalPrecursors[j+=4],&k,sizeof(int));
      for(k=0;k<(int)p.possibleCharges->size();k++){
        memcpy(&scanHeader->additionalPrecursors[j+=k*4],&p.possibleCharges->at(k),sizeof(int));
      }
    }
  }
  
  pFI->bs->getFilterLine(scanHeader->filterLine);
  pFI->bs->getIDString(scanHeader->idString);

  switch(pFI->bs->getActivation()){
    case 1: strcpy(scanHeader->activationMethod,"CID"); break;
    case 2: strcpy(scanHeader->activationMethod,"HCD"); break;
    case 3: strcpy(scanHeader->activationMethod,"ETD"); break;
    case 4: strcpy(scanHeader->activationMethod,"ETD+SA"); break;
    case 5: strcpy(scanHeader->activationMethod,"ECD"); break;
    default: strcpy(scanHeader->activationMethod,""); break;
  }

}

//MH: Indexes in RAMP are stored in an array indexed by scan number, with -1 for the offset
//if the scan number does not exist.
ramp_fileoffset_t* readIndex(RAMPFILE *pFI, ramp_fileoffset_t indexOffset, int *iLastScan){
  vector<cindex>* v;
#ifdef MZP_MZ5
  vector<cMz5Index>* v2;
#endif
  ramp_fileoffset_t* rIndex;
  unsigned int i;
  switch(pFI->fileType){
    case 1:
    case 3:
      v=pFI->mzML->getSpecIndex();
      rIndex = (ramp_fileoffset_t *) malloc((pFI->mzML->highScan()+2)*sizeof(ramp_fileoffset_t));
      memset(rIndex,-1,(pFI->mzML->highScan()+2)*sizeof(ramp_fileoffset_t));
      for(i=0;i<v->size();i++) rIndex[v->at(i).scanNum]=(ramp_fileoffset_t)v->at(i).offset;
      rIndex[v->at(i-1).scanNum+1]=-1;
      *iLastScan=v->at(i-1).scanNum;
      break;
    case 2:
    case 4:
      v=pFI->mzXML->getIndex();
      rIndex = (ramp_fileoffset_t *) malloc((pFI->mzXML->highScan()+2)*sizeof(ramp_fileoffset_t));
      memset(rIndex,-1,(pFI->mzXML->highScan()+2)*sizeof(ramp_fileoffset_t));
      for(i=0;i<v->size();i++) rIndex[v->at(i).scanNum]=(ramp_fileoffset_t)v->at(i).offset;
      rIndex[v->at(i-1).scanNum+1]=-1;
      *iLastScan=v->at(i-1).scanNum;
      break;
#ifdef MZP_MZ5
    case 5:
      v2=pFI->mz5->getSpecIndex();
      rIndex = (ramp_fileoffset_t *) malloc((pFI->mz5->highScan()+2)*sizeof(ramp_fileoffset_t));
      memset(rIndex,-1,(pFI->mz5->highScan()+2)*sizeof(ramp_fileoffset_t));
      for(i=0;i<v2->size();i++) rIndex[v2->at(i).scanNum]=(ramp_fileoffset_t)v2->at(i).offset;
      rIndex[v2->at(i-1).scanNum+1]=-1;
      *iLastScan=v2->at(i-1).scanNum;
      break;
#endif
    default:
      rIndex=NULL;
      *iLastScan=0;
      break;
  }
  v=NULL;
#ifdef MZP_MZ5
  v2=NULL;
#endif
  return rIndex;
}

int readMsLevel(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex){
  vector<cindex>* v;
#ifdef MZP_MZ5
  vector<cMz5Index>* v2;
#endif
  unsigned int i;
  
  if(lScanIndex<0) return 0;

  switch(pFI->fileType){
    case 1:
    case 3:
      v=pFI->mzML->getSpecIndex();
      for(i=0;i<v->size();i++) {
        if(v->at(i).offset==(f_off)lScanIndex) {
          pFI->mzML->readSpectrum(v->at(i).scanNum);
          break;
        }
      }
      break;
    case 2:
    case 4:
      v=pFI->mzXML->getIndex();
      for(i=0;i<v->size();i++) {
        if(v->at(i).offset==(f_off)lScanIndex) {
          pFI->mzXML->readSpectrum(v->at(i).scanNum);
          break;
        }
      }
      break;
#ifdef MZP_MZ5
    case 5:
      v2=pFI->mz5->getSpecIndex();
      for(i=0;i<v2->size();i++) {
        if(v2->at(i).offset==(f_off)lScanIndex) {
          pFI->mz5->readSpectrum(v2->at(i).scanNum);
          break;
        }
      }
      break;
#endif
    default:
      pFI->bs->clear();
      break;
  }
  v=NULL;
#ifdef MZP_MZ5
  v2=NULL;
#endif

  return pFI->bs->getMSLevel();
}

void readMSRun(RAMPFILE *pFI, struct RunHeaderStruct *runHeader){

  vector<cindex>* v;
#ifdef MZP_MZ5
  vector<cMz5Index>* v2;
#endif

  //memset(scanHeader,0,sizeof(struct ScanHeaderStruct));
  runHeader->dEndTime=0.0;
  runHeader->dStartTime=0.0;
  runHeader->endMZ=0.0;
  runHeader->highMZ=0.0;
  runHeader->lowMZ=0.0;
  runHeader->scanCount=0;
  runHeader->startMZ=0.0;
  
  switch(pFI->fileType){
    case 1:
    case 3:
      v=pFI->mzML->getSpecIndex();
      runHeader->scanCount=(int)v->size();
      pFI->mzML->readHeader(v->at(0).scanNum);
      runHeader->dStartTime=pFI->bs->getRTime(false);
      pFI->mzML->readHeader(v->at(v->size()-1).scanNum);
      runHeader->dEndTime=pFI->bs->getRTime(false);
      pFI->bs->clear();
      v=NULL;
      break;
    case 2:
    case 4:
      v=pFI->mzXML->getIndex();
      runHeader->scanCount=(int)v->size();
      pFI->mzXML->readHeader(v->at(0).scanNum);
      runHeader->dStartTime=pFI->bs->getRTime(false);
      pFI->mzXML->readHeader(v->at(v->size()-1).scanNum);
      runHeader->dEndTime=pFI->bs->getRTime(false);
      pFI->bs->clear();
      v=NULL;
      break;
#ifdef MZP_MZ5
    case 5:
      v2=pFI->mz5->getSpecIndex();
      runHeader->scanCount=v2->size();
      pFI->mz5->readHeader(v2->at(0).scanNum);
      runHeader->dStartTime=pFI->bs->getRTime(false);
      pFI->mz5->readHeader(v2->at(v2->size()-1).scanNum);
      runHeader->dEndTime=pFI->bs->getRTime(false);
      pFI->bs->clear();
      v2=NULL;
      break;
#endif
    default:
      break;
  }

}

//MH: Matching the index is very indirect, but requires less code,
//making this wrapper much easier to read
RAMPREAL* readPeaks(RAMPFILE* pFI, ramp_fileoffset_t lScanIndex){
  return  readPeaks(pFI, lScanIndex, 0);
}

RAMPREAL* readPeaks(RAMPFILE* pFI, ramp_fileoffset_t lScanIndex, unsigned long scanI) {
  vector<cindex>* v;
#ifdef MZP_MZ5
  vector<cMz5Index>* v2;
#endif
  unsigned int i;
  RAMPREAL* pPeaks=NULL;

  if(lScanIndex<0) return pPeaks;

  switch(pFI->fileType){
    case 1:
    case 3:
      v=pFI->mzML->getSpecIndex();
      if(!scanI){
        for(i=0;i<v->size();i++) {
          if(v->at(i).offset==(f_off)lScanIndex) {
            pFI->mzML->readSpectrum(v->at(i).scanNum);
            break;
          }
        }
      } else {
        if (v->at(scanI).offset == (f_off)lScanIndex) {
          pFI->mzML->readSpectrum(v->at(scanI).scanNum);
          break;
        }
      }
      break;
    case 2:
    case 4:
      v=pFI->mzXML->getIndex();
      if(!scanI){
        for(i=0;i<v->size();i++) {
          if(v->at(i).offset==(f_off)lScanIndex) {
            pFI->mzXML->readSpectrum(v->at(i).scanNum);
            break;
          }
        }
      } else {
        if (v->at(scanI).offset == (f_off)lScanIndex) {
          pFI->mzXML->readSpectrum(v->at(scanI).scanNum);
          break;
        }
      }
      break;
#ifdef MZP_MZ5
    case 5:
      v2=pFI->mz5->getSpecIndex();
      if(!scanI){
        for(i=0;i<v2->size();i++) {
          if(v2->at(i).offset==(f_off)lScanIndex) {
            pFI->mz5->readSpectrum(v2->at(i).scanNum);
            break;
          }
        }
      } else {
        if (v2->at(scanI).offset == (f_off)lScanIndex) {
          pFI->mz5->readSpectrum(v2->at(scanI).scanNum);
          break;
        }
      }
      break;
#endif
    default:
      pFI->bs->clear();
      break;
  }
  v=NULL;
#ifdef MZP_MZ5
  v2=NULL;
#endif

  unsigned int j=0;
  if(pFI->bs->size()>0){
    pPeaks = (RAMPREAL *) malloc((pFI->bs->size()+1) * 2 * sizeof(RAMPREAL) + 1);
    for(i=0;i<pFI->bs->size();i++){
      pPeaks[j++]=pFI->bs->operator [](i).mz;
      pPeaks[j++]=pFI->bs->operator [](i).intensity;
    }
  } else {
    pPeaks = (RAMPREAL *) malloc(2 * sizeof(RAMPREAL));
  }
  pPeaks[j]=-1;

  return pPeaks;
}

int readPeaksCount(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex){
  return readPeaksCount(pFI, lScanIndex, 0);
}

int readPeaksCount(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex, unsigned long scanI){
  ScanHeaderStruct s;
  readHeader(pFI, lScanIndex, &s, scanI);
  return s.peaksCount;
}

void readRunHeader(RAMPFILE *pFI, ramp_fileoffset_t *pScanIndex, struct RunHeaderStruct *runHeader, int iLastScan){
  vector<cindex>* v;
#ifdef MZP_MZ5
  vector<cMz5Index>* v2;
#endif
  unsigned int i;

  runHeader->scanCount=0;
  runHeader->dEndTime=0.0;
  runHeader->dStartTime=0.0;
  runHeader->endMZ=0.0;
  runHeader->highMZ=0.0;
  runHeader->lowMZ=0.0;
  runHeader->startMZ=0.0;
  
  switch(pFI->fileType){
    case 1:
    case 3:
      v=pFI->mzML->getSpecIndex();
      runHeader->scanCount=(int)v->size();
      
      pFI->mzML->readHeader(v->at(0).scanNum);
      runHeader->dStartTime=(double)pFI->bs->getRTime(false);
      runHeader->lowMZ=pFI->bs->getLowMZ();
      runHeader->highMZ=pFI->bs->getHighMZ();
      runHeader->startMZ=runHeader->lowMZ;
      runHeader->endMZ=runHeader->highMZ;
      
      for(i=1;i<v->size();i++) {
        pFI->mzML->readHeader(v->at(i).scanNum);
        if(pFI->bs->getLowMZ()<runHeader->lowMZ) {
          runHeader->lowMZ=pFI->bs->getLowMZ();
          runHeader->startMZ=runHeader->lowMZ;
        }
        if(pFI->bs->getHighMZ()>runHeader->highMZ){
          runHeader->highMZ=pFI->bs->getHighMZ();
          runHeader->endMZ=runHeader->highMZ;
        }
      }
      pFI->mzML->readHeader(v->at(v->size()-1).scanNum);
      break;

    case 2:
    case 4:
      v=pFI->mzXML->getIndex();
      runHeader->scanCount=(int)v->size();
      
      pFI->mzXML->readHeader(v->at(0).scanNum);
      runHeader->dStartTime=(double)pFI->bs->getRTime(false);
      runHeader->lowMZ=pFI->bs->getLowMZ();
      runHeader->highMZ=pFI->bs->getHighMZ();
      runHeader->startMZ=runHeader->lowMZ;
      runHeader->endMZ=runHeader->highMZ;
      
      for(i=1;i<v->size();i++) {
        pFI->mzXML->readHeader(v->at(i).scanNum);
        if(pFI->bs->getLowMZ()<runHeader->lowMZ) {
          runHeader->lowMZ=pFI->bs->getLowMZ();
          runHeader->startMZ=runHeader->lowMZ;
        }
        if(pFI->bs->getHighMZ()>runHeader->highMZ){
          runHeader->highMZ=pFI->bs->getHighMZ();
          runHeader->endMZ=runHeader->highMZ;
        }
      }
      pFI->mzXML->readHeader(v->at(v->size()-1).scanNum);
      break;

#ifdef MZP_MZ5
    case 5:
      v2=pFI->mz5->getSpecIndex();
      runHeader->scanCount=v2->size();
      
      pFI->mz5->readHeader(v2->at(0).scanNum);
      runHeader->dStartTime=(double)pFI->bs->getRTime(false);
      runHeader->lowMZ=pFI->bs->getLowMZ();
      runHeader->highMZ=pFI->bs->getHighMZ();
      runHeader->startMZ=runHeader->lowMZ;
      runHeader->endMZ=runHeader->highMZ;
      
      for(i=1;i<v2->size();i++) {
        pFI->mz5->readHeader(v2->at(i).scanNum);
        if(pFI->bs->getLowMZ()<runHeader->lowMZ) {
          runHeader->lowMZ=pFI->bs->getLowMZ();
          runHeader->startMZ=runHeader->lowMZ;
        }
        if(pFI->bs->getHighMZ()>runHeader->highMZ){
          runHeader->highMZ=pFI->bs->getHighMZ();
          runHeader->endMZ=runHeader->highMZ;
        }
      }
      pFI->mz5->readHeader(v2->at(v2->size()-1).scanNum);
      break;
#endif

    default:
      pFI->bs->clear();
      v=NULL;
#ifdef MZP_MZ5
      v2=NULL;
#endif
      return;
  }
  v=NULL;
#ifdef MZP_MZ5
  v2=NULL;
#endif

}




//--------------------------------------------------
// CACHED RAMP FUNCTIONS
//--------------------------------------------------
void clearScanCache(struct ScanCacheStruct* cache){
  for (int i=0; i<cache->size; i++) {
    if (cache->peaks[i] == NULL) continue;
    free(cache->peaks[i]);
    cache->peaks[i] = NULL;
  }
  memset(cache->headers, 0, cache->size * sizeof(struct ScanHeaderStruct));
}

void freeScanCache(struct ScanCacheStruct* cache){
  if (cache) {
    for (int i=0; i<cache->size; i++){
      if (cache->peaks[i] != NULL) free(cache->peaks[i]);
    }
    free(cache->peaks);
    free(cache->headers);
    free(cache);
  }
}

int getCacheIndex(struct ScanCacheStruct* cache, int seqNum) {
  int seqNumStart = cache->seqNumStart;
  int size = cache->size;

  // First access, just set the start to seqNum.
  if (seqNumStart == 0) cache->seqNumStart = seqNum;

  // If requested scan is less than cache start, shift cache window
  // left to start at requested scan.
  else if (seqNum < seqNumStart) shiftScanCache(cache, (int) (seqNum - seqNumStart));

  // If requested scan is greater than cache end, shift cache window
  // right so last entry is requested scan.
  else if (seqNum >= seqNumStart + size) shiftScanCache(cache, (int) (seqNum - (seqNumStart + size - 1)));

  return (int) (seqNum - cache->seqNumStart);
}

struct ScanCacheStruct* getScanCache(int size){
  struct ScanCacheStruct* cache = (struct ScanCacheStruct*) malloc(sizeof(struct ScanCacheStruct));
  cache->seqNumStart = 0;
  cache->size = size;
  cache->headers = (struct ScanHeaderStruct*) calloc(size, sizeof(struct ScanHeaderStruct));
  cache->peaks = (RAMPREAL**) calloc(size, sizeof(RAMPREAL*));
  return cache;
}

const struct ScanHeaderStruct* readHeaderCached(struct ScanCacheStruct* cache, int seqNum, RAMPFILE* pFI, ramp_fileoffset_t lScanIndex){
  int i = getCacheIndex(cache, seqNum);
  if (cache->headers[i].msLevel == 0) readHeader(pFI, lScanIndex, cache->headers + i);
  return cache->headers + i;
}

int  readMsLevelCached(struct ScanCacheStruct* cache, int seqNum, RAMPFILE* pFI, ramp_fileoffset_t lScanIndex){
  const struct ScanHeaderStruct* header = readHeaderCached(cache, seqNum, pFI, lScanIndex);
  return header->msLevel;
}

const RAMPREAL* readPeaksCached(struct ScanCacheStruct* cache, int seqNum, RAMPFILE* pFI, ramp_fileoffset_t lScanIndex){
  int i = getCacheIndex(cache, seqNum);
  if (cache->peaks[i] == NULL) cache->peaks[i] = readPeaks(pFI, lScanIndex);
  return cache->peaks[i];
}

void shiftScanCache(struct ScanCacheStruct* cache, int nScans) {
  int i;
  cache->seqNumStart += nScans;
  if (abs(nScans) > cache->size) {
    // If the shift is larger than the size of the cache window,
    // just clear the whole cache.
    clearScanCache(cache);
  } else if (nScans > 0) {
    // Shifting window to the right.  Memory moves right, with new
    // empty scans on the end.

    // Free the peaks that memmove will overwrite.
    for (i = 0; i < nScans; i++) {
      if (cache->peaks[i] != NULL) free(cache->peaks[i]);
    }
    memmove(cache->peaks, cache->peaks + nScans, (cache->size - nScans) * sizeof(RAMPREAL*));
    memset(cache->peaks + cache->size - nScans, 0, nScans * sizeof(RAMPREAL*));
    memmove(cache->headers, cache->headers + nScans,(cache->size - nScans) * sizeof(struct ScanHeaderStruct));
    memset(cache->headers + cache->size - nScans, 0, nScans * sizeof(struct ScanHeaderStruct));
  } else if (nScans < 0) {
    // Shifting window to the left.  Memory moves right, with new
    // empty scans at the beginning.
    nScans = -nScans;

    // Free the peaks that memmove will overwrite.
    for (i = 0; i < nScans; i++) {
      if (cache->peaks[cache->size - 1 - i] != NULL) free(cache->peaks[cache->size - 1 - i]);
    }
    memmove(cache->peaks + nScans, cache->peaks, (cache->size - nScans) * sizeof(RAMPREAL*));
    memset(cache->peaks, 0, nScans * sizeof(RAMPREAL*));
    memmove(cache->headers  + nScans, cache->headers, (cache->size - nScans) * sizeof(struct ScanHeaderStruct));
    memset(cache->headers, 0, nScans * sizeof(struct ScanHeaderStruct));
  }
}



//--------------------------------------------------
// DEAD FUNCTIONS
//--------------------------------------------------
int isScanAveraged(struct ScanHeaderStruct *scanHeader){
  cerr << "call to unsupported function: isScanAveraged(struct ScanHeaderStruct *scanHeader)" << endl;
  return 0;
}

int isScanMergedResult(struct ScanHeaderStruct *scanHeader){
  cerr << "call to unsupported function: isScanMergedResult(struct ScanHeaderStruct *scanHeader)" << endl;
  return 0;
}

int rampSelfTest(char *filename){
  cerr << "call to unsupported function: rampSelfTest(char *filename)" << endl;
  return 0;
}

char* rampTrimBaseName(char *buf){
  cerr << "call to unsupported function: rampTrimBaseName(char *buf)" << endl;
  return buf;
}

int rampValidateOrDeriveInputFilename(char *inbuf, int inbuflen, char *spectrumName){
  cerr << "call to unsupported function: rampValidateOrDeriveInputFilename(char *inbuf, int inbuflen, char *spectrumName)" << endl;
  return 0;
}

double readStartMz(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex){
  cerr << "call to unsupported function: readStartMz(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex)" << endl;
  return 0.0;
}

double readEndMz(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex){
  cerr << "call to unsupported function: readEndMz(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex)" << endl;
  return 0.0;
}

void setRampOption(long option){
  cerr << "call to unsupported function: setRampOption(long option)" << endl;
}
