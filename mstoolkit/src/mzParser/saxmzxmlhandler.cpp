/************************************************************
 *              SAXMzxmlHandler.cpp
 * Adapted from the SAXMzmlHandler.cpp fork at ISB
 * Originally adapted from SAXMzdataHandler.cpp
 * August 2008
 * Ronald Beavis
 *
 * December 2010 - Drastically modified and cannibalized to create
 * robust, yet generic, mzXML parser
 * Mike Hoopmann, Institute for Systems Biology
 *
 * See http://sashimi.sourceforge.net/schema_revision/mzXML_3.2/ for
 * mzXML schema information.
 *
 * Inspired by DtaSAX2Handler.cpp
 * copyright            : (C) 2002 by Pedrioli Patrick, ISB, Proteomics
 * email                : ppatrick@systemsbiology.org
 * Artistic License granted 3/11/2005
 *******************************************************/

#include "mzParser.h"



mzpSAXMzxmlHandler::mzpSAXMzxmlHandler(BasicSpectrum* bs){
  m_bInMsInstrument=false;
  m_bInDataProcessing=false;
  m_bInScan=false;
  m_bInPrecursorMz=false;
  m_bInMsRun=false;
  m_bInIndex=false;
  m_bInPeaks=false;
  m_bCompressedData=false;
  m_bHeaderOnly=false;
  m_bLowPrecision=false;
  m_bNetworkData=true;
  m_bNoIndex=true;
  m_bScanIndex=false;
  m_bIndexSorted = true;
  spec=bs;
  indexOffset=-1;
}

mzpSAXMzxmlHandler::~mzpSAXMzxmlHandler(){
  spec=NULL;
}

void mzpSAXMzxmlHandler::startElement(const XML_Char *el, const XML_Char **attr){

  string s;

  if (isElement("dataProcessing",el)) {
    m_bInDataProcessing=true;

  } else if(isElement("index",el)){
    if(!strcmp(getAttrValue("name", attr),"scan")) m_bScanIndex=true;
    m_vIndex.clear();
    m_bInIndex=true;

  } else if(isElement("msDetector",el)){
    m_instrument.detector=getAttrValue("value", attr);

  } else if(isElement("msInstrument",el)){
    m_bInMsInstrument=true;
    m_instrument.clear();
    m_instrument.id=getAttrValue("id", attr);

  } else if(isElement("msIonisation",el)){
    m_instrument.ionization=getAttrValue("value", attr);

  } else if(isElement("msManufacturer",el)){
    m_instrument.manufacturer=getAttrValue("value", attr);

  } else if(isElement("msMassAnalyzer",el)){
    m_instrument.analyzer=getAttrValue("value", attr);

  } else if(isElement("msModel",el)){
    m_instrument.model=getAttrValue("value", attr);

  } else if(isElement("msRun",el)){
    m_bInMsRun=true;

  } else if(isElement("offset",el) && m_bScanIndex){
    m_strData.clear();
    curIndex.scanNum=atoi(getAttrValue("id", attr));
    curIndex.idRef="";

  } else if(isElement("peaks",el)){
    m_strData.clear();
    m_bInPeaks=true;

    if(!strcmp("network",getAttrValue("byteOrder", attr))) m_bNetworkData=true;
    else m_bNetworkData=true; //make it true anyway.
    if(!strcmp("64",getAttrValue("precision", attr))) m_bLowPrecision=false;
    else m_bLowPrecision=true;

    m_bCompressedData=false;
    s=getAttrValue("compressionType", attr);
    if(!strcmp("zlib",&s[0])) m_bCompressedData=true;
    else if(!strcmp("none",&s[0])) m_bCompressedData=false;
    else if(s.length()>0) {
      cout << "Halting! Unknown compression type: " <<  &s[0] << endl;
      exit(-5);
    }
    s=getAttrValue("compressedLen", attr);
    if(s.length()>0) m_compressLen = (uLong)atoi(&s[0]);
    else m_compressLen=0;

    if(m_bHeaderOnly) stopParser();

  }  
  else if (isElement("precursorMz", el)) {
    m_strData.clear();
    s=getAttrValue("precursorCharge", attr);
    if(s.length()>0) m_precursorIon.charge=atoi(&s[0]);
    else  m_precursorIon.charge=0;
    s=getAttrValue("precursorIntensity", attr);
    if(s.length()>0) m_precursorIon.intensity=atof(&s[0]);
    else  m_precursorIon.intensity=0.0;
    s=getAttrValue("precursorScanNum", attr);
    if(s.length()>0) spec->setPrecursorScanNum(atoi(&s[0]));
    else spec->setPrecursorScanNum(0);
    m_bInPrecursorMz = true;
    
    s=getAttrValue("activationMethod", attr);
    if(s.length()>0){
      if(!strcmp("CID",&s[0])) spec->setActivation(CID);
      else if(!strcmp("ETD",&s[0])) spec->setActivation(ETD);
      else if(!strcmp("HCD",&s[0])) spec->setActivation(HCD);
      else if(!strcmp("ECD",&s[0])) spec->setActivation(ECD);
      else if(!strcmp("ETD+SA",&s[0])) spec->setActivation(ETDSA);
    } else {
      spec->setActivation(none);
    }

  }  
  else if (isElement("scan", el)) {
    if(m_bInScan){
      pushSpectrum();
      stopParser();
    } else {
      m_bInScan=true;
      spec->setScanNum(atoi(getAttrValue("num", attr)));
      spec->setMSLevel(atoi(getAttrValue("msLevel", attr)));
      spec->setBasePeakIntensity(atof(getAttrValue("basePeakIntensity", attr)));
      spec->setBasePeakMZ(atof(getAttrValue("basePeakMz", attr)));
      spec->setCentroid((bool)getAttrValue("centroided",attr));
      spec->setCollisionEnergy(atof(getAttrValue("collisionEnergy", attr)));
      spec->setCompensationVoltage(atof(getAttrValue("CompensationVoltage", attr)));
      s=getAttrValue("filterLine", attr);
      spec->setFilterLine(&s[0]);
      spec->setHighMZ(atof(getAttrValue("highMz", attr)));
      spec->setLowMZ(atof(getAttrValue("lowMz", attr)));
      if(spec->getHighMZ()==0) spec->setHighMZ(atof(getAttrValue("endMz", attr)));
      if(spec->getLowMZ()==0) spec->setLowMZ(atof(getAttrValue("startMz", attr)));
      spec->setTotalIonCurrent(atof(getAttrValue("totIonCurrent",attr)));
      m_peaksCount = atoi(getAttrValue("peaksCount", attr));
      spec->setPeaksCount(m_peaksCount);
      
      s=getAttrValue("retentionTime", attr);
      if(s.length()>0){
        float f=(float)atof(&s[2]);
        if(s[s.length()-1]=='S') spec->setRTime(f/60.0f);
        else spec->setRTime(f);
      } else {
        spec->setRTime(0.0f);
      }
    }

  }
}


void mzpSAXMzxmlHandler::endElement(const XML_Char *el) {

  if(isElement("dataProcessing", el))  {
    m_bInDataProcessing=false;

  } else if(isElement("index",el)){
    m_bInIndex=false;
    posIndex=-1;
    stopParser();
    if (!m_bIndexSorted) {
      qsort(&m_vIndex[0],m_vIndex.size(),sizeof(cindex),cindex::compare);
      m_bIndexSorted=true;
    }

  } else if(isElement("msInstrument",el)){
    m_vInstrument.push_back(m_instrument);
    m_bInMsInstrument=false;

  } else if(isElement("msRun",el)){
    m_bInMsRun=false;

  } else if(isElement("offset",el) && m_bScanIndex){
    if(!m_bInIndex){ //bad index reference
      cout << "Index offset points to wrong location. Please rebuild index." << endl;
      m_bNoIndex=true;
      stopParser();
    }
    curIndex.offset=mzpatoi64(&m_strData[0]);
    m_vIndex.push_back(curIndex);
    if (m_bIndexSorted && m_vIndex.size() > 1) {
      if (m_vIndex[m_vIndex.size()-1].scanNum<m_vIndex[m_vIndex.size()-2].scanNum) {
        m_bIndexSorted = false;
      }
    }

  } else if(isElement("peaks",el)){
    if(m_bLowPrecision && m_bCompressedData) decompress32();
    else if(m_bLowPrecision && !m_bCompressedData) decode32();
    else if(!m_bLowPrecision && m_bCompressedData) decompress64();
    else decode64();
    m_bInPeaks = false;

  } else if(isElement("precursorMz", el)) {
    m_precursorIon.mz=atof(&m_strData[0]);
    m_precursorIon.isoMZ = m_precursorIon.mz;
    spec->setPrecursorIon(m_precursorIon);
    m_bInPrecursorMz=false;
    
  } else if(isElement("scan",el)) {
    m_bInScan = false;
    pushSpectrum();
    stopParser();

  }
}

void mzpSAXMzxmlHandler::characters(const XML_Char *s, int len) {
  m_strData.append(s, len);
}

bool mzpSAXMzxmlHandler::readHeader(int num){
  spec->clear();

  if(m_bNoIndex){
    cout << "Currently only supporting indexed mzXML" << endl;
    return false;
  }

  //if no scan was requested, grab the next one
  if(num<0){
    posIndex++;
    if(posIndex>=(int)m_vIndex.size()) return false;
    m_bHeaderOnly=true;
    parseOffset(m_vIndex[posIndex].offset);
    m_bHeaderOnly=false;
    return true;
  }

  //Assumes scan numbers are in order
  size_t mid=m_vIndex.size()/2;
  size_t upper=m_vIndex.size();
  size_t lower=0;
  while(m_vIndex[mid].scanNum!=num){
    if(lower==upper) break;
    if(m_vIndex[mid].scanNum>num){
      upper=mid-1;
      mid=(lower+upper)/2;
    } else {
      lower=mid+1;
      mid=(lower+upper)/2;
    }
  }

  //need something faster than this, perhaps binary search
  if(m_vIndex[mid].scanNum==num) {
    m_bHeaderOnly=true;
    parseOffset(m_vIndex[mid].offset);
    //force scan number; this was done for files where scan events are not numbered
    if(spec->getScanNum()!=m_vIndex[mid].scanNum) spec->setScanNum(m_vIndex[mid].scanNum);
    spec->setScanIndex((int)mid+1); //set the index, which starts from 1, so offset by 1
    m_bHeaderOnly=false;
    posIndex=(int)mid;
    return true;
  }
  return false;
}

bool mzpSAXMzxmlHandler::readSpectrum(int num){
  spec->clear();

  if(m_bNoIndex){
    cout << "Currently only supporting indexed mzXML" << endl;
    return false;
  }

  //if no scan was requested, grab the next one
  if(num<0){
    posIndex++;
    if(posIndex>=(int)m_vIndex.size()) return false;
    parseOffset(m_vIndex[posIndex].offset);
    return true;
  }

  //Assumes scan numbers are in order
  size_t mid=m_vIndex.size()/2;
  size_t upper=m_vIndex.size();
  size_t lower=0;
  while(m_vIndex[mid].scanNum!=num){
    if(lower==upper) break;
    if(m_vIndex[mid].scanNum>num){
      upper=mid-1;
      mid=(lower+upper)/2;
    } else {
      lower=mid+1;
      mid=(lower+upper)/2;
    }
  }

  //need something faster than this perhaps
  if(m_vIndex[mid].scanNum==num) {
    parseOffset(m_vIndex[mid].offset);
    //force scan number; this was done for files where scan events are not numbered
    if(spec->getScanNum()!=m_vIndex[mid].scanNum) spec->setScanNum(m_vIndex[mid].scanNum);
    spec->setScanIndex((int)mid+1); //set the index, which starts from 1, so offset by 1
    posIndex=(int)mid;
    return true;
  }
  return false;
}

void mzpSAXMzxmlHandler::pushSpectrum(){

  specDP dp;
  for(unsigned int i=0;i<vdM.size();i++)  {
    dp.mz = vdM[i];
    dp.intensity = vdI[i];
    spec->addDP(dp);
  }
  
}

void mzpSAXMzxmlHandler::decompress32(){

  vdM.clear();
  vdI.clear();
  if(m_peaksCount < 1) return;
  
  union udata {
    float f;
    uint32_t i;
  } uData;

  uLong uncomprLen;
  uint32_t* data;
  int length;
  const char* pData = m_strData.data();
  size_t stringSize = m_strData.size();
  
  //Decode base64
  char* pDecoded = (char*) new char[m_compressLen];
  memset(pDecoded, 0, m_compressLen);
  length = b64_decode_mio( (char*) pDecoded , (char*) pData, stringSize );
  pData=NULL;

  //zLib decompression
  data = new uint32_t[m_peaksCount*2];
  uncomprLen = m_peaksCount * 2 * sizeof(uint32_t);
  uncompress((Bytef*)data, &uncomprLen, (const Bytef*)pDecoded, length);
  delete [] pDecoded;

  //write data to arrays
  int n = 0;
  for(int i=0;i<m_peaksCount;i++){
    uData.i = dtohl(data[n++], m_bNetworkData);
    vdM.push_back((double)uData.f);
    uData.i = dtohl(data[n++], m_bNetworkData);
    vdI.push_back((double)uData.f);
  }
  delete [] data;
}

void mzpSAXMzxmlHandler::decompress64(){

  vdM.clear();
  vdI.clear();
  if(m_peaksCount < 1) return;
  
  union udata {
    double d;
    uint64_t i;
  } uData;

  uLong uncomprLen;
  uint64_t* data;
  int length;
  const char* pData = m_strData.data();
  size_t stringSize = m_strData.size();
  
  //Decode base64
  char* pDecoded = (char*) new char[m_compressLen];
  memset(pDecoded, 0, m_compressLen);
  length = b64_decode_mio( (char*) pDecoded , (char*) pData, stringSize );
  pData=NULL;

  //zLib decompression
  data = new uint64_t[m_peaksCount*2];
  uncomprLen = m_peaksCount * 2 * sizeof(uint64_t);
  uncompress((Bytef*)data, &uncomprLen, (const Bytef*)pDecoded, length);
  delete [] pDecoded;

  //write data to arrays
  int n = 0;
  for(int i=0;i<m_peaksCount;i++){
    uData.i = dtohl(data[n++], m_bNetworkData);
    vdM.push_back(uData.d);
    uData.i = dtohl(data[n++], m_bNetworkData);
    vdI.push_back(uData.d);
  }
  delete [] data;

}


void mzpSAXMzxmlHandler::decode32(){
// This code block was revised so that it packs floats correctly
// on both 64 and 32 bit machines, by making use of the uint32_t
// data type. -S. Wiley
  const char* pData = m_strData.data();
  size_t stringSize = m_strData.size();
  
  size_t size = m_peaksCount * 2 * sizeof(uint32_t);
  char* pDecoded = (char *) new char[size];
  memset(pDecoded, 0, size);

  if(m_peaksCount > 0) {
    // Base64 decoding
    // By comparing the size of the unpacked data and the expected size
    // an additional check of the data file integrity can be performed
    int length = b64_decode_mio( (char*) pDecoded , (char*) pData, stringSize );
    if(length != size) {
      cout << " decoded size " << length << " and required size " << (unsigned long)size << " dont match:\n";
      cout << " Cause: possible corrupted file.\n";
      exit(EXIT_FAILURE);
    }
  }

  // And byte order correction
  union udata {
    float fData;
    uint32_t iData;  
  } uData; 

  vdM.clear();
  vdI.clear();
  int n = 0;
  uint32_t* pDecodedInts = (uint32_t*)pDecoded; // cast to uint_32 for reading int sized chunks
  for(int i = 0; i < m_peaksCount; i++) {
    uData.iData = dtohl(pDecodedInts[n++], m_bNetworkData);
    vdM.push_back((double)uData.fData);
    uData.iData = dtohl(pDecodedInts[n++], m_bNetworkData);
    vdI.push_back((double)uData.fData);
  }

  // Free allocated memory
  delete[] pDecoded;
}

void mzpSAXMzxmlHandler::decode64(){

// This code block was revised so that it packs floats correctly
// on both 64 and 32 bit machines, by making use of the uint32_t
// data type. -S. Wiley
  const char* pData = m_strData.data();
  size_t stringSize = m_strData.size();

  size_t size = m_peaksCount * 2 * sizeof(uint64_t);
  char* pDecoded = (char *) new char[size];
  memset(pDecoded, 0, size);

  if(m_peaksCount > 0) {
    // Base64 decoding
    // By comparing the size of the unpacked data and the expected size
    // an additional check of the data file integrity can be performed
    int length = b64_decode_mio( (char*) pDecoded , (char*) pData, stringSize );
    if(length != size) {
      cout << " decoded size " << length << " and required size " << (unsigned long)size << " dont match:\n";
      cout << " Cause: possible corrupted file.\n";
      exit(EXIT_FAILURE);
    }
  }

  // And byte order correction
  union udata {
    double fData;
    uint64_t iData;  
  } uData; 

  vdM.clear();
  vdI.clear();
  int n = 0;
  uint64_t* pDecodedInts = (uint64_t*)pDecoded; // cast to uint_64 for reading int sized chunks
  for(int i = 0; i < m_peaksCount; i++) {
    uData.iData = dtohl(pDecodedInts[n++], m_bNetworkData);
    vdM.push_back(uData.fData);
    uData.iData = dtohl(pDecodedInts[n++], m_bNetworkData);
    vdI.push_back(uData.fData);
  }

  // Free allocated memory
  delete[] pDecoded;
}

unsigned long mzpSAXMzxmlHandler::dtohl(uint32_t l, bool bNet) {

#ifdef OSX
  if (!bNet)
  {
    l = (l << 24) | ((l << 8) & 0xFF0000) |
      (l >> 24) | ((l >> 8) & 0x00FF00);
  }
#else
  if (bNet)
  {
    l = (l << 24) | ((l << 8) & 0xFF0000) |
      (l >> 24) | ((l >> 8) & 0x00FF00);
  }
#endif
  return l;
}

uint64_t mzpSAXMzxmlHandler::dtohl(uint64_t l, bool bNet) {

#ifdef OSX
  if (!bNet)
  {
    l = (l << 56) | ((l << 40) & 0xFF000000000000LL) | ((l << 24) & 0x0000FF0000000000LL) | ((l << 8) & 0x000000FF00000000LL) |
      (l >> 56) | ((l >> 40) & 0x0000000000FF00LL) | ((l >> 24) & 0x0000000000FF0000LL) | ((l >> 8) & 0x00000000FF000000LL) ;
  }
#else
  if (bNet)
  {
    l = (l << 56) | ((l << 40) & 0x00FF000000000000LL) | ((l << 24) & 0x0000FF0000000000LL) | ((l << 8) & 0x000000FF00000000LL) |
      (l >> 56) | ((l >> 40) & 0x000000000000FF00LL) | ((l >> 24) & 0x0000000000FF0000LL) | ((l >> 8) & 0x00000000FF000000LL) ;
  }
#endif
  return l;
}

//Finding the index list offset is done without the xml parser
//to speed things along. This can be problematic if the <indexListOffset>
//tag is placed anywhere other than the end of the mzML file.
f_off mzpSAXMzxmlHandler::readIndexOffset() {

  char buffer[200];
  char chunk[CHUNK];
  char* start;
  char* stop;
  int readBytes;
  size_t sz;

  if(!m_bGZCompression){
    FILE* f=fopen(&m_strFileName[0],"r");
    mzpfseek(f,-200,SEEK_END);
    sz = fread(buffer,1,200,f);
    fclose(f);
    start=strstr(buffer,"<indexOffset>");
    stop=strstr(buffer,"</indexOffset>");
  } else {
    readBytes = gzObj.extract(fptr, gzObj.getfilesize()-200, (unsigned char*)chunk, CHUNK);
    start=strstr(chunk,"<indexOffset>");
    stop=strstr(chunk,"</indexOffset>");
  }

  if(start==NULL || stop==NULL) {
    cerr << "No index list offset found. File will not be read." << endl;
    return 0;
  }

  char offset[64];
  int len=(int)(stop-start-13);
  strncpy(offset,start+13,len);
  offset[len]='\0';
  return mzpatoi64(offset);

}

bool mzpSAXMzxmlHandler::load(const char* fileName){
  if(!open(fileName)) return false;
  indexOffset = readIndexOffset();
  if(indexOffset==0){
    m_bNoIndex=true;
    return false;
  } else {
    m_bNoIndex=false;
    if(!parseOffset(indexOffset)){
      cerr << "Cannot parse index. Make sure index offset is correct or rebuild index." << endl;
      return false;
    }
    posIndex=-1;
  }
  m_vInstrument.clear();
  parseOffset(0);
  return true;
}


void mzpSAXMzxmlHandler::stopParser(){
  m_bStopParse=true;
  XML_StopParser(m_parser,false);

  //reset mzML flags
  m_bInMsInstrument=false;
  m_bInDataProcessing=false;
  m_bInScan=false;
  m_bInPrecursorMz=false;
  m_bInMsRun=false;
  m_bInIndex=false;
  m_bInPeaks=false;
  m_bLowPrecision=false;

  //reset other flags
  m_bScanIndex=false;
}

int mzpSAXMzxmlHandler::highScan() {
  if(m_vIndex.size()==0) return 0;
  return m_vIndex[m_vIndex.size()-1].scanNum;
}

int mzpSAXMzxmlHandler::lowScan() {
  if(m_vIndex.size()==0) return 0;
  return m_vIndex[0].scanNum;
}

vector<cindex>* mzpSAXMzxmlHandler::getIndex(){
  return &m_vIndex;
}

f_off mzpSAXMzxmlHandler::getIndexOffset(){
  return indexOffset;
}

instrumentInfo mzpSAXMzxmlHandler::getInstrument(){
  return m_instrument;
}

int mzpSAXMzxmlHandler::getPeaksCount(){
  return m_peaksCount;
}
