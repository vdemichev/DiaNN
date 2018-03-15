/************************************************************
 *              SAXMzmlHandler.cpp
 * Adapted from SAXMzdataHandler.cpp
 * August 2008
 * Ronald Beavis
 *
 * April 2009 - Support for referenceable param groups and mzML 1.1.0
 * Fredrik Levander
 *
 * December 2010 - Drastically modified and cannibalized to create
 * robust, yet generic, mzML parser
 * Mike Hoopmann, Institute for Systems Biology
 *
 * Premiere version janvier 2005
 * Patrick Lacasse
 * placasse@mat.ulaval.ca
 *
 * 3/11/2005 (Brendan MacLean): Use eXpat SAX parser, and create SAXSpectraHandler
 *
 * November 2005
 * Fredrik Levander 
 * A few changes to handle MzData 1.05.
 *
 * Updated to handle version 1.04 and 1.05. (Rob Craig)
 *
 *
 * See http://psidev.sourceforge.net/ms/#mzdata for
 * mzData schema information.
 *
 * Inspired by DtaSAX2Handler.cpp
 * copyright            : (C) 2002 by Pedrioli Patrick, ISB, Proteomics
 * email                : ppatrick@systemsbiology.org
 * Artistic License granted 3/11/2005
 *******************************************************/

 // THIS FILE HAS BEEN MODIFIED BY VADIM DEMICHEV

#include "mzParser.h"

mzpSAXMzmlHandler::mzpSAXMzmlHandler(BasicSpectrum* bs){
  m_bChromatogramIndex = false;
  m_bInmzArrayBinary = false;
  m_bInintenArrayBinary = false;
  m_bInRefGroup = false;
  m_bNetworkData = false; //always little-endian for mzML
  m_bNumpressLinear = false;
  m_bNumpressPic = false;
  m_bNumpressSlof = false;
  m_bLowPrecision = false;
  m_bInSpectrumList=false;
  m_bInChromatogramList=false;
  m_bInIndexedMzML=false;
  m_bInIndexList=false;
  m_bInProduct=false;
  m_bHeaderOnly=false;
  m_bSpectrumIndex=false;
  m_bNoIndex=true;
  m_bIndexSorted = true;
  m_bZlib=false;
  m_iDataType=0;
  spec=bs;
  indexOffset=-1;
  m_scanPRECCount = 0;
  m_scanSPECCount = 0;
  m_scanIDXCount = 0;
  chromat=NULL;
}

mzpSAXMzmlHandler::mzpSAXMzmlHandler(BasicSpectrum* bs, BasicChromatogram* cs){
  m_bChromatogramIndex = false;
  m_bInmzArrayBinary = false;
  m_bInintenArrayBinary = false;
  m_bInRefGroup = false;
  m_bNetworkData = false; //always little-endian for mzML
  m_bNumpressLinear = false;
  m_bNumpressPic = false;
  m_bNumpressSlof = false;
  m_bLowPrecision = false;
  m_bInSpectrumList=false;
  m_bInChromatogramList=false;
  m_bInIndexedMzML=false;
  m_bInIndexList=false;
  m_bInProduct=false;
  m_bHeaderOnly=false;
  m_bSpectrumIndex=false;
  m_bNoIndex=true;
  m_bIndexSorted = true;
  m_bZlib=false;
  m_iDataType=0;
  spec=bs;
  chromat=cs;
  indexOffset=-1;
  m_scanPRECCount = 0;
  m_scanSPECCount = 0;
  m_scanIDXCount = 0;
}

mzpSAXMzmlHandler::~mzpSAXMzmlHandler(){
  chromat=NULL;
  spec=NULL;
}

void mzpSAXMzmlHandler::startElement(const XML_Char *el, const XML_Char **attr){

  if (isElement("binaryDataArray",el)){
    m_bNumpressLinear=false;
    string s=getAttrValue("encodedLength", attr);
    m_encodedLen=atoi(&s[0]);

  } else if (isElement("binaryDataArrayList",el)) {
    if(m_bHeaderOnly) stopParser();

  } else if (isElement("chromatogram",el)) {
    string s=getAttrValue("id", attr);
    chromat->setIDString(&s[0]);
    m_peaksCount = atoi(getAttrValue("defaultArrayLength", attr));

  } else if (isElement("chromatogramList",el)) {
    m_bInChromatogramList=true;

  } else if(isElement("index",el) && m_bInIndexList){
    if(!strcmp(getAttrValue("name", attr),"spectrum")) m_bSpectrumIndex=true;
    if(!strcmp(getAttrValue("name", attr),"chromatogram")) m_bChromatogramIndex=true;

  } else if (isElement("indexedmzML",el)) {
    m_vIndex.clear();
    m_bInIndexedMzML=true;

  } else if(isElement("indexList",el)) {
    m_bInIndexList=true;

  } else if(isElement("offset",el) && m_bChromatogramIndex){
    m_strData.clear();
    curChromatIndex.idRef=string(getAttrValue("idRef", attr));

  } else if(isElement("offset",el) && m_bSpectrumIndex){
    m_strData.clear();
    curIndex.idRef=string(getAttrValue("idRef", attr));
    if(strstr(&curIndex.idRef[0],"scan=")!=NULL)  {
      curIndex.scanNum=atoi(strstr(&curIndex.idRef[0],"scan=")+5);
    } else if(strstr(&curIndex.idRef[0],"scanId=")!=NULL) {
      curIndex.scanNum=atoi(strstr(&curIndex.idRef[0],"scanId=")+7);
    } else if(strstr(&curIndex.idRef[0],"S")!=NULL) {
      curIndex.scanNum=atoi(strstr(&curIndex.idRef[0],"S")+1);
    } else {
      curIndex.scanNum=++m_scanIDXCount;
      //Suppressing warning.
      //cout << "WARNING: Cannot extract scan number in index offset line: " << &curIndex.idRef[0] << "\tDefaulting to " << m_scanIDXCount << endl;
    }

  } else if(isElement("precursor",el)) {
    string s=getAttrValue("spectrumRef", attr);

    //if spectrumRef is not provided
    if(s.length()<1){
      spec->setPrecursorScanNum(0);
    } else {
      if(strstr(&s[0],"scan=")!=NULL)  {
        spec->setPrecursorScanNum(atoi(strstr(&s[0],"scan=")+5));
      } else if(strstr(&s[0],"scanId=")!=NULL) {
        spec->setPrecursorScanNum(atoi(strstr(&s[0],"scanId=")+7));
      } else if(strstr(&s[0],"S")!=NULL) {
        spec->setPrecursorScanNum(atoi(strstr(&s[0],"S")+1));
      } else {
        spec->setPrecursorScanNum(++m_scanPRECCount);
        //Suppressing warning.
        //cout << "WARNING: Cannot extract precursor scan number spectrum line: " << &s[0] << "\tDefaulting to " << m_scanPRECCount << endl;
      }
    }

  } else if (isElement("product", el)) {
    m_bInProduct=true;

  }  else if (isElement("referenceableParamGroup", el)) {
    const char* groupName = getAttrValue("id", attr);
    m_ccurrentRefGroupName = string(groupName);
    m_bInRefGroup = true;

  } else if (isElement("run", el)){
    stopParser();

  }  else if (isElement("softwareParam", el)) {
    const char* name = getAttrValue("name", attr);
    const char* accession = getAttrValue("accession", attr);
    const char* version = getAttrValue("version", attr);

  }  else if (isElement("spectrum", el)) {
    string s=getAttrValue("id", attr);
    spec->setIDString(&s[0]);
    if(strstr(&s[0],"scan=")!=NULL)  {
      spec->setScanNum(atoi(strstr(&s[0],"scan=")+5));
    } else if(strstr(&s[0],"scanId=")!=NULL) {
      spec->setScanNum(atoi(strstr(&s[0],"scanId=")+7));
    } else if(strstr(&s[0],"S")!=NULL) {
      spec->setScanNum(atoi(strstr(&s[0],"S")+1));
    } else {
      spec->setScanNum(++m_scanSPECCount);
      //Suppressing warning.
      //cout << "WARNING: Cannot extract scan number spectrum line: " << &s[0] << "\tDefaulting to " << m_scanSPECCount << endl;
    }
    m_peaksCount = atoi(getAttrValue("defaultArrayLength", attr));
    spec->setPeaksCount(m_peaksCount);

  } else if (isElement("spectrumList",el)) {
    m_bInSpectrumList=true;

  }  else if (isElement("cvParam", el)) {
    const char* name = getAttrValue("name", attr);
    const char* accession = getAttrValue("accession", attr);
    const char* value = getAttrValue("value", attr);
    const char* unitName = getAttrValue("unitName", attr);
    const char* unitAccession = getAttrValue("unitAccession", attr);
    if (m_bInRefGroup) {
      cvParam m_cvParam;
      m_cvParam.refGroupName = string(m_ccurrentRefGroupName);
      m_cvParam.name = string(name);
      m_cvParam.accession = string(accession);
      m_cvParam.value = string(value);
      m_cvParam.unitName = string(unitName);
      m_cvParam.unitAccession = string(unitAccession);
      m_refGroupCvParams.push_back(m_cvParam);
    }  else {
      processCVParam(name,accession,value,unitName,unitAccession);
    }

  }  else if (isElement("referenceableParamGroupRef", el))  {
    const char* groupName = getAttrValue("ref", attr);
    for (unsigned int i=0;i<m_refGroupCvParams.size();i++)  {
      if( strcmp(groupName,&m_refGroupCvParams[i].refGroupName[0])==0 )  {
        processCVParam(&m_refGroupCvParams[i].name[0], &m_refGroupCvParams[i].accession[0], &m_refGroupCvParams[i].value[0], &m_refGroupCvParams[i].unitName[0], &m_refGroupCvParams[i].unitAccession[0]);
      }
    }
  } else if (isElement("userParam", el)) {
    const char* name = getAttrValue("name", attr);
    const char* dtype = getAttrValue("type", attr);
    const char* value = getAttrValue("value", attr);
    if(strcmp(name,"[Thermo Trailer Extra]Monoisotopic M/Z:")==0){
      m_precursorIon.monoMZ=atof(value);
    }
  }

  if(isElement("binary", el))  {
    m_strData.clear();
  }

}


void mzpSAXMzmlHandler::endElement(const XML_Char *el) {

  if(isElement("binary", el))  {
    processData();
    m_strData.clear();

  } else if(isElement("binaryDataArray", el))  {
    m_bZlib=false;
    m_bInintenArrayBinary = false;
    m_bInmzArrayBinary = false;
    m_bNumpressLinear=false;
    m_bNumpressSlof=false;
    m_bNumpressPic=false;
    m_iDataType=0;

  } else if(isElement("chromatogram",el)) {
    pushChromatogram();
    stopParser();

  } else if(isElement("chromatogramList",el)) {
    m_bInChromatogramList=false;

  } else if(isElement("componentList",el)) {
    m_vInstrument.push_back(m_instrument);

  } else if(isElement("index",el)){
    m_bSpectrumIndex=false;
    m_bChromatogramIndex=false;
    

  } else if(isElement("indexList",el)){
    m_bInIndexList=false;
    stopParser();
    if (!m_bIndexSorted) {
      qsort(&m_vIndex[0],m_vIndex.size(),sizeof(cindex),cindex::compare);
      m_bIndexSorted=true;
    }

  } else if (isElement("isolationWindow",el)){
    if (chromat != NULL && m_bInProduct) chromat->setProduct(m_precursorIon.isoMZ, m_precursorIon.isoLowerMZ, m_precursorIon.isoUpperMZ);

  } else if(isElement("offset",el) && m_bChromatogramIndex){
    curChromatIndex.offset=mzpatoi64(&m_strData[0]);
    m_vChromatIndex.push_back(curChromatIndex);

  } else if(isElement("offset",el) && m_bSpectrumIndex){
    curIndex.offset=mzpatoi64(&m_strData[0]);
    m_vIndex.push_back(curIndex);
    if (m_bIndexSorted && m_vIndex.size() > 1) {
      if (m_vIndex[m_vIndex.size()-1].scanNum < m_vIndex[m_vIndex.size()-2].scanNum) {
        m_bIndexSorted = false;
      }
    }

  } else if(isElement("precursorList",el)){

  } else if (isElement("product", el)){
    m_bInProduct=false;
    
  } else if(isElement("referenceableParamGroup", el)) {
    m_bInRefGroup = false;

  } else if(isElement("selectedIon",el)) {
    spec->setPrecursorIon(m_precursorIon);
    m_precursorIon.clear();

  }  else if(isElement("spectrum", el)) {
    pushSpectrum();
    stopParser();
    
  } else if(isElement("spectrumList",el)) {
    m_bInSpectrumList = false;

  }
}

void mzpSAXMzmlHandler::characters(const XML_Char *s, int len) {
  m_strData.append(s, len);
}

void mzpSAXMzmlHandler::processCVParam(const char* name, const char* accession, const char* value, const char* unitName, const char* unitAccession)
{
  if(!strcmp(name, "32-bit float") || !strcmp(accession,"MS:1000521"))  {
    m_bLowPrecision = true;
    m_iDataType=1;

  } else if(!strcmp(name, "64-bit float") || !strcmp(accession,"MS:1000523"))  {
    m_bLowPrecision = false;
    m_iDataType=2;

  } else if(!strcmp(name, "base peak intensity") || !strcmp(accession,"MS:1000505"))  {
    spec->setBasePeakIntensity(atof(value));

  } else if(!strcmp(name, "base peak m/z") || !strcmp(accession,"MS:1000504"))  {
    spec->setBasePeakMZ(atof(value));
    
  } else if(!strcmp(name, "centroid spectrum") || !strcmp(accession,"MS:1000127")) {
    spec->setCentroid(true);

  } else if(!strcmp(name, "charge state") || !strcmp(accession,"MS:1000041"))  {
    m_precursorIon.charge = atoi(value);

  } else if(!strcmp(name, "collision-induced dissociation") || !strcmp(accession,"MS:1000133"))  {
    if(spec->getActivation()==ETD) spec->setActivation(ETDSA);
    else spec->setActivation(CID);

  } else if(!strcmp(name, "collision energy") || !strcmp(accession,"MS:1000045"))  {
    spec->setCollisionEnergy(atof(value));

  } else if(!strcmp(name,"electron multiplier") || !strcmp(accession,"MS:1000253")) {
    m_instrument.detector=name;

  } else if(!strcmp(name, "electron transfer dissociation") || !strcmp(accession,"MS:1000598"))  {
    if(spec->getActivation()==CID) spec->setActivation(ETDSA);
    else spec->setActivation(ETD);

  } else if(!strcmp(name, "FAIMS compensation voltage") || !strcmp(accession,"MS:1001581"))  {
    spec->setCompensationVoltage(atof(value));
    
  } else if(!strcmp(name, "filter string") || !strcmp(accession,"MS:1000512"))  {
    char str[128];
    strncpy(str,value,127);
    str[127]='\0';
    spec->setFilterLine(str);

  } else if(!strcmp(name, "highest observed m/z") || !strcmp(accession,"MS:1000527"))  {
    spec->setHighMZ(atof(value));

  } else if(!strcmp(name,"inductive detector") || !strcmp(accession,"MS:1000624")) {
    m_instrument.detector=name;

  } else if(!strcmp(name, "intensity array") || !strcmp(accession,"MS:1000515"))  {
    m_bInintenArrayBinary = true;
    m_bInmzArrayBinary = false;

  } else if (!strcmp(name, "ion injection time") || !strcmp(accession, "MS:1000927"))  {
    spec->setIonInjectionTime(atof(value));

  } else if(!strcmp(name,"isolation window target m/z") || !strcmp(accession,"MS:1000827")) {
    m_precursorIon.isoMZ=atof(value);

  } else if (!strcmp(name, "isolation window lower offset") || !strcmp(accession, "MS:1000828")) {
    m_precursorIon.isoLowerMZ = atof(value);

  } else if (!strcmp(name, "isolation window upper offset") || !strcmp(accession, "MS:1000829")) {
    m_precursorIon.isoUpperMZ = atof(value);

  } else if(!strcmp(name,"LTQ Velos") || !strcmp(accession,"MS:1000855")) {
    m_instrument.model=name;

  } else if(!strcmp(name, "lowest observed m/z") || !strcmp(accession,"MS:1000528"))  {
    spec->setLowMZ(atof(value));

  } else if( !strcmp(name, "MS1 spectrum") || !strcmp(accession,"MS:1000579") ){
    spec->setMSLevel(1);

  } else if( !strcmp(name, "ms level") || !strcmp(accession,"MS:1000511") ){
    spec->setMSLevel(atoi(value));

  } else if( !strcmp(name, "MS-Numpress linear prediction compression") || !strcmp(accession,"MS:1002312") ){
    m_bNumpressLinear = true;

  } else if( !strcmp(name, "MS-Numpress positive integer compression") || !strcmp(accession,"MS:1002313") ){
    m_bNumpressPic = true;

  } else if( !strcmp(name, "MS-Numpress short logged float compression") || !strcmp(accession,"MS:1002314") ){
    m_bNumpressSlof = true;

  } else if(!strcmp(name, "m/z array") || !strcmp(accession,"MS:1000514"))  {
    m_bInmzArrayBinary = true;
    m_bInintenArrayBinary = false;

  } else if(!strcmp(name,"nanoelectrospray") || !strcmp(accession,"MS:1000398")) {
    m_instrument.ionization=name;

  } else if(!strcmp(name,"orbitrap") || !strcmp(accession,"MS:1000484")) {
    m_instrument.analyzer=name;

  } else if(!strcmp(name,"peak intensity") || !strcmp(accession,"MS:1000042")) {
    m_precursorIon.intensity=atof(value);

  } else if(!strcmp(name,"positive scan") || !strcmp(accession,"MS:1000130")) {
    spec->setPositiveScan(true);

  } else if(!strcmp(name,"possible charge state") || !strcmp(accession,"MS:1000633")) {
    m_precursorIon.possibleCharges->push_back(atoi(value));

  } else if(!strcmp(name,"profile spectrum") || !strcmp(accession,"MS:1000128")) {
    spec->setCentroid(false);

  } else if(!strcmp(name,"radial ejection linear ion trap") || !strcmp(accession,"MS:1000083")) {
    m_instrument.analyzer=name;

  } else if(!strcmp(name, "scan start time") || !strcmp(accession,"MS:1000016"))  {
    if(!strcmp(unitName, "minute") || !strcmp(unitAccession,"UO:0000031"))  {
      spec->setRTime((float)atof(value));
    } else {
      spec->setRTime((float)atof(value)/60.0f); //assume if not minutes, then seconds
    }

  } else if(!strcmp(name, "scan window lower limit") || !strcmp(accession,"MS:1000501"))    {
    //TODO: should we also check the units???
    spec->setLowMZ(atof(value));
    
  } else if(!strcmp(name, "scan window upper limit") || !strcmp(accession,"MS:1000500"))    {
    //TODO: should we also check the units???
    spec->setHighMZ(atof(value));
    
  } else if(!strcmp(name, "selected ion m/z") || !strcmp(accession,"MS:1000744"))  {
    m_precursorIon.mz=atof(value); //in Thermo instruments this is the monoisotopic peak (if known) or the selected ion peak.
    if(m_precursorIon.monoMZ!=0) m_precursorIon.monoMZ=atof(value); //if the monoisotopic peak was specified earlier, this is a better value to use.
    else if (m_precursorIon.mz<m_precursorIon.isoMZ) m_precursorIon.monoMZ = atof(value); //failsafe? for when thermo trailer info is missing.

  } else if(!strcmp(name, "time array") || !strcmp(accession,"MS:1000595"))  {
    m_bInmzArrayBinary = true; //note that this uses the m/z designation, although it is a time series
    m_bInintenArrayBinary = false;

  } else if(!strcmp(name, "total ion current") || !strcmp(accession,"MS:1000285"))  {
    spec->setTotalIonCurrent(atof(value));

  } else if(!strcmp(name,"Thermo RAW file") || !strcmp(accession,"MS:1000563")) {
    m_instrument.manufacturer="Thermo Scientific";

  } else if(!strcmp(name, "zlib compression") || !strcmp(accession,"MS:1000574"))  {
    m_bZlib=true;
  }
}

void mzpSAXMzmlHandler::processData()
{
  if(m_bInmzArrayBinary) {
    decode(vdM);
    //if(m_bLowPrecision && !m_bCompressedData) decode32(vdM);
    //else if(m_bLowPrecision && m_bCompressedData) decompress32(vdM);
    //else if(!m_bLowPrecision && !m_bCompressedData) decode64(vdM);
    //else decompress64(vdM);
  } else if(m_bInintenArrayBinary) {
    decode(vdI);
    //if(m_bLowPrecision && !m_bCompressedData) decode32(vdI);
    //else if(m_bLowPrecision && m_bCompressedData) decompress32(vdI);
    //else if(!m_bLowPrecision && !m_bCompressedData) decode64(vdI);
    //else decompress64(vdI);
  }
  //m_bCompressedData=false;
}

bool mzpSAXMzmlHandler::readChromatogram(int num){
  if(chromat==NULL) return false;
  chromat->clear();

  if(m_bNoIndex){
    cout << "Currently only supporting indexed mzML" << endl;
    return false;
  }

  //if no scan was requested, grab the next one
  if(num<0) posChromatIndex++;
  else posChromatIndex=num;
  
  if(posChromatIndex>=(int)m_vChromatIndex.size()) return false;
  parseOffset(m_vChromatIndex[posChromatIndex].offset);
  return true;
}

bool mzpSAXMzmlHandler::readHeader(int num){
  spec->clear();

  if(m_bNoIndex){
    cout << "Currently only supporting indexed mzML" << endl;
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

bool mzpSAXMzmlHandler::readSpectrum(int num){
  spec->clear();

  if(m_bNoIndex){
    cout << "Currently only supporting indexed mzML" << endl;
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
  //for(unsigned int i=0;i<m_vIndex.size();i++){
    if(m_vIndex[mid].scanNum==num) {
      parseOffset(m_vIndex[mid].offset);
      //force scan number; this was done for files where scan events are not numbered
      if(spec->getScanNum()!=m_vIndex[mid].scanNum) spec->setScanNum(m_vIndex[mid].scanNum);
      spec->setScanIndex((int)mid+1); //set the index, which starts from 1, so offset by 1
      posIndex=(int)mid;
      return true;
    }
  //}
  return false;
}

void mzpSAXMzmlHandler::pushChromatogram(){
  TimeIntensityPair tip;
  for(unsigned int i=0;i<vdM.size();i++)  {
    tip.time = vdM[i];
    tip.intensity = vdI[i];
    chromat->addTIP(tip);
  }
}

void mzpSAXMzmlHandler::pushSpectrum(){

  specDP dp;
  for(unsigned int i=0;i<vdM.size();i++)  {
    dp.mz = vdM[i];
    dp.intensity = vdI[i];
    spec->addDP(dp);
  }
  
}

void mzpSAXMzmlHandler::decode(vector<double>& d){

  //If there is no data, back out now
  d.clear();
  if(m_peaksCount < 1) return;

  //For byte order correction
  union udata32 {
    float d;
    uint32_t i;  
  } uData32; 

  union udata64 {
    double d;
    uint64_t i;  
  } uData64; 

  const char* pData = m_strData.data();
  size_t stringSize = m_strData.size();

  char* decoded = new char[m_encodedLen];  //array for decoded base64 string
  int decodeLen;
  Bytef* unzipped = NULL;
  uLong unzippedLen;

  int i;

  //Base64 decoding
  decodeLen = b64_decode_mio(decoded,(char*)pData,stringSize);

  //zlib decompression
  if(m_bZlib) {

    if(m_iDataType==1) {
      unzippedLen = m_peaksCount*sizeof(uint32_t);
    } else if(m_iDataType==2) {
      unzippedLen = m_peaksCount*sizeof(uint64_t);
    } else {
      if(!m_bNumpressLinear && !m_bNumpressSlof && !m_bNumpressPic){
        cout << "Unknown data format to unzip. Stopping file read." << endl;
        exit(EXIT_FAILURE);
      }
    //don't know the unzipped size of numpressed data, so assume it to be no larger than unpressed 64-bit data
    unzippedLen = m_peaksCount*sizeof(uint64_t);
    }

    unzipped = new Bytef[unzippedLen];
    uncompress((Bytef*)unzipped, &unzippedLen, (const Bytef*)decoded, (uLong)decodeLen);
    delete [] decoded;

  }

  //Numpress decompression
  if(m_bNumpressLinear || m_bNumpressSlof || m_bNumpressPic){
    double* unpressed=new double[m_peaksCount];
  
    try{
        if(m_bNumpressLinear){
          if(m_bZlib) ms::numpress::MSNumpress::decodeLinear((unsigned char*)unzipped,(const size_t)unzippedLen,unpressed);
          else ms::numpress::MSNumpress::decodeLinear((unsigned char*)decoded,decodeLen,unpressed);
        } else if(m_bNumpressSlof){
          if(m_bZlib) ms::numpress::MSNumpress::decodeSlof((unsigned char*)unzipped,(const size_t)unzippedLen,unpressed);
          else ms::numpress::MSNumpress::decodeSlof((unsigned char*)decoded,decodeLen,unpressed);
        } else if(m_bNumpressPic){
          if(m_bZlib) ms::numpress::MSNumpress::decodePic((unsigned char*)unzipped,(const size_t)unzippedLen,unpressed);
          else ms::numpress::MSNumpress::decodePic((unsigned char*)decoded,decodeLen,unpressed);
        }
    } catch (const char* ch){
      cout << "Exception: " << ch << endl;
      exit(EXIT_FAILURE);
    }

    if(m_bZlib) delete [] unzipped;
    else delete [] decoded;
    for(i=0;i<m_peaksCount;i++) d.push_back(unpressed[i]);
    delete [] unpressed;
    return;
  }

  //Byte order correction
  if(m_bZlib){
    if(m_iDataType==1){
      uint32_t* unzipped32 = (uint32_t*)unzipped;
      for(i=0;i<m_peaksCount;i++){
        uData32.i = dtohl(unzipped32[i], m_bNetworkData);
        d.push_back(uData32.d);
      }
    } else if(m_iDataType==2) {
      uint64_t* unzipped64 = (uint64_t*)unzipped;
      for(i=0;i<m_peaksCount;i++){
        uData64.i = dtohl(unzipped64[i], m_bNetworkData);
        d.push_back(uData64.d);
      }
    }
    delete [] unzipped;
  } else {
    if(m_iDataType==1){
      uint32_t* decoded32 = (uint32_t*)decoded;
      for(i=0;i<m_peaksCount;i++){
        uData32.i = dtohl(decoded32[i], m_bNetworkData);
        d.push_back(uData32.d);
      }
    } else if(m_iDataType==2) {
      uint64_t* decoded64 = (uint64_t*)decoded;
      for(i=0;i<m_peaksCount;i++){
        uData64.i = dtohl(decoded64[i], m_bNetworkData);
        d.push_back(uData64.d);
      }
    }
    delete [] decoded;
  }

}

unsigned long mzpSAXMzmlHandler::dtohl(uint32_t l, bool bNet) {

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

uint64_t mzpSAXMzmlHandler::dtohl(uint64_t l, bool bNet) {

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
f_off mzpSAXMzmlHandler::readIndexOffset() {

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
    start=strstr(buffer,"<indexListOffset>");
    stop=strstr(buffer,"</indexListOffset>");
  } else {
    readBytes = gzObj.extract(fptr, gzObj.getfilesize()-200, (unsigned char*)chunk, CHUNK);
    chunk[200]='\0';
    start=strstr(chunk,"<indexListOffset>");
    stop=strstr(chunk,"</indexListOffset>");
  }

  if(start==NULL || stop==NULL) {
    cerr << "No index list offset found. File will not be read." << endl;
    return 0;
  }

  char offset[64];
  int len=(int)(stop-start-17);
  strncpy(offset,start+17,len);
  offset[len]='\0';
  return mzpatoi64(offset);

}

bool mzpSAXMzmlHandler::load(const char* fileName){
  if(!open(fileName)) return false;
  m_vInstrument.clear();
  m_vIndex.clear();
  m_vChromatIndex.clear();
  parseOffset(0);
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
    posChromatIndex=-1;
  }
  return true;
}


void mzpSAXMzmlHandler::stopParser(){
  m_bStopParse=true;
  XML_StopParser(m_parser,false);

  //reset mzML flags
  m_bInmzArrayBinary = false;
  m_bInintenArrayBinary = false;
  m_bInRefGroup = false;
  m_bInSpectrumList=false;
  m_bInChromatogramList=false;
  m_bInIndexedMzML=false;
  m_bInIndexList=false;

  //reset other flags
  m_bSpectrumIndex=false;
}

int mzpSAXMzmlHandler::highChromat() {
  return (int)m_vChromatIndex.size();
}

int mzpSAXMzmlHandler::highScan() {
  if(m_vIndex.size()==0) return 0;
  return m_vIndex[m_vIndex.size()-1].scanNum;
}

int mzpSAXMzmlHandler::lowScan() {
  if(m_vIndex.size()==0) return 0;
  return m_vIndex[0].scanNum;
}

vector<cindex>* mzpSAXMzmlHandler::getChromatIndex(){
  return &m_vChromatIndex;
}

f_off mzpSAXMzmlHandler::getIndexOffset(){
  return indexOffset;
}

vector<instrumentInfo>* mzpSAXMzmlHandler::getInstrument(){
  return &m_vInstrument;
}

int mzpSAXMzmlHandler::getPeaksCount(){
  return m_peaksCount;
}

vector<cindex>* mzpSAXMzmlHandler::getSpecIndex(){
  return &m_vIndex;
}
