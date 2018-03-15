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

// THIS FILE HAS BEEN MODIFIED BY VADIM DEMICHEV

#include "mzMLWriter.h"

using namespace MSToolkit;
using namespace ms::numpress::MSNumpress;

MzMLWriter::MzMLWriter(){
  iSpecList=0;
  iChromList=0;
  bFileOpen=false;
  bZlib=false;
  bNumpress=false;
  bTabs=false;
  index=0;
  chromIndex=0;
}

MzMLWriter::~MzMLWriter(){
  if(bFileOpen){
    cout << "WARNING: Destroying MzMLWriter before active file was closed. Expect malformed file and/or loss of data." << endl;
  }
}

bool MzMLWriter::closeList(){
  if (iSpecList==1){
    iSpecList = 2;
  } else if (iChromList == 1){
    iChromList=2;
  } else {
    cout << "Error: closeList() - no lists are open." << endl;
    return false;
  }
  if (bTabs) fprintf(fptr, "   ");
  fprintf(fptr, "</spectrumList>\n");
  return true;
}

bool MzMLWriter::closeMzML(){

  if (iSpecList == 1 || iChromList == 1) {
    cout << "Error: closeMzML() - a spectrumList or chromatogramList is still open." << endl;
    return false;
  }

  if(bTabs) fprintf(fptr,"  ");
  fprintf(fptr,"</run>\n");
  if(bTabs) fprintf(fptr," ");
  fprintf(fptr,"</mzML>\n");
  if(!writeIndex()) return false;
  fprintf(fptr,"</indexedmzML>\n");
  if (iSpecList == 2){
    fseek(fptr,fSpecList,0);
    fprintf(fptr,"%.10d",index);
  }
  if (iChromList == 2){
    fseek(fptr, fChromList, 0);
    fprintf(fptr, "%.10d", chromIndex);
  }

  fclose(fptr);
  fptr=NULL;
  bFileOpen=false;
  return true;
}

bool MzMLWriter::createList(bool specList){
  if (specList){
    if (iSpecList>0){
      cout << "Error: createList() - spectrumList already created." << endl;
      return false;
    } else if (iChromList == 1){
      cout << "Error: createList() - chromatogramList not closed." << endl;
      return false;
    }
    iSpecList=1;
    if (bTabs)fprintf(fptr, "   ");
    fprintf(fptr, "<spectrumList count=\"");
    fSpecList = ftell(fptr);
    fprintf(fptr, "%.10d\">\n", 0);
    return true;
  }

  if (iChromList>0){
    cout << "Error: createList() - chromatogramList already created." << endl;
    return false;
  } else if (iSpecList == 1){
    cout << "Error: createList() - spectrumList not closed." << endl;
    return false;
  }
  iChromList = 1;
  if (bTabs)fprintf(fptr, "   ");
  fprintf(fptr, "<chromatogramList count=\"");
  fChromList = ftell(fptr);
  fprintf(fptr, "%.10d\">\n", 0);
  return true;

}

bool MzMLWriter::createMzML(char* fn){
  if(bFileOpen){
    cout << "Error: createMzML() - cannot create new mzML file. Please finalize existing file first." << endl;
    return false;
  }
  fptr = fopen(fn,"wt");
  if(!fptr){
    cout << "Error: createMzML() - cannot create " << fn << endl;
    return false;
  }
  bFileOpen=true;
  chromIndex = 0;
  index=0;
  vIndex.clear();
  vChromIndex.clear();

  iSpecList = 0;
  iChromList = 0;

  fprintf(fptr,"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
  fprintf(fptr,"<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd\">\n");
  if(bTabs)fprintf(fptr," ");
  fprintf(fptr,"<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\">\n");
  if(bTabs)fprintf(fptr,"  ");
  fprintf(fptr,"<run id=\"0\">\n");

  return true;
}

void MzMLWriter::setNumpress(bool b){
  if(bFileOpen){
    cout << "WARNING: Calling setNumpress before active file was closed. Expect inconsistent data indentation." << endl;
  }
  bNumpress=b;
}

void MzMLWriter::setTabs(bool b){
  if(bFileOpen){
    cout << "WARNING: Calling setTabs before active file was closed. Expect inconsistent data indentation." << endl;
  }
  bTabs=b;
}

void MzMLWriter::setZlib(bool b){
  if(bFileOpen){
    cout << "WARNING: Calling setZlib before active file was closed. Expect inconsistent data indentation." << endl;
  }
  bZlib=b;
}

bool MzMLWriter::writeChromatogram(BasicChromatogram& c){
  if (!bFileOpen){
    cout << "Error: writeChromatogram() - cannot write chromatogram, mzML file not open." << endl;
    return false;
  }
  if (iChromList == 0){
    cout << "Error: writeChromatogram() - cannot write chromatogram, chromatogramList not created." << endl;
    return false;
  }
  if (iChromList == 2){
    cout << "Error: writeChromatogram() - cannot write chromatogram, chromatogramList closed." << endl;
    return false;
  }
  if (!exportChromatogram(c, 4)) return false;
  return true;
}

bool MzMLWriter::writeIndex(){
  if(!bFileOpen){
    cout << "Error: writeIndex() - cannot write index, mzML file not open." << endl;
    return false;
  }

  unsigned int i;
  f_off offset;

  if(bTabs) fprintf(fptr," ");

  offset=ftell(fptr);

  int count=0;
  if (vIndex.size()>0)count++;
  if (vChromIndex.size()>0)count++;

  fprintf(fptr,"<indexList count=\"%d\">\n",count);
  
  if (vIndex.size()>0){
    if(bTabs) {
      fprintf(fptr," ");
      fprintf(fptr," ");
    }
    fprintf(fptr,"<index name=\"spectrum\">\n");
    for(i=0;i<vIndex.size();i++) {
      if(!exportOffset(vIndex[i].id,vIndex[i].offset,3)) return false;
    }
    if(bTabs) {
      fprintf(fptr," ");
      fprintf(fptr," ");
    }
    fprintf(fptr,"</index>\n");
  }

  if (vChromIndex.size()>0){
    if (bTabs) {
      fprintf(fptr, " ");
      fprintf(fptr, " ");
    }
    fprintf(fptr, "<index name=\"chromatogram\">\n");
    for (i = 0; i<vChromIndex.size(); i++) {
      if (!exportOffset(vChromIndex[i].id, vChromIndex[i].offset, 3)) return false;
    }
    if (bTabs) {
      fprintf(fptr, " ");
      fprintf(fptr, " ");
    }
    fprintf(fptr, "</index>\n");
  }

  if(bTabs) fprintf(fptr," ");
  fprintf(fptr,"</indexList>\n");
  if(bTabs) fprintf(fptr," ");
  fprintf(fptr,"<indexListOffset>%lld</indexListOffset>\n",offset);

  return true;
}

bool MzMLWriter::writeSpectra(Spectrum& s){
  if(!bFileOpen){
    cout << "Error: writeSpectra() - cannot write Spectra, mzML file not open." << endl;
    return false;
  }
  if (iSpecList == 0){
    cout << "Error: writeSpectra() - cannot write Spectra, spectrumList not created." << endl;
    return false;
  }
  if (iSpecList == 2){
    cout << "Error: writeSpectra() - cannot write Spectra, spectrumList closed." << endl;
    return false;
  }
  if(!exportSpectrum(s,4)) return false;
  return true;
}

bool MzMLWriter::writeSpectra(MSObject& o){
  for (int i = 0; i < o.size(); i++){
    writeSpectra(o.at(i));
  }
  return true;
}

bool MzMLWriter::exportBinary(char* str, int len, int tabs){
  int i;
  string tbs="";
  if(bTabs) {
    tbs.append(tabs,' ');
    fprintf(fptr,"%s",&tbs[0]);
  }
  fprintf(fptr,"<binary>");
  for(i=0;i<len;i++) fprintf(fptr,"%c",str[i]);
  fprintf(fptr,"</binary>\n");
  return true;
}

bool MzMLWriter::exportBinaryDataArray(BasicChromatogram& c, bool bRT, int tabs){

  int i;
  string tbs = "";
  vector<double>  d;
  vector<float>   f;

  //put data in single array
  if (bRT){
    for (i = 0; i<c.size(); i++) d.push_back(c[i].time);
  } else {
    if (bNumpress) {
      for (i = 0; i<c.size(); i++) d.push_back(c[i].intensity);
    } else {
      for (i = 0; i<c.size(); i++) f.push_back((float)c[i].intensity);
    }
  }

  //numpress if requested
  unsigned char* numpr = NULL;
  double fixedPoint;
  size_t numprLen;
  if (bNumpress){
    numpr = new unsigned char[d.size()*sizeof(double) + 8];
    if (bRT){
      fixedPoint = optimalLinearFixedPoint((double*)&d[0], d.size());
      numprLen = encodeLinear(&d[0], d.size(), numpr, fixedPoint);
    } else {
      fixedPoint = optimalSlofFixedPoint((double*)&d[0], d.size());
      numprLen = encodeSlof(&d[0], d.size(), numpr, fixedPoint);
    }
  }

  //zlib if requested
  unsigned char *zCompr = NULL;
  uLong len, zLen;
  if (bZlib){
    if (bNumpress){
      len = (uLong)numprLen;
      zLen = compressBound(len);
      zCompr = (unsigned char*)calloc((uInt)zLen, 1);
      compress(zCompr, &zLen, (const Bytef*)numpr, len);
    } else {
      if (bRT) len = (uLong)d.size()*sizeof(double);
      else len = (uLong)f.size()*sizeof(float);
      zLen = compressBound(len);
      zCompr = (unsigned char*)calloc((uInt)zLen, 1);
      if (bRT) compress(zCompr, &zLen, (const Bytef*)&d[0], len);
      else compress(zCompr, &zLen, (const Bytef*)&f[0], len);
    }
  }

  //convert to base64
  int sz64;
  char* arr64 = NULL;
  if (bNumpress){
    if (bZlib) sz64 = zLen;
    else sz64 = numprLen;
  } else {
    if (bZlib) sz64 = zLen;
    else {
      if (bRT) sz64 = d.size()*sizeof(double);
      else sz64 = f.size()*sizeof(float);
    }
  }
  i = sz64 % 3;
  if (i>0)sz64 += (3 - i);
  sz64 = sz64 * 4 / 3;

  arr64 = new char[sz64];
  if (bNumpress){
    if (bZlib) i = b64_encode(arr64, (char*)zCompr, zLen);
    else i = b64_encode(arr64, (char*)numpr, numprLen);
  } else {
    if (bZlib) i = b64_encode(arr64, (char*)zCompr, zLen);
    else {
      if (bRT) i = b64_encode(arr64, (char*)&d[0], d.size()*sizeof(double));
      else i = b64_encode(arr64, (char*)&f[0], f.size()*sizeof(float));
    }
  }

  //write to file
  if (bTabs) {
    tbs.append(tabs, ' ');
    fprintf(fptr, "%s", &tbs[0]);
  }
  fprintf(fptr, "<binaryDataArray encodedLength=\"%d\">\n", i);
  if (bZlib) exportCvParam("MS:1000574", "MS", "zlib compression", "", "", "", "", tabs + 1);
  if (bNumpress) {
    if (bRT) exportCvParam("MS:1002312", "MS", "MS-Numpress linear prediction compression", "", "", "", "", tabs + 1);
    else exportCvParam("MS:1002314", "MS", "MS-Numpress short logged float compression", "", "", "", "", tabs + 1);
  } else {
    if (bRT) exportCvParam("MS:1000523", "MS", "64-bit float", "", "", "", "", tabs + 1);
    else exportCvParam("MS:1000521", "MS", "32-bit float", "", "", "", "", tabs + 1);
  }
  if (bRT) exportCvParam("MS:1000595", "MS", "time array", "UO:0000010", "MS", "second", "", tabs + 1);
  else exportCvParam("MS:1000515", "MS", "intensity array", "MS:1000131", "MS", "number of detector counts", "", tabs + 1);
  if (!exportBinary(arr64, i, tabs + 1)) return false;
  if (bTabs) fprintf(fptr, "%s", &tbs[0]);
  fprintf(fptr, "</binaryDataArray>\n");

  //clean up memory
  delete[] arr64;
  if (zCompr != NULL) delete[] zCompr;
  if (numpr != NULL) delete[] numpr;

  return true;
}

bool MzMLWriter::exportBinaryDataArray(Spectrum& s, bool bMZ, int tabs){

  int i;
  string tbs="";
  vector<double>  d;
  vector<float>   f;

  //put data in single array
  if(bMZ){
    for(i=0;i<s.size();i++) d.push_back(s[i].mz);
  } else {
    if(bNumpress) {
      for(i=0;i<s.size();i++) d.push_back(s[i].intensity);
    } else {
      for(i=0;i<s.size();i++) f.push_back(s[i].intensity);
    }
  }

  //numpress if requested
  unsigned char* numpr=NULL;
  double fixedPoint;
  size_t numprLen;
  if(bNumpress){
    numpr = new unsigned char[d.size()*sizeof(double)+8];
    if(bMZ){
      fixedPoint = optimalLinearFixedPoint((double*)&d[0],d.size());
      numprLen = encodeLinear(&d[0],d.size(),numpr,fixedPoint);
    } else {
      fixedPoint = optimalSlofFixedPoint((double*)&d[0],d.size());
      numprLen = encodeSlof(&d[0],d.size(),numpr,fixedPoint);
    }
  }

  //zlib if requested
	unsigned char *zCompr=NULL;
  uLong len, zLen;
  if(bZlib){
    if(bNumpress){
      len = (uLong)numprLen;
      zLen = compressBound(len);
      zCompr = (unsigned char*)calloc((uInt)zLen, 1);
      compress(zCompr, &zLen, (const Bytef*)numpr, len);
    } else {
      if(bMZ) len = (uLong)d.size()*sizeof(double);
      else len = (uLong)f.size()*sizeof(float);
      zLen = compressBound(len);
      zCompr = (unsigned char*)calloc((uInt)zLen, 1);
      if(bMZ) compress(zCompr, &zLen, (const Bytef*)&d[0], len);
      else compress(zCompr, &zLen, (const Bytef*)&f[0], len);
    }
  }

  //convert to base64
  int sz64;
  char* arr64=NULL;
  if(bNumpress){
    if(bZlib) sz64=zLen;
    else sz64=numprLen;
  } else {
    if(bZlib) sz64=zLen;
    else {
      if(bMZ) sz64=d.size()*sizeof(double);
      else sz64=f.size()*sizeof(float);
    }
  }
  i=sz64%3;
  if(i>0)sz64+=(3-i);
  sz64=sz64*4/3;

  arr64=new char[sz64];
  if(bNumpress){
    if(bZlib) i=b64_encode(arr64,(char*)zCompr,zLen);
    else i=b64_encode(arr64,(char*)numpr,numprLen);
  } else {
    if(bZlib) i=b64_encode(arr64,(char*)zCompr,zLen);
    else {
      if(bMZ) i=b64_encode(arr64,(char*)&d[0],d.size()*sizeof(double));
      else i=b64_encode(arr64,(char*)&f[0],f.size()*sizeof(float));
    }
  }

  //write to file
  if(bTabs) {
    tbs.append(tabs,' ');
    fprintf(fptr,"%s",&tbs[0]);
  }
  fprintf(fptr,"<binaryDataArray encodedLength=\"%d\">\n",i);
  if(bZlib) exportCvParam("MS:1000574","MS","zlib compression","","","","",tabs+1);
  if(bNumpress) {
    if(bMZ) exportCvParam("MS:1002312","MS","MS-Numpress linear prediction compression","","","","",tabs+1);
    else exportCvParam("MS:1002314","MS","MS-Numpress short logged float compression","","","","",tabs+1);
  } else {
    if(bMZ) exportCvParam("MS:1000523","MS","64-bit float","","","","",tabs+1);
    else exportCvParam("MS:1000521","MS","32-bit float","","","","",tabs+1);
  }
  if(bMZ) exportCvParam("MS:1000514","MS","m/z array","MS:1000040","MS","m/z","",tabs+1);
  else exportCvParam("MS:1000515","MS","intensity array","MS:1000131","MS","number of detector counts","",tabs+1);
  if(!exportBinary(arr64,i,tabs+1)) return false;
  if(bTabs) fprintf(fptr,"%s",&tbs[0]);
  fprintf(fptr,"</binaryDataArray>\n");

  //clean up memory
  delete [] arr64;
  if(zCompr!=NULL) delete [] zCompr;
  if(numpr!=NULL) delete [] numpr;

  return true;
}

bool MzMLWriter::exportBinaryDataArrayList(BasicChromatogram& c, int tabs){
  string tbs = "";
  if (bTabs) {
    tbs.append(tabs, ' ');
    fprintf(fptr, "%s", &tbs[0]);
  }
  fprintf(fptr, "<binaryDataArrayList count=\"2\">\n");
  if (!exportBinaryDataArray(c, true, tabs + 1)) return false;
  if (!exportBinaryDataArray(c, false, tabs + 1)) return false;
  if (bTabs) fprintf(fptr, "%s", &tbs[0]);
  fprintf(fptr, "</binaryDataArrayList>\n");
  return true;
}

bool MzMLWriter::exportBinaryDataArrayList(Spectrum& s, int tabs){
  string tbs="";
  if(bTabs) {
    tbs.append(tabs,' ');
    fprintf(fptr,"%s",&tbs[0]);
  }
  fprintf(fptr,"<binaryDataArrayList count=\"2\">\n");
  if(!exportBinaryDataArray(s,true,tabs+1)) return false;
  if(!exportBinaryDataArray(s,false,tabs+1)) return false;
  if(bTabs) fprintf(fptr,"%s",&tbs[0]);
  fprintf(fptr,"</binaryDataArrayList>\n");
  return true;
}

bool MzMLWriter::exportActivation(Spectrum& s, int tabs){
  char tmp[128];
  string value;
  string tbs="";
  if(bTabs) {
    tbs.append(tabs,' ');
    fprintf(fptr,"%s",&tbs[0]);
  }
  
  fprintf(fptr,"<activation>\n");
  switch(s.getActivationMethod()){
    case mstCID:
      exportCvParam("MS:1000133","MS","collision-induced dissociation","","","","",tabs+1);
      //TODO: get activation energy
      //exportCvParam("MS:1000045","MS","collision energy","","","","",tabs+1);
      break;
    case mstETD:
    case mstHCD:
      break;
    default:
      break;
  }
  if(bTabs) fprintf(fptr,"%s",&tbs[0]);
  fprintf(fptr,"</activation>\n");

  return true;
}

bool MzMLWriter::exportChromatogram(BasicChromatogram& c, int tabs){
  char tmp[128];
  string value;
  string tbs = "";
  sMzMLIndex x;

  if (bTabs) {
    tbs.append(tabs, ' ');
    fprintf(fptr, "%s", &tbs[0]);
  }

  c.getIDString(tmp);
  x.id = tmp;
  x.offset = ftell(fptr);
  vChromIndex.push_back(x);

  fprintf(fptr, "<chromatogram index=\"%d\" id=\"%s\" defaultArrayLength=\"%d\">\n", chromIndex++, tmp, c.size());
  if (!exportPrecursor(c, tabs + 1)) return false;
  if (c.getProdMZ()>0) {
    if (!exportProduct(c,tabs+1)) return false;
  }
  if (!exportBinaryDataArrayList(c, tabs + 1)) return false;

  if (bTabs) fprintf(fptr, "%s", &tbs[0]);
  fprintf(fptr, "</chromatogram>\n");
  return true;
}

bool MzMLWriter::exportCvParam(string ac, string ref, string name, string unitAc, string unitRef, string unitName, string value, int tabs){

  string ex="";
  if(bTabs) ex.append(tabs,' ');
  ex+="<cvParam cvRef=\"" + ref + "\" accession=\"" + ac + "\" name=\"" + name + "\"";
  if(value.size()>0) ex+=" value=\"" + value + "\"";
  if(unitRef.size()>0) ex+=" unitCvRef=\"" + unitRef + "\"";
  if(unitAc.size()>0) ex+=" unitAccession=\"" + unitAc + "\"";
  if(unitName.size()>0) ex+=" unitName=\"" + unitName + "\"";
  ex+="/>";
  fprintf(fptr,"%s\n",&ex[0]);
  return true;
}

bool MzMLWriter::exportIsolationWindow(BasicChromatogram& c, bool bPre, int tabs){
  char tmp[128];
  string value;
  string tbs = "";
  if (bTabs) {
    tbs.append(tabs, ' ');
    fprintf(fptr, "%s", &tbs[0]);
  }

  fprintf(fptr, "<isolationWindow>\n");
  if (bPre) sprintf(tmp, "%.4lf", c.getPreMZ());
  else sprintf(tmp, "%.4lf", c.getProdMZ());
  value = tmp;
  exportCvParam("MS:1000827", "MS", "isolation window target m/z", "MS:1000040", "MS", "m/z", value, tabs + 1);
  if (bPre) sprintf(tmp, "%.4lf", c.getPreOffsetLower());
  else sprintf(tmp, "%.4lf", c.getProdOffsetLower());
  value = tmp;
  exportCvParam("MS:1000828", "MS", "isolation window lower offset", "MS:1000040", "MS", "m/z", value, tabs + 1);
  if (bPre) sprintf(tmp, "%.4lf", c.getPreOffsetUpper());
  else sprintf(tmp, "%.4lf", c.getProdOffsetUpper());
  value = tmp;
  exportCvParam("MS:1000829", "MS", "isolation window upper offset", "MS:1000040", "MS", "m/z", value, tabs + 1);
  if (bTabs) fprintf(fptr, "%s", &tbs[0]);
  fprintf(fptr, "</isolationWindow>\n");

  return true;
}

bool MzMLWriter::exportIsolationWindow(Spectrum& s, int tabs){
  char tmp[128];
  string value;
  string tbs="";
  if(bTabs) {
    tbs.append(tabs,' ');
    fprintf(fptr,"%s",&tbs[0]);
  }
  
  fprintf(fptr,"<isolationWindow>\n");
  sprintf(tmp,"%.2lf",s.getMZ());
  value=tmp;
  exportCvParam("MS:1000827","MS","isolation window target m/z","MS:1000040","MS","m/z",value,tabs+1);
  if(bTabs) fprintf(fptr,"%s",&tbs[0]);
  fprintf(fptr,"</isolationWindow>\n");

  return true;
}

bool MzMLWriter::exportOffset(string idRef, f_off offset, int tabs){

  string tbs="";
  if(bTabs) {
    tbs.append(tabs,' ');
    fprintf(fptr,"%s",&tbs[0]);
  }
  fprintf(fptr,"<offset idRef=\"%s\">%lld</offset>\n",&idRef[0],offset);
  return true;
}

bool MzMLWriter::exportPrecursor(BasicChromatogram& c, int tabs){
  string tbs = "";
  if (bTabs) {
    tbs.append(tabs, ' ');
    fprintf(fptr, "%s", &tbs[0]);
  }

  fprintf(fptr, "<precursor>\n");
  if (!exportIsolationWindow(c, true, tabs + 1)) return false;
  if (c.getPreMZ()>0) {
    if (!exportSelectedIonList(c, tabs + 1)) return false;
  }
  //if (!exportActivation(s, tabs + 1)) return false;

  if (bTabs) fprintf(fptr, "%s", &tbs[0]);
  fprintf(fptr, "</precursor>\n");

  return true;
}

bool MzMLWriter::exportPrecursor(Spectrum& s, int tabs){
  string tbs="";
  if(bTabs) {
    tbs.append(tabs,' ');
    fprintf(fptr,"%s",&tbs[0]);
  }
  
  fprintf(fptr,"<precursor>\n");
  if(!exportIsolationWindow(s,tabs+1)) return false;
  //if(s.getMonoMZ()>0) {
    if(!exportSelectedIonList(s,tabs+1)) return false;
  //}
  if(!exportActivation(s,tabs+1)) return false;

  if(bTabs) fprintf(fptr,"%s",&tbs[0]);
  fprintf(fptr,"</precursor>\n");

  return true;
}

bool MzMLWriter::exportPrecursorList(Spectrum& s, int tabs){
  string tbs="";
  if(bTabs) {
    tbs.append(tabs,' ');
    fprintf(fptr,"%s",&tbs[0]);
  }
  fprintf(fptr,"<precursorList count=\"1\">\n");
  if(!exportPrecursor(s,tabs+1)) return false;
  if(bTabs) fprintf(fptr,"%s",&tbs[0]);
  fprintf(fptr,"</precursorList>\n");
  return true;
}

bool MzMLWriter::exportProduct(BasicChromatogram& c, int tabs){
  if (bTabs) exportTabs(tabs);
  fprintf(fptr, "<product>\n");
  if (!exportIsolationWindow(c, false, tabs + 1)) return false;
  if (bTabs) exportTabs(tabs);
  fprintf(fptr, "</product>\n");

  return true;
}

bool MzMLWriter::exportScan(Spectrum& s, int tabs){
  char tmp[128];
  string value;
  
  if (bTabs) exportTabs(tabs);
  fprintf(fptr,"<scan>\n");

  sprintf(tmp,"%f",s.getRTime());
  value=tmp;
  exportCvParam("MS:1000016","MS","scan start time","UO:0000031","UO","minute",value,tabs+1);
  s.getRawFilter(tmp,128);
  if (strlen(tmp) > 1){
    value=tmp;
    exportCvParam("MS:1000512", "MS", "filter string", "", "", "",value,tabs+1);
  }

  if (!exportScanWindowList(s, tabs + 1)) return false;

  if (bTabs) exportTabs(tabs);
  fprintf(fptr,"</scan>\n");

  return true;
}

bool MzMLWriter::exportScanList(Spectrum& s, int tabs){
  if (bTabs) exportTabs(tabs);
  fprintf(fptr,"<scanList count=\"1\">\n");
  if(!exportScan(s,tabs+1)) return false;
  if (bTabs) exportTabs(tabs);
  fprintf(fptr,"</scanList>\n");
  return true;
}

bool MzMLWriter::exportScanWindow(Spectrum& s, int tabs){
  char tmp[128];
  string value;
 
  if (bTabs) exportTabs(tabs);
  fprintf(fptr, "<scanWindow>\n");

  sprintf(tmp, "%.8lf", s.getScanWindowLower());
  value = tmp;
  exportCvParam("MS:1000501", "MS", "scan window lower limit", "MS:1000040", "MS", "m/z", value, tabs + 1);

  sprintf(tmp, "%.8lf", s.getScanWindowUpper());
  value = tmp;
  exportCvParam("MS:1000500", "MS", "scan window upper limit", "MS:1000040", "MS", "m/z", value, tabs + 1);

  if (bTabs) exportTabs(tabs);
  fprintf(fptr, "</scanWindow>\n");
  return true;
}

bool MzMLWriter::exportScanWindowList(Spectrum& s, int tabs){
  if (bTabs) exportTabs(tabs);
  fprintf(fptr, "<scanWindowList count=\"1\">\n");

  if (!exportScanWindow(s, tabs + 1)) return false;

  if (bTabs) exportTabs(tabs);
  fprintf(fptr, "</scanWindowList>\n");
  return true;
}

bool MzMLWriter::exportSelectedIon(BasicChromatogram& c, int tabs){
  char tmp[128];
  string value;
  
  if (bTabs) exportTabs(tabs);
  fprintf(fptr, "<selectedIon>\n");

  sprintf(tmp, "%.8lf", c.getPreMZ());
  value = tmp;
  exportCvParam("MS:1000744", "MS", "selected ion m/z", "MS:1000040", "MS", "m/z", value, tabs + 1);
  sprintf(tmp, "%d", c.getCharge());
  value = tmp;
  exportCvParam("MS:1000041", "MS", "charge state", "", "", "", value, tabs + 1);
  exportCvParam("MS:1000042", "MS", "peak intensity", "MS:1000132", "MS", "percent of base peak",value,tabs+1);

  if (bTabs) exportTabs(tabs);
  fprintf(fptr, "</selectedIon>\n");

  return true;
}

bool MzMLWriter::exportSelectedIon(Spectrum& s, int tabs){
  int i;
  char tmp[128];
  string value;
  
  if (bTabs) exportTabs(tabs);
  fprintf(fptr,"<selectedIon>\n");
  if(s.getMonoMZ()>0) sprintf(tmp,"%.12lf",s.getMonoMZ());
  else sprintf(tmp, "%.2lf", s.getMZ());
  value=tmp;
  exportCvParam("MS:1000744","MS","selected ion m/z","MS:1000040","MS","m/z",value,tabs+1);
  if(s.sizeZ()==1) {
    sprintf(tmp,"%d",s.atZ(0).z);
    value=tmp;
    exportCvParam("MS:1000041","MS","charge state","","","",value,tabs+1);
  } else {
    for(i=0;i<s.sizeZ();i++){
      sprintf(tmp,"%d",s.atZ(i).z);
      value=tmp;
      exportCvParam("MS:1000633","MS","possible charge state","","","",value,tabs+1);
    }
  }
  if (bTabs) exportTabs(tabs);
  fprintf(fptr,"</selectedIon>\n");

  return true;
}

bool MzMLWriter::exportSelectedIonList(BasicChromatogram& c, int tabs){
  if (bTabs) exportTabs(tabs);
  fprintf(fptr, "<selectedIonList count=\"1\">\n");
  if (!exportSelectedIon(c, tabs + 1)) return false;
  if (bTabs) exportTabs(tabs);
  fprintf(fptr, "</selectedIonList>\n");
  return true;
}

bool MzMLWriter::exportSelectedIonList(Spectrum& s, int tabs){
  if (bTabs) exportTabs(tabs);
  fprintf(fptr,"<selectedIonList count=\"1\">\n");
  if(!exportSelectedIon(s,tabs+1)) return false;
  if (bTabs) exportTabs(tabs);
  fprintf(fptr,"</selectedIonList>\n");
  return true;
}

bool MzMLWriter::exportSpectrum(Spectrum& s, int tabs){
  char tmp[128];
  string value;
  sMzMLIndex x;

  if (bTabs) exportTabs(tabs);

  sprintf(tmp,"scan=%d",s.getScanNumber());
  x.id=tmp;
  x.offset=ftell(fptr);
  vIndex.push_back(x);

  fprintf(fptr,"<spectrum index=\"%d\" id=\"scan=%d\" defaultArrayLength=\"%d\">\n",index++,s.getScanNumber(),s.size());
  sprintf(tmp,"%d",s.getMsLevel());
  value=tmp;
  exportCvParam("MS:1000511","MS","ms level","","","",value,tabs+1);
  sprintf(tmp,"%.12lf",s[0].mz);
  value=tmp;
  exportCvParam("MS:1000528","MS","lowest observed m/z","MS:1000040","MS","m/z",value,tabs+1);
  sprintf(tmp,"%.12lf",s[s.size()-1].mz);
  value=tmp;
  exportCvParam("MS:1000527","MS","highest observed m/z","MS:1000040","MS","m/z",value,tabs+1);
  if(!exportScanList(s,tabs+1)) return false;
  if(s.getMsLevel()>1){
    if(!exportPrecursorList(s,tabs+1)) return false;
  }
  if(!exportBinaryDataArrayList(s,tabs+1)) return false;

  if (bTabs) exportTabs(tabs);
  fprintf(fptr,"</spectrum>\n");
  return true;
}

void MzMLWriter::exportTabs(int tabs){
  string tbs = "";
  tbs.append(tabs, ' ');
  fprintf(fptr, "%s", &tbs[0]);
}
