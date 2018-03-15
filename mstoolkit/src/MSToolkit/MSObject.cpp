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
#include "MSObject.h"
#include <iostream>

using namespace std;
using namespace MSToolkit;

MSObject::MSObject(){
  vSpectrum = new vector<Spectrum>;
  fileName="";
  fileType=Unspecified;
	for(int i=0;i<16;i++)	strcpy(header.header[i],"\0");
}

MSObject::~MSObject(){
  delete vSpectrum;
}

MSObject::MSObject(const MSObject& m){
  unsigned int i;
  vSpectrum = new vector<Spectrum>;

  for(i=0;i<m.vSpectrum->size();i++){
    vSpectrum->push_back(m.vSpectrum->at(i));
  }

  fileType = m.fileType;
  fileName = m.fileName;
	for(i=0;i<16;i++)	strcpy(header.header[i],m.header.header[i]);

}

MSObject& MSObject::operator=(const MSObject& m){
  unsigned int i;
  if (this!=&m){
    delete vSpectrum;
    vSpectrum = new vector<Spectrum>;
    for(i=0;i<m.vSpectrum->size();i++){
      vSpectrum->push_back(m.vSpectrum->at(i));
    }
    fileType = m.fileType;
    fileName = m.fileName;
		for(i=0;i<16;i++)	strcpy(header.header[i],m.header.header[i]);
  }
  return *this;
}

  
void MSObject::add(Spectrum& s){
  vSpectrum->push_back(s);
};

bool MSObject::addToHeader(char* c){
	if(strlen(c)>127) return false;
	for(int i=0;i<16;i++){
		if(header.header[i][0]=='\0'){
			strcpy(header.header[i],c);
			return true;
		};
	};
	return false;
};

bool MSObject::addToHeader(string s){
	if(s.size()>127) return false;
	for(int i=0;i<16;i++){
		if(header.header[i][0]=='\0'){
			strcpy(header.header[i],&s[0]);
			return true;
		};
	};
	return false;
};

Spectrum& MSObject::at(unsigned int i){
  return vSpectrum->at(i);
};

Peak_T& MSObject::at(unsigned int i, unsigned int j){
  return vSpectrum->at(i).at(j);
};

void MSObject::clear(){
  delete vSpectrum;
  vSpectrum = new vector<Spectrum>;
	for(int i=0;i<16;i++) strcpy(header.header[i],"\0");
};

MSHeader& MSObject::getHeader(){
	return header;
};

void MSObject::setHeader(const MSHeader& h){
	for(int i=0;i<16;i++)	strcpy(header.header[i],h.header[i]);
};

int MSObject::size(){
  return (int)vSpectrum->size();
};

