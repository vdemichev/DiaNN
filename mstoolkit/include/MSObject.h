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
#ifndef _MSOBJECT_H
#define _MSOBJECT_H

#include "MSToolkitTypes.h"
#include "Spectrum.h"

#include <string>

using namespace std;

namespace MSToolkit {
class MSObject {
 public:
  //Constructors & Destructors
  MSObject();
  MSObject(const MSObject&);
  ~MSObject(); 

  //Operator Functions
  MSObject& operator=(const MSObject&);
  
  //Functions
  void add(Spectrum&);
  bool addToHeader(char*);
  bool addToHeader(string);
  Spectrum& at(unsigned int);
  Peak_T& at(unsigned int, unsigned int);
  void clear();
  void erase(unsigned int);
  void erase(unsigned int, unsigned int);
  MSHeader& getHeader();
  void setHeader(const MSHeader& h);
  int size();
  
 protected:
 private:
  vector<Spectrum> *vSpectrum;
  string fileName;
	MSHeader header;
  MSSpectrumType fileType;
  
};

}
#endif

