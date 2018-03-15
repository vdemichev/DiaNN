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
#ifndef _PEPXMLWRITER_H
#define _PEPXMLWRITER_H

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>
#include <vector>

using namespace std;

//Some simple data structures for PepXMLWriter (pxw)
typedef struct pxwBasicXMLTag {
  string name;
  string value;
} pxwBasicXMLTag;

typedef struct pxwModAA{
  int position;
  double mass;
} pxwModAA;

typedef struct pxwMSMSRunSummary {
  string base_name;
  string raw_data_type;
  string raw_data;
  string search_engine;
} pxwMSMSRunSummary;

typedef struct pxwProtein {
  string protein;
  char peptide_next_aa;
  char peptide_prev_aa;
  int protein_link_pos_a;
  int protein_link_pos_b;
} pxwProtein;

//Use classes for more complicated structures with dynamic arrays
//and additional fuctions.
class PXWModInfo{
public:
  double mod_cterm_mass;
  double mod_nterm_mass;
  string modified_peptide;
  
  PXWModInfo(){
    mod_cterm_mass=0;
    mod_nterm_mass=0;
    modified_peptide.clear();
    mods=new vector<pxwModAA>;
  }
  PXWModInfo(const PXWModInfo& s){
    mod_cterm_mass=s.mod_cterm_mass;
    mod_nterm_mass=s.mod_nterm_mass;
    modified_peptide = s.modified_peptide;
    mods=new vector<pxwModAA>;
    for(size_t i=0;i<s.mods->size();i++) mods->push_back(s.mods->at(i));
  }
  ~PXWModInfo(){
    delete mods;
  }
  PXWModInfo& operator=(const PXWModInfo& s){
    if(this!=&s){
      mod_cterm_mass = s.mod_cterm_mass;
      mod_nterm_mass = s.mod_nterm_mass;
      modified_peptide = s.modified_peptide;
      delete mods;
      mods=new vector<pxwModAA>;
      for(size_t i=0;i<s.mods->size();i++) mods->push_back(s.mods->at(i));
    }
    return *this;
  }

  void addMod(pxwModAA& p){
    mods->push_back(p);
  }
  void addMod(int pos, double mass){
    pxwModAA p;
    p.position=pos;
    p.mass=mass;
    addMod(p);
  }
  void clear(){
    mod_cterm_mass=0;
    mod_nterm_mass=0;
    modified_peptide.clear();
    mods->clear();
  }
  pxwModAA& getMod(int index){
    return mods->at(index);
  }
  pxwModAA& getMod(size_t index){
	  return mods->at(index);
  }
  size_t sizeMods(){
    return mods->size();
  }

private:
  vector<pxwModAA>* mods;
};

class PXWSearchSummary {
public:
  string base_name;
  string search_database;
  string search_engine;
  string search_engine_version;
  int precursor_mass_type; //0=monoisotopic, 1=average
  int fragment_mass_type; //0=monoisotopic, 1=average
  vector<pxwBasicXMLTag>* parameters;

  PXWSearchSummary(){
    base_name.clear();
    search_database.clear();
    search_engine.clear();
    search_engine_version.clear();
    precursor_mass_type=0;
    fragment_mass_type=0;
    parameters=new vector<pxwBasicXMLTag>;
  }
  PXWSearchSummary(const PXWSearchSummary& s){
    base_name=s.base_name;
    search_database=s.search_database;
    search_engine=s.search_engine;
    search_engine_version=s.search_engine_version;
    precursor_mass_type=s.precursor_mass_type;
    fragment_mass_type=s.fragment_mass_type;
    parameters=new vector<pxwBasicXMLTag>;
    for(size_t i=0;i<s.parameters->size();i++) parameters->push_back(s.parameters->at(i));
  }
  ~PXWSearchSummary(){
    delete parameters;
  }
  PXWSearchSummary& operator=(const PXWSearchSummary& s){
    if(this!=&s){
      base_name=s.base_name;
      search_database = s.search_database;
      search_engine=s.search_engine;
      search_engine_version=s.search_engine_version;
      precursor_mass_type=s.precursor_mass_type;
      fragment_mass_type=s.fragment_mass_type;
      delete parameters;
      parameters=new vector<pxwBasicXMLTag>;
      for(size_t i=0;i<s.parameters->size();i++) parameters->push_back(s.parameters->at(i));
    }
    return *this;
  }
};

class PXWSearchHit {
public:
  int hit_rank;
  string peptide;
  int num_tot_proteins;
  double calc_neutral_pep_mass;
  double calc_neutral_xl_mass;
  double massdiff;
  double xl_massdiff;
  PXWModInfo modInfo;
  string xlink_type; //na,loop,xl

  PXWSearchHit(){
    hit_rank=0;
    peptide.clear();
    num_tot_proteins=0;
    calc_neutral_pep_mass=0;
    calc_neutral_xl_mass=0;
    massdiff=0;
    xl_massdiff=0;
    modInfo.clear();
    xlink_type="na";
    proteins=new vector<pxwProtein>;
    searchScores=new vector<pxwBasicXMLTag>;
    xlScores=new vector<pxwBasicXMLTag>;
  }
  PXWSearchHit(const PXWSearchHit& s){
    size_t i;
    hit_rank=s.hit_rank;
    peptide=s.peptide;
    num_tot_proteins=s.num_tot_proteins;
    calc_neutral_pep_mass=s.calc_neutral_pep_mass;
    calc_neutral_xl_mass=s.calc_neutral_xl_mass;
    massdiff=s.massdiff;
    xl_massdiff=s.xl_massdiff;
    modInfo=s.modInfo;
    xlink_type=s.xlink_type;
    proteins=new vector<pxwProtein>;
    searchScores=new vector<pxwBasicXMLTag>;
    xlScores=new vector<pxwBasicXMLTag>;
    for(i=0;i<s.proteins->size();i++) proteins->push_back(s.proteins->at(i));
    for(i=0;i<s.searchScores->size();i++) searchScores->push_back(s.searchScores->at(i));
    for(i=0;i<s.xlScores->size();i++) xlScores->push_back(s.xlScores->at(i));
  }
  ~PXWSearchHit(){
    delete proteins;
    delete searchScores;
    delete xlScores;
  }
  PXWSearchHit& operator=(const PXWSearchHit& s){
    if(this!=&s){
      size_t i;
      hit_rank=s.hit_rank;
      peptide=s.peptide;
      num_tot_proteins=s.num_tot_proteins;
      calc_neutral_pep_mass=s.calc_neutral_pep_mass;
      calc_neutral_xl_mass=s.calc_neutral_xl_mass;
      massdiff=s.massdiff;
      xl_massdiff=s.xl_massdiff;
      modInfo=s.modInfo;
      xlink_type=s.xlink_type;
      delete proteins;
      delete searchScores;
      delete xlScores;
      proteins=new vector<pxwProtein>;
      searchScores=new vector<pxwBasicXMLTag>;
      xlScores=new vector<pxwBasicXMLTag>;
      for(i=0;i<s.proteins->size();i++) proteins->push_back(s.proteins->at(i));
      for(i=0;i<s.searchScores->size();i++) searchScores->push_back(s.searchScores->at(i));
      for(i=0;i<s.xlScores->size();i++) xlScores->push_back(s.xlScores->at(i));
    }
    return *this;
  }

  void addProtein(pxwProtein& p){
    proteins->push_back(p);
  }
  void addProtein(char* protein, char peptide_next_aa, char peptide_prev_aa, int protein_link_pos_a = 0, int protein_link_pos_b = 0){
    pxwProtein p;
    p.protein=protein;
    p.peptide_next_aa=peptide_next_aa;
    p.peptide_prev_aa=peptide_prev_aa;
    p.protein_link_pos_a = protein_link_pos_a;
    p.protein_link_pos_b = protein_link_pos_b;
    addProtein(p);
  }
  void addProtein(string& protein, char peptide_next_aa, char peptide_prev_aa, int protein_link_pos_a = 0, int protein_link_pos_b = 0){
    pxwProtein p;
    p.protein=protein;
    p.peptide_next_aa=peptide_next_aa;
    p.peptide_prev_aa=peptide_prev_aa;
    p.protein_link_pos_a=protein_link_pos_a;
    p.protein_link_pos_b=protein_link_pos_b;
    addProtein(p);
  }
  void addScore(pxwBasicXMLTag& s){
    searchScores->push_back(s);
  }
  void addScore(const char* name, const char* value){
    pxwBasicXMLTag x;
    x.name=name;
    x.value=value;
    addScore(x);
  }
  void addScore(string& name, string& value){
    pxwBasicXMLTag x;
    x.name=name;
    x.value=value;
    addScore(x);
  }
  void addXLScore(pxwBasicXMLTag& s){
    xlScores->push_back(s);
  }
  void addXLScore(const char* name, const char* value){
    pxwBasicXMLTag x;
    x.name=name;
    x.value=value;
    addXLScore(x);
  }
  void addXLScore(string& name, string& value){
    pxwBasicXMLTag x;
    x.name=name;
    x.value=value;
    addXLScore(x);
  }
  void clear(){
    hit_rank=0;
    peptide.clear();
    num_tot_proteins=0;
    calc_neutral_pep_mass=0;
    calc_neutral_xl_mass=0;
    massdiff=0;
    xl_massdiff=0;
    xlink_type="na";
    modInfo.clear();
    proteins->clear();
    searchScores->clear();
    xlScores->clear();
  }
  pxwProtein& getProtein(int index){
    return proteins->at(index);
  }
  pxwProtein& getProtein(size_t index){
    return proteins->at(index);
  }
  pxwBasicXMLTag& getScore(int index){
    return searchScores->at(index);
  }
  pxwBasicXMLTag& getScore(size_t index){
    return searchScores->at(index);
  }
  pxwBasicXMLTag& getXLScore(int index){
	  return xlScores->at(index);
  }
  pxwBasicXMLTag& getXLScore(size_t index){
    return xlScores->at(index);
  }
  size_t sizeProteins(){
    return proteins->size();
  }
  size_t sizeScores(){
    return searchScores->size();
  }
  size_t sizeXLScores(){
    return xlScores->size();
  }

private:
  vector<pxwProtein>* proteins;
  vector<pxwBasicXMLTag>* searchScores;
  vector<pxwBasicXMLTag>* xlScores;
};

typedef struct pxwSearchHitPair{
  PXWSearchHit* a;
  PXWSearchHit* b;
  string identifier;
  double mass;
  pxwSearchHitPair(){
    a=NULL;
    b=NULL;
    identifier.clear();
    mass=0;
  }
  pxwSearchHitPair(const pxwSearchHitPair& p){
    a=NULL;
    b=NULL;
    identifier=p.identifier;
    mass=p.mass;
    if(p.a!=NULL) {
      a=new PXWSearchHit();
      *a=*p.a;
    }
    if(p.b!=NULL){
      b=new PXWSearchHit();
      *b=*p.b;
    }
  }
  ~pxwSearchHitPair(){
    if(a!=NULL) delete a;
    if(b!=NULL) delete b;
  }
  pxwSearchHitPair& operator=(const pxwSearchHitPair& p){
    if(this!=&p){
      if(a!=NULL) delete a;
      if(b!=NULL) delete b;
      a=NULL;
      b=NULL;
      identifier=p.identifier;
      mass=p.mass;
      if(p.a!=NULL) {
        a=new PXWSearchHit();
        *a=*p.a;
      }
      if(p.b!=NULL){
        b=new PXWSearchHit();
        *b=*p.b;
      }
    }
    return *this;
  }
} pxwSearchHitPair;

typedef struct pxwSampleEnzyme{
  string name;
  string cut;
  string no_cut;
  string sense;
  int maxNumInternalCleavages;
  int minNumTermini;
} pxwSampleEnzyme;

class PXWSpectrumQuery {
public:
  string spectrum;
  int start_scan;
  int end_scan;
  double retention_time_sec;
  double precursor_neutral_mass;
  int assumed_charge;
  
  PXWSpectrumQuery(){
    spectrum.clear();
    start_scan=0;
    end_scan=0;
    retention_time_sec=0;
    precursor_neutral_mass=0;
    assumed_charge=0;
    searchHits=new vector<pxwSearchHitPair>;
  }
  PXWSpectrumQuery(const PXWSpectrumQuery& s){
    spectrum=s.spectrum;
    start_scan=s.start_scan;
    end_scan=s.end_scan;
    retention_time_sec=s.retention_time_sec;
    precursor_neutral_mass=s.precursor_neutral_mass;
    assumed_charge=s.assumed_charge;
    searchHits=new vector<pxwSearchHitPair>;
    for(size_t i=0;i<s.searchHits->size();i++) searchHits->push_back(s.searchHits->at(i));
  }
  ~PXWSpectrumQuery(){
    delete searchHits;
  }
  PXWSpectrumQuery& operator=(const PXWSpectrumQuery& s){
    if(this!=&s){
      spectrum=s.spectrum;
      start_scan=s.start_scan;
      end_scan=s.end_scan;
      retention_time_sec=s.retention_time_sec;
      precursor_neutral_mass=s.precursor_neutral_mass;
      assumed_charge=s.assumed_charge;
      delete searchHits;
      searchHits=new vector<pxwSearchHitPair>;
      for(size_t i=0;i<s.searchHits->size();i++) searchHits->push_back(s.searchHits->at(i));
    }
    return *this;
  }

  void addSearchHit(PXWSearchHit* s, PXWSearchHit* s2=NULL, string* xl=NULL, double* xlMass=NULL){
    pxwSearchHitPair p;
    p.a = new PXWSearchHit(*s);
    if(s2!=NULL) {
      if(xl==NULL || xlMass==NULL){
        printf("PXWSpectrumQuery.addSearchHit(): cross-linked peptides must contain linker and mass. Exiting.\n");
        exit(-4);
      }
      p.b = new PXWSearchHit(*s2);
    }
    if(xl!=NULL) p.identifier=*xl;
    if(xlMass!=NULL) p.mass=*xlMass;
    searchHits->push_back(p);
  }
  void clear(){
    searchHits->clear();
  }
  pxwSearchHitPair& getSearchHit(int index){
    return searchHits->at(index);
  }
  pxwSearchHitPair& getSearchHit(size_t index){
    return searchHits->at(index);
  }
  size_t sizeSearchHits(){
    return searchHits->size();
  }

private:
  vector<pxwSearchHitPair>* searchHits;
};

class PepXMLWriter {
public:

  PepXMLWriter();
  ~PepXMLWriter();

  void  closePepXML         ();
  bool  createPepXML        (char* fn, pxwMSMSRunSummary& run, pxwSampleEnzyme* enzyme=NULL, PXWSearchSummary* search=NULL);
  void  writeSpectrumQuery  (PXWSpectrumQuery& s);

private:

  void addTab               ();
  void deleteTab            ();
  void resetTabs            ();
  void writeAltProtein      (pxwProtein& s);
  void writeModAAMass       (pxwModAA& s);
  void writeModInfo         (PXWModInfo& s);
  void writeLine            (const char* str);
  void writeLinkedPeptide   (PXWSearchHit& s, bool alpha=true);
  void writeSearchHit       (pxwSearchHitPair& s);

  int index;
  int iTabs;
  int spQueryIndex;
  FILE* fptr;

  bool bTabs;
  bool bFileOpen;
  char strTabs[128];

  vector<string> vTagState;

};

#endif
