/*
mzpMz5Structs - The code is
open source under the FreeBSD License, please see LICENSE file
for detailed information.

Copyright (C) 2011, Mike Hoopmann, Institute for Systems Biology
Version 1.0, January 4, 2011.
Version 1.1, March 14, 2012.
*/

#include "mzParser.h"
#ifdef MZP_MZ5
StrType getStringType() {
	StrType stringtype(PredType::C_S1, H5T_VARIABLE);
	return stringtype;
}

StrType getFStringType(const size_t len) {
	StrType stringtype(PredType::C_S1, len);
	return stringtype;
}

CompType BinaryDataMZ5::getType() {
	CompType ret(sizeof(BinaryDataMZ5));
	size_t offset = 0;
	ret.insertMember("xParams", offset, ParamListMZ5::getType());
	offset += ParamListMZ5::getType().getSize();
	ret.insertMember("yParams", offset, ParamListMZ5::getType());
	offset += ParamListMZ5::getType().getSize();
	ret.insertMember("xrefDataProcessing", offset, RefMZ5::getType());
	offset += RefMZ5::getType().getSize();
	ret.insertMember("yrefDataProcessing", offset, RefMZ5::getType());
	return ret;
}

CompType ChromatogramMZ5::getType() {
	CompType ret(sizeof(ChromatogramMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("id", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("params", offset, ParamListMZ5::getType());
	offset += sizeof(ParamListMZ5Data);
	ret.insertMember("precursor", offset, PrecursorMZ5::getType());
	offset += sizeof(PrecursorMZ5);
	ret.insertMember("productIsolationWindow", offset, ParamListMZ5::getType());
	offset += sizeof(ParamListMZ5Data);
	ret.insertMember("refDataProcessing", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5Data);
	ret.insertMember("index", offset, PredType::NATIVE_ULONG);
	offset += sizeof(unsigned long);
	return ret;
}

VarLenType ComponentListMZ5::getType() {
	CompType c = ComponentMZ5::getType();
	VarLenType ret(&c);
	return ret;
}

CompType ComponentMZ5::getType() {
	CompType ret(sizeof(ComponentMZ5));
	size_t offset = 0;
	ret.insertMember("paramList", offset, ParamListMZ5::getType());
	offset += sizeof(ParamListMZ5Data);
	ret.insertMember("order", offset, PredType::NATIVE_ULONG);
	return ret;
}

CompType ComponentsMZ5::getType() {
	CompType ret(sizeof(ComponentsMZ5));
	size_t offset = 0;
	ret.insertMember("sources", offset, ComponentListMZ5::getType());
	offset += sizeof(ComponentListMZ5);
	ret.insertMember("analyzers", offset, ComponentListMZ5::getType());
	offset += sizeof(ComponentListMZ5);
	ret.insertMember("detectors", offset, ComponentListMZ5::getType());
	return ret;
}

CompType ContVocabMZ5::getType() {
	CompType cvtype(sizeof(ContVocabMZ5Data));
	StrType stringtype = getStringType();
	cvtype.insertMember("uri", HOFFSET(ContVocabMZ5Data, uri), stringtype);
	cvtype.insertMember("fullname", HOFFSET(ContVocabMZ5Data, fullname), stringtype);
	cvtype.insertMember("id", HOFFSET(ContVocabMZ5Data, id), stringtype);
	cvtype.insertMember("version", HOFFSET(ContVocabMZ5Data, version), stringtype);
	return cvtype;
}

CVParamMZ5::CVParamMZ5() {
	init(0, ULONG_MAX, ULONG_MAX);
}

CVParamMZ5::CVParamMZ5(const CVParamMZ5& cvparam) {
	init(cvparam.value, cvparam.typeCVRefID, cvparam.unitCVRefID);
}

CVParamMZ5& CVParamMZ5::operator=(const CVParamMZ5& rhs) {
	if (this != &rhs) init(rhs.value, rhs.typeCVRefID, rhs.unitCVRefID);
	return *this;
}

CVParamMZ5::~CVParamMZ5(){}

CompType CVParamMZ5::getType() {
	CompType ret(sizeof(CVParamMZ5Data));
	StrType stringtype = getFStringType(CVL);
	ret.insertMember("value", HOFFSET(CVParamMZ5Data, value), stringtype);
	ret.insertMember("cvRefID", HOFFSET(CVParamMZ5Data, typeCVRefID), PredType::NATIVE_ULONG);
	ret.insertMember("uRefID", HOFFSET(CVParamMZ5Data, unitCVRefID), PredType::NATIVE_ULONG);
	return ret;
}

void CVParamMZ5::init(const char* value, const unsigned long& cvrefid, const unsigned long& urefid) {
	if (value) strcpy(this->value, value);
	else this->value[0] = '\0';
	this->value[CVL - 1] = '\0';
	this->typeCVRefID = cvrefid;
	this->unitCVRefID = urefid;
}

CVRefMZ5::CVRefMZ5() {
	init("","",ULONG_MAX);
}

CVRefMZ5::CVRefMZ5(const CVRefMZ5& cvref) {
	init(cvref.name, cvref.prefix, cvref.accession);
}

CVRefMZ5::~CVRefMZ5() {
	delete [] name;
	delete [] prefix;
}

CVRefMZ5& CVRefMZ5::operator=(const CVRefMZ5& cvp) {
	if (this != &cvp){
		delete [] name;
		delete [] prefix;
		init(cvp.name, cvp.prefix, cvp.accession);
	}
	return *this;
}

void CVRefMZ5::init(const char* name, const char* prefix,	const unsigned long accession){
	if(name) strcpy(this->name,name);
	strcpy(this->prefix,prefix);
	this->accession=accession;
}

CompType CVRefMZ5::getType() {
	CompType ret(sizeof(CVRefMZ5Data));
	StrType stringtype = getStringType();
	ret.insertMember("name", HOFFSET(CVRefMZ5Data, name), stringtype);
	ret.insertMember("prefix", HOFFSET(CVRefMZ5Data, prefix), stringtype);
	ret.insertMember("accession", HOFFSET(CVRefMZ5Data, accession),	PredType::NATIVE_ULONG);
	return ret;
}

CompType DataProcessingMZ5::getType() {
	CompType ret(sizeof(DataProcessingMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("id", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("method", offset, ProcessingMethodListMZ5::getType());
	offset += sizeof(ProcessingMethodListMZ5);
	return ret;
}

FileInformationMZ5::FileInformationMZ5() {
	this->majorVersion = MZ5_FILE_MAJOR_VERSION;
	this->minorVersion = MZ5_FILE_MINOR_VERSION;
	this->didFiltering = 1;
	this->deltaMZ = 1;
	this->translateInten = 1;
}

FileInformationMZ5::FileInformationMZ5(const FileInformationMZ5& rhs) {
	init(rhs.majorVersion, rhs.minorVersion, rhs.didFiltering, rhs.deltaMZ, rhs.translateInten);
}

FileInformationMZ5::FileInformationMZ5(const mzpMz5Config& c) {
	unsigned short didfiltering = c.doFiltering() ? 1 : 0;
	unsigned short deltamz = c.doTranslating() ? 1 : 0;
	unsigned short translateinten = c.doTranslating() ? 1 : 0;
	init(MZ5_FILE_MAJOR_VERSION, MZ5_FILE_MINOR_VERSION, didfiltering, deltamz,	translateinten);
}

FileInformationMZ5::~FileInformationMZ5(){}

FileInformationMZ5& FileInformationMZ5::operator=(const FileInformationMZ5& rhs) {
	if (this != &rhs){
		init(rhs.majorVersion, rhs.minorVersion, rhs.didFiltering, rhs.deltaMZ, rhs.translateInten);
	}
	return *this;
}

void FileInformationMZ5::init(const unsigned short majorVersion, const unsigned short minorVersion, const unsigned didFiltering, const unsigned deltaMZ, const unsigned translateInten) {
	this->majorVersion = majorVersion;
	this->minorVersion = minorVersion;
	this->didFiltering = didFiltering;
	this->deltaMZ = deltaMZ;
	this->translateInten = translateInten;
}

CompType FileInformationMZ5::getType() {
	CompType ret(sizeof(FileInformationMZ5Data));
	ret.insertMember("majorVersion", HOFFSET(FileInformationMZ5Data, majorVersion), PredType::NATIVE_USHORT);
	ret.insertMember("minorVersion", HOFFSET(FileInformationMZ5Data, minorVersion), PredType::NATIVE_USHORT);
	ret.insertMember("didFiltering", HOFFSET(FileInformationMZ5Data, didFiltering), PredType::NATIVE_USHORT);
	ret.insertMember("deltaMZ", HOFFSET(FileInformationMZ5Data, deltaMZ), PredType::NATIVE_USHORT);
	ret.insertMember("translateInten", HOFFSET(FileInformationMZ5Data, translateInten), PredType::NATIVE_USHORT);
	return ret;
}

CompType InstrumentConfigurationMZ5::getType() {
	CompType ret(sizeof(InstrumentConfigurationMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("id", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("params", offset, ParamListMZ5::getType());
	offset += sizeof(ParamListMZ5Data);
	ret.insertMember("components", offset, ComponentsMZ5::getType());
	offset += sizeof(ComponentsMZ5);
	ret.insertMember("refScanSetting", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5Data);
	ret.insertMember("refSoftware", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5Data);
	return ret;
}

CompType ParamGroupMZ5::getType() {
	CompType ret(sizeof(ParamGroupMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("id", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("params", offset, ParamListMZ5::getType());
	return ret;
}

CompType ParamListMZ5::getType() {
	CompType ret(sizeof(ParamListMZ5Data));
	ret.insertMember("cvstart", HOFFSET(ParamListMZ5Data, cvParamStartID), PredType::NATIVE_ULONG);
	ret.insertMember("cvend", HOFFSET(ParamListMZ5Data, cvParamEndID), PredType::NATIVE_ULONG);
	ret.insertMember("usrstart", HOFFSET(ParamListMZ5Data, userParamStartID), PredType::NATIVE_ULONG);
	ret.insertMember("usrend", HOFFSET(ParamListMZ5Data, userParamEndID), PredType::NATIVE_ULONG);
	ret.insertMember("refstart", HOFFSET(ParamListMZ5Data, refParamGroupStartID), PredType::NATIVE_ULONG);
	ret.insertMember("refend", HOFFSET(ParamListMZ5Data, refParamGroupEndID), PredType::NATIVE_ULONG);
	return ret;
}

VarLenType ParamListsMZ5::getType() {
	CompType c(ParamListMZ5::getType());
	VarLenType ret(&c);
	return ret;
}

VarLenType PrecursorListMZ5::getType() {
	CompType c(PrecursorMZ5::getType());
	VarLenType ret(&c);
	return ret;
}

CompType PrecursorMZ5::getType() {
	CompType ret(sizeof(PrecursorMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("externalSpectrumId", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("activation", offset, ParamListMZ5::getType());
	offset += ParamListMZ5::getType().getSize();
	ret.insertMember("isolationWindow", offset, ParamListMZ5::getType());
	offset += ParamListMZ5::getType().getSize();
	ret.insertMember("selectedIonList", offset, ParamListsMZ5::getType());
	offset += ParamListsMZ5::getType().getSize();
	ret.insertMember("refSpectrum", offset, RefMZ5::getType());
	offset += RefMZ5::getType().getSize();
	ret.insertMember("refSourceFile", offset, RefMZ5::getType());
	offset += RefMZ5::getType().getSize();
	return ret;
}

VarLenType ProcessingMethodListMZ5::getType() {
	CompType c = ProcessingMethodMZ5::getType();
	VarLenType ret(&c);
	return ret;
}

CompType ProcessingMethodMZ5::getType() {
	CompType ret(sizeof(ProcessingMethodMZ5));
	size_t offset = 0;
	ret.insertMember("params", offset, ParamListMZ5::getType());
	offset += sizeof(ParamListMZ5Data);
	ret.insertMember("refSoftware", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5Data);
	ret.insertMember("order", offset, PredType::NATIVE_ULONG);
	return ret;
}

VarLenType RefListMZ5::getType() {
	CompType c = RefMZ5::getType();
	VarLenType ret(&c);
	return ret;
}

CompType RefMZ5::getType() {
	CompType ret(sizeof(RefMZ5Data));
	ret.insertMember("refID", HOFFSET(RefMZ5Data, refID), PredType::NATIVE_ULONG);
	return ret;
}

CompType RunMZ5::getType() {
	CompType ret(sizeof(RunMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("id", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("startTimeStamp", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("fid", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("facc", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("params", offset, ParamListMZ5::getType());
	offset += sizeof(ParamListMZ5);
	ret.insertMember("refSpectrumDP", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5);
	ret.insertMember("refChromatogramDP", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5);
	ret.insertMember("refDefaultInstrument", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5);
	ret.insertMember("refSourceFile", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5);
	ret.insertMember("refSample", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5);
	return ret;
}

CompType SampleMZ5::getType() {
	CompType ret(sizeof(SampleMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("id", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("name", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("params", offset, ParamListMZ5::getType());
	return ret;
}

VarLenType ScanListMZ5::getType() {
	CompType c = ScanMZ5::getType();
	VarLenType ret(&c);
	return ret;
}

CompType ScanMZ5::getType() {
	CompType ret(sizeof(ScanMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("externalSpectrumID", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("params", offset, ParamListMZ5::getType());
	offset += sizeof(ParamListMZ5Data);
	ret.insertMember("scanWindowList", offset, ParamListsMZ5::getType());
	offset += sizeof(ParamListsMZ5);
	ret.insertMember("refInstrumentConfiguration", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5Data);
	ret.insertMember("refSourceFile", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5Data);
	ret.insertMember("refSpectrum", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5Data);
	return ret;
}

CompType ScansMZ5::getType() {
	CompType ret(sizeof(ScansMZ5));
	size_t offset = 0;
	ret.insertMember("params", offset, ParamListMZ5::getType());
	offset += sizeof(ParamListMZ5);
	ret.insertMember("scanList", offset, ScanListMZ5::getType());
	return ret;
}

CompType ScanSettingMZ5::getType() {
	CompType ret(sizeof(ScanSettingMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("id", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("params", offset, ParamListMZ5::getType());
	offset += sizeof(ParamListMZ5Data);
	ret.insertMember("refSourceFiles", offset, RefListMZ5::getType());
	offset += sizeof(RefListMZ5Data);
	ret.insertMember("targets", offset, ParamListsMZ5::getType());
	return ret;
}

CompType SoftwareMZ5::getType() {
	CompType ret(sizeof(SoftwareMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("id", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("version", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("params", offset, ParamListMZ5::getType());
	return ret;
}

CompType SourceFileMZ5::getType() {
	CompType ret(sizeof(SourceFileMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("id", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("location", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("name", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("params", offset, ParamListMZ5::getType());
	return ret;
}

CompType SpectrumMZ5::getType() {
	CompType ret(sizeof(SpectrumMZ5));
	StrType stringtype = getStringType();
	size_t offset = 0;
	ret.insertMember("id", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("spotID", offset, stringtype);
	offset += stringtype.getSize();
	ret.insertMember("params", offset, ParamListMZ5::getType());
	offset += sizeof(ParamListMZ5Data);
	ret.insertMember("scanList", offset, ScansMZ5::getType());
	offset += sizeof(ScansMZ5);
	ret.insertMember("precursors", offset, PrecursorListMZ5::getType());
	offset += sizeof(PrecursorListMZ5);
	ret.insertMember("products", offset, ParamListsMZ5::getType());
	offset += sizeof(ParamListsMZ5);
	ret.insertMember("refDataProcessing", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5Data);
	ret.insertMember("refSourceFile", offset, RefMZ5::getType());
	offset += sizeof(RefMZ5Data);
	ret.insertMember("index", offset, PredType::NATIVE_ULONG);
	offset += PredType::NATIVE_ULONG.getSize();
	return ret;
}

CompType UserParamMZ5::getType() {
	CompType ret(sizeof(UserParamMZ5Data));
	StrType namestringtype = getFStringType(USRNL);
	StrType valuestringtype = getFStringType(USRVL);
	StrType typestringtype = getFStringType(USRTL);
	ret.insertMember("name", HOFFSET(UserParamMZ5Data, name), namestringtype);
	ret.insertMember("value", HOFFSET(UserParamMZ5Data, value), valuestringtype);
	ret.insertMember("type", HOFFSET(UserParamMZ5Data, type), typestringtype);
	ret.insertMember("uRefID", HOFFSET(UserParamMZ5Data, unitCVRefID), PredType::NATIVE_ULONG);
	return ret;
}
#endif
