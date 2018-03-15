/*
mzpMz5Config - The code is
open source under the FreeBSD License, please see LICENSE file
for detailed information.

Copyright (C) 2011, Mike Hoopmann, Institute for Systems Biology
Version 1.0, January 4, 2011.
Version 1.1, March 14, 2012.
*/

#include "mzParser.h"
#ifdef MZP_MZ5
mzpMz5Config::mzpMz5Config(){
	init(false, true, true);
}

mzpMz5Config::~mzpMz5Config(){}

const bool mzpMz5Config::doFiltering() const {
	return doFiltering_;
}

const bool mzpMz5Config::doTranslating() const {
	return doTranslating_;
}

const size_t mzpMz5Config::getBufferInB() {
	return (bufferInMB_ * 1024L * 1024L);
}

const DataType& mzpMz5Config::getDataTypeFor(const MZ5DataSets v) {
	if (variableTypes_.find(v) != variableTypes_.end()) {
		return variableTypes_.find(v)->second;
	}
	throw out_of_range("[mzpMz5Config::getDataTypeFor]: out of range");
}

const string& mzpMz5Config::getNameFor(const MZ5DataSets v) {
	if (variableNames_.find(v) != variableNames_.end()) {
		return variableNames_.find(v)->second;
	}
	throw out_of_range("[mzpMz5Config::getNameFor]: out of range");
}

const size_t& mzpMz5Config::getRdccSlots() {
	return rdccSolts_;
}

MZ5DataSets mzpMz5Config::getVariableFor(const std::string& name) {
	if (variableVariables_.find(name) != variableVariables_.end()) {
		return variableVariables_.find(name)->second;
	}
	throw out_of_range("[mzpMz5Config::getVariableFor]: out of range");
}

void mzpMz5Config::init(const bool filter, const bool deltamz, const bool translateinten) {
	variableNames_.insert(pair<MZ5DataSets, string>(ControlledVocabulary, "ControlledVocabulary"));
	variableNames_.insert(pair<MZ5DataSets, string>(CVReference, "CVReference"));
	variableNames_.insert(pair<MZ5DataSets, string>(CVParam, "CVParam"));
	variableNames_.insert(pair<MZ5DataSets, string>(UserParam, "UserParam"));
	variableNames_.insert(pair<MZ5DataSets, string>(RefParam, "RefParam"));
	variableNames_.insert(pair<MZ5DataSets, string>(FileContent, "FileContent"));
	variableNames_.insert(pair<MZ5DataSets, string>(Contact, "Contact"));
	variableNames_.insert(pair<MZ5DataSets, string>(ParamGroups, "ParamGroups"));
	variableNames_.insert(pair<MZ5DataSets, string>(SourceFiles, "SourceFiles"));
	variableNames_.insert(pair<MZ5DataSets, string>(Samples, "Samples"));
	variableNames_.insert(pair<MZ5DataSets, string>(Software, "Software"));
	variableNames_.insert(pair<MZ5DataSets, string>(ScanSetting, "ScanSetting"));
	variableNames_.insert(pair<MZ5DataSets, string>(InstrumentConfiguration, "InstrumentConfiguration"));
	variableNames_.insert(pair<MZ5DataSets, string>(DataProcessing, "DataProcessing"));
	variableNames_.insert(pair<MZ5DataSets, string>(Run, "Run")); 
	variableNames_.insert(pair<MZ5DataSets, string>(SpectrumMetaData, "SpectrumMetaData"));
	variableNames_.insert(pair<MZ5DataSets, string>(SpectrumBinaryMetaData, "SpectrumListBinaryData"));
	variableNames_.insert(pair<MZ5DataSets, string>(ChromatogramMetaData, "ChromatogramList"));
	variableNames_.insert(pair<MZ5DataSets, string>(ChromatogramBinaryMetaData, "ChromatogramListBinaryData"));
	variableNames_.insert(pair<MZ5DataSets, string>(ChromatogramIndex, "ChromatogramIndex"));
	variableNames_.insert(pair<MZ5DataSets, string>(SpectrumIndex, "SpectrumIndex"));
	variableNames_.insert(pair<MZ5DataSets, string>(SpectrumMZ, "SpectrumMZ"));
	variableNames_.insert(pair<MZ5DataSets, string>(SpectrumIntensity, "SpectrumIntensity"));
	variableNames_.insert(pair<MZ5DataSets, string>(ChromatogramTime, "ChomatogramTime"));
	variableNames_.insert(pair<MZ5DataSets, string>(ChromatogramIntensity, "ChromatogramIntensity"));
	variableNames_.insert(pair<MZ5DataSets, string>(FileInformation, "FileInformation"));

	for (map<MZ5DataSets, string>::iterator it = variableNames_.begin(); it != variableNames_.end(); ++it) {
		variableVariables_.insert(pair<string, MZ5DataSets>(it->second, it->first));
	}

	variableTypes_.insert(pair<MZ5DataSets, DataType>(ControlledVocabulary, ContVocabMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(FileContent, ParamListMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(Contact, ParamListMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(CVReference, CVRefMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(ParamGroups, ParamGroupMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(SourceFiles, SourceFileMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(Samples, SampleMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(Software, SoftwareMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(ScanSetting, ScanSettingMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(InstrumentConfiguration, InstrumentConfigurationMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(DataProcessing, DataProcessingMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(Run, RunMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(SpectrumMetaData, SpectrumMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(SpectrumBinaryMetaData, BinaryDataMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(ChromatogramMetaData, ChromatogramMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(ChromatogramBinaryMetaData, BinaryDataMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(ChromatogramIndex, PredType::NATIVE_ULONG));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(SpectrumIndex, PredType::NATIVE_ULONG));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(FileInformation, FileInformationMZ5::getType()));

	variableTypes_.insert(pair<MZ5DataSets, DataType>(SpectrumMZ, PredType::NATIVE_DOUBLE));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(SpectrumIntensity, PredType::NATIVE_DOUBLE));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(ChromatogramIntensity, PredType::NATIVE_DOUBLE));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(ChromatogramTime, PredType::NATIVE_DOUBLE));

	variableTypes_.insert(pair<MZ5DataSets, DataType>(CVParam, CVParamMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(UserParam, UserParamMZ5::getType()));
	variableTypes_.insert(pair<MZ5DataSets, DataType>(RefParam, RefMZ5::getType()));

	hsize_t spectrumChunkSize = 5000L; // 1000=faster random read, 10000=better compression
	hsize_t chromatogramChunkSize = 1000L;
	hsize_t spectrumMetaChunkSize = 2000L; // should be modified in case of on demand access
	// hsize_t chromatogramMetaChunkSize = 10L; // usually one experiment does not contain a lot of chromatograms, so this chunk size is small in order to save storage space
	hsize_t cvparamChunkSize = 5000L;
	hsize_t userparamChunkSize = 100L;

	variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(SpectrumMZ, spectrumChunkSize));
	variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(SpectrumIntensity, spectrumChunkSize));
	variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(ChromatogramTime, chromatogramChunkSize));
	variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(ChromatogramIntensity, chromatogramChunkSize));

	bufferInMB_ = 8L; // this affects all datasets.
	size_t sizeOfDouble = static_cast<size_t> (sizeof(double));
	size_t bufferInByte = bufferInMB_ * 1024L * 1024L;
	size_t numberOfChunksInBuffer = (bufferInByte / sizeOfDouble) / spectrumChunkSize;
	rdccSolts_ = 41957L; // for 32 mb, 10000 chunk size

	hsize_t spectrumBufferSize = spectrumChunkSize * (numberOfChunksInBuffer / 4L);
	hsize_t chromatogramBufferSize = chromatogramChunkSize * 10L;

	variableBufferSizes_.insert(pair<MZ5DataSets, size_t>(SpectrumMZ, spectrumBufferSize));
	variableBufferSizes_.insert(pair<MZ5DataSets, size_t>(SpectrumIntensity, spectrumBufferSize));
	variableBufferSizes_.insert(pair<MZ5DataSets, size_t>(ChromatogramTime, chromatogramBufferSize));
	variableBufferSizes_.insert(pair<MZ5DataSets, size_t>(ChromatogramIntensity, chromatogramBufferSize));

	//if (config_.binaryDataEncoderConfig.compression == pwiz::msdata::BinaryDataEncoder::Compression_Zlib) {
		doTranslating_ = deltamz && translateinten;
		deflateLvl_ = 1;

		variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(SpectrumMetaData, spectrumMetaChunkSize));
		variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(SpectrumBinaryMetaData, spectrumMetaChunkSize));
		variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(SpectrumIndex, spectrumMetaChunkSize));
		// should not affect file size to much, if chromatogram information are not compressed
		// variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(ChromatogramMetaData, chromatogramMetaChunkSize));
		// variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(ChromatogramBinaryMetaData, chromatogramMetaChunkSize));
		// variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(ChromatogramIndex, chromatogramMetaChunkSize));
		variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(CVParam, cvparamChunkSize));
		variableChunkSizes_.insert(pair<MZ5DataSets, hsize_t>(UserParam, userparamChunkSize));
	//} else {
	//	deflateLvl_ = 0;
	//}

	spectrumLoadPolicy_ = SLP_InitializeAllOnFirstCall;
	chromatogramLoadPolicy_ = CLP_InitializeAllOnFirstCall;

	doFiltering_ = filter;
}

void mzpMz5Config::setFiltering(const bool flag) const {
	doFiltering_ = flag;
}

void mzpMz5Config::setTranslating(const bool flag) const {
	doTranslating_ = flag;
}

#endif
