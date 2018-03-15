/*
mz5handler - The code is
open source under the FreeBSD License, please see LICENSE file
for detailed information.

Copyright (C) 2011, Mike Hoopmann, Institute for Systems Biology
Version 1.0, January 4, 2011.
Version 1.1, March 14, 2012.
*/

#include "mzParser.h"
#ifdef MZP_MZ5

mzpMz5Handler::mzpMz5Handler(mzpMz5Config* c, BasicSpectrum* s){
	config_=c;
	spec=s;
}

mzpMz5Handler::mzpMz5Handler(mzpMz5Config* c, BasicSpectrum* s, BasicChromatogram* bc){
	config_=c;
	spec=s;
	chromat=bc;
}

mzpMz5Handler::~mzpMz5Handler(){
	config_=NULL;
	spec=NULL;
	chromat=NULL;
}

void mzpMz5Handler::clean(const MZ5DataSets v, void* data, const size_t dsend) {
    hsize_t dim[1] =
    { static_cast<hsize_t> (dsend) };
    DataSpace dsp(1, dim);
    DataSet::vlenReclaim(data, config_->getDataTypeFor(v), dsp);
    free(data);
    data = 0;
    dsp.close();
}

vector<cMz5Index>* mzpMz5Handler::getChromatIndex(){
	return &m_vChromatIndex;
}

void mzpMz5Handler::getData(vector<double>& data, const MZ5DataSets v, const hsize_t start, const hsize_t end) {
	hsize_t scount = end - start;
	data.resize(scount);
	if (scount > 0) {
		map<MZ5DataSets, DataSet>::iterator it = bufferMap_.find(v);
		if (it == bufferMap_.end()) {
			DataSet ds = file_->openDataSet(config_->getNameFor(v));
			bufferMap_.insert(pair<MZ5DataSets, DataSet>(v, ds));
			it = bufferMap_.find(v);
		}
		DataSet dataset = it->second;
		DataSpace dataspace = dataset.getSpace();
		hsize_t offset[1];
		offset[0] = start;
		hsize_t count[1];
		count[0] = scount;
		dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

		hsize_t dimsm[1];
		dimsm[0] = scount;
		DataSpace memspace(1, dimsm);

		dataset.read(&data[0], PredType::NATIVE_DOUBLE, memspace, dataspace);

		if (v == SpectrumMZ && config_->doTranslating()) {
			size_t ms = data.size();
			double s = 0;
			for (size_t i = 0; i < ms; ++i) {
				data[i] = data[i] + s;
				s = data[i];
			}
		}
		if (v == SpectrumIntensity && config_->doTranslating()) {
			//there is no translating to do...
			//Translator_mz5::reverseTranslateIntensity(data);
		}
		memspace.close();
		dataspace.close();
	}
}

const map<MZ5DataSets, size_t>& mzpMz5Handler::getFields() {
	return fields_;
}

vector<cMz5Index>* mzpMz5Handler::getSpecIndex(){
	return &m_vIndex;
}

int mzpMz5Handler::highChromat() {
	return m_vChromatIndex.size();
}

int mzpMz5Handler::highScan() {
	if(m_vIndex.size()==0) return 0;
	return m_vIndex[m_vIndex.size()-1].scanNum;
}

int mzpMz5Handler::lowScan() {
	if(m_vIndex.size()==0) return 0;
	return m_vIndex[0].scanNum;
}

void mzpMz5Handler::processCVParams(unsigned long index){
	if(cvRef[cvParams_[index].typeCVRefID].group==0){ //MS:
		switch(cvRef[cvParams_[index].typeCVRefID].ref) {
			case 1000016:
				if(cvRef[cvParams_[index].unitCVRefID].ref==31) spec->setRTime((float)atof(cvParams_[index].value));
				else spec->setRTime((float)atof(cvParams_[index].value)/60.0f); //assume seconds if not minutes
				break;
			case 1000041:
				spec->setPrecursorCharge(atoi(cvParams_[index].value));
				break;
			case 1000042:
				spec->setPrecursorIntensity(atof(cvParams_[index].value));
				break;
			case 1000045:
				spec->setCollisionEnergy(atof(cvParams_[index].value));
				break;
			case 1000127:
				spec->setCentroid(true);
				break;
			case 1000285:
				spec->setTotalIonCurrent(atof(cvParams_[index].value));
				break;
			case 1000504:
				spec->setBasePeakMZ(atof(cvParams_[index].value));
				break;
			case 1000505:
				spec->setBasePeakIntensity(atof(cvParams_[index].value));
				break;
			case 1000511:
				spec->setMSLevel(atoi(cvParams_[index].value));
				break;	
			case 1000512:
				spec->setFilterLine(cvParams_[index].value);
				break;
			case 1000527:
				spec->setHighMZ(atof(cvParams_[index].value));
				break;
			case 1000528:
				spec->setLowMZ(atof(cvParams_[index].value));
				break;
			case 1000744:
				spec->setPrecursorMZ(atof(cvParams_[index].value));
				break;
			default:
				//cout << "Unknown/Unparsed CV: " << cvRef[cvParams_[index].typeCVRefID].group << ":" << cvRef[cvParams_[index].typeCVRefID].ref << endl;
				break;
		}
	} else { //unknown:
		//cout << "Unknown/Unparsed CV: " << cvRef[cvParams_[index].typeCVRefID].group << ":" << cvRef[cvParams_[index].typeCVRefID].ref << endl;
	}

}

bool mzpMz5Handler::readChromatogram(int num){
	if(chromat==NULL) return false;
	chromat->clear();

	//if no chromatogram was requested, grab the next one
	if(num<0)	posChromatIndex++;
	else posChromatIndex=num;
	if(posChromatIndex>=(int)m_vChromatIndex.size()) return false;
	
	//Read the chromatogram
	vector<double> time, inten;
	TimeIntensityPair tip;
	if(posChromatIndex==0){
		getData(time, ChromatogramTime, 0, (hsize_t)m_vChromatIndex[posChromatIndex].offset);
		getData(inten, ChromatogramIntensity, 0, (hsize_t)m_vChromatIndex[posChromatIndex].offset);
	} else {
		getData(time, ChromatogramTime, (hsize_t)m_vChromatIndex[posChromatIndex-1].offset, (hsize_t)m_vChromatIndex[posChromatIndex].offset);
		getData(inten, ChromatogramIntensity, (hsize_t)m_vChromatIndex[posChromatIndex-1].offset, (hsize_t)m_vChromatIndex[posChromatIndex].offset);
	}
	for(unsigned int i=0;i<time.size();i++){
		tip.time=time[i];
		tip.intensity=inten[i];
		chromat->addTIP(tip);
	}

	//Read the metadata
	unsigned long start=m_vChromatIndex[posChromatIndex].cvStart;
	unsigned long stop=start+m_vChromatIndex[posChromatIndex].cvLen;
	for(size_t i=start;i<stop;i++) processCVParams(i);
	chromat->setIDString(&m_vChromatIndex[posChromatIndex].idRef[0]);

	return true;
}

void* mzpMz5Handler::readDataSet(const MZ5DataSets v, size_t& dsend, void* ptr) {
	DataSet ds = file_->openDataSet(config_->getNameFor(v));
	DataSpace dsp = ds.getSpace();
	hsize_t start[1], end[1];
	dsp.getSelectBounds(start, end);
	dsend = (static_cast<size_t> (end[0])) + 1;
	DataType dt = config_->getDataTypeFor(v);
	if (ptr == 0) ptr = calloc(dsend, dt.getSize());
	ds.read(ptr, dt);
	dsp.close();
	ds.close();
	return ptr;
}

bool mzpMz5Handler::readFile(const string filename){

	FileCreatPropList fcparm = FileCreatPropList::DEFAULT;
	FileAccPropList faparm = FileAccPropList::DEFAULT;

	int mds_nelemts;
	size_t rdcc_nelmts, rdcc_nbytes;
	double rdcc_w0;
	faparm.getCache(mds_nelemts, rdcc_nelmts, rdcc_nbytes, rdcc_w0);
	//TODO do not set global buffer size, instead set dataset specific buffer size
	rdcc_nbytes = config_->getBufferInB();
	//TODO can be set to 1 if chunks that have been fully read/written will never be read/written again
	//rdcc_w0 = 1.0;
	rdcc_nelmts = config_->getRdccSlots();
	faparm.setCache(mds_nelemts, rdcc_nelmts, rdcc_nbytes, rdcc_w0);

	try {
		file_ = new H5File(filename, H5F_ACC_RDONLY, fcparm, faparm);
	} catch (FileIException&){
		return false;
	}
	closed_ = false;

	hsize_t start[1], end[1];
	size_t dsend = 0;
	DataSet dataset;
	DataSpace dataspace;
	string oname;
	MZ5DataSets v;
	for (hsize_t i = 0; i < file_->getNumObjs(); ++i) {
		oname = file_->getObjnameByIdx(i);
		dataset = file_->openDataSet(oname);
		dataspace = dataset.getSpace();
		dataspace.getSelectBounds(start, end);
		dsend = (static_cast<size_t> (end[0])) + 1;
		try {
			v = config_->getVariableFor(oname);
			fields_.insert(pair<MZ5DataSets, size_t>(v, dsend));
		} catch (out_of_range&) {
		}
		dataspace.close();
		dataset.close();
	}

	map<MZ5DataSets, size_t>::const_iterator it;
	it = fields_.find(FileInformation);
	if (it != fields_.end()) {
		DataSet ds = file_->openDataSet(config_->getNameFor(FileInformation));
		DataSpace dsp = ds.getSpace();
		hsize_t start[1], end[1];
		dsp.getSelectBounds(start, end);
		dsend = (static_cast<size_t> (end[0])) + 1;
		DataType dt = config_->getDataTypeFor(FileInformation);
		FileInformationMZ5* fi = (FileInformationMZ5*) (calloc(dsend, dt.getSize()));
		ds.read(fi, dt);
		dsp.close();
		ds.close();

		if (dsend == 1) {
			if (fi[0].majorVersion == MZ5_FILE_MAJOR_VERSION && fi[0].minorVersion == MZ5_FILE_MINOR_VERSION) {
				config_->setFiltering(fi[0].didFiltering > 0 ? true : false);
				config_->setTranslating(fi[0].deltaMZ && fi[0].translateInten);
			}
		}
		
		hsize_t dim[1] = { static_cast<hsize_t> (dsend) };
		DataSpace dspr(1, dim);
		DataSet::vlenReclaim(fi, config_->getDataTypeFor(FileInformation), dspr);
		free(fi);
		fi = 0;
		dspr.close();
	} else {
		it = fields_.find(Run);
		if (it == fields_.end()) {
			throw runtime_error("mzpMz5Handler::readFile(): given file is not mz5.");
			return false;
		}
	}

	//Read the CV Reference List
	size_t numberOfRef = fields_.find(CVReference)->second;
	CVRefMZ5* cvrl = (CVRefMZ5*) readDataSet(CVReference, dsend);
	cvRef.clear();
	CVRefItem cvr;
	for(int xx=0;xx<numberOfRef;xx++){
		if(strcmp(cvrl[xx].prefix,"MS")==0) cvr.group=0;
		else if(strcmp(cvrl[xx].prefix,"UO")==0) cvr.group=1;
		else cvr.group=2;
		cvr.ref=cvrl[xx].accession;
		cvRef.push_back(cvr);
		//cout << cvrl[xx].name << "\t" << cvrl[xx].prefix << ":" << cvrl[xx].accession << endl;
	}
	clean(CVReference, cvrl, dsend);
	
	//Read the CVParams
	size_t numberOfCV = fields_.find(CVParam)->second;
	cvParams_.resize(numberOfCV);
	readDataSet(CVParam, dsend, &cvParams_[0]);

	//Build the index
	m_vIndex.clear();
	m_scanIDXCount=0;
	size_t numberOfSpectra_ = fields_.find(SpectrumMetaData)->second;
	vector<unsigned long> index;
	index.resize(numberOfSpectra_);
	readDataSet(SpectrumIndex,dsend,&index[0]);
	
	BinaryDataMZ5* binaryParamsData_ = (BinaryDataMZ5*) calloc(numberOfSpectra_, sizeof(BinaryDataMZ5));
	readDataSet(SpectrumBinaryMetaData, dsend, binaryParamsData_);

	SpectrumMZ5* spectrumData_ = (SpectrumMZ5*) calloc(numberOfSpectra_, sizeof(SpectrumMZ5));
	readDataSet(SpectrumMetaData, dsend, spectrumData_);
	
	for(size_t i=0;i<dsend;i++){
		/*
		if(i==0 || i==3237 || i==3238 || i==3239){
			cout << i << "\t" << spectrumData_[i].precursorList.len << endl;
			cout << i << "\t" << spectrumData_[i].paramList.cvParamStartID << "\t" << spectrumData_[i].paramList.cvParamEndID << endl;
			for(size_t j=spectrumData_[i].paramList.cvParamStartID;j<binaryParamsData_[i].xParamList.cvParamEndID;j++){
				cout << cvParams_[j].typeCVRefID << "\t" << cvParams_[j].unitCVRefID << "\t" << cvParams_[j].value << endl;
			}
			cout << i << "\t" << binaryParamsData_[i].xParamList.cvParamStartID << "\t" << binaryParamsData_[i].xParamList.cvParamEndID << endl;
			cout << i << "\t" << binaryParamsData_[i].yParamList.cvParamStartID << "\t" << binaryParamsData_[i].yParamList.cvParamEndID << endl;
			cout << endl;
		}
		*/
		curIndex.cvStart=spectrumData_[i].paramList.cvParamStartID;
		curIndex.cvLen=binaryParamsData_[i].xParamList.cvParamEndID-curIndex.cvStart;
		curIndex.idRef=spectrumData_[i].id;
		curIndex.offset=index[i];
		if(strstr(&curIndex.idRef[0],"scan=")!=NULL)	{
			curIndex.scanNum=atoi(strstr(&curIndex.idRef[0],"scan=")+5);
		} else if(strstr(&curIndex.idRef[0],"scanId=")!=NULL) {
			curIndex.scanNum=atoi(strstr(&curIndex.idRef[0],"scanId=")+7);
		} else if(strstr(&curIndex.idRef[0],"S")!=NULL) {
			curIndex.scanNum=atoi(strstr(&curIndex.idRef[0],"S")+1);
		} else {
			curIndex.scanNum=++m_scanIDXCount;
		}
		m_vIndex.push_back(curIndex);
	}

	//Build the chromatogram index
	m_vChromatIndex.clear();
	size_t numberOfChromats_ = fields_.find(ChromatogramMetaData)->second;
	index.clear();
	index.resize(numberOfChromats_);
	readDataSet(ChromatogramIndex,dsend,&index[0]);
	
	free(binaryParamsData_);
	binaryParamsData_ = (BinaryDataMZ5*) calloc(numberOfChromats_, sizeof(BinaryDataMZ5));
	readDataSet(ChromatogramBinaryMetaData, dsend, binaryParamsData_);

	ChromatogramMZ5* chromatogramData_ = (ChromatogramMZ5*) calloc(numberOfSpectra_, sizeof(ChromatogramMZ5));
	readDataSet(ChromatogramMetaData, dsend, chromatogramData_);
	
	for(size_t i=0;i<dsend;i++){
		/*
		if(i==0 || i==3237 || i==3238 || i==3239){
			cout << i << "\t" << spectrumData_[i].precursorList.len << endl;
			cout << i << "\t" << spectrumData_[i].paramList.cvParamStartID << "\t" << spectrumData_[i].paramList.cvParamEndID << endl;
			for(size_t j=spectrumData_[i].paramList.cvParamStartID;j<binaryParamsData_[i].xParamList.cvParamEndID;j++){
				cout << cvParams_[j].typeCVRefID << "\t" << cvParams_[j].unitCVRefID << "\t" << cvParams_[j].value << endl;
			}
			cout << i << "\t" << binaryParamsData_[i].xParamList.cvParamStartID << "\t" << binaryParamsData_[i].xParamList.cvParamEndID << endl;
			cout << i << "\t" << binaryParamsData_[i].yParamList.cvParamStartID << "\t" << binaryParamsData_[i].yParamList.cvParamEndID << endl;
			cout << endl;
		}
		*/
		curChromatIndex.cvStart=chromatogramData_[i].paramList.cvParamStartID;
		curChromatIndex.cvLen=binaryParamsData_[i].xParamList.cvParamEndID-curChromatIndex.cvStart;
		curChromatIndex.idRef=chromatogramData_[i].id;
		curChromatIndex.offset=index[i];
		m_vChromatIndex.push_back(curChromatIndex);
	}

	free(spectrumData_);
	free(chromatogramData_);
	free(binaryParamsData_);

	posIndex=-1;
	posChromatIndex=-1;

	return true;

}

bool mzpMz5Handler::readHeader(int num){
	spec->clear();

	//if no scan was requested, grab the next one
	if(num<0){
		posIndex++;
		if(posIndex>=(int)m_vIndex.size()) return false;
	} else { //otherwise do binary search for scan number
		//Assumes scan numbers are in order
		int mid=m_vIndex.size()/2;
		int upper=m_vIndex.size();
		int lower=0;
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
		if(m_vIndex[mid].scanNum==num) posIndex=mid;
	}

	//Read the metadata
	unsigned long start=m_vIndex[posIndex].cvStart;
	unsigned long stop=start+m_vIndex[posIndex].cvLen;
	for(size_t i=start;i<stop;i++) processCVParams(i);

	//Get the peak count - this wastes time, though
	vector<double> mz;
	if(posIndex==0) getData(mz, SpectrumMZ, 0, (hsize_t)m_vIndex[posIndex].offset);
	else getData(mz, SpectrumMZ, (hsize_t)m_vIndex[posIndex-1].offset, (hsize_t)m_vIndex[posIndex].offset);
	spec->setPeaksCount((int)mz.size());

	if(spec->getScanNum()!=m_vIndex[posIndex].scanNum) spec->setScanNum(m_vIndex[posIndex].scanNum);
	spec->setScanIndex(posIndex+1); //set the index, which starts from 1, so offset by 1

	return true;

}

bool mzpMz5Handler::readSpectrum(int num){
	spec->clear();

	//if no scan was requested, grab the next one
	if(num<0){
		posIndex++;
		if(posIndex>=(int)m_vIndex.size()) return false;
	} else { //otherwise do binary search for scan number
		//Assumes scan numbers are in order
		int mid=m_vIndex.size()/2;
		int upper=m_vIndex.size();
		int lower=0;
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
		if(m_vIndex[mid].scanNum==num) posIndex=mid;
	}

	//Read the peaks
	vector<double> mz, inten;
	specDP dp;
	if(posIndex==0){
		getData(mz, SpectrumMZ, 0, (hsize_t)m_vIndex[posIndex].offset);
		getData(inten, SpectrumIntensity, 0, (hsize_t)m_vIndex[posIndex].offset);
	} else {
		getData(mz, SpectrumMZ, (hsize_t)m_vIndex[posIndex-1].offset, (hsize_t)m_vIndex[posIndex].offset);
		getData(inten, SpectrumIntensity, (hsize_t)m_vIndex[posIndex-1].offset, (hsize_t)m_vIndex[posIndex].offset);
	}
	for(unsigned int i=0;i<mz.size();i++){
		dp.mz=mz[i];
		dp.intensity=inten[i];
		spec->addDP(dp);
	}
	spec->setPeaksCount((int)mz.size());

	//Read the metadata
	unsigned long start=m_vIndex[posIndex].cvStart;
	unsigned long stop=start+m_vIndex[posIndex].cvLen;
	for(size_t i=start;i<stop;i++) processCVParams(i);

	if(spec->getScanNum()!=m_vIndex[posIndex].scanNum) spec->setScanNum(m_vIndex[posIndex].scanNum);
	spec->setScanIndex(posIndex+1); //set the index, which starts from 1, so offset by 1

	return true;

}

#endif
