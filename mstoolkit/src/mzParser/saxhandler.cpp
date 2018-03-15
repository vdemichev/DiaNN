/*
 Copyright (C) 2003-2004 Ronald C Beavis, all rights reserved
 X! tandem 
 This software is a component of the X! proteomics software
 development project

Use of this software governed by the Artistic license, as reproduced here:

The Artistic License for all X! software, binaries and documentation

Preamble
The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some 
semblance of artistic control over the development of the package, 
while giving the users of the package the right to use and distribute 
the Package in a more-or-less customary fashion, plus the right to 
make reasonable modifications. 

Definitions
"Package" refers to the collection of files distributed by the Copyright 
	Holder, and derivatives of that collection of files created through 
	textual modification. 

"Standard Version" refers to such a Package if it has not been modified, 
	or has been modified in accordance with the wishes of the Copyright 
	Holder as specified below. 

"Copyright Holder" is whoever is named in the copyright or copyrights 
	for the package. 

"You" is you, if you're thinking about copying or distributing this Package. 

"Reasonable copying fee" is whatever you can justify on the basis of 
	media cost, duplication charges, time of people involved, and so on. 
	(You will not be required to justify it to the Copyright Holder, but 
	only to the computing community at large as a market that must bear 
	the fee.) 

"Freely Available" means that no fee is charged for the item itself, 
	though there may be fees involved in handling the item. It also means 
	that recipients of the item may redistribute it under the same
	conditions they received it. 

1. You may make and give away verbatim copies of the source form of the 
Standard Version of this Package without restriction, provided that 
you duplicate all of the original copyright notices and associated 
disclaimers. 

2. You may apply bug fixes, portability fixes and other modifications 
derived from the Public Domain or from the Copyright Holder. A 
Package modified in such a way shall still be considered the Standard 
Version. 

3. You may otherwise modify your copy of this Package in any way, provided 
that you insert a prominent notice in each changed file stating how and 
when you changed that file, and provided that you do at least ONE of the 
following: 

a.	place your modifications in the Public Domain or otherwise make them 
	Freely Available, such as by posting said modifications to Usenet 
	or an equivalent medium, or placing the modifications on a major 
	archive site such as uunet.uu.net, or by allowing the Copyright Holder 
	to include your modifications in the Standard Version of the Package. 
b.	use the modified Package only within your corporation or organization. 
c.	rename any non-standard executables so the names do not conflict 
	with standard executables, which must also be provided, and provide 
	a separate manual page for each non-standard executable that clearly 
	documents how it differs from the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

4. You may distribute the programs of this Package in object code or 
executable form, provided that you do at least ONE of the following: 

a.	distribute a Standard Version of the executables and library files, 
	together with instructions (in the manual page or equivalent) on 
	where to get the Standard Version. 
b.	accompany the distribution with the machine-readable source of the 
	Package with your modifications. 
c.	give non-standard executables non-standard names, and clearly 
	document the differences in manual pages (or equivalent), together 
	with instructions on where to get the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

5. You may charge a reasonable copying fee for any distribution of 
this Package. You may charge any fee you choose for support of 
this Package. You may not charge a fee for this Package itself. 
However, you may distribute this Package in aggregate with other 
(possibly commercial) programs as part of a larger (possibly 
commercial) software distribution provided that you do not a
dvertise this Package as a product of your own. You may embed this 
Package's interpreter within an executable of yours (by linking); 
this shall be construed as a mere form of aggregation, provided that 
the complete Standard Version of the interpreter is so embedded. 

6. The scripts and library files supplied as input to or produced as 
output from the programs of this Package do not automatically fall 
under the copyright of this Package, but belong to whomever generated 
them, and may be sold commercially, and may be aggregated with this 
Package. If such scripts or library files are aggregated with this 
Package via the so-called "undump" or "unexec" methods of producing 
a binary executable image, then distribution of such an image shall 
neither be construed as a distribution of this Package nor shall it 
fall under the restrictions of Paragraphs 3 and 4, provided that you 
do not represent such an executable image as a Standard Version of 
this Package. 

7. C subroutines (or comparably compiled subroutines in other languages) 
supplied by you and linked into this Package in order to emulate 
subroutines and variables of the language defined by this Package 
shall not be considered part of this Package, but are the equivalent 
of input as in Paragraph 6, provided these subroutines do not change 
the language in any way that would cause it to fail the regression 
tests for the language. 

8. Aggregation of this Package with a commercial distribution is always 
permitted provided that the use of this Package is embedded; that is, 
when no overt attempt is made to make this Package's interfaces visible 
to the end user of the commercial distribution. Such use shall not be 
construed as a distribution of this Package. 

9. The name of the Copyright Holder may not be used to endorse or promote 
products derived from this software without specific prior written permission. 

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED 
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE. 

The End 
*/

#include "mzParser.h"

// Static callback handlers
static void mzp_startElementCallback(void *data, const XML_Char *el, const XML_Char **attr)
{
	((mzpSAXHandler*) data)->startElement(el, attr);
}

static void mzp_endElementCallback(void *data, const XML_Char *el)
{
	((mzpSAXHandler*) data)->endElement(el);
}

static void mzp_charactersCallback(void *data, const XML_Char *s, int len)
{
	((mzpSAXHandler*) data)->characters(s, len);
}

mzpSAXHandler::mzpSAXHandler()
{
	fptr = NULL;
	m_bGZCompression = false;
	fptr = NULL;
	m_parser = XML_ParserCreate(NULL);
	XML_SetUserData(m_parser, this);
	XML_SetElementHandler(m_parser, mzp_startElementCallback, mzp_endElementCallback);
	XML_SetCharacterDataHandler(m_parser, mzp_charactersCallback);
}


mzpSAXHandler::~mzpSAXHandler()
{
	if(fptr!=NULL) fclose(fptr);
	fptr = NULL;
	XML_ParserFree(m_parser);
}

void mzpSAXHandler::startElement(const XML_Char *el, const XML_Char **attr)
{
}

void mzpSAXHandler::endElement(const XML_Char *el)
{
}


void mzpSAXHandler::characters(const XML_Char *s, int len)
{
}

bool mzpSAXHandler::open(const char* fileName){
	if(fptr!=NULL) fclose(fptr);
	if(m_bGZCompression) fptr=fopen(fileName,"rb");
	else fptr=fopen(fileName,"r");
	if(fptr==NULL){
		//cerr << "Failed to open input file '" << fileName << "'.\n";
		return false;
	}
	setFileName(fileName);

	//Build the index if gz compressed
	if(m_bGZCompression){
		gzObj.free_index();

		int len;
		len = gzObj.build_index(fptr, SPAN);
    
		if (len < 0) {
        fclose(fptr);
        switch (len) {
        case Z_MEM_ERROR:
            fprintf(stderr, "Error reading .gz file: out of memory\n");
            break;
        case Z_DATA_ERROR:
            fprintf(stderr, "Error reading .gz file: compressed data error in %s\n", fileName);
            break;
        case Z_ERRNO:
            fprintf(stderr, "Error reading .gz file: read error on %s\n", fileName);
            break;
        default:
            fprintf(stderr, "Error reading .gz file: error %d while building index\n", len);
        }
				fptr=NULL;
        return false;
    }
	}

	return true;

}

bool mzpSAXHandler::parse()
{
	if (fptr == NULL){
		cerr << "Error parse(): No open file." << endl;
		return false;
	}

	char buffer[CHUNK];  //CHUNK=16384
	int readBytes = 0;
	bool success = true;
	int chunk=0;

	if(m_bGZCompression){
		while (success && (readBytes = gzObj.extract(fptr, 0+chunk*CHUNK, (unsigned char*)buffer, CHUNK))>0) {
			success = (XML_Parse(m_parser, buffer, readBytes, false) != 0);
			chunk++;
		}
	} else {
		while (success && (readBytes = (int) fread(buffer, 1, sizeof(buffer), fptr)) != 0){
			success = (XML_Parse(m_parser, buffer, readBytes, false) != 0);
		}
	}
	success = success && (XML_Parse(m_parser, buffer, 0, true) != 0);

	if (!success)
	{
		XML_Error error = XML_GetErrorCode(m_parser);

		cerr << m_strFileName
			<< "(" << XML_GetCurrentLineNumber(m_parser) << ")"
			<< " : error " << (int) error << ": ";

		switch (error)
		{
			case XML_ERROR_SYNTAX:
			case XML_ERROR_INVALID_TOKEN:
			case XML_ERROR_UNCLOSED_TOKEN:
				cerr << "Syntax error parsing XML.";
				break;

			// TODO: Add more descriptive text for interesting errors.

			default:
				cerr << "XML Parsing error.";
				break;
		} 
		cerr << "\n";
		return false;
	}
	return true;
}

//This function operates similarly to the parse() function.
//However, it accepts a file offset to begin parsing at a specific point.
//The parser will halt file reading when stop flag is triggered.
bool mzpSAXHandler::parseOffset(f_off offset){

	if (fptr == NULL){
		cerr << "Error parseOffset(): No open file." << endl;
		return false;
	}
	char buffer[CHUNK]; //CHUNK=16384
	int readBytes = 0;
	bool success = true;
	int chunk=0;
	
	XML_ParserReset(m_parser,"ISO-8859-1");
	XML_SetUserData(m_parser, this);
	XML_SetElementHandler(m_parser, mzp_startElementCallback, mzp_endElementCallback);
	XML_SetCharacterDataHandler(m_parser, mzp_charactersCallback);

	mzpfseek(fptr,offset,SEEK_SET);
	m_bStopParse=false;

	if(m_bGZCompression){
		while (success && (readBytes = gzObj.extract(fptr, offset+chunk*CHUNK, (unsigned char*)buffer, CHUNK))>0) {
			success = (XML_Parse(m_parser, buffer, readBytes, false) != 0);
			chunk++;
			if(m_bStopParse) break;
		}
	} else {
		while (success && (readBytes = (int) fread(buffer, 1, sizeof(buffer), fptr)) != 0) {
			success = (XML_Parse(m_parser, buffer, readBytes, false) != 0);
			if(m_bStopParse) break;
		}
	}

	if (!success && !m_bStopParse)
	{
		XML_Error error = XML_GetErrorCode(m_parser);

		cerr << m_strFileName
			<< "(" << XML_GetCurrentLineNumber(m_parser) << ")"
			<< " : parseOffset() " << (int) error << ": ";

		switch (error)
		{
			case XML_ERROR_SYNTAX:
			case XML_ERROR_INVALID_TOKEN:
			case XML_ERROR_UNCLOSED_TOKEN:
				cerr << "Syntax error parsing XML." << endl;
				break;

			// TODO: Add more descriptive text for interesting errors.

			default:
				cerr << "Spectrum XML Parsing error:\n";
				cerr << XML_ErrorString(error) << endl;
				break;
		} 
		exit(-7);
		return false;
	}
	return true;
}

void mzpSAXHandler::setGZCompression(bool b){
	m_bGZCompression=b;
}
