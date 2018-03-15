/* zran.c -- example of zlib/gzip stream indexing and random access
 * Copyright (C) 2005 Mark Adler
 * For conditions of distribution and use, see copyright notice in zlib.h
   Version 1.0  29 May 2005  Mark Adler */

/* Illustrate the use of Z_BLOCK, inflatePrime(), and inflateSetDictionary()
   for random access of a compressed file.  A file containing a zlib or gzip
   stream is provided on the command line.  The compressed stream is decoded in
   its entirety, and an index built with access points about every SPAN bytes
   in the uncompressed output.  The compressed file is left open, and can then
   be read randomly, having to decompress on the average SPAN/2 uncompressed
   bytes before getting to the desired block of data.

   An access point can be created at the start of any deflate block, by saving
   the starting file offset and bit of that block, and the 32K bytes of
   uncompressed data that precede that block.  Also the uncompressed offset of
   that block is saved to provide a referece for locating a desired starting
   point in the uncompressed stream.  build_index() works by decompressing the
   input zlib or gzip stream a block at a time, and at the end of each block
   deciding if enough uncompressed data has gone by to justify the creation of
   a new access point.  If so, that point is saved in a data structure that
   grows as needed to accommodate the points.

   To use the index, an offset in the uncompressed data is provided, for which
   the latest accees point at or preceding that offset is located in the index.
   The input file is positioned to the specified location in the index, and if
   necessary the first few bits of the compressed data is read from the file.
   inflate is initialized with those bits and the 32K of uncompressed data, and
   the decompression then proceeds until the desired offset in the file is
   reached.  Then the decompression continues to read the desired uncompressed
   data from the file.

   Another approach would be to generate the index on demand.  In that case,
   requests for random access reads from the compressed data would try to use
   the index, but if a read far enough past the end of the index is required,
   then further index entries would be generated and added.

   There is some fair bit of overhead to starting inflation for the random
   access, mainly copying the 32K byte dictionary.  So if small pieces of the
   file are being accessed, it would make sense to implement a cache to hold
   some lookahead and avoid many calls to extract() for small lengths.

   Another way to build an index would be to use inflateCopy().  That would
   not be constrained to have access points at block boundaries, but requires
   more memory per access point, and also cannot be saved to file due to the
   use of pointers in the state.  The approach here allows for storage of the
   index in a file.
 */

#include "mzParser.h"


Czran::Czran(){
	index=NULL;
	buffer=NULL;
	lastBuffer=NULL;
	bufferOffset=0;
	bufferLen=0;
	fileSize=0;
	lastBufferOffset=0;
}

Czran::~Czran(){
	if(index!=NULL) free_index();
	if(buffer!=NULL) free(buffer);
	if(lastBuffer!=NULL) free(lastBuffer);
	buffer=NULL;
	lastBuffer=NULL;
}

/* Deallocate an index built by build_index() */
void Czran::free_index(){
    if (index != NULL) {
        free(index->list);
        free(index);
				index=NULL;
    }
}

/* Add an entry to the access point list.  If out of memory, deallocate the
   existing list and return NULL. */
gz_access * Czran::addpoint(int bits,f_off in, f_off out, unsigned left, unsigned char *window) {
    point *next;

    /* if list is empty, create it (start with eight points) */
    if (index == NULL) {
        index = (gz_access*)malloc(sizeof(gz_access));
        if (index == NULL) return NULL;
        index->list = (point*)malloc(sizeof(point) << 3);
        if (index->list == NULL) {
            free(index);
            return NULL;
        }
        index->size = 8;
        index->have = 0;
    }

    /* if list is full, make it bigger */
    else if (index->have == index->size) {
        index->size <<= 1;
        next = (point*)realloc(index->list, sizeof(point) * index->size);
        if (next == NULL) {
            free_index();
            return NULL;
        }
        index->list = next;
    }

    /* fill in entry and increment how many we have */
    next = index->list + index->have;
    next->bits = bits;
    next->in = in;
    next->out = out;
    if (left)
        memcpy(next->window, window + WINSIZE - left, left);
    if (left < WINSIZE)
        memcpy(next->window + left, window, WINSIZE - left);
    index->have++;

    /* return list, possibly reallocated */
    return index;
}

/* Make one entire pass through the compressed stream and build an index, with
   access points about every span bytes of uncompressed output -- span is
   chosen to balance the speed of random access against the memory requirements
   of the list, about 32K bytes per access point.  Note that data after the end
   of the first zlib or gzip stream in the file is ignored.  build_index()
   returns the number of access points on success (>= 1), Z_MEM_ERROR for out
   of memory, Z_DATA_ERROR for an error in the input file, or Z_ERRNO for a
   file read error.  On success, *built points to the resulting index. */
int Czran::build_index(FILE *in, f_off span){
	return build_index(in,span,&index);
}
int Czran::build_index(FILE *in, f_off span, gz_access **built){
    int ret;
    f_off totin, totout;        /* our own total counters to avoid 4GB limit */
    f_off last;                 /* totout value of last access point */
    gz_access *index;							/* access points being generated */
    z_stream strm;
    unsigned char input[READCHUNK];
    unsigned char window[WINSIZE];

    /* initialize inflate */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit2(&strm, 47);      /* automatic zlib or gzip decoding */
    if (ret != Z_OK)
        return ret;

    /* inflate the input, maintain a sliding window, and build an index -- this
       also validates the integrity of the compressed data using the check
       information at the end of the gzip or zlib stream */
    totin = totout = last = 0;
    index = NULL;               /* will be allocated by first addpoint() */
    strm.avail_out = 0;
    do {
        /* get some compressed data from input file */
        strm.avail_in = fread(input, 1, READCHUNK, in);
        if (ferror(in)) {
            ret = Z_ERRNO;
            goto build_index_error;
        }
        if (strm.avail_in == 0) {
            ret = Z_DATA_ERROR;
            goto build_index_error;
        }
        strm.next_in = input;

        /* process all of that, or until end of stream */
        do {
            /* reset sliding window if necessary */
            if (strm.avail_out == 0) {
                strm.avail_out = WINSIZE;
                strm.next_out = window;
            }

            /* inflate until out of input, output, or at end of block --
               update the total input and output counters */
            totin += strm.avail_in;
            totout += strm.avail_out;
            ret = inflate(&strm, Z_BLOCK);      /* return at end of block */
            totin -= strm.avail_in;
            totout -= strm.avail_out;
            if (ret == Z_NEED_DICT)
                ret = Z_DATA_ERROR;
            if (ret == Z_MEM_ERROR || ret == Z_DATA_ERROR)
                goto build_index_error;
            if (ret == Z_STREAM_END)
                break;

            /* if at end of block, consider adding an index entry (note that if
               data_type indicates an end-of-block, then all of the
               uncompressed data from that block has been delivered, and none
               of the compressed data after that block has been consumed,
               except for up to seven bits) -- the totout == 0 provides an
               entry point after the zlib or gzip header, and assures that the
               index always has at least one access point; we avoid creating an
               access point after the last block by checking bit 6 of data_type
             */
            if ((strm.data_type & 128) && !(strm.data_type & 64) &&
                (totout == 0 || totout - last > span)) {
                index = addpoint(strm.data_type & 7, totin,
                                 totout, strm.avail_out, window);
                if (index == NULL) {
                    ret = Z_MEM_ERROR;
                    goto build_index_error;
                }
                last = totout;
            }
        } while (strm.avail_in != 0);
    } while (ret != Z_STREAM_END);

    /* clean up and return index (release unused entries in list) */
		fileSize=strm.total_out;
    (void)inflateEnd(&strm);
    index = (gz_access*)realloc(index, sizeof(point) * index->have);
    index->size = index->have;
    *built = index;
    return index->size;

    /* return error */
  build_index_error:
    (void)inflateEnd(&strm);
		fileSize=0;
    if (index != NULL)
        free_index();
    return ret;
}

/* Use the index to read len bytes from offset into buf, return bytes read or
   negative for error (Z_DATA_ERROR or Z_MEM_ERROR).  If data is requested past
   the end of the uncompressed data, then extract() will return a value less
   than len, indicating how much as actually read into buf.  This function
   should not return a data error unless the file was modified since the index
   was generated.  extract() may also return Z_ERRNO if there is an error on
   reading or seeking the input file. */
int Czran::extract(FILE *in, f_off offset) {

		int ret, len;
    point *here;
		z_stream strm;
    unsigned char input[READCHUNK];
    f_off marker;

		/* find where in stream to start */
		here = index->list;
		ret = index->have;
		while (--ret && here[1].out <= offset)
				here++;
    if(ret>0) marker=here[1].out;
    else marker=0;

		/* initialize file and inflate state to start there */
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		strm.avail_in = 0;
		strm.next_in = Z_NULL;
		ret = inflateInit2(&strm, -15);         /* raw inflate */
		if (ret != Z_OK)
				return ret;
		ret = mzpfseek(in, here->in - (here->bits ? 1 : 0), SEEK_SET);
		if (ret == -1)
				goto extract_ret;
		if (here->bits) {
				ret = getc(in);
				if (ret == -1) {
						ret = ferror(in) ? Z_ERRNO : Z_DATA_ERROR;
						goto extract_ret;
				}
				(void)inflatePrime(&strm, here->bits, ret >> (8 - here->bits));
		}
		(void)inflateSetDictionary(&strm, here->window, WINSIZE);

		if(marker>0) len = (int)(marker-here->out);
		else len = (int)(fileSize-here->out);

		if(buffer!=NULL) free(buffer);
		buffer = (unsigned char*)malloc(len);
		if(buffer==NULL){
			ret = Z_MEM_ERROR;
			goto extract_ret;
		}
		bufferOffset=here->out;

		strm.avail_in = 0;
		strm.avail_out = len;
    strm.next_out = buffer;

    /* uncompress until avail_out filled, or end of stream */
    do {
        if (strm.avail_in == 0) {
            strm.avail_in = fread(input, 1, READCHUNK, in);
            if (ferror(in)) {
                ret = Z_ERRNO;
                goto extract_ret;
            }
            if (strm.avail_in == 0) {
                ret = Z_DATA_ERROR;
                goto extract_ret;
            }
            strm.next_in = input;
        }
        ret = inflate(&strm, Z_NO_FLUSH);       /* normal inflate */
				if (ret == Z_NEED_DICT) 
            ret = Z_DATA_ERROR;
				if (ret == Z_MEM_ERROR || ret == Z_DATA_ERROR)
            goto extract_ret;
				if (ret == Z_STREAM_END)
            break;
    } while (strm.avail_out != 0);

    /* compute number of uncompressed bytes read after offset */
    ret = len - strm.avail_out;

    /* clean up and return bytes read or error */
  extract_ret:
		if(ret<0) bufferLen=0;
		else bufferLen=ret;
    (void)inflateEnd(&strm);
    return ret;

}

int Czran::extract(FILE *in, f_off offset, unsigned char *buf, int len){

	int ret, seg;

	//see if request was to the last buffer
	//last buffer is only stored when it doesn't reside in large buffer
	if(lastBuffer!=NULL && offset==lastBufferOffset){
		memcpy(buf,lastBuffer,lastBufferLen);
		return lastBufferLen;
	}

	//if we don't have the offset, grab it.
	if(buffer==NULL || offset<bufferOffset || offset>=(bufferOffset+bufferLen)){	
		ret=extract(in,offset);
	}

	lastBufferOffset=offset;

	//check if we have all the requested bytes
	if( (offset+len) <= (bufferOffset+bufferLen) ){

		memcpy(buf,buffer+(offset-bufferOffset),len);
		if(lastBuffer!=NULL){
			free(lastBuffer);
			lastBuffer=NULL;
		}
		return len;

	} else {

		if(lastBuffer!=NULL) free(lastBuffer);
		lastBuffer=(unsigned char*)malloc(len);

		//otherwise grab what we can and try extending buffer
		seg=(int)(bufferLen-(offset-bufferOffset));
		memcpy(buf,buffer+(offset-bufferOffset),seg);
		len-=seg;

		//get next block
		ret=extract(in,offset+seg);

		//add remaining bytes
		if(ret<len){
			memcpy(buf+seg,buffer,ret);
			ret+=seg;
		} else {
			memcpy(buf+seg,buffer,len);
			ret=seg+len;
		}

		if(ret!=lastBufferLen){
			if(lastBuffer!=NULL) free(lastBuffer);
			lastBuffer=(unsigned char*)malloc(ret);
			lastBufferLen=ret;
		}
		memcpy(lastBuffer,buf,ret);
		return ret;

	}

}

f_off Czran::getfilesize(){
	return fileSize;
}

