/* downloaded from web */
//#include <cmath>
//#include <stdio.h>
//#include "stdafx.h"
//#include <cstring>

//#include "base64.h"
#include "mzParser.h"

inline int getPosition( char buf );

static const unsigned char *b64_tbl = (const unsigned char*) "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
static const unsigned char b64_pad = '=';

/* base64 encode a group of between 1 and 3 input chars into a group of  4 output chars */
static void encode_group (unsigned char output[],
						  const unsigned char input[],
						  int n)
{
	unsigned char ingrp[3];

	ingrp[0] = n > 0 ? input[0] : 0;
	ingrp[1] = n > 1 ? input[1] : 0;
	ingrp[2] = n > 2 ? input[2] : 0;

	/* upper 6 bits of ingrp[0] */
	output[0] = n > 0 ? b64_tbl[ingrp[0] >> 2] : b64_pad;

	/* lower 2 bits of ingrp[0] | upper 4 bits of ingrp[1] */
	output[1] = n > 0 ? b64_tbl[((ingrp[0] & 0x3) << 4) | (ingrp[1] >> 4)] : b64_pad;

	/* lower 4 bits of ingrp[1] | upper 2 bits of ingrp[2] */
	output[2] = n > 1 ? b64_tbl[((ingrp[1] & 0xf) << 2) | (ingrp[2] >> 6)] : b64_pad;

	/* lower 6 bits of ingrp[2] */
	output[3] = n > 2 ? b64_tbl[ingrp[2] & 0x3f] : b64_pad;

}

int b64_encode (char *dest,
                const char *src,
                int len)
{
  int outsz = 0;

  while (len > 0)
  {
    encode_group ( (unsigned char*) dest + outsz, (const unsigned char*) src, len > 3 ? 3 : len);
    len -= 3;
    src += 3;
    outsz += 4;
  }

  return outsz;

}



/* base64 decode a group of 4 input chars into a group of between 0 and
* 3 output chars */
static void decode_group (unsigned char output[],
						  const unsigned char input[],
						  int *n)
{
	unsigned char *t1, *t2;
	*n = 0;

	if (input[0] == '=')
		return;

	t1 = (unsigned char*) strchr ((const char*)b64_tbl, input[0]);
	t2 = (unsigned char*) strchr ((const char*)b64_tbl, input[1]);

	output[(*n)++] = (unsigned char)(((t1 - b64_tbl) << 2) | ((t2 - b64_tbl) >> 4));

	if (input[2] == '=')
		return;

	t1 = (unsigned char*) strchr ((const char*)b64_tbl, input[2]);

	output[(*n)++] = (unsigned char)(((t2 - b64_tbl) << 4) | ((t1 - b64_tbl) >> 2));

	if (input[3] == '=')
		return;

	t2 = (unsigned char*) strchr ((const char*)b64_tbl, input[3]);

	output[(*n)++] = (unsigned char)(((t1 - b64_tbl) << 6) | (t2 - b64_tbl));

	return;
}



inline int getPosition( char buf )
{

	if( buf > 96 )		// [a-z]
		return (buf - 71);
	else if( buf > 64 )		// [A-Z]
		return (buf - 65);
	else if( buf > 47 )		// [0-9]
		return (buf + 4);
	else if( buf == 43 )
		return 63;
	else				// buf == '/'
		return 64;
}


// Returns the total number of bytes decoded
int b64_decode_mio ( char *dest,  char *src, size_t size )
{
	char *temp = dest;
	char *end = dest + size;

	for (;;)
	{
		int register a;
		int register b;
		int t1,t2,t3,t4;

		if (!(t1 = *src++) || !(t2 = *src++) || !(t3 = *src++) || !(t4 = *src++))
			return (int)(temp-dest);

		if (t1 == 61 || temp >= end)		// if == '='
			return(int)(temp-dest);

		if( t1 > 96 )		// [a-z]
			a = (t1 - 71);
		else if( t1 > 64 )		// [A-Z]
			a = (t1 - 65);
		else if( t1 > 47 )		// [0-9]
			a = (t1 + 4);
		else if( t1 == 43 )
			a = 62;
		else				// src[0] == '/'
			a = 63;     


		if( t2 > 96 )		// [a-z]
			b = (t2 - 71);
		else if( t2 > 64 )		// [A-Z]
			b = (t2 - 65);
		else if( t2 > 47 )		// [0-9]
			b = (t2 + 4);
		else if( t2 == 43 )
			b = 62;
		else				// src[0] == '/'
			b = 63;     

		*temp++ = ( a << 2) | ( b >> 4);

		if (t3 == 61 || temp >= end)
			return (int)(temp-dest);;

		if( t3 > 96 )		// [a-z]
			a = (t3 - 71);
		else if( t3 > 64 )		// [A-Z]
			a = (t3 - 65);
		else if( t3 > 47 )		// [0-9]
			a = (t3 + 4);
		else if( t3 == 43 )
			a = 62;
		else				// src[0] == '/'
			a = 63;     

		*temp++ = ( b << 4) | ( a >> 2);

		if (t4 == 61 || temp >= end)
			return (int)(temp-dest);;

		if( t4 > 96 )		// [a-z]
			b = (t4 - 71);
		else if( t4 > 64 )		// [A-Z]
			b = (t4 - 65);
		else if( t4 > 47 )		// [0-9]
			b = (t4 + 4);
		else if( t4 == 43 )
			b = 62;
		else				// src[0] == '/'
			b = 63;    

		*temp++ = ( a << 6) | ( b );
	}
}
