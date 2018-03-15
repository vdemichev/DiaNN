/*
MSNumpress.hpp
johan.teleman@immun.lth.se

Copyright 2013 Johan Teleman

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
/*
==================== encodeInt ====================
Some of the encodings described below use a integer compression referred to simply as

encodeInt()

This encoding works on a 4 byte integer, by truncating initial zeros or ones.
If the initial (most significant) half byte is 0x0 or 0xf, the number of such
halfbytes starting from the most significant is stored in a halfbyte. This initial
count is then followed by the rest of the ints halfbytes, in little-endian order.
A count halfbyte c of

0 <= c <= 8 		is interpreted as an initial c 		0x0 halfbytes
9 <= c <= 15		is interpreted as an initial (c-8) 	0xf halfbytes

Ex:
int		c		rest
0 	=> 	0x8
-1	=>	0xf		0xf
23	=>	0x6 	0x7	0x1
*/

#ifndef _MSNUMPRESS_HPP_
#define _MSNUMPRESS_HPP_

#include <cstddef>
#include <vector>

// defines whether to throw an exception when a number cannot be encoded safely
// with the given parameters
#ifndef THROW_ON_OVERFLOW
#define THROW_ON_OVERFLOW true
#endif

namespace ms {
  namespace numpress {

    namespace MSNumpress {

      double optimalLinearFixedPoint(
        const double *data,
        size_t dataSize);

      /**
      * Encodes the doubles in data by first using a
      *   - lossy conversion to a 4 byte 5 decimal fixed point representation
      *   - storing the residuals from a linear prediction after first two values
      *   - encoding by encodeInt (see above)
      *
      * The resulting binary is maximally 8 + dataSize * 5 bytes, but much less if the
      * data is reasonably smooth on the first order.
      *
      * This encoding is suitable for typical m/z or retention time binary arrays.
      * On a test set, the encoding was empirically show to be accurate to at least 0.002 ppm.
      *
      * @data		pointer to array of double to be encoded (need memorycont. repr.)
      * @dataSize	number of doubles from *data to encode
      * @result		pointer to where resulting bytes should be stored
      * @fixedPoint	the scaling factor used for getting the fixed point repr.
      * 				This is stored in the binary and automatically extracted
      * 				on decoding.
      * @return		the number of encoded bytes
      */
      size_t encodeLinear(
        const double *data,
        const size_t dataSize,
        unsigned char *result,
        double fixedPoint);

      /**
      * Calls lower level encodeLinear while handling vector sizes appropriately
      *
      * @data		vector of doubles to be encoded
      * @result		vector of resulting bytes (will be resized to the number of bytes)
      */
      void encodeLinear(
        const std::vector<double> &data,
        std::vector<unsigned char> &result,
        double fixedPoint);

      /**
      * Decodes data encoded by encodeLinear.
      *
      * result vector guaranteed to be shorter or equal to (|data| - 8) * 2
      *
      * Note that this method may throw a const char* if it deems the input data to be corrupt, i.e.
      * that the last encoded int does not use the last byte in the data. In addition the last encoded
      * int need to use either the last halfbyte, or the second last followed by a 0x0 halfbyte.
      *
      * @data		pointer to array of bytes to be decoded (need memorycont. repr.)
      * @dataSize	number of bytes from *data to decode
      * @result		pointer to were resulting doubles should be stored
      * @return		the number of decoded doubles, or -1 if dataSize < 4 or 4 < dataSize < 8
      */
      size_t decodeLinear(
        const unsigned char *data,
        const size_t dataSize,
        double *result);

      /**
      * Calls lower level decodeLinear while handling vector sizes appropriately
      *
      * Note that this method may throw a const char* if it deems the input data to be corrupt, i.e..
      * that the last encoded int does not use the last byte in the data. In addition the last encoded
      * int need to use either the last halfbyte, or the second last followed by a 0x0 halfbyte.
      *
      * @data		vector of bytes to be decoded
      * @result		vector of resulting double (will be resized to the number of doubles)
      */
      void decodeLinear(
        const std::vector<unsigned char> &data,
        std::vector<double> &result);

      /////////////////////////////////////////////////////////////


      /**
      * Encodes the doubles in data by storing the residuals from a linear prediction after first two values.
      *
      * The resulting binary is the same size as the input data.
      *
      * This encoding is suitable for typical m/z or retention time binary arrays, and is
      * intended to be used before zlib compression to improve compression.
      *
      * @data		pointer to array of doubles to be encoded (need memorycont. repr.)
      * @dataSize	number of doubles from *data to encode
      * @result		pointer to were resulting bytes should be stored
      */
      size_t encodeSafe(
        const double *data,
        const size_t dataSize,
        unsigned char *result);


      /**
      * Decodes data encoded by encodeSafe.
      *
      * result vector is the same size as the input data.
      *
      * Might throw const char* is something goes wrong during decoding.
      *
      * @data		pointer to array of bytes to be decoded (need memorycont. repr.)
      * @dataSize	number of bytes from *data to decode
      * @result		pointer to were resulting doubles should be stored
      * @return		the number of decoded bytes
      */
      size_t decodeSafe(
        const unsigned char *data,
        const size_t dataSize,
        double *result);

      /////////////////////////////////////////////////////////////

      /**
      * Encodes ion counts by simply rounding to the nearest 4 byte integer,
      * and compressing each integer with encodeInt.
      *
      * The handleable range is therefore 0 -> 4294967294.
      * The resulting binary is maximally dataSize * 5 bytes, but much less if the
      * data is close to 0 on average.
      *
      * @data		pointer to array of double to be encoded (need memorycont. repr.)
      * @dataSize	number of doubles from *data to encode
      * @result		pointer to were resulting bytes should be stored
      * @return		the number of encoded bytes
      */
      size_t encodePic(
        const double *data,
        const size_t dataSize,
        unsigned char *result);

      /**
      * Calls lower level encodePic while handling vector sizes appropriately
      *
      * @data		vector of doubles to be encoded
      * @result		vector of resulting bytes (will be resized to the number of bytes)
      */
      void encodePic(
        const std::vector<double> &data,
        std::vector<unsigned char> &result);

      /**
      * Decodes data encoded by encodePic
      *
      * result vector guaranteed to be shorter of equal to |data| * 2
      *
      * Note that this method may throw a const char* if it deems the input data to be corrupt, i.e.
      * that the last encoded int does not use the last byte in the data. In addition the last encoded
      * int need to use either the last halfbyte, or the second last followed by a 0x0 halfbyte.
      *
      * @data		pointer to array of bytes to be decoded (need memorycont. repr.)
      * @dataSize	number of bytes from *data to decode
      * @result		pointer to were resulting doubles should be stored
      * @return		the number of decoded doubles
      */
      size_t decodePic(
        const unsigned char *data,
        const size_t dataSize,
        double *result);

      /**
      * Calls lower level decodePic while handling vector sizes appropriately
      *
      * Note that this method may throw a const char* if it deems the input data to be corrupt, i.e.
      * that the last encoded int does not use the last byte in the data. In addition the last encoded
      * int need to use either the last halfbyte, or the second last followed by a 0x0 halfbyte.
      *
      * @data		vector of bytes to be decoded
      * @result		vector of resulting double (will be resized to the number of doubles)
      */
      void decodePic(
        const std::vector<unsigned char> &data,
        std::vector<double> &result);

      /////////////////////////////////////////////////////////////


      double optimalSlofFixedPoint(
        const double *data,
        size_t dataSize);

      /**
      * Encodes ion counts by taking the natural logarithm, and storing a
      * fixed point representation of this. This is calculated as
      *
      * unsigned short fp = log(d + 1) * fixedPoint + 0.5
      *
      * the result vector is exactly |data| * 2 + 8 bytes long
      *
      * @data		pointer to array of double to be encoded (need memorycont. repr.)
      * @dataSize	number of doubles from *data to encode
      * @result		pointer to were resulting bytes should be stored
      * @return		the number of encoded bytes
      */
      size_t encodeSlof(
        const double *data,
        const size_t dataSize,
        unsigned char *result,
        double fixedPoint);

      /**
      * Calls lower level encodeSlof while handling vector sizes appropriately
      *
      * @data		vector of doubles to be encoded
      * @result		vector of resulting bytes (will be resized to the number of bytes)
      */
      void encodeSlof(
        const std::vector<double> &data,
        std::vector<unsigned char> &result,
        double fixedPoint);

      /**
      * Decodes data encoded by encodeSlof
      *
      * The return will include exactly (|data| - 8) / 2 doubles.
      *
      * Note that this method may throw a const char* if it deems the input data to be corrupt.
      *
      * @data		pointer to array of bytes to be decoded (need memorycont. repr.)
      * @dataSize	number of bytes from *data to decode
      * @result		pointer to were resulting doubles should be stored
      * @return		the number of decoded doubles
      */
      size_t decodeSlof(
        const unsigned char *data,
        const size_t dataSize,
        double *result);

      /**
      * Calls lower level decodeSlof while handling vector sizes appropriately
      *
      * Note that this method may throw a const char* if it deems the input data to be corrupt.
      *
      * @data		vector of bytes to be decoded
      * @result		vector of resulting double (will be resized to the number of doubles)
      */
      void decodeSlof(
        const std::vector<unsigned char> &data,
        std::vector<double> &result);

    } // namespace MSNumpress
  } // namespace msdata
} // namespace pwiz

#endif // _MSNUMPRESS_HPP_
