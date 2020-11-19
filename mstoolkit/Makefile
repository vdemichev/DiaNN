#Set your paths here.
ZLIB_PATH = ./src/zlib-1.2.8
MZPARSER_PATH = ./src/mzParser
EXPAT_PATH = ./src/expat-2.2.0
SQLITE_PATH = ./src/sqlite-3.7.7.1
MST_PATH = ./src/MSToolkit

HEADER_PATH = ./include

MZPARSER = mzp.MSNumpress.o mzp.mzp_base64.o mzp.BasicSpectrum.o mzp.mzParser.o mzp.RAMPface.o mzp.saxhandler.o mzp.saxmzmlhandler.o \
	mzp.saxmzxmlhandler.o mzp.Czran.o mzp.mz5handler.o mzp.mzpMz5Config.o mzp.mzpMz5Structs.o mzp.BasicChromatogram.o mzp.PWIZface.o
MZPARSERLITE = mzp.MSNumpress.o mzp.mzp_base64_lite.o mzp.BasicSpectrum_lite.o mzp.mzParser_lite.o mzp.RAMPface_lite.o mzp.saxhandler_lite.o mzp.saxmzmlhandler_lite.o \
  mzp.saxmzxmlhandler_lite.o mzp.Czran_lite.o mzp.mz5handler_lite.o mzp.mzpMz5Config_lite.o mzp.mzpMz5Structs_lite.o mzp.BasicChromatogram_lite.o mzp.PWIZface_lite.o
EXPAT = xmlparse.o xmlrole.o xmltok.o
ZLIB = adler32.o compress.o crc32.o deflate.o inffast.o inflate.o infback.o inftrees.o trees.o uncompr.o zutil.o
MSTOOLKIT = Spectrum.o MSObject.o mzMLWriter.o pepXMLWriter.o
READER = MSReader.o
READERLITE = MSReaderLite.o
SQLITE = sqlite3.o 

CC = g++
GCC = gcc
NOSQLITE = -D_NOSQLITE

CFLAGS = -O3 -static -I. -I$(HEADER_PATH) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DHAVE_EXPAT_CONFIG_H
LIBS = -lm -ldl

CMD = $(CC) -o ../diann-linux ../src/diann.cpp -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++17 -lstdc++fs -fpermissive -O3 -mfpmath=sse -msse2 -march=core2 -w -L. -lmstoolkitlite -I. -I$(HEADER_PATH) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DHAVE_EXPAT_CONFIG_H
# MacOS 10.15.7
ifneq ($(OS),Windows_NT)
 ifeq ($(shell uname),Darwin)
  CMD = $(CC) -o ../diann-linux ../src/diann.cpp -lpthread -std=c++14 -fpermissive -O3 -mfpmath=sse -msse2 -march=core2 -w -L. -lmstoolkitlite -I. -I$(HEADER_PATH) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DHAVE_EXPAT_CONFIG_H
 endif
endif

all:  $(ZLIB) $(MZPARSERLITE) $(MSTOOLKIT) $(READERLITE) $(EXPAT)
	ar rcs libmstoolkitlite.a $(ZLIB) $(EXPAT) $(MZPARSERLITE) $(MSTOOLKIT) $(READERLITE)
	$(CMD)
clean:
	rm -f *.o libmstoolkitlite.a libmstoolkit.a msSingleScan diann-linux
	rm -f ../diann-linux

# zLib objects

adler32.o : $(ZLIB_PATH)/adler32.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/adler32.c -c

compress.o : $(ZLIB_PATH)/compress.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/compress.c -c

crc32.o : $(ZLIB_PATH)/crc32.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/crc32.c -c

deflate.o : $(ZLIB_PATH)/deflate.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/deflate.c -c

inffast.o : $(ZLIB_PATH)/inffast.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/inffast.c -c

inflate.o : $(ZLIB_PATH)/inflate.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/inflate.c -c

infback.o : $(ZLIB_PATH)/infback.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/infback.c -c

inftrees.o : $(ZLIB_PATH)/inftrees.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/inftrees.c -c

trees.o : $(ZLIB_PATH)/trees.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/trees.c -c

uncompr.o : $(ZLIB_PATH)/uncompr.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/uncompr.c -c

zutil.o : $(ZLIB_PATH)/zutil.c
	$(GCC) $(CFLAGS) $(ZLIB_PATH)/zutil.c -c



#mzParser objects
mzp.%.o : $(MZPARSER_PATH)/%.cpp
	$(CC) $(CFLAGS) $(MZ5) $< -c -o $@

#mzParserLite objects
mzp.%_lite.o : $(MZPARSER_PATH)/%.cpp
	$(CC) $(CFLAGS) $< -c -o $@


#expat objects
xmlparse.o : $(EXPAT_PATH)/xmlparse.c
	$(GCC) $(CFLAGS) $(EXPAT_PATH)/xmlparse.c -c
xmlrole.o : $(EXPAT_PATH)/xmlrole.c
	$(GCC) $(CFLAGS) $(EXPAT_PATH)/xmlrole.c -c
xmltok.o : $(EXPAT_PATH)/xmltok.c
	$(GCC) $(CFLAGS) $(EXPAT_PATH)/xmltok.c -c



#SQLite object
sqlite3.o : $(SQLITE_PATH)/sqlite3.c
	$(GCC) $(CFLAGS) $(SQLITE_PATH)/sqlite3.c -c




#MSToolkit objects

Spectrum.o : $(MST_PATH)/Spectrum.cpp
	$(CC) $(CFLAGS) $(MST_PATH)/Spectrum.cpp -c

MSReader.o : $(MST_PATH)/MSReader.cpp
	$(CC) $(CFLAGS) $(MZ5) $(MST_PATH)/MSReader.cpp -c

MSReaderLite.o : $(MST_PATH)/MSReader.cpp
	$(CC) $(CFLAGS) $(NOSQLITE) $(MST_PATH)/MSReader.cpp -c -o MSReaderLite.o

MSObject.o : $(MST_PATH)/MSObject.cpp
	$(CC) $(CFLAGS) $(MST_PATH)/MSObject.cpp -c

mzMLWriter.o : $(MST_PATH)/mzMLWriter.cpp
	$(CC) $(CFLAGS) $(MST_PATH)/mzMLWriter.cpp -c
	
pepXMLWriter.o : $(MST_PATH)/pepXMLWriter.cpp
	$(CC) $(CFLAGS) $(MST_PATH)/pepXMLWriter.cpp -c






