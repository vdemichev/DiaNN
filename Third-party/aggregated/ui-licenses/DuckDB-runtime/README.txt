DuckDB compiler-runtime notices

These notices cover support code incorporated into the supplied DuckDB library.
The compilers and their development tools are not included in the distribution.

Linux: GCC 9.3.0 runtime/header code uses the GCC Runtime Library Exception,
reproduced with the applicable copyrights in GCC-9.txt. The base GPL version 3
text is supplied in ../system-libraries/license-texts/GPL-3.txt.
The separately installed libstdc++ and libgcc libraries remain OS components;
their notices are under ../system-libraries/.

Windows: LLVM 22.1.8 libc++, libc++abi, libunwind and compiler-rt use the Apache
License 2.0 with LLVM Exceptions and the retained legacy grants in their four
component notice files. MinGW-w64 runtime support at source revision
c28e9555bb8800c53449f42a465ad9a5676fce88 uses the public-domain, permissive and
Zope Public License grants reproduced in MinGW-runtime-selected.txt.

Windows system services and the Universal C Runtime are supplied by Windows;
no Microsoft DLL is redistributed here. See ../system-libraries/windows/.
