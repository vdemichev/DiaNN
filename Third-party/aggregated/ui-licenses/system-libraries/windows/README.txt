Microsoft Windows system libraries and Universal C Runtime

Copyright Microsoft Corporation. All rights reserved.

Windows supplies the native windowing, graphics, text, input, accessibility,
security, file-system and process services used by DIA-NN. These include the
Win32, COM, Shell, GDI, DirectWrite, Direct3D, DXGI and DirectComposition system
libraries and their operating-system dependencies. They are governed by the
license terms accompanying the installed copy of Microsoft Windows.

Microsoft Windows license information:
https://www.microsoft.com/useterms

DuckDB also uses the Windows Universal C Runtime through its api-ms-win-crt-*
interfaces and the Windows Restart Manager. These are supplied by Windows 10
and later under the applicable Microsoft license terms. The supplied DuckDB
library does not require a separate Microsoft Visual C++ runtime installation.

The portable DIA-NN UI package contains no Microsoft system or runtime DLLs.
This notice does not replace their existing licenses or grant redistribution
rights for those components.
