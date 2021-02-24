#include <iostream>

/// \file
/// Implementation of the assert function for MiniDNN
#pragma once

inline void __M_Assert(const char* expr_str, bool expr, const char* file,
                       int line, const char* msg)
{
    if (!expr)
    {
        std::cerr << "Assert failed:\t" << msg << "\n"
                  << "Expected:\t" << expr_str << "\n"
                  << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}

#ifndef NDEBUG
#define M_Assert(Expr, Msg) \
__M_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
#define M_Assert(Expr, Msg) ;
#endif

