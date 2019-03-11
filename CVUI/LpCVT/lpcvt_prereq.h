#pragma once

// ?
#pragma warning(disable: 4251)
#pragma warning(disable: 4244)
#pragma warning(disable: 4267)

#ifdef LPCVT_EXPORTS
#define LPCVT_API extern "C" __declspec(dllexport)
#define LPCVT_CLASS __declspec(dllexport)
#define LPCVT_TEMPLATE __declspec(dllexport)
#else
#define LPCVT_API extern "C" __declspec(dllimport)
#define LPCVT_CLASS __declspec(dllimport)
#define LPCVT_TEMPLATE
#endif