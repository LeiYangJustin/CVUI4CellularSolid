#pragma once

// ?
#pragma warning(disable: 4251)
#pragma warning(disable: 4244)
#pragma warning(disable: 4267)

#ifdef OPT_EXPORTS
#define OPT_API extern "C" __declspec(dllexport)
#define OPT_CLASS __declspec(dllexport)
#define OPT_TEMPLATE __declspec(dllexport)
#else
#define OPT_API extern "C" __declspec(dllimport)
#define OPT_CLASS __declspec(dllimport)
#define OPT_TEMPLATE
#endif