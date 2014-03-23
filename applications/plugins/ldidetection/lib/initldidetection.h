#ifndef INITldidetection_H
#define INITldidetection_H


#include <sofa/helper/system/config.h>

#ifdef SOFA_BUILD_LDIDETECTION
#define SOFA_LDIDETECTION_API SOFA_EXPORT_DYNAMIC_LIBRARY
#else
#define SOFA_LDIDETECTION_API  SOFA_IMPORT_DYNAMIC_LIBRARY
#endif

/** \mainpage
  This is the ldidetection plugin.
  */

#endif
