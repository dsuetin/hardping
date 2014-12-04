/* include/LHAPDF/Version.h.  Generated from Version.h.in by configure.  */
// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include <string>

/* "LHAPDF version string" */
#define LHAPDF_VERSION "6.0.4"

/* "LHAPDF version as an int" */
#define LHAPDF_VERSION_CODE 60004

// Separate int-valued macro for conditional compilation. Doesn't exist in LHAPDF5.
// Use like "#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6 ..."
#define LHAPDF_MAJOR_VERSION 6

namespace LHAPDF {


  /// Get the LHAPDF library version code (as a string)
  inline std::string version() {
    return LHAPDF_VERSION;
  }


}
