// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include "LHAPDF/Info.h"

namespace LHAPDF {


  /// Class for PDF set metadata and manipulation
  class Config : public Info {
  public:

    /// @name Fetching/creation
    //@{

    /// Get the global configuration object
    ///
    /// The global config is populated by reading from lhapdf.conf if it is
    /// found in the search paths. It is a singleton, hence the 'get' accessor
    /// rather than a constructor.
    ///
    /// @note The LHAPDF system is responsible for deletion of the returned
    /// object. Do NOT delete it yourself!
    static Config& get() {
      static Config _cfg; //< Could we use the Info(path) constructor for automatic init-once behaviour?
      // Test for emptiness and only initialise *once*:
      if (_cfg._metadict.empty()) {
        string confpath = findFile("lhapdf.conf").string();
        if (!confpath.empty()) _cfg.load(confpath);
      }
      return _cfg;
    }

    //@}

  private:

    /// Hide the default constructor
    Config() { }

    //@}


  public:

    /// @name Global config metadata specialisations
    //@{

    /// @brief All available PDF set names
    ///
    /// @note Taken from scanning the directories in LHAPDF_DATA_PATH for
    /// viable PDF sets.
    // std::vector<std::string> availablePDFSets() const {
    //   /// @todo Write!
    // }

    //@}

  };


  /// @name Convenient verbosity control
  //@{

  /// Convenient way to get the current verbosity level
  ///
  /// @note Verbosity, like any other flag, can also be set at lower levels. But who does that, really?!?
  inline int verbosity() {
    return Config::get().get_entry_as<int>("Verbosity", 1);
  }

  /// Convenient way to set the verbosity level
  ///
  /// @note Verbosity, like any other flag, can also be set at lower levels. But who does that, really?!?
  inline void setVerbosity(int v) {
    Config::get().set_entry("Verbosity", v);
  }

  //@}


}
