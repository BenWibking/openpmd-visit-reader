// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ****************************************************************************
//  File: openpmdPluginInfo.C
// ****************************************************************************

#include <openpmdPluginInfo.h>

#include <visit-config.h>
VISIT_PLUGIN_VERSION(openpmd,DBP_EXPORT)

VISIT_DATABASE_PLUGIN_ENTRY(openpmd,General)

// ****************************************************************************
//  Method: openpmdGeneralPluginInfo::GetName
//
//  Purpose:
//    Return the name of the database plugin.
//
//  Returns:    A pointer to the name of the database plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
openpmdGeneralPluginInfo::GetName() const
{
    return "openpmd-api";
}

// ****************************************************************************
//  Method: openpmdGeneralPluginInfo::GetVersion
//
//  Purpose:
//    Return the version of the database plugin.
//
//  Returns:    A pointer to the version of the database plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
openpmdGeneralPluginInfo::GetVersion() const
{
    return "";
}

// ****************************************************************************
//  Method: openpmdGeneralPluginInfo::GetID
//
//  Purpose:
//    Return the id of the database plugin.
//
//  Returns:    A pointer to the id of the database plugin.
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

const char *
openpmdGeneralPluginInfo::GetID() const
{
    return "openpmd_";
}
// ****************************************************************************
//  Method: openpmdGeneralPluginInfo::EnabledByDefault
//
//  Purpose:
//    Return true if this plugin should be enabled by default; false otherwise.
//
//  Returns:    true/false
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

bool
openpmdGeneralPluginInfo::EnabledByDefault() const
{
    return true;
}
// ****************************************************************************
//  Method: openpmdGeneralPluginInfo::HasWriter
//
//  Purpose:
//    Return true if this plugin has a database writer.
//
//  Returns:    true/false
//
//  Programmer: generated by xml2info
//  Creation:   omitted
//
// ****************************************************************************

bool
openpmdGeneralPluginInfo::HasWriter() const
{
    return false;
}
// ****************************************************************************
//  Method:  openpmdGeneralPluginInfo::GetDefaultFilePatterns
//
//  Purpose:
//    Returns the default patterns for a openpmd database.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
std::vector<std::string>
openpmdGeneralPluginInfo::GetDefaultFilePatterns() const
{
    std::vector<std::string> defaultPatterns;
    defaultPatterns.push_back("*.pmd");

    return defaultPatterns;
}

// ****************************************************************************
//  Method:  openpmdGeneralPluginInfo::AreDefaultFilePatternsStrict
//
//  Purpose:
//    Returns if the file patterns for a openpmd database are
//    intended to be interpreted strictly by default.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
bool
openpmdGeneralPluginInfo::AreDefaultFilePatternsStrict() const
{
    return false;
}

// ****************************************************************************
//  Method:  openpmdGeneralPluginInfo::OpensWholeDirectory
//
//  Purpose:
//    Returns if the openpmd plugin opens a whole directory name
//    instead of a single file.
//
//  Programmer:  generated by xml2info
//  Creation:    omitted
//
// ****************************************************************************
bool
openpmdGeneralPluginInfo::OpensWholeDirectory() const
{
    return false;
}
