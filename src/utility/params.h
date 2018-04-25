#ifndef __UTILITY_PARAMS_H__
#define __UTILITY_PARAMS_H__

#include <utility>
#include <string>
#include <vector>

using std::pair;
using std::string;
using std::vector;

namespace Utility {
// Constants for spectrum in .sptxt file.
const string BINARYFILEOFFSET_SPTXT = "BinaryFileOffset";
const string NAME_SPTXT = "Name";
const string NUMPEAKS_SPTXT = "NumPeaks";
const string PRECURSORMZ_SPTXT = "PrecursorMZ";

// Constants for spectrum in .msp file.
const string NAME_MSP = "Name";
const string NUMPEAKS_MSP = "Num peaks";
const string PRECURSORMZ_MSP = "Parent";
// const string BINARYFILEOFFSET_MSP = "BinaryFileOffset";

// Constants for SpectraST search result in .pep.xml file.
const string LIB_FILE_OFFSET_XML = "lib_file_offset";
const string SPECTRUM_XML = "spectrum=\"";
const string VALUE_XML = "value";

// Constants for spectrum in .mgf file.
const string BEGIN_IONS_MGF = "BEGIN IONS";
const string CHARGE_MGF = "CHARGE";
const string END_IONS_MGF = "END IONS";
const string PEPMASS_MGF = "PEPMASS";
const string TITLE_MGF = "TITLE"; 

// Constants for raw peaks size for pre-allocating memory.
const int MAX_PEAK_SIZE = 10000;
}  // namespace Utility
#endif
