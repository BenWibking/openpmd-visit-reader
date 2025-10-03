#include "avtopenpmdFileFormat.h"

#include <avtDatabaseMetaData.h>
#include <DebugStream.h>
#include <InvalidFilesException.h>
#include <VisItException.h>

#include <iostream>
#include <memory>

class accesstester : public avtopenpmdFileFormat {
public:
  using avtopenpmdFileFormat::avtopenpmdFileFormat;
  using avtopenpmdFileFormat::PopulateDatabaseMetaData;
};

int main(int argc, char **argv) {
  const char *path = argc > 1 ? argv[1] : "/Users/benwibking/openpmd-api-plugin/example_data/bp4_3d.pmd";

  try {
    accesstester reader(path);
    auto md = std::make_unique<avtDatabaseMetaData>();
    reader.PopulateDatabaseMetaData(md.get(), 0);
    std::cout << "Successfully populated metadata for " << path << "\n";
  } catch (const VisItException &ex) {
    std::cerr << "VisItException: " << ex.Message() << "\n";
    return 1;
  } catch (const std::exception &ex) {
    std::cerr << "std::exception: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}
