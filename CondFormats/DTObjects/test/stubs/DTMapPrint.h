
/*----------------------------------------------------------------------

Toy EDAnalyzer for testing purposes only.

----------------------------------------------------------------------*/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/DTObjects/interface/DTReadOutMapping.h"
#include "CondFormats/DataRecord/interface/DTReadOutMappingRcd.h"

namespace edmtest {
  class DTMapPrint : public edm::EDAnalyzer {
  public:
    explicit DTMapPrint(edm::ParameterSet const& p);
    explicit DTMapPrint(int i);
    virtual ~DTMapPrint();
    virtual void analyze(const edm::Event& e, const edm::EventSetup& c);

  private:
    edm::ESGetToken<DTReadOutMapping, DTReadOutMappingRcd> es_token;
  };
}  // namespace edmtest
