#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "EventModel/EventHeader.h"
#include "PhotonConverSvc/PhotonConverSvc.h"
#include <algorithm>

using CLHEP::HepLorentzVector;
using namespace Event;
using std::cout;
using std::endl;
StatusCode PhotonConverSvc::initialize() {
    MsgStream log(messageService(), name());
    log << MSG::INFO << "@initialize()" << endreq;

    StatusCode sc = Service::initialize();
    sc = serviceLocator()->service("EventDataSvc", eventSvc_, true);
    return sc;
}

StatusCode PhotonConverSvc::finalize() {
    MsgStream log(messageService(), name());
    log << MSG::INFO << "@initialize()" << endreq;
    StatusCode sc = Service::finalize();
    return sc;
}

StatusCode PhotonConverSvc::queryInterface(const InterfaceID& riid, void** ppvIF) {
    if (IPhotonConverSvc::interfaceID().versionMatch(riid)) {
        *ppvIF = dynamic_cast<IPhotonConverSvc*>(this);
    } else {
        return Service::queryInterface(riid, ppvIF);
    }
    addRef();
    // cout<<"PhotonConverSvc::Inf:queryInterface"<<endl;
    return StatusCode::SUCCESS;
}
