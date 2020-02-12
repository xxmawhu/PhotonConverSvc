#ifndef PhotonConverSvc_PhotonConverSvc_H
#define PhotonConverSvc_PhotonConverSvc_H

#include "GaudiKernel/Service.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "McTruth/McParticle.h"
#include "McTruth/MdcMcHit.h"
#include "MdcRecEvent/RecMdcHit.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "PhotonConverSvc/IPhotonConverSvc.h"
#include "EventNavigator/EventNavigator.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "DatabaseSvc/IDatabaseSvc.h"
#include "BesDChain/CDDecayList.h"
#include "BesDChain/CDPhotonList.h"
#include "BesDChain/CDPi0List.h"
#include "BesStdSelector/Selector.h"
#include <mysql.h>
#include "McDecayModeSvc/McDecayModeSvc.h"

#include "HadronInfo/AvailableInfo.h"
#include "TupleSvc/DecayTree.h"
#include <map>
#include <vector>
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
using std::vector;
using Event::McParticle;
using Event::McParticleCol;
using Event::MdcMcHitCol;

class PhotonConverSvc : public Service,
                        virtual public IPhotonConverSvc,
                        public AvailableInfo {
   public:
    PhotonConverSvc(const std::string& name, ISvcLocator* svcLoc);
    virtual ~PhotonConverSvc();
    virtual StatusCode initialize();
    virtual StatusCode finalize();
    virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvIF);

    virtual void GetInfoI(const std::string&, int&);
    virtual void GetInfoD(const std::string&, double&);
    virtual void GetInfoVd(const std::string&, std::vector<double>&);
    void Feed(const CDCandidate& sig);
    void SetDecayTree(DecayTree);
    void SetEcm(double Ecm);

   private:
    EvtRecTrack* m_tracks[2];
    HepLorentzVector m_p4Beam;
    DecayTree m_decayTree;
    double m_Ecm;
    double m_Rxy, m_Rx, m_Ry;
    double m_mrec;
    int m_Ngamam;
    vector<double> m_EEGList;
    CDPhotonList m_PhotonList;
    bool GetParList();
    IDataProviderSvc* eventSvc_;
    IDatabaseSvc* m_dbsvc;
    // McDecayModeSvc* m_MCDecayModeSvc;
    mutable EventNavigator* m_navigator;
    bool m_usecbE;
    int m_runID, m_eventID, m_status;
    void ReadDb(int run, double&);
    void AnaBeamStatus();
    void UpdateAvialInfo();
};
#endif
