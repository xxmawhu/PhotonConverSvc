#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/PropertyMgr.h"
#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "GammaConv/GammaConv.h"
#include "PhotonConverSvc/PhotonConverSvc.h"
#include <algorithm>
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

using CLHEP::HepLorentzVector;
using namespace Event;
using namespace BesStdSelector;
using std::cout;
using std::endl;
PhotonConverSvc::PhotonConverSvc(const std::string& name, ISvcLocator* svcLoc)
    : Service(name, svcLoc),
      m_Rxy(0),
      m_Rx(0),
      m_Ry(0),
      m_mrec(0),
      m_runID(0),
      m_eventID(0),
      m_Ngamam(0) {
    UpdateAvialInfo();
    m_tracks[0] = NULL;
    m_tracks[1] = NULL;
    // calibrated beam Energy
    SetName("PhotonConv");
    declareProperty("ReadEcm", m_readEcm = true);
    declareProperty("SetEcm", m_Ecm = 4.260);
    declareProperty("UseCbE", m_usecbE = false);
}

PhotonConverSvc::~PhotonConverSvc() {}

void PhotonConverSvc::UpdateAvialInfo() {
    // AvailableInfo::Clear();
    AvailableInfo::Add("Rx", "double");
    AvailableInfo::Add("Ry", "double");
    AvailableInfo::Add("Rxy", "double");
    AvailableInfo::Add("Ngamma", "int");
    AvailableInfo::Add("mRec", "double");
    AvailableInfo::Add("Angle", "double");
    AvailableInfo::Add("mEEL", "double", "Ngamma");
    AvailableInfo::Add("Mee", "double");
}
void PhotonConverSvc::GetInfoI(const std::string& info_name, int& targe) {
    // cout << "Info in PhotonConverSvc::GetInfoI: "
    //    << "info_name = " << info_name << endl;
    if (info_name == "Ngamma") {
        targe = m_EEGList.size();
    }
}
void PhotonConverSvc::GetInfoD(const std::string& info_name, double& targe) {
    if (info_name == "Rx") {
        targe = m_Rx;
    } else if (info_name == "Ry") {
        targe = m_Ry;
    } else if (info_name == "Rxy") {
        targe = m_Rxy;
    } else if (info_name == "mRec") {
        targe = m_mrec;
    } else if (info_name == "Angle") {
        targe = m_angle/3.1415926*180;
    } else if (info_name == "Mee") {
        targe = m_Mee;
    }
}
void PhotonConverSvc::GetInfoVd(const std::string& info_name,
                                std::vector<double>& targe) {
    if (info_name == "mEEL") {
        targe = m_EEGList;
    }
}

void PhotonConverSvc::SetDecayTree(DecayTree decayTree) {
    m_decayTree = decayTree;
}
void PhotonConverSvc::SetEcm(double Ecm) { m_Ecm = Ecm; }

void PhotonConverSvc::ReadDb(int run, double& Ecm) {
    Gaudi::svcLocator()->service("DatabaseSvc", m_dbsvc, true);
    // calibrated beam Energy
    if (m_usecbE) {
        char stmt1[1024];
        snprintf(stmt1, 1024,
                 "select beam_energy, px, py, pz "
                 "from RunParams where run_number = %d",
                 run);
        DatabaseRecordVector res;
        int row_no = m_dbsvc->query("offlinedb", stmt1, res);
        if (row_no == 0) {
            m_status = 401;
            return;
        }
        DatabaseRecord* records = res[0];
        Ecm = 2 * records->GetDouble("beam_energy");
        /// cout << "calibrated beam_energy = " <<
        /// records->GetDouble("beam_energy")
        ///    << endl;
    } else {
        // use online beam Energy
        char stmt1[1024];
        snprintf(stmt1, 1024,
                 "select BER_PRB, BPR_PRB "
                 "from RunParams where run_number = %d",
                 run);
        DatabaseRecordVector res;
        int row_no = m_dbsvc->query("run", stmt1, res);
        if (row_no == 0) {
            return;
        }
        DatabaseRecord* records = res[0];
        double E_E, E_P;
        E_E = records->GetDouble("BER_PRB");
        E_P = records->GetDouble("BPR_PRB");
        Ecm = E_E + E_P;
        /// cout << "beam_energy = " << Ecm;
        /// cout  << endl;
    }
    return;
}
void PhotonConverSvc::AnaBeamStatus() {
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc_,
                                                 "/Event/EventHeader");
    if (m_runID == eventHeader->runNumber() &&
        m_eventID == eventHeader->eventNumber()) {
        //  cout << "Info in PhotonConverSvc::AnaBeamStatus: "
        //      << "#run = " << m_run << ", #id = " << m_event << endl;
        return;
    } else {
        m_runID = eventHeader->runNumber();
        m_eventID = eventHeader->eventNumber();
        //  cout << "Info in PhotonConverSvc::AnaBeamStatus: "
        //      << "#run = " << m_run << ", #id = " << m_event << endl;
    }
    this->ReadDb(abs(m_runID), m_Ecm);
}

void PhotonConverSvc::GetParList() {
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc_,
                                                 "/Event/EventHeader");
    if (m_runID == eventHeader->runNumber() &&
        m_eventID == eventHeader->eventNumber()) {
        return;
    } else {
        m_runID = eventHeader->runNumber();
        m_eventID = eventHeader->eventNumber();
    }
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc_,
                                          "/Event/EvtRec/EvtRecEvent");
    SmartDataPtr<EvtRecTrackCol> evtRecTrackCol(eventSvc_,
                                                "/Event/EvtRec/EvtRecTrackCol");
    if (soloPhotonSelector.vetoPi0()) {
        SmartDataPtr<EvtRecPi0Col> recPi0Col(eventSvc_,
                                             "/Event/EvtRec/EvtRecPi0Col");
        CDPi0List Pi0List(pi0Selector);
        dc_fill(Pi0List, recPi0Col->begin(), recPi0Col->end());
        /// cout << "Info in <PhotonConverSvc::GetParList>: ";
        /// cout << "Pi0List.size = " << Pi0List.size() << endl;

        // fill into  the vector
        vector<const EvtRecPi0*> _pi0s;
        for (CDPi0List::iterator itr = Pi0List.particle_begin();
             itr != Pi0List.particle_end(); ++itr) {
            const EvtRecPi0* navPi0 = (*itr).particle().navPi0();
            _pi0s.push_back(navPi0);
        }
        BesStdSelector::soloPhotonSelector.setPi0s(_pi0s);
    }

    int nCharged = evtRecEvent->totalCharged();
    EvtRecTrackIterator neu_begin = evtRecTrackCol->begin() + nCharged;
    EvtRecTrackIterator neu_end = evtRecTrackCol->end();
    ///  cout << "Info in <PhotonConverSvc::GetParList>: ";
    ///  cout << "#showers = " << neu_end - neu_begin << endl;
    m_PhotonList = CDPhotonList(neu_begin, neu_end, soloPhotonSelector);
    /// cout << "Info in <PhotonConverSvc::GetParList>: ";
    /// cout << "#solo photon = " << m_PhotonList.size() << endl;
    return;
}

void PhotonConverSvc::Feed(const CDCandidate& sig) {
    m_EEGList.clear();
    HepLorentzVector p4Ep, p4Em;
    for (int i = 0; i < m_decayTree.size(); ++i) {
        /// cout << "PID[" << i << "] = " << m_decayTree.PID(i) << endl;
        if (m_decayTree.PID(i) == 11) {
            const CDCandidate& trk = sig.decay().child(i);
            m_tracks[0] =
                const_cast<EvtRecTrack*>(trk.finalChildren().first[0]);
            p4Ep = trk.p4();
            continue;
        }
        if (m_decayTree.PID(i) == -11) {
            const CDCandidate& trk = sig.decay().child(i);
            p4Em = trk.p4();
            m_tracks[1] =
                const_cast<EvtRecTrack*>(trk.finalChildren().first[0]);
            continue;
        }
    }
    if (m_tracks[0] == NULL || m_tracks[1] == NULL) {
        return;
    }
    ///  cout << "Info in <PhotonConverSvc::Feed>: ";
    ///  cout << "get two tracks successful!" << endl;
    HepVector helix1 = m_tracks[0]->mdcKalTrack()->getZHelixE();
    HepVector helix2 = m_tracks[1]->mdcKalTrack()->getZHelixE();
    HepPoint3D IP(0, 0, 0);
    GammaConv gammaConv(helix1, helix2, IP);
    m_Rxy = gammaConv.getRxy();
    m_Rx = gammaConv.getRx();
    m_Ry = gammaConv.getRy();
    if (m_readEcm) {
        this->AnaBeamStatus();
    }
    /// cout << "Info in <PhotonConverSvc::Feed>: ";
    /// cout << "Ecm = " << m_Ecm << endl;
    m_angle = p4Em.angle(p4Ep.vect());
    m_Mee = (p4Ep + p4Em).m();
    m_p4Beam = HepLorentzVector(0, 0, 0, m_Ecm);
    p4Em.boost(-0.011, 0, 0);
    p4Ep.boost(-0.011, 0, 0);
    m_mrec = (m_p4Beam - p4Ep - p4Em).m();
    /// cout << "Info in <PhotonConverSvc::Feed>: ";
    /// cout << "mRec = " << m_mrec << endl;
    GetParList();
    m_Ngamam = m_PhotonList.size();
    for (CDPhotonList::iterator itr = m_PhotonList.particle_begin();
         itr != m_PhotonList.particle_end(); ++itr) {
        const CDCandidate& photon = (*itr).particle();
        double mass = (photon.p4() + p4Ep + p4Em).m();
        m_EEGList.push_back(mass);
    }
}
