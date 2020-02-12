/* ====================================================
#   Copyright (C)2020 All rights reserved.
#
#   Author        : Xin-Xin MA
#   Email         : xxmawhu@163.com
#   File Name     : IPhotonConverSvc.h
#   Create Time   : 2020-02-08 15:07
#   Last Modified : 2020-02-08 15:07
#   Describe      :
#
# ====================================================*/
#ifndef PhotonConverSvc__IPhotonConverSvc_H
#define PhotonConverSvc__IPhotonConverSvc_H
#include "GaudiKernel/IService.h"
#include "McTruth/McParticle.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "CLHEP/Vector/LorentzVector.h"

/* Decaration of the interface ID */
static const InterfaceID IID_IPhotonConverSvc("IPhotonConverSvc", 1, 0);

class IPhotonConverSvc : virtual public IService {
   public:
    virtual ~IPhotonConverSvc() {}
    static const InterfaceID& interfaceID() { return IID_IPhotonConverSvc; }
};
#endif
