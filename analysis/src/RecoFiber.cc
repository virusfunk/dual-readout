#include "TVector3.h"

#include "RecoFiber.h"

#include <cmath>

RecoFiber::RecoFiber() {
  fCalibS = 0.;
  fCalibC = 0.;

  // [L] = mm, [t] = ns
  fSpeed = 189.5;
  fEffSpeedInv = 1./fSpeed - 1./300.;
  fDepthEM = 149.49;
  fAbsLen = 5677.;
  fCThres = 30.;
}

void RecoFiber::reconstruct(const DRsimInterface::DRsimSiPMData& sipm, RecoInterface::RecoTowerData& recoTower) {
  RecoInterface::RecoFiberData recoFiber(sipm);

  if (fSeg->IsCerenkov(recoFiber.fiberNum)) {
    recoFiber.n = cutXtalk(sipm);
    recoFiber.E = (float)recoFiber.n / fCalibC;
    recoFiber.Ecorr = recoFiber.E;
    recoFiber.t = setTmax(sipm);

    recoFiber.depth = setDepth(recoFiber.t, recoTower, sipm);

    addFjInputs(recoFiber);
  } else {
    recoFiber.E = (float)recoFiber.n / fCalibS;
    recoFiber.t = setTmax(sipm);

    recoFiber.Ecorr = recoFiber.E*std::exp((-recoFiber.depth+fDepthEM)/fAbsLen);
    addFjInputs(recoFiber);
  }

  fData = recoFiber;
  recoTower.fibers.push_back(fData);
}

float RecoFiber::setTmax(const DRsimInterface::DRsimSiPMData& sipm) {
  std::pair<DRsimInterface::hitRange,int> maxima = std::make_pair(std::make_pair(0.,0.),0);
  for (auto timeItr = sipm.timeStruct.begin(); timeItr != sipm.timeStruct.end(); ++timeItr) {
    auto timeObj = *timeItr;

    if (timeObj.second > maxima.second) maxima = timeObj;
  }

  return maxima.first.first;
}

float RecoFiber::setDepth(const float tmax, const RecoInterface::RecoTowerData& recoTower, const DRsimInterface::DRsimSiPMData& sipm) {
  auto position = fSeg->towerposition_from_cell_ID(sipm.SiPMnum); // cm
  TVector3 temp_vec; temp_vec.SetXYZ(position.x(), position.y(), position.z());
  double towerH = fSeg->towerheight_from_cell_ID(sipm.SiPMnum); // cm

  float distance = (float) ( (temp_vec.Mag() - towerH / 2.0) * 10.0 ); // mm

  float depth = ( distance/300. + 2000./fSpeed - tmax )/( fEffSpeedInv ); // 300mm/ns
  if (depth < 0.) return 0.;
  else if (depth > 2000.) return 2000.;
  else return depth;
}

int RecoFiber::cutXtalk(const DRsimInterface::DRsimSiPMData& sipm) {
  int sum = 0;
  for (auto timeItr = sipm.timeStruct.begin(); timeItr != sipm.timeStruct.end(); ++timeItr) {
    auto timeObj = *timeItr;

    if (timeObj.first.first < fCThres) sum += timeObj.second;
  }

  return sum;
}

void RecoFiber::addFjInputs(const RecoInterface::RecoFiberData& recoFiber) {
  auto global = fSeg->position(recoFiber.fiberNum);
  TVector3 vec(global.x(),global.y(),global.z());
  TVector3 p = recoFiber.E*vec.Unit();

  if (fSeg->IsCerenkov(recoFiber.fiberNum)) {
    fFjInputs_C.push_back( fastjet::PseudoJet(p.x(),p.y(),p.z(),recoFiber.E) );
    fFjInputs_depth.push_back(recoFiber.depth);
    fFjInputs_fibernum_C.push_back(recoFiber.fiberNum);
  } else {
    fFjInputs_S.push_back( fastjet::PseudoJet(p.x(),p.y(),p.z(),recoFiber.E) );

    TVector3 p_corr = recoFiber.Ecorr*vec.Unit();
    fFjInputs_Scorr.push_back( fastjet::PseudoJet(p_corr.x(),p_corr.y(),p_corr.z(),recoFiber.Ecorr) );

    fFjInputs_fibernum_S.push_back(recoFiber.fiberNum);
    fFjInputs_fibernum_Scorr.push_back(recoFiber.fiberNum);
  }
}

void RecoFiber::clear() {
  fFjInputs_S.clear();
  fFjInputs_Scorr.clear();
  fFjInputs_C.clear();
  fFjInputs_depth.clear();
  fFjInputs_fibernum_S.clear();
  fFjInputs_fibernum_Scorr.clear();
  fFjInputs_fibernum_C.clear();
}
