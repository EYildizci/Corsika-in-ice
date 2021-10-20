/****************************************************************************
 *                                                                          *
 *  Copyright and any other appropriate legal protection of these           *
 *  computer programs and associated documentation are reserved in          *
 *  all countries of the world.                                             *
 *                                                                          *
 *  These programs or documentation may not be reproduced by any method     *
 *  without prior written consent of Karlsruhe Institute of Technology (KIT)*
 *  ot its delegate. Commercial and military use are explicitly forbidden.  *
 *                                                                          *
 *  The Karlsruhe Institute of Technology welcomes comments concerning the  *
 *  COAST code but undertakes no obligation for maintenance of the programs,*
 *  nor responsibility for their correctness, and accepts no liability      *
 *  whatsoever resulting from the use of its programs.                      *
 *                                                                          *
 ****************************************************************************/

#ifndef __INCLUDE_TPLOTTER_H__
#define __INCLUDE_TPLOTTER_H__

#include <string>
#include <map>
#include "threevector.h"

class TCorsika;
class GroundArea;
class ScenarioParams;

namespace crs {
  class CParticle;
  class CInteraction;
  class MEventHeader;
  class MEventEnd;
  class MRunHeader;
};

class TPlotter {

public:
  TPlotter ();
  ~TPlotter ();

  void Welcome() const;
  void InitializeRadioSimulation();
  void SetDirectory (std::string dir) {fDirectory=dir;}
  void SetFileName (std::string file) {fFileName=file;}
  void SetThinning (bool thinn) {fThinning=thinn;}
  void SetSlant (bool slant) {fSlant=slant;}
  void SetCurved (bool flag) {fCurved=flag;}
  void SetStackInput(bool flag) {fStackInput=flag;}
  void SetPreshower(bool flag) {fPreshower=flag;}
    
  void SetRunHeader (const crs::MRunHeader &header);
  void SetShowerHeader (const crs::MEventHeader &header);
  void SetShowerTrailer (const crs::MEventEnd &trailer);
    
  void SetEvent (int no) {fEventNo=no;}
  void SetRun (int no) {fRunNo=no;}
  void SetShowerZenith (float zenith);
  void SetShowerAzimuth (float azimuth);

  void SetFirstInteraction(const double x, const double y, const double z, 
			   const double X, const double t);
  void SetShowerAxis(const double zenith, const double azimuth);
  
  bool IsThinned() const {return fThinning;}
    
  void Init();
  void Write();

  void AddTrack (const crs::CParticle &pre, const crs::CParticle &post);
  void AddInteraction (const crs::CInteraction& interaction);

private:
  void Clear();

  double GetEffectiveRefractiveIndexBetween(const ThreeVector& p1, const ThreeVector& p2) const;
  double GetHeightAtPosition(const ThreeVector& p1) const;
  
  void Rotate(double x, double y, double z,
	      double &sx, double &sy, double &sz,
	      int inverse) const;
  
  int fVerbosityLevel;

  std::string fDirectory;
  std::string fFileName;
  
  bool fThinning;
  bool fSlant;
  bool fCurved;
  bool fStackInput;
  bool fPreshower;
  
  // to figure out the primary track
  bool  fPrimaryTrack;
  bool  fFirstInteraction;
  
  // point of the first interaction
  double fFirstInteractionX;
  double fFirstInteractionY;
  double fFirstInteractionZ;
  double fFirstInteractionTime;
  double fFirstInteractionDist;
  double fFirstInteractionSlantDepth; // slant depth

  // intermediate storage of magnetic field data
  double fBx;
  double fBz;

  // catesian shower axis
  double fAxisX;
  double fAxisY;
  double fAxisZ;
  
  // for rotations
  double fCosZenith;
  double fSinZenith;
  double fCosAzimuth;
  double fSinAzimuth;
  
  // information from CORSIKA
  double fHeightFirstInt; 
  float fXmax;
  float fSkimmingAltitude;
  bool fSkimming;
  int fEventNo;
  int fRunNo;
  double fObservationLevel;
  double fMaxShowerRange;
  float fZenith;
  float fAzimuth;
  float fShowerEnergy;
  int fShowerPrimary;
  
  struct ParticleDef {
      std::string name;
      int color;
  };
  
  std::map<int, ParticleDef> fParticles;
  
  // to keep track of versions
  std::string fCorsikaVersion;
  std::string fCoastVersion;
  std::string fHistVersion;

  TCorsika* fCORSIKA;

  // for directRadio
  GroundArea* fGroundArea;
  double fCoreHitTime;
  ThreeVector fRadioCore;
  double fGroundLevelRefractivity;
  ScenarioParams* itsScenarioParams;
};



#endif
