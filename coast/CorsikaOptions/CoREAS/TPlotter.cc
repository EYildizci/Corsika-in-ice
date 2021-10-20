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

#include <TPlotter.h>
#include <TCorsika.h>

#include <crs/CorsikaConsts.h>
#include <crs/CParticle.h>
#include <crs/CInteraction.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MRunHeader.h>
using namespace crs;

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include "groundelement.h"
#include "groundarea.h"
#include "scenarioparams.h"
using namespace std;



/* -------------------------------------------------------------------
   Define the particle names/types and CORSIKA ID here !!!!!!!!!!!!!!!!
*/
const unsigned int nParticles = 2;
const char* gParticle_name [nParticles] = {"electron", "positron"};
const int gParticle_id     [nParticles] = {3, 2};
// -------------------------------------------------------------------



/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   - first find the point of the first interaction, and define the range
     where the shower should get histogramed

   - rotate all particles in shower coordinate system, with (0,0,0) at 
     the first interaction

   - define planes in equal slant depth bins 

   - fill histograms    
 */

TPlotter::TPlotter() {
  
  //fVerbosityLevel = 100; // all includning TCorsika
  //fVerbosityLevel = 40; // all from TPlotter
  fVerbosityLevel = 0; // none
  
  fPrimaryTrack = true;
  fFirstInteraction = false;
  
  fSlant = false;
  fCurved = false;
  fStackInput = false;
  fPreshower = false;
  
  fCORSIKA = 0;
  
  Clear ();
  fDirectory = "";
    
  for (unsigned int j=0; j<nParticles; j++) {
    fParticles[gParticle_id[j]].name = gParticle_name[j];
  }

  // for direct radio
  fGroundArea = 0;
  fCoreHitTime = 0.0;
}



TPlotter::~TPlotter () {
  
  fParticles.clear();
  Clear();
  if (fGroundArea)
    delete fGroundArea;
}


void 
TPlotter::Welcome() 
   const 
{
  cout << "\n"
       << " *******************************************************\n"
       << " **                                               \n"
       << " ** You are using the CoREAS V" << setprecision(2) << ProgramVersion << " radio simulation code\n"
       << " **                                               \n"
       << " **  THIN/CURVED/SLANT/STACKIN/PRESHOWER: "
       << fThinning << "/" << fCurved << "/" << fSlant << "/" << fStackInput << "/" << fPreshower << "\n"
       << " **                                               \n"
       << " ** Please cite the following article:            \n"
       << " ** T. Huege, M. Ludwig, C.W. James, AIP Conf. Proc.\n"
       << " ** 1535, 128-132 (2013), doi:10.1063/1.4807534     \n"
       << " **                                               \n"
       << " *******************************************************\n\n"
       << endl;
}

void
TPlotter::InitializeRadioSimulation()
{
  // create SIM file names from run number
  ostringstream fname;
  fname << "SIM" << setw (6) << setfill('0') << fRunNo;
  itsScenarioParams = new ScenarioParams(fDirectory+fname.str(), fDirectory+fname.str());
  
  bool readok = itsScenarioParams->ReadFromFiles();

  if (readok)
  {
    // overwrite values with the correct ones from CORSIKA
    itsScenarioParams->itsShowerZenithAngle = fZenith*180./Pi;
    itsScenarioParams->itsShowerAzimuthAngle = fAzimuth*180./Pi;
    itsScenarioParams->itsPrimaryParticleType = fShowerPrimary;
    itsScenarioParams->itsPrimaryParticleEnergy = fShowerEnergy*1.e9;
    itsScenarioParams->itsMagneticFieldStrength = sqrt(fBx*fBx+fBz*fBz)/100.0;       // convert from microTesla to Gauss before writing out
    itsScenarioParams->itsMagneticFieldInclinationAngle = atan2(fBz,fBx)*180.0/M_PI; // convert from radians to degrees before writing out
    // temporarily delete values for depth of shower maximum and distance of shower maximum
    itsScenarioParams->itsDepthOfShowerMaximum = -1.0;
    itsScenarioParams->itsDistanceOfShowerMaximum = -1.0;

    cout << " CoREAS: The following parameters are used for the radio calculation\n\n";
    cout << (*itsScenarioParams); // write used parameter values into CORSIKA output log
    cout << "\n";

    // now calculate the (radio) core position at the height of itsCoreCoordinateVertical
    const ThreeVector showerAxis(cos(fAzimuth)*sin(fZenith),sin(fAzimuth)*sin(fZenith),-cos(fZenith));
    const double p = (itsScenarioParams->itsCoreCoordinateVertical-fObservationLevel)/cos(fZenith);
    fRadioCore = ThreeVector(0.0,0.0,fObservationLevel)-p*showerAxis;
    
//  cout << " CoREAS radioCore: " << fRadioCore << "\n" << endl;
     
    const ThreeVector firstInteraction(fFirstInteractionX, fFirstInteractionY, fFirstInteractionZ);
    const double distFirstIntToCore = (firstInteraction-fRadioCore).GetLength();
    fCoreHitTime = fFirstInteractionTime+distFirstIntToCore/SpeedOfLight; // SpeedOfLight is in cgs units

    // initialize ground area
    fGroundArea = new GroundArea(itsScenarioParams->itsTimeResolution, 
	 itsScenarioParams->itsTimeLowerBoundary, 
	 itsScenarioParams->itsTimeUpperBoundary, 
	 itsScenarioParams->itsAutomaticTimeBoundaries, 
	 itsScenarioParams->itsResolutionReductionScale, 
	 firstInteraction, 
	 fFirstInteractionTime-fCoreHitTime,
	 itsScenarioParams->itsAntennaPositions,
         fRadioCore, showerAxis);

    if (fObservationLevel > fGroundArea->GetLowestHeight())
      cout << " CoREAS warning: lowest observation level is at " << fObservationLevel << " cm and thus above lowest antenna height at " << fGroundArea->GetLowestHeight() << " cm!\n\n";

    fGroundLevelRefractivity = itsScenarioParams->itsGroundLevelRefractiveIndex-1.0;

  }
  else
  {    
    // error message is being printed out by ScenarioParams class
    exit (11);
  }

  // debugging
/*
  for (int i=0; i<20; ++i)
  {
    cout << "Between " << i*100 << " m and " << i*i*100 << "m is average refractivity noncurved ";
    fCurved = false;
    cout << GetEffectiveRefractiveIndexBetween(ThreeVector(i*100*100*100,0,i*100*100), ThreeVector(0,0,i*i*100*100))-1;
    fCurved = true;
    cout << " and curved " << GetEffectiveRefractiveIndexBetween(ThreeVector(i*100*100*100,0,i*100*100), ThreeVector(0,0,i*i*100*100))-1 << "\n";
  }
*/

}


void 
TPlotter::Clear () 
{
  fPrimaryTrack = true;
  fFirstInteraction = false;
  fFirstInteractionX = 0;
  fFirstInteractionY = 0;
  fFirstInteractionZ = 0;
  
  fAxisX = 0;
  fAxisY = 0;
  fAxisZ = 0;
  
  fCosZenith = 0;
  fSinZenith = 0;
  fCosAzimuth = 0;
  fSinAzimuth = 0;
  
  fXmax = 0;
  fEventNo = 0;
  fObservationLevel = 0;
  fHeightFirstInt = 0;
  fZenith = 0;
  fAzimuth = 0;
  fShowerEnergy = 0;
  fShowerPrimary = 0;

  fCORSIKA = 0;
}


void TPlotter::SetRunHeader(const crs::MRunHeader &header) 
{
  ostringstream ssTmp; 
  ssTmp.str(""); ssTmp << header.GetVersion(); fCorsikaVersion = ssTmp.str();
  ssTmp.str(""); ssTmp << setprecision(2) << ProgramVersion; fCoastVersion = ssTmp.str();
  
  SetRun ((int)header.GetRunID ());
  
  // create filename
  ostringstream fname;
  fname << "DAT" << setw (6) << setfill('0') << fRunNo;
  
  SetFileName(fname.str ());
  
  if ((int)header.GetNumberOfShowers() != 1) {
    cout << "\n\n\nCoREAS: Error, NSHOW is set to a value of " << (int)header.GetNumberOfShowers() <<". CoREAS only supports NSHOW = 1.\n\n\n" << endl;
    exit(11);
  }
  
}

void 
TPlotter::SetShowerHeader(const crs::MEventHeader &header) 
{
  if (header.GetSkimmingIncidence()) {
    fSkimming = true;
    fSkimmingAltitude = header.GetSkimmingAltitude() * cm;
    cout << " CoREAS: Detected skimming geometry with impact=" 
         << fSkimmingAltitude/km
         << " km" << endl;
  }
  
  // to flag the track of the primary particle
  fPrimaryTrack = true;
  fFirstInteraction = false;
  fFirstInteractionX = 0;
  fFirstInteractionY = 0;
  fFirstInteractionZ = 0;
  
  SetEvent(header.GetEventNumber());
  fShowerPrimary = (int)header.GetParticleId();
  
  SetShowerAxis(header.GetTheta()*rad, header.GetPhi()*rad);
  fShowerEnergy = header.GetEnergy()*GeV;
  
  fHeightFirstInt = header.GetZFirst()*cm;
  fObservationLevel = header.GetObservationHeight(header.GetNObservationLevels()-1)*cm;
  
  if (fHeightFirstInt<0) {
    if (fVerbosityLevel>=10) {
      cout << " CoREAS: height of first interaction is NEGATIVE, corrected." 
           << endl;
    }
    fHeightFirstInt = fabs(fHeightFirstInt);
  }

  // store magnetic field configuration for later access at writeout of .reas file
  fBx = header.GetBx();
  fBz = header.GetBz();

  cout  << " CoREAS: firstint " << fHeightFirstInt/m << "m\n"
        << " CoREAS: obslev " << fObservationLevel/m 
        << endl;

}


void 
TPlotter::SetShowerTrailer(const crs::MEventEnd &trailer) 
{
  fXmax = trailer.GetXmax()*g/cm/cm;
  if (!fSlant && fCosZenith!=0) 
    fXmax /= fCosZenith;
  
  cout << " CoREAS: --------- shower end ------ \n"
       << " CoREAS: xmax: " << fXmax/g*cm*cm << "\n";
  cout << endl;
}


void 
TPlotter::Init() 
{
  if (!fSkimming) 
    fCORSIKA = new TCorsika(fZenith, 
          fObservationLevel,
          fSlant, fCurved);
  else 
    fCORSIKA = new TCorsika(fSkimmingAltitude,
          fObservationLevel);
  
  if (fVerbosityLevel > 0) 
    fCORSIKA->SetVerbosityLevel(fVerbosityLevel);

  if( fCORSIKA == 0 )
    cout << "Warning fCORSIKA not defined" << endl;
  else
    fCORSIKA->SetHeightOfFirstInteraction(fHeightFirstInt);

  if (fCurved && (not fSlant)) { // curved without slant is not supported in REAS
    cout << "\n\n"
         << " ###################################\n"
         << "   CURVED without SLANT is not supported   \n"
         << "    please switch on the SLANT option \n"
         << " ###################################\n\n"
         << endl;
    exit (11);
  }
  
  if (fStackInput || fPreshower) { // STACKIN OPTION NEEDS SPECIAL TREATMENT !!!!
    fPrimaryTrack = false;
    fFirstInteraction = true;
  }
  
  if (fVerbosityLevel > 1 && fCurved) {
    fCORSIKA->DumpAtmosphereTable();
  }
}





void 
TPlotter::Write () 
{
  if (fVerbosityLevel >= 2) {
    cout << " CoREAS: Write " << endl;
    /*
    ostringstream gif_dir;
    gif_dir << fFileName 
    << "_"
      << fEventNo;
    */
  }

  // Write some values into the ScenarioParams and then output the file
  itsScenarioParams->itsDepthOfShowerMaximum = fXmax/g*cm*cm;
  itsScenarioParams->itsDistanceOfShowerMaximum = (fCORSIKA->GetDistanceOfPoint(fRadioCore.GetX(), fRadioCore.GetY(), fRadioCore.GetZ())-fCORSIKA->GetDistanceOfSlantDepth(fXmax))/cm;
  itsScenarioParams->WriteToFiles();                                                                                             

  // output direct radio results
  ostringstream fname;
  fname << fDirectory << "SIM" << setw (6) << setfill('0') << fRunNo << "_coreas";
  ostringstream fsimName;
  fsimName << "SIM" << setw (6) << setfill('0') << fRunNo << "_coreas";
  fGroundArea->WriteResults(fname.str(), fsimName.str());
  
  // clean up
  Clear ();
  delete itsScenarioParams;   
}




// #define Rotate(x, y, sx, sy, cosa, sina) { 
inline
void 
TPlotter::Rotate(double x, double y, double z,
     double &sx, double &sy, double &sz,
     int inverse) const
{

  // rotate around z by azimuth
  double sx_ =             x*fCosAzimuth + inverse * y*fSinAzimuth;
  double sy_ = - inverse * x*fSinAzimuth +           y*fCosAzimuth; 
  double sz_ =   z;
  
  // rotate around y` by zenith
  double sx__ =             sx_*fCosZenith - inverse * sz_*fSinZenith;
  double sy__ = sy_;
  double sz__ = + inverse * sx_*fSinZenith +           sz_*fCosZenith; 

  // rotate back around z`` by -azimuth 
  sx =             sx__*fCosAzimuth - inverse * sy__*fSinAzimuth;
  sy = + inverse * sx__*fSinAzimuth +           sy__*fCosAzimuth; 
  sz =   sz__;

}


void 
TPlotter::AddTrack(const crs::CParticle &pre, 
       const crs::CParticle &post) 
{  
  if (fVerbosityLevel >= 9) {
    /*  double dX = pre.x-fFirstInteractionX; 
  double dY = pre.y-fFirstInteractionY; 
  double dZ = pre.z-fFirstInteractionZ; 
  double length = sqrt (dX*dX + dY*dY + dZ*dZ);*/
    cout << " CoREAS:   PRE> "
   << " id=" << (int)pre.particleId
   << " E=" << pre.energy
   << " w=" << pre.weight
   << " x=" << pre.x/km << "km"
   << " y=" << pre.y/km << "km"
   << " z=" << pre.z/km << "km"
   << " t=" << pre.time*s/ns << "ns"
   << " X=" << pre.depth << "g/cm^2"
      //         << " v= " <<  length/(pre.time*s/ns-fFirstInteractionTime*s/ns)/crs::cSpeedOfLight 
   << endl;
  }
  
  if (fVerbosityLevel >= 10) {
    /*  double dX = post.x-fFirstInteractionX; 
  double dY = post.y-fFirstInteractionY; 
  double dZ = post.z-fFirstInteractionZ; 
  double length = sqrt (dX*dX + dY*dY + dZ*dZ);*/
    cout << " CoREAS:  POST> "
   << " id=" << (int)post.particleId
   << " E=" << post.energy
   << " w=" << post.weight
   << " x=" << post.x/km << "km"
   << " y=" << post.y/km << "km"
   << " z=" << post.z/km << "km"
   << " t=" << post.time*s/ns << "ns"
   << " X=" << post.depth << "g/cm^2"
      //        << " v= " <<  length/(post.time*s/ns-fFirstInteractionTime*s/ns)/crs::cSpeedOfLight
   << endl;
  }
    
  /*
    Skip the track of the primary particle, which is in a different 
    reference system as the shower !!!
  */
  //static double cosZenith_track = fCosZenith;
  //static double azimuth_track = fAzimuth;
  if (fPrimaryTrack) {
    if (fVerbosityLevel >= 2) {
      cout << " CoREAS: Primary track  " << endl;
    }    
    return;
  }
  
  /* Skip the particle if pre and post times are identical */
  if (pre.time-post.time == 0.0)
    return;  
  
  if (fFirstInteraction) {
    
    fFirstInteraction = false;
    SetFirstInteraction(pre.x, pre.y, pre.z, pre.depth*g/cm/cm, pre.time);
    
    if (fVerbosityLevel >= 2) {
      cout << " CoREAS: Primary track end-pos: x=" << fFirstInteractionX/m 
     << " y=" << fFirstInteractionY/m 
     << " z=" << fFirstInteractionZ/m << endl;
    }
  }
  
  const int particleId = abs((int)pre.particleId);
  
  // check if particles to be counted at all
  if (!fParticles.count(particleId)) 
    return;
  
  // calculate direct radio
  // find a better way to specify charge signs
  double particleCharge = gParticleCharge[abs((int)pre.particleId)-1];
      
  ThreeVector prePosition(pre.x, pre.y, pre.z);
  ThreeVector postPosition(post.x, post.y, post.z);

  double preIndex = GetEffectiveRefractiveIndexBetween(prePosition,prePosition);
  double postIndex = GetEffectiveRefractiveIndexBetween(postPosition,postPosition);
      
/*double preGamma = pre.energy/gParticleMass[abs((int)pre.particleId)-1];
  double postGamma = post.energy/gParticleMass[abs((int)post.particleId)-1];
      
  double preBetaValue = sqrt(1.-1./(preGamma*preGamma));
  double postBetaValue = sqrt(1.-1./(postGamma*postGamma));
*/    
  ThreeVector currDirection(post.x-pre.x, post.y-pre.y, post.z-pre.z);
  currDirection*=1./currDirection.GetLength();
      
//ThreeVector preBeta = preBetaValue*currDirection;     // neglects curvature along track from pre to post!
//ThreeVector postBeta = postBetaValue*currDirection;   // neglects curvature along track from pre to post!

  double corrBetaValue = (postPosition-prePosition).GetLength()/(SpeedOfLight*(post.time-pre.time));
  ThreeVector corrBeta = corrBetaValue * currDirection;

/*cout << "--------------------\n";
  cout << "pre.particleId: " << pre.particleId << "\n";
  cout << "post.particleId: " << post.particleId << "\n";
  cout << "particleCharge: " << particleCharge << "\n";
  cout << "pre.energy: " << pre.energy << "\n";
  cout << "post.energy: " << post.energy << "\n";
  cout << "preGamma: " << preGamma << "\n";
  cout << "postGamma: " << postGamma << "\n";
  cout << "preBetaValue: " << preBetaValue << "\n";
  cout << "postBetaValue: " << postBetaValue << "\n";
  cout << "prePosition: " << prePosition.GetX() << ", " << prePosition.GetY() << ", " << prePosition.GetZ() << "\n";	   
  cout << "postPosition: " << postPosition.GetX() << ", " << postPosition.GetY() << ", " << postPosition.GetZ() << "\n";      
  cout << "preBeta: " << preBeta.GetX() << ", " << preBeta.GetY() << ", " << preBeta.GetZ() << "\n";	   
  cout << "postBeta: " << postBeta.GetX() << ", " << postBeta.GetY() << ", " << postBeta.GetZ() << "\n";      
  cout << "corrBeta: " << corrBeta.GetX() << ", " << corrBeta.GetY() << ", " << corrBeta.GetZ() << "\n";      
*/

  const double approxthreshold = 1.0e-3; // set threshold for application of ZHS-like approximation
  ThreeVector startE;
  ThreeVector endE;
  
  const vector<GroundElement*>& groundElements = fGroundArea->GetActiveGroundElementPointers();
  for (vector<GroundElement*>::const_iterator currGroundElement = groundElements.begin();
      currGroundElement != groundElements.end(); ++currGroundElement)  
  {
    ThreeVector obsPosition = (*currGroundElement)->GetPosition();

    // check if particle contribution should be skipped due to special mode of observer bin
    if ((*currGroundElement)->GetMode() != Normal)
    {
      // check whether to apply cuts in particle gamma
      if ((*currGroundElement)->GetMode() == Gamma)
      {
        double preGamma = pre.energy/gParticleMass[abs((int)pre.particleId)-1];
        if ((preGamma < (*currGroundElement)->GetModeDataLow()) || (preGamma >= (*currGroundElement)->GetModeDataHigh())) {
          continue;
        }
      }
      // check whether to apply cuts in particle height asl
      else if ((*currGroundElement)->GetMode() == Height)
      {
        double height = GetHeightAtPosition(ThreeVector(pre.x, pre.y, pre.z));
        if ((height < (*currGroundElement)->GetModeDataLow())
           || (height >= (*currGroundElement)->GetModeDataHigh())) {
          continue;
        }
      }
    }    
        
    ThreeVector preN = obsPosition-prePosition;
    double preDistance = preN.GetLength();
    preN *= 1./preDistance;
    double preT = pre.time-fCoreHitTime;
    double preResT = preT+GetEffectiveRefractiveIndexBetween(prePosition,obsPosition)*preDistance/SpeedOfLight;
        
    ThreeVector postN = obsPosition-postPosition;
    double postDistance = postN.GetLength();
    postN *= 1./postDistance;
    double postT = post.time-fCoreHitTime;
    double postResT = postT+GetEffectiveRefractiveIndexBetween(postPosition,obsPosition)*postDistance/SpeedOfLight;
        
    double preDoppler = (1.0 - preIndex*corrBeta.DottedWith(preN));
    double postDoppler = (1.0 - postIndex*corrBeta.DottedWith(postN));

    // apply special treatment for refractive index, this is motivated by a finite observation time resolution and a ZHS-like limit calculation
    if ( (fGroundLevelRefractivity > 0.0) && ((fabs(preDoppler)<approxthreshold) || (fabs(postDoppler)<approxthreshold)) )
    {
      ThreeVector midPosition = 0.5*(prePosition+postPosition);
      ThreeVector midN = obsPosition-midPosition;
      double midDistance = midN.GetLength();
      midN *= 1./midDistance;
      double midIndex = GetEffectiveRefractiveIndexBetween(midPosition,midPosition);
      double midDoppler = (1.0 - midIndex*corrBeta.DottedWith(midN));

      // check if midDoppler has become zero because of numerical limitations
      if (midDoppler == 0)
      {
        // redo calculation with higher precision
        long double indexL = midIndex;
        long double betaX = corrBeta.GetX();
        long double betaY = corrBeta.GetY();
        long double betaZ = corrBeta.GetZ();
        long double nX = midN.GetX();
        long double nY = midN.GetY();
        long double nZ = midN.GetZ();
        long double doppler = 1.0l - indexL * (betaX*nX+betaY*nY+betaZ*nZ);
        midDoppler = doppler;
      }

      startE = pre.weight * particleCharge * UnitCharge * ( midN.CrossedWith( midN.CrossedWith(corrBeta))) /
               (SpeedOfLight * midDistance * midDoppler);      
                               
      endE = -1.0 * startE;
        
      // force deltaT to be proportional to the Doppler factor (i.e. contribution = constant * track length), analogous to ZHS
      const double midResT = 0.5*(preT+postT)+GetEffectiveRefractiveIndexBetween(midPosition,obsPosition)*midDistance/SpeedOfLight;
      const double segmentlength = (postPosition-prePosition).GetLength();
      double deltaT = segmentlength/(SpeedOfLight*corrBetaValue)*fabs(midDoppler);
        
      if (preResT < postResT)
      { // startE arrives earlier
        preResT = midResT-0.5*deltaT;
        postResT = midResT+0.5*deltaT;
      }
      else
      { // endE arrives earlier
        postResT = midResT-0.5*deltaT;
        preResT = midResT+0.5*deltaT;
      }

      const long double gridresolution = (*currGroundElement)->GetTimeResolution();
      deltaT = postResT-preResT; //have to recalculate, because deltaT can be negative
        
      // redistribute contributions over time scale defined by the observation time resolution
      if (fabs(deltaT) < gridresolution)
      {
        startE *= fabs(deltaT/gridresolution);
        endE *= fabs(deltaT/gridresolution);
          
        const long startBin = static_cast<long>(floor(preResT/gridresolution+0.5l));
        const long endBin = static_cast<long>(floor(postResT/gridresolution+0.5l));
        const double startBinFraction = (preResT/gridresolution)-floor(preResT/gridresolution);
        const double endBinFraction = (postResT/gridresolution)-floor(postResT/gridresolution);

        // only do timing modification if contributions would land in same bin
        if (startBin == endBin)
        {
          // if startE arrives before endE
          if (deltaT >= 0)
          { 
            if ((startBinFraction >= 0.5) && (endBinFraction >= 0.5)) // both points left of bin center
            {
              preResT -= gridresolution; // shift startE to previous gridpoint
            }
            else if ((startBinFraction < 0.5) && (endBinFraction < 0.5)) // both points right of bin center
            {
              postResT += gridresolution; // shift endE to next gridpoint
            }
            else // points on both sides of bin center
            {
              const double leftDist = 1.0-startBinFraction;
              const double rightDist = endBinFraction;
              // check if asymmetry to right or left
              if (rightDist >= leftDist)
              {
                postResT += gridresolution; // shift endE to next gridpoint
              }
              else
              {
                preResT -= gridresolution;  // shift startE to previous gridpoint
              }
            }          
          }
          else // if endE arrives before startE
          {
            if ((startBinFraction >= 0.5) && (endBinFraction >= 0.5)) // both points left of bin center
            {
              postResT -= gridresolution; // shift endE to previous gridpoint
            }
            else if ((startBinFraction < 0.5) && (endBinFraction < 0.5)) // both points right of bin center
            {
              preResT += gridresolution; // shift startE to next gridpoint
            }
            else // points on both sides of bin center
            {
              const double leftDist = 1.0-endBinFraction;
              const double rightDist = startBinFraction;
              // check if asymmetry to right or left
              if (rightDist >= leftDist)
              {
                preResT += gridresolution; // shift startE to next gridpoint
              }
              else
              {
                postResT -= gridresolution;  // shift startE to previous gridpoint
              }
            }
          }
        } 
      }
    }
    else
    {
      // refractive index is unity or not near Cherenkov angle

      // check if preDoppler has become zero in case of refractive index of unity because of numerical limitations
      if (preDoppler == 0)
      {
        // redo calculation with higher precision
        long double indexL = preIndex;
        long double betaX = corrBeta.GetX();
        long double betaY = corrBeta.GetY();
        long double betaZ = corrBeta.GetZ();
        long double nX = preN.GetX();
        long double nY = preN.GetY();
        long double nZ = preN.GetZ();
        long double doppler = 1.0l - indexL * (betaX*nX+betaY*nY+betaZ*nZ);
        preDoppler = doppler;
      }

      startE = pre.weight * particleCharge * UnitCharge * ( preN.CrossedWith( preN.CrossedWith(corrBeta))) / 
               (SpeedOfLight * preDistance * preDoppler);

      // check if postDoppler has become zero in case of refractive index of unity because of numerical limitations
      if (postDoppler == 0)
      {
        // redo calculation with higher precision
        long double indexL = postIndex;
        long double betaX = corrBeta.GetX();
        long double betaY = corrBeta.GetY();
        long double betaZ = corrBeta.GetZ();
        long double nX = postN.GetX();
        long double nY = postN.GetY();
        long double nZ = postN.GetZ();
        long double doppler = 1.0l - indexL * (betaX*nX+betaY*nY+betaZ*nZ);
        postDoppler = doppler;
      }
      
      endE = -1.0 * post.weight * particleCharge * UnitCharge * ( postN.CrossedWith( postN.CrossedWith(corrBeta))) / 
             (SpeedOfLight * postDistance * postDoppler);

      // if preDoppler or postDoppler are below a certain threshold, redistribute contributions over two consecutive bins (take into account finite detector time resolution)
      if ((preDoppler<1.e-9) || (postDoppler<1.e-9))
      {
        const long double gridresolution = (*currGroundElement)->GetTimeResolution();
        double deltaT = postResT-preResT;
        if (fabs(deltaT) < gridresolution)
        {
          startE *= fabs(deltaT/gridresolution);
          endE *= fabs(deltaT/gridresolution);
          
          const long startBin = static_cast<long>(floor(preResT/gridresolution+0.5l));
          const long endBin = static_cast<long>(floor(postResT/gridresolution+0.5l));
          const double startBinFraction = (preResT/gridresolution)-floor(preResT/gridresolution);
          const double endBinFraction = (postResT/gridresolution)-floor(postResT/gridresolution);

          // only do timing modification if contributions would land in same bin
          if (startBin == endBin)
          {
            if ((startBinFraction >= 0.5) && (endBinFraction >= 0.5)) // both points left of bin center
            {
              preResT -= gridresolution; // shift startE to previous gridpoint
            }
            else if ((startBinFraction < 0.5) && (endBinFraction < 0.5)) // both points right of bin center
            {
              postResT += gridresolution; // shift endE to next gridpoint
            }
            else // points on both sides of bin center
            {
              const double leftDist = 1.0-startBinFraction;
              const double rightDist = endBinFraction;
              // check if asymmetry to right or left
              if (rightDist >= leftDist)
              {
                postResT += gridresolution; // shift endE to next gridpoint
              }
              else
              {
                preResT -= gridresolution;  // shift startE to previous gridpoint
              }
            }
          }
        }          
      }
    }
  
    if ((*currGroundElement)->GetMode() == Pattern)
    {
      //get angles in ResponseTable interface coordinate system (CCW from east, CORSIKA: CCW from North)
      const double startTheta = atan2(sqrt(preN.GetX()*preN.GetX()+preN.GetY()*preN.GetY()), -preN.GetZ() );
      const double startPhi = atan2(preN.GetY(), preN.GetX()) + Pi/2.;
      const double endTheta = atan2(sqrt(postN.GetX()*postN.GetX()+postN.GetY()*postN.GetY()), -postN.GetZ() );
      const double endPhi = atan2(postN.GetY(), postN.GetX()) + Pi/2.;

      //get weight factors from antenna response pattern 
      //if no response pattern is used (*currGroundElement)->GetEffectiveAntennaArea() returns 1
      const double weightStartE = sqrt((*currGroundElement)->GetEffectiveAntennaArea(startTheta, startPhi));
      const double weightEndE = sqrt((*currGroundElement)->GetEffectiveAntennaArea(endTheta, endPhi));
      startE *= weightStartE;
      endE *= weightEndE;
    }

    (*currGroundElement)->CollectValuePair(EFieldDataPoint(preResT,startE),EFieldDataPoint(postResT,endE));
  }
}

void 
TPlotter::AddInteraction(const crs::CInteraction& interaction) 
{
  if (fVerbosityLevel >= 10) {
    cout << " CoREAS: ";
    interaction.Dump();
  }
  
  if (fPrimaryTrack) 
    fFirstInteraction = true;
  else 
    fFirstInteraction = false;
  
  fPrimaryTrack = false;
}


void 
TPlotter::SetFirstInteraction(const double x, const double y, const double z, 
            const double X, const double t) 
{
  fFirstInteractionX = x;
  fFirstInteractionY = y;
  fFirstInteractionZ = z;
  fFirstInteractionTime = t;
  fFirstInteractionSlantDepth = X;
  if (X<0) {
    cout << " CoREAS: WARNING FirstInteractionDepth smaller than 0 (" << X/g*cm*cm << " g/cm2) " << endl;
    fFirstInteractionSlantDepth = 0;
  }
  //fFirstInteractionDist = fCORSIKA->GetDistanceOfSlantDepth(fFirstInteractionSlantDepth);
  fFirstInteractionDist = fCORSIKA->GetDistanceOfPoint(x, y, z);
  
  // now that first interaction point is known, set up radio simultion
  InitializeRadioSimulation();                                  
  
  cout << " CoREAS: SetFirstInteraction Dist=" << fFirstInteractionDist/km << "km "
       << ", x=" << fFirstInteractionX/km << "km "
       << ", y=" << fFirstInteractionY/km << "km "
       << ", z=" << fFirstInteractionZ/km << "km "
       << ", SlantDepth=" << fFirstInteractionSlantDepth/g*cm*cm << "g/cm2 "
       << ", Time=" << fFirstInteractionTime << "s "
       << "\n\n";
}


void 
TPlotter::SetShowerAxis(const double zenith, const double azimuth) 
{  
  cout << " \n\n CoREAS: SetShowerAxis zenith=" << zenith/deg << "deg"  
       << "   azimuth=" << azimuth/deg << "deg"  
       << endl;
  fZenith = zenith;
  fAzimuth = azimuth;
  fCosZenith = cos(zenith);
  fSinZenith = sin(zenith);
  fCosAzimuth = cos(azimuth);
  fSinAzimuth = sin(azimuth);
  fAxisZ = fCosZenith;
  fAxisX = fSinZenith * fCosAzimuth;;
  fAxisY = fSinZenith * fSinAzimuth;
}


double TPlotter::GetEffectiveRefractiveIndexBetween(const ThreeVector& p1, const ThreeVector& p2) const
{
  if (fGroundLevelRefractivity == 0.0)
    return 1.0; // more efficient than averaging (although same result)

  if ((fCurved) && (fZenith>=75.0))
  {
    // check if singular point
    if (p1 == p2)
      return 1.0 + fGroundLevelRefractivity * rhof_(GetHeightAtPosition(p1))/rhof_(0.0);

    // integrate the density of the atmosphere along the line of sight between p1 and p2 to calculate averaged refractive index
    const double totalDistance = (p2-p1).GetLength();
    const ThreeVector travelDirection = (p2-p1).GetDirection();
    int numSteps = static_cast<int>(totalDistance/500000.)+1; // one step every 5000 metres
    if (numSteps<10)
      numSteps=10; // do at least 10 steps
    const double stepSize = totalDistance/numSteps;
    double traversedAtmosphericDepth = 0.0;
    for (int jj = 1; jj<=numSteps; ++jj)
    {
      double localHeight = GetHeightAtPosition(p1+stepSize*jj*travelDirection);
      traversedAtmosphericDepth += stepSize * rhof_(localHeight);
    }

    double srcDensity = rhof_(GetHeightAtPosition(p1));
    double srcRefractivity = fGroundLevelRefractivity * srcDensity/rhof_(0.0);
    return 1.0 + srcRefractivity/(srcDensity * totalDistance) * traversedAtmosphericDepth; // averaged refractive index
  }
  else
  {
    // average over refractive index between these two points, return the effective value
    double h1 = p1.GetZ();
    double h2 = p2.GetZ();
  
    if (h1 == h2) // singular point
      return 1.0 + fGroundLevelRefractivity * rhof_(h1)/rhof_(0.0);

    return 1.0 + fGroundLevelRefractivity/(rhof_(0.0) * (h2-h1)) * (thick_(h1)-thick_(h2)); // averaged refractive index
  }
}

double TPlotter::GetHeightAtPosition(const ThreeVector& p1) const
{
  #warning Cross-check and unify TPlotter::GetHeightAtPosition() and TCorsika::GetHeightOfPoint()
  if (fCurved)
    return sqrt(p1.GetX()*p1.GetX() + p1.GetY()*p1.GetY() + p1.GetZ()*p1.GetZ() + 2*cRearth*p1.GetZ() + cRearth*cRearth) - cRearth;
  else
    return p1.GetZ();
}
