/* $Id: TParticleBlockEntry.cc,v 1.1.1.1 2007-07-31 07:00:48 rulrich Exp $   */

#include <crs/TParticleBlockEntry.h>
using namespace crs;


TParticleBlockEntry::TParticleBlockEntry (const CREAL *data, bool thinn) {

  fData = data;
  fThinned = thinn;
}



TParticleBlockEntry::TParticleBlockEntry (const TParticleBlockEntry& p) {

  fData = p.fData;
  fThinned = p.fThinned;
}

/*
  #define IsParticle(p) (((p) > 0) && ((p) < 100000))
  #define IsNucleus(p) (((p) >= 100 000) && ((p) < 9900000))
  #define IsCherenkov(p) ((p) >= 9900000)
*/




bool TParticleBlockEntry::IsParticle () const {
  int p = GetParticleId();
  if (p<0) 
    p = -p;
  return (p<200 && p!=75 && p!=76 && p!=85 && p!=86 && p!=95 && p!=96);
}

bool TParticleBlockEntry::IsNucleus () const {
  int p = GetParticleId();
  if (p<0) 
    p = -p;
  return (p>=200);
}

bool TParticleBlockEntry::IsCherenkov () const {
  int p = GetParticleId();
  if (p<0) 
    p = -p;
  return (p>=9900);
}

bool TParticleBlockEntry::IsMuonProductionInfo () const {
  int p = GetParticleId ();
  if (p<0) 
    p = -p;
  return (p == 75 || p == 76 || p == 85 || p == 86 || p == 95 || p == 96);
}

bool TParticleBlockEntry::IsEmpty () const {
  return (0 == (int)fData [0]);
}


ParticleType TParticleBlockEntry::GetType () const {
    
  if (IsParticle ())
    return eParticle;

  if (IsCherenkov ())
    return eCherenkov;

  if (IsNucleus ())
    return eNucleus;

  if (IsMuonProductionInfo ())
    return eMuonProductionInfo;
    
  if (IsEmpty ())
    return eEmpty;
    
  return eUnknown;
}



std::string TParticleBlockEntry::GetParticleName () const {

  return std::string ("undefined");
}




