#ifndef TRACKING_H
#define TRACKING_H

//#include "./functions.h"
#include "./geometry.h"
#include "./utils.h"
#include <memory>
#include <cmath>

#include "TROOT.h"
#include "Math/Vector3D.h" // XYZVector
#include "Math/Vector4D.h" // LorentzVector

namespace tracking {
  enum TrackMode { Euler=0, RungeKutta4, NMODES };
}

////////////////////////////////////////////
//TParticle has too much info
////////////////////////////////////////////
struct Particle {
public:
  ROOT::Math::XYZVector pos;
  ROOT::Math::XYZVector mom;
  double energy;
  double mass;
  int charge;
};

////////////////////////////////////////////
//simple structure to store track information
////////////////////////////////////////////
class Track {
public:
  template <typename T> 
  using Vec = std::vector<T>;  
  using XYZ = ROOT::Math::XYZVector;

 Track(): mNstepsUsed(0), mEnergies(0), mPositions(0), mMomenta(0) {};
  
  Track(unsigned pNstepsUsed,
	Vec<double> pEnergies, Vec<XYZ> pPositions, Vec<XYZ> pMomenta)
    : mNstepsUsed(pNstepsUsed),
    mEnergies(pEnergies), mPositions(pPositions), mMomenta(pMomenta) {};

  unsigned steps_used() const { return mNstepsUsed; }
  Vec<double> energies() const { return mEnergies; }
  Vec<XYZ> positions() const { return mPositions; }
  Vec<XYZ> momenta() const { return mMomenta; }
  
private:
  unsigned mNstepsUsed;
  Vec<double> mEnergies;
  Vec<XYZ> mPositions;
  Vec<XYZ> mMomenta;
};

////////////////////////////////////////////
//simulates the particle trajectory
////////////////////////////////////////////
class SimParticle {
public:
  template <typename T> 
  using Vec = std::vector<T>;
  using XYZ = ROOT::Math::XYZVector;
  using Ltz = ROOT::Math::PxPyPzMVector;

  SimParticle(Particle pParticle)
    : mParticle(pParticle),
      mTracks(tracking::TrackMode::NMODES), mTrackCheck(tracking::TrackMode::NMODES, false) {}
  
  SimParticle(Particle pParticle, unsigned pNsteps, double pStepSize)
    : mParticle(pParticle),
      mTracks(tracking::TrackMode::NMODES), mTrackCheck(tracking::TrackMode::NMODES, false),
      mNsteps(pNsteps), mStepSize(pStepSize) {};

  const Track& track(const MagnetSystem&, tracking::TrackMode, double, float) &;
  //no copies of big objects, so forbid calling 'track()' on a temporary object
  Track track(const MagnetSystem&, tracking::TrackMode, double) && = delete;
    
private:
  Particle mParticle;
  Vec<Track> mTracks;
  Vec<bool> mTrackCheck;
  
  static constexpr double mSpeedOfLight = 29979245800.0; // (cm/s)
  static constexpr double mEcharge = 1.602176565E-19; // C = A*s
  unsigned mNsteps = 3000;
  double mStepSize = 0.;

  XYZ calc_relativistic_velocity(const XYZ&, double, double) const;
  XYZ calc_lorentz_force(double, const XYZ&, const XYZ&) const;
  Track track_euler( const MagnetSystem&, double, float );
  Track track_rungekutta4( const MagnetSystem&, double );
};

#endif // TRACKING_H
