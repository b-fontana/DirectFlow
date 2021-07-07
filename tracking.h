#ifndef TRACKING_H
#define TRACKING_H

#include "./functions.h"
#include <memory>

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
  
 Track(unsigned pNstepsUsed,
       Vec<double> pEnergies, Vec<XYZ> pPositions, Vec<XYZ> pDirections)
   : mNstepsUsed(pNstepsUsed),
    mEnergies(pEnergies), mPositions(mPositions), mDirections(pDirections) {};
    
  Vec<double> energies() const { return mEnergies; }
  Vec<XYZ> positions() const { return mPositions; }
  Vec<XYZ> directions() const { return mDirections; }
  unsigned steps_used() const { return mNstepsUsed; }
  
 private:
  Vec<double> mEnergies;
  Vec<XYZ> mPositions;
  Vec<XYZ> mDirections;
  unsigned mNstepsUsed;
};

////////////////////////////////////////////
//simulates the particle trajectory
////////////////////////////////////////////
class SimParticle {
 public:
  template <typename T> 
    using Vec = std::vector<T>;
  using XYZ = ROOT::Math::XYZVector;
  
  enum class TrackMode { Euler=0, NMODES };
  
 SimParticle(Particle pParticle)
   : mParticle(pParticle) {};
 SimParticle(Particle pParticle, unsigned pNsteps, double pStepSize)
   : mParticle(pParticle), mNsteps(pNsteps), mStepSize(pStepSize) {};

  const Track& track(const Magnets&, SimParticle::TrackMode);
    
 private:
  Vec<Track> mTracks(TrackMode::NMODES);
  Vec<bool> mTrackCheck(TrackMode::NMODES);
  
  Particle mParticle;

  static constexpr double mCvelocity = 29979245800.0; // (cm/s)
  static constexpr double mEcharge = 1.602176565E-19; // C = A*s
  unsigned mNsteps = 3000;
  float mStepSize = .1f;

  const Track& track_euler(const Magnets& );
  // const Track& track_rungekutta4(const Magnets& );
};

#endif // TRACKING_H
