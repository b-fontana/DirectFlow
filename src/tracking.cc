#include "include/tracking.h"

const Track& SimParticle::track(const Magnets& magnets, tracking::TrackMode mode ) {

  const Track* t = nullptr;
  
  if(mode == tracking::TrackMode::Euler) {
      
    if(mTrackCheck[tracking::TrackMode::Euler] == false) {
      mTracks[tracking::TrackMode::Euler] = Track();
      mTracks[tracking::TrackMode::Euler] = track_euler(magnets);
      mTrackCheck[tracking::TrackMode::Euler] = true;
    }
    t = &( mTracks[tracking::TrackMode::Euler] );
      
  }
  /* else if(mode == TrackMode::RungeKutta4) { */
      
  /*   if(mTrackCheck[TrackMode::RungeKutta4] == false) { */
  /* 	mTracks[TrackMode::RungeKutta4] = track_rungekutta4(); */
  /* 	mTrackCheck[TrackMode::RungeKutta4] = true; */
  /*   } */
  /*   const Track& t = mTracks[TrackMode::RungeKutta4]; */
      
  /* } */
  else {
      
    throw std::invalid_argument("The tracking mode specified is not supported.");

  }
    
  return *t;
}

Track SimParticle::track_euler(const Magnets& magnets)
{
  /*
  Track a particle through the magnetic field with
  pseudo-rapidity eta
  transverse momentum pT (GeV/c)
  azimuthal angle phi (rad)
  charge q (0 = +, 1 = -)
  z-vertex position z (cm)
  particle mass m (GeV/c^2)
  number of tracking steps N_steps
  step size step_size (cm)
  */
  
  double charge = mParticle.charge * mEcharge; // C = A*s

  ROOT::Math::PxPyPzMVector lvector( mParticle.mom.X(), mParticle.mom.Y(), mParticle.mom.Z(), mParticle.mass );
  XYZ partPos = mParticle.pos; //cm
  XYZ partMom = lvector.Vect(); // Gev/c
  // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) -> (cm/s)
  XYZ partVel = partMom * mSpeedOfLight / ( lvector.Gamma() * mParticle.mass );

  double delta_t = mStepSize/(mSpeedOfLight*lvector.Beta()); // s

  Vec<double> energies;
  Vec<XYZ> positions;
  Vec<XYZ> momenta;
  energies.reserve(mNsteps);
  positions.reserve(mNsteps);
  momenta.reserve(mNsteps);

  unsigned nStepsUsed = 0;
  while(nStepsUsed<mNsteps)
    {
      positions.push_back( partPos );
      momenta.push_back( partMom );

      // direction to move without magnetic field in cm
      XYZ posIncr = partMom * ( mStepSize / TMath::Sqrt(partMom.Mag2()) );

      // new position after delta_t without magnetic field
      XYZ partPosNext = partPos + posIncr;
      // center of begin and stop vector without magnetic field
      XYZ partPosMid = (partPos + partPosNext) * 0.5;     

      XYZ Bfield = magnets.field(partPosMid, 350.);

      if(Bfield.Mag2() == 0.0)
	  partPos = partPosNext;
      else
	{
	  // F = q*v X B
	  XYZ force = charge*partVel.Cross(Bfield); // (A*s)*(cm/s)*(kg/(A*s*s)) = (cm*kg)/(s*s)

	  // F = (dp/dt)
	  XYZ forceDelta = force * delta_t * 1.8708026E16; // (cm*kg)/(s*s) * delta_p (cm*kg)/s -> (GeV/c)
	  XYZ partMomNext = partMom + forceDelta;
	  
	  double mag0 = TMath::Sqrt(partMom.Mag2());
	  double mag1 = TMath::Sqrt(partMomNext.Mag2());
	  partMomNext *= mag0 / mag1; // make sure that total momentum doesn't change

	  XYZ momDelta = partMomNext * ( mStepSize / mag1); // direction to move with magnetic field in cm
	  partPos += momDelta; // new position after delta_t with magnetic field

	  partMom = partMomNext;

	  // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) (cm/s)
	  partVel = partMom * mSpeedOfLight / ( lvector.Gamma() * mParticle.mass );
	}

      energies.push_back( TMath::Sqrt(partMom.Mag2() + 0.938*0.938) );
      
      ++nStepsUsed;
      
      if(fabs(partPos.Z()) > 9000.0) break;
    }

  return Track(nStepsUsed, energies, positions, momenta);
}
