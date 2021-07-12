#include "include/tracking.h"

SimParticle::XYZ SimParticle::calc_relativistic_velocity(const XYZ& mom, double gamma, double mass) const {
  return mom * mSpeedOfLight / ( gamma * mass );
}

SimParticle::XYZ SimParticle::calc_lorentz_force(double charge, const XYZ& vel, const XYZ& b) const {
  return charge*vel.Cross(b); // (A*s)*(cm/s)*(kg/(A*s*s)) = (cm*kg)/(s*s)
}

const Track& SimParticle::track(const Magnets& magnets, tracking::TrackMode mode, double scale ) {
  const Track* t = nullptr;
  using m = tracking::TrackMode;
  
  if(mode == m::Euler) {
      
    if(mTrackCheck[m::Euler] == false) {
      mTracks[m::Euler] = track_euler(magnets, scale);
      mTrackCheck[m::Euler] = true;
    }
    t = &( mTracks[m::Euler] );
      
  }
  else if(mode == m::RungeKutta4) {
      
    if(mTrackCheck[m::RungeKutta4] == false) {
      mTracks[m::RungeKutta4] = track_rungekutta4(magnets, scale);
  	mTrackCheck[m::RungeKutta4] = true;
    }
    t = &( mTracks[m::RungeKutta4] );
      
  }
  else   
    throw std::invalid_argument("The tracking mode specified is not supported.");

  return *t;
}

Track SimParticle::track_euler(const Magnets& magnets, double scale)
{ 
  double charge = mParticle.charge * mEcharge; // C = A*s

  Ltz initLorentzVec( mParticle.mom.X(), mParticle.mom.Y(), mParticle.mom.Z(), mParticle.mass );
  XYZ partPos = mParticle.pos; //cm
  XYZ partMom = initLorentzVec.Vect(); // Gev/c
  // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) -> (cm/s)
  XYZ partVel = calc_relativistic_velocity(partMom, initLorentzVec.Gamma(), mParticle.mass);

  double deltaT = mStepSize / ( mSpeedOfLight * initLorentzVec.Beta() ); // s

  Vec<double> energies;
  Vec<XYZ> positions;
  Vec<XYZ> momenta;
  energies.reserve(mNsteps);
  positions.reserve(mNsteps);
  momenta.reserve(mNsteps);

  unsigned nStepsUsed = 0;
  while(nStepsUsed<mNsteps)
    {
      positions.push_back( partPos);
      momenta.push_back( partMom );

      // direction to move without magnetic field in cm
      XYZ posIncr = partMom * ( mStepSize / TMath::Sqrt(partMom.Mag2()) );
      //XYZ posIncr = partMom * mStepSize;

      // new position after deltaT without magnetic field
      XYZ partPosNext = partPos + posIncr;
      // center of begin and stop vector without magnetic field
      XYZ partPosMid = (partPos + partPosNext) * 0.5;     

      XYZ Bfield = magnets.field(partPosMid, scale);

      if(Bfield.Mag2() == 0.0)
	  partPos = partPosNext;
      else
	{
	  // F = q*v X B
	  XYZ force = calc_lorentz_force(charge, partVel, Bfield); // (A*s)*(cm/s)*(kg/(A*s*s)) = (cm*kg)/(s*s)

	  // F = (dp/dt)
	  XYZ forceDelta = force * deltaT * 1.8708026E16; // (cm*kg)/(s*s) * delta_p (cm*kg)/s -> (GeV/c)
	  XYZ partMomNext = partMom + forceDelta;
	  
	  double mag0 = TMath::Sqrt(partMom.Mag2());
	  double mag1 = TMath::Sqrt(partMomNext.Mag2());
	  partMomNext *= mag0 / mag1; // make sure that total momentum doesn't change

	  XYZ momDelta = partMomNext * ( mStepSize / mag1); // direction to move with magnetic field in cm
	  //XYZ momDelta = partMomNext * mStepSize; // direction to move with magnetic field in cm
	  partPos += momDelta; // new position after deltaT with magnetic field

	  partMom = partMomNext;

	  // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) (cm/s)
	  Ltz tmpLorentz( partMom.X(), partMom.Y(), partMom.Z(), mParticle.mass );
	  partVel = calc_relativistic_velocity(partMom, tmpLorentz.Gamma(), mParticle.mass);
	}

      energies.push_back( TMath::Sqrt(partMom.Mag2() + 0.938*0.938) );
      
      ++nStepsUsed;
      
      if(fabs(partPos.Z()) > 9000.0) break;
    }

  return Track(nStepsUsed, energies, positions, momenta);
}

Track SimParticle::track_rungekutta4(const Magnets& magnets, double scale)
{
  double charge = mParticle.charge * mEcharge; // C = A*s

  Ltz initLorentzVec( mParticle.mom.X(), mParticle.mom.Y(), mParticle.mom.Z(), mParticle.mass );
  double deltaT = mStepSize / ( mSpeedOfLight * initLorentzVec.Beta() ); // s
  XYZ partPos = mParticle.pos; //cm
  XYZ partMom = initLorentzVec.Vect(); // Gev/c
  // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) -> (cm/s)

  XYZ partVel = calc_relativistic_velocity(partMom, initLorentzVec.Gamma(), mParticle.mass);

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
      //XYZ posIncr = partMom * ( mStepSize / TMath::Sqrt(partMom.Mag2()) );
      XYZ posIncr = partVel * mStepSize / TMath::Sqrt(partVel.Mag2());

      // new position after delta_t without magnetic field
      XYZ partPosNext = partPos + posIncr;
      // center of begin and stop vector without magnetic field
      XYZ partPosMid = (partPos + partPosNext) * 0.5;

      XYZ Bfield = magnets.field(partPosMid, scale);

      if(Bfield.Mag2() == 0.0)
	  partPos = partPosNext;
      else
	{
	  // F = q*v X B (runge-kutta k1)
	  XYZ k1_p = calc_lorentz_force(charge, partVel, Bfield); // (A*s)*(cm/s)*(kg/(A*s*s)) = (cm*kg)/(s*s))

	  //runge-kutta k2 term
	  XYZ p2 = partMom + mStepSize * k1_p * 0.5;
	  Ltz ltz2( p2.X(), p2.Y(), p2.Z(), mParticle.mass );
	  XYZ vel2 = calc_relativistic_velocity(p2, ltz2.Gamma(), mParticle.mass);
	  XYZ k2_p = calc_lorentz_force(charge, vel2, Bfield); // (A*s)*(cm/s)*(kg/(A*s*s)) = (cm*kg)/(s*s)

	  //runge-kutta k3 term
	  XYZ p3 = partMom + mStepSize * k2_p * 0.5;
	  Ltz ltz3( p3.X(), p3.Y(), p3.Z(), mParticle.mass );
	  XYZ vel3 = calc_relativistic_velocity(p3, ltz3.Gamma(), mParticle.mass);
	  XYZ k3_p = calc_lorentz_force(charge, vel3, Bfield); // (A*s)*(cm/s)*(kg/(A*s*s)) = (cm*kg)/(s*s)

	  //runge-kutta k4 term
	  XYZ p4 = partMom + mStepSize * k3_p;
	  Ltz ltz4( p4.X(), p4.Y(), p4.Z(), mParticle.mass );
	  XYZ vel4 = calc_relativistic_velocity(p4, ltz4.Gamma(), mParticle.mass);
	  XYZ k4_p = calc_lorentz_force(charge, vel4, Bfield); // (A*s)*(cm/s)*(kg/(A*s*s)) = (cm*kg)/(s*s)

	  // (cm*kg)/(s*s) * delta_p (cm*kg)/s -> (GeV/c)
	  XYZ forceDelta = ( deltaT * 1.8708026E16 * 0.16666666 * (k1_p + 2*k2_p + 2*k3_p + k4_p) );
	  XYZ partMomNext = partMom + forceDelta;

	  // std::cerr << nStepsUsed << std::endl;
	  // std::cerr << "K1: " << k1_p.X() << ", " << k1_p.Y() << ", " << k1_p.Z() << std::endl;
	  // std::cerr << "K2: " << k2_p.X() << ", " << k2_p.Y() << ", " << k2_p.Z() << std::endl;
	  // std::cerr << "K3: " << k3_p.X() << ", " << k3_p.Y() << ", " << k3_p.Z() << std::endl;
	  // std::cerr << "K4: " << k4_p.X() << ", " << k4_p.Y() << ", " << k4_p.Z() << std::endl;

	  // std::cerr << "P2: " << p2.X() << ", " << p2.Y() << ", " << p2.Z() << std::endl;
	  // std::cerr << "P3: " << p3.X() << ", " << p3.Y() << ", " << p3.Z() << std::endl;
	  // std::cerr << "P4: " << p4.X() << ", " << p4.Y() << ", " << p4.Z() << std::endl;

	  // std::cerr << "VEL2: " << vel2.X() << ", " << vel2.Y() << ", " << vel2.Z() << std::endl;
	  // std::cerr << "VEL3: " << vel3.X() << ", " << vel3.Y() << ", " << vel3.Z() << std::endl;
	  // std::cerr << "VEL4: " << vel4.X() << ", " << vel4.Y() << ", " << vel4.Z() << std::endl;
	  
	  // std::cerr << "Mom before: " << partMom.X() << ", " << partMom.Y() << ", " << partMom.Z() << std::endl;
	  // std::cerr << "Mom after: " << partMomNext.X() << ", " << partMomNext.Y() << ", " << partMomNext.Z() << std::endl;
	  // std::cerr << "Bfield: " << Bfield.X() << ", " << Bfield.Y() << ", " << Bfield.Z() << std::endl;

	  // v = dx/dt = p * c / (gamma * m) (runge-kutta k1 term)
	  XYZ partPosNext = partPos + ( mStepSize * 0.16666666 * (partVel + 2*vel2 + 2*vel3 + vel4) / TMath::Sqrt(partVel.Mag2()) );
	  partPos = partPosNext;

	  // momentum normalization
	  double mag2 = TMath::Sqrt(partMom.Mag2());
	  double mag3 = TMath::Sqrt(partMomNext.Mag2());
	  partMomNext *= mag2/mag3;
	  partMom = partMomNext;
	  
	  Ltz ltz_for_next_vel( partMom.X(), partMom.Y(), partMom.Z(), mParticle.mass );
	  partVel = calc_relativistic_velocity(partMom, ltz_for_next_vel.Gamma(), mParticle.mass);

	  // std::cout << "Positions: " << nStepsUsed << ", " << partPos.X() << ", " << partPos.Y() << ", " << partPos.Z() << std::endl;
	  // std::cout << "Velocities: " << nStepsUsed << ", " << partVel.X() << ", " << partVel.Y() << ", " << partVel.Z() << std::endl;
	  // std::cout << std::endl;

	}

      energies.push_back( TMath::Sqrt(partMom.Mag2() + 0.938*0.938) );
      
      ++nStepsUsed;

      if(fabs(partPos.Z()) > 9000.0) break;
    }

  return Track(nStepsUsed, energies, positions, momenta);
}
