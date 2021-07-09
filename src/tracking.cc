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
  XYZ part_pos[4];
  XYZ part_dir[3];
  XYZ part_vel[3];
  XYZ part_mom[3];
  XYZ B_field;
  XYZ F_field;

  // part_pos[0].SetXYZ( mParticle.pos.X(), mParticle.pos.Y(), mParticle.pos.Z() );
  part_pos[0] = mParticle.pos;
  //  part_mom[0].SetXYZ( lvector.Px(), lvector.Py(), lvector.Pz() ); // (GeV/c)
  part_mom[0] = lvector.Vect(); // (GeV/c)
  part_dir[0] = part_mom[0]; // (GeV/c)
  part_vel[0] = part_dir[0]; // (GeV/c)
  part_vel[0] *= (mCvelocity/(lvector.Gamma()*mParticle.mass)); // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) -> (cm/s)

  //cout << "gamma = " << lvector.Gamma() << endl;
  //double original_phi = part_mom[0].Phi();

  //double init_momentum = part_mom[0].Mag();
  //double init_pz       = part_mom[0].Pz();

  //cout << "eta = " << eta << ", pT = " << pT << ", p = " << init_momentum << ", pz = " << init_pz << endl;

  double delta_t = mStepSize/(mCvelocity*lvector.Beta()); // s

  // 1.8708026E16; // delta_p  (cm*kg)/s -> (GeV/c)

  Vec<double> energies;
  Vec<XYZ> positions;
  Vec<XYZ> directions;
  energies.reserve(mNsteps);
  positions.reserve(mNsteps);
  directions.reserve(mNsteps);

  unsigned nStepsUsed = 0;
  while(nStepsUsed<mNsteps)
    {
      XYZ tmppos; tmppos.SetXYZ( part_pos[0].X(), part_pos[0].Y(), part_pos[0].Z() );
      positions.push_back( tmppos );
      XYZ tmpdir; tmpdir.SetXYZ( part_dir[0].X(), part_dir[0].Y(), part_dir[0].Z() );
      directions.push_back( tmpdir );

      //printf("istep: %d, pos: {%4.8f, %4.8f, %4.8f} \n",istep,track_pos_dir[istep][0],track_pos_dir[istep][1],track_pos_dir[istep][2]);

      // p = beta*gamma*m0
      // beta = (ds/dt)/c

      //cout << "" << endl;
      //cout << "-----------------------------------------------------------------------------------" << endl;
      //cout << "istep = " << istep << ", old velocity = (" << part_vel[0].Px() << ", " << part_vel[0].Py() << ", " << part_vel[0].Pz()
      //    << "), |v| = " << part_vel[0].Mag() << ", pos = ("
      //    << part_pos[0].X() << ", " << part_pos[0].Y() << ", " << part_pos[0].Z() << ")" << endl;


      //cout << "delta_t = " << delta_t << endl;
      part_dir[0] *= (mStepSize/TMath::Sqrt(part_dir[0].Mag2())); // direction to move without magnetic field in cm
      //cout << "part_dir[0] = (" << part_dir[0].X() << ", " << part_dir[0].Y() << ", " << part_dir[0].Z() << ")" << endl;
      part_pos[1] =  part_pos[0];
      part_pos[1] += part_dir[0]; // new position after delta_t without magnetic field
      part_pos[2] =  part_pos[0];
      part_pos[2] += part_pos[1];
      part_pos[2] *= 0.5; // center of begin and stop vector without magnetic field
      //cout << "Center of MF, part_pos[2] = (" << part_pos[2].X() << ", " << part_pos[2].Y() << ", " << part_pos[2].Z() << ")" << endl;


      B_field = magnets.field(part_pos[2], 350.);
      //printf("B_field: {%f,%f,%f} \n",B_field.X(),B_field.Y(),B_field.Z());

      double Btot = B_field.Mag2();

      if(Btot == 0.0)
	  part_pos[0] = part_pos[1];
      else
	{
	  // F = q*v X B
	  F_field = charge*part_vel[0].Cross(B_field); // (A*s)*(cm/s)*(kg/(A*s*s)) = (cm*kg)/(s*s)

	  // F = (dp/dt)
	  part_dir[1] =  F_field; // (cm*kg)/(s*s)
	  part_dir[1] *= delta_t*1.8708026E16; // delta_p  (cm*kg)/s -> (GeV/c)
	  //cout << "delta_p, part_dir[1] = (" << part_dir[1].X() << ", " << part_dir[1].Y() << ", " << part_dir[1].Z() << "), step_size = " << step_size << endl;
	  part_mom[1] =  part_mom[0];
	  part_mom[1] += part_dir[1];
	  part_mom[1] *= TMath::Sqrt(part_mom[0].Mag2()) / TMath::Sqrt(part_mom[1].Mag2()); // make sure that total momentum doesn't change

	  part_dir[2] =  part_mom[1];
	  part_dir[2] *= (mStepSize/TMath::Sqrt(part_mom[1].Mag2())); // direction to move with magnetic field in cm
	  //cout << "MF dir vector, part_dir[2] = (" << part_dir[2].X() << ", " << part_dir[2].Y() << ", " << part_dir[2].Z() << "), step_size = " << step_size << endl;
	  part_pos[3] =  part_pos[0];
	  part_pos[3] += part_dir[2]; // new position after delta_t with magnetic field

	  part_pos[0] = part_pos[3];
	  part_mom[0] = part_mom[1];
	  part_dir[0] = part_mom[0]; // (GeV/c)
	  part_vel[0] = part_dir[0]; // (GeV/c)
	  part_vel[0] *= (mCvelocity/(lvector.Gamma()*mParticle.mass)); // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) (cm/s)
	}

      energies.push_back( TMath::Sqrt(part_mom[0].Mag2() + 0.938*0.938) );
      //cout << "istep = " << istep << ", z = " << part_pos[0].Z() << ", phi = " << part_mom[0].Phi() << endl;


      ++nStepsUsed;
      
      //if(part_pos[0].Perp() > 200.0) break;
      if(fabs(part_pos[0].Z()) > 9000.0) break;
      //cout << "Total momentum = " << part_mom[0].Mag() << endl;
      //cout << "New momentum = (" << part_mom[0].Px() << ", " << part_mom[0].Py() << ", " << part_mom[0].Pz() << "), pos = ("
      //    << part_pos[0].X() << ", " << part_pos[0].Y() << ", " << part_pos[0].Z() << ")" << endl;
    }

  return Track(nStepsUsed, energies, positions, directions);
}
