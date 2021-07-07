const Track& SimParticle::track(const Magnets& magnets, TrackMode mode ) {
    
  if(mode == TrackMode::Euler) {
      
    if(mTrackCheck[TrackMode::Euler] == false) {
      mTracks[TrackMode::Euler] = track_euler(magnets);
      mTrackCheck[TrackMode::Euler] = true;
    }
    const Track& t = mTracks[TrackMode::Euler];
      
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
    
  return t;
}

Track& SimParticle::track_euler(const Magnets& magnets)
{
  // Track a particle through the STAR magnetic field with
  // pseudo-rapidity eta
  // transverse momentum pT (GeV/c)
  // azimuthal angle phi (rad)
  // charge q (0 = +, 1 = -)
  // z-vertex position z (cm)
  // particle mass m (GeV/c^2)
  // number of tracking steps N_steps
  // step size step_size (cm)
  
  double charge = mParticle.charge * mEcharge; // C = A*s

  TLorentzVector lvector;
  XYZ part_pos[4];
  XYZ part_dir[3];
  XYZ part_vel[3];
  XYZ part_mom[3];
  XYZ B_field;
  XYZ F_field;


  lvector.SetXYZM( mParticle.mom.X(), mParticle.mom.Y(), mParticle.mom.Z(), mParticle.mass);
  part_pos[0].SetXYZ( mParticle.pos.X(), mParticle.pos.Y(), mParticle.pos.Z());
  part_mom[0].SetXYZ(lvector.Px(),lvector.Py(),lvector.Pz()); // (GeV/c)
  part_dir[0] = part_mom[0]; // (GeV/c)
  part_vel[0] = part_dir[0]; // (GeV/c)
  part_vel[0] *= (mCvelocity/(lvector.Gamma()*mParticle.mass)); // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) -> (cm/s)

  //cout << "gamma = " << lvector.Gamma() << endl;
  //double original_phi = part_mom[0].Phi();

  //double init_momentum = part_mom[0].Mag();
  //double init_pz       = part_mom[0].Pz();

  //cout << "eta = " << eta << ", pT = " << pT << ", p = " << init_momentum << ", pz = " << init_pz << endl;

  double delta_t = step_size/(mCvelocity*lvector.Beta()); // s
  double alpha   = ( TMath::Sqrt(part_mom[0].Mag2())/TMath::Sqrt(part_vel[0].Mag2()) ) / 1.8708026E16; // (GeV/c)/(cm/s) -> ((cm*kg)/s)/(cm/s) = kg

  // 1.8708026E16; // delta_p  (cm*kg)/s -> (GeV/c)

  Vec<double> energies;
  Vec<XYZ> positions;
  Vec<XYZ> directions;
  energies.reserve(mNsteps);
  positions.reserve(mNsteps);
  directions.reserve(mNsteps);
	
  for(int istep=0; istep<mNsteps; istep++)
    {
      XYZ tmppos; tmpos.SetXYZ( part_pos[0].X(), part_pos[0].Y(), part_pos[0].Z() );
      positions.push_back( tmppos );
      XYZ tmpdir; tmdir.SetXYZ( part_dir[0].X(), part_dir[0].Y(), part_dir[0].Z() );
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
      part_dir[0] *= (step_size/TMath::Sqrt(part_dir[0].Mag2())); // direction to move without magnetic field in cm
      //cout << "part_dir[0] = (" << part_dir[0].X() << ", " << part_dir[0].Y() << ", " << part_dir[0].Z() << ")" << endl;
      part_pos[1] =  part_pos[0];
      part_pos[1] += part_dir[0]; // new position after delta_t without magnetic field
      part_pos[2] =  part_pos[0];
      part_pos[2] += part_pos[1];
      part_pos[2] *= 0.5; // center of begin and stop vector without magnetic field
      //cout << "Center of MF, part_pos[2] = (" << part_pos[2].X() << ", " << part_pos[2].Y() << ", " << part_pos[2].Z() << ")" << endl;


      double PosXYZ[3] = {part_pos[2].X(),part_pos[2].Y(),part_pos[2].Z()};
      double BXYZ[3];
      //StarMag->Interpolate3DBfield(TMath::Sqrt(TMath::Power(part_pos[2].X(),2)+TMath::Power(part_pos[2].Y(),2)),part_pos[2].Z(),part_pos[2].Phi(),Br,Bz,Bphi); // kGauss, cm
      //StarMag->BField(PosXYZ,BXYZ); // <-
      //B_field.SetMagThetaPhi(TMath::Sqrt(TMath::Power(Br,2)+TMath::Power(Bz,2)),TMath::ATan2(Br,Bz),Bphi); // in kGauss
      //B_field.SetXYZ(BXYZ[0],BXYZ[1],BXYZ[2]); // <-
      //B_field *= 0.1; // kg/(A*s*s) // <-

      //B_field = Ali_forward_B_field(part_pos[2]);
      B_field = magnets.field(part_pos[2], 350.);

      //printf("B_field: {%f,%f,%f} \n",B_field.X(),B_field.Y(),B_field.Z());


      double Btot   = B_field.Mag2();

      if(Btot == 0.0)
        {
	  part_pos[0] = part_pos[1];
	  //printf("istep: %d, pos: {%f,%f,%f} \n",istep,part_pos[0].X(),part_pos[0].Y(),part_pos[0].Z());
        }

      //cout << "B_field = (" << B_field.X() << ", " << B_field.Y() << ", " << B_field.Z() << ")" << endl;
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
      part_dir[2] *= (step_size/TMath::Sqrt(part_mom[1].Mag2())); // direction to move with magnetic field in cm
      //cout << "MF dir vector, part_dir[2] = (" << part_dir[2].X() << ", " << part_dir[2].Y() << ", " << part_dir[2].Z() << "), step_size = " << step_size << endl;
      part_pos[3] =  part_pos[0];
      part_pos[3] += part_dir[2]; // new position after delta_t with magnetic field

      part_pos[0] = part_pos[3];
      part_mom[0] = part_mom[1];
      part_dir[0] = part_mom[0]; // (GeV/c)
      part_vel[0] = part_dir[0]; // (GeV/c)
      part_vel[0] *= (mCvelocity/(lvector.Gamma()*mParticle.mass)); // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) (cm/s)

      energies.push_back( TMath::Sqrt(part_mom[0].Mag2() + 0.938*0.938) );
      //cout << "istep = " << istep << ", z = " << part_pos[0].Z() << ", phi = " << part_mom[0].Phi() << endl;


      mNstepsUsed[TrackMode::Euler] = istep;
      //if(part_pos[0].Perp() > 200.0) break;
      if(fabs(part_pos[0].Z()) > 9000.0) break;
      //cout << "Total momentum = " << part_mom[0].Mag() << endl;
      //cout << "New momentum = (" << part_mom[0].Px() << ", " << part_mom[0].Py() << ", " << part_mom[0].Pz() << "), pos = ("
      //    << part_pos[0].X() << ", " << part_pos[0].Y() << ", " << part_pos[0].Z() << ")" << endl;
    }

  return Track(energies, positions, directions);
}
