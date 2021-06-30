#ifndef FUNCTIONS_PROTON_SPEC_H
#define FUNCTIONS_PROTON_SPEC_H


static const Double_t part_mass   = 0.938; // GeV/c2
static const Double_t delta_s     = 1.25; // step size for tracking in cm, usually 2.5
static const Int_t N_total_steps = 3000; // number of steps for particle tracking, usually 300
static const Int_t step_size     = 10; // step size in cm
static const Double_t c_light     = 29979245800.0; // (cm/s)
static const Double_t e_charge    = 1.602176565E-19; // C = A*s
static Double_t track_pos_dir[N_total_steps][6]; // [step][x,y,z,px,py,pz]

static const Int_t N_magnets = 7;
static const Double_t Magnet_pos_z[2][N_magnets] =
{
    {750.0,1920.0,2300.0,3180.0,3830.0,4730.0,5840.0},
    {750.0+430.0,1920.0+190.0,2300.0+630.0,3180.0+550.0,3830.0+550.0,4730.0+630.0,5840.0+945.0}
};

static const Int_t N_rabbit = 4;
static const Double_t rabbit_pos_z[4][N_rabbit] =  // xstart, xstop, ystart, ystop
{
    {10.0,10.0,10.0,10.0},
    {35.0,35.0,35.0,35.0},
    {11500.0-50.0,11500.0-80.0,11500.0-100.0,11500.0-150.0},
    {11500.0-40.0,11500.0-70.0,11500.0-90.0,11500.0-140.0}
};

static const TString Label_magnets[N_magnets] = {"D(#mu)","D(c)","Q1","Q2","Q3","Q4","D1"};


// Scanned Fermi momentum distribution in Pb from https://arxiv.org/pdf/1405.0583.pdf
static const Double_t pT_fermi[38]   = {0.0125,0.025,0.0375,0.0527778,0.0611111,0.0708333,0.0805556,0.0902778,0.102778,0.115278,0.125,0.136111,0.156944,0.176389,0.194444,0.211111,0.229167,0.248611,0.265278,0.281944,0.295833,0.309722,0.323611,0.336111,0.347222,0.359722,0.376389,0.390278,0.406944,0.425,0.454167,0.483333,0.516667,0.552778,0.583333,0.613889,0.638889,0.648611};
static const Double_t prob_fermi[38] = {0.043824,0.346614,0.665339,1.06375,1.43028,1.81275,2.0996,2.40239,2.64143,2.88048,2.96016,3.00797,2.96016,2.83267,2.70518,2.59363,2.49801,2.37052,2.29084,2.14741,2.01992,1.89243,1.73307,1.62151,1.49402,1.41434,1.30279,1.20717,1.11155,1.06375,0.936255,0.87251,0.792829,0.713147,0.649402,0.585657,0.553785};



//------------------------------------------------------------------------------------------------------------------
void Get_XY_hit_points(Double_t det_z_pos, Int_t N_steps_used, Double_t &posx, Double_t &posy,
                       Double_t &alpha, Double_t &beta, Double_t &gamma)
{
    posx   = -1000.0;
    posy   = -1000.0;
    alpha  = -1000.0;
    beta   = -1000.0;
    for(Int_t i_step = 1; i_step < N_steps_used; i_step++) // loop over all track step positions
    {
        //printf("i_pz: %d, i_px: %d, pos: {%f,%f,%f} \n",i_pz,i_px,track_pos_dir[i_step][0],track_pos_dir[i_step][1],track_pos_dir[i_step][2]);
        Double_t x_trackA = track_pos_dir[i_step][0];
        Double_t y_trackA = track_pos_dir[i_step][1];
        Double_t z_trackA = track_pos_dir[i_step][2];

        Double_t x_trackB = track_pos_dir[i_step-1][0];
        Double_t y_trackB = track_pos_dir[i_step-1][1];
        Double_t z_trackB = track_pos_dir[i_step-1][2];

        if(fabs(z_trackA) >= fabs(det_z_pos))
        {
            posx   = x_trackA;
            posy   = y_trackA;
            alpha  = TMath::ATan2(z_trackB-z_trackA,x_trackB-x_trackA)*TMath::RadToDeg();
            beta   = TMath::ATan2(z_trackB-z_trackA,y_trackB-y_trackA)*TMath::RadToDeg();
            TVector3 vec_particle, vec_plane;
            vec_particle.SetXYZ(x_trackA-x_trackB,y_trackA-y_trackB,z_trackA-z_trackB);
            vec_plane.SetXYZ(0.0,0.0,1.0);
            gamma = vec_particle.Angle(vec_plane)*TMath::RadToDeg();
            //printf("posx in: %f, alphax: %f \n",posx,alphax);
            break;
        }
    }
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
TH1D* Get_fermi_distribution()
{
    TGraph* tg_fermi_prob = new TGraph();
    TH1D* h_fermi_prob = new TH1D("h_fermi_prob","h_fermi_prob",100,0,0.65);
    tg_fermi_prob ->Set(38);
    for(Int_t i_p = 0; i_p  < 38; i_p++)
    {
        tg_fermi_prob ->SetPoint(i_p,pT_fermi[i_p],prob_fermi[i_p]);
    }

    for(Int_t i_bin = 1; i_bin <= h_fermi_prob->GetNbinsX(); i_bin++)
    {
        Double_t pt_val = h_fermi_prob ->GetBinCenter(i_bin);
        Double_t value  = tg_fermi_prob->Eval(pt_val);
        if(value < 0.0) value = 0.0;
        h_fermi_prob    ->SetBinContent(i_bin,value);
    }
    return h_fermi_prob;
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
TVector3 Ali_forward_B_field(TVector3 pos)
{
    // pos: space vector in cm

    TVector3 B_vec;
    B_vec.SetXYZ(0.0,0.0,0.0);
    Double_t scaling_factor  = 1.0; // correct value = 1.0
    Double_t field_direction = 1.0;
    Double_t sign_dir = 1.0;

    if(pos.Z() < 0.0) // Muon dipole and corrector is only on C-side (= negative z-direction)
    {
        sign_dir = -1.0;
        if(pos.Z() < -750.0 && pos.Z() > (-750.0 - 430.0)) // Muon dipole
        {
            //B_vec.SetXYZ(0.75*field_direction*scaling_factor,0.0,0.0);
            B_vec.SetXYZ(0.67*field_direction*scaling_factor,0.0,0.0);
        }
        if(pos.Z() < -1920.0 && pos.Z() > (-1920.0 - 190.0)) // Dipole corrector
        {
            B_vec.SetXYZ(-1.1716*field_direction*scaling_factor,0.0,0.0);
        }


        if(pos.Z() < -2300.0 && pos.Z() > (-2300.0 - 630.0))  // First quadrupole
        {
            B_vec.SetXYZ(200.34*field_direction*scaling_factor*pos.Y()*sign_dir/100.0,-200.34*field_direction*scaling_factor*pos.X()*sign_dir/100.0,0.0);
        }

        if(pos.Z() < -3180.0 && pos.Z() > (-3180.0 - 550.0))  // Second quadrupole
        {
            B_vec.SetXYZ(-200.34*field_direction*scaling_factor*pos.Y()*sign_dir/100.0,200.34*field_direction*scaling_factor*pos.X()*sign_dir/100.0,0.0);
        }
        if(pos.Z() < -3830.0 && pos.Z() > (-3830.0 - 550.0))  // Third quadrupole
        {
            B_vec.SetXYZ(-200.34*field_direction*scaling_factor*pos.Y()*sign_dir/100.0,200.34*field_direction*scaling_factor*pos.X()*sign_dir/100.0,0.0);
        }
        if(pos.Z() < -4730.0 && pos.Z() > (-4730.0 - 630.0))  // Fourth quadrupole
        {
            B_vec.SetXYZ(200.34*field_direction*scaling_factor*pos.Y()*sign_dir/100.0,-200.34*field_direction*scaling_factor*pos.X()*sign_dir/100.0,0.0);
        }

        if(pos.Z() < -5840.0 && pos.Z() > (-5840.0 - 945.0)) // D1
        {
            B_vec.SetXYZ(0.0,-3.529*scaling_factor*sign_dir,0.0);
        }
    }
    else
    {
        sign_dir = +1.0;
        if(pos.Z() > 2300.0 && pos.Z() < (2300.0 + 630.0))  // First quadrupole
        {
            B_vec.SetXYZ(200.34*field_direction*scaling_factor*pos.Y()*sign_dir/100.0,-200.34*field_direction*scaling_factor*pos.X()*sign_dir/100.0,0.0);
        }

        if(pos.Z() > 3180.0 && pos.Z() < (3180.0 + 550.0))  // Second quadrupole
        {
            B_vec.SetXYZ(-200.34*field_direction*scaling_factor*pos.Y()*sign_dir/100.0,200.34*field_direction*scaling_factor*pos.X()*sign_dir/100.0,0.0);
        }
        if(pos.Z() > 3830.0 && pos.Z() < (3830.0 + 550.0))  // Third quadrupole
        {
            B_vec.SetXYZ(-200.34*field_direction*scaling_factor*pos.Y()*sign_dir/100.0,200.34*field_direction*scaling_factor*pos.X()*sign_dir/100.0,0.0);
        }
        if(pos.Z() > 4730.0 && pos.Z() < (4730.0 + 630.0))  // Fourth quadrupole
        {
            B_vec.SetXYZ(200.34*field_direction*scaling_factor*pos.Y()*sign_dir/100.0,-200.34*field_direction*scaling_factor*pos.X()*sign_dir/100.0,0.0);
        }

        if(pos.Z() > 5840.0 && pos.Z() < (5840.0 + 945.0)) // D1
        {
            B_vec.SetXYZ(0.0,-3.529*scaling_factor*sign_dir,0.0);
        }
    }

    return B_vec;
}
//------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------
void Track_eta_pT_phi_q_xyz_m(Double_t px, Double_t py, Double_t pz, Int_t q, Double_t pos_x, Double_t pos_y, Double_t pos_z, Double_t mass,
                                                    Int_t N_steps, Double_t step_size, Int_t &N_steps_used)
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

    Double_t charge = 1.0;
    if(q == 1) charge = -1.0;

    charge *= e_charge; // C = A*s

    TLorentzVector lvector;
    TVector3 part_pos[4];
    TVector3 part_dir[3];
    TVector3 part_vel[3];
    TVector3 part_mom[3];
    TVector3 B_field;
    TVector3 F_field;


    lvector.SetXYZM(px,py,pz,mass);
    part_pos[0].SetXYZ(pos_x,pos_y,pos_z);
    part_mom[0].SetXYZ(lvector.Px(),lvector.Py(),lvector.Pz()); // (GeV/c)
    part_dir[0] = part_mom[0]; // (GeV/c)
    part_vel[0] = part_dir[0]; // (GeV/c)
    part_vel[0] *= (c_light/(lvector.Gamma()*mass)); // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) -> (cm/s)

    //cout << "gamma = " << lvector.Gamma() << endl;
    //Double_t original_phi = part_mom[0].Phi();

    //Double_t init_momentum = part_mom[0].Mag();
    //Double_t init_pz       = part_mom[0].Pz();

    //cout << "eta = " << eta << ", pT = " << pT << ", p = " << init_momentum << ", pz = " << init_pz << endl;

    Double_t delta_t = step_size/(c_light*lvector.Beta()); // s
    Double_t alpha   = (part_mom[0].Mag()/part_vel[0].Mag())/1.8708026E16; // (GeV/c)/(cm/s) -> ((cm*kg)/s)/(cm/s) = kg

    // 1.8708026E16; // delta_p  (cm*kg)/s -> (GeV/c)

    for(Int_t istep = 0; istep < N_steps; istep++)
    {
        track_pos_dir[istep][0] = part_pos[0].X();
        track_pos_dir[istep][1] = part_pos[0].Y();
        track_pos_dir[istep][2] = part_pos[0].Z();
        track_pos_dir[istep][3] = part_dir[0].X();
        track_pos_dir[istep][4] = part_dir[0].Y();
        track_pos_dir[istep][5] = part_dir[0].Z();

        //printf("istep: %d, pos: {%4.8f, %4.8f, %4.8f} \n",istep,track_pos_dir[istep][0],track_pos_dir[istep][1],track_pos_dir[istep][2]);

        // p = beta*gamma*m0
        // beta = (ds/dt)/c

        //cout << "" << endl;
        //cout << "-----------------------------------------------------------------------------------" << endl;
        //cout << "istep = " << istep << ", old velocity = (" << part_vel[0].Px() << ", " << part_vel[0].Py() << ", " << part_vel[0].Pz()
        //    << "), |v| = " << part_vel[0].Mag() << ", pos = ("
        //    << part_pos[0].X() << ", " << part_pos[0].Y() << ", " << part_pos[0].Z() << ")" << endl;


        //cout << "delta_t = " << delta_t << endl;
        part_dir[0] *= (step_size/part_dir[0].Mag()); // direction to move without magnetic field in cm
        //cout << "part_dir[0] = (" << part_dir[0].X() << ", " << part_dir[0].Y() << ", " << part_dir[0].Z() << ")" << endl;
        part_pos[1] =  part_pos[0];
        part_pos[1] += part_dir[0]; // new position after delta_t without magnetic field
        part_pos[2] =  part_pos[0];
        part_pos[2] += part_pos[1];
        part_pos[2] *= 0.5; // center of begin and stop vector without magnetic field
        //cout << "Center of MF, part_pos[2] = (" << part_pos[2].X() << ", " << part_pos[2].Y() << ", " << part_pos[2].Z() << ")" << endl;


        Double_t PosXYZ[3] = {part_pos[2].X(),part_pos[2].Y(),part_pos[2].Z()};
        Double_t BXYZ[3];
        //StarMag->Interpolate3DBfield(TMath::Sqrt(TMath::Power(part_pos[2].X(),2)+TMath::Power(part_pos[2].Y(),2)),part_pos[2].Z(),part_pos[2].Phi(),Br,Bz,Bphi); // kGauss, cm
        //StarMag->BField(PosXYZ,BXYZ); // <-
        //B_field.SetMagThetaPhi(TMath::Sqrt(TMath::Power(Br,2)+TMath::Power(Bz,2)),TMath::ATan2(Br,Bz),Bphi); // in kGauss
        //B_field.SetXYZ(BXYZ[0],BXYZ[1],BXYZ[2]); // <-
        //B_field *= 0.1; // kg/(A*s*s) // <-

        B_field = Ali_forward_B_field(part_pos[2]);

        //printf("B_field: {%f,%f,%f} \n",B_field.X(),B_field.Y(),B_field.Z());


        Double_t Btot   = TMath::Power(B_field.Mag(),2);

        if(Btot == 0.0)
        {
            part_pos[0] = part_pos[1];
            //printf("istep: %d, pos: {%f,%f,%f} \n",istep,part_pos[0].X(),part_pos[0].Y(),part_pos[0].Z());
        }

        if(Btot != 0.0)
        {
            Double_t q_c    = charge;
            Double_t SqBtot = TMath::Sqrt(Btot);
            Double_t SinB   = TMath::Sin(SqBtot*q_c*delta_t/alpha);
            Double_t CosB   = TMath::Cos(SqBtot*q_c*delta_t/alpha);

            Double_t Bx  = B_field.X();
            Double_t By  = B_field.Y();
            Double_t Bz  = B_field.Z();

            Double_t vx0 = part_vel[0].X();
            Double_t vy0 = part_vel[0].Y();
            Double_t vz0 = part_vel[0].Z();

            Double_t term1  = (Bx*vx0 + By*vy0 + Bz*vz0);

            part_vel[1].SetX((Bx*term1 + (By*By*vx0 - Bx*By*vy0 + Bz*(Bz*vx0 - Bx*vz0))*CosB + SqBtot*(Bz*vy0 - By*vz0)*SinB)/Btot);   // (cm/s)
            part_vel[1].SetY((By*term1 + (-Bx*By*vx0 + Bx*Bx*vy0 + Bz*(Bz*vy0 - By*vz0))*CosB + SqBtot*(-Bz*vx0 + Bx*vz0)*SinB)/Btot); // (cm/s)
            part_vel[1].SetZ((Bz*term1 + (-Bz*(Bx*vx0 + By*vy0) + (Bx*Bx + By*By)*vz0)*CosB + SqBtot*(By*vx0 - Bx*vy0)*SinB)/Btot);    // (cm/s)

            //printf("vel: {%f,%f,%f} \n",part_vel[1].X(),part_vel[1].Y(),part_vel[1].Z());

            //part_pos[3].SetX((Bx*delta_t*term1 + (alpha*(-Bz*vy0 + By*vz0)*CosB/q) + (alpha*(By*By*vx0 + Bz*Bz*vx0 - Bx*By*vy0 - Bx*Bz*vz0)*SinB/(SqBtot*q)))/Btot);   // cm
            //part_pos[3].SetY((By*delta_t*term1 + (alpha*(-Bz*vx0 + Bx*vz0)*CosB/q) + (alpha*(-Bx*By*vx0 + Bx*Bx*vy0 + Bz*Bz*vy0 - By*Bz*vz0)*SinB/(SqBtot*q)))/Btot);  // cm
            //part_pos[3].SetZ((Bz*delta_t*term1 + (alpha*(-By*vx0 + Bx*vy0)*CosB/q) + (alpha*(-Bx*Bz*vx0 - By*Bz*vy0 + Bx*Bx*vz0 + By*By*vz0)*SinB/(SqBtot*q)))/Btot);  // cm

            part_pos[3].SetX(((Bx*Bx*delta_t*vx0 + Bx*By*delta_t*vy0 + Bx*Bz*delta_t*vz0) + (((Bz*vy0 - By*vz0)*(alpha-alpha*CosB))/q_c) + ((alpha*(By*By*vx0 - Bx*By*vy0 + Bz*(Bz*vx0 - Bx*vz0))*SinB)/(SqBtot*q_c)))/Btot);
            part_pos[3].SetY(((Bx*By*delta_t*vx0 + By*By*delta_t*vy0 + By*Bz*delta_t*vz0) + (((Bz*vx0 - Bx*vz0)*(alpha-alpha*CosB))/q_c) + ((alpha*(-Bx*By*vx0 + Bx*Bx*vy0 + Bz*(Bz*vy0 - By*vz0))*SinB)/(SqBtot*q_c)))/Btot);
            part_pos[3].SetZ(((Bx*Bz*delta_t*vx0 + By*Bz*delta_t*vy0 + Bz*Bz*delta_t*vz0) + (((By*vx0 - Bx*vy0)*(alpha-alpha*CosB))/q_c) + ((alpha*(-Bz*(Bx*vx0 + By*vy0) + (Bx*Bx + By*By)*vz0)*SinB)/(SqBtot*q_c)))/Btot);

            //printf("term1: %f, delta_t: %f, alpha: %f, pos: {%f,%f,%f} \n",term1,delta_t,alpha,part_pos[3].X(),part_pos[3].Y(),part_pos[3].Z());

            part_pos[2] =  part_pos[0];
            part_pos[2] += part_pos[3];
            //cout << "delta_t = " << delta_t << ", vz0 = " << vz0 << ", pos_z_old = " << part_pos[0].Z() << ", pos_z_new = " << part_pos[2].Z() << endl;
            part_pos[0] = part_pos[2];
            part_vel[0] = part_vel[1];
            part_dir[0] = part_vel[0]; // new

            printf("istep: %d, pos: {%4.8f, %4.8f, %4.8f}, B-field: {%4.5f, %4.5f, %4.5f}  \n",istep,part_pos[0].X(),part_pos[0].Y(),part_pos[0].Z(),Bx,By,Bz);

        }

        if(0)
        {

            //cout << "B_field = (" << B_field.X() << ", " << B_field.Y() << ", " << B_field.Z() << ")" << endl;
            // F = q*v X B
            F_field = charge*part_vel[0].Cross(B_field); // (A*s)*(cm/s)*(kg/(A*s*s)) = (cm*kg)/(s*s)

            // F = (dp/dt)
            part_dir[1] =  F_field; // (cm*kg)/(s*s)
            part_dir[1] *= delta_t*1.8708026E16; // delta_p  (cm*kg)/s -> (GeV/c)
            //cout << "delta_p, part_dir[1] = (" << part_dir[1].X() << ", " << part_dir[1].Y() << ", " << part_dir[1].Z() << "), step_size = " << step_size << endl;
            part_mom[1] =  part_mom[0];
            part_mom[1] += part_dir[1];
            part_mom[1] *= part_mom[0].Mag()/part_mom[1].Mag(); // make sure that total momentum doesn't change

            part_dir[2] =  part_mom[1];
            part_dir[2] *= (step_size/part_mom[1].Mag()); // direction to move with magnetic field in cm
            //cout << "MF dir vector, part_dir[2] = (" << part_dir[2].X() << ", " << part_dir[2].Y() << ", " << part_dir[2].Z() << "), step_size = " << step_size << endl;
            part_pos[3] =  part_pos[0];
            part_pos[3] += part_dir[2]; // new position after delta_t with magnetic field

            part_pos[0] = part_pos[3];
            part_mom[0] = part_mom[1];
            part_dir[0] = part_mom[0]; // (GeV/c)
            part_vel[0] = part_dir[0]; // (GeV/c)
            part_vel[0] *= (c_light/(lvector.Gamma()*part_mass)); // p (GeV/c) = beta (c) *gamma*m0 (GeV/c2) (cm/s)

            //cout << "istep = " << istep << ", z = " << part_pos[0].Z() << ", phi = " << part_mom[0].Phi() << endl;
        }


        N_steps_used = istep;
        //if(part_pos[0].Perp() > 200.0) break;
        if(fabs(part_pos[0].Z()) > 7000.0) break;
        //cout << "Total momentum = " << part_mom[0].Mag() << endl;
        //cout << "New momentum = (" << part_mom[0].Px() << ", " << part_mom[0].Py() << ", " << part_mom[0].Pz() << "), pos = ("
        //    << part_pos[0].X() << ", " << part_pos[0].Y() << ", " << part_pos[0].Z() << ")" << endl;
    }
}
//------------------------------------------------------------------------------------------------------------------

#endif // FUNCTIONS_PROTON_SPEC_H
