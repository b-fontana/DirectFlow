#include "./functions.h"
#include "./geometry.h"
#include "./functions_proton_spec.h"

void v1_beam()
{
  gRandom->SetSeed(0);
  gStyle->SetPalette(56); // 53 = black body radiation, 56 = inverted black body radiator, 103 = sunset, 87 == light temperature

  ROOT::Math::XYZVector vec_p(0.0,0.0,1500.0); // 5000
  double pos_xyz[3] = {0.0,-0.0,-6900.0};

  TLorentzVector* TLV_particle = new TLorentzVector();

  double particle_energy = TMath::Sqrt(vec_p.Mag2() + 0.938*0.938);
  TLV_particle->SetPxPyPzE(vec_p.x(),vec_p.y(),vec_p.z(),particle_energy); // nucleon momentum in cm frame, Fermi momentum + v1 momentum

  // std::vector<Magnets::Magnet> magnetInfo{
  //    {Magnets::DipoleY,    "D1_neg", kBlue,    std::make_pair(0.,-3.529),  std::make_pair(-5840.0-945.0, -5840.0), 100., 100.},
  //    {Magnets::Quadrupole, "Q4_neg", kYellow,  std::make_pair(2.34,-2.34), std::make_pair(-4730.0-630.0, -4730.0), 100., 100.},
  //    {Magnets::Quadrupole, "Q3_neg", kYellow,  std::make_pair(-2.34,2.34), std::make_pair(-3830.0-550.0, -3830.0), 100., 100.},
  //    {Magnets::Quadrupole, "Q2_neg", kYellow,  std::make_pair(-2.34,2.34), std::make_pair(-3180.0-550.0, -3180.0), 100., 100.},
  //    {Magnets::Quadrupole, "Q1_neg", kYellow,  std::make_pair(2.34,-2.34), std::make_pair(-2300.0-630.0, -2300.0), 100., 100.},
  //    {Magnets::DipoleX,    "D_corr", kBlue+1,  std::make_pair(-1.1716,0.), std::make_pair(-1920.0-190.0, -1920.0), 100., 100.},
  //    {Magnets::DipoleX,    "Muon"  , kMagenta, std::make_pair(0.67,0.),    std::make_pair(-750.0-430.0,  -750.0),  100., 100.},
  //    {Magnets::Quadrupole, "Q1_pos", kYellow,  std::make_pair(2.34,-2.34), std::make_pair(2300.0, 2300.0+630.0),   100., 100.},
  //    {Magnets::Quadrupole, "Q2_pos", kYellow,  std::make_pair(-2.34,2.34), std::make_pair(3180.0, 3180.0+550.0),   100., 100.},
  //    {Magnets::Quadrupole, "Q3_pos", kYellow,  std::make_pair(-2.34,2.34), std::make_pair(3830.0, 3830.0+550.0),   100., 100.},
  //    {Magnets::Quadrupole, "Q4_pos", kYellow,  std::make_pair(2.34,-2.34), std::make_pair(4730.0, 4730.0+630.0),   100., 100.},
  //    {Magnets::DipoleY,    "D1_pos", kBlue,    std::make_pair(0.,-3.529), std::make_pair(5840.0, 5840.0+945.0),   100., 100.}
  // };

  // std::vector<Magnets::Magnet> magnetInfo{
  //    {Magnets::DipoleY,    "D1_neg", kGreen,    std::make_pair(0.,-3.529),  std::make_pair(-5840.0-945.0, -5840.0), 100., 100.},
  // };

  std::vector<Magnets::Magnet> magnetInfo{
     {Magnets::DipoleY,    "D1_neg", kGreen,    std::make_pair(0.,-3.529),  std::make_pair(pos_xyz[2]-1500, pos_xyz[2]+3000), 600., 600.},
  };

  Magnets magnets(magnetInfo);
  magnets.draw();
  
  vector<TEveLine*> vec_particle_track;
  vec_particle_track.resize(1);
  vec_particle_track[0] = new TEveLine();


  //-------------------------------------------------------------------------------------------
  // Scanned sigma x (m) as a function of z (m) of the LHC beam  (John Jowett)
  Double_t z_LHC[37]       = {-67.643,-58.738,-54.063,-50.500,-47.161,-44.267,-40.259,-38.701,-32.244,-30.241,-26.679,-23.784,-20.222,-17.105,-13.098,-9.3135,-4.6382,-0.4081,3.37662,8.05195,13.6178,17.6252,20.2968,22.7458,26.0853,28.9796,31.8738,34.7681,38.1076,40.7792,44.1187,46.3451,49.2393,53.0241,58.3673,63.4879,69.0538};
  Double_t sigma_x_LHC[37] = {0.0013124,0.0014062,0.0014562,0.0013687,0.0012687,0.0010499,0.0008562,0.0007999,0.0007562,0.0007812,0.0007624,0.0007437,0.0006249,0.0005249,0.0004062,0.0002812,0.0001437,4.37485e-05,0.0001312,0.0002687,0.0004374,0.0005624,0.0006562,0.0007187,0.0009124,0.0010812,0.0013124,0.0014312,0.0015624,0.0014874,0.0014249,0.0012374,0.0011374,0.0010187,0.0010062,0.0009812,0.0009687};
  //-------------------------------------------------------------------------------------------



  Int_t N_steps_used = 0;
  Track_eta_pT_phi_q_xyz_m(magnets,
                           TLV_particle->Px(),TLV_particle->Py(),TLV_particle->Pz(),
			   0,
			   pos_xyz[0],pos_xyz[1],pos_xyz[2],
			   0.938,
			   N_total_steps,
			   step_size,
			   N_steps_used); // fills track_pos_dir

  std::string filename("track_pos.csv");
  std::fstream file;
  file.open(filename, std::ios_base::out);
    
  for(Int_t i_step = 0; i_step < N_steps_used; i_step++)
    {
      vec_particle_track[0]        ->SetNextPoint(track_pos_dir[i_step][0],track_pos_dir[i_step][1],track_pos_dir[i_step][2]);

      if (!file.is_open()) 
	std::cout << "failed to open " << filename << '\n';
      else {
	if(i_step==0)
	  file << "x" << "," << "z" << std::endl;
	//if(i_step%100==0)
	file << std::to_string(track_pos_dir[i_step][0]) << "," << std::to_string(track_pos_dir[i_step][2]) << std::endl;
      }
      //f >> std::to_string(track_pos_dir[i_step][0]) >> std::endl;
      //printf("i_step: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_step,track_pos_dir[i_step][0],track_pos_dir[i_step][1],track_pos_dir[i_step][2]);

    }

  HistName = "track ";
  HistName += 0;
  std::cout << HistName.Data() << std::endl;
  std::exit(0);
  vec_particle_track[0] ->SetName(HistName.Data());
  vec_particle_track[0] ->SetLineStyle(1);
  vec_particle_track[0] ->SetLineWidth(5);
  vec_particle_track[0] ->SetMainColor(kRed);
  vec_particle_track[0] ->SetMainAlpha(0.7);
  gEve->AddElement(vec_particle_track[0]);
  gEve->Redraw3D(kTRUE);
}

