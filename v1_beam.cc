#include "include/geometry.h"
#include "include/tracking.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <boost/program_options.hpp>

#include <TApplication.h>
#include "Math/Vector3D.h" // XYZVector

#include <TEveManager.h>
#include "TEveLine.h"

#include "TStyle.h"
#include "TRandom.h"

void run(tracking::TrackMode mode)
{
  gRandom->SetSeed(0);
  gStyle->SetPalette(56); // 53 = black body radiation, 56 = inverted black body radiator, 103 = sunset, 87 == light temperature
 
  //define the initial properties of the incident particle
  Particle particle;
  particle.pos = ROOT::Math::XYZVector(0.0, 0.0, -7000.0); // cm
  particle.mom = ROOT::Math::XYZVector(0.0, 0.0, 3500.0); // GeV/c
  particle.mass = 0.938; // GeV/c^2
  particle.energy = TMath::Sqrt(particle.mom.Mag2() + particle.mass*particle.mass);
  particle.charge = +1;
  
  std::vector<Magnets::Magnet> magnetInfo{
     {Magnets::DipoleY,    "D1_neg", kBlue,    std::make_pair(0.,-3.529),  std::make_pair(-5840.0-945.0, -5840.0), 100., 100.},
     {Magnets::Quadrupole, "Q4_neg", kYellow,  std::make_pair(2.34,-2.34), std::make_pair(-4730.0-630.0, -4730.0), 100., 100.},
     {Magnets::Quadrupole, "Q3_neg", kYellow,  std::make_pair(-2.34,2.34), std::make_pair(-3830.0-550.0, -3830.0), 100., 100.},
     {Magnets::Quadrupole, "Q2_neg", kYellow,  std::make_pair(-2.34,2.34), std::make_pair(-3180.0-550.0, -3180.0), 100., 100.},
     {Magnets::Quadrupole, "Q1_neg", kYellow,  std::make_pair(2.34,-2.34), std::make_pair(-2300.0-630.0, -2300.0), 100., 100.},
     {Magnets::DipoleX,    "D_corr", kBlue+1,  std::make_pair(-1.1716,0.), std::make_pair(-1920.0-190.0, -1920.0), 100., 100.},
     {Magnets::DipoleX,    "Muon"  , kMagenta, std::make_pair(0.67,0.),    std::make_pair(-750.0-430.0,  -750.0),  100., 100.},
     {Magnets::Quadrupole, "Q1_pos", kYellow,  std::make_pair(2.34,-2.34), std::make_pair(2300.0, 2300.0+630.0),   100., 100.},
     {Magnets::Quadrupole, "Q2_pos", kYellow,  std::make_pair(-2.34,2.34), std::make_pair(3180.0, 3180.0+550.0),   100., 100.},
     {Magnets::Quadrupole, "Q3_pos", kYellow,  std::make_pair(-2.34,2.34), std::make_pair(3830.0, 3830.0+550.0),   100., 100.},
     {Magnets::Quadrupole, "Q4_pos", kYellow,  std::make_pair(2.34,-2.34), std::make_pair(4730.0, 4730.0+630.0),   100., 100.},
     {Magnets::DipoleY,    "D1_pos", kBlue,    std::make_pair(0.,-3.529), std::make_pair(5840.0, 5840.0+945.0),   100., 100.}
  };

  // std::vector<Magnets::Magnet> magnetInfo{ {Magnets::DipoleY,
  // 					    "D1_neg",
  // 					    kGreen,
  // 					    std::make_pair(0.,-3.529),
  // 					    std::make_pair(-5840.0-945.0, -5840.0),
  // 					    100., 100.},
  // };

  // std::vector<Magnets::Magnet> magnetInfo{ {Magnets::DipoleY,
  // 					    "D1_neg",
  // 					    kGreen,
  // 					    std::make_pair(0.,-0.5),
  // 					    std::make_pair(particle.pos.Z()-1500, particle.pos.Z()+3000), 600., 600.},
  // };

  Magnets magnets(magnetInfo);
  magnets.draw();
  std::vector<TEveLine*> particleTrackViz;
  particleTrackViz.resize(1);
  particleTrackViz[0] = new TEveLine();

  //-------------------------------------------------------------------------------------------
  // Scanned sigma x (m) as a function of z (m) of the LHC beam  (John Jowett)
  // Double_t z_LHC[37]       = {-67.643,-58.738,-54.063,-50.500,-47.161,-44.267,-40.259,-38.701,-32.244,-30.241,-26.679,-23.784,-20.222,-17.105,-13.098,-9.3135,-4.6382,-0.4081,3.37662,8.05195,13.6178,17.6252,20.2968,22.7458,26.0853,28.9796,31.8738,34.7681,38.1076,40.7792,44.1187,46.3451,49.2393,53.0241,58.3673,63.4879,69.0538};
  // Double_t sigma_x_LHC[37] = {0.0013124,0.0014062,0.0014562,0.0013687,0.0012687,0.0010499,0.0008562,0.0007999,0.0007562,0.0007812,0.0007624,0.0007437,0.0006249,0.0005249,0.0004062,0.0002812,0.0001437,4.37485e-05,0.0001312,0.0002687,0.0004374,0.0005624,0.0006562,0.0007187,0.0009124,0.0010812,0.0013124,0.0014312,0.0015624,0.0014874,0.0014249,0.0012374,0.0011374,0.0010187,0.0010062,0.0009812,0.0009687};
  //-------------------------------------------------------------------------------------------

  std::array<unsigned, tracking::TrackMode::NMODES> nsteps = {{ 18000, 18000 }};
  std::array<double, tracking::TrackMode::NMODES> stepsize = {{ 1., 1. }};
  std::array<std::string, tracking::TrackMode::NMODES> suf = {{ "_euler", "_rk4" }};

  double Bscale = 1.;
  SimParticle simp(particle, nsteps[mode], stepsize[mode]);

  const Track& track = simp.track(magnets, mode, Bscale );
  std::vector<double> itEnergies = track.energies();
  std::vector<ROOT::Math::XYZVector> itPositions = track.positions();
  std::vector<ROOT::Math::XYZVector> itMomenta = track.momenta();
  unsigned nStepsUsed = track.steps_used();
    
  std::string filename("data/track" + suf[mode] + ".csv");
  std::fstream file;
  file.open(filename, std::ios_base::out);
  for(unsigned i_step = 0; i_step < nStepsUsed; i_step++)
    {
      particleTrackViz[0]->SetNextPoint( itPositions[i_step].X(), itPositions[i_step].Y(), itPositions[i_step].Z() );

      if (!file.is_open()) 
 	std::cerr << "failed to open " << filename << '\n';
      else {
	if(i_step==0)
	  file << "x,y,z,energy" << std::endl;
	file << std::to_string( itPositions[i_step].X() ) << ","
	     << std::to_string( itPositions[i_step].Y() ) << ","
	     << std::to_string( itPositions[i_step].Z() ) << ","
	     << std::to_string( itEnergies[i_step] )
	     << std::endl;
      }
    }

  std::string histname = "track 0";
  particleTrackViz[0]->SetName( histname.c_str() );
  particleTrackViz[0]->SetLineStyle(1);
  particleTrackViz[0]->SetLineWidth(5);
  particleTrackViz[0]->SetMainColor(kRed);
  particleTrackViz[0]->SetMainAlpha(0.7);
  gEve->AddElement(particleTrackViz[0]);
  gEve->Redraw3D(kTRUE);
}

int main(int argc, char **argv) {
  TApplication myapp("myapp", &argc, argv);

  tracking::TrackMode mode = tracking::TrackMode::Euler;
  if(argc<2)
    std::cerr << "No arguments passed." << std::endl;
  else
    {
      namespace po = boost::program_options;
      po::options_description desc("Options");
      desc.add_options()
	("mode", po::value<std::string>()->default_value("euler"), "numerical solver");
      
      po::variables_map vm;
      po::store(po::parse_command_line(argc,argv,desc), vm);
      po::notify(vm);

      if(vm.count("mode")) {
	std::string m_ = boost::any_cast<std::string>(vm["mode"].value());
	if(m_ == "euler") mode = tracking::TrackMode::Euler;
	else if(m_ == "rk4") mode = tracking::TrackMode::RungeKutta4;
	else throw std::invalid_argument("This mode is not supported.");
      }
      else
	throw std::invalid_argument("Please specify a mode");

      std::cout << "--- Executable options ---" << std::endl;
      for (const auto& it : vm) {
	std::cout << it.first.c_str() << ": ";
	auto& value = it.second.value();
	if (auto v = boost::any_cast<uint32_t>(&value))
	  std::cout << *v << std::endl;
	else if (auto v = boost::any_cast<std::string>(&value))
	  std::cout << *v << std::endl;
	else
	  std::cerr << "type missing" << std::endl;
      }
    }

  //run simulation
  run(mode);

  myapp.Run();
  
  return 0;
}
