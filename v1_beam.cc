#include "include/geometry.h"
#include "include/tracking.h"
#include "include/tqdm.h"

#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <boost/program_options.hpp>

#include <TApplication.h>
#include "Math/Vector3D.h" // XYZVector

#include <TEveManager.h>
#include "TEveLine.h"

#include "TStyle.h"
#include "TRandom.h"

struct InputArgs {
public:
  bool draw;
  float x;
  float y;
  float energy;
  unsigned nparticles;
};

unsigned size_last_batch(unsigned nbatches, unsigned nelems, unsigned batchSize) {
  return nelems-(nbatches-1)*batchSize;
}

void run(tracking::TrackMode mode, const InputArgs& args)
{
  gRandom->SetSeed(0);
  gStyle->SetPalette(56); // 53 = black body radiation, 56 = inverted black body radiator, 103 = sunset, 87 == light temperature

  //set global parameters
  std::fstream file, file2;
  std::string histname1, histname2;
  double Bscale = 1.;
  std::array<unsigned, tracking::TrackMode::NMODES> nsteps = {{ 15000, 13000 }};
  std::array<double, tracking::TrackMode::NMODES> stepsize = {{ .1, 1. }};
  std::array<std::string, tracking::TrackMode::NMODES> suf = {{ "_euler", "_rk4" }};
    
  // generate random positions around input positions
  std::random_device seeder;
  std::mt19937 rng( seeder() );
  std::normal_distribution<double> xdist( args.x, 0.1 ); //beam width of 1 millimeter
  std::normal_distribution<double> ydist( args.y, 0.1 ); //beam width of 1 millimeter
  
  // std::vector<Magnets::Magnet> magnetInfo{
  //    {Magnets::DipoleY,    "D1_neg", kBlue,    std::make_pair(0.,-3.529),
  //     Geometry::Dimensions{-10., 10., -10., 10., -5840.0-945.0, -5840.0}
  //    },
     
  //    {Magnets::Quadrupole, "Q4_neg", kYellow,  std::make_pair(200.34,-200.34),
  //     Geometry::Dimensions{-10., 10., -10., 10., -4730.0-630.0, -4730.0}
  //    },
     
  //    {Magnets::Quadrupole, "Q3_neg", kYellow,  std::make_pair(-200.34,200.34),
  //     Geometry::Dimensions{-10., 10., -10., 10., -3830.0-550.0, -3830.0}
  //    },
     
  //    {Magnets::Quadrupole, "Q2_neg", kYellow,  std::make_pair(-200.34,200.34),
  //     Geometry::Dimensions{-10., 10., -10., 10., -3180.0-550.0, -3180.0}
  //    },
     
  //    {Magnets::Quadrupole, "Q1_neg", kYellow,  std::make_pair(200.34,-200.34),
  //     Geometry::Dimensions{-10., 10., -10., 10., -2300.0-630.0, -2300.0}
  //    },
      
  //    {Magnets::DipoleX,    "D_corr", kBlue+1,  std::make_pair(-1.1716,0.),
  //     Geometry::Dimensions{-10., 10., -10., 10., -1920.0-190.0, -1920.0}
  //    },
      
  //    {Magnets::DipoleX,    "Muon"  , kMagenta, std::make_pair(0.67,0.),
  //     Geometry::Dimensions{-10., 10., -10., 10., -750.0-430.0,  -750.0}
  //    },
      
  //    {Magnets::Quadrupole, "Q1_pos", kYellow,  std::make_pair(200.34,-200.34),
  //     Geometry::Dimensions{-10., 10., -10., 10., 2300.0, 2300.0+630.0}
  //    },
      
  //    {Magnets::Quadrupole, "Q2_pos", kYellow,  std::make_pair(-200.34,200.34),
  //     Geometry::Dimensions{-10., 10., -10., 10., 3180.0, 3180.0+550.0}
  //    },
      
  //    {Magnets::Quadrupole, "Q3_pos", kYellow,  std::make_pair(-200.34,200.34),
  //     Geometry::Dimensions{-10., 10., -10., 10., 3830.0, 3830.0+550.0}
  //    },
      
  //    {Magnets::Quadrupole, "Q4_pos", kYellow,  std::make_pair(200.34,-200.34),
  //     Geometry::Dimensions{-10., 10., -10., 10., 4730.0, 4730.0+630.0}
  //    },
      
  //    {Magnets::DipoleY,    "D1_pos", kBlue,    std::make_pair(0.,-3.529),
  //     Geometry::Dimensions{-10., 10., -10., 10., 5840.0, 5840.0+945.0}
  //    }
  // };

  // std::vector<Magnets::Magnet> magnetInfo{ {Magnets::DipoleY,
  // 					    "D1_neg",
  // 					    kGreen,
  // 					    std::make_pair(0.,-0.5),
  // 					    std::make_pair(particle.pos.Z()-1500, particle.pos.Z()+3000), 60., 60.},
  // };

  std::vector<Magnets::Magnet> magnetInfo{};
  Magnets magnets(magnetInfo, args.draw);

  //figure 3.3 in ALICE ZDC TDR (which does not agree perfectly with the text: see dimensions in Chapters 3.4 and 3.5)
  //available on July 15th 2021 here: https://cds.cern.ch/record/381433/files/Alice-TDR.pdf

  std::vector<Calorimeters::Calorimeter> caloInfo{
      // {Calorimeters::Neutron, "NeutronZDC", kCyan-3, Geometry::Dimensions{-8/2., 8./2., -8/2., 8/2., -11613, -11613+100}},
      // {Calorimeters::Proton,  "ProtonZDC",  kCyan+3, Geometry::Dimensions{10.82, 10.82+22., -13./2., 13./2., -11563, -11563+150}}
  };
  Calorimeters calos(caloInfo, args.draw);
  
  constexpr float batchSize = 50.f; //each batch will hold 100 particles (positive z + negative z sides)
  const unsigned nbatches = ceil(args.nparticles/batchSize);
  std::cout << " --- Simulation Informations --- " << std::endl;
  std::cout << "Batch Size: " << batchSize << " (last batch: " << size_last_batch(nbatches, args.nparticles, batchSize) << ")" << std::endl;
  std::cout << "Number of batches: " << nbatches << std::endl;
  std::cout << "--------------------------" << std::endl;
  unsigned batchSize_;

  std::string roundx = std::to_string(args.x).substr(0,8);
  std::string roundy = std::to_string(args.y).substr(0,8);
  std::string rounden = std::to_string(args.energy).substr(0,10);
  std::replace( roundx.begin(), roundx.end(), '.', 'p');
  std::replace( roundy.begin(), roundy.end(), '.', 'p');
  std::replace( rounden.begin(), rounden.end(), '.', 'p');
  std::string str_initpos = "_" + roundx + "X_" + roundy + "Y_" + rounden + "En";

  std::string filename("data/track" + suf[mode] + str_initpos + ".csv");
  std::string filename2("data/histo" + suf[mode] + str_initpos + ".csv");
  file2.open(filename2, std::ios_base::out);
      
  for (unsigned ibatch : tq::trange(nbatches))
    //for(unsigned ibatch=0; ibatch<nbatches; ++ibatch)
    {
      batchSize_ = ibatch==nbatches-1 ? size_last_batch(nbatches, args.nparticles, batchSize) : batchSize;
      
      //define the initial properties of the incident particle
      std::vector<Particle> p1(batchSize_);
      std::vector<Particle> p2(batchSize_);
      for(unsigned i=0; i<batchSize_; ++i) {
	//negative z side
	p1[i].pos = ROOT::Math::XYZVector( xdist(rng), ydist(rng), -500.0 ); // cm //-7000
	p1[i].mom = ROOT::Math::XYZVector(0.0, 0.0, args.energy); // GeV/c
	p1[i].mass = 0.938; // GeV/c^2
	p1[i].energy = args.energy; //TMath::Sqrt(particle.mom.Mag2() + particle.mass*particle.mass);
	p1[i].charge = +1;
	//positive z side
	p2[i].pos = ROOT::Math::XYZVector( xdist(rng), ydist(rng), 500.0 ); // cm //7000
	p2[i].mom = ROOT::Math::XYZVector(0.0, 0.0, -args.energy); // GeV/c
	p2[i].mass = 0.938; // GeV/c^2
	p2[i].energy = args.energy;
	p2[i].charge = +1;
      }

      std::vector<TEveLine*> particleTrackViz1(batchSize_);
      std::vector<TEveLine*> particleTrackViz2(batchSize_);
      if(args.draw) {
	for(unsigned i=0; i<batchSize_; ++i) {
	  particleTrackViz1[i] = new TEveLine();
	  particleTrackViz2[i] = new TEveLine();
	}
      }

      //-------------------------------------------------------------------------------------------
      // Scanned sigma x (m) as a function of z (m) of the LHC beam  (John Jowett)
      // Double_t z_LHC[37]       = {-67.643,-58.738,-54.063,-50.500,-47.161,-44.267,-40.259,-38.701,-32.244,-30.241,-26.679,-23.784,-20.222,-17.105,-13.098,-9.3135,-4.6382,-0.4081,3.37662,8.05195,13.6178,17.6252,20.2968,22.7458,26.0853,28.9796,31.8738,34.7681,38.1076,40.7792,44.1187,46.3451,49.2393,53.0241,58.3673,63.4879,69.0538};
      // Double_t sigma_x_LHC[37] = {0.0013124,0.0014062,0.0014562,0.0013687,0.0012687,0.0010499,0.0008562,0.0007999,0.0007562,0.0007812,0.0007624,0.0007437,0.0006249,0.0005249,0.0004062,0.0002812,0.0001437,4.37485e-05,0.0001312,0.0002687,0.0004374,0.0005624,0.0006562,0.0007187,0.0009124,0.0010812,0.0013124,0.0014312,0.0015624,0.0014874,0.0014249,0.0012374,0.0011374,0.0010187,0.0010062,0.0009812,0.0009687};
      //-------------------------------------------------------------------------------------------

      std::vector<SimParticle> simp1;
      std::vector<SimParticle> simp2;
      for(unsigned i=0; i<batchSize_; ++i) {
	simp1.push_back( SimParticle(p1[i], nsteps[mode], stepsize[mode]) );
	simp2.push_back( SimParticle(p2[i], nsteps[mode], stepsize[mode]) );
      }

      std::vector<const Track*> tracks1(batchSize_);
      std::vector<const Track*> tracks2(batchSize_);
      for(unsigned i=0; i<batchSize_; ++i) {
	tracks1[i] = &( simp1[i].track( magnets, mode, Bscale ));
	tracks2[i] = &( simp2[i].track( magnets, mode, Bscale ));
      }

      std::vector< std::vector<double> > itEnergies1(batchSize_), itEnergies2(batchSize_);
      std::vector<std::vector<ROOT::Math::XYZVector>> itPositions1(batchSize_), itPositions2(batchSize_);
      std::vector<std::vector<ROOT::Math::XYZVector>> itMomenta1(batchSize_), itMomenta2(batchSize_);
      std::vector<unsigned> nStepsUsed1(batchSize_), nStepsUsed2(batchSize_);

      for(unsigned i=0; i<batchSize_; ++i) {
	//negative z side
	itEnergies1[i]  = tracks1[i]->energies();
	itPositions1[i] = tracks1[i]->positions();
	itMomenta1[i]   = tracks1[i]->momenta();
	nStepsUsed1[i]  = tracks1[i]->steps_used();
	//positive z side
	itEnergies2[i]  = tracks2[i]->energies();
	itPositions2[i] = tracks2[i]->positions();
	itMomenta2[i]   = tracks2[i]->momenta();
	nStepsUsed2[i]  = tracks2[i]->steps_used();
      }

      unsigned minelem = *std::min_element(std::begin(nStepsUsed1), std::end(nStepsUsed1));
      file.open(filename, std::ios_base::out);
      for(unsigned i_step = 0; i_step<minelem; i_step++)
	{
	  if(args.draw) {
	    for(unsigned ix=0; ix<batchSize_; ix++) {
	      particleTrackViz1[ix]->SetNextPoint(itPositions1[ix][i_step].X(),
						  itPositions1[ix][i_step].Y(),
						  itPositions1[ix][i_step].Z() );

	      particleTrackViz2[ix]->SetNextPoint(itPositions2[ix][i_step].X(),
						  itPositions2[ix][i_step].Y(),
						  itPositions2[ix][i_step].Z() );
	    }
	  }
  
	  if (!file.is_open()) 
	    std::cerr << "failed to open " << filename << '\n';
	  else {
	    if(i_step==0)
	      file << "x1,y1,z1,energy1,x2,y2,z2,energy2" << std::endl;
	    file << std::to_string( itPositions1[0][i_step].X() ) << ","
		 << std::to_string( itPositions1[0][i_step].Y() ) << ","
		 << std::to_string( itPositions1[0][i_step].Z() ) << ","
		 << std::to_string( itEnergies1[0][i_step] ) << ","
		 << std::to_string( itPositions2[0][i_step].X() ) << ","
		 << std::to_string( itPositions2[0][i_step].Y() ) << ","
		 << std::to_string( itPositions2[0][i_step].Z() ) << ","
		 << std::to_string( itEnergies2[0][i_step] ) << ","
		 << std::endl;
	  }

	}
      file.close();

      for(unsigned ix=0; ix<batchSize_; ix++) {
	bool origin_found = false;
	for(unsigned i_step = 0; i_step<minelem; i_step++) {

	  if( TMath::Sqrt( (itPositions1[ix][i_step]-itPositions2[ix][i_step]).Mag2() ) < .1 )
	    {
	      origin_found = true;
	      if(ix==0 and ibatch==0)
		file2 << "iBatch,Idx,mom1X,mom1Y,mom1Z,mom2X,mom2Y,mom2Z,sumMomX,sumMomY,sumMomZ" << std::endl;
	      file2 << std::to_string( ibatch ) << ","
		    << std::to_string( ix ) << ","
		    << std::to_string( itMomenta1[ix][i_step].X() ) << ","
		    << std::to_string( itMomenta1[ix][i_step].Y() ) << ","
		    << std::to_string( itMomenta1[ix][i_step].Z() ) << ","
		    << std::to_string( itMomenta2[ix][i_step].X() ) << ","
		    << std::to_string( itMomenta2[ix][i_step].Y() ) << ","
		    << std::to_string( itMomenta2[ix][i_step].Z() ) << ","
		    << std::to_string( itMomenta1[ix][i_step].X() + itMomenta2[ix][i_step].X() ) << ","
		    << std::to_string( itMomenta1[ix][i_step].Y() + itMomenta2[ix][i_step].Y() ) << ","
		    << std::to_string( itMomenta1[ix][i_step].Z() + itMomenta2[ix][i_step].Z() )
		    << std::endl;
	      break; //only at most one output row per particle pair
	    }
	}
	if(!origin_found)
	  std::cout << "ERROR: Batch " << ibatch << " Idx " << ix << " no origin found." << std::endl;
      }
	    
      if(args.draw) {
	for(unsigned ix=0; ix<batchSize_; ix++) {
	  histname1 = "track_zpos_ " + std::to_string(ix);
	  particleTrackViz1[ix]->SetName( histname1.c_str() );
	  particleTrackViz1[ix]->SetLineStyle(1);
	  particleTrackViz1[ix]->SetLineWidth(2);
	  particleTrackViz1[ix]->SetMainAlpha(0.7);
	  particleTrackViz1[ix]->SetMainColor(kRed+3);

	  histname2 = "track_zneg_ " + std::to_string(ix);
	  particleTrackViz2[ix]->SetName( histname2.c_str() );
	  particleTrackViz2[ix]->SetLineStyle(1);
	  particleTrackViz2[ix]->SetLineWidth(2);
	  particleTrackViz2[ix]->SetMainAlpha(0.7);
	  particleTrackViz2[ix]->SetMainColor(kRed-7);
	  gEve->AddElement(particleTrackViz1[ix]);
	  gEve->AddElement(particleTrackViz2[ix]);
	}
      }

    } // for ibatch

  file2.close();
  
  if(args.draw)
    gEve->Redraw3D(kTRUE);
}

// run example: ./v1_beam.exe --mode euler --x 0.01 --y 0.01 --energy 1380 --nparticles 1
int main(int argc, char **argv) {
  TApplication myapp("myapp", &argc, argv);
  
  tracking::TrackMode mode = tracking::TrackMode::Euler;
  bool flag_draw = false;
 
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "produce this help message")
    ("mode", po::value<std::string>()->default_value("euler"), "numerical solver")
    ("draw", po::bool_switch(&flag_draw), "whether to draw the geometry with ROOT's Event Display")
    ("x", po::value<float>()->default_value(0.f), "initial beam x position")
    ("y", po::value<float>()->default_value(0.f), "initial beam y position")
    ("energy", po::value<float>()->default_value(1380.f), "beam energy position")
    ("nparticles", po::value<unsigned>()->default_value(1), "number of particles to generate on each beam");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc), vm);
  po::notify(vm);

  if(vm.count("help") or argc<2) {
    std::cerr << desc << std::endl;
    std::exit(0);
  }
      
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
    if (auto v = boost::any_cast<float>(&value))
      std::cout << *v << std::endl;
    else if (auto v = boost::any_cast<bool>(&value)) {
      std::string str_ = *v==1 ? "true" : "false";
      std::cout << *v << std::endl;
    }
    else if (auto v = boost::any_cast<std::string>(&value))
      std::cout << *v << std::endl;
    else if (auto v = boost::any_cast<unsigned>(&value))
      std::cout << *v << std::endl;
    else
      std::cerr << "type missing" << std::endl;
  }

  //run simulation   
  InputArgs info;
  info.draw = flag_draw;
  info.x = boost::any_cast<float>(vm["x"].value());
  info.y = boost::any_cast<float>(vm["y"].value());
  info.energy = boost::any_cast<float>(vm["energy"].value());
  info.nparticles = boost::any_cast<unsigned>(vm["nparticles"].value());
  run(mode, info);

  if(flag_draw)
    myapp.Run();

  std::cout << std::endl;
  return 0;
}
