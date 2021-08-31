#include "include/geometry.h"
#include "include/tracking.h"
#include "include/generator.h"
#include "include/tqdm.h"
#include "include/utils.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <boost/program_options.hpp>

#include <TApplication.h>
#include "Math/Vector3D.h" // XYZVector
#include "TVector3.h"

#include <TEveManager.h>
#include "TEveLine.h"

#include "TStyle.h"
#include "TRandom.h"

template <typename T> using Vec = std::vector<T>;
template <typename T> using Vec2 = std::vector<std::vector<T>>;

struct InputArgs {
public:
  bool draw;
  float x;
  float y;
  float energy;
  float mass;
  unsigned nparticles;
  float zcutoff;
};

void double_rotation(const TVector3& v, std::pair<float,float> angles) {
}

unsigned size_last_batch(unsigned nbatches, unsigned nelems, unsigned batchSize) {
  return nelems-(nbatches-1)*batchSize;
}

void run(tracking::TrackMode mode, const InputArgs& args)
{
  using XYZ = ROOT::Math::XYZVector;
  gStyle->SetPalette(56); // 53 = black body radiation, 56 = inverted black body radiator, 103 = sunset, 87 == light temperature

  //set global variables
  std::fstream file, file2;
  std::string histname1, histname2;

  //set global parameters
  // XYZ line_perp_plane1(args.x, args.y, -args.zcutoff);
  // XYZ line_perp_plane2(args.x, args.y, args.zcutoff);

  double Bscale = 1.;
  XYZ origin(0.f, 0.f, 0.f);
  std::array<unsigned, tracking::TrackMode::NMODES> nsteps = {{30000, 13000 }};
  std::array<double, tracking::TrackMode::NMODES> stepsize = {{ args.zcutoff/500, 1. }};
  std::array<std::string, tracking::TrackMode::NMODES> suf = {{ "_euler", "_rk4" }};
    
  //generate random positions around input positions
  NormalDistribution<double> xdist(args.x, 0.1); //beam width of 1 millimeter
  NormalDistribution<double> ydist(args.y, 0.1); //beam width of 1 millimeter
  BoltzmannDistribution<float> boltzdist(1.f, 0.15, 4, 0.138);  //boltzdist.test("data/boltz.csv");
  UniformDistribution<float> phidist(-M_PI, M_PI);
  UniformDistribution<float> etadist(-2.f, 2.f);
  
  Vec<Magnet> magnetInfo{
     // {Magnet::DipoleY,    "D1_neg", kBlue,    std::make_pair(0.,-3.529),
     //  Geometry::Dimensions{-10., 10., -10., 10., -5840.0-945.0, -5840.0}
     // },
     
     // {Magnet::Quadrupole, "Q4_neg", kYellow,  std::make_pair(200.34,-200.34),
     //  Geometry::Dimensions{-10., 10., -10., 10., -4730.0-630.0, -4730.0}
     // },
     
     // {Magnet::Quadrupole, "Q3_neg", kYellow,  std::make_pair(-200.34,200.34),
     //  Geometry::Dimensions{-10., 10., -10., 10., -3830.0-550.0, -3830.0}
     // },
     
     // {Magnet::Quadrupole, "Q2_neg", kYellow,  std::make_pair(-200.34,200.34),
     //  Geometry::Dimensions{-10., 10., -10., 10., -3180.0-550.0, -3180.0}
     // },
     
     // {Magnet::Quadrupole, "Q1_neg", kYellow,  std::make_pair(200.34,-200.34),
     //  Geometry::Dimensions{-10., 10., -10., 10., -2300.0-630.0, -2300.0}
     // },
      
     // {Magnet::DipoleX,    "D_corr", kBlue+1,  std::make_pair(-1.1716,0.),
     //  Geometry::Dimensions{-10., 10., -10., 10., -1920.0-190.0, -1920.0}
     // },
      
     // {Magnet::DipoleX,    "Muon"  , kMagenta, std::make_pair(0.67,0.),
     //  Geometry::Dimensions{-10., 10., -10., 10., -750.0-430.0,  -750.0}
     // },
      
     // {Magnet::Quadrupole, "Q1_pos", kYellow,  std::make_pair(200.34,-200.34),
     //  Geometry::Dimensions{-10., 10., -10., 10., 2300.0, 2300.0+630.0}
     // },
      
     // {Magnet::Quadrupole, "Q2_pos", kYellow,  std::make_pair(-200.34,200.34),
     //  Geometry::Dimensions{-10., 10., -10., 10., 3180.0, 3180.0+550.0}
     // },
      
     // {Magnet::Quadrupole, "Q3_pos", kYellow,  std::make_pair(-200.34,200.34),
     //  Geometry::Dimensions{-10., 10., -10., 10., 3830.0, 3830.0+550.0}
     // },
      
     // {Magnet::Quadrupole, "Q4_pos", kYellow,  std::make_pair(200.34,-200.34),
     //  Geometry::Dimensions{-10., 10., -10., 10., 4730.0, 4730.0+630.0}
     // },
      
     // {Magnet::DipoleY,    "D1_pos", kBlue,    std::make_pair(0.,-3.529),
     //  Geometry::Dimensions{-10., 10., -10., 10., 5840.0, 5840.0+945.0}
     // }
  };

  // Vec<Magnets::Magnet> magnetInfo{ {Magnets::DipoleY,
  // 					    "D1_neg",
  // 					    kGreen,
  // 					    std::make_pair(0.,-0.5),
  // 					    std::make_pair(particle.pos.Z()-1500, particle.pos.Z()+3000), 60., 60.},
  // };

  //Vec<Magnets::Magnet> magnetInfo{};
  MagnetSystem magnets(magnetInfo);

  //figure 3.3 in ALICE ZDC TDR (which does not agree perfectly with the text: see dimensions in Chapters 3.4 and 3.5)
  //available on July 15th 2021 here: https://cds.cern.ch/record/381433/files/Alice-TDR.pdf

  Vec<Calo> caloInfo{
      // {Calo::Neutron, "NeutronZDC", kCyan-3, Dimensions{-8/2., 8./2., -8/2., 8/2., -11613, -11613+100}},
      // {Calo::Proton,  "ProtonZDC",  kCyan+3, Dimensions{10.82, 10.82+22., -13./2., 13./2., -11563, -11563+150}}
  };
  CaloSystem calos(caloInfo);

  if(args.draw) {
    float beamcap = args.zcutoff+100;
    BuildGeom(Dimensions{0., 0., 0., 0., -beamcap, beamcap}, //beamline coordinates
	      magnets, calos);
  }
  
  constexpr float batchSize = 1500.f;
  const unsigned nbatches = ceil(args.nparticles/batchSize);
  std::cout << " --- Simulation Information --- " << std::endl;
  std::cout << "Batch Size: " << batchSize << " (last batch: " << size_last_batch(nbatches, args.nparticles, batchSize) << ")" << std::endl;
  std::cout << "Number of batches: " << nbatches << std::endl;
  std::cout << "Step Size: " << stepsize[mode] << std::endl;
  std::cout << "--------------------------" << std::endl;
  unsigned batchSize_;

  std::string roundx = std::to_string(args.x).substr(0,8);
  std::string roundy = std::to_string(args.y).substr(0,8);
  std::string rounden = std::to_string(args.energy).substr(0,10);
  std::string roundsz = std::to_string(stepsize[mode]);
  std::replace( roundx.begin(), roundx.end(), '.', 'p');
  std::replace( roundy.begin(), roundy.end(), '.', 'p');
  std::replace( rounden.begin(), rounden.end(), '.', 'p');
  std::replace( roundsz.begin(), roundsz.end(), '.', 'p');
  std::string str_initpos = "_" + roundx + "X_" + roundy + "Y_" + rounden + "En_";
  std::string sz = roundsz + "SZ";

  std::string filename("data/track" + suf[mode] + str_initpos + sz + ".csv");
  std::string filename2("data/histo" + suf[mode] + str_initpos + sz + ".csv");
  file2.open(filename2, std::ios_base::out);
      
  for (unsigned ibatch : tq::trange(nbatches))
    //for(unsigned ibatch=0; ibatch<nbatches; ++ibatch)
    {
      batchSize_ = ibatch==nbatches-1 ? size_last_batch(nbatches, args.nparticles, batchSize) : batchSize;
      
      //define the initial properties of the incident particle
      Vec<Particle> p1(batchSize_);
      Vec<Particle> p2(batchSize_);
      for(unsigned i=0; i<batchSize_; ++i) {
	//negative z side
	p1[i].pos = XYZ( xdist.generate(), ydist.generate(), -args.zcutoff-50 ); // cm //-7000
	p1[i].mom = XYZ(0.0, 0.0, args.energy); // GeV/c
	p1[i].mass = args.mass; // GeV/c^2
	p1[i].energy = args.energy; //TMath::Sqrt(particle.mom.Mag2() + particle.mass*particle.mass);
	p1[i].charge = +1;
	//positive z side
	p2[i].pos = XYZ( xdist.generate(), ydist.generate(), args.zcutoff+50 ); // cm //7000
	p2[i].mom = XYZ(0.0, 0.0, -args.energy); // GeV/c
	p2[i].mass = args.mass; // GeV/c^2
	p2[i].energy = args.energy;
	p2[i].charge = +1;
      }

      Vec<TEveLine*> particleTrackViz1(batchSize_);
      Vec<TEveLine*> particleTrackViz2(batchSize_);
      if(args.draw) {
	for(unsigned i=0; i<batchSize_; ++i) {
	  particleTrackViz1[i] = new TEveLine();
	  particleTrackViz2[i] = new TEveLine();
	}
      }

      Vec<SimParticle> simp1;
      Vec<SimParticle> simp2;
      for(unsigned i=0; i<batchSize_; ++i) {
	simp1.push_back( SimParticle(p1[i], nsteps[mode], stepsize[mode]) );
	simp2.push_back( SimParticle(p2[i], nsteps[mode], stepsize[mode]) );
      }

      Vec<const Track*> tracks1(batchSize_);
      Vec<const Track*> tracks2(batchSize_);
      for(unsigned i=0; i<batchSize_; ++i) {
	tracks1[i] = &( simp1[i].track( magnets, mode, Bscale, args.zcutoff ));
	tracks2[i] = &( simp2[i].track( magnets, mode, Bscale, args.zcutoff ));
      }

      Vec2<double> itEnergies1(batchSize_), itEnergies2(batchSize_);
      Vec2<XYZ> itPositions1(batchSize_), itPositions2(batchSize_);
      Vec2<XYZ> itMomenta1(batchSize_), itMomenta2(batchSize_);
      Vec<unsigned> nStepsUsed1(batchSize_), nStepsUsed2(batchSize_);

      Vec<float> psi1(batchSize_);
      Vec<float> psi2(batchSize_);
      Vec<float> psi_angles(batchSize_);
      Vec<float> corr(batchSize_);

      std::pair<float,float> nomAngles = calculate_angles_to_beamline(args.x, args.y, args.zcutoff);

      //unit vectors
      TVector3 uvX_1(1,0,0), uvY_1(0,1,0);
      TVector3 uvX_2(1,0,0), uvY_2(0,1,0);

      uvX_1.RotateY( -nomAngles.second );
      uvY_1.Rotate( nomAngles.first, uvX_1);

      uvX_2.RotateY( nomAngles.second );
      uvY_2.Rotate( -nomAngles.first, uvX_2 );

      // std::cout << TMath::Sqrt(uvX_1.Mag2()) << std::endl;
      // print_pos<TVector3>("X", uvX_1);
      // print_pos<TVector3>("Y", uvY_1);
      // std::cout << TMath::Sqrt(uvX_1.Mag2()) << std::endl;
      // std::cout << TMath::Sqrt(uvY_1.Mag2()) << std::endl;
      // std::cout << "---" << std::endl;
      // std::exit(0);

      if(uvX_1.Dot(uvY_1) > 1e-15 or uvX_2.Dot(uvY_2) > 1e-15) {
	std::cout << "The vectors must be perpendicular!" << std::endl;
	std::cout << uvX_1.Dot(uvY_1) << ", " << uvX_2.Dot(uvY_2) << std::endl;
	std::exit(0);
      }
            
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

	XYZ last1_ = itPositions1[i].back();
	TVector3 last1V(last1_.X(), last1_.Y(), last1_.Z());
	last1V.Rotate( -nomAngles.first, uvX_1);
	last1V.RotateY( nomAngles.second );
	
	double last1X_ = last1V.Dot(uvX_1);
	double last1Y_ = last1V.Dot(uvY_1);
	
	XYZ last2_ = itPositions2[i].back();
	TVector3 last2V(last2_.X(), last2_.Y(), last2_.Z());
	last2V.Rotate( nomAngles.first, uvX_2);
	last2V.RotateY( -nomAngles.second );

	double last2X_ = last2V.Dot(uvX_2);
	double last2Y_ = last2V.Dot(uvY_2);

	psi1[i] = std::atan2( last1Y_, last1X_ ) + M_PI;
	psi2[i] = std::atan2( last2Y_, last2X_ ) + M_PI;

	//check if the two particles "crossed"
	//this catches number of iterations that are too small
	assert(last1_.Z() > last2_.Z());
	
	// // XYZ point_check(0., 0., 1.);
	
	// //calculate intersection between particle trajectory and plane
	// XYZ is1 = intersect_plane_with_line(line_perp_plane1, origin,
	// 				    origin, last1V, last1V);
	// XYZ is2 = intersect_plane_with_line(line_perp_plane2, origin,
	// 				    origin, last2V, last2V);
	// // XYZ is1_check = intersect_plane_with_line(point_check, origin,
	// // 					  origin, last1V, last1V);
	// // XYZ is2_check = intersect_plane_with_line(point_check, origin,
	// // 					  origin, last2V, last2V);
	
	// //calculate intersection between beam line and plane
	// XYZ is1_origin = intersect_plane_with_line(line_perp_plane1, origin,
	// 					   line_perp_plane1, origin,
	// 					   last1V);
	// XYZ is2_origin = intersect_plane_with_line(line_perp_plane2, origin,
	// 					   line_perp_plane2, origin,
	// 					   last2V);
	
	// //cordinate transformation	
	// is1 = rotate_coordinates(is1, nomAngles.first, nomAngles.second);
	// is2 = rotate_coordinates(is2, -1*(nomAngles.first), -1*(nomAngles.second));

	// is1_origin = rotate_coordinates(is1_origin, nomAngles.first, nomAngles.second);
	// is2_origin = rotate_coordinates(is2_origin, -1*(nomAngles.first), -1*(nomAngles.second));
	
	// //cordinate translation
	// is1 = translate_coordinates(is1, is1_origin);
	// is2 = translate_coordinates(is2, is2_origin);
	
	// //psi angles do not depend on Z
	// psi1[i] = std::atan2( is1.Y(), is1.X() ) + M_PI;
	// psi2[i] = std::atan2( is2.Y(), is2.X() ) + M_PI;
	

	float psi2_tmp = psi2[i]+M_PI>2*M_PI ? psi2[i]-M_PI : psi2[i]+M_PI;
	psi_angles[i] = distance_two_angles(psi1[i], psi2_tmp);
	psi_angles[i] /= 2.;
      }

      unsigned minelem = *std::min_element(std::begin(nStepsUsed1), std::end(nStepsUsed1));
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
	}
	  
      // file.open(filename, std::ios_base::out);
      // for(unsigned i_step = 0; i_step<minelem; i_step++)
      // 	{
      // 	  if (!file.is_open()) 
      // 	    std::cerr << "failed to open " << filename << '\n';
      // 	  else {
      // 	    if(i_step==0)
      // 	      file << "x1,y1,z1,energy1,x2,y2,z2,energy2" << std::endl;
      // 	    file << std::to_string( itPositions1[0][i_step].X() ) << ","
      // 		 << std::to_string( itPositions1[0][i_step].Y() ) << ","
      // 		 << std::to_string( itPositions1[0][i_step].Z() ) << ","
      // 		 << std::to_string( itEnergies1[0][i_step] ) << ","
      // 		 << std::to_string( itPositions2[0][i_step].X() ) << ","
      // 		 << std::to_string( itPositions2[0][i_step].Y() ) << ","
      // 		 << std::to_string( itPositions2[0][i_step].Z() ) << ","
      // 		 << std::to_string( itEnergies2[0][i_step] ) << ","
      // 		 << std::endl;
      // 	  }

      // 	}
      // file.close();

      for(unsigned ix=0; ix<batchSize_; ix++) {
  
	unsigned id1 = get_index_closer_to_origin(itPositions1[ix], minelem);
	unsigned id2 = get_index_closer_to_origin(itPositions2[ix], minelem);
	
	if(ix==0 and ibatch==0)
	  file2 << "iBatch,Idx,sumMomX,sumMomY,sumMomZ,PsiA,PsiB,Psi,Phi,Cos" << std::endl;

	XYZ momSum = itMomenta1[ix][id1] + itMomenta2[ix][id2];
	//pT_boltzmann->Boost(Li->BoostVector())
	// where Li is a LorentzVector (pTi, 2GeV)

	ROOT::Math::PtEtaPhiMVector boltzPT(boltzdist.generate(),
					    etadist.generate(), phidist.generate(),
					    args.mass);

	float nNucleons = 200.f;
	ROOT::Math::PxPyPzEVector kickPT(momSum.X()/nNucleons,
					 momSum.Y()/nNucleons,
					 momSum.Z()/nNucleons,
					 args.mass);
	ROOT::Math::PtEtaPhiMVector totalPT = boltzPT + kickPT;
	
	float totalPhi = totalPT.Phi();
	if(totalPhi<0)
	  totalPhi = 2*M_PI + totalPhi; //convert from [-Pi;Pi[ to [0;2Pi[
		
	corr[ix] = std::cos( distance_two_angles(totalPhi, psi_angles[ix]) );
	
	file2 << std::to_string( ibatch ) << ","
	      << std::to_string( ix ) << ","
	      << std::to_string( momSum.X() ) << ","
	      << std::to_string( momSum.Y() ) << ","
	      << std::to_string( momSum.Z() ) << ","
	      << std::to_string( psi1[ix] ) << ","
	      << std::to_string( psi2[ix] ) << ","
	      << std::to_string( psi_angles[ix] ) << ","
	      << std::to_string( totalPhi ) << ","
	      << std::to_string( corr[ix] ) 
	      << std::endl;	
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

// run example: ./v1_beam.exe --mode euler --x 0.08 --y 0.08 --energy 1380 --nparticles 1 --zcutoff 5000.
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
    ("nparticles", po::value<unsigned>()->default_value(1), "number of particles to generate on each beam")
    ("zcutoff", po::value<float>()->default_value(5000.f), "cutoff at which to apply the fake deflection");
      
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
  info.mass = 0.938; //GeV
  info.nparticles = boost::any_cast<unsigned>(vm["nparticles"].value());
  info.zcutoff = boost::any_cast<float>(vm["zcutoff"].value());
  assert(info.zcutoff > 0);
  
  run(mode, info);

  if(flag_draw)
    myapp.Run();

  std::cout << std::endl;
  return 0;
}
