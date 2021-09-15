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
#include "TLorentzVector.h"

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
  float energy_scale;
  float width_scale;
  float mass;
  float mass_interaction;
  unsigned npartons;
  unsigned nparticles;
  float zcutoff;
};

struct Globals {
  static constexpr float distanceToDetector = 11600; //cm
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

  double Bscale = 1.;
  XYZ origin(0.f, 0.f, 0.f);
  std::array<unsigned, tracking::TrackMode::NMODES> nsteps = {{30000, 13000 }};
  std::array<double, tracking::TrackMode::NMODES> stepsize = {{ static_cast<double>(args.zcutoff)/500., 1. }};
  std::array<std::string, tracking::TrackMode::NMODES> suf = {{ "_euler", "_rk4" }};
    
  //generate random positions around input positions
  NormalDistribution<double> xdist(args.x, args.width_scale * 0.1); //beam width of 1 millimeter
  NormalDistribution<double> ydist(args.y, args.width_scale * 0.1); //beam width of 1 millimeter
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
  std::string roundens = std::to_string(args.energy_scale).substr(0,5);
  std::string roundws = std::to_string(args.width_scale).substr(0,5);
  std::string roundsz = std::to_string(stepsize[mode]);
  std::string roundmint = std::to_string(args.mass_interaction).substr(0,8);
  std::string roundnpartons = std::to_string(args.npartons).substr(0,8);
  std::replace( roundx.begin(), roundx.end(), '.', 'p');
  std::replace( roundy.begin(), roundy.end(), '.', 'p');
  std::replace( rounden.begin(), rounden.end(), '.', 'p');
  std::replace( roundens.begin(), roundens.end(), '.', 'p');
  std::replace( roundws.begin(), roundws.end(), '.', 'p');
  std::replace( roundsz.begin(), roundsz.end(), '.', 'p');
  std::replace( roundmint.begin(), roundmint.end(), '.', 'p');
  std::replace( roundnpartons.begin(), roundnpartons.end(), '.', 'p');
  std::string str_initpos = "_" + roundx + "X_" + roundy + "Y_" + rounden + "En_" + roundens + "EnScale_" + roundws + "WScale_";
  std::string extra = roundsz + "SZ_" + roundmint + "MInt_" + roundnpartons + "NP" ;

  std::string filename("data/track" + suf[mode] + str_initpos + extra + ".csv");
  std::string filename2("data/histo" + suf[mode] + str_initpos + extra + ".csv");
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
	p1[i].mom = XYZ(0.0, 0.0, args.energy / static_cast<float>(args.npartons)); // GeV/c
	p1[i].mass = args.mass; // GeV/c^2
	p1[i].energy = args.energy / static_cast<float>(args.npartons);
	p1[i].charge = +1;
	//positive z side
	p2[i].pos = XYZ( xdist.generate(), ydist.generate(), args.zcutoff+50 ); // cm //7000
	p2[i].mom = XYZ(0.0, 0.0, args.energy_scale * -args.energy / static_cast<float>(args.npartons)); // GeV/c
	p2[i].mass = args.mass; // GeV/c^2
	p2[i].energy = args.energy / static_cast<float>(args.npartons);
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

      Vec<float> xHit(batchSize_);
      Vec<float> yHit(batchSize_);
      Vec<float> psi1(batchSize_);
      Vec<float> psi2(batchSize_);
      Vec<unsigned> cat(batchSize_, 99);
      Vec<float> psi_angles(batchSize_);
      Vec<float> corr(batchSize_);

      //std::pair<float,float> nomAngles = calculate_angles_to_beamline(args.x, args.y, args.zcutoff);

      //unit vectors
      TVector3 uZ1(-args.x, -args.y, args.zcutoff);
      uZ1 = uZ1.Unit();
      TVector3 uX1 = uZ1.Orthogonal();
      TVector3 uY1 = uX1;
      uY1.Rotate( M_PI/2, uZ1 );
      assert( (M_PI/2) - uX1.Angle(uY1) < 1e-15 );
      
      TVector3 uZ2(-args.x, -args.y, -args.zcutoff);
      uZ2 = uZ2.Unit();
      TVector3 uX2 = uZ2.Orthogonal();
      TVector3 uY2 = uX2;
      uY2.Rotate( M_PI/2, uZ2 );
      assert( (M_PI/2) - uX2.Angle(uY2) < 1e-15 );

      if(uX1.Dot(uY1) > 1e-15 or uX2.Dot(uY2) > 1e-15) {
	std::cout << "The vectors must be perpendicular!" << std::endl;
	std::cout << uX1.Dot(uY1) << ", " << uX2.Dot(uY2) << std::endl;
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
	TVector3 check1(-p1[i].pos.X(), -p1[i].pos.Y(), args.zcutoff);
	if( check1.Angle(last1V) > 1e-7 ) {
	  std::cout << "The trajectory is not as it should!" << std::endl;
	  std::cout << "Angle1: " << check1.Angle(last1V) << std::endl;
	  std::exit(0);
	}
	double last1X_ = last1_.Dot( uX1 );
	double last1Y_ = last1_.Dot( uY1 );
	
	XYZ last2_ = itPositions2[i].back();
	TVector3 last2V(last2_.X(), last2_.Y(), last2_.Z());
	TVector3 check2(-p2[i].pos.X(), -p2[i].pos.Y(), -args.zcutoff);
	if( check2.Angle(last2V) > 1e-7 ) {
	  std::cout << "The trajectory is not as it should!" << std::endl;
	  std::cout << "Angle2: " << check2.Angle(last2V) << std::endl;
	  std::exit(0);
	}
	double last2X_ = last2_.Dot( uX2 );
	double last2Y_ = last2_.Dot( uY2 );

	TVector3 last1Det = Globals::distanceToDetector * last1V.Unit();
	xHit[i] = last1Det.Dot( uX1 );
	yHit[i] = last1Det.Dot( uY1 );
	psi1[i] = std::atan2( last1Y_, last1X_ ) + M_PI;
	psi2[i] = std::atan2( last2Y_, last2X_ ) + M_PI;

	//define categories according to relative angular difference
	float category_bound = M_PI/6;
	float diff = std::abs(psi1[i]-psi2[i]);
	if( diff < category_bound or diff > 2*M_PI-category_bound )
	  cat[i] = 1;
	else if(diff < M_PI+category_bound and diff > M_PI-category_bound)
	  cat[i] = 2;
	else
	  cat[i] = 0;
	
	//check if the two particles "crossed"
	//this catches number of iterations that are too small
	assert(last1_.Z() > last2_.Z());

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
	  file2 << "iBatch,Idx,sumMomX,sumMomY,sumMomZ,XHit,YHit,PsiA,PsiB,cat1,Psi,Phi,Eta,Cos" << std::endl;
	
	TLorentzVector momLorentz1;
	momLorentz1.SetXYZM(itMomenta1[ix][id1].X(),
			    itMomenta1[ix][id1].Y(),
			    itMomenta1[ix][id1].Z(), args.mass);

	TLorentzVector momLorentz2;
	momLorentz2.SetXYZM(itMomenta2[ix][id2].X(),
			    itMomenta2[ix][id2].Y(),
			    itMomenta2[ix][id2].Z(), args.mass);

	
	TLorentzVector momSum = momLorentz1 + momLorentz2;
	
	//momSum.SetE( TMath::Sqrt( sq(momSum.X()) + sq(momSum.Y()) + sq(momSum.Z()) + sq(args.mass_interaction) ) );

	//print_pos_4D("momLor1", momLorentz1);
	//print_pos_4D("momLor2", momLorentz2);
	//print_pos_4D("momSum", momSum);

	TLorentzVector boltzPT;
	float mass_pion = 0.139;
	float boltzgen = boltzdist.generate();
	float etagen = etadist.generate();
	float phigen = phidist.generate();
	// std::cout << boltzgen << std::endl;
	// std::cout << etagen << std::endl;
	// std::cout << phigen << std::endl;
	// std::exit(0);

	boltzPT.SetPtEtaPhiM(boltzgen,
			     etagen, phigen,
			     mass_pion);

	//std::cout << std::endl;
	// print_pos_4D("boltz before boost", boltzPT);
	// std::cout << "boltz pt: " << TMath::Sqrt(boltzPT.Px()*boltzPT.Px() + boltzPT.Py()*boltzPT.Py()) << std::endl;

	TVector3 bVector = momSum.BoostVector();
	boltzPT.Boost(bVector);
	// print_pos("bVector", bVector);
	// print_pos_4D("boltz after boost", boltzPT);
	// std::exit(0);
	// float nNucleons = 200.f;
	// TLorentzVector kickPT;
	// kickPT.SetPxPyPzE(momSum.Px()/nNucleons,
	// 		  momSum.Py()/nNucleons,
	// 		  momSum.Pz()/nNucleons,
	// 		  args.energy / static_cast<float>(args.npartons));
	// TLorentzVector totalPT = boltzPT + kickPT;
	// std::cout << std::endl;
	// print_pos_4D("kick", kickPT);
	// print_pos_4D("boltz", boltzPT);
	// print_pos_4D("total", totalPT);

	// std::cout << std::endl;
	// std::cout << "kick pt: " << TMath::Sqrt(kickPT.Px()*kickPT.Px() + kickPT.Py()*kickPT.Py()) << ", " << TMath::Sqrt(kickPT.X()*kickPT.X() + kickPT.Y()*kickPT.Y()) << std::endl;
	// std::cout << "boltz pt: " << TMath::Sqrt(boltzPT.Px()*boltzPT.Px() + boltzPT.Py()*boltzPT.Py()) << ", " << TMath::Sqrt(boltzPT.X()*boltzPT.X() + boltzPT.Y()*boltzPT.Y()) << std::endl;
	// std::exit(0);
	float totalPhi = boltzPT.Phi();
	float totalEta = boltzPT.Eta();
	
	if(totalPhi<0)
	  totalPhi = 2*M_PI + totalPhi; //convert from [-Pi;Pi[ to [0;2Pi[
		
	corr[ix] = std::cos( distance_two_angles(totalPhi, psi_angles[ix]) );
	
	file2 << std::to_string( ibatch ) << ","
	      << std::to_string( ix ) << ","
	      << std::to_string( momSum.Px() ) << ","
	      << std::to_string( momSum.Py() ) << ","
	      << std::to_string( momSum.Pz() ) << ","
	      << std::to_string( xHit[ix] ) << ","
	      << std::to_string( yHit[ix] ) << ","
	      << std::to_string( psi1[ix] ) << ","
	      << std::to_string( psi2[ix] ) << ","
	      << std::to_string( cat[ix] ) << ","
	      << std::to_string( psi_angles[ix] ) << ","
	      << std::to_string( totalPhi ) << ","
	      << std::to_string( totalEta ) << ","
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
  //https://www.boost.org/doc/libs/1_45_0/doc/html/boost/program_options/typed_value.html#id903171-bb
  desc.add_options()
    ("help,h", "produce this help message")
    ("mode", po::value<std::string>()->default_value("euler"), "numerical solver")
    ("draw", po::bool_switch(&flag_draw), "whether to draw the geometry with ROOT's Event Display")
    ("x", po::value<float>()->required(), "initial beam x position")
    ("y", po::value<float>()->required(), "initial beam y position")
    ("energy", po::value<float>()->required(), "beam energy position")
    ("energy_scale", po::value<float>()->default_value(1.f), "factor to scale the energy of the right beam")
    ("width_scale", po::value<float>()->default_value(1.f), "factor to scale the width of both beams")
    ("mass_interaction", po::value<float>()->default_value(0.938), "modelled interaction mass [GeV]")
    ("npartons", po::value<unsigned>()->default_value(1), "number of partons in a proton colliding")
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
  info.energy_scale = boost::any_cast<float>(vm["energy_scale"].value());
  info.width_scale = boost::any_cast<float>(vm["width_scale"].value());
  info.mass = 0.938; //GeV
  info.mass_interaction = boost::any_cast<float>(vm["mass_interaction"].value()); //GeV
  info.npartons = boost::any_cast<unsigned>(vm["npartons"].value()); //GeV
  info.nparticles = boost::any_cast<unsigned>(vm["nparticles"].value());
  info.zcutoff = boost::any_cast<float>(vm["zcutoff"].value());
  assert(info.zcutoff > 0);
  
  run(mode, info);

  if(flag_draw)
    myapp.Run();

  std::cout << std::endl;
  return 0;
}
