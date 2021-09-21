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

struct InputArgs {
public:
  bool draw;
  double energy;
  double mass;
  unsigned nparticles;
};

double calculate_intersection_distance(const Vec<XYZ>& p1, const Vec<XYZ>& p2,
				       unsigned s, double t, double stepSize) {

  unsigned nMatches(0);
  for(unsigned i(0); i<s; ++i) 
    {
      if( local_distance(p1[i], p2[i]) < t )
	++nMatches;
    }
  return static_cast<double>(nMatches)*stepSize;
}

unsigned size_last_batch(unsigned nbatches, unsigned nelems, unsigned batchSize) {
  return nelems-(nbatches-1)*batchSize;
}

void run(const InputArgs& args)
{
  using XYZ = ROOT::Math::XYZVector;

  gStyle->SetPalette(56); // 53 = black body radiation, 56 = inverted black body radiator, 103 = sunset, 87 == light temperature

  //set global variables
  std::fstream file;
  XYZ origin(0., 0., 0.);
  unsigned nsteps = 30000;
  double stepsize = 0.001;

  double dist = 5.;
  double energy = args.energy;
  double mass = args.mass;
    
  //generate random positions around input positions
  UniformDistribution<double> thetadist(0, M_PI/2);
  UniformDistribution<double> phidist(-M_PI, M_PI);
  
  //Vec<Magnets::Magnet> magnetInfo{};
  MagnetSystem magnets({});

  //Vec<Calo> caloInfo{};
  CaloSystem calos({});

  if(args.draw) {
    float beamcap = dist;
    BuildGeom(Dimensions{0., 0., 0., 0., -beamcap, beamcap}, //beamline coordinates
	      magnets, calos);
  }
  
  constexpr float batchSize = 500.f;
  const unsigned nbatches = ceil(args.nparticles/batchSize);
  
  std::cout << " --- Simulation Information --- " << std::endl;
  std::cout << "Batch Size: " << batchSize << " (last batch: " << size_last_batch(nbatches, args.nparticles, batchSize) << ")" << std::endl;
  std::cout << "Number of batches: " << nbatches << std::endl;
  std::cout << "Step Size: " << stepsize << std::endl;
  std::cout << "--------------------------" << std::endl;
  unsigned batchSize_;

  std::string filename("data/collision_prob.csv");
  file.open(filename, std::ios_base::out);
      
  for (unsigned ibatch : tq::trange(nbatches))
    {
      batchSize_ = ibatch==nbatches-1 ? size_last_batch(nbatches, args.nparticles, batchSize) : batchSize;
      
      //define the initial properties of the incident particle
      Vec<Particle> p_fixed(batchSize_);
      Vec<Particle> p_moving(batchSize_);

      Vec<double> Thetas(batchSize_);
      Vec<double> Phis(batchSize_);
      
      for(unsigned i=0; i<batchSize_; ++i) {
	//negative z side
	p_fixed[i].pos = XYZ( 0., 0., -dist ); //cm
	p_fixed[i].mom = XYZ(0.0, 0.0, TMath::Sqrt(energy*energy - mass*mass) ); // GeV/c
	p_fixed[i].mass = args.mass; // GeV/c^2
	p_fixed[i].energy = energy;
	p_fixed[i].charge = +1;

	//positive z side
        Thetas[i] = thetadist.generate();
	Phis[i]   = phidist.generate();
	p_moving[i].pos = Polar(dist, Thetas[i], Phis[i]); //cm

	double mom_magnitude = TMath::Sqrt(energy*energy - mass*mass);
	p_moving[i].mom = -1*p_moving[i].pos.Unit()*mom_magnitude; // GeV/c

	p_moving[i].mass = args.mass; // GeV/c^2
	p_moving[i].energy = energy;
	p_moving[i].charge = +1;
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
	simp1.push_back( SimParticle(p_fixed[i], nsteps, stepsize) );
	simp2.push_back( SimParticle(p_moving[i], nsteps, stepsize) );
      }

      Vec<Track> tracks1(batchSize_);
      Vec<Track> tracks2(batchSize_);
      for(unsigned i=0; i<batchSize_; ++i) {
	tracks1[i] = simp1[i].track_straight();
	tracks2[i] = simp2[i].track_straight();
      }

      Vec<double> intDist(batchSize_);
            
      for(unsigned i=0; i<batchSize_; ++i) {
	//negative z side
	Vec<double> en1 = tracks1[i].energies();
	Vec<XYZ> pos1   = tracks1[i].positions();
	Vec<XYZ> mom1   = tracks1[i].momenta();
	unsigned nStepsUsed1  = tracks1[i].steps_used();

	//positive z side
	Vec<double> en2 = tracks2[i].energies();
	Vec<XYZ> pos2   = tracks2[i].positions();
	Vec<XYZ> mom2   = tracks2[i].momenta();
	unsigned nStepsUsed2  = tracks2[i].steps_used();

	unsigned minelem = std::min(nStepsUsed1, nStepsUsed2);

	for(unsigned i_step = 0; i_step<minelem; i_step++)
	  {
	    
	    if(args.draw) {
	      particleTrackViz1[i]->SetNextPoint(pos1[i_step].X(),
						 pos1[i_step].Y(),
						 pos1[i_step].Z() );
	      
	      particleTrackViz2[i]->SetNextPoint(pos2[i_step].X(),
						 pos2[i_step].Y(),
						 pos2[i_step].Z() );
	    }
	  }

	intDist[i] = calculate_intersection_distance(pos1, pos2, minelem, 0.01, stepsize);
  	
	if(i==0 and ibatch==0)
	  file << "iBatch,Idx,intDist,Theta,Phi" << std::endl;
		
	file << std::to_string( ibatch ) << ","
	     << std::to_string( i ) << ","
	     << std::to_string( intDist[i] ) << ","
	     << std::to_string( Thetas[i] ) << ","
	     << std::to_string( Phis[i] )
	     << std::endl;	
      }
      

      if(args.draw) {
	for(unsigned ix=0; ix<batchSize_; ix++) {
	  std::string histname1 = "track_fixed_ " + std::to_string(ix);
	  particleTrackViz1[ix]->SetName( histname1.c_str() );
	  particleTrackViz1[ix]->SetLineStyle(1);
	  particleTrackViz1[ix]->SetLineWidth(2);
	  particleTrackViz1[ix]->SetMainAlpha(0.7);
	  particleTrackViz1[ix]->SetMainColor(kRed+3);

	  std::string histname2 = "track_moving_ " + std::to_string(ix);
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

  file.close();

  if(args.draw)
    gEve->Redraw3D(kTRUE);
}

// run example: ./collision_prob --energy 1380 --nparticles 1
int main(int argc, char **argv) {
  TApplication myapp("myapp", &argc, argv);
  
  bool flag_draw = false;
 
  namespace po = boost::program_options;
  po::options_description desc("Options");
  //https://www.boost.org/doc/libs/1_45_0/doc/html/boost/program_options/typed_value.html#id903171-bb
  desc.add_options()
    ("help,h", "produce this help message")
    ("draw", po::bool_switch(&flag_draw), "whether to draw the geometry with ROOT's Event Display")
    ("nparticles", po::value<unsigned>()->default_value(1), "number of particles to generate on each beam")
    ("energy", po::value<float>()->required(), "beam energy position");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc), vm);
  po::notify(vm);

  if(vm.count("help") or argc<2) {
    std::cerr << desc << std::endl;
    std::exit(0);
  }
      
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
  info.nparticles = boost::any_cast<unsigned>(vm["nparticles"].value());
  info.energy = boost::any_cast<float>(vm["energy"].value());
  info.mass = 0.938; //GeV
  
  run(info);

  if(flag_draw)
    myapp.Run();

  std::cout << std::endl;
  return 0;
}
