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

#include <stdio.h>

struct InputArgs {
public:
  long double energy;
  long double mass;
  unsigned nparticles;
};

long double calculate_intersection_distance(const Vec<XYZ>& p1, const Vec<XYZ>& p2,
				       unsigned s, long double t, long double stepSize) {

  unsigned nMatches(0);
  for(unsigned i(0); i<s; ++i) 
    {
      if( local_distance(p1[i], p2[i]) < t )
	++nMatches;
    }
  return static_cast<long double>(nMatches);
}

unsigned size_last_batch(unsigned nbatches, unsigned nelems, unsigned batchSize) {
  return nelems-(nbatches-1)*batchSize;
}

long double track_straight(XYZ pos1, XYZ mom1,
			XYZ pos2, XYZ mom2,
			long double step, unsigned nSteps)
{
  XYZ partPos1 = pos1; // cm
  XYZ partMom1 = mom1; // Gev/c
  XYZ partPos2 = pos2; // cm
  XYZ partMom2 = mom2; // Gev/c

  long double threshold = .1;
  long double thisStepSize = 0.01;
  long double distClose = 0.;
  
  Vec<long double> dists;

  unsigned nStepsUsed(0);
  while(nStepsUsed<nSteps)
    {
      //long double energy = TMath::Sqrt(partMom.Mag2() + mParticle.mass*mParticle.mass);
      
      XYZ posIncr1 = partMom1.Unit() * thisStepSize;
      XYZ posIncr2 = partMom2.Unit() * thisStepSize;
      if( TMath::Sqrt(posIncr1.Mag2())-thisStepSize>1e-5 or
	  TMath::Sqrt(posIncr2.Mag2())-thisStepSize>1e-5 )
	std::cout << "ERROR: " <<  TMath::Sqrt(posIncr1.Mag2()) << ", " <<  thisStepSize << ", " << TMath::Sqrt(posIncr2.Mag2()) << ", " <<  thisStepSize << std::endl;
	
      partPos1 += posIncr1;
      partPos2 += posIncr2;

      long double dist = local_distance(partPos1, partPos2);
      if(dist < threshold)
	distClose += thisStepSize;
      
      if( dist > 2*1.1 )
      	break;
      else if( dist < threshold + step )
	thisStepSize = step / 10000.;
      else
	thisStepSize = step;

      if(nStepsUsed%100==0)
	std::cout << dist << ", " << thisStepSize << " || ";
      ++nStepsUsed;
    }
  std::cout << std::endl;
  return distClose;
}

//assumes the velocity of oth particles is the same
//assumes particle1 flies along the z direction (y=0 always)
long double solve_order2_distance_intersection(long double z1, long double z2, long double y2,
					  long double velocity, long double angle,
					  long double threshold) {
  
  auto afunc = [velocity, angle]() {
		 long double vel2 = velocity*velocity;
		 return 2*(1+std::cos(angle))*vel2;
	       };
s
  auto bfunc = [z1, z2, y2, velocity, angle]() {
		 long double zdiff = z1-z2;
		 return 2*velocity*(zdiff*(1+std::cos(angle))-y2*std::sin(angle));
	       };

  auto cfunc = [z1, z2, y2, threshold]() {
		 long double zdiff = z1-z2;
		 return zdiff*zdiff + y2*y2 - threshold*threshold;
	       };

  long double aval = afunc();
  long double bval = bfunc();
  long double cval = cfunc();
  long double sqrtval = std::sqrt(bval*bval - 4*aval*cval);

  long double time1 = ( -bval + sqrtval ) / ( 2 * aval );
  long double time2 = ( -bval - sqrtval ) / ( 2 * aval );

  return std::abs(time1 - time2);
}

void run(const InputArgs& args)
{
  using XYZ = ROOT::Math::XYZVector;

  //set global variables
  std::fstream file;
  XYZ origin(0., 0., 0.);

  long double dist = 1.1;
  long double energy = args.energy;
  long double mass = args.mass;
  //unsigned nsteps = 500000;
  long double stepsize = 0.001;
    
  //generate random positions around input positions
  long double nomAngle = 2 * TMath::ASin(0.1/5000);
  //UniformDistribution<long double> thetadist(nomAngle-nomAngle/10, nomAngle+nomAngle/10);
  UniformDistribution<long double> thetadist(nomAngle, nomAngle);

  const long double mom_magnitude = TMath::Sqrt(energy*energy - mass*mass);
  
  constexpr long double batchSize = 1500.;
  const unsigned nbatches = ceil(args.nparticles/batchSize);
  
  std::cout << " --- Simulation Information --- " << std::endl;
  std::cout << "Batch Size: " << batchSize << " (last batch: "
	    << size_last_batch(nbatches, args.nparticles, batchSize)
	    << ")" << std::endl;
  std::cout << "Number of batches: " << nbatches << std::endl;
  std::cout << "Step Size: " << stepsize << std::endl;
  std::cout << "--------------------------" << std::endl;
  unsigned batchSize_;

  std::string filename("data/collision_prob.csv");
  file.open(filename, std::ios_base::out);

  long double nomDistProxy = solve_order2_distance_intersection(-dist,
								dist*std::cos(nomAngle),
								dist*std::sin(nomAngle),
								mom_magnitude/mass,
								nomAngle,
								0.01);
  
  for (unsigned ibatch : tq::trange(nbatches))
    {
      batchSize_ = ibatch==nbatches-1 ? size_last_batch(nbatches, args.nparticles, batchSize) : batchSize;
      
      //define the initial properties of the incident particle
      Vec<Particle> p_fixed(batchSize_);
      Vec<Particle> p_moving(batchSize_);

      Vec<long double> Thetas(batchSize_);
	
      for(unsigned i=0; i<batchSize_; ++i) {
	//negative z side
	p_fixed[i].pos = XYZ( 0., 0., -dist ); //cm
	p_fixed[i].mom = XYZ(0.0, 0.0, mom_magnitude ); // GeV/c
	p_fixed[i].mass = args.mass; // GeV/c^2
	p_fixed[i].energy = energy;
	p_fixed[i].charge = +1;

	//positive z side
        Thetas[i] = thetadist.generate();
	p_moving[i].pos = XYZ(0., dist*std::sin(Thetas[i]), dist*std::cos(Thetas[i])); //cm
	p_moving[i].mom = -1*p_moving[i].pos.Unit()*mom_magnitude; // GeV/c

	p_moving[i].mass = args.mass; // GeV/c^2
	p_moving[i].energy = energy;
	p_moving[i].charge = +1;

	// long double distProxy = track_straight(p_fixed[i].pos, p_fixed[i].mom,
	// 				  p_moving[i].pos, p_moving[i].mom,
	// 				  stepsize, nsteps);

	long double distProxy = solve_order2_distance_intersection(p_fixed[i].pos.Z(),
							      p_moving[i].pos.Z(),
							      p_moving[i].pos.Y(),
							      mom_magnitude/mass,
							      Thetas[i],
							      0.01); //threshold [cm]


	printf("%.20Lf\n", distProxy);
	distProxy = nomDistProxy-distProxy;
	//distProxy *= 10000;
	
	//std::cout << distProxy << ", " << std::abs(nomDistProxy-distProxy) << std::endl;
	
	if(i==0 and ibatch==0)
	  file << "iBatch,Idx,distProxy,Theta" << std::endl;

	char str1[40], str2[40];
	sprintf(str1, "%.20Lf", distProxy);
	sprintf(str2, "%.20Lf", Thetas[i]);
	
	file << std::to_string( ibatch ) << ","
	     << std::to_string( i ) << ","
	  //	     << std::to_string( distProxy ) << ","
	     << str1 << ","
	     << str2
	     << std::endl;	
      }
      

    } // for ibatch

  file.close();
}

// run example: ./collision_prob --energy 1380 --nparticles 1
int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description desc("Options");
  //https://www.boost.org/doc/libs/1_45_0/doc/html/boost/program_options/typed_value.html#idsthetadist903171-bb
  desc.add_options()
    ("help,h", "produce this help message")
    ("nparticles", po::value<unsigned>()->default_value(1), "number of particles to generate on each beam")
    ("energy", po::value<long double>()->required(), "beam energy position");
      
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
    if (auto v = boost::any_cast<long double>(&value))
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
  info.nparticles = boost::any_cast<unsigned>(vm["nparticles"].value());
  info.energy = boost::any_cast<long double>(vm["energy"].value());
  info.mass = 0.938; //GeV
  
  run(info);

  std::cout << std::endl;
  return 0;
}
