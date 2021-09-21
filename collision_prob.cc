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
};

unsigned size_last_batch(unsigned nbatches, unsigned nelems, unsigned batchSize) {
  return nelems-(nbatches-1)*batchSize;
}

void run(tracking::TrackMode mode, const InputArgs& args)
{
  using XYZ = ROOT::Math::XYZVector;
  using Polar = ROOT::Math::Polar3DVector;
  gStyle->SetPalette(56); // 53 = black body radiation, 56 = inverted black body radiator, 103 = sunset, 87 == light temperature

  //set global variables
  std::fstream file;
  std::string histname;
  XYZ origin(0., 0., 0.);
  unsigned nsteps = 30000;
  double stepsize = 0.1;

  double distance = 10.;
  double energy = 1380.;
  double mass = 0.938;
    
  //generate random positions around input positions
  UniformDistribution<double> thetadist(0, M_PI);
  UniformDistribution<double> phidist(-M_PI, M_PI);
  
  //Vec<Magnets::Magnet> magnetInfo{};
  MagnetSystem magnets({});

  //Vec<Calo> caloInfo{};
  CaloSystem calos({});

  if(args.draw) {
    float beamcap = distance;
    BuildGeom(Dimensions{0., 0., 0., 0., -beamcap, beamcap}, //beamline coordinates
	      magnets, calos);
  }
  
  constexpr float batchSize = 1500.f;
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
      for(unsigned i=0; i<batchSize_; ++i) {
	//negative z side
	p_fixed[i].pos = XYZ( 0., 0., -distance ); //cm
	p_fixed[i].mom = XYZ(0.0, 0.0, TMath::Sqrt(energy*energy - mass*mass) ); // GeV/c
	p_fixed[i].mass = 0.938; // GeV/c^2
	p_fixed.energy = energy;
	p_fixed[i].charge = +1;

	//positive z side
	double theta = thetadist.generate();
	double phi = phidist.generate();
	p_moving[i].pos = Polar(distance, theta, phi); //cm

	double mom_magnitude = TMath::Sqrt(energy*energy - mass*mass);
	p_moving[i].mom = -1*p_moving[i].pos.Unit()*mom_magnitude; // GeV/c

	p_moving[i].mass = 0.938; // GeV/c^2
	p_moving.energy = energy;
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
	simp1.push_back( SimParticle(p_fixed[i], nsteps[mode], stepsize[mode]) );
	simp2.push_back( SimParticle(p_moving[i], nsteps[mode], stepsize[mode]) );
      }

      Vec<const Track*> tracks1(batchSize_);
      Vec<const Track*> tracks2(batchSize_);
      for(unsigned i=0; i<batchSize_; ++i) {
	tracks1[i] = &( simp1[i].track_straight());
	tracks2[i] = &( simp2[i].track_straight());
      }

      Vec2<double> itEnergies1(batchSize_), itEnergies2(batchSize_);
      Vec2<XYZ> itPositions1(batchSize_), itPositions2(batchSize_);
      Vec2<XYZ> itMomenta1(batchSize_), itMomenta2(batchSize_);
            
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

	XYZ last1Pos_ = itPositions1[i].back();
	TVector3 last1PosV_(last1Pos_.X(), last1Pos_.Y(), last1Pos_.Z());
	TVector3 check1(-p1[i].pos.X(), -p1[i].pos.Y(), args.zcutoff);
	if( check1.Angle(last1PosV_) > 1e-7 ) {
	  std::cout << "The trajectory is not as it should!" << std::endl;
	  std::cout << "Angle1: " << check1.Angle(last1PosV_) << std::endl;
	  std::exit(0);
	}

	double last1X_ = last1Pos_.Dot( uX1 );
	double last1Y_ = last1Pos_.Dot( uY1 );
	
	XYZ last2Pos_ = itPositions2[i].back();
	TVector3 last2V(last2Pos_.X(), last2Pos_.Y(), last2Pos_.Z());
	TVector3 check2(-p2[i].pos.X(), -p2[i].pos.Y(), -args.zcutoff);
	if( check2.Angle(last2V) > 1e-7 ) {
	  std::cout << "The trajectory is not as it should!" << std::endl;
	  std::cout << "Angle2: " << check2.Angle(last2V) << std::endl;
	  std::exit(0);
	}
	double last2X_ = last2Pos_.Dot( uX2 );
	double last2Y_ = last2Pos_.Dot( uY2 );

	//hits distribution without fermi boost
	TVector3 last1Det = Globals::distanceToDetector * last1PosV_.Unit();
	xHitNoBoost[i] = last1Det.Dot( uX1 );
	yHitNoBoost[i] = last1Det.Dot( uY1 );

	//fermi momentum correction
	float fermiMom = fermidist.generate();
	Double_t fermiPhi = phidist.generate();
	Double_t fermiTheta = thetadist.generate();
	TVector3 fermiVec;
	fermiVec.SetPtThetaPhi(1.0, fermiTheta, fermiPhi);
	fermiVec *= fermiMom/fermiVec.Mag();
	fermiVec.SetY(fermiVec.Y() + args.fermi_shift);

	XYZ last1Mom_ = itMomenta1[i].back();
	TLorentzVector last1MomLtz_;
	last1MomLtz_.SetPxPyPzE(last1Mom_.X(), last1Mom_.Y(), last1Mom_.Z(), args.energy);

	TLorentzVector fermiVecLtz(fermiVec, TMath::Sqrt(fermiVec.Mag2() + args.mass*args.mass));
	//std::cout << std::endl;
	//print_pos_4D("fermiVecLtz before boost", fermiVecLtz);
	//print_pos_4D("last1MomLtz before boost", last1MomLtz_);
	TVector3 fermiBoost = last1MomLtz_.BoostVector();
	//print_pos("fermiBoost", fermiBoost);
	fermiPzBeforeBoost[i] = fermiVecLtz.Pz();
	fermiVecLtz.Boost( fermiBoost );
	fermiPzAfterBoost[i] = fermiVecLtz.Pz();
	
	last1Det = Globals::distanceToDetector * (fermiVecLtz.Vect()).Unit();
	xHit[i] = last1Det.Dot( uX1 );
	yHit[i] = last1Det.Dot( uY1 );
	//std::cout << "Generated numbers: FermiPT " << fermiMom << ", Theta " << fermiTheta << ", Phi " << fermiPhi << std::endl;
	// print_pos_4D("fermiVecLtz after boost", fermiVecLtz);
	// print_pos("fermiLtzVector coordinates at detector", last1Det);
	//std::exit(0);
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
	assert(last1Pos_.Z() > last2Pos_.Z());

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
	  file2 << "iBatch,Idx,sumMomX,sumMomY,sumMomZ,FermiPzBeforeBoost,FermiPzAfterBoost,XHitNoBoost,YHitNoBoost,XHit,YHit,PsiA,PsiB,cat1,Psi,Phi,Eta,Cos" << std::endl;
	
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
	      << std::to_string( fermiPzBeforeBoost[ix] ) << ","
	      << std::to_string( fermiPzAfterBoost[ix] ) << ","
	      << std::to_string( xHitNoBoost[ix] ) << ","
	      << std::to_string( yHitNoBoost[ix] ) << ","
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
  info.x = boost::any_cast<float>(vm["x"].value());
  
  run(mode, info);

  if(flag_draw)
    myapp.Run();

  std::cout << std::endl;
  return 0;
}
