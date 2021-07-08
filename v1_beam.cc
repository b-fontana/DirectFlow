//#include "./functions.h"
#include "include/geometry.h"
#include "include/tracking.h"
//#include "./functions_proton_spec.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#include <TApplication.h>
#include <TEveManager.h>
#include "TEveLine.h"
#include "TEveBox.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"

#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLCamera.h"
#include "TGLPerspectiveCamera.h"
#include "TGFrame.h"
#include "TGLUtil.h"
#include "TGLLightSet.h"
#include "TGLCameraOverlay.h"

#include "Math/Vector3D.h" // XYZVector

#include "THnSparse.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"
#include "TVector.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFitResult.h"
#include "TList.h"
#include "TChain.h"
//#include "Math/GSLMinimizer.h"
//#include "Math/GSLSimAnMinimizer.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "Math/Functor.h"
#include "TMinuitMinimizer.h"
//#include <GSLMultiMinimizer.h>
//#include "Math/GSLMultiMinimizer.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "TRotation.h"
#include "TSVDUnfold.h"
#include "TSystemDirectory.h"
#include "TMatrixT.h"
#include "TVector2.h"
#include "TEllipse.h"
#include "TPaveLabel.h"
#include "TBox.h"
#include "TSystemDirectory.h"
#include "TAttImage.h"
#include "TImage.h"
#include "TPaletteAxis.h"
#include "THelix.h"
#include "TView.h"
#include "TSystem.h"
#include "TASImage.h"
#include "TImage.h"
#include "TArrayD.h"
//#include "TPython.h"
#include "TArrow.h"
#include "TKey.h"
#include "TSpectrum.h"
#include "TNtuple.h"
#include "TDatime.h"

void run()
{
  gRandom->SetSeed(0);
  gStyle->SetPalette(56); // 53 = black body radiation, 56 = inverted black body radiator, 103 = sunset, 87 == light temperature

  //define the initial properties of the incident particle
  Particle particle;
  particle.pos = ROOT::Math::XYZVector(0.0, 0.0, -6900.0);
  particle.mom = ROOT::Math::XYZVector(0.0, 0.0, 1500.0); // GeV
  particle.mass = 0.938; // GeV
  particle.energy = TMath::Sqrt(particle.mom.Mag2() + particle.mass*particle.mass);
  particle.charge = +1;
  
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

  //define the magnetic field
  std::vector<Magnets::Magnet> magnetInfo{ {Magnets::DipoleY,
					    "D1_neg",
					    kGreen,
					    std::make_pair(0.,-3.529),
					    std::make_pair(particle.pos.Z()-1500, particle.pos.Z()+3000), 600., 600.},
  };

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

  SimParticle simp(particle, 3000000, 1.f);

  std::vector<double> itEnergies = simp.track(magnets, tracking::TrackMode::Euler ).energies();

  std::vector<ROOT::Math::XYZVector> itPositions = simp.track(magnets, tracking::TrackMode::Euler ).positions();
  std::vector<ROOT::Math::XYZVector> itDirections = simp.track(magnets, tracking::TrackMode::Euler ).directions();
  unsigned nStepsUsed = simp.track(magnets, tracking::TrackMode::Euler ).steps_used();
  std::cout << nStepsUsed << std::endl;
    
  std::string filename("track_pos.csv");
  std::fstream file;
  file.open(filename, std::ios_base::out);
  for(unsigned i_step = 0; i_step < nStepsUsed; i_step++)
    {
      particleTrackViz[0]->SetNextPoint( itPositions[i_step].X(), itPositions[i_step].Y(), itPositions[i_step].Z() );

      if (!file.is_open()) 
	std::cout << "failed to open " << filename << '\n';
      else {
	if(i_step==0)
	  file << "x,y,z,energy" << std::endl;
	//if(i_step%100==0)
	file << std::to_string( itPositions[i_step].X() ) << ","
	     << std::to_string( itPositions[i_step].Y() ) << ","
	     << std::to_string( itPositions[i_step].Z() ) << ","
	     << std::to_string( itEnergies[i_step] )
	     << std::endl;
      }
      //f >> std::to_string(track_pos_dir[i_step][0]) >> std::endl;
      //printf("i_step: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_step,track_pos_dir[i_step][0],track_pos_dir[i_step][1],track_pos_dir[i_step][2]);

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
  run();
  myapp.Run();
  return 0;
}

