#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <string>
#include <vector>

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

#include "TROOT.h"
#include "TMath.h"
#include "Math/Vector3D.h" // XYZVector

//#include "./functions.h"

class Geometry {
 public:
  Geometry() {
    TEveManager::Create();
    this->draw_beam_axis();
  };

 private:
  void draw_beam_axis() {

    TEveLine* TEveLine_beam_axis = NULL;
    TEveLine_beam_axis = new TEveLine();
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,-7650.0);
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,7650.0);
    TEveLine_beam_axis ->SetName("beam axis");
    TEveLine_beam_axis ->SetLineStyle(1);
    TEveLine_beam_axis ->SetLineWidth(4);
    TEveLine_beam_axis ->SetMainAlpha(0.7);
    TEveLine_beam_axis ->SetMainColor(kBlue);
    gEve->AddElement(TEveLine_beam_axis);
    
  }
};

class Magnets: public Geometry {
 public:
  enum Type { DipoleX, DipoleY, Quadrupole, NTYPES };

  struct Magnet {
    Magnets::Type type;
    std::string label;
    int color;
    std::pair<double,double> intensity; //x and y intensity (including sign)
    std::pair<double,double> z; //start and stop in z
    double x;
    double y;
  };
    
 Magnets(const std::vector<Magnet>& pMagnetsInfo)
   : Geometry(), mMagnetsInfo(pMagnetsInfo) {};
  
  void draw() {

    std::vector<TEveBox*> magnets( mMagnetsInfo.size() );

    
    for(unsigned im=0; im<mMagnetsInfo.size(); im++)
      {
	magnets[im] = new TEveBox;
	magnets[im]->SetName( mMagnetsInfo[im].label.c_str() );

	magnets[im]->SetVertex(0, -mMagnetsInfo[im].x, -mMagnetsInfo[im].y, mMagnetsInfo[im].z.first);
	magnets[im]->SetVertex(1,  mMagnetsInfo[im].x, -mMagnetsInfo[im].y, mMagnetsInfo[im].z.first);
	magnets[im]->SetVertex(2,  mMagnetsInfo[im].x,  mMagnetsInfo[im].y, mMagnetsInfo[im].z.first);
	magnets[im]->SetVertex(3, -mMagnetsInfo[im].x,  mMagnetsInfo[im].y, mMagnetsInfo[im].z.first);
	magnets[im]->SetVertex(4, -mMagnetsInfo[im].x, -mMagnetsInfo[im].y, mMagnetsInfo[im].z.second);
	magnets[im]->SetVertex(5,  mMagnetsInfo[im].x, -mMagnetsInfo[im].y, mMagnetsInfo[im].z.second);
	magnets[im]->SetVertex(6,  mMagnetsInfo[im].x,  mMagnetsInfo[im].y, mMagnetsInfo[im].z.second);
	magnets[im]->SetVertex(7, -mMagnetsInfo[im].x,  mMagnetsInfo[im].y, mMagnetsInfo[im].z.second);

	magnets[im]->SetMainColor(mMagnetsInfo[im].color);
	magnets[im]->SetMainTransparency(75); // the higher the value the more transparent
	
	gEve->AddElement(magnets[im]);
      }

  }

  ROOT::Math::XYZVector field(ROOT::Math::XYZVector pos, double scale=1.0) const {
    ROOT::Math::XYZVector Bfield(0.,0.,0.);
    double fieldDir = 1.0;
    
    for(auto && info : mMagnetsInfo)
      {
	if(info.type == Type::DipoleX)
	  {
	    if(pos.Z() < info.z.second and pos.Z() > info.z.first) {
	      if(info.intensity.second != 0)
		std::cout << "Are you sure this is a dipole along x?"<< std::endl;
	      double fieldIntensity = info.intensity.first*fieldDir*scale;
	      Bfield.SetXYZ(fieldIntensity, 0., 0.);
	    }
	  }
	
	else if(info.type == Type::DipoleY)
	  {
	    if(pos.Z() < info.z.second and pos.Z() > info.z.first) {
	      if(info.intensity.first != 0)
		std::cout << "Are you sure this is a dipole along y?"<< std::endl;
	      double fieldIntensity = info.intensity.second*fieldDir*scale;
	      fieldIntensity *= pos.Z()<0 ? -1.f : 1.f;
	      Bfield.SetXYZ(0., fieldIntensity, 0.);
	    }
	  }

	else if(info.type == Type::Quadrupole)
	  {
	    if(pos.Z() < info.z.second and pos.Z() > info.z.first) {
	      if(info.intensity.second == 0 or info.intensity.first == 0)
		std::cout << "Are you sure this is a quadrupole?"<< std::endl;
	      double fieldIntensityX = info.intensity.first*fieldDir*scale*pos.Y();
	      double fieldIntensityY = info.intensity.first*fieldDir*scale*pos.X();
	      const float sign = get_sign_direction(pos.Z());
	      fieldIntensityX *= sign;
	      fieldIntensityY *= sign;
	      Bfield.SetXYZ(fieldIntensityX, fieldIntensityY, 0.);
	    }
	  }
      }
    return Bfield;
  }

 private:
  static constexpr std::pair<double,double> magnetSize = std::make_pair(100.0,100.0);
  std::vector<Magnet> mMagnetsInfo;

  const float get_sign_direction(double z) const {
    return z<0 ? -1.f : 1.f;
  }
};

#endif //GEOMETRY_H
