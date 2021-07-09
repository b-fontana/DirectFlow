#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <string>
#include <vector>

#include <TEveManager.h>
#include "TEveLine.h"
#include "TEveBox.h"

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
  using XYZ = ROOT::Math::XYZVector;
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
  
  void draw();
  XYZ field(XYZ, double) const;

private:
  static constexpr std::pair<double,double> magnetSize = std::make_pair(100.0,100.0);
  std::vector<Magnet> mMagnetsInfo;

  const float get_sign_direction(double z) const {
    return z<0 ? -1.f : 1.f;
  }
};

#endif //GEOMETRY_H
