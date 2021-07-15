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
  class Dimensions {
  public:
    using pair = std::pair<double, double>;
    ROOT::Math::XYZVector beg;
    ROOT::Math::XYZVector end;

    Dimensions(double x1, double x2,
	       double y1, double y2,
	       double z1, double z2) {
      beg.SetXYZ(x1, y1, z1);
      end.SetXYZ(x2, y2, z2);
      
      mX = std::make_pair( beg.X(), end.X() );
      mY = std::make_pair( beg.Y(), end.Y() );
      mZ = std::make_pair( beg.Z(), end.Z() );
    }
      
    const pair X() const { return mX; }
    const pair Y() const { return mY; }
    const pair Z() const { return mZ; }

  private:
    pair mX, mY, mZ;
  };

  Geometry() {
    TEveManager::Create();
    this->draw_beam_axis();
    this->draw_beams();
  };

  //virtual void draw();
  
private:
  void draw_beam_axis() {

    TEveLine* TEveLine_beam_axis = NULL;
    TEveLine_beam_axis = new TEveLine();
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,-6300-5840.);
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,1000);
    TEveLine_beam_axis ->SetName("beam axis");
    TEveLine_beam_axis ->SetLineStyle(9);
    TEveLine_beam_axis ->SetLineWidth(1);
    TEveLine_beam_axis ->SetMainAlpha(0.7);
    TEveLine_beam_axis ->SetMainColor(kBlue);
    gEve->AddElement(TEveLine_beam_axis);
    
  }

  void draw_beams() {

    const unsigned nbeams = 2;
    std::array<std::string,nbeams> names_ = { {"incoming beam",
					       "outgoing beam"} };
    std::array<Dimensions,nbeams> coords_{{ Dimensions{9.3,  0., 0., 0., -6300-5840., 0.},
					    Dimensions{-9.3, 0., 0., 0., -6300-5840., 0.} }};
    std::array<unsigned,nbeams> colors_ = {{kGreen-9,kGreen-1}};
    for(unsigned i=0; i<2; ++i) {
      TEveLine* TEveLine_beam_axis = NULL;
      TEveLine_beam_axis = new TEveLine();
      TEveLine_beam_axis ->SetNextPoint( coords_[i].X().first,  coords_[i].Y().first,  coords_[i].Z().first  );
      TEveLine_beam_axis ->SetNextPoint( coords_[i].X().second, coords_[i].Y().second, coords_[i].Z().second );
      TEveLine_beam_axis ->SetName( names_[i].c_str() );
      TEveLine_beam_axis ->SetLineStyle(1);
      TEveLine_beam_axis ->SetLineWidth(1);
      TEveLine_beam_axis ->SetMainAlpha(0.8);
      TEveLine_beam_axis ->SetMainColor( colors_[i] );
      gEve->AddElement(TEveLine_beam_axis);
    }
  }
};

class Magnets: public Geometry {
public:
  using XYZ = ROOT::Math::XYZVector;
  enum Type { DipoleX, DipoleY, Quadrupole, NTYPES };

  struct Magnet {
    Magnets::Type type; // one of the possibilities in the enum above
    std::string label; // name to be shown in the event display
    int color;
    std::pair<double,double> intensity; //B field intensity along x and y (with sign) [T]
    Geometry::Dimensions dims; //beginning and end coordinates (x, y and z) [cm]
  };
    
  Magnets(const std::vector<Magnet>& pMagnetsInfo)
    : Geometry(), mMagnetsInfo(pMagnetsInfo) {};
  
  void draw();
  XYZ field(XYZ, double) const;

private:
  std::vector<Magnet> mMagnetsInfo;

  const float get_sign_direction(double z) const {
    return z<0 ? -1.f : 1.f;
  }
};


class Calorimeters: public Geometry {
public:
  enum Type { Neutron, Proton, NTYPES };

  struct Calorimeter {
    Calorimeters::Type type; // one of the possibilities in the enum above
    std::string label; // name to be shown in the event display
    int color;
    Geometry::Dimensions dims; //beginning and end coordinates (x, y and z) [cm]
  };
    
  Calorimeters(const std::vector<Calorimeter>& pCalosInfo)
    : Geometry(), mCalosInfo(pCalosInfo) {};
  
  void draw() const;

private:
  std::vector<Calorimeter> mCalosInfo;
};

#endif //GEOMETRY_H
