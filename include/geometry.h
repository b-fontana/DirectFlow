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

////////////////////////////////
/////Dimensions/////////////////
////////////////////////////////
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

  const float X1() const {return mX.first;}
  const float X2() const {return mX.second;}
  const float Y1() const {return mY.first;}
  const float Y2() const {return mY.second;}
  const float Z1() const {return mZ.first;}
  const float Z2() const {return mZ.second;}

private:
  pair mX, mY, mZ;
};

////////////////////////////////
/////Calorimeter////////////////
////////////////////////////////
class Calo {
public:
  enum Type { Neutron, Proton, NTYPES };
  
  Type type; // one of the possibilities in the enum above
  std::string label; // name to be shown in the event display
  int color;
  Dimensions dims; //beginning and end coordinates (x, y and z) [cm]
};

////////////////////////////////
/////Magnet/////////////////////
////////////////////////////////
class Magnet {
public:
  enum Type { DipoleX, DipoleY, Quadrupole, NTYPES };
  
  Type type; // one of the possibilities in the enum above
  std::string label; // name to be shown in the event display
  int color;
  std::pair<double,double> intensity; //B field intensity along x and y (with sign) [T]
  Dimensions dims; //beginning and end coordinates (x, y and z) [cm]
};
  
////////////////////////////////
/////Group of Magnets///////////
////////////////////////////////
class MagnetSystem {
public:
  using XYZ = ROOT::Math::XYZVector;
    
  MagnetSystem(const std::vector<Magnet>& pMagnets)
    : mMagnets(pMagnets) {};
  
  void draw() const;
  XYZ field(XYZ, double) const;

  const float get_sign_direction(double z) const {
    return z<0 ? -1.f : 1.f;
  }
  
private:
  std::vector<Magnet> mMagnets;
};

////////////////////////////////
/////Group of Calorimeters//////
////////////////////////////////
class CaloSystem {
public:
    
  CaloSystem(const std::vector<Calo>& pCalos)
    : mCalos(pCalos) {};
  
  void draw() const;

private:
  std::vector<Calo> mCalos;
};


class BuildGeom {
public:

  BuildGeom(Dimensions ax,
	    const MagnetSystem& pMagnetSyst,
	    const CaloSystem& pCaloSyst)
    : mMSyst(new MagnetSystem(pMagnetSyst)),
      mCSyst(new CaloSystem(pCaloSyst))
  {
    init_(ax);
    mMSyst->draw();
    mCSyst->draw();
  };

  BuildGeom(Dimensions ax,
	    const MagnetSystem& pMagnetSyst)
    : mMSyst(new MagnetSystem(pMagnetSyst))
  {
    init_(ax);
    mMSyst->draw();
  };
  
  BuildGeom(Dimensions ax,
	    const CaloSystem& pCaloSyst)
    : mCSyst(new CaloSystem(pCaloSyst))
  {
    init_(ax);
    mCSyst->draw();
  };

  void init_(const Dimensions& ax) {
    TEveManager::Create();
    draw_beam_axis_(ax);
  }
  
  void draw_beam_axis_(Dimensions ax) {
    TEveLine* bax = new TEveLine();
    bax->SetNextPoint( ax.X1(), ax.Y1(), ax.Z1() );
    bax->SetNextPoint( ax.X2(), ax.Y2(), ax.Z2() );
    bax->SetName("beam_axis");
    bax->SetLineStyle(9);
    bax->SetLineWidth(1);
    bax->SetMainAlpha(0.7);
    bax->SetMainColor(kBlue);
    gEve->AddElement(bax);
  }
  
private:
  std::unique_ptr<MagnetSystem> mMSyst;
  std::unique_ptr<CaloSystem> mCSyst;
};

#endif //GEOMETRY_H

  // void draw_beams_() {

  //   const unsigned nbeams = 2;
  //   std::array<std::string,nbeams> names_ = { {"incoming beam",
  // 					       "outgoing beam"} };
  //   std::array<Dimensions,nbeams> coords_{{ Dimensions{9.3,  0., 0., 0., -6300-5840., 0.},
  // 					    Dimensions{-9.3, 0., 0., 0., -6300-5840., 0.} }};
  //   std::array<unsigned,nbeams> colors_ = {{kGreen-9,kGreen-1}};
  //   for(unsigned i=0; i<2; ++i) {
  //     TEveLine* TEveLine_beam_axis = NULL;
  //     TEveLine_beam_axis = new TEveLine();
  //     TEveLine_beam_axis ->SetNextPoint( coords_[i].X().first,  coords_[i].Y().first,  coords_[i].Z().first  );
  //     TEveLine_beam_axis ->SetNextPoint( coords_[i].X().second, coords_[i].Y().second, coords_[i].Z().second );
  //     TEveLine_beam_axis ->SetName( names_[i].c_str() );
  //     TEveLine_beam_axis ->SetLineStyle(1);
  //     TEveLine_beam_axis ->SetLineWidth(1);
  //     TEveLine_beam_axis ->SetMainAlpha(0.8);
  //     TEveLine_beam_axis ->SetMainColor( colors_[i] );
  //     gEve->AddElement(TEveLine_beam_axis);
  //   }
  // }
