#include "include/geometry.h"

void Magnets::draw() {

  std::vector<TEveBox*> magnets( mMagnetsInfo.size() );

    
  for(unsigned im=0; im<mMagnetsInfo.size(); im++)
    {
      magnets[im] = new TEveBox;
      magnets[im]->SetName( mMagnetsInfo[im].label.c_str() );

      Geometry::Dimensions d = mMagnetsInfo[im].dims;
      magnets[im]->SetVertex(0, d.X().first,  d.Y().first,  d.Z().first);
      magnets[im]->SetVertex(1, d.X().second, d.Y().first,  d.Z().first);
      magnets[im]->SetVertex(2, d.X().second, d.Y().second, d.Z().first);
      magnets[im]->SetVertex(3, d.X().first,  d.Y().second, d.Z().first);
      magnets[im]->SetVertex(4, d.X().first,  d.Y().first,  d.Z().second);
      magnets[im]->SetVertex(5, d.X().second, d.Y().first,  d.Z().second);
      magnets[im]->SetVertex(6, d.X().second, d.Y().second, d.Z().second);
      magnets[im]->SetVertex(7, d.X().first,  d.Y().second, d.Z().second);

      magnets[im]->SetMainColor(mMagnetsInfo[im].color);
      magnets[im]->SetMainTransparency(75); // the higher the value the more transparent
	
      gEve->AddElement(magnets[im]);
    }

}

Magnets::XYZ Magnets::field(XYZ pos, double scale=1.0) const {
  
  XYZ Bfield(0.,0.,0.);
  const double fieldDir = 1.0;
    
  for(auto && info : mMagnetsInfo)
    {
      if(info.type == Type::DipoleX)
	{
	  if(pos.Z() < info.dims.Z().second and pos.Z() > info.dims.Z().first) {
	    if(info.intensity.second != 0)
	      std::cout << "Are you sure this is a dipole along x?"<< std::endl;
	    double fieldIntensity = info.intensity.first*fieldDir*scale;
	    Bfield.SetXYZ(fieldIntensity, 0., 0.);
	  }
	}
	
      else if(info.type == Type::DipoleY)
	{
	  if(pos.Z() < info.dims.Z().second and pos.Z() > info.dims.Z().first) {
	    if(info.intensity.first != 0)
	      std::cout << "Are you sure this is a dipole along y?"<< std::endl;
	    double fieldIntensity = info.intensity.second*fieldDir*scale;
	    fieldIntensity *= pos.Z()<0 ? -1.f : 1.f;
	    Bfield.SetXYZ(0., fieldIntensity, 0.);
	  }
	}

      else if(info.type == Type::Quadrupole)
	{
	  if(pos.Z() < info.dims.Z().first and pos.Z() > info.dims.Z().second) {
	    if(info.intensity.second == 0 or info.intensity.first == 0)
	      std::cout << "Are you sure this is a quadrupole?"<< std::endl;

	    // dividing by 100 to convert from centimeters to meters (assuming coordinates were given in cm)
	    double fieldIntensityX = info.intensity.first*fieldDir*scale*pos.Y() / 100.;
	    double fieldIntensityY = info.intensity.second*fieldDir*scale*pos.X() / 100.;
	    
	    const float sign = get_sign_direction(pos.Z());
	    fieldIntensityX *= sign;
	    fieldIntensityY *= sign;
	    Bfield.SetXYZ(fieldIntensityX, fieldIntensityY, 0.);
	  }
	}
    }
  return Bfield;
}

void Calorimeters::draw() const {

  std::vector<TEveBox*> calos( mCalosInfo.size() );

    
  for(unsigned im=0; im<mCalosInfo.size(); im++)
    {
      calos[im] = new TEveBox;
      calos[im]->SetName( mCalosInfo[im].label.c_str() );

      Geometry::Dimensions d = mCalosInfo[im].dims;
      calos[im]->SetVertex(0, d.X().first,  d.Y().first,  d.Z().first);
      calos[im]->SetVertex(1, d.X().second, d.Y().first,  d.Z().first);
      calos[im]->SetVertex(2, d.X().second, d.Y().second, d.Z().first);
      calos[im]->SetVertex(3, d.X().first,  d.Y().second, d.Z().first);
      calos[im]->SetVertex(4, d.X().first,  d.Y().first,  d.Z().second);
      calos[im]->SetVertex(5, d.X().second, d.Y().first,  d.Z().second);
      calos[im]->SetVertex(6, d.X().second, d.Y().second, d.Z().second);
      calos[im]->SetVertex(7, d.X().first,  d.Y().second, d.Z().second);

      calos[im]->SetMainColor(mCalosInfo[im].color);
      calos[im]->SetMainTransparency(75); // the higher the value the more transparent
	
      gEve->AddElement(calos[im]);
    }

}
