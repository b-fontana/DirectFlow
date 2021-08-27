#include "include/geometry.h"

void MagnetSystem::draw() const {

  std::vector<TEveBox*> magnets( mMagnets.size() );

    
  for(unsigned im=0; im<mMagnets.size(); im++)
    {
      magnets[im] = new TEveBox;
      magnets[im]->SetName( mMagnets[im].label.c_str() );

      Dimensions d = mMagnets[im].dims;
      magnets[im]->SetVertex(0, d.X().first,  d.Y().first,  d.Z().first);
      magnets[im]->SetVertex(1, d.X().second, d.Y().first,  d.Z().first);
      magnets[im]->SetVertex(2, d.X().second, d.Y().second, d.Z().first);
      magnets[im]->SetVertex(3, d.X().first,  d.Y().second, d.Z().first);
      magnets[im]->SetVertex(4, d.X().first,  d.Y().first,  d.Z().second);
      magnets[im]->SetVertex(5, d.X().second, d.Y().first,  d.Z().second);
      magnets[im]->SetVertex(6, d.X().second, d.Y().second, d.Z().second);
      magnets[im]->SetVertex(7, d.X().first,  d.Y().second, d.Z().second);

      magnets[im]->SetMainColor(mMagnets[im].color);
      magnets[im]->SetMainTransparency(75); // the higher the value the more transparent
	
      gEve->AddElement(magnets[im]);
    }

}

MagnetSystem::XYZ MagnetSystem::field(XYZ pos, double scale=1.0) const {
  
  XYZ Bfield(0.,0.,0.);
  const double fieldDir = 1.0;
    
  for(auto && info : mMagnets)
    {
      if(info.type == Magnet::DipoleX)
	{
	  if(pos.Z() < info.dims.Z().second and pos.Z() > info.dims.Z().first) {
	    if(info.intensity.second != 0)
	      std::cout << "Are you sure this is a dipole along x?"<< std::endl;
	    double fieldIntensity = info.intensity.first*fieldDir*scale;
	    Bfield.SetXYZ(fieldIntensity, 0., 0.);
	  }
	}
	
      else if(info.type == Magnet::DipoleY)
	{
	  if(pos.Z() < info.dims.Z().second and pos.Z() > info.dims.Z().first) {
	    if(info.intensity.first != 0)
	      std::cout << "Are you sure this is a dipole along y?"<< std::endl;
	    double fieldIntensity = info.intensity.second*fieldDir*scale;
	    fieldIntensity *= pos.Z()<0 ? -1.f : 1.f;
	    Bfield.SetXYZ(0., fieldIntensity, 0.);
	  }
	}

      else if(info.type == Magnet::Quadrupole)
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

void CaloSystem::draw() const {

  std::vector<TEveBox*> calos( mCalos.size() );

    
  for(unsigned im=0; im<mCalos.size(); im++)
    {
      calos[im] = new TEveBox;
      calos[im]->SetName( mCalos[im].label.c_str() );

      Dimensions d = mCalos[im].dims;
      calos[im]->SetVertex(0, d.X().first,  d.Y().first,  d.Z().first);
      calos[im]->SetVertex(1, d.X().second, d.Y().first,  d.Z().first);
      calos[im]->SetVertex(2, d.X().second, d.Y().second, d.Z().first);
      calos[im]->SetVertex(3, d.X().first,  d.Y().second, d.Z().first);
      calos[im]->SetVertex(4, d.X().first,  d.Y().first,  d.Z().second);
      calos[im]->SetVertex(5, d.X().second, d.Y().first,  d.Z().second);
      calos[im]->SetVertex(6, d.X().second, d.Y().second, d.Z().second);
      calos[im]->SetVertex(7, d.X().first,  d.Y().second, d.Z().second);

      calos[im]->SetMainColor(mCalos[im].color);
      calos[im]->SetMainTransparency(75); // the higher the value the more transparent
	
      gEve->AddElement(calos[im]);
    }

}
