---
layout: post
title:  "General Informations"
date:   2021-07-20 10:00:18 +0200
categories: plots
---

The protons are fired in the Z direction with opposited momenta (same magnitude) at a distance of 70 meters with 1.38GeV.

The magnets are not perfectly symmetric, hence some small misalignment:

```c++
  std::vector<Magnets::Magnet> magnetInfo{
     {Magnets::DipoleY,    "D1_neg", kBlue,    std::make_pair(0.,-3.529),
      Geometry::Dimensions{-10., 10., -10., 10., -5840.0-945.0, -5840.0}
     },
     
     {Magnets::Quadrupole, "Q4_neg", kYellow,  std::make_pair(200.34,-200.34),
      Geometry::Dimensions{-10., 10., -10., 10., -4730.0-630.0, -4730.0}
     },
     
     {Magnets::Quadrupole, "Q3_neg", kYellow,  std::make_pair(-200.34,200.34),
      Geometry::Dimensions{-10., 10., -10., 10., -3830.0-550.0, -3830.0}
     },
     
     {Magnets::Quadrupole, "Q2_neg", kYellow,  std::make_pair(-200.34,200.34),
      Geometry::Dimensions{-10., 10., -10., 10., -3180.0-550.0, -3180.0}
     },
     
     {Magnets::Quadrupole, "Q1_neg", kYellow,  std::make_pair(200.34,-200.34),
      Geometry::Dimensions{-10., 10., -10., 10., -2300.0-630.0, -2300.0}
     },
      
     {Magnets::DipoleX,    "D_corr", kBlue+1,  std::make_pair(-1.1716,0.),
      Geometry::Dimensions{-10., 10., -10., 10., -1920.0-190.0, -1920.0}
     },
      
     {Magnets::DipoleX,    "Muon"  , kMagenta, std::make_pair(0.67,0.),
      Geometry::Dimensions{-10., 10., -10., 10., -750.0-430.0,  -750.0}
     },
      
     {Magnets::Quadrupole, "Q1_pos", kYellow,  std::make_pair(200.34,-200.34),
      Geometry::Dimensions{-10., 10., -10., 10., 2300.0, 2300.0+630.0}
     },
      
     {Magnets::Quadrupole, "Q2_pos", kYellow,  std::make_pair(-200.34,200.34),
      Geometry::Dimensions{-10., 10., -10., 10., 3180.0, 3180.0+550.0}
     },
      
     {Magnets::Quadrupole, "Q3_pos", kYellow,  std::make_pair(-200.34,200.34),
      Geometry::Dimensions{-10., 10., -10., 10., 3830.0, 3830.0+550.0}
     },
      
     {Magnets::Quadrupole, "Q4_pos", kYellow,  std::make_pair(200.34,-200.34),
      Geometry::Dimensions{-10., 10., -10., 10., 4730.0, 4730.0+630.0}
     },
      
     {Magnets::DipoleY,    "D1_pos", kBlue,    std::make_pair(0.,-3.529),
      Geometry::Dimensions{-10., 10., -10., 10., 5840.0, 5840.0+945.0}
     }
  };
```

-----------------------

**Apparatus**

<figure>
<img src="/DirectFlow/images/Apparatus.jpg" alt="Apparatus" style="width:100%">
<figcaption>Apparatus</figcaption>
</figure> 

---------------------------------

You can find the code used for plotting [here][plotcode].

[plotcode]: https://github.com/b-fontana/DirectFlow/blob/master/python/trajectory.py
