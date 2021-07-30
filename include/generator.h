#ifndef GENERATOR_H
#define GENERATOR_H

#include "TF1.h"
#include <random>
#include <fstream>

template <class T>
class Generator {
public:
  Generator() { mRng = std::mt19937( mSeeder() ); }

  virtual T generate() = 0;

protected:
  std::mt19937 mRng;
  
private:
  std::random_device mSeeder;
};

template <class T>
class UniformDistribution final : public Generator<T> {
public:

  UniformDistribution(T left, T right) : Generator<T>() {
    mDist = std::uniform_real_distribution<T>(left, right);
  }
  
  T generate() { return mDist(Generator<T>::mRng); }

private:
  std::uniform_real_distribution<T> mDist;
};

template <class T>
class NormalDistribution final : public Generator<T> {
public:

  NormalDistribution(T mean, T sigma) : Generator<T>() {
    mDist = std::normal_distribution<T>(mean, sigma);
  }
  
  T generate() { return mDist(Generator<T>::mRng); }

private:
  std::normal_distribution<T> mDist;
};

template <class T>
class BoltzmannDistribution final : public Generator<T> {
public:

  BoltzmannDistribution(T pB, T pTemp, T pN, T pM0)
    : Generator<T>(), mB(pB), mTemp(pTemp), mN(pN), mM0(pM0) {    
    mDist = new TF1("BoltzmannDistribution_generator",
		    boltzmann_signature, 0.0, 100, 4); //GeV
    mDist->SetNpx(5000);
    mDist->SetParameters(mB, mTemp, mN, mM0);
    // mDist->SetLineColor(kBlack);
    // mDist->SetLineStyle(1);
    // mDist->SetLineWidth(1);
    //mDist->Draw();
  }

  T generate() { return mDist->GetRandom(); }

  void test(std::string filename, unsigned niters=100000) {
    std::fstream filetest;
    filetest.open(filename, std::ios_base::out);
    filetest << "data" << std::endl;
    for(unsigned i = 0; i<niters; ++i)
      filetest << generate() << std::endl;
    filetest.close();
  }

private:
  TF1* mDist;
  T mB, mTemp, mN, mM0;

  static double boltzmann_signature(double* x_val, double* par) {
    // One over pT term is removed -> original pT distribution
    // Fit function for d2N/(2pi dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    double pT, y, B, Temp, n, m0;
    pT      = x_val[0];
    B       = par[0];
    Temp    = par[1];
    n       = par[2];
    m0      = par[3];
    double mT = TMath::Sqrt(pT*pT+m0*m0);
    y = pT*B/TMath::Power(1.0+(mT-m0)/(n*Temp),n);
    return y;
  }
};

#endif // GENERATOR_H
