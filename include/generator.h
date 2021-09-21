#ifndef GENERATOR_H
#define GENERATOR_H

#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TRandom.h"
#include <random>
#include <fstream>

template <class T>
class Generator {
public:
  Generator() {
    mRng = std::mt19937( mSeeder() );
    gRandom->SetSeed(0);
  }

  virtual T generate() = 0;

  void test(std::string filename, unsigned niters=100000) {
    std::fstream filetest;
    filetest.open(filename, std::ios_base::out);
    filetest << "data" << std::findendl;
    for(unsigned i = 0; i<niters; ++i)
      filetest << generate() << std::endl;
    filetest.close();
  }

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

template <class T>
class FermiDistribution final : public Generator<T> {
public:

  FermiDistribution() : Generator<T>() {    

    mFermiProbGraph = new TGraph();
    mFermiProbGraph->Set(mNPoints);
    
    mDist = new TH1D("h_fermi_prob","h_fermi_prob",
		     static_cast<int>(mNPoints*2.5), 0, 0.65);
    
    for(int i(0); i<mNPoints; i++)
      mFermiProbGraph->SetPoint(i, mPtFermi[i], mProbFermi[i]);

    for(int i(1); i<=mDist->GetNbinsX(); i++)
      {
	double pt_val = mDist->GetBinCenter(i);
        double value  = mFermiProbGraph->Eval(pt_val);
        if(value < 0.0) value = 0.0;
        mDist->SetBinContent(i,value);
      }
  }

  T generate() { return mDist->GetRandom(); }

private:
  TGraph* mFermiProbGraph;
  TH1D* mDist;
  
  static constexpr int mNPoints = 71;
  static constexpr double mPtFermi[mNPoints] = {0.0206,0.0272,0.0322,0.0361,0.0419,0.0485,0.0551,0.0614,0.0664,0.0707,0.0773,0.0831,0.0901,0.0959,0.104,0.114,0.122,0.131,0.135,0.141,0.149,0.16,0.168,0.178,0.185,0.194,0.203,0.209,0.213,0.218,0.225,0.231,0.245,0.252,0.259,0.267,0.275,0.282,0.289,0.295,0.304,0.309,0.316,0.323,0.33,0.338,0.346,0.353,0.363,0.371,0.377,0.384,0.394,0.403,0.411,0.425,0.438,0.448,0.461,0.476,0.49,0.505,0.519,0.533,0.551,0.565,0.579,0.592,0.609,0.624,0.638};
  static constexpr double mProbFermi[mNPoints] = {0.195,0.356,0.461,0.572,0.761,0.963,1.18,1.4,1.58,1.73,1.94,2.14,2.33,2.51,2.67,2.83,2.92,2.99,3,3.01,2.98,2.93,2.88,2.83,2.79,2.72,2.67,2.62,2.6,2.58,2.53,2.5,2.42,2.37,2.32,2.28,2.2,2.14,2.07,2.01,1.95,1.88,1.82,1.75,1.68,1.59,1.54,1.46,1.4,1.33,1.28,1.24,1.19,1.14,1.1,1.05,1.02,0.991,0.97,0.928,0.921,0.901,0.88,0.852,0.796,0.768,0.733,0.698,0.663,0.628,0.579};
};

#endif // GENERATOR_H
