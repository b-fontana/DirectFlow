#include <iostream>
#include <vector>
#include <fstream>
#include <boost/program_options.hpp>

#include "TRandom.h"

#include "include/tqdm.h"
#include "include/generator.h"

struct InputArgs {
public:
  unsigned niterations;
};

template <typename T>
class Point {
public:
  Point(T x, T y): x_(x), y_(y) {}
  T x() { return x_; }
  T y() { return y_; }
  
private:
  T x_;
  T y_;
};
  
//calculates area of polygon with four non-intersecting vertices
template <typename T>
T polygon_area(Point<T> v1, Point<T> v2, Point<T> v3, Point<T> v4) {
  T sum = 0.;

  auto determinant = [](Point<T> p1, Point<T> p2) -> T {
		       return p1.x()*p2.y() - p2.x()*p1.y();
		     };
    
  sum += determinant(v1, v2);
  sum += determinant(v2, v3);
  sum += determinant(v3, v4);
  sum += determinant(v4, v1);

  return 0.5 * sum;
}

unsigned size_last_batch(unsigned nbatches, unsigned nelems, unsigned batchSize) {
  return nelems-(nbatches-1)*batchSize;
}


void print_startup_info(const InputArgs& args) {  
  std::cout << " --- Simulation Information --- " << std::endl;
  // std::cout << "Batch Size: " << batchSize << " (last batch: "
  // 	    << size_last_batch(nbatches, args.nparticles, batchSize)
  // 	    << ")" << std::endl;
  // std::cout << "Number of batches: " << nbatches << std::endl;
  std::cout << "--------------------------" << std::endl;
}

//Calculate the four intersecting vertices of four lines, each pair representing
//a beam with distance equal to their beam width.
//The left beam is fixed, the right beam position follows theta:
//y1_left  = 0
//y2_left  = beam_width
//y1_right = slope * x
//y2_right = slope * x + b2
//where the slope is obtained from theta.
//The four vertices are encoded in a clockwise order as v1, v2, v3 and v4,
//starting from the lower left, where v4 always coincides with the origin.
//'kbig' is reserved for the distance between v1 and x=0 along the x axis
//'ksmall' is reserved for the distance between v2 and x=0 along the x axis
//The function calculates and stores the theta-dependent area of the 4-vertex polygon.
double calculate_area_with_output(double theta, double beamWidth2) {

  const double slope = std::tan(theta);
  const double b2 = beamWidth2 / std::cos(theta);
  const double kbig = -beamWidth2 / std::sin(theta);
  const double ksmall = (beamWidth2-b2) / slope;
      
  Point<double> v1(kbig, 0.);
  Point<double> v2(ksmall, beamWidth2);
  Point<double> v3(beamWidth2/slope, beamWidth2);
  Point<double> v4(0, 0); //by definition

  //the area will be negative since the vertices are arranged in clockwise order
  //see https://mathworld.wolfram.com/PolygonArea.html
  return 0.01 * std::abs( polygon_area<double>(v1, v2, v3, v4) ); //convert from fm^2 to barn
 }

void run(const InputArgs& args)
{
  print_startup_info(args);
  
  //set output-related variables
  std::fstream file;
  std::string filename("data/area_of_intersection.csv");
  file.open(filename, std::ios_base::out);
  file << "Idx,Area0,Area1,Area2,Area3,Area4,Theta" << std::endl;
  
  //set global variables
  double leadR = 6.68; // [fm]
  std::vector<double> beamWidth2{leadR, 2*leadR, 3*leadR, 4*leadR, 5*leadR}; // [fm]

  UniformDistribution<double> thetadist(2.5e-4, 4e-4);
  //UniformDistribution<double> thetadist(M_PI/100, M_PI/2);

  //Monte Carlo
  for (unsigned iter : tq::trange(args.niterations))
    {
      double theta = thetadist.generate();
      
      double area0 = calculate_area_with_output(theta, beamWidth2[0]);
      double area1 = calculate_area_with_output(theta, beamWidth2[1]);
      double area2 = calculate_area_with_output(theta, beamWidth2[2]);
      double area3 = calculate_area_with_output(theta, beamWidth2[3]);
      double area4 = calculate_area_with_output(theta, beamWidth2[4]);

      char astr0[40], astr1[40], astr2[40], astr3[40], astr4[40];
      char theta_str[40];
      sprintf(astr0, "%.20f", area0);
      sprintf(astr1, "%.20f", area1);
      sprintf(astr2, "%.20f", area2);
      sprintf(astr3, "%.20f", area3);
      sprintf(astr4, "%.20f", area4);
      sprintf(theta_str, "%.20f", theta);
      
      file << std::to_string( iter ) << ","
	   << astr0 << ","
	   << astr1 << ","
	   << astr2 << ","
	   << astr3 << ","
	   << astr4 << ","
	   << theta_str
	   << std::endl;	
    }

  file.close();
}

// run example: ./area_of_interaction.exe --niterations 10000
int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "produce this help message")
    ("niterations", po::value<unsigned>()->required(), "number of loop iterations");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc), vm);
  po::notify(vm);

  if(vm.count("help") or argc<2) {
    std::cerr << desc << std::endl;
    std::exit(0);
  }
      
  std::cout << "--- Executable options ---" << std::endl;
  for (const auto& it : vm) {
    std::cout << it.first.c_str() << ": ";
    auto& value = it.second.value();
    if (auto v = boost::any_cast<double>(&value))
      std::cout << *v << std::endl;
    else if (auto v = boost::any_cast<bool>(&value)) {
      std::string str_ = *v==1 ? "true" : "false";
      std::cout << *v << std::endl;
    }
    else if (auto v = boost::any_cast<std::string>(&value))
      std::cout << *v << std::endl;
    else if (auto v = boost::any_cast<unsigned>(&value))
      std::cout << *v << std::endl;
    else
      std::cerr << "type missing" << std::endl;
  }
  
  //run simulation   
  InputArgs info;
  info.niterations = boost::any_cast<unsigned>(vm["niterations"].value());
  
  run(info);

  std::cout << std::endl;
  return 0;
}
