#ifndef DUMMY
#define DUMMY

#include <iostream>

class MyClass {
public:
  MyClass() { std::cout << "constructor" << std::endl; }
  void print();
};

#endif //DUMMY
