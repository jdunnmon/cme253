#include "atom.hpp"
#include "protein.hpp"
#include <fstream>
#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage:" << std::endl;
    std::cout << "  " << argv[0] << " <input file> " << std::endl;
    return 0;
  }

  std::string inputfile   = argv[1];

  Protein prot;
  std::vector<double> distances;

  int status = prot.Setup(inputfile);
  // if (status)
  // {
  //   std::cerr << "ERROR: System setup was unsuccessful!" << std::endl;
  //   return 1;
  // }

  distances = prot.ContactFeaturizer();
  for(unsigned int k = 0; k < distances.size(); k++){
                        std::cout << distances[k] << std::endl;
                    }


  return 0;

}