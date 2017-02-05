#include "atom.hpp"
#include <fstream>
#include <iostream>
#include <fstream>
#include <vector>

std::vector<double> Atom::GetXyzPosition(double x, double y, double z){
    std::vector<double> pos;
    pos.push_back(x);
    pos.push_back(y);
    pos.push_back(z);
    return(pos);
}