#ifndef ATOM_HPP
#define ATOM_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "atom.hpp"

class Atom
{
    private:

    public:
        int atomnum;
        int resnum;
        std::string code;
        double x;
        double y;
        double z;
        //std::vector<double> pos;

        std::vector<double> GetXyzPosition(double x, double y, double z);
};

#endif /* ATOM_HPP */

