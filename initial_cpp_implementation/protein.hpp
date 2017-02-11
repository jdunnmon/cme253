#ifndef PROTEIN_HPP
#define PROTEIN_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "atom.hpp"
#include "protein.hpp"

class Protein
{
    private:

    public:
        std::vector<Atom> residues;
        int Setup(std::string inputfile);

        double ComputeDistance(Atom a1, Atom a2);
        std::vector<double> ContactFeaturizer();

};

#endif /* PROTEIN_HPP */