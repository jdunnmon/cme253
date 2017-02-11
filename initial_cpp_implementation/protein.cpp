#include "atom.hpp"
#include "protein.hpp"
#include <fstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

int Protein::Setup(std::string inputfile){
  std::ifstream f(inputfile.c_str());
  if (f.is_open()) {
        std::string klass, r_code, r_resname, chain;
        int r_atomnum, r_resnum;
        double r_x, r_y, r_z, occ, temp;
        std::string elem;

        while (f >> klass >> r_atomnum >> r_code >> r_resname
                 >> chain >> r_resnum >> r_x >> r_y >> r_z
                 >> occ >> temp >> elem){
            Atom temp_atom;
            temp_atom.atomnum = r_atomnum;
            temp_atom.code = r_code;
            temp_atom.resnum = r_resnum;
            temp_atom.x = r_x;
            temp_atom.y = r_y;
            temp_atom.z = r_z;
            residues.push_back(temp_atom); 
            std::cout << r_resnum << std::endl;
        }

        for(unsigned int k = 0; k < residues.size(); k++){
                        std::cout << residues[k].x << std::endl;
                    }
  }

    return 0;
}

double Protein::ComputeDistance(Atom a1, Atom a2){
    double sum, ans;
    sum = pow((a1.x-a2.x),2) + pow((a1.y-a2.y),2) + pow((a1.z-a2.z),2);
    ans = sqrt(sum)/10.;
    return ans;
}

std::vector<double> Protein::ContactFeaturizer(){
    std::vector<double> distances;
    for(unsigned int i = 0; i < residues.size() - 2; i++){
        for(unsigned int j = i+3; j < residues.size(); j++){
            distances.push_back(ComputeDistance(residues[i],residues[j]));
        }
    }
    // for(a1 = residues[0]; a1 < residues.size() - 3; a1++){
    //     for(a2 = residues[a1.resnum+2], a2 < residues.size(), a2++){


    return distances;

}




