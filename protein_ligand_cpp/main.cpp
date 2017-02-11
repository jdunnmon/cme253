#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <string>
#include <vector>

double GetRadius(char elem){
    if (elem == 'N'){
        return 0.56;
    }
    else if (elem == 'H'){
        return 0.53;
    }
    else if (elem == 'C'){
        return 0.67;
    }
    else if (elem == 'O'){
        return 0.48;
    }
    else if (elem == 'S'){
        return 0.88;
    }
    else {
        return 0.;
    }

}

void ProteinSetup(std::string protein_inputfile,
                  std::vector<int>& prot_atomnums,
                  std::vector<int>& prot_resnums,
                  std::vector<std::vector<double>>& prot_xyz_coords,
                  std::vector<double>& prot_radii){
    std::ifstream f(protein_inputfile.c_str());
      if (f.is_open()) {
            std::string klass, code, resname, chain;
            int atomnum, resnum;
            double x, y, z, occ, temp;

            while (f >> klass >> atomnum >> code >> resname
                     >> chain >> resnum >> x >> y >> z
                     >> occ >> temp){

                std::vector<double> temp_coord;
                temp_coord.push_back(x);
                temp_coord.push_back(y);
                temp_coord.push_back(z); 

                prot_atomnums.push_back(atomnum);
                prot_resnums.push_back(resnum);
                prot_xyz_coords.push_back(temp_coord);

                double rad = GetRadius(code.at(0));
                prot_radii.push_back(rad);

            }

            // some checks
            if(prot_atomnums.size() != prot_resnums.size()){
                std::cerr << "ERROR: Problem in protein file" << std::endl;
            }

            if(prot_atomnums.size() != prot_xyz_coords.size()){
                std::cerr << "ERROR: Problem in protein file" << std::endl;
            }

            // for(unsigned int k = 0; k < prot_xyz_coords.size(); k++){
            //                 std::cout << prot_xyz_coords[k][0] << std::endl;
            // }

    }
    std::cout << "Lines in protein file : " << prot_atomnums.size() << std::endl;
}

void LigandTrajSetup(std::string ligand_inputfile,
                     std::vector<int>& lig_trajnums,
                     std::vector<int>& lig_atomnums,
                     std::vector<int>& lig_resnums,
                     std::vector<std::vector<double>>& lig_xyz_coords,
                     std::vector<double>& lig_radii){
    std::ifstream f(ligand_inputfile.c_str());
    if (f.is_open()) {
        std::string klass, code, resname, chain;
        int trajnum, atomnum, resnum;
        double x, y, z, occ, temp;

        while (f >> trajnum
                 >> klass >> atomnum >> code >> resname
                 >> chain >> resnum >> x >> y >> z
                 >> occ >> temp){

            std::vector<double> temp_coord;
            temp_coord.push_back(x);
            temp_coord.push_back(y);
            temp_coord.push_back(z); 

            lig_trajnums.push_back(trajnum);
            lig_atomnums.push_back(atomnum);
            lig_resnums.push_back(resnum);
            lig_xyz_coords.push_back(temp_coord);

            double rad = GetRadius(code.at(0));
            lig_radii.push_back(rad);

        }

        // some checks
        if(lig_atomnums.size() != lig_trajnums.size()){
            std::cerr << "ERROR: Problem in ligand file" << std::endl;
        }

        if(lig_atomnums.size() != lig_resnums.size()){
            std::cerr << "ERROR: Problem in ligand file" << std::endl;
        }

        if(lig_atomnums.size() != lig_xyz_coords.size()){
            std::cerr << "ERROR: Problem in ligand file" << std::endl;
        }
    }

    std::cout << "Lines in ligand file : " << lig_atomnums.size() << std::endl;
}

double ComputeSquaredDistance(std::vector<double> v1, std::vector<double> v2){
    double dist_squared;
    dist_squared = {  (v1[0]-v2[0])*(v1[0]-v2[0])
                    + (v1[1]-v2[1])*(v1[1]-v2[1])
                    + (v1[2]-v2[2])*(v1[2]-v2[2]) };
    return dist_squared;
}

std::vector<double> LPContactFeaturizer(std::vector<int>& prot_atomnums,
                                        std::vector<std::vector<double>>& prot_xyz_coords,
                                        std::vector<int>& lig_trajnums,
                                        std::vector<std::vector<double>>& lig_xyz_coords){

    std::vector<double> all_distances;

    // for every line in the ligand traj file
    for (unsigned int ii = 0; ii < lig_trajnums.size(); ii++){

        // for every atom in the protein i.e. line in the protein file
        for (unsigned int jj =0; jj < prot_atomnums.size(); jj++){ 

            double temp_dist = ComputeSquaredDistance(lig_xyz_coords[ii],
                                                      prot_xyz_coords[jj]);

            temp_dist = sqrt(temp_dist)/10.;
            all_distances.push_back(temp_dist);
        }
    }

    return all_distances;
}

std::vector<double> CudaLPContactFeaturizer(std::vector<int>& prot_atomnums,
                                        std::vector<std::vector<double>>& prot_xyz_coords,
                                        std::vector<int>& lig_trajnums,
                                        std::vector<std::vector<double>>& lig_xyz_coords){

    std::vector<double> all_distances;

    int n_frames = (lig_trajnums[lig_trajnums.size() - 1]) + 1;
    std::cout << "Number of frames in trajectory : " << n_frames << std::endl;


    // call a kernel on each frame (18 x 4500 distances)
    // optimize concurrency to accomplish this

    // transfer chunks of protein and entire trajectory
    // put protein in shared memory

    // shared memory optimization - optimize memory access patterns
    // so that we don't have bank conflicts that slow us down

    return all_distances;
}


int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    std::cout << "Usage:" << std::endl;
    {std::cout << "  " << argv[0] << " <protein input file> "
                       << " <ligand input file> " << std::endl;}
    return 0;
  }

  std::string protein_inputfile = argv[1];
  std::string ligand_inputfile = argv[2];

  std::vector<int> prot_atomnums;
  std::vector<int> prot_resnums;
  std::vector<std::vector<double>> prot_xyz_coords;
  std::vector<double> prot_radii;

  std::vector<int> lig_trajnums;
  std::vector<int> lig_atomnums;
  std::vector<int> lig_resnums;
  std::vector<std::vector<double>> lig_xyz_coords;
  std::vector<double> lig_radii;

  ProteinSetup(protein_inputfile,
               prot_atomnums,
               prot_resnums,
               prot_xyz_coords,
               prot_radii);

  LigandTrajSetup(ligand_inputfile,
                  lig_trajnums,
                  lig_atomnums,
                  lig_resnums,
                  lig_xyz_coords,
                  lig_radii);

  // num_protein_atoms = 
  // n_frames = 

  std::vector<double> distances = LPContactFeaturizer(prot_atomnums,
                                                      prot_xyz_coords,
                                                      lig_trajnums,
                                                      lig_xyz_coords);

  std::vector<double> cudaDistances = CudaLPContactFeaturizer(prot_atomnums,
                                                              prot_xyz_coords,
                                                              lig_trajnums,
                                                              lig_xyz_coords);

  // for(unsigned int k = 0; k < distances.size(); k++){
  //                       std::cout << distances[k] << std::endl;
  //                   }

  std::ofstream f("cpp_distances.txt");
  if(f.is_open()){
    for(unsigned int k = 0; k < distances.size(); k++){
                        f << distances[k] << std::endl;
                    }
  }
  f.close();

  return 0;

}

