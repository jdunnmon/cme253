#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

//cuda error checking macros
#ifdef DEBUG
#define CUDA_CALL(F)  if( (F) != cudaSuccess ) \
  {printf("Error %s at %s:%d\n", cudaGetErrorString(cudaGetLastError()), \
   __FILE__,__LINE__); exit(-1);} 
#define CUDA_CHECK()  if( (cudaPeekAtLastError()) != cudaSuccess ) \
  {printf("Error %s at %s:%d\n", cudaGetErrorString(cudaGetLastError()), \
   __FILE__,__LINE__-1); exit(-1);} 
#else
#define CUDA_CALL(F) (F)
#define CUDA_CHECK() 
#endif

#define INDX( row, col ld) ( ( (col)* (ld) )+ (row) )

/*
Atom class and methods
*/

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

std::vector<double> Atom::GetXyzPosition(double x, double y, double z){
    std::vector<double> pos;
    pos.push_back(x);
    pos.push_back(y);
    pos.push_back(z);
    return(pos);
}

/*
Protein class and methods
*/

class Protein
{
    private:

    public:
        std::vector<Atom> residues;
        
    //    std::int length;
        void Setup(std::string inputfile);

        double ComputeDistance(Atom a1, Atom a2);
        std::vector<double> ContactFeaturizer();

        //double ComputeDistanceCuda(Atom a1, Atom a2);
      //  std::vector<double> ContactFeaturizerCuda();

};

void Protein::Setup(std::string inputfile){
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

}

double Protein::ComputeDistance(Atom a1, Atom a2){
    double sum, ans;
    sum = pow((a1.x-a2.x),2) + pow((a1.y-a2.y),2) + pow((a1.z-a2.z),2);
    ans = sqrt(sum)/10.; // Divide by 10 to get Angstroms
    return ans;
}

std::vector<double> Protein::ContactFeaturizer(){
    std::vector<double> distances;
    
    for(unsigned int i = 0; i < residues.size() - 2; i++){
        for(unsigned int j = i+3; j < residues.size(); j++){
            distances.push_back(ComputeDistance(residues[i],residues[j]));
        }
    }

    return distances;

}

/*__global__ std::vector<double> Protein::ContactFeaturizerCuda(int * ){
    std::vector<double> distances;
    
    int i = threadIdx.x+blockIdx.x*blockDim.x;
    int j = threadIdx.y+blockIdx.y*blockDim.y;
      
    if j>(i+2) && j<residues.size()
    
}
    


    return distances;
}*/

/*
Execute program
*/

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

  prot.Setup(inputfile);
  distances = prot.ContactFeaturizer();
  for(unsigned int k = 0; k < distances.size(); k++){
                        std::cout << distances[k] << std::endl;
                    }


  int npoints = 5;
  int size = npoints*3*sizeof(double);
  int distances = npoints*npoints*sizeof(double);
  

  return 0;

}


