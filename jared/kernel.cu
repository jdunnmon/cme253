/*
 *  Copyright 2016 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

//#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include "./debug.h"

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

                //double rad = GetRadius(code.at(0));
                double rad = 0.0;
                prot_radii.push_back(rad);

            }

            // some checks
            if(prot_atomnums.size() != prot_resnums.size()){
                std::cerr << "ERROR: Problem in protein file" << std::endl;
            }

            if(prot_atomnums.size() != prot_xyz_coords.size()){
                std::cerr << "ERROR: Problem in protein file" << std::endl;
            }
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

            //double rad = GetRadius(code.at(0));
            double rad = 0.0;
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
    for (unsigned int ii = 0; ii < lig_trajnums.size(); ii++){
        for (unsigned int jj =0; jj < prot_atomnums.size(); jj++){ 
            double temp_dist = ComputeSquaredDistance(lig_xyz_coords[ii],
                                                      prot_xyz_coords[jj]);
            temp_dist = sqrt(temp_dist)/10.;
            all_distances.push_back(temp_dist);
        }
    }
    return all_distances;
}


__global__ void cuContacts(double *pxyz, double *lxyz, double *cudists, int *plength, int *llength)
{
  int pidx = threadIdx.x + blockIdx.x * blockDim.x;
  int lidx = threadIdx.y + blockIdx.y * blockDim.y;
  //cudists[index] = 20.;
  
 /* if(pidx % 3 == 0 && lidx % 3 == 0 ){
     cudists[pidx/3 + plength[0]*(lidx/3)] = ( sqrt(
               (pxyz[pidx]-lxyz[lidx])*(pxyz[pidx]-lxyz[lidx])
             + (pxyz[pidx+1]-lxyz[lidx+1])*(pxyz[pidx+1]-lxyz[lidx+1])
             + (pxyz[pidx+2]-lxyz[lidx+2])*(pxyz[pidx+2]-lxyz[lidx+2])  )/10. ); 
  }*/

  if ( (pidx < plength[0]) && (lidx< llength[0])){
    cudists[pidx+plength[0]*lidx] = ( sqrt(
               (pxyz[pidx*3]-lxyz[lidx*3])*(pxyz[pidx*3]-lxyz[lidx*3])
             + (pxyz[pidx*3+1]-lxyz[lidx*3+1])*(pxyz[pidx*3+1]-lxyz[lidx*3+1])
             + (pxyz[pidx*3+2]-lxyz[lidx*3+2])*(pxyz[pidx*3+2]-lxyz[lidx*3+2])  )/10. );

  }  

  __syncthreads();

}

//#define N (2048*2048)
//#define THREADS_PER_BLOCK 512

#define THREADS_PER_BLOCK_X 32
#define THREADS_PER_BLOCK_Y 32

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

  /* compute distanes using cpp*/
  std::vector<double> distances = LPContactFeaturizer(prot_atomnums,
                                                      prot_xyz_coords,
                                                      lig_trajnums,
                                                      lig_xyz_coords);

  /* print out cpp distances to a file */
  std::cout << "Number of cpp distances computed : " << distances.size() << std::endl;
  std::cout << "First cpp distance computed : " << distances[0] << std::endl;
  if(prot_atomnums.size()*lig_trajnums.size() > 50){
    std::cout << "50th cpp distance computed : " << distances[50] << std::endl; }
  if(prot_atomnums.size()*lig_trajnums.size() > 100000){
    std::cout << "100kth cpp distance computed : " << distances[100000] << std::endl; }
  /*
  std::ofstream f("cpp_distances.txt");
  if(f.is_open()){
    for(unsigned int k = 0; k < distances.size(); k++){
                        f << distances[k] << std::endl;
                    }
  }
  f.close();
  */
  
  double *pxyz, *lxyz, *cudists;
  double *d_pxyz, *d_lxyz, *d_cudists;
  int *plength, *d_plength;
  int *llength, *d_llength;
  int protein_size = prot_atomnums.size()*3;
  int ligand_traj_size = lig_trajnums.size()*3;
  int cudists_size = protein_size/3 * ligand_traj_size/3;

/* get GPU device number and name */
  int dev;
  cudaDeviceProp deviceProp;
  checkCUDA( cudaGetDevice( &dev ) );
  checkCUDA( cudaGetDeviceProperties( &deviceProp, dev ) );
  printf("Using GPU %d: %s\n", dev, deviceProp.name );

/* allocate space for device copies of a, b, c */
  checkCUDA( cudaMalloc( (void **) &d_pxyz, protein_size*sizeof(double)) );
  checkCUDA( cudaMalloc( (void **) &d_lxyz, ligand_traj_size*sizeof(double)) );
  checkCUDA( cudaMalloc( (void **) &d_cudists, cudists_size*sizeof(double) ));
  checkCUDA( cudaMalloc( (void **) &d_plength, sizeof(int) ));
  checkCUDA( cudaMalloc( (void **) &d_llength, sizeof(int) ));

/* allocate space for host copies of a, b, c and setup input values */
  pxyz = (double *)malloc( protein_size *sizeof(double));
  lxyz = (double *)malloc( ligand_traj_size *sizeof(double));
  cudists = (double *)malloc( cudists_size *sizeof(double));
  plength = (int *)malloc( sizeof(int));
  llength = (int *)malloc( sizeof(int));

  for(unsigned int pp = 0; pp < prot_atomnums.size(); pp++){
    pxyz[pp*3] = prot_xyz_coords[pp][0];
    pxyz[pp*3+1] = prot_xyz_coords[pp][1];
    pxyz[pp*3+2] = prot_xyz_coords[pp][2];
    //std::cout << "Last index of pxyz was " << pp*3+2 << std::endl;
  }

  for(unsigned int ll = 0; ll < lig_trajnums.size(); ll++){
    lxyz[ll*3] = lig_xyz_coords[ll][0];
    lxyz[ll*3+1] = lig_xyz_coords[ll][1];
    lxyz[ll*3+2] = lig_xyz_coords[ll][2];
    //std::cout << "Last index of lxyz was " << ll*3+2 << std::endl;
  }
  // for(unsigned int xx = 0; xx < lig_trajnums.size()*3; xx++){
  //   std::cout << lxyz[xx] << std::endl;
  // }

  plength[0] = prot_atomnums.size();
  llength[0] = lig_trajnums.size();

/* copy inputs to device */
  checkCUDA( cudaMemcpy( d_pxyz, pxyz, protein_size*sizeof(double), cudaMemcpyHostToDevice ) );
  checkCUDA( cudaMemcpy( d_lxyz, lxyz, ligand_traj_size*sizeof(double), cudaMemcpyHostToDevice ) );
  checkCUDA( cudaMemcpy( d_plength, plength, sizeof(int), cudaMemcpyHostToDevice) );  
  checkCUDA( cudaMemcpy( d_llength, llength, sizeof(int), cudaMemcpyHostToDevice) );


/* zero out the C array */
  checkCUDA( cudaMemset( d_cudists, 0, cudists_size*sizeof(double) ) );

/* setup threadblock size and grid sizes*/

  dim3 threads(plength[0], llength[0], 1);
  dim3 blocks(cudists_size/plength[0]+1,
              cudists_size/llength[0]+1, 1);

/* launch the kernel on the GPU */
  cuContacts<<< blocks, threads >>>( d_pxyz, d_lxyz, d_cudists, d_plength, d_llength );
  checkKERNEL();

/* copy result back to host */
  checkCUDA( cudaMemcpy( cudists, d_cudists, cudists_size*sizeof(double), cudaMemcpyDeviceToHost ) );

  //std::cout << "Number of cuda distances computed : " << sizeof cudists << std::endl; 
  std::cout << "First cuda distance computed : " << cudists[0] << std::endl;
  if(prot_atomnums.size()*lig_trajnums.size() > 50){
    std::cout << "50th cuda distance computed : " << cudists[50] << std::endl; }
  if(prot_atomnums.size()*lig_trajnums.size() > 100000){
   std::cout << "100kth cuda distance computed : " << cudists[100000] << std::endl; }


  std::ofstream f("distances.txt");
  if(f.is_open()){
    for(unsigned int k = 0; k < distances.size(); k++){
                        f << distances[k] << "  " << cudists[k] << std::endl;
                    }
  }
  f.close();


  free(pxyz);
  free(lxyz);
  free(cudists);
  free(plength);
  checkCUDA( cudaFree( d_pxyz ) );
  checkCUDA( cudaFree( d_lxyz ) );
  checkCUDA( cudaFree( d_cudists ) );
  checkCUDA( cudaFree( d_plength ) );

  checkCUDA( cudaDeviceReset () ); 
    
  return 0;
} /* end main */
