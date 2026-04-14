//###########################################################################
//# Copyright (c), The PANNAdevs group. All rights reserved.                #
//# This file is part of the PANNA code.                                    #
//#                                                                         #
//# The code is hosted on GitLab at https://gitlab.com/PANNAdevs/panna      #
//# For further information on the license, see the LICENSE.txt file        #
//###########################################################################

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_panna.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <omp.h>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */
// Helper funcion for fast integer power
double ipow(double b, int e){
    double r = 1.0;
    while (e>0){
        if (e & 1) r *= b;
        e >>= 1;
        b *= b;
    }
    return r;
}

// Helper function for matrix multiplication
void matmul(double* __restrict__ a, const double* b, const double* c, int iM, int jM, int kM){
  for(int i=0; i<iM; i++)
    for(int k=0; k<kM; k++)
      for(int j=0;j<jM; j++)
        a[i*jM+j] += b[i*kM+k]*c[k*jM+j];
}


// ########################################################
//                       Constructor
// ########################################################
//

PairPANNA::PairPANNA(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  single_enable = 0;
}

// ########################################################
// ########################################################


// ########################################################
//                       Destructor
// ########################################################
//

PairPANNA::~PairPANNA()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

// ########################################################
// ########################################################

// Function computing gvect and its derivative
void PairPANNA::compute_gvect(int ind1, double **x, int* type,
                              int* neighs, int num_neigh,
                              double *G, double* dGdx){
  const double epscorr = 0.001;
  double posx = x[ind1][0];
  double posy = x[ind1][1];
  double posz = x[ind1][2];
  // Elements to store neigh list for angular part
  // We allocate max possible size, so we don't need to reallocate
  int nan = 0;
  int ang_neigh[num_neigh], ang_type[num_neigh];
  double dists[num_neigh], diffx[num_neigh], diffy[num_neigh], diffz[num_neigh];
  double iRij[num_neigh], fcij[num_neigh], sindij[num_neigh];
  // Other variables to be reused in loops
  double dx, dy, dz, Rij, coeff, cent, gauss, fc, dgp, derx, dery, derz;
  double cos_ijk, fact0, fact1, fact2, fact3, Gterm, dg0, dg1, dg2;
  double dgdxj, dgdyj, dgdzj, dgdxk, dgdyk, dgdzk;
  double Rcent[par.RsN_ang], Rexp[par.RsN_ang];
  double sin_ijk[par.ThetasN], fcrad[par.ThetasN], mod_norm[par.ThetasN], rad_mod[par.ThetasN];
  int nind, indsh, indsh2, indsh3;

  // Loop on neighbours, compute radial part, store quantities for angular
  for(int n=0; n<num_neigh; n++){
    nind = neighs[n];
    dx = x[nind][0]-posx;
    dy = x[nind][1]-posy;
    dz = x[nind][2]-posz;
    Rij = sqrt(dx*dx+dy*dy+dz*dz);
    indsh = (type[nind]-1)*par.RsN_rad;
    //Exclude atoms that are not present in the default group all
    //This is mostly needed for GCMC
    int igroup = group->find("all");
    int groupbit = group->bitmask[igroup];

    if (Rij < par.Rc_rad and atom->mask[nind] & groupbit){

      for(int indr=0; indr<par.RsN_rad; indr++){
        cent = Rij - par.Rsi_rad[indr];
        gauss = exp( - par.eta_rad[indr] * cent * cent);
        fc = 0.5 * ( 1.0 + cos(Rij * par.iRc_rad) );
        dgp = ( par.iRc_rad_half * sin(Rij * par.iRc_rad) +
                par.twoeta_rad[indr] * fc * cent ) * gauss / Rij;
        G[indsh+indr] += gauss * fc;
        indsh2 = (indsh+indr)*(num_neigh+1)*3;
        derx = dgp*dx;
        dery = dgp*dy;
        derz = dgp*dz;
        dGdx[indsh2 + num_neigh*3     ] += derx;
        dGdx[indsh2 + num_neigh*3 + 1 ] += dery;
        dGdx[indsh2 + num_neigh*3 + 2 ] += derz;
        dGdx[indsh2 + n*3     ] -= derx;
        dGdx[indsh2 + n*3 + 1 ] -= dery;
        dGdx[indsh2 + n*3 + 2 ] -= derz;
      }
    }
    // If within radial cutoff, store quantities
    if (Rij < par.Rc_ang and atom->mask[nind] & groupbit){
      ang_neigh[nan] = n;
      ang_type[nan] = type[nind];
      dists[nan] = Rij;
      diffx[nan] = dx;
      diffy[nan] = dy;
      diffz[nan] = dz;
      iRij[nan] = 1.0/Rij;
      fcij[nan] = 0.5 * ( 1.0 + cos(Rij * par.iRc_ang) );
      sindij[nan] = sin(Rij * par.iRc_ang);
      nan++;
    }
  }

  // Loop on angular neighbours and fill angular part
  for(int n=0; n<nan-1; n++){
    for(int m=n+1; m<nan; m++){
      // Compute cosine
      cos_ijk = (diffx[n]*diffx[m] + diffy[n]*diffy[m] + diffz[n]*diffz[m])*(iRij[n]*iRij[m]);
      // Clamping for numerical stability
      cos_ijk = std::max(-0.999999999, std::min(cos_ijk, 0.999999999));
      // Gvect shift due to species
      indsh = par.typsh[ang_type[n]-1][ang_type[m]-1];
      // Loop over radial and angular separately to get independent terms
      for(int Rsi=0; Rsi<par.RsN_ang; Rsi++){
        Rcent[Rsi] = 0.5 * (dists[n] + dists[m]) - par.Rsi_ang[Rsi];
        Rexp[Rsi] = 2.0 * exp( - par.eta_ang[Rsi] * Rcent[Rsi] * Rcent[Rsi]);
      }
      for(int Thi=0; Thi<par.ThetasN; Thi++){
        sin_ijk[Thi] = sqrt(1.0 - cos_ijk * cos_ijk + 
                           epscorr * par.Thi_sin[Thi] * par.Thi_sin[Thi]);
        fcrad[Thi] = 0.5 * ( 1.0 + par.Thi_cos[Thi] * cos_ijk + par.Thi_sin[Thi] * sin_ijk[Thi] );
        mod_norm[Thi] = ipow( 0.5 * (1.0 + sqrt(1.0 + epscorr * par.Thi_sin[Thi] * par.Thi_sin[Thi] ) ), 
                        par.zeta[Thi]);
        rad_mod[Thi] = ipow(fcrad[Thi], par.zeta[Thi]-1) / mod_norm[Thi];
      }
      // Loop over all bins
      for(int Rsi=0; Rsi<par.RsN_ang; Rsi++){
        for(int Thi=0; Thi<par.ThetasN; Thi++){
          indsh2 = Rsi * par.ThetasN + Thi;
          fact0 = Rexp[Rsi] * rad_mod[Thi];
          fact1 = fact0 * fcij[n] * fcij[m];
          fact2 = par.zeta_half[Thi] * fact1 * ( par.Thi_cos[Thi] - par.Thi_sin[Thi] * 
                               cos_ijk / sin_ijk[Thi] );
          fact3 = par.iRc_ang_half * fact0 * fcrad[Thi];
          Gterm = fact1 * fcrad[Thi];
          // Filling G contribution
          G[indsh+indsh2] += Gterm;
          dg0 = -iRij[n] * ( par.eta_ang[Rsi] * Rcent[Rsi] * Gterm
                    + fact2 * cos_ijk * iRij[n]
                    + fact3 * fcij[m] * sindij[n] );
          dg1 = fact2 * iRij[n] * iRij[m];
          dg2 = -iRij[m] * ( par.eta_ang[Rsi] * Rcent[Rsi] * Gterm
                    + fact2 * cos_ijk * iRij[m]
                    + fact3 * fcij[n] * sindij[m] );
          // Computing the derivative contributions
          dgdxj = dg0*diffx[n] + dg1*diffx[m];
          dgdyj = dg0*diffy[n] + dg1*diffy[m];
          dgdzj = dg0*diffz[n] + dg1*diffz[m];
          dgdxk = dg1*diffx[n] + dg2*diffx[m];
          dgdyk = dg1*diffy[n] + dg2*diffy[m];
          dgdzk = dg1*diffz[n] + dg2*diffz[m];
          // Filling all the interested terms
          indsh3 = (indsh+indsh2)*(num_neigh+1)*3;
          dGdx[indsh3 + ang_neigh[n]*3     ] += dgdxj;
          dGdx[indsh3 + ang_neigh[n]*3 + 1 ] += dgdyj;
          dGdx[indsh3 + ang_neigh[n]*3 + 2 ] += dgdzj;
          dGdx[indsh3 + ang_neigh[m]*3     ] += dgdxk;
          dGdx[indsh3 + ang_neigh[m]*3 + 1 ] += dgdyk;
          dGdx[indsh3 + ang_neigh[m]*3 + 2 ] += dgdzk;
          dGdx[indsh3 + num_neigh*3     ] -= dgdxj + dgdxk;
          dGdx[indsh3 + num_neigh*3 + 1 ] -= dgdyj + dgdyk;
          dGdx[indsh3 + num_neigh*3 + 2 ] -= dgdzj + dgdzk;
        }
      }
    }
  }

}

double PairPANNA::compute_network(double *G, double *dEdG, int type){
  // *1 layer input
  // *2 layer output
  double *lay1, *lay2, *dlay1, *dlay2;
  dlay1 = new double[par.layers_size[type][0]*par.gsize];
  lay1 = G;

  std::fill_n(dlay1,par.layers_size[type][0]*par.gsize,0.0);
  // dG_i/dG_i = 1
  for(int i=0; i<par.gsize; i++) dlay1[i*par.gsize+i] = 1.0;
  // Loop over layers
  for(int l=0; l<par.Nlayers[type]; l++){
    int size1 = par.layers_size[type][l];
    int size2 = par.layers_size[type][l+1];
    lay2 = new double[size2];
    dlay2 = new double[size2*par.gsize];
    std::fill_n(dlay2,size2*par.gsize,0.0);
    // Matrix vector multiplication done by hand for now...
    // We compute W.x+b and W.(dx/dg)
    for(int i=0; i<size2; i++){
      // a_i = b_i
      lay2[i] = network[type][2*l+1][i];
      for(int j=0;j<size1; j++){
        // a_i += w_ij * x_j
        lay2[i] += network[type][2*l][i*size1+j]*lay1[j];
      }
    }
    // da_i/dg_k += w_ij * dx_j/dg_k
    matmul(dlay2,network[type][2*l],dlay1,size2,par.gsize,size1);

    // Apply appropriate activation
    // Gaussian
    if(par.layers_activation[type][l]==1){
      for(int i=0; i<size2; i++){
        double tmp = exp(-lay2[i]*lay2[i]);
        for(int k=0; k<par.gsize; k++)
          dlay2[i*par.gsize+k] *= -2.0*lay2[i]*tmp;
        lay2[i] = tmp;
      }
    }
    // ReLU
    else if(par.layers_activation[type][l]==3){
      for(int i=0; i<size2; i++){
        if(lay2[i]<0){
          lay2[i] = 0.0;
          for(int k=0; k<par.gsize; k++) dlay2[i*par.gsize+k] = 0.0;
        }
      }
    }
    // Tanh
    else if(par.layers_activation[type][l]==4){
      for(int i=0; i<size2; i++){
        double tmp = tanh(lay2[i]);
        for(int k=0; k<par.gsize; k++)
          dlay2[i*par.gsize+k] *= (1 - tmp * tmp);
        lay2[i] = tmp;
      }
    }
    // Otherwise it's linear and nothing needs to be done

    if(l!=0) delete[] lay1;
    delete[] dlay1;
    lay1 = lay2;
    dlay1 = dlay2;
  }
  // myfile.close();
  for(int i=0;i<par.gsize;i++) dEdG[i]=dlay1[i];
  double E = lay1[0];
  delete[] lay1;
  delete[] dlay1;
  return E;
}

// ########################################################
//                       COMPUTE
// ########################################################
// Determine the energy and forces for the current structure.

void PairPANNA::compute(int eflag, int vflag)
{
  // if (eflag || vflag) ev_setup(eflag,vflag);
  // else evflag = vflag_fdotr = 0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  // I'll assume the order is the same.. we'll need to create a mapping if not the case
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int igroup = group->find("all");
  int groupbit = group->bitmask[igroup];

  double time_feat = 0.0;
  double time_fwbw = 0.0;
  double time_featalloc = 0.0;
  double time_featzero = 0.0;
  double time_chain = 0.0;
  double time_sum = 0.0;
  double time_virial = 0.0;
  double start;

  // Array to store energy contribution from each atom
  double Etmp[inum];
  std::fill_n(Etmp,inum,0.0);
  // Getting the max possible number of neighbors
  int nnmax = 0;
  for(int a=0; a<inum; a++) nnmax = std::max(nnmax, numneigh[ilist[a]]);
  // Array to store force contributions from each atom
  double ftmp[inum*(nnmax+1)*3];
  std::fill_n(ftmp,inum*(nnmax+1)*3,0.0);

  // Looping on local atoms
  #pragma omp parallel for
  for(int a=0; a<inum; a++){
    int myind = ilist[a];
    if (!(atom->mask[myind] & groupbit)) continue;

    start = omp_get_wtime();
    int nn = numneigh[myind];
    // Allocate this gvect and dG/dx and zero them
    double G[par.gsize];
    double dEdG[par.gsize];
    // dGdx has (numn+1)*3 derivs per elem: neigh first, then the atom itself
    double dGdx[par.gsize*(nn+1)*3];
    time_featalloc += (omp_get_wtime() - start);
    start = omp_get_wtime();
    std::fill_n(G,par.gsize,0.0);
    std::fill_n(dGdx,par.gsize*(nn+1)*3,0.0);
    time_featzero += (omp_get_wtime() - start);
    start = omp_get_wtime();
    // Calculate Gvect and derivatives
    compute_gvect(myind, x, type, firstneigh[myind], nn, G, dGdx);
    time_feat += (omp_get_wtime() - start);
    start = omp_get_wtime();
  
    // Apply network
    Etmp[a] = compute_network(G,dEdG,type[myind]-1);
    time_fwbw += (omp_get_wtime() - start);
    start = omp_get_wtime();

    // Calculate forces
    int shift = (nn+1)*3;
    int tshift = a*(nnmax+1)*3;
    for(int n=0; n<nn; n++){
      for(int j=0; j<par.gsize; j++){
        ftmp[tshift + 3*n    ] -= dEdG[j] * dGdx[j*shift + 3*n    ];
        ftmp[tshift + 3*n + 1] -= dEdG[j] * dGdx[j*shift + 3*n + 1];
        ftmp[tshift + 3*n + 2] -= dEdG[j] * dGdx[j*shift + 3*n + 2];
      }
    }
    for(int j=0; j<par.gsize; j++){
      ftmp[tshift + 3*nn    ] -= dEdG[j] * dGdx[j*shift + 3*nn    ];
      ftmp[tshift + 3*nn + 1] -= dEdG[j] * dGdx[j*shift + 3*nn + 1];
      ftmp[tshift + 3*nn + 2] -= dEdG[j] * dGdx[j*shift + 3*nn + 2];
    }
  }

    time_chain += (omp_get_wtime() - start);
    start = omp_get_wtime();

  // Summing all contributions to global values
  for(int a=0; a<inum; a++){
    int myind = ilist[a];
    if (!(atom->mask[myind] & groupbit)) continue;

    int nn = numneigh[myind];
    int tshift = a*(nnmax+1)*3;
    for(int n=0; n<nn; n++){
      int nind = firstneigh[myind][n];
      f[nind][0] += ftmp[tshift + 3*n    ];
      f[nind][1] += ftmp[tshift + 3*n + 1];
      f[nind][2] += ftmp[tshift + 3*n + 2];
    }
    f[myind][0] += ftmp[tshift + 3*nn    ];
    f[myind][1] += ftmp[tshift + 3*nn + 1];
    f[myind][2] += ftmp[tshift + 3*nn + 2];

    if (eflag_global) eng_vdwl += Etmp[a];
    if (eflag_atom) eatom[myind] += Etmp[a];
  }
    time_sum += (omp_get_wtime() - start);
    start = omp_get_wtime();
  if (vflag_fdotr) {
    virial_fdotr_compute();
  }
    time_virial += (omp_get_wtime() - start);

  // if (comm->me == 0) { // Only print from the master MPI rank
    printf("Step Time - Feat: %f s | fwbw: %f s | featalloc: %f s | featzero: %f s | chain: %f s | sum: %f s | virial: %f s\n",
    time_feat, time_fwbw, time_featalloc, time_featzero, time_chain, time_sum, time_virial);
  // }

}

// ########################################################
// ########################################################

// Get a new line skipping comments or empty lines
// Set value=... if [...], return 1
// Fill key,value if 'key=value', return 2
// Set value=... if ..., return 3
// Return 0 if eof, <0 if error, >0 if okay
int PairPANNA::get_input_line(std::ifstream* file, std::string* key, std::string* value){
  std::string line;
  int parsed = 0; int vc = 1;
  while(!parsed){
    std::getline(*file,line);
    // Exit on EOF
    if(file->eof()) return 0;
    // Exit on bad read
    if(file->bad()) return -1;
    // Remove spaces
    line.erase (std::remove(line.begin(), line.end(), ' '), line.end());
    // Skip empty line
    if(line.length()==0) continue;
    // Skip comments
    if(line.at(0)=='#') continue;
    // Parse headers
    if(line.at(0)=='['){
      *value = line.substr(1,line.length()-2);
      return 1;
    }
    // Check if we have version information:
    if(line.at(0)=='!') { vc=0 ;}
    // Look for equal sign
    std::string eq = "=";
    size_t eqpos = line.find(eq);
    // Parse key-value pair
    if(eqpos != std::string::npos){
      *key = line.substr(0,eqpos);
      *value = line.substr(eqpos+1,line.length()-1);
      if (vc == 0) { vc = 1 ; return 3; }
      return 2;
    }
    std::cout << line << std::endl;
    parsed = 1;
  }
  return -1;
}

int PairPANNA::get_parameters(char* directory, char* filename)
{
  //const double panna_pi = 3.14159265358979323846;
  // Parsing the potential parameters
  std::ifstream params_file;
  std::ifstream weights_file;
  std::string key, value;
  std::string dir_string(directory);
  std::string param_string(filename);
  std::string file_string(dir_string+"/"+param_string);
  std::string wfile_string;

  // Initializing some parameters before reading:
  par.Nspecies = -1;
  // Flags to keep track of set parameters
  int Npars = 17;
  int parset[Npars];
  for(int i=0;i<Npars;i++) parset[i]=0;
  int *spset;
  std::string version = "v0"; int gversion=0;
  double tmp_eta_rad; double tmp_eta_ang ; int tmp_zeta ;
  //
  params_file.open(file_string.c_str());
  // section keeps track of input file sections
  // -1 in the beginning
  // 0 for gvect params
  // i for species i (1 based)
  int section = -1;
  // parseint checks the status of input parsing
  int parseint = get_input_line(&params_file,&key,&value);
  while(parseint>0){
    // Parse line
    if(parseint==1){
      // Gvect param section
      if(value=="GVECT_PARAMETERS"){ section = 0; }
      // For now other sections are just species networks
      else {
        // First time after params are read: do checks
        if(section==0){
          if (gversion==0){
            if(parset[5]==0){
            // Set steps if they were omitted
            par.Rsst_rad = (par.Rc_rad - par.Rs0_rad) / par.RsN_rad; parset[5]=1;}
            if(parset[10]==0){
            par.Rsst_ang = (par.Rc_ang - par.Rs0_ang) / par.RsN_ang; parset[10]=1;}}
          else if (gversion==1){parset[5]=1 ; parset[10]=1;}
          // Check that all parameters have been set
          for(int p=0;p<Npars;p++){
            if(parset[p]==0){
              std::cout << "Parameter " << p << " not set!" << std::endl;  return -1; } }
          // Calculate Gsize
          par.gsize = par.Nspecies * par.RsN_rad + (par.Nspecies*(par.Nspecies+1))/2 * par.RsN_ang * par.ThetasN;
        } //section 0 ended
        int match = 0;
        for(int s=0;s<par.Nspecies;s++){
          // If species matches the list, change section
          if(value==par.species[s]){
            section = s+1;
            match = 1;
          }
        }
        if(match==0){
          std::cout << "Species " << value << " not found in species list." << std::endl;
          return -2;
        }
      }
    }// A header is parsed
    else if(parseint==2){
      // Parse param section
      if(section==0){
	std::string comma = ",";
        if(key=="Nspecies"){
          par.Nspecies = std::atoi(value.c_str());
          // Small check
          if(par.Nspecies<1){
            std::cout << "Nspecies needs to be >0." << std::endl;
            return -2; }
          parset[0] = 1;
          // Allocate species list
          par.species = new std::string[par.Nspecies];
          // Allocate network quantities
          par.Nlayers = new int[par.Nspecies];
          par.layers_size = new int*[par.Nspecies];
          par.layers_activation = new int*[par.Nspecies];
          network = new double**[par.Nspecies];
          // Keep track of set species
          spset = new int[par.Nspecies];
          for(int s=0;s<par.Nspecies;s++) {
            par.Nlayers[s] = -1;
            spset[s]=0; } }
        else if(key=="species"){
          //std::string comma = ",";
          size_t pos = 0;
          int s = 0;
          // Parse species list
          while ((pos = value.find(comma)) != std::string::npos) {
            if(s>par.Nspecies-2){
              std::cout << "Species list longer than Nspecies." << std::endl;
              return -2; }
            par.species[s] = value.substr(0, pos);
            value.erase(0, pos+1);  s++; }
          if(value.length()>0){
            par.species[s] = value; s++; };
          if(s<par.Nspecies){
            std::cout << "Species list shorter than Nspecies." << std::endl;
            return -2; }
          parset[1] = 1;
          std::cout << "Done reading." << std::endl;
          std::cout << par.species[1] << std::endl;
           }
        // Common features are read.
        // From here on what will be read depends on the gversion
        if(gversion == 0){ // Potentials compatible with OPENKIM
	  std::cout << "G Version is " << gversion << std::endl;
          if(key=="eta_rad"){
            tmp_eta_rad = std::atof(value.c_str()); parset[2] = 0; }
          else if(key=="Rc_rad"){
            par.Rc_rad = std::atof(value.c_str()); parset[3] = 1; }
          else if(key=="Rs0_rad"){
            par.Rs0_rad = std::atof(value.c_str());  parset[4] = 1; }
          else if(key=="Rsst_rad"){
            par.Rsst_rad = std::atof(value.c_str()); parset[5] = 1; }
          else if(key=="RsN_rad"){
            par.RsN_rad = std::atoi(value.c_str()); parset[6] = 1;
            par.eta_rad = new double[par.RsN_rad];
	    par.twoeta_rad = new double[par.RsN_rad];
            par.Rs_rad  = new double[par.RsN_rad];
            for(int i=0;i<par.RsN_rad;i++) par.eta_rad[i]=tmp_eta_rad;
            for(int i=0;i<par.RsN_rad;i++) par.Rs_rad[i]= par.Rs0_rad + i *(par.Rc_rad - par.Rs0_rad) / par.RsN_rad ;
            parset[14]=1; parset[2]=1;}
          else if(key=="eta_ang"){
            tmp_eta_ang = std::atof(value.c_str()); parset[7] = 0; }
          else if(key=="Rc_ang"){
            par.Rc_ang = std::atof(value.c_str()); parset[8] = 1; }
          else if(key=="Rs0_ang"){
            par.Rs0_ang = std::atof(value.c_str()); parset[9] = 1; }
          else if(key=="Rsst_ang"){
            par.Rsst_ang = std::atof(value.c_str()); parset[10] = 1; }
          else if(key=="RsN_ang"){
            par.RsN_ang = std::atoi(value.c_str()); parset[11] = 1;
            par.eta_ang = new double[par.RsN_ang];
            par.Rs_ang  = new double[par.RsN_ang];
            for(int i=0;i<par.RsN_ang;i++) par.eta_ang[i]=tmp_eta_ang;
            for(int i=0;i<par.RsN_ang;i++) par.Rs_ang[i]= par.Rs0_ang + i *(par.Rc_ang - par.Rs0_ang) / par.RsN_ang ;
            parset[15]=1; parset[7]=1;}
          else if(key=="zeta"){
            tmp_zeta = std::atof(value.c_str()); parset[12] = 0; }
          else if(key=="ThetasN"){
            par.ThetasN = std::atoi(value.c_str()); parset[13] = 1;
            par.zeta = new int[par.ThetasN];
	    par.zeta_half = new double[par.ThetasN];
            par.Thetas = new double[par.ThetasN];
            for(int i=0;i<par.ThetasN;i++) par.zeta[i]=tmp_zeta; parset[12]=1;
            for(int i=0;i<par.ThetasN;i++) par.Thetas[i]= (0.5f+ i)*(M_PI/par.ThetasN); parset[16]=1;}
        }//gversion = 0
        else if(gversion ==1 ){
          //First read allocation sizes
          if(key=="RsN_rad"){
            par.RsN_rad = std::atoi(value.c_str());
            par.eta_rad = new double[par.RsN_rad];
	    par.twoeta_rad = new double[par.RsN_rad];
            par.Rs_rad = new double[par.RsN_rad]; parset[6]=1;}
            else if(key=="RsN_ang"){
            par.RsN_ang = std::atoi(value.c_str());
            par.eta_ang = new double[par.RsN_ang];
            par.Rs_ang = new double[par.RsN_ang];
            
            // MOVED FROM ThetasN: Allocate to RsN_ang instead
            par.zeta = new int[par.RsN_ang];
            par.zeta_half = new double[par.RsN_ang];
            par.Thetas = new double[par.RsN_ang]; 
            
            parset[11]=1;
          }
          else if(key=="ThetasN"){
            par.ThetasN = std::atoi(value.c_str());
            // REMOVED allocations from here
            parset[13]=1;
          }
          // Then cutoffs
          else if(key=="Rc_rad"){
            par.Rc_rad = std::atof(value.c_str()); parset[3] = 1; }
          else if(key=="Rc_ang"){
            par.Rc_ang = std::atof(value.c_str()); parset[8] = 1; }
          // Then param arrays
          else if(key=="eta_rad"){
            size_t pos = 0; int s = 0;
            value=value.substr(1, value.size() - 2); //get rid of [ ]
            while ((pos = value.find(comma)) != std::string::npos) {
              par.eta_rad[s] = std::atof(value.substr(0, pos).c_str());
              value.erase(0, pos+1);  s++; }
            if(value.length()>0){ par.eta_rad[s] = std::atof(value.c_str()); s++; }; parset[2] = 1; }
          else if(key=="eta_ang") {
            size_t pos = 0; int s = 0;
            value=value.substr(1, value.size() - 2); //get rid of [ ]
            while ((pos = value.find(comma)) != std::string::npos) {
              par.eta_ang[s] = std::atof(value.substr(0, pos).c_str());
              value.erase(0, pos+1);  s++; }
            if(value.length()>0){ par.eta_ang[s] = std::atof(value.c_str()); s++; }; parset[7] = 1; }
          else if(key=="zeta") {
            size_t pos = 0; int s = 0;
            value=value.substr(1, value.size() - 2); //get rid of [ ]
            while ((pos = value.find(comma)) != std::string::npos) {
              par.zeta[s] = std::atoi(value.substr(0, pos).c_str());
              value.erase(0, pos+1);  s++; }
            if(value.length()>0){ par.zeta[s] = std::atoi(value.c_str()); s++; }; parset[12] = 1; }
          // Then the bin center arrays
          else if(key=="Rs_rad"){
            size_t pos = 0; int s = 0;
            value=value.substr(1, value.size() - 2); //get rid of [ ]
            while ((pos = value.find(comma)) != std::string::npos) {
              par.Rs_rad[s] = std::atof(value.substr(0, pos).c_str());
              value.erase(0, pos+1);  s++; }
            if(value.length()>0){ par.Rs_rad[s] = std::atof(value.c_str()); s++; }; parset[14] = 1;
            par.Rs0_rad=par.Rs_rad[0];                           parset[4]=1;}
          else if(key=="Rs_ang") {
            size_t pos = 0; int s = 0;
            value=value.substr(1, value.size() - 2); //get rid of [ ]
            while ((pos = value.find(comma)) != std::string::npos) {
              par.Rs_ang[s] = std::atof(value.substr(0, pos).c_str());
              value.erase(0, pos+1);  s++; }
            if(value.length()>0){ par.Rs_ang[s] = std::atof(value.c_str()); s++; }; parset[15] = 1;
            par.Rs0_ang=par.Rs_ang[0];                           parset[9]=1;}
          else if(key=="Thetas") {
            size_t pos = 0; int s = 0;
            value=value.substr(1, value.size() - 2); //get rid of [ ]
            while ((pos = value.find(comma)) != std::string::npos) {
              par.Thetas[s] = std::atof(value.substr(0, pos).c_str());
              value.erase(0, pos+1);  s++; }
            if(value.length()>0){ par.Thetas[s] = std::atof(value.c_str()); s++; };  parset[16] = 1; }
        } //gversion = 1
      } //Section 0 (Parameter parsing) is finished.
      // Parse species network
      else if(section<par.Nspecies+1){
        int s=section-1;
        // Read species network
        if(key=="Nlayers"){
          par.Nlayers[s] = std::atoi(value.c_str());
          // This has the extra gvect size
          par.layers_size[s] = new int[par.Nlayers[s]+1];
          par.layers_size[s][0] = par.gsize;
          par.layers_size[s][1] = 0;
          par.layers_activation[s] = new int[par.Nlayers[s]];
          for(int i=0;i<par.Nlayers[s]-1;i++) par.layers_activation[s][i]=1;
          par.layers_activation[s][par.Nlayers[s]-1]=0;
          network[s] = new double*[2*par.Nlayers[s]];
        }
        else if(key=="sizes"){
          if(par.Nlayers[s]==-1){
            std::cout << "Sizes cannot be set before Nlayers." << std::endl;
            return -3;
          }
          std::string comma = ",";
          size_t pos = 0;
          int l = 0;
          // Parse layers list
          while ((pos = value.find(comma)) != std::string::npos) {
            if(l>par.Nlayers[s]-2){
              std::cout << "Layers list longer than Nlayers." << std::endl;
              return -3;
            }
            std::string lsize = value.substr(0, pos);
            par.layers_size[s][l+1] = std::atoi(lsize.c_str());
            value.erase(0, pos+1);
            l++;
          }
          if(value.length()>0){
            par.layers_size[s][l+1] = std::atoi(value.c_str());
            l++;
          };
          if(l<par.Nlayers[s]){
            std::cout << "Layers list shorter than Nlayers." << std::endl;
            return -3;
          }
        }
        else if(key=="activations"){
          if(par.Nlayers[s]==-1){
            std::cout << "Activations cannot be set before Nlayers." << std::endl;
            return -3;
          }
          std::string comma = ",";
          size_t pos = 0;
          int l = 0;
          // Parse layers list
          while ((pos = value.find(comma)) != std::string::npos) {
            if(l>par.Nlayers[s]-2){
              std::cout << "Activations list longer than Nlayers." << std::endl;
              return -3;
            }
            std::string lact = value.substr(0, pos);
            int actnum = std::atoi(lact.c_str());
            if (actnum!=0 && actnum!=1 && actnum!=3 && actnum!=4 ){
              std::cout << "Activations unsupported: " << actnum << std::endl;
              return -3;
            }
            par.layers_activation[s][l] = actnum;
            value.erase(0, pos+1);
            l++;
          }
          if(value.length()>0){
            int actnum = std::atoi(value.c_str());
            if (actnum!=0 && actnum!=1 && actnum!=3 && actnum!=4){
              std::cout << "Activations unsupported: " << actnum << std::endl;
              return -3;
            }
            par.layers_activation[s][l] = actnum;
            l++;
          };
          if(l<par.Nlayers[s]){
            std::cout << "Activations list shorter than Nlayers." << std::endl;
            return -3;
          }
        }
        else if(key=="file"){
          if(par.layers_size[s][1]==0){
            std::cout << "Layers sizes unset before filename for species " << par.species[s] << std::endl;
            return -3;
          }
          // Read filename and load weights
          wfile_string = dir_string+"/"+value;
          weights_file.open(wfile_string.c_str(), std::ios::binary);
          if(!weights_file.is_open()){
            std::cout << "Error reading weights file for " << par.species[s] << std::endl;
            return -3;
          }
          for(int l=0; l<par.Nlayers[s]; l++){
            // Allocate and read the right amount of data
            // Weights
            network[s][2*l] = new double[par.layers_size[s][l]*par.layers_size[s][l+1]];
            for(int i=0; i<par.layers_size[s][l]; i++) {
              for(int j=0; j<par.layers_size[s][l+1]; j++) {
                float num;
                weights_file.read(reinterpret_cast<char*>(&num), sizeof(float));
                if(weights_file.eof()){
                  std::cout << "Weights file " << wfile_string << " is too small." << std::endl;
                  return -3;
                }
                network[s][2*l][j*par.layers_size[s][l]+i] = (double)num;
              }
            }
            // Biases
            network[s][2*l+1] = new double[par.layers_size[s][l+1]];
            for(int d=0; d<par.layers_size[s][l+1]; d++) {
              float num;
              weights_file.read(reinterpret_cast<char*>(&num), sizeof(float));
              if(weights_file.eof()){
                std::cout << "Weights file " << wfile_string << " is too small." << std::endl;
                return -3;
              }
              network[s][2*l+1][d] = (double)num;
            }
          }
          // Check if we're not at the end
          std::ifstream::pos_type fpos = weights_file.tellg();
          weights_file.seekg(0, std::ios::end);
          std::ifstream::pos_type epos = weights_file.tellg();
          if(fpos!=epos){
            std::cout << "Weights file " << wfile_string << " is too big." << std::endl;
            return -3;
          }
          weights_file.close();
          spset[section-1] = 1;
        }
      }
      else{
        return -3;
      }
    }
    else if(parseint==3){
      // Version information is read:
      if(key == "!version") {
        version = value ;
        std::cout << "Network version " << value << std::endl; }
      else if(key == "!gversion") {
        gversion = std::atoi(value.c_str()) ;
        std::cout << "Gvector version " << value << std::endl;}
    }
    // Get new line
    parseint = get_input_line(&params_file,&key,&value);
  }

  // Derived params - for both gvect types done here
  par.cutmax = par.Rc_rad>par.Rc_ang ? par.Rc_rad : par.Rc_ang;
  for(int i=0; i<par.RsN_rad; i++) {
    par.twoeta_rad[i] = 2.0*par.eta_rad[i];}
  for(int i=0; i<par.ThetasN; i++) {
    par.zeta_half[i] = 0.5f*par.zeta[i];}
  par.iRc_rad = M_PI/par.Rc_rad;
  par.iRc_rad_half = 0.5*par.iRc_rad;
  par.iRc_ang = M_PI/par.Rc_ang;
  par.iRc_ang_half = 0.5*par.iRc_ang;
  //par.Rsi_rad = new float[par.RsN_rad];
  //for(int indr=0; indr<par.RsN_rad; indr++) par.Rsi_rad[indr] = par.Rs0_rad + indr * par.Rsst_rad;
  par.Rsi_rad = par.Rs_rad;
  //par.Rsi_ang = new float[par.RsN_ang];
  //for(int indr=0; indr<par.RsN_ang; indr++) par.Rsi_ang[indr] = par.Rs0_ang + indr * par.Rsst_ang;
  par.Rsi_ang = par.Rs_ang;
  par.Thi_cos = new double[par.ThetasN];
  par.Thi_sin = new double[par.ThetasN];
  for(int indr=0; indr<par.ThetasN; indr++)  {
    //double ti = (indr + 0.5f) * M_PI / par.ThetasN;
    double ti = par.Thetas[indr];
    par.Thi_cos[indr] = cos(ti);
    par.Thi_sin[indr] = sin(ti);
  }
  for(int s=0;s<par.Nspecies;s++){
    if(spset[s]!=1){
      std::cout << "Species network undefined for " << par.species[s] << std::endl;
      return -4;
    }
  }

  // Precalculate gvect shifts for any species pair
  par.typsh = new int*[par.Nspecies];
  for(int s=0; s<par.Nspecies; s++){
    par.typsh[s] = new int[par.Nspecies];
    for(int ss=0; ss<par.Nspecies; ss++){
      if(s<ss) par.typsh[s][ss] = par.Nspecies*par.RsN_rad +
                  (s*par.Nspecies - (s*(s+1))/2 + ss) *
                  par.RsN_ang * par.ThetasN;
      else par.typsh[s][ss] = par.Nspecies*par.RsN_rad +
                  (ss*par.Nspecies - (ss*(ss+1))/2 + s) *
                  par.RsN_ang * par.ThetasN;
    }
  }
  params_file.close();
  delete[] spset;
  return(0);
}

// ########################################################
//                       ALLOCATE
// ########################################################
// Allocates all necessary arrays.

void PairPANNA::allocate()
{

  allocated = 1;
  int n = atom->ntypes;
  
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 1;
    }
  }
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

// ########################################################
// ########################################################

// ########################################################
//                       COEFF
// ########################################################
// Load all the gvectors and NN parameters

void PairPANNA::coeff(int narg, char **arg)
{

  if (!allocated) {
    allocate();
  }

  // We now expect a directory and the parameters file name (inside the directory) with all params
  if (narg != 4) {
    error->all(FLERR,"Format of pair_coeff command is\npair_coeff * *  network_directory parameter_file\n");
  }

  std::cout << "Loading PANNA pair parameters from " << arg[2] << "/" << arg[3] << std::endl;
  int gpout = get_parameters(arg[2], arg[3]);
  if(gpout==0){
    std::cout << "Network loaded!" << std::endl;
  }
  else{
    std::cout << "Error " << gpout << " while loading network!" << std::endl;
    exit(1);
  }

  for (int i=1; i<=atom->ntypes; i++) {
    for (int j=1; j<=atom->ntypes; j++) {
      cutsq[i][j] = par.cutmax * par.cutmax;
    }
  }
}

// ########################################################
// ########################################################

// ########################################################
//                       INIT_STYLE
// ########################################################
// Set up the pair style to be a NN potential.

void PairPANNA::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style PANNA requires newton pair on");

  // request FULL neighbor list
  int list_style = NeighConst::REQ_FULL;
  neighbor->add_request(this, list_style);

}


// ########################################################
// ########################################################

// ########################################################
//                       init_one
// ########################################################
// Initilize 1 pair interaction.  Needed by LAMMPS but not
// used in this style.

double PairPANNA::init_one(int i, int j)
{
  return sqrt(cutsq[i][j]);
}

// ########################################################
// ########################################################



// ########################################################
//                       WRITE_RESTART
// ########################################################
// Writes restart file. Not implemented.

void PairPANNA::write_restart(FILE *fp)
{

}

// ########################################################


// ########################################################
//                       READ_RESTART
// ########################################################
// Reads from restart file. Not implemented.

void PairPANNA::read_restart(FILE *fp)
{

}

// ########################################################


// ########################################################
//                       WRITE_RESTART_SETTINGS
// ########################################################
// Writes settings to restart file. Not implemented.

void PairPANNA::write_restart_settings(FILE *fp)
{

}

// ########################################################
// ########################################################



// ########################################################
//                       READ_RESTART_SETTINGS
// ########################################################
// Reads settings from restart file. Not implemented.

void PairPANNA::read_restart_settings(FILE *fp)
{

}

// ########################################################
// ########################################################

// Not implemented.
void PairPANNA::write_data(FILE *fp)
{

}

// Not implemented.
void PairPANNA::write_data_all(FILE *fp)
{

}

// Not implemented.
double PairPANNA::single(int i, int j, int itype, int jtype, double rsq,
                      double factor_coul, double factor_lj,
                      double &fforce)
{
  return 1;
}

/* ---------------------------------------------------------------------- */



// ########################################################
//                       Settings
// ########################################################
// Initializes settings. No setting needed.

void PairPANNA::settings(int narg, char* argv[])
{
  if (narg != 0) {
    error->all(FLERR,"pair_panna requires no arguments.\n");
  }

}

// ########################################################
// ########################################################
