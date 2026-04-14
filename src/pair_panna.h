//###########################################################################
//# Copyright (c), The PANNAdevs group. All rights reserved.                #
//# This file is part of the PANNA code.                                    #
//#                                                                         #
//# The code is hosted on GitLab at https://gitlab.com/PANNAdevs/panna      #
//# For further information on the license, see the LICENSE.txt file        #
//###########################################################################
#ifdef PAIR_CLASS

PairStyle(panna, PairPANNA)

#else

#ifndef LMP_PAIR_PANNA_H
#define LMP_PAIR_PANNA_H

#include "pair.h"
#include <string>
#include <fstream>

namespace LAMMPS_NS {

class PairPANNA : public Pair {
public:
  PairPANNA(class LAMMPS *);
  virtual ~PairPANNA();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);

  int get_parameters(char*, char *);
  int get_input_line(std::ifstream*, std::string*, std::string*);

  void compute_gvect(int, double**, int*, int*, int, double*, double*);
  // double Gradial_d(double, int, double*);
  // void Grad_fill(double, double, double, double , int, int  , int , double*, double*);
  // double Gangular_d(double, double, double, int, int, double*);
  double compute_network(double*, double*, int);
  // double compute_network(double*, double*, int, int);

  struct parameters{
    int Nspecies;
    // Gvector parameters
    double *eta_rad;
    double Rc_rad;
    double Rs0_rad;
    double Rsst_rad;
    int RsN_rad;
    double *eta_ang;
    double Rc_ang;
    double Rs0_ang;
    double Rsst_ang;
    int RsN_ang;
    int *zeta;
    int ThetasN;
    std::string* species;
    int gsize;
    double * Rs_rad;
    double * Rs_ang;
    double * Thetas;

    // Network parameters
    int *Nlayers;
    int **layers_size;
    int **layers_activation;

    // Useful precalculated quantities
    double cutmax;
    double *twoeta_rad;
    double *zeta_half;
    double iRc_rad;
    double iRc_rad_half;
    double iRc_ang;
    double iRc_ang_half;
    double *Rsi_rad;
    double *Rsi_ang;
    double *Thi_cos;
    double *Thi_sin;
    int **typsh;
  };


protected:
  // Gvect and NN parameters
  struct parameters par;
  // The network [species, layers, array]
  double ***network;

  void allocate();

};

}

#endif
#endif
