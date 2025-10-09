/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(runner,PairRuNNer);
// clang-format on
#else

#ifndef LMP_PAIR_RUNNER_H
#define LMP_PAIR_RUNNER_H

#include "pair.h"

extern "C" {
int runner_lammps_api_version();
void runner_lammps_interface_init(const char *path, int *npath, double *cutoff, double *cfenergy,
                                  double *cflength, int *nnp_generation, int *num_committee_members,
                                  bool *l_hirshfeld_vdw, bool *ltwo_body, bool *lcheck_extrap,
                                  int *rank, int *size);

void runner_lammps_interface_transfer_atoms_and_neighbor_lists(
    int *nlocal, int *nghost, int *atomic_numbers, int *inum, int *sum_num_neigh, int *ilist,
    int *num_neigh, int *first_neigh, int *neigh, double *lattice, double *xyz, bool *lperiodic);

void runner_interface_short_range(int *nlocal, int *nghost, int *inum, int *nmax, int *ilist,
                                  double *energy, double *forces, double *d_energy_d_strain,
                                  double *hirsh_volumes, double *atomic_charges,
                                  double *elec_negativities, double *hardness);
// void runner_lammps_interface_hirshfeld_vdw(int*, int*, int*, int*, double*,
// double*, double*, double*, double*);

void runner_interface_evaluate_electrostatics_3g_part_1(int *num_atoms, double *xyz,
                                                        double *total_charge, double *lattice,
                                                        bool *lperiodic, double *atomic_charges,
                                                        double *energy, double *force_outer,
                                                        double *d_energy_d_q,
                                                        double *d_energy_d_strain);

void runner_interface_reinitialize_electrostatics(int *num_atoms_global, double *pos_global,
                                                  int *atomic_numbers_global);

void runner_interface_calc_screening(int *nlocal, int *nghost, double *atomic_charges,
                                     double *energy, double *forces, double *de_dq,
                                     double *d_energy_d_strain);

void runner_interface_evaluate_electrostatics_3g_part_2(int *nlocal, int *nghost, int *nglobal,
                                                        int *committee_member_idx, double *energy,
                                                        double *forces, double *d_energy_d_q,
                                                        double *d_energy_d_q_sum_global,
                                                        double *d_energy_d_strain);

void runner_interface_compute_charges_4g(int *num_atoms, double *total_charge,
                                         double *atomic_electronegativities,
                                         double *atomic_hardness, double *atomic_charges,
                                         bool *luse_prev_q, int icomm);

void runner_interface_finalize_step();

void runner_interface_short_range_4g(int *nlocal, int *nghost, int *inum, int *nmax, int *ilist,
                                     double *, double *energy, double *forces,
                                     double *d_energy_d_strain, double *d_energy_d_q);

void runner_interface_evaluate_electrostatics_4g_part_1(int *nglobal, double *d_energy_dq,
                                                        double *energy, double *forces,
                                                        double *d_energy_d_strain,
                                                        double *lagrange_charges, int icomm);

void runner_interface_evaluate_electrostatics_4g_part_2(int *nlocal, int *nghost,
                                                        int *committee_member_idx,
                                                        double *lagrange_charges, double *forces,
                                                        double *d_energy_d_strain);

void runner_interface_hirshfeld_vdw(int *nlocal, int *nghost, int *inum, int *ilist,
                                    int *committee_member_idx, double *hirsh_volumes,
                                    double *energy, double *forces, double *d_energy_d_strain);

void runner_interface_two_body(int *nlocal, int *nghost, double *energy, double *forces,
                               double *d_energy_d_strain);

void runner_interface_extrapolation_warnings(char **c_ptr_extrap_msg, long *len_extrap_msg);

void runner_interface_dealloc_extrapolation_warnings();

void runner_interface_extrapolation_count(long *extraplation_count, long *total_extrapolation_count,
                                          bool *lreset);
}

namespace LAMMPS_NS {

class PairRuNNer : public Pair {
 public:
  /**
   * Constructor of the pair style.
   */
  PairRuNNer(class LAMMPS *);

  /**
   * // Destructor of the pair style.
   */
  ~PairRuNNer() override;

  /**
   * The main computation routine for the pair style.
   * @param eflag Whether the energy needs to be tallied, i.e. summed
   *   across processes. This can often be skipped as the potential
   *   energy is not necessary to drive an MD simulation and often only
   *   printed for diagnostic purposes.
   * @param vflag Whether the virial needs to be tallied, i.e. summed across
   *   processes. The virial is mostly used for pressure calculation, and
   *   therefore only needed when diagnostic output has been requested in
   *   a timestep.
   */
  void compute(int eflag, int vflag) override;

  /**
   * Parse global settings from the pair_style line of the LAMMPS input file.
   * @param narg The number of arguments. This includes parameters and their
   *   values, so a line like `pair_style runner dir "."` gives narg = 2.
   * @param arg A list of the arguments in the line. Total of `narg`.
   */
  void settings(int narg, char **arg) override;

  /**
   * Parse atom type coefficients from the pair_coeff line of the LAMMPS input file.
   * @param narg The number of arguments. This includes parameters and their
   *   values, so a line like `pair_style runner dir "."` gives narg = 2.
   * @param arg A list of the arguments in the line. Total of `narg`.
   */
  void coeff(int narg, char **arg) override;

  /**
   * Init specific to this pair style.
   */
  void init_style() override;

  /**
   * Init for one type pair i,j and corresponding j,i
   * This function is called in the init phase of the simulation
   * It returns the cutoff, which is then used by LAMMPS for the
   * neighborlist calculation
   */
  double init_one(int, int) override;

  /**
   * Set coefficients for one or more type pairs.
   */
  void allocate();

  /**
   * Pack local atom information into communication buffer.
   * Before calling this routine, the global `commstyle` variable
   * should be adjusted to the desired property (see further below).
   */
  int pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc) override;

  /**
   * Unpack local atom information from buffer into ghost atom storage.
   * Before calling this routine, the global `commstyle` variable
   * should be adjusted to the desired property (see further below).
   */
  void unpack_forward_comm(int n, int first, double *buf) override;

  /**
   * Pack ghost atom contributions into communication buffer.
   */
  int pack_reverse_comm(int n, int first, double *buf) override;

  /**
   * Add ghost atom contributions from communication buffer to local atom.
   */
  void unpack_reverse_comm(int n, int *list, double *buf) override;

  /**
   * Pack the xyz coordinates from the local atoms on each process into one
   * global structure on root.
   * This is required for electrostatics calculations.
   * @param rank The MPI rank of this process.
   * @param size The number of MPI processes.
   * @param natoms The total number of atoms in the simulation box.
   * @param inum The number of local atoms on this process for which a neighbor
   *   list has been built.
   * @param ilist The indices of those atoms owned by this process for which a
   *   neighbor list has been built. Has a total of `inum` members.
   * @param x The Cartesian coordinates of the `inum` local atoms on this
   *   process for which a neighbor list has been built.
   * @param runner_types
   * @param xyz_global The Cartesian coordinates of all `natoms` atoms. Only
   *   filled on the process with MPI rank 0.
   * @param z_global The atomic numbers of all `natoms` atoms. Only populated on
   *   the root process.
   */
  void pack_structure(int rank, int size, int natoms, int inum, int *ilist, tagint *tag, double **x,
                      int *runner_types, double *&xyz_global, int *&z_global);

  /**
   * Pack the values of a single atomic property from the local atoms
   * on each process into a single array on the root MPI process.
   * @param rank The MPI rank of this process.
   * @param size The number of MPI processes.
   * @param natoms The total number of atoms in the simulation box.
   * @param inum The number of local atoms on this process for which a neighbor
   *   list has been built.
   * @param ilist The indices of those atoms owned by this process for which a
   *   neighbor list has been built. Has a total of `inum` members.
   * @param local_property The property value for each of the `inum` local atoms
   *   on this process for which a neighbor list has been built.
   * @param global_property Collected property values for all `natoms` atoms in
   *   the simulation box. Only populated on the root MPI process.
   */
  void pack_atomic_property(int rank, int size, int natoms, int inum, int *ilist, tagint *tag,
                            double *local_property, double *&global_property);

  /**
   * Distribute the values of a single atomic property from the root process to
   * all local atoms owned by the child processes.
   * @param rank The MPI rank of this process.
   * @param size The number of MPI processes.
   * @param natoms The total number of atoms in the simulation box.
   * @param inum The number of local atoms on this process for which a neighbor
   *   list has been built.
   * @param ilist The indices of those atoms owned by this process for which a
   *   neighbor list has been built. Has a total of `inum` members.
   * @param nprop The number of property values per atom.
   * @param global_property The property values for all `natoms` atoms in
   *   the simulation box. Only needs to be populated on the root MPI process.
   * @param local_property The distributed property value for each of the `inum` local atoms
   *   on this process for which a neighbor list has been built.
   */
  void unpack_local_atomic_properties(int rank, int size, int natoms, int inum, int *ilist,
                                      tagint *tag, int nprop, double *&global_properties,
                                      double *local_properties);

 private:
  double cflength;    // Length conversion factor.
  double cfenergy;    // Energy conversion factor.
  bool
      luse_prev_q;    // Use charges from previous timestep as initial guess for iterative qeq solvers.
  bool lwrite_f_comm;    // Write committee forces into f_comm array
  bool lcheck_extrap;    // Flag enabling checks for feature extrapolation
  long
      max_extrap;    // Maximal number of allowed timesteps with feature extrapolations during the MD simulation
  bool lshow_ew;         // Flag enabling output of extrapolation warnings to log file
  long sum_ew_freq;      // Frequency where extrapolation warning summary is printed to log file
  long reset_ew_freq;    // Frequency where extrapolation count is reseted to 0
  long
      local_extrap_sum;    // Sum of recorded extrapolations per process (being reseted at every extrapolation summary)
  double cutoff;    // Max feature map cutoff.
  double
      total_charge;    // The total charge of the structure being simulated. Must be 0 for periodic systems.
  char *directory;    // directory containing RuNNer potential files
  int *map;           // Mapping from atom types to elements
  int nmax;           // Allocated size of per-atom arrays. This is usually
      // max(nlocal + nghost) across all processes, with some padding that can become quite large.

  // Additional per-atom arrays
  double *atomic_charge, *hirshfeld_volume, *electronegativity, *lagrange_charges, *de_dq,
      *screening_de_dq, *hardness, *committee_storage;
  bool lhirshfeld_vdw;
  bool ltwo_body;
  int nnp_generation;
  int num_committee_members;    // specified in input.nn
  int commstyle;                // communication flag for forward and reverse
                                // communication
  const int COMMATCHARGE = 1;
  const int COMMELECNEGATIVITY = 2;
  const int COMMHIRSHVOLUME = 3;
  const int COMMDEDQ = 4;
  const int COMMSCREENINGDEDQ = 5;
  const int COMMLAMBDACHARGE = 6;
  const int COMMCOMMITTEESTORAGE = 7;
};

}    // namespace LAMMPS_NS

#endif
#endif
