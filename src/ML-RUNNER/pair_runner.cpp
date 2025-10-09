/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Knut Nikolas Lausch, Gunnar Schmitz,
                        Alexander Knoll
------------------------------------------------------------------------- */

#include "pair_runner.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;

PairRuNNer::PairRuNNer(LAMMPS *lmp) : Pair(lmp)
{
  // HDNNP is not pairwise additive, due to three body terms
  single_enable = 0;
  // Do not write binary restart files
  restartinfo = 0;
  // Parameters are read from input.nn, therefore pair_coeff only has single
  // command.
  one_coeff = 1;
  // Many-body potential flag.
  manybody_flag = 1;
  // Currently no unit conversion necessary as RuNNer is unit-agnostic.
  unit_convert_flag = 0;
  // We generally calculate the virial ourselves
  // and do not need to call virial_fdotr().
  // In case of 2G HDNNP, we could call virial_fdotr()
  // for performace reasons.
  // flag must then be overwritten in init_style().
  no_virial_fdotr_compute = 1;
  map = nullptr;

  // Additional per-atom arrays for communication
  nmax = 0;
  atomic_charge = nullptr;
  hirshfeld_volume = nullptr;
  electronegativity = nullptr;
  lagrange_charges = nullptr;
  screening_de_dq = nullptr;
  de_dq = nullptr;
  committee_storage = nullptr;

  // Set commsize needed by this pairstyle.
  comm_forward = 1;    // Forward communication (1 double per atom)
  comm_reverse = 1;    // Reverse communication (1 double per atom)
  commstyle = 0;

  // Sum of extrapolation is zero upon initialization of this pair style
  local_extrap_sum = 0;

  // variable for output using pair compute
  nextra = 0;
  pvector = nullptr;
}

PairRuNNer::~PairRuNNer()
{
  // Deallocate member variables
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(atomic_charge);
    memory->destroy(hirshfeld_volume);
    memory->destroy(electronegativity);
    memory->destroy(lagrange_charges);
    memory->destroy(de_dq);
    memory->destroy(screening_de_dq);
    memory->destroy(committee_storage);
    delete[] map;
    delete[] directory;
    delete[] pvector;
  }
}

void PairRuNNer::compute(int eflag, int vflag)
{
  const bool debug = false;

  int inum, jnum, ii, jj, i, j;
  int *ilist;
  int *jlist;
  int *numneigh, **firstneigh;
  int *runner_num_neigh, *runner_first_neighbor, *runner_jlist, *runner_types;

  // Initializes flags, which signal if energy and virial need to be tallied.
  // In a typical MD simulation, these two properties are only required for
  // output purposes (reporting potential energy and pressure), so their
  // collection across all processes (=tallying) can be skipped in all timesteps
  // where no reporting has been requested.
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  // Number of atoms owned by this process.
  int nlocal = atom->nlocal;
  // Number of ghost atoms on this process. Ghost atoms are atoms owned by another
  // process which are inside the neighbor list cutoff for local atoms on this
  // process. We have to calculate the force contribution from local atoms to
  // these ghost atoms.
  int nghost = atom->nghost;
  int ntotal = nlocal + nghost;
  // Total number of atoms in the simulation box. This is mpi_sum(nlocal) across
  // all processes.
  int natoms = static_cast<int>(atom->natoms);
  int *type = atom->type;
  tagint *tag = atom->tag;

  // Interface variables
  bool lperiodic;
  double *committee_energy, *committee_force, *committee_d_energy_d_strain, *lattice;
  double *committee_atomic_charge, *committee_hirshfeld_volume, *committee_electronegativity,
      *committee_hardness;

  lattice = new double[9];

  // MPI
  int size, rank;
  rank = comm->me;
  size = comm->nprocs;

  if (nlocal == 0) {
    std::cout << "Error in PairRuNNer. No local atoms on process " << rank << "." << std::endl;
    std::cout << "Try adjusting simulation box partitioning with the "
                 "`balance` command"
              << std::endl;
    std::cout << "or restart the simulation using fewer processors." << std::endl;
    MPI_Abort(world, -1);
  }

  MPI_Barrier(world);
  if (debug) std::cout << "Entered PairRuNNer::compute" << std::endl;

  // Allocate additional per-atom arrays
  if (atom->nmax > nmax) {
    memory->destroy(atomic_charge);
    memory->destroy(hirshfeld_volume);
    memory->destroy(electronegativity);
    memory->destroy(lagrange_charges);
    memory->destroy(de_dq);
    memory->destroy(screening_de_dq);
    memory->destroy(committee_storage);
    nmax = atom->nmax;

    memory->create(atomic_charge, nmax, "pair:atomic_charge");
    memory->create(hirshfeld_volume, nmax, "pair:hirshfeld_volume");
    memory->create(electronegativity, nmax, "pair:electronegativity");
    memory->create(lagrange_charges, nmax, "pair:lagrange_charges");
    memory->create(de_dq, nmax, "pair:de_dq");
    memory->create(screening_de_dq, nmax, "pair:screening_de_dq");
    memory->create(committee_storage, nmax, "pair:committee_storage");
  }

  // Set additional per-atom arrays to zero
  memset(atomic_charge, 0.0, nmax * (sizeof *atomic_charge));
  memset(hirshfeld_volume, 0.0, nmax * (sizeof *atomic_charge));
  memset(electronegativity, 0.0, nmax * (sizeof *electronegativity));
  memset(lagrange_charges, 0.0, nmax * (sizeof *lagrange_charges));
  memset(de_dq, 0.0, nmax * (sizeof *de_dq));
  memset(screening_de_dq, 0.0, nmax * (sizeof *screening_de_dq));
  memset(committee_storage, 0.0, nmax * (sizeof *committee_storage));

  committee_energy = new double[num_committee_members];
  committee_force = new double[ntotal * 3 * num_committee_members];
  committee_d_energy_d_strain = new double[9 * num_committee_members];
  committee_atomic_charge = new double[nmax * num_committee_members];
  committee_hirshfeld_volume = new double[nmax * num_committee_members];
  committee_electronegativity = new double[nmax * num_committee_members];
  committee_hardness = new double[nmax * num_committee_members];

  // Set interface variables to zero
  memset(committee_energy, 0.0, num_committee_members * (sizeof *committee_energy));
  memset(committee_force, 0.0, (ntotal * 3 * num_committee_members) * (sizeof *committee_force));
  memset(committee_d_energy_d_strain, 0.0,
         9 * num_committee_members * (sizeof *committee_d_energy_d_strain));
  memset(committee_atomic_charge, 0.0,
         nmax * num_committee_members * (sizeof *committee_atomic_charge));
  memset(committee_hirshfeld_volume, 0.0,
         nmax * num_committee_members * (sizeof *committee_hirshfeld_volume));
  memset(committee_electronegativity, 0.0,
         nmax * num_committee_members * (sizeof *committee_electronegativity));
  memset(committee_hardness, 0.0, nmax * num_committee_members * (sizeof *committee_hardness));

  // Neighborlist information
  // Number of local atoms for which a neighbor list has been built
  // This is usually equal to nlocal, but may also be smaller when some
  // atoms have been removed from the neighbor list calculation by the user.
  inum = list->inum;
  // Local index of those inum atoms for which a neighbor list was built.
  ilist = list->ilist;
  // The number of neighbors for each atom in ilist (dimension inum).
  numneigh = list->numneigh;
  // Pointer array to the first neighbor of each atom in ilist (dimension inum).
  firstneigh = list->firstneigh;

  // Calculate total number of neighbors
  int num_neigh_sum = 0;
  for (ii = 0; ii < inum; ii++) num_neigh_sum += numneigh[ii];

  // create  arrays for interface
  runner_num_neigh = new int[inum];
  runner_first_neighbor = new int[inum];
  runner_jlist = new int[num_neigh_sum];
  runner_types = new int[ntotal];

  // Collect neighbor data
  int irunner = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];            // get local index of atom ii
    jlist = firstneigh[i];    // pointer to first neighbor of atom
    jnum = numneigh[i];       // number of neighbors jj of atom
    runner_num_neigh[ii] = jnum;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];     // index of neighbor jj
      j &= NEIGHMASK;    // masks bits, which encode pair information
      if (jj == 0) runner_first_neighbor[ii] = irunner + 1;    // plus one due to Fortran
      runner_jlist[irunner] = j;
      irunner++;    // runs till num_neigh_sum
    }
  }
  // collect atomic numbers of all atoms by converting element id to atomic
  // number using map
  for (ii = 0; ii < ntotal; ii++) runner_types[ii] = map[type[ii]];

  // get lattice parameters
  lattice[0] = domain->xprd;
  lattice[1] = 0.0;
  lattice[2] = 0.0;
  lattice[3] = domain->xy;
  lattice[4] = domain->yprd;
  lattice[5] = 0.0;
  lattice[6] = domain->xz;
  lattice[7] = domain->yz;
  lattice[8] = domain->zprd;

  // get periodicity information
  switch (domain->periodicity[0] + domain->periodicity[1] + domain->periodicity[2]) {
    case 3:
      lperiodic = true;
      break;
    case 0:
      lperiodic = false;
      break;
    default:
      if (nnp_generation == 2) {
        lperiodic = false;
        break;
      } else {
        error->all(
            FLERR,
            "Periodicity in only one or two dimensions not supported when using 3G or 4G HDNNPs,");
      }
  }

  if (lperiodic && total_charge != 0.0)
    error->all(FLERR,
               "Periodic systems must be charge neutral (total_charge = 0.0) when using pair_style "
               "runner.");

  if (debug) std::cout << "Transfer atoms and neighbor lists to RuNNer interface" << std::endl;

  runner_lammps_interface_transfer_atoms_and_neighbor_lists(
      &nlocal, &nghost, runner_types, &inum, &num_neigh_sum, ilist, runner_num_neigh,
      runner_first_neighbor, runner_jlist, lattice, &x[0][0], &lperiodic);

  if (debug) std::cout << "RuNNer short-range predicition" << std::endl;

  runner_interface_short_range(&nlocal, &nghost, &inum, &nmax, ilist, committee_energy,
                               committee_force, committee_d_energy_d_strain,
                               committee_hirshfeld_volume, committee_atomic_charge,
                               committee_electronegativity, committee_hardness);

  if (lhirshfeld_vdw) {
    if (debug) std::cout << "RuNNer long-range vdW interactions" << std::endl;

    for (int i = 0; i < num_committee_members; i++) {

      double vdw_energy = 0.0;

      double *vdw_forces = new double[ntotal * 3];
      memset(vdw_forces, 0.0, ntotal * 3 * (sizeof *vdw_forces));

      double *vdw_d_energy_d_strain = new double[9];
      memset(vdw_d_energy_d_strain, 0.0, 9 * (sizeof *vdw_d_energy_d_strain));

      for (int ii = 0; ii < nmax; ii++)
        hirshfeld_volume[ii] = committee_hirshfeld_volume[nmax * i + ii];

      // Communicate Hirshfeld volumes from local atoms to ghost atoms.
      commstyle = COMMHIRSHVOLUME;
      comm->forward_comm(this);

      // Calculate dispersion energies and forces using Hirshfeld volumes
      // and volume gradients (stored on runner side)
      runner_interface_hirshfeld_vdw(&nlocal, &nghost, &inum, ilist, &i, hirshfeld_volume,
                                     &vdw_energy, vdw_forces, vdw_d_energy_d_strain);

      // Add electrostatic interactions to short-range results
      committee_energy[i] += vdw_energy;

      for (ii = 0; ii < ntotal * 3; ii++) committee_force[ntotal * 3 * i + ii] += vdw_forces[ii];

      for (ii = 0; ii < 9; ii++)
        committee_d_energy_d_strain[i * 9 + ii] += vdw_d_energy_d_strain[ii];

      delete[] vdw_forces;
      delete[] vdw_d_energy_d_strain;
    }
  }

  if (ltwo_body) {
    if (debug) std::cout << "RuNNer two-body potential" << std::endl;

    double two_body_energy = 0.0;

    double *two_body_forces = new double[ntotal * 3];
    memset(two_body_forces, 0.0, ntotal * 3 * (sizeof *two_body_forces));

    double *two_body_d_energy_d_strain = new double[9];
    memset(two_body_d_energy_d_strain, 0.0, 9 * (sizeof *two_body_d_energy_d_strain));

    // Calculate two-body energies and forces
    runner_interface_two_body(&nlocal, &nghost, &two_body_energy, two_body_forces,
                              two_body_d_energy_d_strain);

    for (int i = 0; i < num_committee_members; i++) {
      // Add two-body interactions to short-range results
      committee_energy[i] += two_body_energy;

      for (ii = 0; ii < ntotal * 3; ii++)
        committee_force[ntotal * 3 * i + ii] += two_body_forces[ii];

      for (ii = 0; ii < 9; ii++)
        committee_d_energy_d_strain[i * 9 + ii] += two_body_d_energy_d_strain[ii];

      delete[] two_body_forces;
      delete[] two_body_d_energy_d_strain;
    }
  }

  if (nnp_generation >= 3) {
    // Collect local atoms into a global structure. This needs to be done only
    // once for all committee members.
    double *xyz_global = new double[natoms * 3];
    int *z_global = new int[natoms];

    pack_structure(rank, size, natoms, inum, ilist, tag, x, runner_types, xyz_global, z_global);

    if (rank == 0) {
      // Reinitialize the electrostatics calculator on the root
      // process. It is completely reallocated if the global
      // number of atoms changed. Otherwise, only the ordering
      // of atoms is updated.
      runner_interface_reinitialize_electrostatics(&natoms, &xyz_global[0], z_global);
    }

    if (nnp_generation == 3) {
      if (debug) std::cout << "RuNNer 3G long-range electrostatics" << std::endl;

      for (int i = 0; i < num_committee_members; i++) {

        // Collect the charge predicted by the committee members into one global array.
        double *q_global = new double[natoms];
        pack_atomic_property(rank, size, natoms, inum, ilist, tag,
                             &committee_atomic_charge[nmax * i], q_global);

        // long-range electrostatics variables
        double runner_elec_energy = 0.0;
        double *elec_force_global = new double[natoms * 3];
        memset(elec_force_global, 0.0, natoms * 3 * (sizeof *elec_force_global));

        double *de_dq_global = new double[natoms];
        memset(de_dq_global, 0.0, natoms * (sizeof *de_dq_global));

        double *runner_elec_d_energy_d_strain = new double[9];
        memset(runner_elec_d_energy_d_strain, 0.0, 9 * (sizeof *runner_elec_d_energy_d_strain));

        if (rank == 0) {
          // Calculate long-range electrostatics on root using the global
          // structure.
          runner_interface_evaluate_electrostatics_3g_part_1(
              &natoms, &xyz_global[0], &total_charge, lattice, &lperiodic, &q_global[0],
              &runner_elec_energy, &elec_force_global[0], &de_dq_global[0],
              &runner_elec_d_energy_d_strain[0]);
        }

        MPI_Barrier(world);

        // Broadcast and unpack electrostatic results back to each local
        // process.
        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 1, q_global,
                                       atomic_charge);

        memset(de_dq, 0.0, ntotal * (sizeof *de_dq));
        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 1, de_dq_global,
                                       de_dq);

        double *runner_elec_forces = new double[ntotal * 3];
        memset(runner_elec_forces, 0.0, ntotal * 3 * (sizeof *runner_elec_forces));
        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 3, elec_force_global,
                                       runner_elec_forces);

        // Communicate scaled atomic charges from local atoms to ghost
        // atoms for screening calculation.
        commstyle = COMMATCHARGE;
        comm->forward_comm(this);

        double screening_energy = 0.0;
        double *screening_forces = new double[ntotal * 3];
        memset(screening_forces, 0.0, ntotal * 3 * (sizeof *screening_forces));

        memset(screening_de_dq, 0.0, ntotal * (sizeof *screening_de_dq));

        double *screening_d_energy_d_strain = new double[9];
        memset(screening_d_energy_d_strain, 0.0, 9 * (sizeof *screening_d_energy_d_strain));

        // Apply screening
        runner_interface_calc_screening(&nlocal, &nghost, atomic_charge, &screening_energy,
                                        screening_forces, screening_de_dq,
                                        screening_d_energy_d_strain);

        // Communicate screening de_dq from ghost atoms to local atoms
        commstyle = COMMSCREENINGDEDQ;
        comm->reverse_comm(this);

        // Determine sum of de_dq of local atoms across all procs.
        double de_dq_sum_local = 0.0;
        for (j = 0; j < inum; j++) {
          ii = ilist[j];
          de_dq[ii] -= screening_de_dq[ii];
          de_dq_sum_local += de_dq[ii];
        }
        double de_dq_sum_global = 0.0;
        MPI_Allreduce(&de_dq_sum_local, &de_dq_sum_global, 1, MPI_DOUBLE, MPI_SUM, world);

        runner_interface_evaluate_electrostatics_3g_part_2(
            &nlocal, &nghost, &natoms, &i, &runner_elec_energy, runner_elec_forces, de_dq,
            &de_dq_sum_global, runner_elec_d_energy_d_strain);

        // Add electrostatic interactions to short-range results
        committee_energy[i] += runner_elec_energy - screening_energy;

        for (ii = 0; ii < ntotal * 3; ii++)
          committee_force[ntotal * 3 * i + ii] += runner_elec_forces[ii] - screening_forces[ii];

        for (ii = 0; ii < 9; ii++)
          committee_d_energy_d_strain[i * 9 + ii] +=
              runner_elec_d_energy_d_strain[ii] - screening_d_energy_d_strain[ii];

        delete[] screening_forces;
        delete[] screening_d_energy_d_strain;
        delete[] elec_force_global;
        delete[] de_dq_global;
        delete[] runner_elec_forces;
        delete[] runner_elec_d_energy_d_strain;
        delete[] q_global;
      }

    } else if (nnp_generation == 4) {
      if (debug) std::cout << "RuNNer 4G non-local" << std::endl;

      // Part 1. For each committee member:
      //   - collect predicted electronegativities and hardness values on root.
      //   - perform QeQ to obtain charges for all atoms
      //   - transfer charges back to each process (local and ghost atoms)
      for (int i = 0; i < num_committee_members; i++) {

        // Collect electronegativities and hardness values from each process into
        // one structure on the root process for the computation of charges via
        // QeQ.
        double *electronegativity_global = new double[natoms];
        double *hardness_global = new double[natoms];
        double *q_global = new double[natoms];

        pack_atomic_property(rank, size, natoms, inum, ilist, tag,
                             &committee_electronegativity[nmax * i], electronegativity_global);

        pack_atomic_property(rank, size, natoms, inum, ilist, tag, &committee_hardness[nmax * i],
                             hardness_global);

        if (rank == 0) {
          // compute charges using qeq on root using global structure
          runner_interface_compute_charges_4g(&natoms, &total_charge, &electronegativity_global[0],
                                              &hardness_global[0], &q_global[0], &luse_prev_q,
                                              i + 1);
        }
        MPI_Barrier(world);

        // unpack global charges from root to respective processes
        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 1, q_global,
                                       atomic_charge);

        // Communicate qeq charges from local atoms to ghost atoms
        // for screening calculation
        commstyle = COMMATCHARGE;
        comm->forward_comm(this);

        for (ii = 0; ii < nmax; ii++) committee_atomic_charge[nmax * i + ii] = atomic_charge[ii];

        delete[] electronegativity_global;
        delete[] hardness_global;
        delete[] q_global;
      }

      double *committee_d_energy_d_q = new double[ntotal * num_committee_members];

      // Perform short-range prediction for all committee members at once.
      runner_interface_short_range_4g(&nlocal, &nghost, &inum, &nmax, ilist,
                                      committee_atomic_charge, committee_energy, committee_force,
                                      committee_d_energy_d_strain, committee_d_energy_d_q);
      MPI_Barrier(world);

      // Part 2. For each committee member:
      //   - Calculate electrostatic screening
      //   - sum up short-range and screening dE/dQ
      //   - determine lagrange charges, evaluate electrostatics, and sum up
      //     contributions
      //   - calculate force contributions of dChi/dxyz and dHardness/dxyz.
      for (int i = 0; i < num_committee_members; i++) {

        double screening_energy = 0.0;
        double *screening_forces = new double[ntotal * 3];
        memset(screening_forces, 0.0, ntotal * 3 * (sizeof *screening_forces));

        memset(screening_de_dq, 0.0, ntotal * (sizeof *screening_de_dq));

        double *screening_d_energy_d_strain = new double[9];
        memset(screening_d_energy_d_strain, 0.0, 9 * (sizeof *screening_d_energy_d_strain));

        // Apply screening
        runner_interface_calc_screening(&nlocal, &nghost, &committee_atomic_charge[nmax * i],
                                        &screening_energy, screening_forces, screening_de_dq,
                                        screening_d_energy_d_strain);

        // Communicate screening de_dq from ghost atoms to local atoms
        commstyle = COMMSCREENINGDEDQ;
        comm->reverse_comm(this);

        // Determine sum of de_dq of local atoms across all procs.
        memset(de_dq, 0.0, ntotal * (sizeof *de_dq));
        for (j = 0; j < inum; j++) {
          ii = ilist[j];
          de_dq[ii] += committee_d_energy_d_q[ntotal * i + ii] - screening_de_dq[ii];
        }

        double runner_elec_energy = 0.0;
        double *elec_force_global = new double[natoms * 3];
        memset(elec_force_global, 0.0, natoms * 3 * (sizeof *elec_force_global));

        double *de_dq_global = new double[natoms];
        memset(de_dq_global, 0.0, natoms * (sizeof *de_dq_global));

        double *runner_elec_d_energy_d_strain = new double[9];
        memset(runner_elec_d_energy_d_strain, 0.0, 9 * (sizeof *runner_elec_d_energy_d_strain));

        double *lagrange_global = new double[natoms];
        memset(lagrange_global, 0.0, natoms * (sizeof *lagrange_global));

        // Communicate de_dq and pack into global structure for
        // determination of lagrange charges.

        // communicate ghost atom contributions to local atoms
        commstyle = COMMDEDQ;
        comm->reverse_comm(this);

        pack_atomic_property(rank, size, natoms, inum, ilist, tag, de_dq, de_dq_global);

        if (rank == 0) {
          // serial step determining lagrange charges and
          // electrostatic contribution from global de_dq
          runner_interface_evaluate_electrostatics_4g_part_1(
              &natoms, de_dq_global, &runner_elec_energy, elec_force_global,
              runner_elec_d_energy_d_strain, lagrange_global, i + 1);
        }

        MPI_Barrier(world);

        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 1, lagrange_global,
                                       lagrange_charges);

        // communicate lagrange charges from local atoms to ghost
        // atoms for calculation of force trick part 2
        commstyle = COMMLAMBDACHARGE;
        comm->forward_comm(this);

        double *runner_elec_forces = new double[ntotal * 3];
        memset(runner_elec_forces, 0.0, ntotal * 3 * (sizeof *runner_elec_forces));
        unpack_local_atomic_properties(rank, size, natoms, inum, ilist, tag, 3, elec_force_global,
                                       runner_elec_forces);

        // Apply remaining force contributions from predicited
        // electronegativities and lagrange charges to
        // electrostatic forces.
        runner_interface_evaluate_electrostatics_4g_part_2(&nlocal, &nghost, &i, lagrange_charges,
                                                           runner_elec_forces,
                                                           runner_elec_d_energy_d_strain);

        // Add electrostatic interactions to short-range results
        committee_energy[i] += runner_elec_energy - screening_energy;

        for (ii = 0; ii < ntotal * 3; ii++)
          committee_force[ntotal * 3 * i + ii] += runner_elec_forces[ii] - screening_forces[ii];

        for (ii = 0; ii < 9; ii++)
          committee_d_energy_d_strain[i * 9 + ii] +=
              runner_elec_d_energy_d_strain[ii] - screening_d_energy_d_strain[ii];

        delete[] runner_elec_forces;
        delete[] screening_forces;
        delete[] screening_d_energy_d_strain;
        delete[] elec_force_global;
        delete[] runner_elec_d_energy_d_strain;
        delete[] lagrange_global;
        delete[] de_dq_global;
      }

      delete[] committee_d_energy_d_q;
    }

    delete[] xyz_global;
    delete[] z_global;
  }

  // Copy results from RuNNer back into LAMMPS atom array
  if (debug) std::cout << "Transferring results from RuNNer to LAMMPS" << std::endl;

  // Forces
  irunner = 0;
  for (i = 0; i < num_committee_members; i++) {
    for (ii = 0; ii < ntotal; ii++) {
      for (jj = 0; jj < 3; jj++) {
        f[ii][jj] += committee_force[irunner] / cfenergy * cflength /
            num_committee_members;    // runner_force is a vector
        irunner++;
      }
    }
  }

  // write individual committee forces into f_comm array
  if (lwrite_f_comm) {
    int flag, cols, ghost;
    int idx = atom->find_custom_ghost("f_comm", flag, cols, ghost);
    double **f_comm = atom->darray[idx];
    for (jj = 0; jj < 3; jj++) {                       // xyz components
      for (i = 0; i < num_committee_members; i++) {    // committee members
        irunner = jj + i * ntotal * 3;
        memset(committee_storage, 0.0, nmax * (sizeof *committee_storage));
        for (ii = 0; ii < ntotal; ii++) {    // nlocal + nghost atoms
          committee_storage[ii] = committee_force[irunner];
          irunner += 3;
        }
        commstyle = COMMCOMMITTEESTORAGE;
        comm->reverse_comm(
            this);    // communicate forces so that local atoms have full contribution
        for (ii = 0; ii < ntotal; ii++)
          f_comm[ii][jj + i * 3] =
              committee_storage[ii] / cfenergy * cflength;    // copy into customd d2_array
      }
    }
  }

  // Potential energy
  if (eflag_global) eng_vdwl = 0.0;
  for (i = 0; i < num_committee_members; i++) {
    if (eflag_global) eng_vdwl = eng_vdwl + committee_energy[i] / cfenergy / num_committee_members;
  }

  // Write committee energies into pair compute vector
  memset(pvector, 0.0, num_committee_members);
  for (i = 0; i < num_committee_members; i++) pvector[i] = committee_energy[i] / cfenergy;

  // Charges if charge atom style is used
  if (q != NULL) {
    memset(
        q, 0.0,
        nmax *
            (sizeof *q));    // This array does not seem to get reset to zero every timestep by LAMMPS
    for (ii = 0; ii < ntotal; ii++) {
      for (jj = 0; jj < num_committee_members; jj++) {
        q[ii] += committee_atomic_charge[ii + jj * nmax] / num_committee_members;
      }
    }
  }

  // In case of 2G HDNNP, virial could be calculated via F dot r.
  if (vflag_fdotr) virial_fdotr_compute();

  // Virial is -1.0 * d_energy_d_strain
  if (vflag_global) {
    memset(virial, 0.0, 9);
    for (i = 0; i < num_committee_members; i++) {
      virial[0] -= committee_d_energy_d_strain[0 + 9 * i] / cfenergy / num_committee_members;
      virial[1] -= committee_d_energy_d_strain[4 + 9 * i] / cfenergy / num_committee_members;
      virial[2] -= committee_d_energy_d_strain[8 + 9 * i] / cfenergy / num_committee_members;
      virial[3] -= committee_d_energy_d_strain[0 + 9 * i] / cfenergy / num_committee_members;
      virial[4] -= committee_d_energy_d_strain[6 + 9 * i] / cfenergy / num_committee_members;
      virial[5] -= committee_d_energy_d_strain[7 + 9 * i] / cfenergy / num_committee_members;
    }
  }

  // Handling of feature extrapolations
  // - Adding extrapolation warnings to log and screen if `show_extrap` is set to yes.
  // - Stops the MD simulation if the number of extrapolations exceeds
  // the threshold defined by `max_extrap`.
  // - Resets the number of total number of recorded extrapolations during a simulation if
  // the timestep is a multiple of `reset_ew_freq` and larger than 0.
  // - Prints a summary of the recorded extrapolations at every intervall until the timestep is
  // a multiple of `sum_ew_freq` and larger than 0.
  if (lcheck_extrap) {

    long timestep = update->ntimestep;

    // Add extrapolation warnings to log
    if (lshow_ew) {
      char *c_ptr_extrap_msg = nullptr;
      long len_extrap_msg = 0;

      runner_interface_extrapolation_warnings(&c_ptr_extrap_msg, &len_extrap_msg);

      long *len_extrap_msg_array = nullptr;
      if (rank == 0) {
        len_extrap_msg_array = new long[size];

        // Write extrapolation message owned by rank 0
        if (len_extrap_msg > 0) {

          // Write header for extrapolation warnings
          utils::logmesg(lmp, "RuNNer2: Feature extrapolations on process 0 at timestep {:d}\n",
                         timestep);
          std::string extrap_msg = std::string(c_ptr_extrap_msg);
          utils::logmesg(lmp, extrap_msg);
        }
      }

      // Gather length of the extrapolation messages owned by the other processes.
      MPI_Gather(&len_extrap_msg, 1, MPI_LONG, len_extrap_msg_array, 1, MPI_LONG, 0, world);

      // Write extrapolation message owned by the other processes.
      if (rank == 0) {

        // Loop over the processes
        for (int irank = 1; irank < size; irank++) {

          if (len_extrap_msg_array[irank] > 0) {

            // Receive extrapolation message from irank if there is one
            char *c_extrap_msg = nullptr;
            c_extrap_msg = new char[len_extrap_msg_array[irank]];
            MPI_Status mpi_stat;
            MPI_Recv(c_extrap_msg, len_extrap_msg_array[irank], MPI_BYTE, irank, 0, world,
                     &mpi_stat);

            // Write header for extrapolation warnings
            utils::logmesg(lmp,
                           "RuNNer2: Feature extrapolations on process {:d} at timestep {:d}\n",
                           irank, timestep);

            // Convert C char array to string and write extrapolation warning to screen and log
            std::string extrap_msg = std::string(c_extrap_msg);
            utils::logmesg(lmp, extrap_msg);
            delete[] c_extrap_msg;
          }
        }

      } else if (len_extrap_msg > 0) {

        MPI_Send(c_ptr_extrap_msg, len_extrap_msg, MPI_BYTE, 0, 0, world);
      }

      c_ptr_extrap_msg = nullptr;
      if (rank == 0) delete[] len_extrap_msg_array;
    }

    // Total number of extrapolation accumulated on this process during the simulation
    long local_extrap_count_total = 0;
    // Number of extrapolation accumulated on this process during this this timestep
    long extrap_count_timestep = 0;

    // Sets the flag `lreset` to reset the total extrapolation count if the timestep
    // is a multiple of reset_ew_freq and larger than zero
    bool lreset = false;
    if (reset_ew_freq > 0 && timestep % reset_ew_freq == 0 && timestep > 0) { lreset = true; }

    // Retrive the number of extrapolations during this timestep and during the whole simulation
    // on each process 1and reset the latter if `lreset` is true.
    runner_interface_extrapolation_count(&extrap_count_timestep, &local_extrap_count_total,
                                         &lreset);

    // Number of extrapolations recorded on this process (being reseted at every summary)
    local_extrap_sum = local_extrap_sum + extrap_count_timestep;

    // Total number of extrapolation accumulated across all processes "..."
    long global_extrap_count_total = 0;
    MPI_Reduce(&local_extrap_count_total, &global_extrap_count_total, 1, MPI_LONG, MPI_SUM, 0,
               world);

    // Abort simulation if the total number of extrapolations across all processes exceeds `max_extrap`
    if (max_extrap > -1 && global_extrap_count_total > max_extrap && rank == 0) {
      error->one(
          FLERR,
          "Maximal number of allowed extrapolations have been exceeded during the simulation!\n"
          "Current extrapolation count: {:10.3e}",
          double(global_extrap_count_total));
    }

    // Prints a summary of the recorded extrapolations at every intervall until the timestep is
    // a multiple of `sum_ew_freq` and larger than 0.
    if (sum_ew_freq > 0 && timestep % sum_ew_freq == 0 && timestep > 0) {

      // Number of extrapolations recorded on all process for the extrapolation summary
      long global_extrap_sum = 0;
      MPI_Reduce(&local_extrap_sum, &global_extrap_sum, 1, MPI_LONG, MPI_SUM, 0, world);

      if (rank == 0) {
        utils::logmesg(lmp,
                       "RuNNer2: HDNNP Extrapolation Summary\n"
                       "Timestep: {:10d} Sum: {:10.3e} Frequency: {:10.3e} per timestep\n",
                       timestep, double(global_extrap_sum),
                       double(global_extrap_sum) / double(sum_ew_freq));
      }

      local_extrap_sum = 0;
    }

    // Deallocates the character array containing the extrapolation message on the Fortran side
    // and frees up the internal memory of the `ExtrapolationHandler` (see RuNNer 2 documentation)
    runner_interface_dealloc_extrapolation_warnings();
  }

  MPI_Barrier(world);

  // Deallocate internal arrays
  delete[] runner_types;
  delete[] runner_num_neigh;
  delete[] runner_first_neighbor;
  delete[] runner_jlist;
  delete[] committee_energy;
  delete[] committee_force;
  delete[] committee_d_energy_d_strain;
  delete[] committee_atomic_charge;
  delete[] committee_hirshfeld_volume;
  delete[] committee_electronegativity;
  delete[] committee_hardness;
  delete[] lattice;

  runner_interface_finalize_step();

  if (debug) std::cout << "compute done" << std::endl;
}

void PairRuNNer::settings(int narg, char **arg)
{
  int iarg = 0;

  // default settings
  directory = utils::strdup("./");
  cflength = 1.0;
  cfenergy = 1.0;
  luse_prev_q = false;
  lwrite_f_comm = false;
  total_charge = 0.0;
  lcheck_extrap = false;
  max_extrap = 100;
  lshow_ew = false;
  sum_ew_freq = 0;
  reset_ew_freq = 0;
  nextra = 1;    // defaults to committee size of 1

  while (iarg < narg) {
    // set RuNNer potential directory
    if (strcmp(arg[iarg], "dir") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      delete[] directory;
      directory = utils::strdup(arg[iarg + 1]);
      iarg += 2;
      // length unit conversion factor
    } else if (strcmp(arg[iarg], "cflength") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      cflength = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      // energy unit conversion factor
    } else if (strcmp(arg[iarg], "cfenergy") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      cfenergy = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "f_comm") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      lwrite_f_comm = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "committee_size") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      nextra = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "total_charge") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      total_charge = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "use_prev_q") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      luse_prev_q = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "check_extrap") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      lcheck_extrap = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "max_extrap") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      max_extrap = utils::bnumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "show_ew") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      lshow_ew = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "sum_ew_freq") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      sum_ew_freq = utils::bnumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "reset_ew_freq") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal pair_style command");
      reset_ew_freq = utils::bnumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal pair_style command");
  }

  // check if linked to the correct RuNNer library API version
  if (runner_lammps_api_version() != 2)
    error->all(FLERR,
               "RuNNer LAMMPS wrapper API version is not compatible "
               "with this version of LAMMPS");
}

void PairRuNNer::allocate()
{
  allocated = 1;                 // sets allocated flag, checked in coeff function
  int np1 = atom->ntypes + 1;    // +1 because typeID start at 1 not 0

  setflag = memory->create(setflag, np1, np1, "pair:setflag");
  cutsq = memory->create(cutsq, np1, np1, "pair:cutsq");
  map = new int[np1];
}

void PairRuNNer::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ntypes = atom->ntypes;

  // We only have the mapping of the element to type in the pair_coeff
  // line narg 0 and 1 are two wildcards * * for I,J, so mapping
  // declaration starts at 2
  if (narg != (2 + ntypes))
    error->all(FLERR, "Number of arguments {} is not correct, it should be {}", narg, 2 + ntypes);

  // iarg - 1, because we start at 2, so first entry is 1.
  for (int iarg = 2; iarg < narg; iarg++)
    map[iarg - 1] = utils::inumeric(FLERR, arg[iarg], false, lmp);

  // clear setflag since coeff() might be called once with I,J = * *
  for (int iat = 1; iat <= ntypes; iat++)
    for (int jat = iat; jat <= ntypes; jat++) setflag[iat][jat] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  int count = 0;
  for (int iat = 1; iat <= ntypes; iat++)
    for (int jat = iat; jat <= ntypes; jat++)
      if (map[iat] >= 0 && map[jat] >= 0) {
        setflag[iat][jat] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

void PairRuNNer::init_style()
{
  // Require newton pair on
  // Switches reverse communication on on every time step, which adds
  // data from ghost atoms to corresponding local atoms Required for
  // RuNNer forces
  if (force->newton_pair != 1) error->all(FLERR, "Pair style runner requires newton pair on");

  // request full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);

  // Read model coefficients and do initialization on RuNNer side.
  // Returns max cutoff for LAMMPS neighborlist.
  // Also returns booleans if additional atomic properties are
  // predicted, which need to be communicated between local and ghost
  // atoms.
  int n_directory_len = strlen(directory);
  int rank = comm->me;
  int size = comm->nprocs;

  runner_lammps_interface_init(directory, &n_directory_len, &cutoff, &cfenergy, &cflength,
                               &nnp_generation, &num_committee_members, &lhirshfeld_vdw, &ltwo_body,
                               &lcheck_extrap, &rank, &size);

  // In the 2G case, this has a performance benefit when multiple MPI processes
  // are used.
  if (nnp_generation == 2) no_virial_fdotr_compute = 0;    // Overwrite default flag

  // Error checking for output by compute pair command
  if (nextra == num_committee_members) {
    // array for storing committee energies for output by compute pair command
    if (pvector) delete[] pvector;
    pvector = new double[nextra];
  } else {
    error->all(FLERR,
               "committee_size specification for pair style runner "
               "inconsistent between lammps input and input.nn.");
  }
  // Error checking for custom f_comm per-atom array
  if (lwrite_f_comm) {
    int flag, cols, ghost;
    int idx = atom->find_custom_ghost("f_comm", flag, cols, ghost);
    if (idx < 0 || ghost == 0)
      error->all(FLERR, "f_comm array not correctly specified in fix for pair_style runner.");
    if (flag != 1)
      error->all(FLERR, "f_comm array needs to be of type double for pair_style runner.");
    if (cols != num_committee_members * 3)
      error->all(FLERR,
                 "Not enough columns specified for f_comm array in fix for pair_style runner.");
  }
}

double PairRuNNer::init_one(int /*i*/, int /*j*/)
{
  return cutoff;
}

/* ----------------------------------------------------------------------
communication between local and ghost atoms
------------------------------------------------------------------------- */

int PairRuNNer::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m;

  if (commstyle == COMMHIRSHVOLUME) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = hirshfeld_volume[j];
    }
  } else if (commstyle == COMMATCHARGE) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = atomic_charge[j];
    }
  } else if (commstyle == COMMELECNEGATIVITY) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = electronegativity[j];
    }
  } else if (commstyle == COMMDEDQ) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = de_dq[j];
    }
  } else if (commstyle == COMMSCREENINGDEDQ) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = screening_de_dq[j];
    }
  } else if (commstyle == COMMLAMBDACHARGE) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = lagrange_charges[j];
    }
  } else if (commstyle == COMMCOMMITTEESTORAGE) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = committee_storage[j];
    }
  }

  return m;
}

void PairRuNNer::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  if (commstyle == COMMHIRSHVOLUME) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) hirshfeld_volume[i] = buf[m++];
  } else if (commstyle == COMMATCHARGE) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) atomic_charge[i] = buf[m++];
  } else if (commstyle == COMMELECNEGATIVITY) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) electronegativity[i] = buf[m++];
  } else if (commstyle == COMMDEDQ) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) de_dq[i] = buf[m++];
  } else if (commstyle == COMMSCREENINGDEDQ) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) screening_de_dq[i] = buf[m++];
  } else if (commstyle == COMMLAMBDACHARGE) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) lagrange_charges[i] = buf[m++];
  } else if (commstyle == COMMCOMMITTEESTORAGE) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) committee_storage[i] = buf[m++];
  }
}

int PairRuNNer::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  if (commstyle == COMMHIRSHVOLUME) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) buf[m++] = hirshfeld_volume[i];
  } else if (commstyle == COMMATCHARGE) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) buf[m++] = atomic_charge[i];
  } else if (commstyle == COMMELECNEGATIVITY) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) buf[m++] = electronegativity[i];
  } else if (commstyle == COMMDEDQ) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) buf[m++] = de_dq[i];
  } else if (commstyle == COMMSCREENINGDEDQ) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) buf[m++] = screening_de_dq[i];
  } else if (commstyle == COMMLAMBDACHARGE) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) buf[m++] = lagrange_charges[i];
  } else if (commstyle == COMMCOMMITTEESTORAGE) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) buf[m++] = committee_storage[i];
  }

  return m;
}

void PairRuNNer::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  if (commstyle == COMMHIRSHVOLUME) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      hirshfeld_volume[j] += buf[m++];
    }
  } else if (commstyle == COMMATCHARGE) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      atomic_charge[j] += buf[m++];
    }
  } else if (commstyle == COMMELECNEGATIVITY) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      electronegativity[j] += buf[m++];
    }
  } else if (commstyle == COMMDEDQ) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      de_dq[j] += buf[m++];
    }
  } else if (commstyle == COMMSCREENINGDEDQ) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      screening_de_dq[j] += buf[m++];
    }
  } else if (commstyle == COMMLAMBDACHARGE) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      lagrange_charges[j] += buf[m++];
    }
  } else if (commstyle == COMMCOMMITTEESTORAGE) {
    m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      committee_storage[j] += buf[m++];
    }
  }
}

void PairRuNNer::pack_structure(int rank, int size, int natoms, int inum, int *ilist, tagint *tag,
                                double **x, int *runner_types, double *&xyz_global, int *&z_global)
{
  int i, ii;
  int start;
  int tag_minus_one;

  // First collect all Cartesian coordinates.
  double *xyz_local = new double[natoms * 3];
  memset(xyz_local, 0, (natoms * 3) * (sizeof *xyz_local));
  memset(xyz_global, 0, (natoms * 3) * (sizeof *xyz_global));

  double xtmp, ytmp, ztmp;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    // positions of atom i
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // get atom ID of atom i
    tag_minus_one = static_cast<int>(tag[i]) - 1;    // tags start at 1
    start = tag_minus_one * 3;

    // store at atom ID position of local array to get consistent order between timesteps
    xyz_local[start] = xtmp;
    xyz_local[start + 1] = ytmp;
    xyz_local[start + 2] = ztmp;
  }
  MPI_Barrier(world);

  // Communicate local positions to root process.
  MPI_Reduce(xyz_local, xyz_global, natoms * 3, MPI_DOUBLE, MPI_SUM, 0, world);

  // Second collect the atomic numbers z.
  int *z_local = new int[natoms];

  memset(z_local, 0.0, natoms * (sizeof *z_local));
  memset(z_global, 0.0, natoms * (sizeof *z_global));

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // get atom ID of atom i
    tag_minus_one = static_cast<int>(tag[i]) - 1;
    z_local[tag_minus_one] = runner_types[i];
  }
  MPI_Barrier(world);

  // Communicate local atomic numbers to root process.
  MPI_Reduce(z_local, z_global, natoms, MPI_INT, MPI_SUM, 0, world);

  // Deallocation of local arrays.
  delete[] xyz_local;
  delete[] z_local;
}

void PairRuNNer::pack_atomic_property(int rank, int size, int natoms, int inum, int *ilist,
                                      tagint *tag, double *local_property, double *&global_property)
{
  int i, ii;
  int tag_minus_one;

  // Prepare an array of the local properties that is natoms long,
  // where the local atoms of each process are sorted according to their unique atom ID
  // to achieve a consistent order between timesteps.
  double *local_property_sorted = new double[natoms];

  memset(global_property, 0.0, natoms * (sizeof *global_property));
  memset(local_property_sorted, 0.0, natoms * (sizeof *local_property_sorted));

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    tag_minus_one = static_cast<int>(tag[i]) - 1;
    local_property_sorted[tag_minus_one] = local_property[i];
  }
  MPI_Barrier(world);

  // Communicate local charges to root process.
  MPI_Reduce(local_property_sorted, global_property, natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  // Deallocation of local arrays.
  delete[] local_property_sorted;
}

void PairRuNNer::unpack_local_atomic_properties(int rank, int size, int natoms, int inum,
                                                int *ilist, tagint *tag, int nprop,
                                                double *&global_properties,
                                                double *local_properties)
{
  int i, ii, iprop;
  int start;
  int tag_minus_one;

  // Broadcast global property results
  MPI_Bcast(global_properties, natoms * nprop, MPI_DOUBLE, 0, world);

  // sort back from atom ID position in global array to local order
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    tag_minus_one = static_cast<int>(tag[i]) - 1;
    start = tag_minus_one * nprop;

    for (iprop = 0; iprop < nprop; iprop++) {
      local_properties[i * nprop + iprop] = global_properties[start + iprop];
    }
  }
}
