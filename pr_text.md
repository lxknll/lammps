### Summary

We have implemented a new pair style that provides an interface to our Fortran
library RuNNer (Ruhr University Neural Network Energy Representation). RuNNer is
a program for training and predicting machine learning potentials (MPLs),
specifically high-dimensional neural network potentials (HDNNPs). See
"Implementation Details" for a short overview of HDNNPs.

- the pair style provides access to several different architectures of machine
- learning models. These are categorized into generations:
  - second-generation (2G) HDNNPs are short-ranged many-body potentials. These
    are very similar to the already existing ML-HDNNP package, that is an
    interface to the C++ code n2p2.
  - third-generation (3G) HDNNPs extend upon 2G-HDNNPs by including long-range
    electrostatic interactions.
  - fourth-generation (4G) HDNNPs also include long-range electrostatics, but
    with a global charge equilibration (QEq) step.
  - All generations of HDNNPs can optionally be augmented with long-range
    dispersion interactions and a pairwise repulsive potential.
  - All generations support running multiple potentials at once and using the
    averaged forces (energy, virial) for propagating the simulation
    (committee approach).
- we use `fix/atom_property` for optionally writing the forces of each committee
  member to a file.
- we use `compute/energy` for optionally writing the energies of each committee
  member to a file.

### Related Issue(s)
None

### Author(s)
Knut Nikolas Lausch (1)
Alexander L. M. Knoll (1)
Moritz R. Schäfer (1)
Gunnar Schmitz (1)
Jörg Behler, joerg.behler@rub.de (1)

(1) Lehrstuhl für Theoretische Chemie II, Research Center Chemical Sciences and
Sustainability, Ruhr Universität Bochum, Gesundheitscampus Süd 25, 44801 Bochum,
Germany

### Licensing

By submitting this pull request, I agree that my contribution will be included
in LAMMPS and redistributed under either the GNU General Public License version
2 (GPL v2) or the GNU Lesser General Public License version 2.1 (LGPL v2.1).

### Backward Compatibility
None

### Implementation Notes

#### HDNNPs
In our MLPs, we construct a feature vector for each atom
based on the Cartesian coordinates of it's neighbors. This feature vector is
passed through a machine-learning model, which yields atomic output quantities
for further processing:

- In a 2G-HDNNP, the output quantities are atomic energies. They are summed up to
  yield the total energy of the system. We label this contribution "short-range",
  because it only depends on the local chemical environment of each atom.
- In a 3G-HDNNP, we keep the prediction of atomic energies and add a second model
  which predicts atomic charges (the "electrostatic model"). We use
  these charges to calculate an additional long-range electrostatic energy
  contribution.
- In a 4G-HDNNP, the electrostatic model predicts environment-dependent atomic
  electronegativities. Optionally, users can also predict element-specific or
  environment-dependent atomic hardness values. These values are run through
  a charge equilibration scheme, yielding atomic charges again. In a second step,
  these atomic charges are fed back into a 2G-HDNNP model, so that the predicted
  atomic energies contain a global descriptor in the form of the atomic charges.

#### Initialization
The pair style implements the required initialization routines`init_style`, 
`settings`, `allocate`, and `coeff`. In `init_style` we perform a call to our
library setup routine `runner_lammps_interface_init`. Apart from some general
setup specific to our program, it is notable that this routine also parses
several files (program settings, potential-specific input files) on the root MPI
process and broadcasts the contents. We found this to be very portable and
stable in all our tests.

In case the user requested a 2G-HDNNP, we calculate the virial via
`fdotr`. For potentials with long-range electrostatic contributions (3G, 4G) this
is not possible, and we handle the virial calculation ourselves.

We currently initialize a global interface object on the library side, so it is
not possible to initialize our pair style multiple times without serious problems.

#### `compute` Routine - General Workflow
We implement a single `compute` routine that manages the calculations for all
implemented generations. The program sequence runs like this:
1. `runner_lammps_interface_transfer_atoms_and_neighbor_lists`: Transfer atoms
  and neighbor list from LAMMPS to the format required by our library.
2. `runner_interface_short_range`: Predict all properties that only depend on
  the local environment of an atom. This includes atomic energies (2G), atomic
  charges (3G/4G), atomic electronegativities (4G), atomic hardnesses (4G),
  atomic Hirshfeld volumes (vdW), and short-range force and virial
  contributions. In the case of a 2G-HDNNP, this is already the end of the
  routine.
3. Based on the predicted atomic properties, we may calculate vdW contributions
  (`runner_interface_hirshfeld_vdw`).
4. If requested, a repulsive two-body potential may be calculated through
  `runner_interface_two_body`.
5. For a 3G-HDNNP:
  5.1 Collect the atomic charges on the root process.
  5.2 Evaluate long-range electrostatics on the root process in
    `runner_interface_evaluate_electrostatics_3g_part_1`
  5.3 Unpack the electrostatic results back to each process.
  5.4 Evaluate additional local electrostatic screening and gradient contributions
    on each process.
6. For a 4G-HDNNP:
  6.1 Collect the atomic electronegativities and hardnesses on the root process.
  6.2 Perform charge equilibration with `runner_interface_compute_charges_4g`.
  6.3 Unpack the globally calculated charges pack to local and ghost atoms on 
     each MPI process.
  6.4 On each process, calculate the short-range (2G) atomic energies with
    `runner_interface_short_range_4g`.
  6.5 Optionally calculate local electrostatic screening.
  6.6 In a second serial step, we calculate any globally-dependent electrostatic
    contributions (`runner_interface_evaluate_electrostatics_4g_part_1`).
  6.7 The global contributions are unpacked again and final gradient contributions
    are added in `runner_interface_evaluate_electrostatics_4g_part_2`.

#### Committee treatment
All implemented methods support the prediction of multiple MLPs with different
parameters. This so-called committee approach is a very intuitive error metric
and helpful for active learning. It is very efficient, as committee members must
always share the same feature vector for each atom.

Many of our library routines exploit this fact. However, some calculation steps
must happen serially with respect to the committee members. These are:

- the Hirshfeld vdW calculation
- the 3G-HDNNP long-range electrostatic calculation
- the 4G-HDNNP QEq step and the long-range electrostatic calculation. Please note
  that the short-range prediction step in between happens for all committee members
  at once again, which is why there are two separate loops over committee members.

#### Computational Scaling and Parallelization
All implemented methods scale linearly with the number of atoms. Depending on the
system, we spend approximately 1% to 2% of serial CPU time transfering atoms
and neighbor lists to our code. The remaining time is mostly spend in the calculation
of the feature vectors (80%) and the forward and backward passes through our
machine learning models (15%).

##### OpenMP
Our library uses OpenMP where possible and useful. We tested the scaling
for systems up to a few thousand atoms and it is near linear in most cases.

##### MPI
2G-HDNNPs scale linearly with the number of MPI tasks because there is no
additional communication necessary (apart from what LAMMPS does).

As is evident from the "General Workflow" section, 3G-HDNNPs and 4G-HDNNPs
require a lot of MPI communication because we need global information for
the long-range electrostatics calculation and charge equilibration steps. This
results in a decreased MPI scalability.

### Post Submission Checklist

<!--Please check the fields below as they are completed after the pull request has been submitted. Delete lines that don't apply-->

- [ ] The feature or features in this pull request is complete
- [x] Licensing information is complete
- [x] Corresponding author information is complete
- [x] The source code follows the LAMMPS formatting guidelines
- [ ] Suitable new documentation files and/or updates to the existing docs are included
- [ ] The added/updated documentation is integrated and tested with the documentation build system
- [x] The feature has been verified to work with the CMake based build system
- [ ] Suitable tests have been added to the unittest tree.
- [ ] A package specific README file has been included or updated
- [ ] One or more example input decks are included

### Further Information, Files, and Links

- The n2p2 documentation: https://compphysvienna.github.io/n2p2/index.html
