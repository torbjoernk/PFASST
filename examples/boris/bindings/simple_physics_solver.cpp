#include "simple_physics_solver.hpp"

#include <cassert>
#include <iostream>
#include <cmath>
using namespace std;

#define UNUSED(expr) (void)(expr)


namespace simple_physics_solver
{
  SimplePhysicsSolverConfig::SimplePhysicsSolverConfig(const double omega_e, const double omega_b,
                                                       const double epsilon, const double sigma)
    :   omega_e(omega_e)
      , omega_b(omega_b)
      , epsilon(epsilon)
      , sigma(sigma)
      , sigma2(sigma * sigma)
  {
    this->external_e_field_matrix[0][0] = double(1.0);
    this->external_e_field_matrix[0][1] = double(0.0);
    this->external_e_field_matrix[0][2] = double(0.0);
    this->external_e_field_matrix[1][0] = double(0.0);
    this->external_e_field_matrix[1][1] = double(1.0);
    this->external_e_field_matrix[1][2] = double(0.0);
    this->external_e_field_matrix[2][0] = double(0.0);
    this->external_e_field_matrix[2][1] = double(0.0);
    this->external_e_field_matrix[2][2] = double(-2.0);

    this->b_field_matrix[0][0] = double(0.0);
    this->b_field_matrix[0][1] = double(1.0);
    this->b_field_matrix[0][2] = double(0.0);
    this->b_field_matrix[1][0] = double(-1.0);
    this->b_field_matrix[1][1] = double(0.0);
    this->b_field_matrix[1][2] = double(0.0);
    this->b_field_matrix[2][0] = double(0.0);
    this->b_field_matrix[2][1] = double(0.0);
    this->b_field_matrix[2][2] = double(0.0);
  }

  SimplePhysicsSolverConfig::~SimplePhysicsSolverConfig()
  {}


  void evaluate_external_e_field(const double* positions, const double* charges, const double* masses,
                                 const size_t num_particles, const double t,
                                 const SimplePhysicsSolverConfig* config,
                                 double* forces)
  {
    UNUSED(t);
    double pre_factor = (- config->epsilon) * (config->omega_e * config->omega_e);
    double factor = double(0.0);
    for (size_t i = 0; i < num_particles; ++i) {
      factor = pre_factor / (charges[i] / masses[i]);
      internal::scale_mat_mul_vec(config->external_e_field_matrix, positions + (i * DIM), factor, forces + (i * DIM));
    }
  }

  void evaluate_internal_e_field(const double* positions, const double* charges, const double* masses,
                                 const size_t num_particles, const double t,
                                 const SimplePhysicsSolverConfig* config,
                                 double* exyz, double* phis)
  {
    UNUSED(masses); UNUSED(t);
    double r = double(0.0),
           r3 = double(0.0),
           dist2 = double(0.0);
    double* dist = new double[DIM];

    // copy this ranks particles; these copies will get reused as send-recv buffers later on
    double* other_positions = new double[num_particles * DIM];
    double* other_charges = new double[num_particles];
    double* other_masses = new double[num_particles];
    for (size_t i = 0; i < num_particles; ++i) {
//       cout << "[" << config->mpi_rank << "] p[" << i << "] = ";
//       internal::print_vec(positions + i * DIM);
//       cout << endl;
      // null result values
      phis[i] = double(0.0);

      for (size_t d = 0 ; d < DIM; ++d) {
        other_positions[i * DIM + d] = positions[i * DIM + d];
        exyz[i * DIM + d] = double(0.0);
      }
      other_charges[i] = charges[i];
      other_masses[i] = masses[i];
    }

    // rank_id where the other particles origin from
    int current_paired = config->mpi_rank;
    // rank_id where other particles in the next round will origin from
    int next_recv = (config->mpi_rank == 0) ? config->mpi_size - 1 : config->mpi_rank - 1;
    // rank_id where other particles will be send to after the loop
    int next_send = (config->mpi_rank == config->mpi_size - 1) ? 0 : config->mpi_rank + 1;

    const int prev_rank = next_recv;
    const int next_rank = next_send;

//     cout << "[" << config->mpi_rank << "] start hot-potato" << endl;
    // The Hot-Potato-Ring-Parallelization-Loop
    for (size_t proc = 0; proc < config->mpi_size; ++proc) {
      // computing forces on particle i
      for (size_t i = 0; i < num_particles; ++i) {
        // effects of particle j on particle i
        for (size_t j = 0; j < num_particles; ++j) {
          // skip pairing the very same particle
          if (proc == 0 && i == j) {
            continue;
          }
//           cout << "[" << config->mpi_rank << "] pairing " << i << " with " << j << endl;

          // compute distance between the two particles
          dist2 = double(0.0);
          for (size_t d = 0; d < DIM; ++d) {
            dist[d] = positions[i * DIM + d] - other_positions[j * DIM + d];
            dist2 += dist[d] * dist[d];
          }
//           cout << "[" << config->mpi_rank << "] dist(p[" << i << "], p[" << j << "]) = ";
//           internal::print_vec(dist);
//           cout << endl;

          // compute the phi
          r = sqrt(dist2 + config->sigma2);
          phis[i] += other_charges[j] / r;

          // compute the Coulomb force
          r3 = r * r * r;
          for (size_t d = 0; d < DIM; ++d) {
            exyz[i * DIM + d] += dist[d] / r3 * other_charges[j];
          }
//           cout << "[" << config->mpi_rank << "] updated force for p[" << i << "] = ";
//           internal::print_vec(exyz + i * DIM);
//           cout << endl;
        }
      }
#ifdef WITH_MPI
      // send-recv
//       cout << "[" << config->mpi_rank << "] exchaning data" << endl;
      MPI_Status stat_pos, stat_charge, stat_mass;
      // NOTE: this will probably only work for spacial parallelization up to 10,000 ranks
      MPI_Sendrecv_replace(other_positions, num_particles * DIM, MPI_DOUBLE,
                           next_rank, 10000 + current_paired,
                           prev_rank, 10000 + next_recv,
                           config->space_comm, &stat_pos);
//       cout << "[" << config->mpi_rank << "] got positions: ";
//       for (size_t p = 0; p < num_particles; ++p) {
//         internal::print_vec(other_positions + p * DIM);
//         cout << "\t";
//       }
//       cout << endl;
      MPI_Sendrecv_replace(other_charges, num_particles, MPI_DOUBLE,
                           next_rank, 20000 + current_paired,
                           prev_rank, 20000 + next_recv,
                           config->space_comm, &stat_charge);
      MPI_Sendrecv_replace(other_masses, num_particles, MPI_DOUBLE,
                           next_rank, 30000 + current_paired,
                           prev_rank, 30000 + next_recv,
                           config->space_comm, &stat_mass);

      // update loop
      current_paired = next_recv;
      next_recv = (next_recv == 0) ? config->mpi_size - 1 : next_recv - 1;
      next_send = (next_send == config->mpi_size - 1) ? 0 : next_send + 1;
#endif
    }
//     cout << "[" << config->mpi_rank << "] done" << endl;

    delete[] other_positions;
    delete[] other_charges;
    delete[] other_masses;
    delete[] dist;
  }


  void evaluate_e_field(const double* positions, const double* charges, const double* masses,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config,
                        double* forces)
  {
    double* external = new double[num_particles * DIM];
    double* internal = new double[num_particles * DIM];
    double* phis = new double[num_particles];

    evaluate_external_e_field(positions, charges, masses, num_particles, t, config, external);
    evaluate_internal_e_field(positions, charges, masses, num_particles, t, config, internal, phis);
    for (size_t i = 0; i < num_particles; ++i) {
      for (size_t d = 0; d < DIM; ++d) {
        forces[i * DIM + d] = external[i * DIM + d] + internal[i * DIM + d];
      }
    }

    delete[] external;
    delete[] internal;
    delete[] phis;
  }


  void get_b_field_vector(const SimplePhysicsSolverConfig* config, double* b_field_vector)
  {
    for (size_t i = 0; i < DIM; ++i) {
      b_field_vector[i] = double(0.0);
    }
    b_field_vector[2] = double(1.0);
    for (size_t i = 0; i < DIM; ++i) {
      b_field_vector[i] *= config->omega_b;
    }
  }

  void evaluate_b_field(const double* velocities, const double* charges, const double* masses,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config, double* forces)
  {
    UNUSED(t);
    for (size_t i = 0; i < num_particles; ++i) {
      internal::scale_mat_mul_vec(config->b_field_matrix, velocities + (i * DIM),
                                  config->omega_b / (charges[i] / masses[i]),
                                  forces + (i * DIM));
    }
  }


  double compute_energy(const double* positions, const double* velocities,
                        const double* charges, const double* masses,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config)
  {
    double e_kin = double(0.0);
    double e_pot = double(0.0);

    double* exyz = new double[num_particles * DIM];
    double* phis = new double[num_particles];
    double* temp = new double[DIM];
    double v2 = double(0.0);

    evaluate_internal_e_field(positions, charges, masses, num_particles, t, config, exyz, phis);

    for (size_t i = 0; i < num_particles; ++i) {
      // potential energy (induced by external electric field, position and internal electric field (i.e. phis))
      internal::scale_mat_mul_vec(config->external_e_field_matrix, positions + (i * DIM),
                                  (- config->epsilon * config->omega_e * config->omega_e / double(2.0) * (charges[i] / masses[i])),
                                  temp);
      e_pot += (charges[i] * phis[i]) - internal::scalar_prod(temp, positions + (i * DIM));

      // kinetic energy (induced by velocity)
      v2 = internal::scalar_prod(velocities + (i * DIM), velocities + (i * DIM));
      e_kin += masses[i] / double(2.0) * v2;
    }

    delete[] exyz;
    delete[] phis;
    delete[] temp;
    return e_kin + e_pot;
  }

  namespace internal
  {
    inline void cross_prod(const double first[DIM], const double second[DIM], double cross_prod[DIM])
    {
      cross_prod[0] = first[1] * second[2] - first[2] * second[1];
      cross_prod[1] = first[2] * second[0] - first[0] * second[2];
      cross_prod[2] = first[0] * second[1] - first[1] * second[0];
    }

    inline double scalar_prod(const double first[DIM], const double second[DIM])
    {
      double prod = 0.0;
      for (size_t m = 0; m < DIM; ++m) {
        prod += first[m] * second[m];
      }
      return prod;
    }

    inline void scale_mat_mul_vec(const double mat[DIM][DIM], const double vec[DIM],
                                  const double factor, double prod[DIM])
    {
      for (size_t i = 0; i < DIM; ++i) {
        prod[i] = double(0.0);
        for (size_t j = 0; j < DIM; ++j) {
          prod[i] += factor * mat[i][j] * vec[j];
        }
      }
    }

    inline void print_vec(const double vec[DIM])
    {
      cout << "[";
      for (size_t i = 0; i < DIM; ++i) {
        cout << vec[i];
        if (i < DIM - 1) {
          cout << " , ";
        }
      }
      cout << "]";
    }
  }  // ::simple_physics_solver::internal
}  // ::simple_physics_solver
