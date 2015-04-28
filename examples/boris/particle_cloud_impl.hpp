#include "particle_cloud.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <random>
#include <string>
using namespace std;

#include <boost/format.hpp>

#include <pfasst/globals.hpp>
#include <pfasst/site_config.hpp>
#include <pfasst/logging.hpp>
#ifdef WITH_MPI
  #include <mpi.h>
#endif

#include "particle_util.hpp"

#define PFASST_RANDOM_SEED 42


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
#ifdef WITH_MPI
      template<typename precision>
      inline mpi::MPICommunicator& ParticleCloud<precision>::as_mpi(ICommunicator* comm)
      {
        auto mpi = dynamic_cast<mpi::MPICommunicator*>(comm);
        assert(mpi);
        return *mpi;
      }
#endif


      template<typename precision>
      ParticleCloud<precision>::ParticleCloud(const size_t num_particles, const size_t dim,
                                              const precision default_charge, const precision default_mass)
        :   _dim(dim)
          , _num_particles(num_particles)
          , _positions(num_particles * dim)
          , _velocities(num_particles * dim)
          , _charges(num_particles, default_charge)
          , _masses(num_particles, default_mass)
          , _default_charge(default_charge)
          , _default_mass(default_mass)
#ifdef WITH_MPI
          , recv_request(2)
          , send_request(2)
#endif
      {
        this->zero();
      }

      template<typename precision>
      ParticleCloud<precision>::~ParticleCloud()
      {}

      template<typename precision>
      void ParticleCloud<precision>::zero()
      {
        fill(this->_positions.begin(), this->_positions.end(), precision(0.0));
        fill(this->_velocities.begin(), this->_velocities.end(), precision(0.0));
        fill(this->_charges.begin(), this->_charges.end(), this->_default_charge);
        fill(this->_masses.begin(), this->_masses.end(), this->_default_mass);
#ifdef WITH_MPI
        fill(this->recv_request.begin(), this->recv_request.end(), MPI_REQUEST_NULL);
        fill(this->send_request.begin(), this->send_request.end(), MPI_REQUEST_NULL);
#endif
      }

      template<typename precision>
      void ParticleCloud<precision>::copy(shared_ptr<const encap::Encapsulation<precision>> other)
      {
        shared_ptr<const ParticleCloud<precision>> other_c = dynamic_pointer_cast<const ParticleCloud<precision>>(other);
        assert(other_c);
//         this->copy(other_c);
//       }

//       template<typename precision>
//       void ParticleCloud<precision>::copy(shared_ptr<const ParticleCloud<precision>> other)
//       {
        this->_dim = other_c->dim();
        this->_num_particles = other_c->size();
        this->_positions = other_c->positions();
        this->_velocities = other_c->velocities();
        this->_charges = other_c->charges();
        this->_masses = other_c->masses();
        this->_default_charge = other_c->_default_charge;
        this->_default_mass = other_c->_default_mass;
      }

      template<typename precision>
      inline size_t ParticleCloud<precision>::size() const
      {
        return this->_num_particles;
      }

      template<typename precision>
      inline size_t ParticleCloud<precision>::dim() const
      {
        return this->_dim;
      }

      template<typename precision>
      ParticleCloudComponent<precision>& ParticleCloud<precision>::positions()
      {
        return this->_positions;
      }

      template<typename precision>
      const ParticleCloudComponent<precision>& ParticleCloud<precision>::positions() const
      {
        return this->_positions;
      }

      template<typename precision>
      ParticleCloudComponent<precision>& ParticleCloud<precision>::velocities()
      {
        return this->_velocities;
      }

      template<typename precision>
      const ParticleCloudComponent<precision>& ParticleCloud<precision>::velocities() const
      {
        return this->_velocities;
      }

      template<typename precision>
      vector<precision>& ParticleCloud<precision>::charges()
      {
        return this->_charges;
      }

      template<typename precision>
      const vector<precision>& ParticleCloud<precision>::charges() const 
      {
        return this->_charges;
      }

      template<typename precision>
      vector<precision>& ParticleCloud<precision>::masses()
      {
        return this->_masses;
      }

      template<typename precision>
      const vector<precision>& ParticleCloud<precision>::masses() const
      {
        return this->_masses;
      }

      template<typename precision>
      vector<shared_ptr<Particle<precision>>> ParticleCloud<precision>::particles() const
      {
        vector<shared_ptr<Particle<precision>>> particles(this->size());
        for (size_t index = 0; index < this->size(); ++index) {
          particles[index] = this->operator[](index);
        }
        return particles;
      }

      /**
       * @todo take particles of other ranks into account!
       */
      template<typename precision>
      ParticleComponent<precision> ParticleCloud<precision>::center_of_mass() const
      {
        ParticleComponent<precision> center(this->dim());
        for (size_t p = 0; p < this->size(); ++p) {
          // conter += position[p]
          std::transform(center.cbegin(), center.cend(), this->positions().cbegin() + (p * this->dim()),
                         center.begin(), std::plus<precision>());
        }
        return center / this->size();
      }

      template<typename precision>
      shared_ptr<Particle<precision>> ParticleCloud<precision>::operator[](const size_t index) const
      {
        shared_ptr<Particle<precision>> particle = make_shared<Particle<precision>>(this->dim());
        copy_n(this->positions().cbegin() + (index * this->dim()), this->dim(), particle->pos().begin());
        copy_n(this->velocities().cbegin() + (index * this->dim()), this->dim(), particle->vel().begin());
        particle->set_charge(this->_charges[index]);
        particle->set_mass(this->_masses[index]);
        return particle;
      }

      template<typename precision>
      shared_ptr<Particle<precision>> ParticleCloud<precision>::at(const size_t index) const
      {
        assert(this->size() > index);
        return this->operator[](index);
      }

      template<typename precision>
      void ParticleCloud<precision>::set_at(const size_t index, const shared_ptr<Particle<precision>>& particle)
      {
        assert(this->size() > index);
        assert(particle->dim() == this->dim());
        copy_n(particle->pos().cbegin(), this->dim(), this->positions().begin() + (index * this->dim()));
        copy_n(particle->vel().cbegin(), this->dim(), this->velocities().begin() + (index * this->dim()));
        this->_masses[index] = particle->mass();
        this->_charges[index] = particle->charge();
      }

      template<typename precision>
#ifdef WITH_MPI
      void ParticleCloud<precision>::distribute_around_center(const shared_ptr<Particle<precision>>& center,
                                                              pfasst::mpi::MPICommunicator space_comm)
#else
      void ParticleCloud<precision>::distribute_around_center(const shared_ptr<Particle<precision>>& center)
#endif
      {
        CVLOG(3, "Boris") << LOG_INDENT << "distributing " << this->size()
                                        << " particles around center " << center;
        assert(this->size() > 0);

        precision scale = 1000.0;

        default_random_engine rd_gen(PFASST_RANDOM_SEED);
        precision max_pos = max(center->pos());
        precision max_vel = max(center->vel());
        uniform_real_distribution<precision> dist_pos(- max_pos / scale, max_pos / scale);
        uniform_real_distribution<precision> dist_vel(- max_vel / scale, max_vel / scale);
        CVLOG(4, "Boris") << LOG_INDENT << "random displacement range for";
        CVLOG(4, "Boris") << LOG_INDENT << " ... position: "
                                        << boost::format("[%.4f, %.4f]") % dist_pos.min() % dist_pos.max();
        CVLOG(4, "Boris") << LOG_INDENT << " ... velocity: "
                                        << boost::format("[%.4f, %.4f]") % dist_vel.min() % dist_vel.max();

#ifdef WITH_MPI
        // make sure the same random numbers are drawn independent of the degree of spatial parallelization
        if (space_comm.size() > 1) {
          CVLOG(4, "Boris") << LOG_INDENT << "skipping a total of "
                                          << space_comm.rank() * this->size() * this->dim()
                                          << " random numbers each for positions and velocities";
          for (size_t skip = 0; skip < space_comm.rank() * this->size() * this->dim(); ++skip) {
            dist_pos(rd_gen);
            dist_vel(rd_gen);
          }
        }
#endif

        for (size_t p = 0; p < this->size(); ++p) {
          ParticleComponent<precision> pos_rand(this->dim());
          ParticleComponent<precision> vel_rand(this->dim());
          std::transform(center->pos().cbegin(), center->pos().cend(), pos_rand.begin(),
                         [&](const precision& c) { return c + dist_pos(rd_gen); });
          std::transform(center->pos().cbegin(), center->pos().cend(), vel_rand.begin(),
                         [&](const precision& c) { return c + dist_vel(rd_gen); });
          std::copy(pos_rand.cbegin(), pos_rand.cend(), this->positions().begin() + (p * this->dim()));
          std::copy(vel_rand.cbegin(), vel_rand.cend(), this->velocities().begin() + (p * this->dim()));
          CVLOG(5, "Boris") << LOG_INDENT << "p=" << (p+1) << ": " << this->at(p);
        }
        CVLOG(3, "Boris") << LOG_INDENT << "center after distribute: " << this->center_of_mass();
      }

      template<typename precision>
      precision ParticleCloud<precision>::norm0() const
      {
        return std::max(max_abs(this->positions()), max_abs(this->velocities()));
      }

#ifdef WITH_MPI
      template<typename precision>
      void ParticleCloud<precision>::post(ICommunicator* comm, int tag)
      {
        auto& mpi = as_mpi(comm);
        if (mpi.size() == 1) {
          CLOG(DEBUG, "Communicator") << "no posted data to receive as only one rank in communicator '" << mpi.name() << "'";
          return;
        }
        if (mpi.rank() == 0) {
          CLOG(DEBUG, "Communicator") << "first rank will not receive posted data";
          return;
        }

        int src_rank = (mpi.rank() - 1) % mpi.size();
        CVLOG(9, "Communicator") << "receiving posted data with tag " << tag
                                 << " from rank " << src_rank << " of communicator " << mpi.name();
        int err = MPI_Irecv(this->positions().data(),
                            sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                            src_rank, tag, mpi.comm, &(this->recv_request[0]));
        if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
        CVLOG(9, "Communicator") << "received posted positions: " << this->positions();

        err = MPI_Irecv(this->velocities().data(),
                        sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                        src_rank, tag + 1, mpi.comm, &(this->recv_request[1]));
        if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
        CVLOG(9, "Communicator") << "received posted velocities: " << this->velocities();
      }

      template<typename precision>
      void ParticleCloud<precision>::recv(ICommunicator* comm, int tag, bool blocking)
      {
        auto& mpi = as_mpi(comm);
        if (mpi.size() == 1) {
          CLOG(DEBUG, "Communicator") << "nothing to receive as only one rank in communicator '" << mpi.name() << "'";
          return;
        }
        if (mpi.rank() == 0) {
          CLOG(DEBUG, "Communicator") << "first rank will not receive";
          return;
        }

        int err;
        MPI_Status stat;
        if (blocking) {
          int src_rank = (mpi.rank() - 1) % mpi.size();
          CVLOG(9, "Communicator") << "blocked receive with tag " << tag
                                   << " from rank " << src_rank << " of communicator " << mpi.name();
          err = MPI_Recv(this->positions().data(),
                         sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                         src_rank, tag, mpi.comm, &stat);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          CVLOG(9, "Communicator") << "received positions: " << this->positions();

          err = MPI_Recv(this->velocities().data(),
                         sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                         src_rank, tag + 1, mpi.comm, &stat);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          CVLOG(9, "Communicator") << "received velocities: " << this->positions();
        } else {
          CVLOG(9, "Communicator") << "immediate receiving with tag " << tag
                                   << " in communicator " << mpi.name();
          for (auto req : this->recv_request) {
            CVLOG(9, "Communicator") << "waiting for previous request";
            err = MPI_Wait(&req, &stat);
            if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          }
        }
      }

      template<typename precision>
      void ParticleCloud<precision>::send(ICommunicator* comm, int tag, bool blocking)
      {
        auto& mpi = as_mpi(comm);
        if (mpi.size() == 1) {
          CLOG(DEBUG, "Communicator") << "nothing to send as only one rank in communicator '" << mpi.name() << "'";
          return;
        }
        if (mpi.rank() == mpi.size() - 1) {
          CLOG(DEBUG, "Communicator") << "last rank will not send";
          return;
        }

        int dest_rank = (mpi.rank() + 1) % mpi.size();
        int err = MPI_SUCCESS;
        if (blocking) {
          CVLOG(9, "Communicator") << "blocked sending with tag " << tag
                                   << " to rank " << dest_rank << " of communicator " << mpi.name();
          err = MPI_Send(this->positions().data(),
                         sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                         dest_rank, tag, mpi.comm);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          CVLOG(9, "Communicator") << "sent positions: " << this->positions();

          err = MPI_Send(this->velocities().data(),
                         sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                         dest_rank, tag + 1, mpi.comm);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          CVLOG(9, "Communicator") << "sent velocities: " << this->velocities();
        } else {
          CVLOG(9, "Communicator") << "immediate sending with tag " << tag
                                   << " to rank " << dest_rank << " of communicator " << mpi.name();
          MPI_Status stat;
          for (auto req : this->send_request) {
            CVLOG(9, "Communicator") << "waiting for previous request";
            err = MPI_Wait(&req, &stat);
            if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          }

          err = MPI_Isend(this->positions().data(),
                          sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                          dest_rank, tag, mpi.comm, &send_request[0]);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          CVLOG(9, "Communicator") << "sent positions: " << this->positions();

          err = MPI_Isend(this->velocities().data(),
                          sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                          dest_rank, tag + 1, mpi.comm, &send_request[1]);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          CVLOG(9, "Communicator") << "sent velocities: " << this->velocities();
        }
      }

      template<typename precision>
      void ParticleCloud<precision>::broadcast(ICommunicator* comm)
      {
        auto& mpi = as_mpi(comm);

        int root_rank = comm->size() - 1;
        CVLOG(9, "Communicator") << "broadcasting with root " << root_rank << " of '" << mpi.name() << "'";
        int err = MPI_Bcast(this->positions().data(),
                            sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                            root_rank, mpi.comm);
        if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
        CVLOG(9, "Communicator") << "broadcasted positions: " << this->positions();

        err = MPI_Bcast(this->velocities().data(),
                        sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                        root_rank, mpi.comm);
        if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
        CVLOG(9, "Communicator") << "broadcasted velocities: " << this->velocities();
      }
#endif

      template<typename precision>
      void ParticleCloud<precision>::log(el::base::type::ostream_t& os) const
      {
        os << fixed << setprecision(LOG_PRECISION);
        os << "ParticleCloud(n=" << this->size() << ", d=" << this->dim() << ", particles=" << this->particles() << ")";
        os.unsetf(ios_base::floatfield);
      }


      template<typename precision>
      static precision distance(const Particle<precision>& first, const Particle<precision>& second)
      {
        assert(first.dim() == second.dim());
        vector<precision> diff = first.pos() - second.pos();
        return norm0(diff);
      }

      template<typename precision>
      static precision distance(const shared_ptr<Particle<precision>> first, const shared_ptr<Particle<precision>> second)
      {
        return distance(*(first.get()), *(second.get()));
      }

      template<typename precision>
      static vector<precision> distance_to_reference(const ParticleCloud<precision>& cloud,
                                                     const Particle<precision>&      reference)
      {
        vector<precision> distances(cloud.size());
        for (size_t i = 0; i < distances.size(); ++i) {
          distances[i] = distance(*cloud[i], reference);
        }
        return distances;
      }

      template<typename precision>
      static vector<precision> distance_to_reference(const shared_ptr<ParticleCloud<precision>>& cloud,
                                                     const shared_ptr<Particle<precision>>&      reference)
      {
        return distance_to_reference(*cloud, *reference);
      }


      template<typename precision>
      inline MAKE_LOGGABLE(shared_ptr<ParticleCloud<precision>>, sp_cloud, os)
      {
        os << "<" << addressof(sp_cloud) << ">";
        sp_cloud->log(os);
        return os;
      }

      template<typename precision>
      inline MAKE_LOGGABLE(shared_ptr<const ParticleCloud<precision>>, sp_cloud, os)
      {
        os << "<" << addressof(sp_cloud) << ">";
        sp_cloud->log(os);
        return os;
      }


      template<typename precision>
      ParticleCloudFactory<precision>::ParticleCloudFactory(const size_t num_particles, const size_t dim,
                                                            const precision default_charge,
                                                            const precision default_mass)
        :   _num_particles(num_particles)
          , _dim(dim)
          , _default_charge(default_charge)
          , _default_mass(default_mass)
      {}

      template<typename precision>
      inline size_t ParticleCloudFactory<precision>::num_particles() const
      {
        return this->_num_particles;
      }

      template<typename precision>
      inline size_t ParticleCloudFactory<precision>::dim() const
      {
        return this->_dim;
      }

      template<typename precision>
      shared_ptr<encap::Encapsulation<precision>>
      ParticleCloudFactory<precision>::create(const encap::EncapType)
      {
        return make_shared<ParticleCloud<precision>>(this->num_particles(), this->dim(),
                                                     this->_default_charge, this->_default_mass);
      }
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst
