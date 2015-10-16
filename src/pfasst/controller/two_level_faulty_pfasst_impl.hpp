#include "pfasst/controller/two_level_faulty_pfasst.hpp"

#include <cassert>
using namespace std;

#include "pfasst/util.hpp"
#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/config.hpp"


namespace pfasst
{
  template<class TransferT, class CommT>
  void
  TwoLevelFaultyPfasst<TransferT, CommT>::init_loggers()
  {
    TwoLevelPfasst<TransferT, CommT>::init_loggers();
    log::add_custom_logger("FT_PFASST");
  }

  template<class TransferT, class CommT>
  void
  TwoLevelFaultyPfasst<TransferT, CommT>::init_opts()
  {
    config::options::add_option<string>("Faulty PFASST", "reset",
                                        "comma separated list of '<rank>-<step>-<iter>' tuples specifying sweeper resets");
  }

  template<class TransferT, class CommT>
  TwoLevelFaultyPfasst<TransferT, CommT>::TwoLevelFaultyPfasst()
    : TwoLevelPfasst<TransferT, CommT>()
  {
    TwoLevelFaultyPfasst<TransferT, CommT>::init_loggers();
    this->set_logger_id("FT_PFASST");
  }

  template<class TransferT, class CommT>
  void
  TwoLevelFaultyPfasst<TransferT, CommT>::set_options()
  {
    TwoLevelPfasst<TransferT, CommT>::set_options();

    if (config::has_value("reset")) {
      auto resets = split(config::get_value<string>("reset"), ',');

      // loop over reset points
      for (auto r : resets) {
        auto p = split(r, '-');
        if (p.size() != 3) {
          ML_CLOG(FATAL, this->get_logger_id(),
                  "The faulty point is not specified as a dash-delimited '<rank>-<step>-<iter>' tuple: " << r);
          throw runtime_error("faulty-points must be specified as dash-delimited '<rank>-<step>-<iter>' tuples");
        }

        size_t rank = strtoul(p[0].c_str(), nullptr, 10);
        if (rank >= this->get_communicator()->get_size()) {
          ML_CLOG(FATAL, this->get_logger_id(),
                  "Specified faulty rank is out of range [0,"
                  << this->get_communicator()->get_size() << ") :" << rank);
          throw runtime_error("invalid rank specified");
        }

        size_t step = strtoul(p[1].c_str(), nullptr, 10);
        size_t iter = strtoul(p[2].c_str(), nullptr, 10);

        if (this->_faults.find(rank) == this->_faults.end()) {
          map<size_t, size_t> step_iter = {{step, iter}};
          this->_faults[rank] = step_iter;
        } else {
          if (this->_faults[rank].find(step) == this->_faults[rank].end()) {
            this->_faults[rank][step] = iter;
          } else {
            ML_CLOG(WARNING, this->get_logger_id(),
                    "Only one reset per rank and time step possible.");
          }
        }

        ML_CLOG(INFO, this->get_logger_id(),
                "will simulate a fault at rank " << rank << " on step " << step << " in iteration " << iter << ".");
      }
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelFaultyPfasst<TransferT, CommT>::setup()
  {
    TwoLevelPfasst<TransferT, CommT>::setup();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelFaultyPfasst<TransferT, CommT>::run()
  {
    Controller<TransferT, CommT>::run();

    assert(this->get_communicator() != nullptr);

    if (this->get_status()->get_num_steps() % this->get_communicator()->get_size() != 0) {
      ML_CLOG(ERROR, this->get_logger_id(), "Number of time steps (" << this->get_status()->get_num_steps()
                                         << ") must be a multiple of the number of processors ("
                                         << this->get_communicator()->get_size() << ").");
      throw logic_error("number time steps must be multiple of number processors");
    }

    const size_t num_blocks = this->get_status()->get_num_steps() / this->get_communicator()->get_size();

    if (num_blocks == 0) {
      ML_CLOG(ERROR, this->get_logger_id(), "Invalid Duration: There are more time processes ("
                                         << this->get_communicator()->get_size() << ") than time steps ("
                                         << this->get_status()->get_num_steps() << ").");
      throw logic_error("invalid duration: too many time processes for given time steps");
    }

    // iterate over time blocks (i.e. time-parallel blocks)
    do {
      this->status()->step() = this->_time_block * this->get_communicator()->get_size() \
                               + this->get_communicator()->get_rank();
      if (this->_time_block == 0) {
        this->status()->time() += this->get_status()->get_dt() * this->status()->get_step();
      }

      ML_CLOG(INFO, this->get_logger_id(), "");
      ML_CLOG(INFO, this->get_logger_id(), "Time Step " << (this->get_status()->get_step() + 1)
                                                        << " of " << this->get_status()->get_num_steps()
                                                        << " (i.e. t0=" << this->get_status()->get_time() << ")");

      // XXX: required here?
      this->_prev_status->clear();
      this->_prev_status_temp->clear();

      this->status()->set_primary_state(PrimaryState::PREDICTING);

      // iterate on each time step (i.e. iterations on single time step)
      do {
        if (this->get_status()->get_primary_state() == (+PrimaryState::PREDICTING)) {
          this->predictor();

        } else if (this->get_status()->get_primary_state() == (+PrimaryState::ITERATING)) {
          ML_CLOG(INFO, this->get_logger_id(), "");
          ML_CLOG(INFO, this->get_logger_id(), "Iteration " << this->get_status()->get_iteration());

          this->cycle_down();

          this->recv_coarse();
          this->sweep_coarse();
          this->send_coarse();

          this->cycle_up();

          if (be_faulty()) {
            // backup step and iteration index
            const auto step = this->get_status()->get_step();
            const auto iter = this->get_status()->get_iteration();
            this->get_coarse()->reset();
            this->get_fine()->reset();
            this->status()->step() = step;
            this->status()->iteration() = iter;
            this->status()->set_primary_state(PrimaryState::PREDICTING);
          }

          this->sweep_fine();
          this->send_fine();

        } else {
          ML_CLOG(FATAL, this->get_logger_id(), "Something went severly wrong with the states.");
          ML_CLOG(FATAL, this->get_logger_id(), "Expected state: PREDICTING or ITERATING, got: "
                                                << (+this->get_status()->get_primary_state())._to_string());
          throw runtime_error("something went severly wrong");
        }

        // convergence check
      } while(this->advance_iteration());

      ML_CLOG(INFO, this->get_logger_id(), "");
      ML_CLOG(INFO, this->get_logger_id(), "Time Step done.");
    } while(this->advance_time(this->get_communicator()->get_size()));
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelFaultyPfasst<TransferT, CommT>::advance_iteration()
  {
    this->status()->set_primary_state(PrimaryState::INTER_ITER);

    this->get_check_prev_status();

    const bool fine_converged = this->get_fine()->converged(true);
    const bool previous_done = (this->get_communicator()->is_first())
                               ? true
                               : this->_prev_status->get_primary_state() <= (+PrimaryState::FAILED);
    ML_CLOG(DEBUG, this->get_logger_id(), "this status: " << to_string(this->status()));
    ML_CLOG(DEBUG, this->get_logger_id(), "prev status: " << to_string(this->_prev_status));

    if (previous_done && fine_converged) {
      ML_CLOG(INFO, this->get_logger_id(), "FINE sweeper has converged as well as previous process.");

      // receive potentially pending fine data of previous process
      this->recv_fine();

      this->status()->set_primary_state(PrimaryState::CONVERGED);
      this->_prev_status->set_primary_state(PrimaryState::UNKNOWN_PRIMARY);

    } else {
      ML_CLOG_IF(previous_done && !fine_converged, INFO, this->get_logger_id(),
        "previous process has converged but FINE sweeper not yet.");

      if (Controller<TransferT, CommT>::advance_iteration()) {
        ML_CLOG(INFO, this->get_logger_id(), "FINE sweeper has not yet converged and additional iterations to do.");
        this->get_fine()->save();
        this->get_coarse()->save();
        this->status()->set_primary_state(PrimaryState::ITERATING);

      } else {
        ML_CLOG(WARNING, this->get_logger_id(), "FINE sweeper has not yet converged and iterations threshold reached.");

        // receive potentially pending fine data of previous process
        this->recv_fine(true);

        this->status()->set_primary_state(PrimaryState::FAILED);
      }
    }

    this->send_status();

    this->get_fine()->converged(false);

    return (this->get_status()->get_primary_state() > (+PrimaryState::FAILED));
  }

  template<class TransferT, class CommT>
  void
  TwoLevelFaultyPfasst<TransferT, CommT>::recv_coarse()
  {
    if (!this->get_communicator()->is_first()) {
      if (this->_prev_status->get_primary_state() > (+PrimaryState::FAILED)) {
        ML_CVLOG(2, this->get_logger_id(), "looking for coarse data");
        this->get_coarse()
            ->initial_state()
            ->recv(this->get_communicator(),
                   this->get_communicator()->get_rank() - 1,
                   this->compute_tag(TagType::DATA, TagLevel::COARSE, TagModifier::PREV_STEP),
                   true);
      } else {
        ML_CLOG(WARNING, this->get_logger_id(), "previous process doesn't send any coarse data any more");
      }
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelFaultyPfasst<TransferT, CommT>::get_check_prev_status()
  {
    if (!this->get_communicator()->is_first()) {
      if (this->_prev_status->get_primary_state() > (+PrimaryState::FAILED)) {
        ML_CLOG(DEBUG, this->get_logger_id(), "prev known status: " << to_string(this->_prev_status));
        this->_prev_status_temp->clear();

        ML_CVLOG(1, this->get_logger_id(), "looking for updated state of previous process");
        this->_prev_status_temp->recv(this->get_communicator(),
                                      this->get_communicator()->get_rank() - 1,
                                      this->compute_tag(TagType::STATUS, TagLevel::FINE, TagModifier::PREV_STEP), true);
        // copy latest received status to the place where we use it from
        *(this->_prev_status) = *(this->_prev_status_temp);
        ML_CLOG(DEBUG, this->get_logger_id(), "Status received: " << to_string(this->_prev_status));

        if (this->_prev_status->get_primary_state() == (+PrimaryState::FAILED)) {
          ML_CLOG(WARNING, this->get_logger_id(), "previous process failed");

        } else if (this->_prev_status->get_primary_state() == (+PrimaryState::CONVERGED)) {
          ML_CLOG(WARNING, this->get_logger_id(), "previous process has converged; this process not");

        } else {
          ML_CVLOG(1, this->get_logger_id(), "previous process not finished");
        }
      }
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelFaultyPfasst<TransferT, CommT>::predictor()
  {
    assert(this->get_status()->get_iteration() == 0);

    ML_CLOG(INFO, this->get_logger_id(), "");
    ML_CLOG(INFO, this->get_logger_id(), "PFASST Prediction step");

    // restrict fine initial condition ...
    this->get_transfer()->restrict_initial(this->get_fine(), this->get_coarse());
    // ... and spread it to all nodes on the coarse level
    this->get_coarse()->spread();
    this->get_coarse()->save();

    // perform PFASST prediction sweeps on coarse level
    for (size_t predict_step = 0;
         predict_step <= this->get_communicator()->get_rank();
         ++predict_step) {
      // do the sweeper's prediction once ...
      if (predict_step == 0) {
        this->predict_coarse();
      } else {
        // and default sweeps for subsequent processes
        this->recv_coarse();
        this->sweep_coarse();
      }

      this->send_coarse();
    }

    // return to fine level
    ML_CVLOG(1, this->get_logger_id(), "cycle up onto fine level");
    this->get_transfer()->interpolate(this->get_coarse(), this->get_fine(), true);
    this->sweep_fine();

    this->send_fine();

    // finalize prediction step
    this->get_coarse()->save();
    this->get_fine()->save();
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelFaultyPfasst<TransferT, CommT>::be_faulty()
  {
    const size_t rank = this->get_communicator()->get_rank();
    const size_t step = this->get_status()->get_step();
    const size_t iter = this->get_status()->get_iteration();

    if (this->_faults.find(rank) != this->_faults.end()) {
      if (this->_faults[rank].find(step) != this->_faults[rank].end()) {
        if (this->_faults[rank][step] == iter) {
          ML_CLOG(WARNING, this->get_logger_id(), "I'AM FAULTY NOW!");
          this->_faults[rank].erase(step);
          return true;
        }
      }
    }
    ML_CLOG(DEBUG, this->get_logger_id(), "no fault on step " << step << " in iteration " << iter);
    return false;
  }

  template<class TransferT, class CommT>
  int
  TwoLevelFaultyPfasst<TransferT, CommT>::compute_tag(const TagType type, const TagLevel level,
                                                      const TagModifier mod) const
  {
    int tag = (type == (+TagType::DATA)) ? 1 : 0;

    const size_t step = this->get_status()->get_step()
                        - ((   mod == (+TagModifier::PREV_STEP)
                            || mod == (+TagModifier::PREV_ITER_PREV_STEP)) ? 1 : 0);
    tag += (step + 1) * 1000000;

    tag += ((+level)._to_integral() + 1) * 100;

    ML_CLOG(DEBUG, this->get_logger_id(),
            "computing tag for " << (+type)._to_string() << " communication "
            << "on " << (+level)._to_string() << " level "
            << (mod == (+TagModifier::UNMOD) ? "without modifier" : "with modifier ")
            << (mod != (+TagModifier::UNMOD) ? (+mod)._to_string() : "")
            << " --> " << tag);

    return tag;
  }
}  // ::pfasst
