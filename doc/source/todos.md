# ToDos                                                                                {#page_todos}

- Transfer
  - think about: bind each transfer to two sweepers (maybe with backreferences)?

- Communicator
  - MPI-speciallization for non-blocking point-to-point and one-sided communication

- check test coverage

- write API documentation

- coding scheme
  - write up guide
    - non-public member variables with leading underscore
    - accessor methods for non-public variables
      - no-logic-setter methods just return reference
      - consistency-checking setter methods prepend with `set_`
      - read-only accessors of the form `const T get_variable() const`
  - cleanup quadrature stuff to comply with new scheme

- Logging
  - reformat logging output (especially DEBUG level)
  
- CMake
  - update GTest/GMock to point to github.com/google/googletest and deal with new layout
    (or switch to Catch completely)
