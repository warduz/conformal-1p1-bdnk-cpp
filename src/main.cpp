#include "Evolver.h"

int main(int argc, char** argv) {
  if (argc > 2) {
    throw std::runtime_error("Invalid number of command line arguments.");
  }

  Evolver evolver { argv[1] };
}