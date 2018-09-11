// This simulation uses third version of numerical recipe.
// In the future, the solution can be replaced by spectro-method.
#include "kinetic_simulation.h"

int main(int argc, char* argv[]) {
  make_dir();
  read_para();
  simulation();
  return 0;
}