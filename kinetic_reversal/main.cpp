#include "kinetic_simulation.h"

int main(int argc, char* argv[]) {
  make_dir();
  read_para();
  simulation();
  return 0;
}