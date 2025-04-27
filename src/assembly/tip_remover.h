//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_TIP_REMOVER_H
#define MEGAHIT_TIP_REMOVER_H

#include <cstdint>
#include "mpienv/mpienv.hpp"

class UnitigGraph;

uint32_t RemoveTips(UnitigGraph &graph, uint32_t max_tip_len, MPIEnviroment &mpienv);

#endif  // MEGAHIT_TIP_REMOVER_H
