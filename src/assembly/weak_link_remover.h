//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_WEAK_LINK_REMOVER_H
#define MEGAHIT_WEAK_LINK_REMOVER_H

#include <cstdint>
#include "mpienv/mpienv.hpp"
class UnitigGraph;

uint32_t DisconnectWeakLinks(UnitigGraph&, MPIEnviroment &mpienv, double local_ratio);

#endif  // MEGAHIT_WEAK_LINK_REMOVER_H
