#ifndef HCLIB_BALE_ACTOR_H
#define HCLIB_BALE_ACTOR_H

#include "hclib-module.h"
#include "hclib_cpp.h"

extern int nic_locale_id;
extern hclib::locale_t *nic;

namespace hclib {

HCLIB_MODULE_INITIALIZATION_FUNC(bale_actor_pre_initialize);
HCLIB_MODULE_INITIALIZATION_FUNC(bale_actor_post_initialize);
HCLIB_MODULE_INITIALIZATION_FUNC(bale_actor_finalize);

int shmem_my_pe();
int shmem_n_pes();

}; //namespace hclib

#endif
