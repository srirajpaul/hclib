/*
 * hcpp-asyncStruct.h
 *  
 *      Author: Vivek Kumar (vivekk@rice.edu)
 */

#ifndef HCPP_ASYNCSTRUCT_H_
#define HCPP_ASYNCSTRUCT_H_

#include <string.h>

#include "hcpp-place.h"
#include "hcpp-task.h"

#ifdef __cplusplus
extern "C" {
#endif

inline hclib_ddf_t ** get_ddf_list(hclib_task_t *t) {
    return t->ddf_list;
}

inline void mark_as_asyncAnyTask(hclib_task_t *t) {
    t->is_asyncAnyType = 1;
}

inline int is_asyncAnyTask(hclib_task_t *t) {
    return t->is_asyncAnyType;
}

void spawn(hclib_task_t * task);
void spawn_at_hpt(place_t* pl, hclib_task_t * task);
void spawn_await_at(hclib_task_t * task, hclib_ddf_t** ddf_list, place_t *pl);
void spawn_await(hclib_task_t * task, hclib_ddf_t** ddf_list);
void spawn_commTask(hclib_task_t * task);
void spawn_gpu_task(hclib_task_t *task);

#ifdef __cplusplus
}
#endif

#endif /* HCPP_ASYNCSTRUCT_H_ */
