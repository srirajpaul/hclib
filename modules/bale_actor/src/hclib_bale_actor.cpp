
#include <shmem.h>
#include "hclib_bale_actor.h"
#include "hclib-locality-graph.h"

#ifdef USE_OFFLOAD
    #define START_IS_OFFLOAD \
        hclib::finish([=] () {\
        hclib::async_nb_at([=] () {
#else
    #define START_IS_OFFLOAD
#endif

#ifdef USE_OFFLOAD
    #define END_IS_OFFLOAD \
        }, nic);\
        });
#else
    #define END_IS_OFFLOAD
#endif

int nic_locale_id = -1;
hclib::locale_t *nic = NULL;

HCLIB_MODULE_INITIALIZATION_FUNC(bale_actor_pre_initialize) {
#ifdef USE_LOCK
    nic_locale_id = hclib_add_known_locale_type("Interconnect");
    HASSERT(nic_locale_id > -1);
#endif
}

int shmem_init_thread(int, int*);

HCLIB_MODULE_INITIALIZATION_FUNC(bale_actor_post_initialize) {

#ifdef USE_LOCK
#ifdef USE_CRAY_SHMEM_7
    //printf("USE_CRAY_SHMEM_7\n");
    int ret = ::shmem_init_thread(SHMEM_THREAD_MULTIPLE);
    assert(ret == SHMEM_THREAD_MULTIPLE);
#else // USE_CRAY_SHMEM_7
    //printf("NO USE_CRAY_SHMEM_7\n");
    int major, minor;
    shmem_info_get_version(&major, &minor);
    if(major>=1 && minor>=4) {
        int provided;
        int ret = ::shmem_init_thread(SHMEM_THREAD_MULTIPLE, &provided);
        assert(ret == 0);
        assert(provided == SHMEM_THREAD_MULTIPLE);
    }
    else {
        printf("WARNING: SHMEM 1.4 or above required\n");
        ::shmem_init();
    }
#endif // USE_CRAY_SHMEM_7
    int n_nics;
    hclib::locale_t **nics = hclib::get_all_locales_of_type(nic_locale_id,
            &n_nics);
    HASSERT(n_nics == 1);
    HASSERT(nics);
    HASSERT(nic == NULL);
    nic = nics[0];

    hclib_locale_mark_special(nic, "COMM");

#else // USE_LOCK
    //printf("No USE_LOCK\n");
    ::shmem_init();
#endif // USE_LOCK
}

HCLIB_MODULE_INITIALIZATION_FUNC(bale_actor_finalize) {
    ::shmem_finalize();
}

int hclib::shmem_my_pe() {
    return ::shmem_my_pe();
}

int hclib::shmem_n_pes() {
    return ::shmem_n_pes();
}

HCLIB_REGISTER_MODULE("bale_actor", bale_actor_pre_initialize, bale_actor_post_initialize, bale_actor_finalize)

