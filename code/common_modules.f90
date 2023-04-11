#include "cppdefs.h"
module mod_common
  !----------------------------------------------------------------
  ! Common modules used by all other modules
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use utils
#ifdef USE_OMP
  use omp_lib
#endif
end module mod_common
