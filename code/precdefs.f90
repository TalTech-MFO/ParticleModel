module mod_precdefs
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  integer, parameter :: rk = selected_real_kind(8)
  integer, parameter :: LEN_CHAR_S = 64
  integer, parameter :: LEN_CHAR_L = 256
end module mod_precdefs
