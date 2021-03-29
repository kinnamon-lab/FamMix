#include <TMB.hpp>
#include "sing_asc_lmm.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model);
  if(model == "sing_asc_lmm") {
    return sing_asc_lmm(this);
  } else {
    error("Unknown model");
  }

  return Type(0);
}
