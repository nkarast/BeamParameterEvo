#pragma once

#include <gsl/gsl_integration.h>

// Simple RAII wrapper 
class gsl_integration_workspace_pp {
  gsl_integration_workspace * wsp;

  public:
  gsl_integration_workspace_pp(const size_t n=1024): wsp(gsl_integration_workspace_alloc(n)) {}
  ~gsl_integration_workspace_pp() { gsl_integration_workspace_free(wsp); }

  operator gsl_integration_workspace*() { return wsp; }
};

