#include "TMath.h"

double ResolutionFull(double *x_val, double *par)
{
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  double resFull = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return resFull;
}
