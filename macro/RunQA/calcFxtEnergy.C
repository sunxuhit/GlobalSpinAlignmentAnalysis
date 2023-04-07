#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "../../StRoot/StAnalysisUtils/StAnalysisCons.h"

using namespace std;

void calcFxtEnergy(double beamEnergy = 3.85)
{
  // beam travels -z
  double massProton = anaUtils::mMassProton;
  // int numProtons  = 79;
  // int numNucleons = 197;
  int numNucleons = 1;

  // beam
  double massBeam = massProton*numNucleons;
  double eLabBeam = beamEnergy*numNucleons;
  double pLabBeam = -1.0*TMath::Sqrt(eLabBeam*eLabBeam - massBeam*massBeam); // p^2 = E^2 - m^2
  double bLabBeam = pLabBeam/eLabBeam;
  double yLabBeam = 0.5*TMath::Log((eLabBeam+pLabBeam)/(eLabBeam-pLabBeam));
  cout << "eLabBeam = " << eLabBeam << "GeV, pLabBeam = " << pLabBeam << "GeV, bLabBeam = " << bLabBeam << ", yLabBeam = " << yLabBeam << endl;

  // target
  double massTarget = massProton*numNucleons;
  double eLabTarget = massTarget;
  double pLabTarget = TMath::Sqrt(eLabTarget*eLabTarget - massTarget*massTarget); // p^2 = E^2 - m^2 & in opposite direction as beam
  double bLabTarget = pLabTarget/eLabTarget;
  double yLabTarget = 0.5*TMath::Log((eLabTarget+pLabTarget)/(eLabTarget-pLabTarget));
  cout << "eLabTarget = " << eLabTarget << "GeV, pLabTarget = " << pLabTarget << "GeV, bLabTarget = " << bLabTarget << ", yLabTarget = " << yLabTarget << endl;

  // center of mass
  double eCtrMassNN = TMath::Sqrt(massBeam*massBeam + massTarget*massTarget + 2.0*eLabBeam*eLabTarget - 2.0*pLabBeam*pLabTarget);
  double eCtrMass   = TMath::Sqrt(massBeam*massBeam + massTarget*massTarget + 2.0*eLabBeam*eLabTarget - 2.0*pLabBeam*pLabTarget)*TMath::Sqrt(1.0/numNucleons/numNucleons);
  double pCtrMass = pLabBeam*massTarget/eCtrMassNN;

  double eCtrBeam = TMath::Sqrt(pCtrMass*pCtrMass + massBeam*massBeam);
  double pCtrBeam = pCtrMass;
  double bCtrBeam = pCtrMass/eCtrBeam;
  double yCtrBeam = 0.5*TMath::Log((eCtrBeam+pCtrBeam)/(eCtrBeam-pCtrBeam));

  double eCtrTarget = TMath::Sqrt(pCtrMass*pCtrMass + massTarget*massTarget);
  double pCtrTarget = -1.0*pCtrMass;
  double bCtrTarget = pCtrTarget/eCtrTarget;
  double yCtrTarget = 0.5*TMath::Log((eCtrTarget+pCtrTarget)/(eCtrTarget-pCtrTarget));

  double yCtrMass = yLabBeam - yCtrBeam;
  cout << "eCtrBeam = " << eCtrBeam << "GeV, pCtrBeam = " << pCtrBeam << "GeV, bCtrBeam = " << bCtrBeam << ", yCtrBeam = " << yCtrBeam << endl;
  cout << "eCtrTarget = " << eCtrTarget << "GeV, pCtrTarget = " << pCtrTarget << "GeV, bCtrTarget = " << bCtrTarget << ", yCtrTarget = " << yCtrTarget << endl;
  cout << "eCtrMass = " << eCtrMass << "GeV, pCtrMass = " << pCtrMass << "GeV, yCtrMass = " << yCtrMass << endl;
}
