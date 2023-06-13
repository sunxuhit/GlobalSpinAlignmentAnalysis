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

void calcLabEnergy(double cmsEnergy = 3.0)
{
  const double massProton = anaUtils::mMassProton;
  const double eLabBeam = (cmsEnergy*cmsEnergy-2.0*massProton*massProton)/(2.0*massProton);
  const double bLabBeam = TMath::Sqrt(eLabBeam*eLabBeam-massProton*massProton)/(eLabBeam+massProton);

  cout << "cmsEnergy = " << cmsEnergy << ", eLabBeam = " << eLabBeam << ", bLabBeam = " << bLabBeam << endl;
}
