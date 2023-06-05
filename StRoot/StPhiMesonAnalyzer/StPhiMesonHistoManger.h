#ifndef StPhiMesonHistoManger_h
#define StPhiMesonHistoManger_h

#include <map>
#include <string>

class TH1F;
class TH2F;

class StPhiMesonHistoManger
{
  public:
    StPhiMesonHistoManger(const int beamType, const int flagME);
    virtual ~StPhiMesonHistoManger();

    // phi QA
    void initPhiQA();
    void fillPhiQA(int cent9, double pt, double yLab, double yCms, double invMass, double refWgt);
    void writePhiQA();
    int getPtBinPhiQA(double pt);
    int getRapBinPhiQA(double y);

    // phi Flow
    bool is0080(int cent9);
    bool is2060(int cent9);

    void initIsoPhiFlow(); // IsoBar
    void fillIsoPhiV2(int cent9, double pt, double yCms, double phi, double Psi2, double invMass, double res2, double refWgt);
    void fillIsoPhiV3(int cent9, double pt, double yCms, double phi, double Psi3, double invMass, double res3, double refWgt);
    void fillIsoPhiYields(int cent9, double pt, double yCms, double invMass, double refWgt);
    void writeIsoPhiFlow();
    int getCentBinIsoPhiFlow(int cent9);
    int getPtBinIsoPhiFlow(double pt);
    int getPsi2BinIsoPhiFlow(double phi, double Psi2);
    int getPsi3BinIsoPhiFlow(double phi, double Psi3);
    double transPsi2(double phi, double Psi2);
    double transPsi3(double phi, double Psi3);
    bool isPsi2InRange(double phiPsi2);
    bool isPsi3InRange(double phiPsi3);

    void initFxtPhiFlow(); // Fxt3p85_2018
    void fillFxtPhiV2(int cent9, double pt, double yCms, double phi, double Psi1, double invMass, double res12, double refWgt);
    void fillFxtPhiYields(int cent9, double pt, double yCms, double invMass, double refWgt);
    void writeFxtPhiFlow();
    int getCentBinFxtPhiFlow(int cent9);
    int getPtBinFxtPhiFlow(double pt);
    int getPsi12BinFxtPhiFlow(double phi, double Psi1);
    double transPsi12(double phi, double Psi1);
    bool isPsi12InRange(double phiPsi12);

  private:
    static const int mNumCentBinQA = 10; // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%, 9: 20-60%
    static const int mNumPtBinQA   = 25; // 25 bins from 0 to 5.0 GeV/c
    static const int mNumRapBinQA  = 25; // 25 bins from -1.25 to 1.25

    static const int mNumCentBinIsoPhiFlow   = 4;  // 0: 0%-80%, 1: 0-10%, 2: 10-40%, 3: 40-80%
    static const int mNumPtBinIsoPhiFlow     = 20; // 0.4 - 7.2 GeV/c
    static const int mNumPsiBinIsoPhiFlow    = 7;  // [0,pi/2] for v2 | [0,pi/3] for v3
    static const int mNumCentBinIsoPhiYileds = 9;  // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%

    static const int mNumCentBinFxtPhiFlow   = 4;  // 0: 0%-60%, 1: 0-10%, 2: 10-40%, 3: 40-60%
    static const int mNumPtBinFxtPhiFlow     = 10; // 0.4 - 3.0 GeV/c
    static const int mNumPsiBinFxtPhiFlow    = 5;  // [0,pi/2] for v2
    static const int mNumCentBinFxtPhiYileds = 9;  // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%

    // QA histograms
    TH2F* h_mInvMassPhiQA[mNumCentBinQA]; // pT vs. invMass
    // Acceptance Histograms
    // 0 = centrality: 9 = 20%-60%, 0-8 from StRefMultCorr 
    // 1 = pt bin: 25 bins from 0 to 5.0 GeV/c
    // 2 = rapidity bin: 25 bins from -1.25 to 1.25
    std::map<string,TH1F*> h_mAcptPhiLab;
    std::map<string,TH1F*> h_mAcptPhiCms;

    // Isobar phi Flow histograms
    // 0 = centrality: 0: 0%-80%, 1: 0-10%, 2: 10-40%, 3: 40-80% for Isobar
    // 1 = pt bin
    // 2 = phi-Psi: 7 bins for Isobar v2 & v3
    std::map<string,TH1F*> h_mInvMassIsoPhiV2;
    std::map<string,TH1F*> h_mInvMassIsoPhiV3;
    // Isobar phi Yields histograms
    // 0 = centrality: 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
    // 1 = pt bin
    std::map<string,TH1F*> h_mInvMassIsoPhiYields;

    // FXT phi Flow histograms
    // 0 = centrality: 0: 0-60%, 1: 0-10%, 2: 10-40%, 3: 40-60%
    // 1 = pt bin
    // 2 = phi-Psi: 5 bins for FXT v2
    std::map<string,TH1F*> h_mInvMassFxtPhiV2;
    // Isobar phi Yields histograms
    // 0 = centrality: 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
    // 1 = pt bin
    std::map<string,TH1F*> h_mInvMassFxtPhiYields;

    const int mType;
    const int mFlagME; // 0 for Same Event, 1 for Mixed Event

    std::string str_mMixEvt[2] = {"SE","ME"};

  ClassDef(StPhiMesonHistoManger,1)
};
#endif
