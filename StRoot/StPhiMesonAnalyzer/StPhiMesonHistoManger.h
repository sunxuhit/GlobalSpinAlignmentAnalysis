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

    void initPhiQA();
    void fillPhiQA(int cent9, double pt, double yLab, double yCms, double invMass, double refWgt);
    void writePhiQA();

    int getPtBinQA(double pt);
    int getRapBinQA(double y);
    bool is2060(int cent9);

  private:
    static const int mNumCentrality = 10; // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%, 9: 20-60%
    static const int mNumPtBinQA    = 25; // 25 bins from 0 to 5.0 GeV/c
    static const int mNumRapBinQA   = 25; // 25 bins from -1.25 to 1.25

    // QA histograms
    TH2F* h_mInvMassPhi[mNumCentrality]; // pT vs. invMass
    // Acceptance Histograms
    // 0 = centrality: 9 = 20%-60%, 0-8 from StRefMultCorr 
    // 1 = pt bin: 25 bins from 0 to 5.0 GeV/c
    // 2 = rapidity bin: 25 bins from -1.25 to 1.25
    std::map<string,TH1F*> h_mAcptPhiLab;
    std::map<string,TH1F*> h_mAcptPhiCms;

    // Flow histograms
    // 0 = centrality: 9 = 20%-60%, 0-8 from StRefMultCorr 
    // 1 = pt bin
    // 2 = phi-Psi1
    std::map<string,TH1F*> h_mInvMassPhiFlow;

    const int mType;
    const int mFlagME; // 0 for Same Event, 1 for Mixed Event

    std::string str_mMixEvt[2] = {"SE","ME"};

  ClassDef(StPhiMesonHistoManger,1)
};
#endif
