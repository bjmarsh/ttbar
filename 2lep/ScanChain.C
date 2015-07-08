// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TCanvas.h"

// CORE
#include "/home/users/bmarsh/CORE/CMS3.h"
#include "/home/users/bmarsh/CORE/ElectronSelections.h"
#include "/home/users/bmarsh/CORE/MuonSelections.h"
#include "/home/users/bmarsh/CORE/JetSelections.h"


//nice plots
#include "/home/users/bmarsh/Software/dataMCplotMaker/dataMCplotMaker.h"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1) {

    // Benchmark
    TBenchmark *bmark = new TBenchmark();
    bmark->Start("benchmark");

    TCanvas *c = new TCanvas("c","",800,600);
    c->cd();

    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

    // mLL
    TH1F *mLL = new TH1F("mLL", "", 300,0,300);
    mLL->SetDirectory(rootdir);

    // pfmet
    TH1F *pfmet = new TH1F("pfmet", "", 300,0,300);
    pfmet->SetDirectory(rootdir);

    // njets
    TH1F *njets = new TH1F("njets", "", 13,-0.5,12.5);
    njets->SetDirectory(rootdir);

    // nbtags
    TH1F *nbtags = new TH1F("nbtags", "", 7,-0.5,6.5);
    nbtags->SetDirectory(rootdir);

    // jetpt
    TH1F *jetpt = new TH1F("jetpt", "", 300,0,300);
    jetpt->SetDirectory(rootdir);

    // bjetpt
    TH1F *bjetpt = new TH1F("bjetpt", "", 300,0,300);
    bjetpt->SetDirectory(rootdir);

    // H_T
    TH1F *HT = new TH1F("HT","", 300,0,900);
    HT->SetDirectory(rootdir);

    // leppt
    TH1F *leppt = new TH1F("leppt","", 300,0,300);
    leppt->SetDirectory(rootdir);


    // Loop over events to Analyze
    unsigned int nEventsTotal = 0;
    unsigned int nEventsChain = chain->GetEntries();
    if( nEvents >= 0 ) nEventsChain = nEvents;
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    int totalEvents = 4992231;

    // File Loop
    while ( (currentFile = (TFile*)fileIter.Next()) ) {

        // Get File Content
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");
        if(fast) TTreeCache::SetLearnEntries(10);
        if(fast) tree->SetCacheSize(128*1024*1024);
        cms3.Init(tree);
    

        // Loop over Events in current file
        if( nEventsTotal >= nEventsChain ) continue;
        unsigned int nEventsTree = tree->GetEntriesFast();
        for( unsigned int event = 0; event < nEventsTree; ++event) {
    
            // Get Event Content
            if( nEventsTotal >= nEventsChain ) continue;
            if(fast) tree->LoadTree(event);
            cms3.GetEntry(event);
            ++nEventsTotal;

            // get the proper scale factors
            double scale;
            scale = evt_scale1fb()*10*totalEvents/nEventsChain;
    
            // Progress
            CMS3::progress( nEventsTotal, nEventsChain );            

            /******** Analysis Code *********/
            
            bool isEnoughMET = (evt_pfmet()>30);

            int nhyps = hyp_ll_p4().size();
            int ngoodhyps = 0;
            int besthyp = -1;
            double bestPt = 0;
            for(int i=0; i<nhyps; i++){
                if(hyp_ll_charge()[i] + hyp_lt_charge()[i] != 0)
                    continue;
                if(hyp_ll_p4()[i].Pt() < 20 || hyp_lt_p4()[i].Pt() < 20)
                    continue;
                if(fabs(hyp_ll_p4()[i].Eta()) > 2.4 || fabs(hyp_lt_p4()[i].Eta()) > 2.4)
                    continue;
                
                ngoodhyps++;
                
                double totalPt = hyp_ll_p4()[i].Pt() + hyp_lt_p4()[i].Pt();
                if(totalPt>bestPt){
                    bestPt = totalPt;
                    besthyp = i;
                }
                
            }
            if(ngoodhyps==0)
                continue;

            int njetstotal = pfjets_p4().size();
            int ngoodjets = 0;
            int ngoodbtags = 0;
            double h_t = 0;
            for(int i=0; i<njetstotal; i++){

                bool isBtag = (pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag()[i] > 0.814);
                
                double jetEta = pfjets_p4()[i].Eta();
                double jetPhi = pfjets_p4()[i].Phi();
                double l1Eta = hyp_ll_p4()[besthyp].Eta();
                double l2Eta = hyp_lt_p4()[besthyp].Eta();
                double l1Phi = hyp_ll_p4()[besthyp].Phi();
                double l2Phi = hyp_lt_p4()[besthyp].Phi();
                double dR1 = sqrt((jetEta-l1Eta)*(jetEta-l1Eta)+(jetPhi-l1Phi)*(jetPhi-l1Phi));
                double dR2 = sqrt((jetEta-l2Eta)*(jetEta-l2Eta)+(jetPhi-l2Phi)*(jetPhi-l2Phi));
                if(dR1 < 0.4 || dR2 < 0.4)
                    continue;

                jetpt->Fill(pfjets_p4()[i].Pt(), scale);
                if(isBtag)
                    bjetpt->Fill(pfjets_p4()[i].Pt(), scale);                

                if(pfjets_p4()[i].Pt() < 40)
                    continue;
                if(fabs(pfjets_p4()[i].Eta()) > 2.4)
                    continue;

                ngoodjets++;
                if(isBtag)
                    ngoodbtags++;

                h_t += pfjets_p4()[i].Pt();
            }

            if(isEnoughMET){
                njets->Fill(ngoodjets, scale);
                nbtags->Fill(ngoodbtags, scale);
            }

            if(ngoodjets < 2)
                continue;
            if(ngoodbtags < 1)
                continue;

            pfmet->Fill(evt_pfmet(), scale);

            if(!isEnoughMET)
                continue;

            int nElecTotal = els_p4().size();
            int nMuonTotal = mus_p4().size();

            for(int i=0; i<nElecTotal; i++){
                leppt->Fill(els_p4()[i].Pt());
            }
            for(int i=0; i<nMuonTotal; i++){
                leppt->Fill(mus_p4()[i].Pt());
            }


            mLL->Fill((hyp_ll_p4()[besthyp]+hyp_lt_p4()[besthyp]).M(), scale);          
            HT->Fill(h_t, scale);
            
            for(int i=0; i<njetstotal; i++){

                bool isBtag = (pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag()[i] > 0.814);
                
                double jetEta = pfjets_p4()[i].Eta();
                double jetPhi = pfjets_p4()[i].Phi();
                double l1Eta = hyp_ll_p4()[besthyp].Eta();
                double l2Eta = hyp_lt_p4()[besthyp].Eta();
                double l1Phi = hyp_ll_p4()[besthyp].Phi();
                double l2Phi = hyp_lt_p4()[besthyp].Phi();
                double dR1 = sqrt((jetEta-l1Eta)*(jetEta-l1Eta)+(jetPhi-l1Phi)*(jetPhi-l1Phi));
                double dR2 = sqrt((jetEta-l2Eta)*(jetEta-l2Eta)+(jetPhi-l2Phi)*(jetPhi-l2Phi));
                if(dR1 < 0.4 || dR2 < 0.4)
                    continue;

                if(fabs(pfjets_p4()[i].Eta()) > 2.4)
                    continue;

                jetpt->Fill(pfjets_p4()[i].Pt(), scale);
                if(isBtag)
                    bjetpt->Fill(pfjets_p4()[i].Pt(), scale);                

            }


        }
  
        // Clean Up
        delete tree;
        file->Close();
        delete file;
    }
    if ( nEventsChain != nEventsTotal ) {
        cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
    }
  

    // mLL
    std::vector<TH1F*> bkg;
    bkg.push_back(mLL);
    std::vector<string> str;
    str.push_back("TTJets");
    TH1F *null = new TH1F("","",1,0,1);
    dataMCplotMaker(null, bkg, str, "mLL", "", "--outputName tt_mLL.pdf --xAxisLabel m_{LL}  --setMinimum 100" );

    // pfmet
    bkg.pop_back();
    bkg.push_back(pfmet);
    str.pop_back();
    str.push_back("TTJets");
    dataMCplotMaker(null, bkg, str, "pfmet", "", "--outputName tt_pfmet.pdf --xAxisLabel MET  --setMinimum 100" );

    // njets
    bkg.pop_back();
    bkg.push_back(njets);
    str.pop_back();
    str.push_back("TTJets");
    dataMCplotMaker(null, bkg, str, "njets", "", "--outputName tt_njets.pdf --xAxisLabel njets  --nDivisions 13 --noXaxisUnit --setMinimum 10" );

    // nbtags
    bkg.pop_back();
    bkg.push_back(nbtags);
    str.pop_back();
    str.push_back("TTJets");
    dataMCplotMaker(null, bkg, str, "nbtags", "", "--outputName tt_nbtags.pdf --xAxisLabel nbtags  --nDivisions 7 --noXaxisUnit" );

    // jetpt
    bkg.pop_back();
    bkg.push_back(jetpt);
    str.pop_back();
    str.push_back("TTJets");
    dataMCplotMaker(null, bkg, str, "jet pT", "", "--outputName tt_jetpt.pdf --xAxisLabel p_{T} --setMinimum 100" );

    // bjetpt
    bkg.pop_back();
    bkg.push_back(bjetpt);
    str.pop_back();
    str.push_back("TTJets");
    dataMCplotMaker(null, bkg, str, "b-jet pT", "", "--outputName tt_bjetpt.pdf --xAxisLabel p_{T}  --setMinimum 10" );

    // H_T
    bkg.pop_back();
    bkg.push_back(HT);
    str.pop_back();
    str.push_back("TTJets");
    dataMCplotMaker(null, bkg, str, "H_{T}", "", "--outputName tt_H_T.pdf --xAxisLabel H_{T}  --setMinimum 10" );

    // leppt
    bkg.pop_back();
    bkg.push_back(leppt);
    str.pop_back();
    str.push_back("TTJets");
    dataMCplotMaker(null, bkg, str, "lepton p_{T}", "", "--outputName tt_leppt.pdf --xAxisLabel p_{T} --setMinimum 10 " );


    // return
    bmark->Stop("benchmark");
    cout << endl;
    cout << nEventsTotal << " Events Processed" << endl;
    cout << "------------------------------" << endl;
    cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
    cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
    cout << endl;
    delete bmark;
    return 0;
}
