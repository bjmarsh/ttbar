{
    gSystem->Load("/home/users/bmarsh/CORE/CMS3_CORE.so");
    gROOT->ProcessLine(".L /home/users/bmarsh/Software/dataMCplotMaker/dataMCplotMaker.cc++");
    gROOT->ProcessLine(".L ScanChain.C++");
    

    TChain *ch = new TChain("Events"); 

    ch->Add("/hadoop/cms/store/group/snt/run2_50ns/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/V07-04-03/*1*.root");
    
    ScanChain(ch); 
}
