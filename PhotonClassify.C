#define PhotonClassify_cxx
#include "PhotonClassify.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void PhotonClassify::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PhotonClassify.C
//      root> PhotonClassify t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

	TFile file_in("fit_function.root");
	file_in.ls();
	TF1 *p1eb_fit2 = (TF1*) file_in.Get("p1eb_fit2");
	p1eb_fit2->Print();

	int p1eb_fit2_good_count = 0;
	int p1eb_fit2_bad_count = 0;
    int p1eb_fit2_bad_high_count = 0;
	int p1eb_fit2_bad_low_count = 0;

	TFile file_out("PhotonClassifyLoopOutput.root","RECREATE");

        TH1D *Photon1_EB_fit2_good_res_hist = new TH1D("Photon1_EB_fit2_good_res_hist","Photon 1 Barrel Resolution Good Fit2",100,-1,1);
        Photon1_EB_fit2_good_res_hist->GetXaxis()->SetTitle("(Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt");
        Photon1_EB_fit2_good_res_hist->GetXaxis()->CenterTitle();
        Photon1_EB_fit2_good_res_hist->GetYaxis()->SetTitle("Number of events");
        Photon1_EB_fit2_good_res_hist->GetYaxis()->CenterTitle();

        TH1D *Photon1_EB_fit2_bad_res_hist = new TH1D("Photon1_EB_fit2_bad_res_hist","Photon 1 Barrel Resolution Bad Fit2",100,-1,1);
        Photon1_EB_fit2_bad_res_hist->GetXaxis()->SetTitle("(Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt");
        Photon1_EB_fit2_bad_res_hist->GetXaxis()->CenterTitle();
        Photon1_EB_fit2_bad_res_hist->GetYaxis()->SetTitle("Number of events");
        Photon1_EB_fit2_bad_res_hist->GetYaxis()->CenterTitle();

        TH1D *Photon1_EB_fit2_bad_high_res_hist = new TH1D("Photon1_EB_fit2_bad_high_res_hist","Photon 1 Barrel Resolution Bad High Fit2",100,-1,1);
        Photon1_EB_fit2_bad_high_res_hist->GetXaxis()->SetTitle("(Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt");
        Photon1_EB_fit2_bad_high_res_hist->GetXaxis()->CenterTitle();
        Photon1_EB_fit2_bad_high_res_hist->GetYaxis()->SetTitle("Number of events");
        Photon1_EB_fit2_bad_high_res_hist->GetYaxis()->CenterTitle();

		TH1D *Photon1_EB_fit2_bad_low_res_hist = new TH1D("Photon1_EB_fit2_bad_low_res_hist","Photon 1 Barrel Resolution Bad Low Fit2",100,-1,1);
        Photon1_EB_fit2_bad_low_res_hist->GetXaxis()->SetTitle("(Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt");
        Photon1_EB_fit2_bad_low_res_hist->GetXaxis()->CenterTitle();
        Photon1_EB_fit2_bad_low_res_hist->GetYaxis()->SetTitle("Number of events");
        Photon1_EB_fit2_bad_low_res_hist->GetYaxis()->CenterTitle();

        TH1D *Photon1_EB_res_hist = new TH1D("Photon1_EB_res_hist","Photon 1 Barrel Resolution",100,-1,1);
        Photon1_EB_res_hist->GetXaxis()->SetTitle("(Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt");
        Photon1_EB_res_hist->GetXaxis()->CenterTitle();
        Photon1_EB_res_hist->GetYaxis()->SetTitle("Number of events");
        Photon1_EB_res_hist->GetYaxis()->CenterTitle();

        TH2D *Photon1_EB_fit2_bad_eta_phi = new TH2D("Photon1_EB_fit2_bad_eta_phi","Photon 1 Barrel Fit2 Bad Eta Phi",100,-3,3,100,-3.15,3.15);
        Photon1_EB_fit2_bad_eta_phi->GetXaxis()->SetTitle("scEta");
        Photon1_EB_fit2_bad_eta_phi->GetXaxis()->CenterTitle();
        Photon1_EB_fit2_bad_eta_phi->GetYaxis()->SetTitle("scPhi");
        Photon1_EB_fit2_bad_eta_phi->GetYaxis()->CenterTitle();

        TH2D *Photon1_EB_fit2_bad_low_eta_phi = new TH2D("Photon1_EB_fit2_bad_low_eta_phi","Photon 1 Barrel Fit2 Bad Low Eta Phi",100,-3,3,100,-3.15,3.15);
        Photon1_EB_fit2_bad_low_eta_phi->GetXaxis()->SetTitle("scEta");
        Photon1_EB_fit2_bad_low_eta_phi->GetXaxis()->CenterTitle();
        Photon1_EB_fit2_bad_low_eta_phi->GetYaxis()->SetTitle("scPhi");
        Photon1_EB_fit2_bad_low_eta_phi->GetYaxis()->CenterTitle();

        TH2D *Photon1_EB_fit2_bad_high_eta_phi = new TH2D("Photon1_EB_fit2_bad_high_eta_phi","Photon 1 Barrel Fit2 Bad High Eta Phi",100,-3,3,100,-3.15,3.15);
        Photon1_EB_fit2_bad_high_eta_phi->GetXaxis()->SetTitle("scEta");
        Photon1_EB_fit2_bad_high_eta_phi->GetXaxis()->CenterTitle();
        Photon1_EB_fit2_bad_high_eta_phi->GetYaxis()->SetTitle("scPhi");
        Photon1_EB_fit2_bad_high_eta_phi->GetYaxis()->CenterTitle();

        TH2D *Photon1_EB_fit2_good_eta_phi = new TH2D("Photon1_EB_fit2_good_eta_phi","Photon 1 Barrel Fit2 Good Eta Phi",100,-3,3,100,-3.15,3.15);
        Photon1_EB_fit2_good_eta_phi->GetXaxis()->SetTitle("scEta");
        Photon1_EB_fit2_good_eta_phi->GetXaxis()->CenterTitle();
        Photon1_EB_fit2_good_eta_phi->GetYaxis()->SetTitle("scPhi");
        Photon1_EB_fit2_good_eta_phi->GetYaxis()->CenterTitle();

	TTree *PhotonClassifyTree = new TTree("PhotonClassifyTree","Photon Classify Tree");
	// if p1eb_bad_classify = 1 then it is a low bad photon
	// if p1eb_classify = 1 then it is bad in general
	// we only want photons with bad classify = 1 and classify = 0 (so no high bad photons)
    Bool_t p1eb_fit2_classify, p1eb_bad_classify;
    Double_t Photon1_res, p1eb_fit2sigma;
	
   Double_t        out_Photon1_pt;
   Double_t        out_Photon1_eta;
   Double_t        out_Photon1_phi;
   Double_t        out_Photon1_scEta;
   Double_t        out_Photon1_scPhi;
   Double_t        out_Photon1_rho;
   Double_t        out_Photon1_chargedHadIso03;
   Double_t        out_Photon1_neutralHadIso03;
   Double_t        out_Photon1_photonIso03;
   Double_t        out_Photon1_rhoCorChargedHadIso03;
   Double_t        out_Photon1_rhoCorNeutralHadIso03;
   Double_t        out_Photon1_rhoCorPhotonIso03;
   Double_t        out_Photon1_corPhotonIso03;
   Double_t        out_Photon1_hadTowerOverEm;
   Double_t        out_Photon1_hadronicOverEm;
   Double_t        out_Photon1_r9;
   Double_t        out_Photon1_r9_5x5;
   Double_t        out_Photon1_sigmaIetaIeta;
   Double_t        out_Photon1_sigmaIetaIeta5x5;
   Double_t        out_Photon1_sigmaEtaEta;
   Double_t        out_Photon1_sigmaIphiIphi;
   Double_t        out_Photon1_sigmaIphiIphi5x5;
   Double_t        out_Photon1_sigmaIetaIphi;
   Double_t        out_Photon1_sigmaIetaIphi5x5;
   Double_t        out_Photon1_maxEnergyXtal;
   Double_t        out_Photon1_iEta;
   Double_t        out_Photon1_iPhi;
//   Double_t        out_Photon1_alphaHighPtID;
//   Double_t        out_Photon1_kappaHighPtID;
//   Double_t        out_Photon1_phoEAHighPtID;
//   Double_t        out_Photon1_chEAegmID;
//   Double_t        out_Photon1_nhEAegmID;
//   Double_t        out_Photon1_phoEAegmID;
   Bool_t          out_Photon1_passEGMLooseID;
   Bool_t          out_Photon1_passEGMMediumID;
   Bool_t          out_Photon1_passEGMTightID;
   Bool_t          out_Photon1_isEB;
   //Bool_t          Photon1_isEE;
   Bool_t          out_Photon1_isEBEtaGap;
   Bool_t          out_Photon1_isEBPhiGap;
   //Bool_t          Photon1_isEERingGap;
   //Bool_t          Photon1_isEEDeeGap;
   Bool_t          out_Photon1_isEBEEGap;
   Bool_t          out_Photon1_passElectronVeto;
   Bool_t          out_Photon1_passHTowOverE;
   Bool_t          out_Photon1_passChIso;
   Bool_t          out_Photon1_passCorPhoIso;
   Bool_t          out_Photon1_passSieie;
   Bool_t          out_Photon1_passHighPtID;
//   Bool_t          out_Photon1_passChIsoDenom;
//   Bool_t          out_Photon1_passCorPhoIsoDenom;
//   Bool_t          out_Photon1_isFakeable;
//   Bool_t          out_Photon1_isNumeratorObjCand;
//   Bool_t          out_Photon1_isDenominatorObj;
   Bool_t          out_Photon1_isSaturated;
   Int_t           out_Event_npv_true;
   Int_t           out_nPV;

	PhotonClassifyTree->Branch("Photon1_pt",&out_Photon1_pt,"Photon1_pt/D");
	PhotonClassifyTree->Branch("Photon1_eta",&out_Photon1_eta,"Photon1_eta/D");
	PhotonClassifyTree->Branch("Photon1_phi",&out_Photon1_phi,"Photon1_phi/D");
	PhotonClassifyTree->Branch("Photon1_scEta",&out_Photon1_scEta,"Photon1_scEta/D");
	PhotonClassifyTree->Branch("Photon1_scPhi",&out_Photon1_scPhi,"Photon1_scPhi/D");
	PhotonClassifyTree->Branch("Photon1_rho",&out_Photon1_rho,"Photon1_rho/D");
	PhotonClassifyTree->Branch("Photon1_chargedHadIso03",&out_Photon1_chargedHadIso03,"Photon1_chargedHadIso03/D");
	PhotonClassifyTree->Branch("Photon1_neutralHadIso03",&out_Photon1_neutralHadIso03,"Photon1_neutralHadIso03/D");
	PhotonClassifyTree->Branch("Photon1_photonIso03",&out_Photon1_photonIso03,"Photon1_photonIso03/D");
	PhotonClassifyTree->Branch("Photon1_rhoCorChargedHadIso03",&out_Photon1_rhoCorChargedHadIso03,"Photon1_rhoCorChargedHadIso03/D");
	PhotonClassifyTree->Branch("Photon1_rhoCorNeutralHadIso03",&out_Photon1_rhoCorNeutralHadIso03,"Photon1_rhoCorNeutralHadIso03/D");
	PhotonClassifyTree->Branch("Photon1_rhoCorPhotonIso03",&out_Photon1_rhoCorPhotonIso03,"Photon1_rhoCorPhotonIso03/D");
	PhotonClassifyTree->Branch("Photon1_corPhotonIso03",&out_Photon1_corPhotonIso03,"Photon1_corPHotonIso03/D");
	PhotonClassifyTree->Branch("Photon1_hadTowerOverEm",&out_Photon1_hadTowerOverEm,"Photon1_hadTowerOverEm/D");
	PhotonClassifyTree->Branch("Photon1_hadronicOverEm",&out_Photon1_hadronicOverEm,"Photon1_hadronicOverEm/D");
	PhotonClassifyTree->Branch("Photon1_r9",&out_Photon1_r9,"Photon1_r9/D");
	PhotonClassifyTree->Branch("Photon1_r9_5x5",&out_Photon1_r9_5x5,"Photon1_r9_5x5/D");
	PhotonClassifyTree->Branch("Photon1_sigmaIetaIeta",&out_Photon1_sigmaIetaIeta,"Photon1_sigmaIetaIeta/D");
	PhotonClassifyTree->Branch("Photon1_sigmaIetaIeta5x5",&out_Photon1_sigmaIetaIeta5x5,"Photon1_sigmaIetaIeta5x5/D");
	PhotonClassifyTree->Branch("Photon1_sigmaEtaEta",&out_Photon1_sigmaEtaEta,"Photon1_sigmaEtaEta/D");
	PhotonClassifyTree->Branch("Photon1_sigmaIphiIphi",&out_Photon1_sigmaIphiIphi,"Photon1_sigmaIphiIphi/D");
	PhotonClassifyTree->Branch("Photon1_sigmaIphiIphi5x5",&out_Photon1_sigmaIphiIphi5x5,"Photon1_sigmaIphiIphi5x5/D");
	PhotonClassifyTree->Branch("Photon1_sigmaIetaIphi",&out_Photon1_sigmaIetaIphi,"Photon1_sigmaIetaIphi/D");
	PhotonClassifyTree->Branch("Photon1_sigmaIetaIphi5x5",&out_Photon1_sigmaIetaIphi5x5,"Photon1_sigmaIetaIphi5x5/D");
	PhotonClassifyTree->Branch("Photon1_maxEnergyXtal",&out_Photon1_maxEnergyXtal,"Photon1_maxEnergyXtal/D");
	PhotonClassifyTree->Branch("Photon1_iEta",&out_Photon1_iEta,"Photon1_iEta/D");
	PhotonClassifyTree->Branch("Photon1_iPhi",&out_Photon1_iPhi,"Photon1_iPhi/D");
//	PhotonClassifyTree->Branch("Photon1_alphaHighPtID",&out_Photon1_alphaHighPtID,"Photon1_alphaHighPtID/D");
//	PhotonClassifyTree->Branch("Photon1_kappaHighPtID",&out_Photon1_kappaHighPtID,"Photon1_kappaHighPtID/D");
//	PhotonClassifyTree->Branch("Photon1_phoEAHighPtID",&out_Photon1_phoEAHighPtID,"Photon1_phoEAHighPtID/D");
//	PhotonClassifyTree->Branch("Photon1_chEAegmID",&out_Photon1_chEAegmID,"Photon1_chEAegmID/D");
//	PhotonClassifyTree->Branch("Photon1_nhEAegmID",&out_Photon1_nhEAegmID,"Photon1_nhEAegmID/D");
//	PhotonClassifyTree->Branch("Photon1_phoEAegmID",&out_Photon1_phoEAegmID,"Photon1_phoEAegmID/D");
	PhotonClassifyTree->Branch("Photon1_passEGMLooseID",&out_Photon1_passEGMLooseID,"Photon1_passEGMLooseID/O");
	PhotonClassifyTree->Branch("Photon1_passEGMMediumID",&out_Photon1_passEGMMediumID,"Photon1_passEGMMediumID/O");
	PhotonClassifyTree->Branch("Photon1_passEGMTightID",&out_Photon1_passEGMTightID,"Photon1_passEGMTightID/O");
	PhotonClassifyTree->Branch("Photon1_isEB",&out_Photon1_isEB,"Photon1_isEB/O");
	PhotonClassifyTree->Branch("Photon1_isEBEtaGap",&out_Photon1_isEBEtaGap,"Photon1_isEBEtaGap/O");
	PhotonClassifyTree->Branch("Photon1_isEBPhiGap",&out_Photon1_isEBPhiGap,"Photon1_isEBPhiGap/O");
	PhotonClassifyTree->Branch("Photon1_isEBEEGap",&out_Photon1_isEBEEGap,"Photon1_isEBEEGap/O");
	PhotonClassifyTree->Branch("Photon1_passElectronVeto",&out_Photon1_passElectronVeto,"Photon1_passElectronVeto/O");
	PhotonClassifyTree->Branch("Photon1_passHTowOverE",&out_Photon1_passHTowOverE,"Photon1_passHTowOverE/O");
	PhotonClassifyTree->Branch("Photon1_passChIso",&out_Photon1_passChIso,"Photon1_passChIso/O");
	PhotonClassifyTree->Branch("Photon1_passCorPhoIso",&out_Photon1_passCorPhoIso,"Photon1_passCorPhoIso/O");
	PhotonClassifyTree->Branch("Photon1_passSieie",&out_Photon1_passSieie,"Photon1_passSieie/O");
	PhotonClassifyTree->Branch("Photon1_passHighPtID",&out_Photon1_passHighPtID,"Photon1_passHighPtID/O");
//	PhotonClassifyTree->Branch("Photon1_passChIsoDenom",&out_Photon1_passChIsoDenom,"Photon1_passChIsoDenom/O");
//	PhotonClassifyTree->Branch("Photon1_passCorPhoIsoDenom",&out_Photon1_passCorPhoIsoDenom,"Photon1_passCorPhoIsoDenom/O");
//	PhotonClassifyTree->Branch("Photon1_isFakeable",&out_Photon1_isFakeable,"Photon1_isFakeable/O");
//	PhotonClassifyTree->Branch("Photon1_isNumeratorObjCand",&out_Photon1_isNumeratorObjCand,"Photon1_isNumeratorObjCand/O");
//	PhotonClassifyTree->Branch("Photon1_isDenominatorObj",&out_Photon1_isDenominatorObj,"Photon1_isDenominatorObj/O");
	PhotonClassifyTree->Branch("Photon1_isSaturated",&out_Photon1_isSaturated,"Photon1_isSaturated/O");
	PhotonClassifyTree->Branch("Event_npv_true",&out_Event_npv_true,"Event_npv_true/I");
	PhotonClassifyTree->Branch("nPV",&out_nPV,"nPV/I");
	PhotonClassifyTree->Branch("Photon1_res",&Photon1_res,"Photon1_res/D");
//	PhotonClassifyTree->Branch("p1eb_classify",&p1eb_classify,"p1eb_classify/O");
        PhotonClassifyTree->Branch("p1eb_fit2_classify",&p1eb_fit2_classify,"p1eb_fit2_classify/O");
//	PhotonClassifyTree->Branch("p1eb_fitsigma",&p1eb_fitsigma,"p1eb_fitsigma/D");
        PhotonClassifyTree->Branch("p1eb_fit2sigma",&p1eb_fit2sigma,"p1eb_fit2sigma/D");
	PhotonClassifyTree->Branch("p1eb_bad_classify",&p1eb_bad_classify,"p1eb_bad_classify/O");


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
	if (jentry%10000 == 0) cout << "Entry " << jentry << endl;
	
	Photon1_res = ((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt);
        p1eb_fit2sigma = p1eb_fit2->Eval(GenPhoton1_pt);

	out_Photon1_pt = Photon1_pt;
	out_Photon1_eta = Photon1_eta;
	out_Photon1_phi = Photon1_phi;
	out_Photon1_scEta = Photon1_scEta;
	out_Photon1_scPhi = Photon1_scPhi;
	out_Photon1_rho = Photon1_rho;
	out_Photon1_chargedHadIso03 = Photon1_chargedHadIso03;
	out_Photon1_neutralHadIso03 = Photon1_neutralHadIso03;
	out_Photon1_photonIso03 = Photon1_photonIso03;
	out_Photon1_rhoCorChargedHadIso03 = Photon1_rhoCorChargedHadIso03;
	out_Photon1_rhoCorNeutralHadIso03 = Photon1_rhoCorNeutralHadIso03;
        out_Photon1_rhoCorPhotonIso03 = Photon1_rhoCorPhotonIso03;
        out_Photon1_corPhotonIso03 = Photon1_corPhotonIso03;
        out_Photon1_hadTowerOverEm = Photon1_hadTowerOverEm;
        out_Photon1_hadronicOverEm = Photon1_hadronicOverEm;
	out_Photon1_r9 = Photon1_r9;
	out_Photon1_r9_5x5 = Photon1_r9_5x5;
	out_Photon1_sigmaIetaIeta = Photon1_sigmaIetaIeta;
	out_Photon1_sigmaIetaIeta5x5 = Photon1_sigmaIetaIeta5x5;
	out_Photon1_sigmaEtaEta = Photon1_sigmaEtaEta;
	out_Photon1_sigmaIphiIphi = Photon1_sigmaIphiIphi;
	out_Photon1_sigmaIphiIphi5x5 = Photon1_sigmaIphiIphi5x5;
	out_Photon1_sigmaIetaIphi = Photon1_sigmaIetaIphi;
	out_Photon1_sigmaIetaIphi5x5 = Photon1_sigmaIetaIphi5x5;
	out_Photon1_maxEnergyXtal = Photon1_maxEnergyXtal;
	out_Photon1_iEta = Photon1_iEta;
	out_Photon1_iPhi = Photon1_iPhi;
//	out_Photon1_alphaHighPtID = Photon1_alphaHighPtID;
//	out_Photon1_kappaHighPtID = Photon1_kappaHighPtID;
//	out_Photon1_phoEAHighPtID = Photon1_phoEAHighPtID;
//	out_Photon1_chEAegmID = Photon1_chEAegmID;
//	out_Photon1_nhEAegmID = Photon1_nhEAegmID;
//	out_Photon1_phoEAegmID = Photon1_phoEAegmID;
	out_Photon1_passEGMLooseID = Photon1_passEGMLooseID;
	out_Photon1_passEGMMediumID = Photon1_passEGMMediumID;
	out_Photon1_passEGMTightID = Photon1_passEGMTightID;
	out_Photon1_isEB = Photon1_isEB;
	out_Photon1_isEBEtaGap = Photon1_isEBEtaGap;
	out_Photon1_isEBPhiGap = Photon1_isEBPhiGap;
	out_Photon1_isEBEEGap = Photon1_isEBEEGap;
	out_Photon1_passElectronVeto = Photon1_passElectronVeto;
	out_Photon1_passHTowOverE = Photon1_passHTowOverE;
	out_Photon1_passChIso = Photon1_passChIso;
	out_Photon1_passCorPhoIso = Photon1_passCorPhoIso;
	out_Photon1_passSieie = Photon1_passSieie;
	out_Photon1_passHighPtID = Photon1_passHighPtID;
//	out_Photon1_passChIsoDenom = Photon1_passChIsoDenom;
//	out_Photon1_passCorPhoIsoDenom = Photon1_passCorPhoIsoDenom;
//	out_Photon1_isFakeable = Photon1_isFakeable;
//	out_Photon1_isNumeratorObjCand = Photon1_isNumeratorObjCand;
//	out_Photon1_isDenominatorObj = Photon1_isDenominatorObj;
	out_Photon1_isSaturated = Photon1_isSaturated;
	out_Event_npv_true = Event_npv_true;
	out_nPV = nPV;

	
	if (isGood && GenPhoton1_pt > 100.) {
		if (fabs(Photon1_scEta) < 1.4442) {
			Double_t p1eb_resval = ((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt);
			Photon1_EB_res_hist->Fill((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt);

			// if good section
        	if (fabs(p1eb_resval) <= 2*p1eb_fit2sigma) {
            	p1eb_fit2_good_count = p1eb_fit2_good_count + 1;
            	p1eb_fit2_classify = 0;
				p1eb_bad_classify = 0;
            	Photon1_EB_fit2_good_res_hist->Fill((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt);
            	Photon1_EB_fit2_good_eta_phi->Fill(Photon1_scEta,Photon1_phi);
		 	}// end if good

			// if bad section
        	if (fabs(p1eb_resval) > 2*p1eb_fit2sigma) {
        		p1eb_fit2_bad_count = p1eb_fit2_bad_count + 1;
        		p1eb_fit2_classify = 1;
            	Photon1_EB_fit2_bad_res_hist->Fill((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt);
            	Photon1_EB_fit2_bad_eta_phi->Fill(Photon1_scEta,Photon1_phi);

				if (p1eb_resval > 2*p1eb_fit2sigma) {
            		p1eb_fit2_bad_high_count++;
					Photon1_EB_fit2_bad_high_res_hist->Fill((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt);
            		Photon1_EB_fit2_bad_high_eta_phi->Fill(Photon1_scEta,Photon1_phi);
            		p1eb_bad_classify = 0;
         		}
                                
         		if (p1eb_resval < -2*p1eb_fit2sigma) {
            		p1eb_fit2_bad_low_count++;
					Photon1_EB_fit2_bad_low_res_hist->Fill((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt);
            		Photon1_EB_fit2_bad_low_eta_phi->Fill(Photon1_scEta,Photon1_phi);
            		p1eb_bad_classify = 1;
        		}
			}// end if bad
		}// end if Photon1_scEta < 1.4442
	PhotonClassifyTree->Fill();
	}// end if isGood && GenPhoton1_pt
}// end Long64_t nentries

cout << "p1eb fit2 Good: " << p1eb_fit2_good_count << endl;
cout << "p1eb fit2 Bad: " << p1eb_fit2_bad_count << endl;
cout << "p1eb fit2 Bad Low: " << p1eb_fit2_bad_low_count << endl;
cout << "p1eb fit2 Bad High: " << p1eb_fit2_bad_high_count << endl;

Photon1_EB_fit2_good_res_hist->Write();
Photon1_EB_fit2_bad_res_hist->Write();
Photon1_EB_res_hist->Write();
Photon1_EB_fit2_bad_eta_phi->Write();
Photon1_EB_fit2_bad_low_eta_phi->Write();
Photon1_EB_fit2_bad_high_eta_phi->Write();
Photon1_EB_fit2_good_eta_phi->Write();

PhotonClassifyTree->Write();

file_out.ls();
file_out.Close();
}// end PhotonClassify Loop

