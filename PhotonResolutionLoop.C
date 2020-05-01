// this is the file that gives PhotonResolutionLoopOutput.root
#define PhotonResolutionLoop_cxx
#include "PhotonResolutionLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void PhotonResolutionLoop::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PhotonResolutionLoop.C
//      root> PhotonResolutionLoop t
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

        std::vector<Double_t> Photon1_EB_binedge{100,150,200,250,300,350,400,450,500,550,600,700,800,900,1000,2100};
        std::vector<TH1D*> Photon1_EB_res_vect_02;
	std::vector<TH1D*> Photon1_EB_pt_vect;

                // 1d hist of all p1eb resolutions
                TH1D *Photon1_EB_res_hist = new TH1D("Photon1_EB_res_hist","Photon 1 Barrel Resolution",100,-1,1);
                Photon1_EB_res_hist->GetXaxis()->SetTitle("(Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt");
                Photon1_EB_res_hist->GetXaxis()->CenterTitle();
                Photon1_EB_res_hist->GetYaxis()->SetTitle("Number of events");
                Photon1_EB_res_hist->GetYaxis()->CenterTitle();

                // 1d hist of all p1eb pts
                TH1D *Photon1_EB_pt_hist = new TH1D("Photon1_EB_pt_hist","Photon 1 Barrel Pt",100,100,2100);
                Photon1_EB_pt_hist->GetXaxis()->SetTitle("GenPhoton1_pt");
                Photon1_EB_pt_hist->GetXaxis()->CenterTitle();
                Photon1_EB_pt_hist->GetYaxis()->SetTitle("Number of events");
                Photon1_EB_pt_hist->GetYaxis()->CenterTitle();

                // 1d hist of just p1eb resolutions within 0.2 resolution
                TH1D *Photon1_EB_res_02_hist = new TH1D("Photon1_EB_res_02_hist","Photon 1 Barrel 02 Resolution",100,-0.2,0.2);
                Photon1_EB_res_02_hist->GetXaxis()->SetTitle("(Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt");
                Photon1_EB_res_02_hist->GetXaxis()->CenterTitle();
                Photon1_EB_res_02_hist->GetYaxis()->SetTitle("Number of events");
                Photon1_EB_res_02_hist->GetYaxis()->CenterTitle();

                // 1d hist of just p1eb pts within 0.2 resolution
                TH1D *Photon1_EB_pt_02_hist = new TH1D("Photon1_EB_pt_02_hist","Photon 1 Barrel 02 Pt",100,100,2100);
                Photon1_EB_pt_02_hist->GetXaxis()->SetTitle("GenPhoton1_pt");
                Photon1_EB_pt_02_hist->GetXaxis()->CenterTitle();
                Photon1_EB_pt_02_hist->GetYaxis()->SetTitle("Number of events");
                Photon1_EB_pt_02_hist->GetYaxis()->CenterTitle();

	for (int i=0;i<(Photon1_EB_binedge.size()-1);i++) { 
                TString Photon1_EB_res_02_name = Form("Photon1_EB_res_02_bin%d",i);
		TString Photon1_EB_pt_name = Form("Photon1_EB_pt_bin%d", i);

                TString p1eb_pt_title;
                Double_t ymin = Photon1_EB_binedge[i];
                Double_t ymax = Photon1_EB_binedge[i+1];
                p1eb_pt_title = Form("Photon1 EB pt bin%d ",i);
                p1eb_pt_title += ymin;
                p1eb_pt_title += " to ";
                p1eb_pt_title += ymax;

                TString p1eb_res_title;
                p1eb_res_title = Form("Photon1 EB res bin%d pt ",i);
                p1eb_res_title += ymin;
                p1eb_res_title += " to ";
                p1eb_res_title += ymax;

		TH1D *Photon1_EB_pt = new TH1D(Photon1_EB_pt_name,p1eb_pt_title,100,Photon1_EB_binedge[i],Photon1_EB_binedge[i+1]); 
                Photon1_EB_pt->GetXaxis()->SetTitle("GenPhoton1_pt");
                Photon1_EB_pt->GetXaxis()->CenterTitle();
                Photon1_EB_pt->GetYaxis()->SetTitle("Number of events");
                Photon1_EB_pt->GetYaxis()->CenterTitle();

                TH1D *Photon1_EB_res_02 = new TH1D(Photon1_EB_res_02_name,p1eb_res_title,100,-0.2,0.2);
                Photon1_EB_res_02->GetXaxis()->SetTitle("(Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt");
                Photon1_EB_res_02->GetXaxis()->CenterTitle();
                Photon1_EB_res_02->GetYaxis()->SetTitle("Number of events");
                Photon1_EB_res_02->GetYaxis()->CenterTitle();

		Photon1_EB_pt_vect.push_back(Photon1_EB_pt); 
                Photon1_EB_res_vect_02.push_back(Photon1_EB_res_02);

	} // end p1 eb for

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

    if (jentry%10000 == 0) cout << "Entry " << jentry << endl;

// Photon 1 EB binned histograms
if (isGood && GenPhoton1_pt > 100.) {
        // in barrel
	if (fabs(Photon1_scEta) < 1.4442) {
		int this_bin = -1;
		Photon1_EB_res_hist->Fill((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt);
		Photon1_EB_pt_hist->Fill(GenPhoton1_pt);
                // if resolution is between 0.2
		if (-0.2 < ((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt) && ((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt) < 0.2) {
			Photon1_EB_pt_02_hist->Fill(GenPhoton1_pt);
	                Photon1_EB_res_02_hist->Fill((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt);
		}
                // loop thru bins
		for (Int_t i=0; i<(Photon1_EB_binedge.size()-1); i++) {
			if (GenPhoton1_pt>= Photon1_EB_binedge[i] & GenPhoton1_pt< Photon1_EB_binedge[i+1]) {
				this_bin = i;
			} //end if
		} // end for

		if (this_bin >= 0) {
                        Photon1_EB_res_vect_02.at(this_bin)->Fill((Photon1_pt-GenPhoton1_pt)/GenPhoton1_pt);
			Photon1_EB_pt_vect.at(this_bin)->Fill(GenPhoton1_pt);
		} //end if
	} // end if barrel

} //end isgood
} //end of entry loop

TFile file_out("PhotonResolutionLoopOutput.root","RECREATE");

//Photon1_EB_res_hist->Write();
//Photon1_EB_pt_hist->Write();

//only withing range -0.2 to 0.2
Photon1_EB_res_02_hist->Write();
Photon1_EB_pt_02_hist->Write();

for (Int_t i=0; i<(Photon1_EB_binedge.size()-1); i++) {
        Photon1_EB_res_vect_02.at(i)->Write();
	Photon1_EB_pt_vect.at(i)->Write();
} // end p1 eb

file_out.ls();
file_out.Close();
}

