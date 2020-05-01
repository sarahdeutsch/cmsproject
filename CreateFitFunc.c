{
    // this is the file that gives you the fit function while not plotting a ton of stuff!
gStyle->SetOptStat("ourmen");   //adding over and underflow bins
TFile file_in("PhotonResolutionLoopOutput.root");
file_in.ls();

TFile file_out("fit_function.root","RECREATE");

//TH1D *Photon1_EB_res_hist = (TH1D *) file_in.Get("Photon1_EB_res_hist");
//TH1D *Photon1_EB_pt_hist = (TH1D *) file_in.Get("Photon1_EB_pt_hist");
TH1D *Photon1_EB_res_02_hist = (TH1D *) file_in.Get("Photon1_EB_res_02_hist");
TH1D *Photon1_EB_pt_02_hist = (TH1D *) file_in.Get("Photon1_EB_pt_02_hist");

std::vector<TH1D*> Photon1_EB_res_vect_02;
std::vector<TH1D*> Photon1_EB_pt_vect;

TCanvas *p1_EB_res_02[15];
TCanvas *p1_EB_pt[15];
TCanvas *p1_EB_res_02_gaus[15];
TCanvas *p1_EB_res_02_gaus_log[15];

std::vector<Double_t> Photon1_EB_binedge{100,150,200,250,300,350,400,450,500,550,600,700,800,900,1000,2100};
int p1_eb_bins =15;

for (int i=0;i<p1_eb_bins;i++) {
        TString Photon1_EB_res_02_name = Form("Photon1_EB_res_02_bin%d",i);
	    TString Photon1_EB_pt_name = Form("Photon1_EB_pt_bin%d",i);

        TH1D *Photon1_EB_res_02 = (TH1D *) file_in.Get(Photon1_EB_res_02_name);
	    TH1D *Photon1_EB_pt = (TH1D *) file_in.Get(Photon1_EB_pt_name);

        Photon1_EB_res_vect_02.push_back(Photon1_EB_res_02);
	    Photon1_EB_pt_vect.push_back(Photon1_EB_pt);
}

for (int i=0;i<p1_eb_bins;i++) {
    //Photon1_EB_res_vect_02[i]->Write();
    p1_EB_res_02_gaus_log[i] = new TCanvas(Form("p1_EB_res_02_gaus_log%d",i));
        gPad->SetLogy();
        Photon1_EB_res_vect_02[i]->Fit("gaus");
        Photon1_EB_res_vect_02[i]->Draw();
        Double_t p1eb_XcutValue  = 2*Photon1_EB_res_vect_02[i]->GetFunction("gaus")->GetParameter(2);
        TLine *p1eb_cutObj = new TLine(-p1eb_XcutValue, Photon1_EB_res_vect_02[i]->GetMinimum(), -p1eb_XcutValue, Photon1_EB_res_vect_02[i]->GetMaximum());
        TLine *p1eb_cutObj2 = new TLine(p1eb_XcutValue, Photon1_EB_res_vect_02[i]->GetMinimum(), p1eb_XcutValue, Photon1_EB_res_vect_02[i]->GetMaximum());
        p1eb_cutObj->SetLineWidth(2);
        p1eb_cutObj2->SetLineWidth(2);
        p1eb_cutObj->SetLineStyle(9);
        p1eb_cutObj2->SetLineStyle(9);
        p1eb_cutObj->Draw("lsames");
        p1eb_cutObj2->Draw("lsames");
        p1_EB_res_02_gaus_log[i]->SaveAs(Form("/Users/sarahdeutsch/Box Sync/diphoton_GGJets_studies/final_project/PhotonResolutionLoopImages/p1eb/res02/p1_EB_res_02_gaus_log%d",i));
//	TFile file_out("fit_function.root","RECREATE");
	    Photon1_EB_res_vect_02[i]->Write();
        gPad->Update();
}

Double_t Photon1_EB_meanval[15];
for (int i=0;i<p1_eb_bins;i++) {
        Photon1_EB_meanval[i] = Photon1_EB_pt_vect[i]->GetMean();
}

Double_t Photon1_EB_gaus_sigma_02[15], Photon1_EB_gaus_sigma_error_02[15];
Double_t Photon1_EB_gaus_mean_02[15], Photon1_EB_gaus_mean_error_02[15];
for (int i=0;i<p1_eb_bins;i++) {
        Photon1_EB_gaus_sigma_02[i] = Photon1_EB_res_vect_02[i]->GetFunction("gaus")->GetParameter(2);
	    Photon1_EB_gaus_sigma_error_02[i] = Photon1_EB_res_vect_02[i]->GetFunction("gaus")->GetParError(2);
        Photon1_EB_gaus_mean_02[i] = Photon1_EB_res_vect_02[i]->GetFunction("gaus")->GetParameter(1);
        Photon1_EB_gaus_mean_error_02[i] = Photon1_EB_res_vect_02[i]->GetFunction("gaus")->GetParError(1);
}

TF1 *p1eb_myfit = new TF1("p1eb_myfit","sqrt(sq([0]/sqrt(x))+sq([1]/x)+sq([2]))",100,2100);
p1eb_myfit->SetParameters(0,0,0);
p1eb_myfit->SetParNames("a","b","c");
TGraphErrors *Photon1_EB_meanpt_sigma_error_fit_02 = new TGraphErrors(p1_eb_bins, Photon1_EB_meanval, Photon1_EB_gaus_sigma_02, 0, Photon1_EB_gaus_sigma_error_02);
p1_EB_meanpt_sigma_error_fit_02 = new TCanvas("p1_EB_meanpt_sigma_error_fit_02");
Photon1_EB_meanpt_sigma_error_fit_02->Fit("p1eb_myfit");
Photon1_EB_meanpt_sigma_error_fit_02->Draw("AC*");
Photon1_EB_meanpt_sigma_error_fit_02->SetMarkerStyle(7);
Photon1_EB_meanpt_sigma_error_fit_02->SetTitle("Photon 1 Barrel sigma vs Mean pt with error 0.2 Fit");
Photon1_EB_meanpt_sigma_error_fit_02->GetXaxis()->SetTitle("Mean pt");
Photon1_EB_meanpt_sigma_error_fit_02->GetXaxis()->CenterTitle();
Photon1_EB_meanpt_sigma_error_fit_02->GetYaxis()->SetTitle("gaus sigma");
Photon1_EB_meanpt_sigma_error_fit_02->GetYaxis()->CenterTitle();
gStyle->SetOptFit();
p1_EB_meanpt_sigma_error_fit_02->SaveAs("/Users/sarahdeutsch/Box Sync/diphoton_GGJets_studies/final_project/PhotonResolutionLoopImages/p1eb/res02/p1_EB_meanpt_sigma_error_fit_02");
p1eb_myfit->Write();
gPad->Update();

TF1 *p1eb_fit2 = new TF1("p1eb_fit2","([0]+ ([1]/(TMath::Power(x,[2])))+([3]/(TMath::Power(x,[4]))))",100,2100);
p1eb_fit2->SetParameters(p1eb_myfit->GetParameter(2),p1eb_myfit->GetParameter(1),1,p1eb_myfit->GetParameter(0),0.5);
p1eb_fit2->SetParNames("a","b","c","d","e");
TGraphErrors *Photon1_EB_meanpt_sigma_error_fit2_02 = new TGraphErrors(p1_eb_bins, Photon1_EB_meanval, Photon1_EB_gaus_sigma_02, 0, Photon1_EB_gaus_sigma_error_02);
p1_EB_meanpt_sigma_error_fit2_02 = new TCanvas("p1_EB_meanpt_sigma_error_fit2_02");
Photon1_EB_meanpt_sigma_error_fit2_02->Fit("p1eb_fit2");

Photon1_EB_meanpt_sigma_error_fit2_02->Draw("AC*");
Photon1_EB_meanpt_sigma_error_fit2_02->SetMarkerStyle(7);
Photon1_EB_meanpt_sigma_error_fit2_02->SetTitle("Photon 1 Barrel sigma vs Mean pt with error  0.2 Fit2");
Photon1_EB_meanpt_sigma_error_fit2_02->GetXaxis()->SetTitle("Mean pt");
Photon1_EB_meanpt_sigma_error_fit2_02->GetXaxis()->CenterTitle();
Photon1_EB_meanpt_sigma_error_fit2_02->GetYaxis()->SetTitle("gaus sigma");
Photon1_EB_meanpt_sigma_error_fit2_02->GetYaxis()->CenterTitle();
gStyle->SetOptFit();
p1_EB_meanpt_sigma_error_fit2_02->SaveAs("/Users/sarahdeutsch/Box Sync/diphoton_GGJets_studies/final_project/PhotonResolutionLoopImages/p1eb/res02/p1_EB_meanpt_sigma_error_fit2_02");

p1eb_fit2->Write();
gPad->Update();
}





