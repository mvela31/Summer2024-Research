#include "histogram.h"

void DrawSurf() {
    TFile *f = new TFile("OO_DS_10k_0.5GeV_Npart.root");
    TTree *PairTree = (TTree*)f->Get("PairTree"); // retreive the TTree from the file
    TTree *MixedPairTree = (TTree*)f->Get("MixedPairTree"); // retreive the TTree from the file


    const double del_phi_eta_bins = 55;
    const double del_Xmax = 3*M_PI/2;
    const double del_Ymax = 2.5;
    const double del_Xmin = -M_PI/2;
    const double del_Ymin = -2.5;

    TAxis* del_phi_axis = getLinAxis(del_Xmin, del_Xmax, del_phi_eta_bins, "#Delta#phi", "#Delta#phi");
    TAxis* del_eta_axis = getLinAxis(del_Ymin, del_Ymax, del_phi_eta_bins, "#Delta#eta", "#Delta#eta");

    TH2D* del_hist_same = makeTH2D("del_hist_same", del_phi_axis, del_eta_axis);
    TH2D* del_hist_mixed = makeTH2D("del_hist_mixed", del_phi_axis, del_eta_axis);


    Double_t del_phi, del_eta;
    Int_t Npart_same;
    PairTree->SetBranchAddress("del_phi", &del_phi);
    PairTree->SetBranchAddress("del_eta", &del_eta);
    PairTree->SetBranchAddress("Npart",   &Npart_same);

    Double_t mdel_phi, mdel_eta;
    Int_t Npart_mixed;
    MixedPairTree->SetBranchAddress("mdel_phi", &mdel_phi);
    MixedPairTree->SetBranchAddress("mdel_eta", &mdel_eta);
    MixedPairTree->SetBranchAddress("Npart",   &Npart_mixed);


    /* if (del_phi > 2 * M_PI) {
        del_phi -= 2*M_PI;
    } else if (del_phi < 2 * M_PI) {
        del_phi += 2*M_PI;
    }
    */
    
    Long64_t nentries1 = PairTree->GetEntries();
    for (Long64_t i=0;i<nentries1;i++) {
        PairTree->GetEntry(i);
        if (Npart_same < 34.3636 && Npart_same > 19.0909) del_hist_same->Fill(del_phi, del_eta);
    }
    Long64_t nentries2 = MixedPairTree->GetEntries();
    for (Long64_t i=0;i<nentries2;i++) {
        MixedPairTree->GetEntry(i);
        if (Npart_mixed < 34.3636 && Npart_mixed > 19.0909) del_hist_mixed->Fill(mdel_phi, mdel_eta);
    }
    
    TFile *file0 = new TFile("OO_Collisions/OO_DS_1k_0.5GeV_output_corr_020.root","recreate");
    file0->cd();
    del_hist_same->SetOption("SURF1");
    del_hist_mixed->SetOption("SURF1");
    del_hist_same->Scale(1/del_hist_same->GetEntries());
    del_hist_mixed->Scale(1/del_hist_mixed->GetEntries());

    TH2D* rp = (TH2D*)del_hist_same->Clone();
    rp->Divide(del_hist_mixed);
    rp->SetName("ratio_plot");
    rp->SetOption("SURF1");
    rp->GetZaxis()->SetTitle("C(#Delta#phi, #Delta#eta)");

    TProfile* projXZ = new TProfile("projectionXZ", 
                                    "Projection C(#Delta#phi,#Delta#eta) into#Delta#phi;#Delta#phi;C(#Delta#phi)",
                                    rp->GetNbinsX(), rp->GetXaxis()->GetXmin(), rp->GetXaxis()->GetXmax());

    // Fill the profile with the average Z values
    for (int ix = 1; ix <= rp->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= rp->GetNbinsY(); ++iy) {
            double zValue = rp->GetBinContent(ix, iy);
            double xValue = rp->GetXaxis()->GetBinCenter(ix);
            projXZ->Fill(xValue, zValue);
        }
    }

    

    del_hist_same->Write("del_hist_same");
    del_hist_mixed->Write("del_hist_mixed");
    rp->Write("ratio_plot");
    projXZ->Write("projectionXZ");

    TCanvas *c1 = new TCanvas("c1", "Projection of Z vs X", 800, 600);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptStat(0);
    projXZ->SetMarkerStyle(20);
    projXZ->SetMarkerColor(30);
    projXZ->SetMarkerSize(1.5);
    projXZ->Draw();
    c1->SaveAs("OO_Collisions/projectionXZ_10k.pdf");
    
    delete c1;
    file0->Close();
    f->Close();

}