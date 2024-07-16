#include "histogram.h"

//=====================
void Tree_reader_Angantyr(){
  //bins and ranges
const int nCollBins = 55;
const double nCollxMin = 0;
const double nCollxMax = 35;

const int NpartBins = 55;
const double NpartxMin = 0;
const double NpartxMax = 35;

const int Npart_2DBins = 100;
const double Npart_2DxMin = 0;
const double Npart_2DxMax = 35;

const int epsilon2Bins = 100;
const double epsilon2xMin = 0;
const double epsilon2xMax = 1.1;

const int epsilon3Bins = 100;
const double epsilon3xMin = 0;
const double epsilon3xMax = 1.1;

const int epsilon2_1DBins = 55;
const double epsilon2_1DxMin = 0;
const double epsilon2_1DxMax = 1.1;

const int epsilon3_1DBins = 55;
const double epsilon3_1DxMin = 0;
const double epsilon3_1DxMax = 1.1;


//axes
TAxis* nCollAxis = getLinAxis(nCollxMin, nCollxMax, nCollBins, "N_{coll}", "N_{coll}");
TAxis* NpartAxis = getLinAxis(NpartxMin, NpartxMax, NpartBins, "N_{part}", "N_{part}");
TAxis* Npart_2DAxis = getLinAxis(Npart_2DxMin, Npart_2DxMax, Npart_2DBins, "N_{part}", "N_{part}");
TAxis* epsilon2Axis = getLinAxis(epsilon2xMin, epsilon2xMax, epsilon2Bins, "#varepsilon_{2}", "#varepsilon_{2}");
TAxis* epsilon3Axis = getLinAxis(epsilon3xMin, epsilon3xMax, epsilon3Bins, "#varepsilon_{3}", "#varepsilon_{3}");
TAxis* epsilon2_1DAxis = getLinAxis(epsilon2_1DxMin, epsilon2_1DxMax, epsilon2_1DBins, "#varepsilon_{2}", "#varepsilon_{2}");
TAxis* epsilon3_1DAxis = getLinAxis(epsilon3_1DxMin, epsilon3_1DxMax, epsilon3_1DBins, "#varepsilon_{3}", "#varepsilon_{3}");

//declaring histograms
TH1D* nColl_hist = makeTH1D("nColl_hist", nCollAxis);
TH1D* Npart_hist = makeTH1D("Npart_hist", NpartAxis);
TH2D* epsilon2_hist = makeTH2D("epsilon2_hist", Npart_2DAxis, epsilon2Axis);
TH2D* epsilon3_hist = makeTH2D("epsilon3_hist", Npart_2DAxis, epsilon3Axis);
TH1D* epsilon2_1D_hist = makeTH1D("epsilon2_1D_hist", epsilon2_1DAxis);
TH1D* epsilon3_1D_hist = makeTH1D("epsilon3_1D_hist", epsilon3_1DAxis); 

Int_t Npart, nCollTot;
Double_t epsilon2;
Double_t epsilon3;
DOuble_t v_2, v_3;

TFile *f = new TFile("OO_DS_10k.root"); //source file

TTree *theTree = (TTree*)f->Get("MBTree"); // retreive the TTree from the file

//Set the branch adress
theTree->SetBranchAddress("Npart", &Npart);
theTree->SetBranchAddress("nCollTot", &nCollTot);
theTree->SetBranchAddress("epsilon_2", &epsilon2);
theTree->SetBranchAddress("epsilon_3", &epsilon3);
theTree->SetBranchAddress("v_2", &v_2);
theTree->SetBranchAddress("v_3", &v_3);

//==============================================================
    // Start the event loop
    //==============================================================
    Long64_t nentries = theTree->GetEntries();
    for (Long64_t i=0;i<nentries;i++) {
      theTree->GetEntry(i);
      nColl_hist->Fill(nCollTot);
      Npart_hist->Fill(Npart);
      epsilon2_hist->Fill(Npart,epsilon2);
      epsilon3_hist->Fill(Npart,epsilon3);
      epsilon2_1D_hist->Fill(epsilon2);
      epsilon3_1D_hist->Fill(epsilon3);
    }

//Write histograms to a ROOT file
TFile *file0 = new TFile("OO_Collisions/OO_DS_10k_HR_output.root","recreate");
file0->cd();

//write histograms
nColl_hist->Write("nColl_hist");
Npart_hist->Write("Npart_hist");
epsilon2_hist->Write("epsilon2_hist");
epsilon3_hist->Write("epsilon3_hist");
epsilon2_1D_hist->Write("epsilon2_1D_hist");
epsilon3_1D_hist->Write("epsilon3_1D_hist");
file0->Close();
f->Close();


  //==============================================================
  // Draw
  //==============================================================
// create a style for your histograms, draw them, and save them

  plot1DHist(nColl_hist, "DS0.1_nColl_hist", "N_{coll}", "N_{coll}", 0, nCollxMin, nCollxMax, -1, -1);
  plot1DHist(Npart_hist, "DS0.1_Npart_hist", "N_{part}", "N_{part}", 0, NpartxMin, NpartxMax, -1, -1);
  plot2DHist(epsilon2_hist, "DS0.1_epsilon2_hist", "N_{part}", "#varepsilon_{2}", 0, NpartxMin, NpartxMax, epsilon2xMin, epsilon2xMax);
  plot2DHist(epsilon3_hist, "DS0.1_epsilon3_hist", "N_{part}", "#varepsilon_{3}", 0, NpartxMin, NpartxMax, epsilon3xMin, epsilon3xMax);
  plot1DHist(epsilon2_1D_hist, "DS0.1_epsilon2_1D_hist", "#varepsilon_{2}", "#varepsilon_{2}", 0, epsilon2_1DxMin, epsilon2_1DxMax, -1, -1);
  plot1DHist(epsilon3_1D_hist, "DS0.1_epsilon3_1D_hist", "#varepsilon_{3}", "#varepsilon_{3}", 0, epsilon3_1DxMin, epsilon3_1DxMax, -1, -1);
 
/*
  //==============================================================
  // Draw
  //==============================================================
// create a style for your histograms, draw them, and save them


  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);

   auto* c1    = new TCanvas("c1", "c1", 2400, 1800);
   c1->cd();
   nColl_hist->SetMarkerStyle(20);
   nColl_hist->SetMarkerColor(9);
   nColl_hist->SetTitle("Angantyr - n_{Coll}");
   nColl_hist->GetXaxis()->SetTitle("n_{Coll}");
   nColl_hist->GetYaxis()->SetTitle("#frac{dN}{dn_{Coll}}");
   nColl_hist->Draw("P");
   c1->Draw();
   c1->SaveAs("Summer/A_nColl.pdf");
   c1->Close();

   auto* c2    = new TCanvas("c2", "c2", 2400, 1800);
   c2->cd();
   Npart_hist->SetMarkerStyle(20);
   Npart_hist->SetMarkerColor(46);
   Npart_hist->SetTitle("Angantyr - N_{part}");
   Npart_hist->GetXaxis()->SetTitle("N_{part}");
   Npart_hist->GetYaxis()->SetTitle("#frac{dN}{dN_{part}}");
   Npart_hist->Draw("P");
   c2->Draw();
   c2->SaveAs("Summer/A_Npart.pdf");
   c2->Close();

   auto* c3    = new TCanvas("c3", "c3", 2400, 1800);
   c3->cd();
   gStyle->SetPalette(112);
   epsilon2_hist->SetTitle("Angantyr - #varepsilon_{2} vs N_{part}");
   epsilon2_hist->GetXaxis()->SetTitle("N_{part}");
   epsilon2_hist->GetYaxis()->SetTitle("#varepsilon_{2}");
   epsilon2_hist->SetTitleSize(0.05, "Y");
   epsilon2_hist->Draw("colz");
   c3->Draw();
   c3->SaveAs("Summer/A_epsilon2.pdf");
   c3->Close();

   auto* c4    = new TCanvas("c4", "c4", 2400, 1800);
   c4->cd();
   gStyle->SetPalette(103);
   epsilon3_hist->SetTitle("Angantyr - #varepsilon_{3} vs N_{part}");
   epsilon3_hist->GetXaxis()->SetTitle("N_{part}");
   epsilon3_hist->GetYaxis()->SetTitle("#varepsilon_{3}");
   epsilon3_hist->SetTitleSize(0.05, "Y");
   epsilon3_hist->Draw("colz");
   c4->Draw();
   c4->SaveAs("Summer/A_epsilon3.pdf");
   c4->Close();

   auto* c5    = new TCanvas("c2", "c2", 2400, 1800);
   c5->cd();
   epsilon2_1D->SetMarkerStyle(20);
   epsilon2_1D->SetMarkerColor(6);
   epsilon2_1D->SetTitle("Angantyr - epsilon2 1D");
   epsilon2_1D->GetXaxis()->SetTitle("#varepsilon_{2}");
   epsilon3_hist->SetTitleSize(0.05, "X");
   epsilon2_1D->GetYaxis()->SetTitle("#frac{dN}{d#varepsilon_{2}}");
   epsilon2_1D->Draw("P");
   c5->Draw();
   c5->SaveAs("Summer/A_epsilon2_1D.pdf");
   c5->Close();

   auto* c6    = new TCanvas("c2", "c2", 2400, 1800);
   c6->cd();
   epsilon3_1D->SetMarkerStyle(20);
   epsilon3_1D->SetMarkerColor(8);
   epsilon3_1D->SetTitle("Angantyr - epsilon3 1D");
   epsilon3_1D->GetXaxis()->SetTitle("#varepsilon_{3}");
   epsilon3_hist->SetTitleSize(0.05, "X");
   epsilon3_1D->GetYaxis()->SetTitle("#frac{dN}{d#varepsilon_{3}}");
   epsilon3_1D->Draw("P");
   c6->Draw();
   c6->SaveAs("Summer/A_epsilon3_1D.pdf");
   c6->Close();
   */
}







