#include <iomanip>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include "math.h"
#include "string.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TUnixSystem.h"

using namespace std;

inline TAxis* getLinAxis(const Float_t min, const Float_t max, const Int_t nBins, std::string axisName, std::string axisTitle)
  {
    double xBins[nBins+1];
    xBins[0] = min; xBins[nBins] = max;
    Float_t interval = (xBins[nBins] - xBins[0])/nBins;
    for(Int_t iter = 1; iter < nBins; iter++){
      xBins[iter] = xBins[0] + iter*interval;
    }
    TAxis* myAx = new TAxis(nBins, &xBins[0]);
    myAx->SetNameTitle(axisName.c_str(), axisTitle.c_str());
    return myAx;
  }
inline TH1D* makeTH1D(std::string name, TAxis* xaxis, bool sumw = true)
{
    std::stringstream ss;
    ss << ";" << xaxis->GetTitle() << ";Entries";
    TH1D* myH = new TH1D(name.c_str(),ss.str().c_str(),xaxis->GetXbins()->GetSize()-1,xaxis->GetXbins()->GetArray());
    if(sumw) myH->Sumw2();
    return myH;
}
inline TH2D* makeTH2D(std::string name, TAxis* xaxis, TAxis* yaxis)
  {
    std::stringstream ss;
    ss << ";" << xaxis->GetTitle()
       << ";" << yaxis->GetTitle()
       << ";Entries";
    return new TH2D(name.c_str(),ss.str().c_str(),
        xaxis->GetXbins()->GetSize()-1,xaxis->GetXbins()->GetArray(),
        yaxis->GetXbins()->GetSize()-1,yaxis->GetXbins()->GetArray());
}

inline void plot1DHist ( TH1D* h1, std::string hname, std::string xtitle, std::string kinem, int opt,
                        float xmin, float xmax, float ymin, float ymax){
  auto* c1    = new TCanvas("c1", "c1", 2400, 1800);
  c1->cd();
  // auto* leg1  = new TLegend(0.67, 0.65, 0.9, 0.89);
  // leg1->SetNColumns(1);
  // leg1->SetTextSize(0.03);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  gPad->SetBottomMargin(0.15);
  printf("nbins[%d] \n",h1->GetXaxis()->GetNbins());
  // - - - - - - - - -
  for(int bin=1; bin<h1->GetXaxis()->GetNbins()+1;bin++){
    float value = h1->GetBinContent(bin);
    float error = h1->GetBinError(bin);
    float binL  = h1->GetXaxis()->GetBinLowEdge(bin);
    float binH  = h1->GetXaxis()->GetBinUpEdge(bin);
    value /= (binH-binL);
    error /= (binH-binL);
    h1->SetBinContent(bin,value);
    h1->SetBinError(bin,error);
  }
  if(opt==1) {
    gPad->SetLogx();
    gPad->SetLogy();
    h1->GetXaxis()->SetMoreLogLabels();
  }
  if (ymax > 0) h1->GetYaxis()->SetRangeUser(ymin, ymax);
  h1->GetXaxis()->SetRangeUser(xmin, xmax);
  h1->SetMarkerStyle(20);
  h1->SetMarkerColor(30);
  h1->SetMarkerSize(1.5);
  h1->GetXaxis()->SetTitle(xtitle.c_str());
  h1->GetXaxis()->SetTitleOffset(1.5);
  h1->GetYaxis()->SetTitle(Form("dN/d%s",kinem.c_str()));
  h1->Draw();
  c1->Draw();
  c1->SaveAs(Form("OO_Collisions/%s.pdf",hname.c_str()));
  c1->Close();
}

inline void plot2DHist( TH2D* h2, std::string hname, std::string xtitle, std::string ytitle, int opt,
                        float xmin, float xmax, float ymin, float ymax){ //ARIC IS HERE
  auto* c1    = new TCanvas("c1", "c1", 2400, 1800);
  c1->cd();
  auto* leg1  = new TLegend(0.67, 0.65, 0.9, 0.89);
  leg1->SetNColumns(1);
  leg1->SetTextSize(0.03);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.1);

  // - - - - - - - - -
  int binsx = h2->GetNbinsX();
  int binsy = h2->GetNbinsY();
  int totbin = binsx*binsy;
  for(int bin=0; bin<totbin;bin++){
    float value = h2->GetBinContent(bin);
    float binL = h2->GetXaxis()->GetBinLowEdge(bin);
    float binH = h2->GetXaxis()->GetBinUpEdge(bin);
    value /= (binH-binL);
    h2->SetBinContent(bin,value);
  }
  if(opt==1) {
    gPad->SetLogx();
    gPad->SetLogy();
    h2->GetXaxis()->SetMoreLogLabels();
    h2->GetYaxis()->SetMoreLogLabels();
  }
  h2->SetAxisRange(ymin, ymax, "Y");
  h2->SetAxisRange(xmin, xmax, "X");
  h2->GetXaxis()->SetTitle(xtitle.c_str());
  h2->GetYaxis()->SetTitle(ytitle.c_str());
  h2->GetXaxis()->SetTitleOffset(1.6);
  h2->GetYaxis()->SetTitleOffset(1.1);
  // if(opt==2) leg1->AddEntry(h2, “2D”  , “h”);
  h2->Draw("colz");
  if(opt==2){
    TProfile* profx = h2->ProfileX("profx");
    TProfile* profy = h2->ProfileY("profy");
    profx->SetLineColor(kRed);
    profx->SetLineStyle(2);
    profx->SetLineWidth(2);
    profy->SetLineColor(kCyan);
    profy->SetLineStyle(2);
    profy->SetLineWidth(2);
    profx->Draw("same");
    profy->Draw("same");
    leg1->AddEntry(profx, "profX"  , "l");
    leg1->AddEntry(profy, "profY"  , "l");
    c1->Update();
  }
  gPad->Modified();
  gPad->Update();
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->Draw();
  c1->Draw();
  c1->SaveAs(Form("OO_Collisions/%s.pdf",hname.c_str()));
  c1->Close();
}

inline void plotVecHist(  std::vector<TH1D*> v_h1, std::vector<double> cent_lim,
                          std::string hname, int opt, float ymin, float ymax){                       
  auto* c1    = new TCanvas("c1", "c1", 2400, 1800);
  c1->cd();
  auto* leg1  = new TLegend(0.6, 0.65, 0.9, 0.89);
  leg1->SetNColumns(2);
  leg1->SetTextSize(0.03);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  // - - - - - - - - -
  int centBins = cent_lim.size()-1;
  TH1D* h_sum = (TH1D*)v_h1[0]->Clone("h_sum");
  h_sum->Reset();
  // - - - - - - - - -
  for (int cent = 0; cent < centBins; ++cent) {
    float norm = cent_lim[cent+1]-cent_lim[cent];
    v_h1[cent]->Scale(1.0 / (0 - norm));

    h_sum->Add(v_h1[cent]);
    // v_h1[cent]->SetAxisRange(0, 5, “Y”);

    if(ymax>0) v_h1[cent]->SetAxisRange(ymin, ymax, "Y");

    v_h1[cent]->SetMarkerStyle(20+cent);
    v_h1[cent]->SetMarkerColor(kViolet-cent);
    v_h1[cent]->SetMarkerSize(1.5);
    std::ostringstream oss;
    oss << "#varepsilon_{" << opt << "}";
    std::string title = oss.str();
    std::string title2 = "dN/d" + title;
    v_h1[cent]->GetXaxis()->SetTitle(title.c_str());
    v_h1[cent]->GetYaxis()->SetTitle(title2.c_str());
    std::string centLabel = Form("%.0f - %.0f %%",cent_lim[cent+1]*100,cent_lim[cent]*100);
    leg1->AddEntry(v_h1[cent], Form("%s",centLabel.c_str()) , "P");
    if(cent==0) v_h1[cent]->Draw();
    else        v_h1[cent]->Draw("same");
    if(cent==centBins-1){
      leg1->AddEntry(h_sum, Form("Inclusive [0-90%%]") , "lep");
      h_sum->Scale(1.0 / 9.0);
      h_sum->SetMarkerStyle(8);
      h_sum->SetMarkerColor(kBlack);
      h_sum->SetMarkerSize(1.5);
      h_sum->GetXaxis()->SetNdivisions(506);
      h_sum->Draw("same");
  }
  } //end cent loop
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->Draw();
  c1->Draw();
  c1->SaveAs(Form("OO_Collisions/%s.pdf",hname.c_str()));
  c1->Close();
}

// inline void plotCompare(  std::vector<TH1D*> v_h1a, std::vector<TH1D*> v_h1b,
//                           std::vector<double> v_weight_a, std::vector<double> v_weight_b,
//                           std::vector<double> cent_lim, std::string hname, int opt, float ymin, float ymax){

inline void plotCompare(  std::vector<TH1D*> v_h1a, std::vector<TH1D*> v_h1b,
                          std::vector<double> cent_lim, std::string hname, int opt, float ymin, float ymax){
  auto* c1    = new TCanvas("c1", "c1", 2400, 1800);
  c1->cd();
  auto* leg1  = new TLegend(0.67, 0.65, 0.9, 0.89);
  leg1->SetNColumns(2);
  leg1->SetTextSize(0.03);
  auto* leg2  = new TLegend(0.4, 0.65, 0.65, 0.89);
  leg2->SetNColumns(2);
  leg2->SetTextSize(0.03);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  // - - - - - - - - -
  int centBins = cent_lim.size()-1;
  // - - - - - - - - -
  for (int cent = 0; cent < centBins; ++cent) {
    float norm = cent_lim[cent+1]-cent_lim[cent];
    v_h1a[cent]->Scale(1.0 / (v_h1a[cent]->GetBinWidth(1)));
    v_h1b[cent]->Scale(1.0 / (v_h1b[cent]->GetBinWidth(1)));
    if(ymax>0) v_h1a[cent]->SetAxisRange(ymin, ymax, "Y");
    v_h1a[cent]->SetMarkerStyle(20);
    v_h1b[cent]->SetMarkerStyle(22);
    v_h1a[cent]->SetMarkerColor(kViolet-cent);
    v_h1b[cent]->SetMarkerColor(kOrange+cent);
    v_h1a[cent]->SetMarkerSize(1.5);
    v_h1b[cent]->SetMarkerSize(1.5);
    std::ostringstream oss;
    oss << "#varepsilon_{" << opt << "}";
    std::string title = oss.str();
    std::string title2 = "dN/d" + title;
    v_h1a[cent]->GetXaxis()->SetTitle(title.c_str());
    v_h1a[cent]->GetYaxis()->SetTitle(title2.c_str());
    std::string centLabel = Form("%.0f - %.0f %%",cent_lim[cent]*100,cent_lim[cent+1]*100);
    leg1->AddEntry(v_h1a[cent], Form("%s",centLabel.c_str()) , "p");
    leg2->AddEntry(v_h1b[cent], Form("%s",centLabel.c_str()) , "p");
    if(cent==0) { v_h1a[cent]->Draw(""),      v_h1b[cent]->Draw("same"); }
    else        { v_h1a[cent]->Draw("same"),  v_h1b[cent]->Draw("same"); }
  } //end cent loop
  //gPad->SetLogy();
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->Draw();
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->Draw();
  c1->Draw();
  c1->SaveAs(Form("OO_Collisions/%s.pdf",hname.c_str()));
  c1->Close();
}

inline void plotCompare3Gens(  std::vector<TH1D*> v_h1,
                            std::string hname, int opt, float ymin, float ymax) {          
  auto* c1    = new TCanvas("c1", "c1", 2400, 1800);
  c1->cd();
  auto* leg1  = new TLegend(0.3, 0.8, 0.7, 0.9);
  leg1->SetNColumns(3);
  leg1->SetTextSize(0.03); 
  //int result = v_h1[0]->GetMaximum();
  //if (v_h1[1]->GetMaximum() > result) result = v_h1[1]->GetMaximum();
  //if (v_h1[2]->GetMaximum() > result) result = v_h1[2]->GetMaximum();
  //if(ymax>0) v_h1[0]->SetAxisRange(0, result+2000, "Y");
  //if(ymax>0) v_h1[0]->SetAxisRange(0, ymax, "Y");
  std::vector<std::string> gens = {"Angantyr", "Glauber", "Trento"};
  for (int i = 0; i < 3; i++) {
    //v_h1[i]->Rebin(2);
    v_h1[i]->Scale(1.0 / (v_h1[i]->GetBinWidth(1)));
    if(ymax>0) v_h1[i]->SetAxisRange(ymin, ymax, "Y");
    v_h1[i]->SetStats(0);

    v_h1[i]->SetMarkerStyle(20+i);
    v_h1[i]->SetMarkerColor(kPink-(i*20));
    v_h1[i]->SetMarkerSize(2);
    std::ostringstream oss;
    oss << "#varepsilon_{" << opt << "}";
    std::string title = oss.str();
    std::string title2 = "dN/d" + title;
    v_h1[i]->GetXaxis()->SetTitle(title.c_str());
    v_h1[i]->GetYaxis()->SetTitle(title2.c_str());
    leg1->AddEntry(v_h1[i], gens[i].c_str(), "P");
    if(i==0) v_h1[i]->Draw();
    else        v_h1[i]->Draw("same");
  }        
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->Draw();
  c1->Draw();
  c1->SaveAs(Form("OO_Collisions/%s.pdf",hname.c_str()));
  c1->Close();            
}