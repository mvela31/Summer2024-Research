// main423.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Leif Lonnblad <leif.lonnblad@fysik.lu.se>

// Keywords: heavy ions; charged multiplicity; centrality;
//           angantyr

// This test program will generate Pb-Pb collisions at
// sqrt(S_NN)=2.76TeV using the Angantyr model for Heavy Ion
// collisions. The analysis will divide the event in centrality
// classes using the same observable as was used for p-Pb in the ATLAS
// analysis in arXiv:1508.00848 [hep-ex] (see main422.cc). The
// centrality classes are same as in the ALICE analysis in
// arXiv:1012.1657 [nucl-ex] although the actual observable used is
// not the same. Histograms of multiplicity distributions are measured
// for each centrality percentile.

// Note that heavy ion collisions are computationally quite CPU
// intensive and generating a single event will take around a second
// on a reasonable desktop. To get reasonable statistics, this program
// will take a couple of hours to run.

#include "Pythia8/Pythia.h"
#include <cmath>

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TRandom2.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TTree.h"
#include "TVector2.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TApplication.h"
#include "TClass.h"
#include <vector>

using namespace Pythia8;

inline void getLogBins(const Float_t lower, const Float_t higher, const UInt_t nBins, std::vector<double>& bins)
{
    
    bins.assign(nBins+1,lower);
    Double_t scaleFactor=std::pow(higher/lower,1./static_cast<double>(nBins));
    for(UInt_t iter = 1; iter <= nBins; iter++){
        bins[iter]=scaleFactor*bins[iter-1];
        //std::cout << " Bins " << iter << " : " <<  bins[iter] << std::endl;
    }
    return;
}

inline TAxis* constructAxis_vD ( std::vector < double > bins, std::string name, std::string title){
      TAxis* theAxis=new TAxis(bins.size()-1,&(bins[0]));
      theAxis->SetNameTitle(name.c_str(), title.c_str());
      return theAxis;
}

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

//std::vector<TH1D*> getLogBins(int nbins, double xmin, double xmax, int nhistograms, double* pT_avgLog) {
inline std::vector<TH1D*> histoVectoMaker(std::string histoID, TAxis* xAxis, int nHistograms){

        std::vector<TH1D*> histograms;
        for (int i = 0; i < nHistograms; ++i) {
            histograms.push_back(new TH1D(Form("h%s_%d", histoID.c_str(), i), Form("%s - bin %d", histoID.c_str(), i), xAxis->GetXbins()->GetSize()-1,xAxis->GetXbins()->GetArray()));
        }
        return histograms;
}

int main() {

  Pythia pythia;

  // Setup the beams.
  pythia.readString("Beams:idA = 1000080160");
  pythia.readString("Beams:idB = 1000080160"); // The lead ions.
  pythia.readString("Beams:eCM = 5020.0");
  pythia.readString("Beams:frameType = 1");
  //pythia.readString("HadronLevel:Rescatter = on");

  //pythia.readString("SoftQCD:all = on");
  // Initialize the Angantyr model to fit the total and semi-inclusive
  // cross sections in Pythia within some tolerance.
  pythia.readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
  // These parameters are typically suitable for sqrt(S_NN)=5TeV
  pythia.readString("HeavyIon:SigFitDefPar = 2.15,17.24,0.33");
  // A simple genetic algorithm is run for 20 generations to fit the
  // parameters.
  pythia.readString("HeavyIon:SigFitNGen = 20");

  // There will be nine centrality bins based on the sum transverse
  // emergy in a rapidity interval between 3.2 and 4.9 obtained from
  // the borders from the generated transverse energy spectrum. The
  // default settings should give approximately the following:
  double genlim[] = {0.0, 215.4, 134.3, 96.0, 67.8,
                     48.2, 33.4, 22.9, 13.4};
  // If you change any parameters these should also be changed.

  // The upper edge of the correponding percentiles:
  double pclim[] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  // Book the pseudorapidity and multiplicity histograms and get
  // counters for number of events and sum of event weights:
  typedef map<double,int,std::greater<double> > MapIdx;
  MapIdx genetaidx;
  vector<TH1D*> etadist(9), lmult(9), hmult(9);
  string etaname("EtadistC"), mname("MultC");
  vector<double> gensumw(9, 0.0), gensumn(9, 0.0);
  TAxis* etaAxis = getLinAxis(-2.7, 2.7, 54, "eta", "eta");
  TAxis* hmultAxis = getLinAxis(-0.5, 2999.5, 75, "hmult", "hmult");
  TAxis* lmultAxis = getLinAxis(-0.5, 299.5, 75, "lmult", "lmult");


  //std::vector<TH1D*> etadist = histoVectoMaker(etaname, etaAxis, 9);
  //std::vector<TH1D*> hmult = histoVectoMaker("MultCH", hmultAxis, 9);
  //std::vector<TH1D*> hmult = histoVectoMaker("MultCL", lmultAxis, 9);


  for ( int i = 0; i < 9; ++i ) {
    genetaidx[genlim[i]] = i;
  }

  // Book histogram for the centrality measure.
  //TAxis* sumetAxis = getLinAxis(0, 400, 200, "SumETfwd", "SumETfwd");
  //TH1D* sumet = makeTH1D("SumETfwd", sumetAxis);

  // Also make a map of all weight to check the generated centrality
  // classes.
  multimap<double,double> gencent;

  // Book a histogram for the distribution of number of wounded
  // nucleons.
  //TAxis* woundedAxis = getLinAxis(-0.5, 417.5, 209, "Nwounded", "Nwounded");
  //TH1D* wounded = makeTH1D("Nwounded", woundedAxis);

  // Profile for average central multiplicity and number of wounded
  // nucleons as a function of centrality (with errors).
  vector<double> cmult(9, 0.0), cmult2(9, 0.0);
  vector<double> wound(9, 0.0), wound2(9, 0.0);

  // Sum up the weights of all generated events.
  std::vector < double > *eta = 0;
  std::vector < double > *pt = 0;
  double negFCalEt = 0;
  double impact_parameter;
  double weight = 0;
  double sumw = 0;
  int nCollTot, wn, wn_diff, nColl, nCollNDTot, Npart;
  double epsilon_2, epsilon_3, v_2, v_3, del_phi, del_eta;

  TTree* theTree = new TTree("MBTree", "MBTree");

  theTree->Branch("negFCalEt",        &negFCalEt); 
  theTree->Branch("eta",              &eta);
  theTree->Branch("pt",               &pt);
  theTree->Branch("weight",           &weight);
  theTree->Branch("impact_parameter", &impact_parameter);
  theTree->Branch("nCollTot",         &nCollTot);
  theTree->Branch("wn",               &wn);
  theTree->Branch("wn_diff",          &wn_diff);
  theTree->Branch("Npart",            &Npart);
  theTree->Branch("nColl",            &nColl);
  theTree->Branch("nCollNDTot",       &nCollNDTot);
  theTree->Branch("epsilon_2",        &epsilon_2);
  theTree->Branch("epsilon_3",        &epsilon_3);
  theTree->Branch("v_2",              &v_2);
  theTree->Branch("v_3",              &v_3);

  TTree* PairTree = new TTree("PairTree", "PairTree");
  PairTree->Branch("del_phi",          &del_phi);
  PairTree->Branch("del_eta",          &del_eta);


  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Loop over events.
  int nEvents = 1000;
  int nDiffp = 0;
  int nAbsp = 0;
  int nDiffmc = 0;
  int nAbsmc = 0;
  int nDiffc = 0;
  int nAbsc = 0;
  for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {
    if (iEvent % 100 == 0) {
      std::cout << "Event: " << iEvent << "/" << nEvents << std::endl;
    }
    if ( !pythia.next() ) continue;
    negFCalEt = 0;
        
    weight           = pythia.info.weight();
    impact_parameter = pythia.info.hiInfo->b();
    nCollNDTot       = pythia.info.hiInfo->nCollNDTot();
    nCollTot         = pythia.info.hiInfo->nCollTot();
    wn               = pythia.info.hiInfo->nAbsTarg() +
                        pythia.info.hiInfo->nDiffTarg() +
                        pythia.info.hiInfo->nAbsProj() +
                        pythia.info.hiInfo->nDiffProj() +
                        pythia.info.hiInfo->nElProj() +
                        pythia.info.hiInfo->nElTarg();
    wn_diff          = pythia.info.hiInfo->nDiffTarg();
    nColl            = wn;
    Npart            = pythia.info.hiInfo->nPartProj() +
                        pythia.info.hiInfo->nPartTarg();

    sumw += weight;
    /* if (Npart < 6) {
      if (pythia.info.hiInfo->nAbsTarg() != 0) {
        nAbsp += pythia.info.hiInfo->nAbsTarg();
      } else if (pythia.info.hiInfo->nDiffProj() != 0 || pythia.info.hiInfo->nDiffTarg() != 0) {
        nDiffp += (pythia.info.hiInfo->nDiffProj() + pythia.info.hiInfo->nDiffTarg());
      }
    }
    if (Npart > 87 && Npart < 130) {
      if (pythia.info.hiInfo->nAbsTarg() != 0) {
        nAbsmc += pythia.info.hiInfo->nAbsTarg();
      } else if (pythia.info.hiInfo->nDiffProj() != 0 || pythia.info.hiInfo->nDiffTarg() != 0) {
        nDiffmc += (pythia.info.hiInfo->nDiffProj() + pythia.info.hiInfo->nDiffTarg());
      }
    }
    if (Npart > 321) {
      if (pythia.info.hiInfo->nAbsTarg() != 0) {
        nAbsc += pythia.info.hiInfo->nAbsTarg();
      } else if (pythia.info.hiInfo->nDiffProj() != 0 || pythia.info.hiInfo->nDiffTarg() != 0) {
        nDiffc += (pythia.info.hiInfo->nDiffProj() + pythia.info.hiInfo->nDiffTarg());
      }
    } */
    
    // First sum up transverse energy for centrality measure and also
    // check that the trigger requiring ar least one charged particle
    // forward and backward.
    double etfwd = 0.0;
    bool trigfwd = false;
    bool trigbwd = false;
    int nc = 0;
    double x = 0;
    double y = 0;
    double z = 0;
    double px = 0;
    double py = 0;
    double et = 0;
    double phi = 0;
    double r_squared_cos_2phi = 0;
    double r_squared_sin_2phi = 0;
    double r_squared_cos_3phi = 0;
    double r_squared_sin_3phi = 0; 
    double r_squared = 0;
    int numParticles = 0;
    int numPairs = 0;
    v_2 = 0;
    v_3 = 0;
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() ) {
        double eta = p.eta();
        if ( p.isCharged() && p.pT() > 0.1 && eta < -2.09 && eta > -3.84 )
          trigfwd = true;
        if ( p.isCharged() && p.pT() > 0.1 && eta > 2.09 && eta < 3.84 )
          trigbwd = true;
        if ( p.pT() > 0.1 && abs(eta) > 3.2 && abs(eta) < 4.9 )
          etfwd += p.eT();
        if ( p.isCharged() && p.pT() > 0.1 && abs(eta) < 0.5 ) ++nc;
      }
      if (abs(pythia.event[i].status()) == 13) {
        x = pythia.event[i].xProd();
        y = pythia.event[i].yProd();
        z = pythia.event[i].zProd();
        px = pythia.event[i].px();
        py = pythia.event[i].py();
        phi = TMath::ATan2(y,x);
        et += pythia.event[i].e()/TMath::CosH(pythia.event[i].eta());
        
        r_squared += (x*x + y*y);
        r_squared_cos_2phi += (x*x + y*y) * (cos(2*phi));
        r_squared_sin_2phi += (x*x + y*y) * (sin(2*phi));
        r_squared_cos_3phi += (x*x + y*y) * (cos(3*phi));
        r_squared_sin_3phi += (x*x + y*y) * (sin(3*phi));
        
        numParticles++;
      }

      double phi_final1 = 0;
      double phi_final2 = 0;
      
      double eta_1 = 0;
      double eta_2 = 0;
      


      if (p.isFinal() && p.isCharged()) { 
        phi_final1 = atan2(p.py(), p.px());
        eta_1 = p.eta();
        for (int j = i + 1; j < pythia.event.size(); ++j) {
          Particle & p2 = pythia.event[j];
          eta_2 = p2.eta();
          if (p2.isFinal() && p2.isCharged()) {
            phi_final2 = atan2(p2.py(), p2.px());
            del_phi = phi_final1 - phi_final2;
            del_eta = eta_1 - eta_2;
            v_2 += cos(2*del_phi);
            v_3 += cos(3*del_phi);
            numPairs++;
            PairTree->Fill();
          }
        }
      } 
    }
    if (numPairs > 0) {
      v_2 /= numPairs;
      v_3 /= numPairs;
    } else {
      v_2 = 0;
      v_3 = 0;
    }
    if (numParticles > 0){
      r_squared_cos_2phi /= numParticles; 
      r_squared_sin_2phi /= numParticles;
      r_squared_cos_3phi /= numParticles;
      r_squared_sin_3phi /= numParticles;
      r_squared /= numParticles;
    }
    if (r_squared > 0) {
      epsilon_2 = pow((pow(r_squared_cos_2phi, 2) + pow(r_squared_sin_2phi, 2)), 0.5) / (r_squared);
      epsilon_3 = pow((pow(r_squared_cos_3phi, 2) + pow(r_squared_sin_3phi, 2)), 0.5) / (r_squared);
    } else {
      epsilon_2 = 0;
      epsilon_3 = 0;
    }
    theTree->Fill();
    // Skip if not triggered


    if ( !(trigfwd && trigbwd) ) continue;

    // Histogram and save the summed Et.
    //sumet->Fill(etfwd, weight);
    gencent.insert(make_pair(etfwd, weight));

    // Also fill the number of (absorptively and diffractively)
    // wounded nucleaons.
    /* int nw = pythia.info.hiInfo->nAbsTarg() +
      pythia.info.hiInfo->nDiffTarg() +
      pythia.info.hiInfo->nAbsProj() +
      pythia.info.hiInfo->nDiffProj();
    wounded->Fill(nw, weight);
 */
    // Find the correct centrality histograms.
    MapIdx::iterator genit = genetaidx.upper_bound(etfwd);
    int genidx = genit== genetaidx.end()? -1: genit->second;

    // Sum the weights in the centrality classes, skip if not in a class.
   /*  if ( genidx < 0 ) continue;
    gensumw[genidx] += weight;
    hmult[genidx]->fill(nc, weight);
    lmult[genidx]->fill(nc, weight);
    gensumn[genidx] += 1.0;
    cmult[genidx] += nc*weight;
    cmult2[genidx] += nc*nc*weight;
    wound[genidx] += nw*weight;
    wound2[genidx] += nw*nw*weight;

    // Go through the event again and fill the eta distributions.
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() && p.isCharged() &&
           abs(p.eta()) < 2.7 && p.pT() > 0.1 ) {
         etadist[genidx]->Fill(p.eta(), weight);
      }
    }
   */
  }
  // The run is over, so we write out some statistics.


  // Now, we just have to normalize and prtint out the histograms. We
  // choose to print the histograms to a file that can be read by
  // eg. gnuplot.
  TFile *output = new TFile("OO_DS_1k.root", "RECREATE");

  theTree->Write();

  const double del_phi_eta_bins = 55;
  const double del_Xmax = 4.5;
  const double del_Ymax = 1.5;
  const double del_Xmin = -0.5;
  const double del_Ymin = -1.5;

  TAxis* del_phi_axis = getLinAxis(del_Xmin, del_Xmax, del_phi_eta_bins, "#delta#phi", "#delta#phi");
  TAxis* del_eta_axis = getLinAxis(del_Ymin, del_Ymax, del_phi_eta_bins, "#delta#eta", "#delta#eta");

  TH2D* del_hist = makeTH2D("del_hist", del_phi_axis, del_eta_axis);

  Long64_t nentries = PairTree->GetEntries();
    for (Long64_t i=0;i<nentries;i++) {
      PairTree->GetEntry(i);
      del_hist->Fill(del_phi, del_eta);
    }
  del_hist->Write("del_hist");

  output->Close();

  /* sumet /= sumw*2.0;
  ofs << "# " << sumet.getTitle() << endl;
  sumet.table(ofs);

  wounded /= sumw*2.0;
  ofs << "\n# " << wounded.getTitle() << endl;
  wounded.table(ofs); */

  // Print out the centrality binned eta distributions and delete the
  // heap-allocate histograms.
  /* for ( int idx = 0; idx < 9; ++idx ) {
    *hmult[idx] /= gensumw[idx]*40.0;
    ofs << "\n# " << hmult[idx]->getTitle() << endl;
    hmult[idx]->table(ofs);
    delete hmult[idx];
    *lmult[idx] /= gensumw[idx]*4.0;
    ofs << "\n# " << lmult[idx]->getTitle() << endl;
    lmult[idx]->table(ofs);
    delete lmult[idx];
    *etadist[idx] /= gensumw[idx]*0.1;
    ofs << "\n# " << etadist[idx]->getTitle() << endl;
    etadist[idx]->table(ofs);
    delete etadist[idx];
  }

  // Print out average central charged multiplicity as a function of
  // centrality.
  ofs << "\n# Nch0\n";
  for ( int idx = 0; idx < 9; ++idx ) {
    double Nch = cmult[idx]/gensumw[idx];
    cmult2[idx] = (cmult2[idx]/gensumw[idx] - pow2(Nch))/gensumn[idx];
    ofs << setprecision(2) << setw(4) << int(pclim[idx]*100.0 + 0.5)
        << setw(10) << Nch << setw(10) << sqrt(cmult2[idx]) <<endl;
  }
  ofs << "\n# Nwc\n";
  for ( int idx = 0; idx < 9; ++idx ) {
    double Nw = wound[idx]/gensumw[idx];
    wound2[idx] = (wound2[idx]/gensumw[idx] - pow2(Nw))/gensumn[idx];
    ofs << setprecision(2) << setw(4) << int(pclim[idx]*100.0 + 0.5)
        << setw(10) << Nw << setw(10) << sqrt(wound2[idx]) <<endl;
  } */

  // Befor we end we print out some statistics. Also, we want to check
  // that our generated centrality classes were the same as we
  // guessed.
  cout << "nDiffp: " << nDiffp << "\n"
       << "nAbsp: " << nAbsp << "\n"
       << "nDiffmc: " << nDiffmc << "\n"
       << "nAbsmc: " << nAbsmc << "\n"
       << "nDiffc: " << nDiffc << "\n"
       << "nAbsc: " << nAbsc << "\n";
  pythia.stat();
  double curr = 0.0;
  double prev = 0.0;
  double acc = 0.0;
  int idxa = 8;
  double lim = sumw*(1.0 - pclim[idxa]);
  vector<double> newlim(9);
  for ( multimap<double, double>::iterator it = gencent.begin();
        it != gencent.end(); ++it ) {
    prev = curr;
    curr = it->first;
    double w = it->second;
    if ( acc < lim && acc + w >= lim ) {
      newlim[idxa--] = prev + (curr - prev)*(lim - acc)/w;
      if ( idxa < 0 ) break;
      lim = sumw*(1.0 - pclim[idxa]);
    }
    acc += w;
  }

  cout << "The generated limits between centrality classes in this run:\n"
       << "   %   assumed    actual\n";
  for ( int idx = 0; idx < 9; ++idx )
    cout << setw(4) << int(pclim[idx]*100.0 + 0.5)
         << setw(10) << fixed << setprecision(1) << genlim[idx]
         << setw(10) << fixed << setprecision(1) << newlim[idx] << endl;

  // And we're done!
  return 0;
}
