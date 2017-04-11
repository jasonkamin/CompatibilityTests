#include <math.h>

TH1D *xjz_raw_stat_pp;
TH1D *xjz_raw_stat_pb;
TH1D *xjz_bkg_stat_pb;

TH1D *dphi_raw_stat_pp;
TH1D *dphi_raw_stat_pb;
TH1D *dphi_bkg_stat_pb;

TH1D *dphi_sum_stat_pp;
TH1D *xjz_sum_stat_pp;

TH1D *xjz_raw_syst_pp;
TH1D *xjz_raw_syst_pb;
TH1D *xjz_bkg_syst_pb;

TH1D *dphi_raw_syst_pp;
TH1D *dphi_raw_syst_pb;
TH1D *dphi_bkg_syst_pb;

TH1D *h_Chisq;
TH1D *h_KolSm;
TH1D *h_AndDr;

void MakeDataCardFromHistos(TH1D *h1, TH1D *h2, TH1D *h3, int startbin, int nBins);
void PrintSystUncertainty(double *percent_uncert, TH1 *h_data, int FirstBin, int LastBin, char *category, char *type);
void PrintSystUncertainty2(double *percent_uncert_lo, double *percent_uncert_hi, TH1 *h_data, int FirstBin, int LastBin, char *category, char *type);

// 21 March 2017
// this should make the fully combined datacard !
// now can accept asymmetric, correlated uncertainties.  
// however, they are input inside of MakeDataCardFromHistos(), 
// which makes it actually pretty inconvenient !! 
void MakeCombineDatacard()
{

  double BlowUpErrors = 1.0;
  char saythis[500];
  int pTbin = 0;

  //TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_5Jul2016.root");
  //TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSumPP_5Jul2016.root");
  //TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_HI_PP_28Jul2016.root");
  //TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_HI_PP_28Jul2016.root");

  //TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_HI_PP_31Aug2016.root");
  //TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_HI_PP_31Aug2016.root");
  TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/Jan23_2017/zJetHistogramSum_HI_PP_20170118.root");
  TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/Jan23_2017/zJetHistogramSum_HI_PP_20170118.root");

  //TFile *f_pb = TFile::Open("$CODE/ZJet/Frankenstein.root");
  //TFile *f_pp = TFile::Open("$CODE/ZJet/Frankenstein.root");

  //xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_ptBin0_hiBin1_jetRAW_final_norm");
  //xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_ptBin0_hiBin1_jetBKG_final_norm");
  //xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get("h1D_xjz_ptBin0_hiBin0_jetRAW_final_norm");
  sprintf(saythis,"h1D_xjz_binJER_ptBin%d_hiBin1_jetRAW_final_norm",pTbin); xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get(saythis);
  sprintf(saythis,"h1D_xjz_binJER_ptBin%d_hiBin1_jetBKG_final_norm",pTbin); xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get(saythis);
  sprintf(saythis,"h1D_xjz_binJER_ptBin%d_hiBin0_jetRAW_final_norm",pTbin); xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get(saythis);

  //xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("xjz_raw_stat_pb");
  //xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("xjz_bkg_stat_pb");
  //xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get("xjz_raw_stat_pp");

  sprintf(saythis,"h1D_dphi_rebin_ptBin%d_hiBin1_jetRAW_final_norm",pTbin); dphi_raw_stat_pb = (TH1D*)f_pb->GetDirectory("HI")->Get(saythis);
  sprintf(saythis,"h1D_dphi_rebin_ptBin%d_hiBin1_jetBKG_final_norm",pTbin); dphi_bkg_stat_pb = (TH1D*)f_pb->GetDirectory("HI")->Get(saythis);
  sprintf(saythis,"h1D_dphi_rebin_ptBin%d_hiBin0_jetRAW_final_norm",pTbin); dphi_raw_stat_pp = (TH1D*)f_pp->GetDirectory("PP")->Get(saythis);

  for(int i=1; i<xjz_raw_stat_pp->GetNbinsX()+1; i++){
    xjz_raw_stat_pp->SetBinError(i,BlowUpErrors*xjz_raw_stat_pp->GetBinError(i));
    xjz_bkg_stat_pb->SetBinError(i,BlowUpErrors*xjz_bkg_stat_pb->GetBinError(i));
  }
  for(int i=1; i<dphi_raw_stat_pp->GetNbinsX()+1; i++){
    dphi_raw_stat_pp->SetBinError(i,BlowUpErrors*dphi_raw_stat_pp->GetBinError(i));
    dphi_bkg_stat_pb->SetBinError(i,BlowUpErrors*dphi_bkg_stat_pb->GetBinError(i));
  }

  //dphi_raw_stat_pb ->Scale(301.562);
  //dphi_bkg_stat_pb ->Scale(301.562);
  //dphi_raw_stat_pp ->Scale(673);
  //dphi_raw_stat_pp ->Scale(301.562);

  double QuickError = 0.0;
  double QuickIntegral_pp = dphi_raw_stat_pp->IntegralAndError(-1,-1,QuickError,"width");
  cout << " integral rawpp: " << QuickIntegral_pp << " +/- " << QuickError << endl;
  double QuickIntegral_pb = dphi_raw_stat_pb->IntegralAndError(-1,-1,QuickError,"width");
  double TotalEqvStatsPbPb = QuickIntegral_pb*QuickIntegral_pb/(QuickError*QuickError);
  cout << " integral rawpb: " << QuickIntegral_pb << " +/- " << QuickError << endl;
  double QuickIntegral_bg = dphi_bkg_stat_pb->IntegralAndError(-1,-1,QuickError,"width");
  cout << " integral bkgpb: " << QuickIntegral_bg << " +/- " << QuickError << endl;

  double RatioTotalIntegrals_pb_over_pp = QuickIntegral_pb/(QuickIntegral_pp+QuickIntegral_bg);

  gStyle->SetOptStat(0);
  TH1D *h_frame = new TH1D("h_frame"," ",1000,0,10);
  TCanvas *c1 = new TCanvas("c1","c1",900,600);
  c1->Divide(3,2);
  c1->cd(1);
  xjz_raw_stat_pb  ->SetMarkerStyle(20);
  xjz_bkg_stat_pb  ->SetMarkerStyle(24);
  xjz_raw_stat_pp  ->SetMarkerStyle(20);
  xjz_raw_stat_pb  ->SetMarkerColor(1);
  xjz_bkg_stat_pb  ->SetMarkerColor(2);
  xjz_raw_stat_pp  ->SetMarkerColor(4);
  xjz_raw_stat_pb  ->SetLineColor  (1);
  xjz_bkg_stat_pb  ->SetLineColor  (2);
  xjz_raw_stat_pp  ->SetLineColor  (4);
  h_frame  ->GetYaxis()->SetRangeUser(0,1.5);
  h_frame  ->GetXaxis()->SetRangeUser(0,2.2);
  h_frame  ->GetXaxis()->SetTitle("XJZ");
  h_frame->DrawCopy();
  xjz_raw_stat_pb  ->DrawCopy("same");
  xjz_bkg_stat_pb  ->DrawCopy("same");
  xjz_raw_stat_pp  ->DrawCopy("same");

  TLegend *c1a = new TLegend(0.4,0.65,0.89,0.89);
  c1a->AddEntry(xjz_raw_stat_pb,"PbPb non-sub","P");
  c1a->AddEntry(xjz_bkg_stat_pb,"PbPb mix evts","P");
  c1a->AddEntry(xjz_raw_stat_pp,"pp","P");
  c1a->Draw();

  c1->cd(4);
  dphi_raw_stat_pb  ->SetMarkerStyle(20);
  dphi_bkg_stat_pb  ->SetMarkerStyle(24);
  dphi_raw_stat_pp  ->SetMarkerStyle(20);
  dphi_raw_stat_pb  ->SetMarkerColor(1);
  dphi_bkg_stat_pb  ->SetMarkerColor(2);
  dphi_raw_stat_pp  ->SetMarkerColor(4);
  dphi_raw_stat_pb  ->SetLineColor  (1);
  dphi_bkg_stat_pb  ->SetLineColor  (2);
  dphi_raw_stat_pp  ->SetLineColor  (4);
  h_frame  ->GetYaxis()->SetRangeUser(0,4);
  h_frame  ->GetXaxis()->SetRangeUser(0,3.2);
  h_frame  ->GetXaxis()->SetTitle("dPhi");
  h_frame->DrawCopy();
  dphi_raw_stat_pb ->DrawCopy("same");
  dphi_bkg_stat_pb ->DrawCopy("same");
  dphi_raw_stat_pp ->DrawCopy("same");


  dphi_sum_stat_pp = (TH1D*)dphi_raw_stat_pp->Clone("dphi_sum_stat_pp");
  dphi_sum_stat_pp ->Add(dphi_bkg_stat_pb);
  //dphi_sum_stat_pp ->Scale(dphi_raw_stat_pb->Integral(-1,-1)/dphi_sum_stat_pp->Integral(-1,-1));

  xjz_sum_stat_pp = (TH1D*)xjz_raw_stat_pp->Clone("xjz_sum_stat_pp");
  xjz_sum_stat_pp ->Add(xjz_bkg_stat_pb);
  //xjz_sum_stat_pp ->Scale(xjz_raw_stat_pb->Integral(-1,-1)/xjz_sum_stat_pp->Integral(-1,-1));


  cout << endl << endl;

  cout << "=====================================================" << endl << endl;
  //MakeDataCardFromHistos(dphi_raw_stat_pp, dphi_raw_stat_pb, dphi_bkg_stat_pb, -1,-1);
  cout << "=====================================================" << endl << endl;
  //MakeDataCardFromHistos(xjz_raw_stat_pp, xjz_raw_stat_pb, xjz_bkg_stat_pb,1,8);// using 10 bin scheme "binJER"
  //MakeDataCardFromHistos(xjz_raw_stat_pp, xjz_raw_stat_pb, xjz_bkg_stat_pb,2,12);//using 16 bin scheme (non JER)
  MakeDataCardFromHistos(xjz_raw_stat_pp, xjz_raw_stat_pb, xjz_bkg_stat_pb,-1,-1);
  cout << "=====================================================" << endl << endl;

  return;





}



void MakeDataCardFromHistos(TH1D *h_pp, TH1D *h_pb, TH1D *h_bg, int startbin, int nBins)
{

  char saythis[500];
  //int startbin = 2;
  //int nBins = 13;
  //int nBins = h_pp->GetNbinsX();

  if(startbin<0)
    startbin = 1;
  if(nBins<0)
    nBins = h_pp->GetNbinsX();

  cout << "# Made with MakeCombineDatacard.C" << endl;
  cout << "# using histograms" << endl
       << "#   " << h_pp->GetName() << endl
       << "#   " << h_pb->GetName() << endl
       << "#   " << h_bg->GetName() << endl;
  cout << "# " << endl;

  cout << "# -----------INPUT-----------" << endl;
  cout << "#   " << h_pp->GetName() << endl;
  cout << "#    bin" << "   Content" << "     Error" << endl;
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"# %6d%10.3f%10.3f",i,h_pp->GetBinContent(i),h_pp->GetBinError(i));
    cout << saythis << endl;
  }
  cout << "# ---------------------------" << endl << endl;
  cout << "#   " << h_pb->GetName() << endl;
  cout << "#    bin" << "   Content" << "     Error" << endl;
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"# %6d%10.3f%10.3f",i,h_pb->GetBinContent(i),h_pb->GetBinError(i));
    cout << saythis << endl;
  }
  cout << "# ---------------------------" << endl << endl;
  cout << "#   " << h_bg->GetName() << endl;
  cout << "#    bin" << "   Content" << "     Error" << endl;
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"# %6d%10.3f%10.3f",i,h_bg->GetBinContent(i),h_bg->GetBinError(i));
    cout << saythis << endl;
  }
  cout << "# ---------------------------" << endl << endl;


  cout << "imax " << (nBins+1 - startbin)*2 << "  number of bins"                << endl;
  cout << "jmax " << "2"   << "  number of processes minus 1"         << endl;
  cout << "kmax " << "*"   << "  number of nuisance parameters" << endl;
  cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;


  cout << endl;
  cout << "# observation is the 'effective counts' in pp and PbPb for each bin. " << endl;

  cout << "bin        ";
  for(int j=0; j<2; j++){
    for(int i=startbin; i<nBins+1; i++){
      if(i<10)
        sprintf(saythis,"  ch%d_b%d",j+1,i);
      else
        sprintf(saythis," ch%d_b%d",j+1,i);
      cout << saythis;
    }
  }
  cout << endl;

  cout << "observation";
  for(int i=startbin; i<nBins+1; i++){
    if(h_pp->GetBinContent(i)>0)
      sprintf(saythis,"%8d",int(floor( h_pp->GetBinContent(i)*h_pp->GetBinContent(i)/(h_pp->GetBinError(i)*h_pp->GetBinError(i))+0.5 )));
    else
      sprintf(saythis,"%8d",0);
    cout << saythis;
  }
  for(int i=startbin; i<nBins+1; i++){
    if(h_pb->GetBinContent(i)>0)
      sprintf(saythis,"%8d",int(floor( h_pb->GetBinContent(i)*h_pb->GetBinContent(i)/(h_pb->GetBinError(i)*h_pb->GetBinError(i))+0.5 )));
    else
      sprintf(saythis,"%8d",0);
    cout << saythis;
  }
  cout << endl;
  cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;


  cout << endl;
  cout << "# the 'rate' is the factor by which we multiply the 'physics value' to get units of effective counts" << endl;
  cout << "# for the bkg's, this something weird... it's the bkg phys_val * PbPb_eff (ie. the PbPb rate below)" << endl;

  cout << "bin        ";
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"   ch1_b%d",i);
    cout << saythis;
  }
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"   ch2_b%d",i);
    cout << saythis;
  }
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"   ch2_b%d",i);
    cout << saythis;
  }
  cout << endl;

  cout << "process    ";
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"%9s","pp");
    cout << saythis;
  }
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"%9s","PbPb");
    cout << saythis;
  }
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"%9s","bkg");
    cout << saythis;
  }
  cout << endl;

  cout << "process    ";
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"%9d",1);
    cout << saythis;
  }
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"%9d",3);
    cout << saythis;
  }
  for(int i=startbin; i<nBins+1; i++){
    sprintf(saythis,"%9d",2);
    cout << saythis;
  }
  cout << endl;

  cout << "rate       ";
  for(int i=startbin; i<nBins+1; i++){
    //sprintf(saythis,"%10.6f", h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i) );
    if(h_pp->GetBinContent(i)>0)
      sprintf(saythis,"%9.4f", h_pp->GetBinContent(i)/(h_pp->GetBinError(i)*h_pp->GetBinError(i)) );
    else
      sprintf(saythis,"  0.0001");
    cout << saythis;
  }
  for(int i=startbin; i<nBins+1; i++){
    //sprintf(saythis,"%10.6f", h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i) );
    if(h_pb->GetBinContent(i)>0)
      sprintf(saythis,"%9.4f", h_pb->GetBinContent(i)/(h_pb->GetBinError(i)*h_pb->GetBinError(i)) );
    else
      sprintf(saythis,"  0.0001");
    cout << saythis;
  }
  for(int i=startbin; i<nBins+1; i++){
    //sprintf(saythis,"%10.6f", h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i) );
    if(h_pb->GetBinContent(i)>0 && h_bg->GetBinContent(i)>0)
      sprintf(saythis,"%9.4f", h_bg->GetBinContent(i)*h_pb->GetBinContent(i)/(h_pb->GetBinError(i)*h_pb->GetBinError(i)) );
    else
      sprintf(saythis,"  0.0001");
    cout << saythis;
  }
  cout << endl;
  cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;

  cout << endl;
  cout << "# these are the PbPb 'physics values' divided by the pp" << endl;
  cout << "# 'physics values' means 'what is in the final fully-normalized plot'" << endl;

  for(int i=startbin; i<nBins+1; i++){
    if(h_pp->GetBinContent(i)>0){
      if(i<10)
        sprintf(saythis,"dx%d           param  %4.2f  100  [0,100]",i,h_pb->GetBinContent(i)/h_pp->GetBinContent(i));
      else
        sprintf(saythis,"dx%d          param  %4.2f  100  [0,100]",i,h_pb->GetBinContent(i)/h_pp->GetBinContent(i));
    }
    else{
      if(i<10)
        sprintf(saythis,"dx%d           param  0.00  100  [0,100]",i);
      else
        sprintf(saythis,"dx%d          param  0.00  100  [0,100]",i);
    }
    cout << saythis << endl;
  }

  cout << endl;
  cout << "# pp 'physics values'" << endl;

  for(int i=startbin; i<nBins+1; i++){
    if(i<10)
      sprintf(saythis,"x%d            rateParam  ch1_b%d   pp    %4.3f  [0,%2.2f]",i,i,h_pp->GetBinContent(i),10.0*h_pp->GetBinContent(h_pp->GetMaximumBin()));
    else
      sprintf(saythis,"x%d           rateParam  ch1_b%d  pp    %4.3f  [0,%2.2f]",i,i,h_pp->GetBinContent(i),10.0*h_pp->GetBinContent(h_pp->GetMaximumBin()));
    cout << saythis << endl;
  }
  for(int i=startbin; i<nBins+1; i++){
    if(i<10)
      sprintf(saythis,"y%d            rateParam  ch2_b%d   PbPb  @0*@1  x%d,dx%d",i,i,i,i);
    else
      sprintf(saythis,"y%d           rateParam  ch2_b%d  PbPb  @0*@1  x%d,dx%d",i,i,i,i);
    cout << saythis << endl;
  }


  for(int j=startbin; j<nBins+1; j++){
    sprintf(saythis,"statbkg%d  lnN    ",j);
    cout << saythis;
    for(int i=startbin; i<nBins+1; i++){
      cout << "  -     ";
    }
    for(int i=startbin; i<nBins+1; i++){
      cout << "  -     ";
    }
    for(int i=startbin; i<nBins+1; i++){
      if(i==j){
        if(h_bg->GetBinContent(i)>0)
          sprintf(saythis,"  %4.4f ",1.0+h_bg->GetBinError(i)/h_bg->GetBinContent(i));
        else
          sprintf(saythis,"  -     ");
        cout << saythis;
      }
      else
        cout << "  -     ";
    }
    cout << endl;
  }


  char namepp[5]; sprintf(namepp,"pp");
  char namepb[5]; sprintf(namepb,"Pb");
  char namebg[5]; sprintf(namebg,"Bg");
  char namecat[50];

  double xjz_syst_pp[10]   = {59.8,  4.9,  3.1,  3.1,  4.1,  7.1,  7.1,  7.1,  7.1,  7.1};
  double xjz_syst_PbPb[10] = {44.6, 11.4,  9.8,  9.8, 12.4, 18.6, 18.6, 18.6, 18.6, 18.6};
  double xjz_syst_bkg[10]  = {19.1,  8.8,  9.1,  9.1,  7.8, 12.4,  8.6,  8.6,  0.1,  0.1};

  //sprintf(namecat,"Tot");
  //PrintSystUncertainty(xjz_syst_pp,   h_pp, startbin, nBins, namecat, namepp);
  //PrintSystUncertainty(xjz_syst_PbPb, h_pb, startbin, nBins, namecat, namepb);
  //PrintSystUncertainty(xjz_syst_bkg,  h_bg, startbin, nBins, namecat, namebg);


  //############################### xjz SYSTEMATICS ZpT>60 #####################################
  //PB, SIG correlation : (currently not used; we separate SIG into RAW and BKG)
  //double xjz_syst_JES_pb[10] = {10.5   10.5   4.5    4.5    7.2    15.7   15.7   15.7   15.7 15.7   jet energy scale
  //double xjz_syst_EES_pb[10] = {0.9    0.9    1.1    1.1    1.1    1.7    1.7    1.7    1.7  1.7     electron energy scale
  //double xjz_syst_JER_pb[10] = {0.3    0.3    6.8    6.8    9.8    2.5    2.5    2.5    2.5  2.5     jet energy resolution
  //double xjz_syst_ZER_pb[10] = {3.4    3.4    1.4    1.4    1.2    6.9    6.9    6.9    6.9  6.9     Z energy resolution
  //double xjz_syst_ZEF_pb[10] = {5.1    5.1    6.4    6.4    4.0    6.2    6.2    6.2    6.2  6.2     Z reconstruction efficiency correction
  //double xjz_syst_TOT_pb[10] = {12.2   12.2   10.5   10.5   12.9   18.5   18.5   18.5   18.5 18.5   total


// ########## (sent from Kaya on March 17th, 2017)
// 1. For pp, you will need to feed the JES uncertainties twice.
// 2. For PbPb, you will need to feed the JES1 once and JES2 twice.
// 
// xjz (zPt > 60) JES up and down systematics
// #####
// pp
// JES UP   :  2.7,  2.7,  1.2,  1.2,  1.7, -5.3, -5.3, -5.3, -5.3, -5.3
// JES DOWN : -2.0, -2.0, -1.9, -1.9, -1.5,  5.5,  5.5,  5.5,  5.5, 5.5
// 
// #####
// PbPb (zpT > 60), SIG 
// JES1 UP   :  8.7,   8.7,  3.2, 3.2, -6.6, -14.2, -14.2, -14.2, -14.2, -14.2
// JES1 DOWN : -10.1, -10.1, 2.7, 2.7,  4.8,   9.8,   9.8,   9.8, 9.8,   9.8
// JES2 UP   :  1.8,   1.8,  1.2, 1.2, -1.3,  -1.4,  -1.4,  -1.4, -1.4,  -1.4
// JES2 DOWN : -2.8,  -2.8,  2.9, 2.9, -0.2,   7.6,   7.6,   7.6, 7.6,   7.6
// 
// ##########
// xjz (zPt > 60), RAW correlation
// #####
// PbPb
// JES1 UP   : -7.2, -7.2, -3.5, -3.5, -7.5, -12.0, -12.0, -12.0, -12.0, -12.0
// JES1 DOWN :  2.0, 2.0, 4.5, 4.5, 5.3, 7.5, 7.5, 7.5, 7.5, 7.5
// JES2 UP   : -3.5, -3.5, -1.3, -1.3, -3.5, 1.2, 1.2, 1.2, 1.2, 1.2
// JES2 DOWN :  3.0, 3.0, 3.6, 3.6, 3.2, 7.8, 7.8, 7.8, 7.8, 7.8
// 
// ##########
// xjz (zPt > 60), BKG correlation
// #####
// PbPb
// JES1 UP   : -24.8, -24.8, -27.5, -27.5, -43.8, -20.1, -20.1, -20.1, -20.1, -20.1
// JES1 DOWN :  16.3, 16.3, 16.4, 16.4, 20.8, 20.5, 20.5, 20.5, 20.5, 20.5
// JES2 UP   :  -9.9, -9.9, -10.6, -10.6, -13.7, -10.8, -10.8, -10.8, -10.8, -10.8
// JES2 DOWN :   9.7, 9.7, 9.4, 9.4, 11.6, 11.2, 11.2, 11.2, 11.2, 11.2
// 
// 
//
// ##########
// rjz JES up and down systematics
// #####
// pp
// rjz bins   pt>40, pt>50, pt >60, 40-50, 50-60, 60-80, pt>80
// JES UP   :  -1.6, -1.5, -1.3, -2.4, -1.6, -1.6, -1.3
// JES DOWN :   1.6,  1.5,  1.2,  2.5,  1.6,  1.7,  1.2
// 
// #####
// PbPb
// rjz bins    pt>40, pt>50, pt >60, 40-50, 50-60, 60-80, pt>80
// JES1 UP   :  -3.3, -3.1, -1.9, -10.1, -8.1, -2.6, -2.0
// JES1 DOWN :   1.9,  2.0,  0.9,   7.5,  4.9,  2.3,  2.3
// JES2 UP   :  -2.1, -1.6, -1.4, -5.3, -4.7, -2.3, -1.6
// JES2 DOWN :   1.9,  1.8,  0.7,  5.5,  4.1,  1.4,  0.4
// 
// ##########
// rjz JES up and down systematics, RAW correlation
// #####
// PbPb
// rjz bins    pt>40, pt>50, pt >60, 40-50, 50-60, 60-80, pt>80
// JES1 UP   : -12.1, -10.8, -8.4, -20.4, -18.2, -9.8, -5.9
// JES1 DOWN :   7.1, 6.5, 4.8, 13.5, 11.3, 5.4, 3.0
// JES2 UP   :  -2.2, -1.7, -1.5, -3.5, -3.3, -3.0, 1.2
// JES2 DOWN :   4.6, 4.2, 3.1, 7.2, 6.5, 3.4, 2.1
// 
// ##########
// rjz JES up and down systematics, BKG correlation
// #####
// PbPb
// rjz bins    pt>40, pt>50, pt >60, 40-50, 50-60, 60-80, pt>80
// JES1 UP   : -29.2, -28.8, -27.5, -30.1, -30.8, -23.9, -25.2
// JES1 DOWN :  17.2, 17.2, 16.0, 19.2, 19.8, 15.6, 16.3
// JES2 UP   : -10.5, -10.3, -10.0, -11.7, -12.3, -9.8, -10.5
// JES2 DOWN :   9.1, 9.9, 9.9, 11.4, 11.6, 9.0, 9.1
// ##########



  //PP, SIG correlation :
  double xjz_syst_JES_pp[10] = {3.4   , 3.4   , 2.8   , 2.8   , 2.6   , 6.8   , 6.8   , 6.8   , 6.8  , 6.8 }; //  jet energy scale
  // if we use the un-symmetrized values, insert them twice, independently.
  double xjz_syst_JES_pp_lo[10] = {2.7,    2.7,    1.2,    1.2,    1.7,   -5.3,   -5.3,   -5.3,   -5.3,  -5.3 }; //  jet energy scale
  double xjz_syst_JES_pp_hi[10] = {-2.0,  -2.0,   -1.9,   -1.9,   -1.5,    5.5,    5.5,    5.5,    5.5,   5.5 }; //  jet energy scale

  double xjz_syst_EES_pp[10] = {2.4   , 2.4   , 1.1   , 1.1   , 1.6   , 2.1   , 2.1   , 2.1   , 2.1  , 2.1 }; //  electron energy scale
  double xjz_syst_JER_pp[10] = {2.6   , 2.6   , 0.8   , 0.8   , 2.6   , 0.3   , 0.3   , 0.3   , 0.3  , 0.3 }; //  jet energy resolution
  double xjz_syst_JAR_pp[10] = {0.4   , 0.4   , 0.4   , 0.4   , 0.1   , 0.1   , 0.1   , 0.1   , 0.1  , 0.1 }; //  jet angular resolution
  double xjz_syst_TOT_pp[10] = {4.9   , 4.9   , 3.1   , 3.1   , 4.1   , 7.1   , 7.1   , 7.1   , 7.1  , 7.1 }; //  total


  //PB RAW correlation :
  double xjz_syst_JES_pb[10] = {8.1   , 8.1   , 6.3   , 6.3   , 8.5   , 15.8  , 15.8  , 15.8  , 15.8 , 15.8}; //     jet energy scale
  // if we use the un-symmetrized values, insert JES1 once and JES2 twice, independently.
  double xjz_syst_JES1_pb_lo[10] = {-7.2, -7.2, -3.5, -3.5, -7.5, -12.0, -12.0, -12.0, -12.0, -12.0}; //     jet energy scale 1
  double xjz_syst_JES1_pb_hi[10] = { 2.0, 2.0, 4.5, 4.5, 5.3, 7.5, 7.5, 7.5, 7.5, 7.5              }; //     jet energy scale 1
  double xjz_syst_JES2_pb_lo[10] = {-3.5, -3.5, -1.3, -1.3, -3.5, 1.2, 1.2, 1.2, 1.2, 1.2          }; //     jet energy scale 2
  double xjz_syst_JES2_pb_hi[10] = { 3.0, 3.0, 3.6, 3.6, 3.2, 7.8, 7.8, 7.8, 7.8, 7.8              }; //     jet energy scale 2

  double xjz_syst_EES_pb[10] = {0.6   , 0.6   , 1.0   , 1.0   , 1.1   , 1.7   , 1.7   , 1.7   , 1.7  , 1.7 }; //    electron energy scale
  double xjz_syst_JER_pb[10] = {3.5   , 3.5   , 2.2   , 2.2   , 8.0   , 2.8   , 2.8   , 2.8   , 2.8  , 2.8 }; //    jet energy resolution
  double xjz_syst_ZER_pb[10] = {2.2   , 2.2   , 0.9   , 0.9   , 1.0   , 6.8   , 6.8   , 6.8   , 6.8  , 6.8 }; //    Z energy resolution
  double xjz_syst_ZEF_pb[10] = {6.8   , 6.8   , 7.0   , 7.0   , 4.1   , 6.2   , 6.2   , 6.2   , 6.2  , 6.2 }; //    Z reconstruction efficiency correction
  double xjz_syst_TOT_pb[10] = {11.4  , 11.4  , 9.8   , 9.8   , 12.4  , 18.6  , 18.6  , 18.6  , 18.6 , 18.6}; //   total


  //PB BKG correlation :
  double xjz_syst_JES_bg[10] = {26.9  , 26.9  , 30.6  , 30.6  ,  47.8 ,  23.5 ,  10.1 ,  10.1 , 0.1  , 0.1 }; //  jet energy scale
  // if we use the un-symmetrized values, insert JES1 once and JES2 twice, independently.
  double xjz_syst_JES1_bg_lo[10] = { -24.8, -24.8, -27.5, -27.5, -43.8, -20.1, -20.1, -20.1, -20.1, -20.1}; //     jet energy scale 1
  double xjz_syst_JES1_bg_hi[10] = {  16.3, 16.3, 16.4, 16.4, 20.8, 20.5, 20.5, 20.5, 20.5, 20.5         }; //     jet energy scale 1
  double xjz_syst_JES2_bg_lo[10] = {  -9.9, -9.9, -10.6, -10.6, -13.7, -10.8, -10.8, -10.8, -10.8, -10.8 }; //     jet energy scale 2
  double xjz_syst_JES2_bg_hi[10] = {   9.7, 9.7, 9.4, 9.4, 11.6, 11.2, 11.2, 11.2, 11.2, 11.2            }; //     jet energy scale 2

  double xjz_syst_EES_bg[10] = {0.4   , 0.4   , 0.8   , 0.8   ,  4.6  ,  1.0  ,  1.0  ,  1.0  , 0.1  , 0.1 }; //  electron energy scale
  double xjz_syst_JER_bg[10] = {7.8   , 7.8   , 15.4  , 15.4  ,  46.8 ,  22.0 ,  1.0  ,  2.0  , 0.1  , 0.1 }; //  jet energy resolution
  double xjz_syst_ZER_bg[10] = {0.8   , 0.8   , 1.0   , 1.0   ,  6.4  ,  5.2  ,  26.0 ,  26.0 , 0.1  , 0.1 }; //  Z energy resolution
  double xjz_syst_ZEF_bg[10] = {8.8   , 8.8   , 9.1   , 9.1   ,  7.8  ,  12.4 ,  8.6  ,  8.6  , 0.1  , 0.1 }; //  Z reconstruction efficiency correction
  double xjz_syst_TOT_bg[10] = {29.4  , 29.4  , 35.5  , 35.5  ,  67.8 ,  34.9 ,  29.3 ,  29.3 , 0.1  , 0.1 }; //  total

  //sprintf(namecat,"JES");  PrintSystUncertainty(xjz_syst_JES_pp, h_pp, startbin, nBins, namecat, namepp);
  sprintf(namecat,"JESa");  PrintSystUncertainty2(xjz_syst_JES_pp_lo, xjz_syst_JES_pp_hi, h_pp, startbin, nBins, namecat, namepp);
  sprintf(namecat,"JESb");  PrintSystUncertainty2(xjz_syst_JES_pp_lo, xjz_syst_JES_pp_hi, h_pp, startbin, nBins, namecat, namepp);
  sprintf(namecat,"EES");  PrintSystUncertainty(xjz_syst_EES_pp, h_pp, startbin, nBins, namecat, namepp);
  sprintf(namecat,"JER");  PrintSystUncertainty(xjz_syst_JER_pp, h_pp, startbin, nBins, namecat, namepp);
  sprintf(namecat,"JAR");  PrintSystUncertainty(xjz_syst_JAR_pp, h_pp, startbin, nBins, namecat, namepp);
                      
  //sprintf(namecat,"JES");  PrintSystUncertainty(xjz_syst_JES_pb, h_pb, startbin, nBins, namecat, namepb);
  sprintf(namecat,"JES1");   PrintSystUncertainty2(xjz_syst_JES1_pb_lo, xjz_syst_JES1_pb_hi, h_pb, startbin, nBins, namecat, namepb);
  sprintf(namecat,"JES2a");  PrintSystUncertainty2(xjz_syst_JES2_pb_lo, xjz_syst_JES2_pb_hi, h_pb, startbin, nBins, namecat, namepb);
  sprintf(namecat,"JES2b");  PrintSystUncertainty2(xjz_syst_JES2_pb_lo, xjz_syst_JES2_pb_hi, h_pb, startbin, nBins, namecat, namepb);
  sprintf(namecat,"EES");  PrintSystUncertainty(xjz_syst_EES_pb, h_pb, startbin, nBins, namecat, namepb);
  sprintf(namecat,"JER");  PrintSystUncertainty(xjz_syst_JER_pb, h_pb, startbin, nBins, namecat, namepb);
  sprintf(namecat,"ZER");  PrintSystUncertainty(xjz_syst_ZER_pb, h_pb, startbin, nBins, namecat, namepb);
  sprintf(namecat,"ZEF");  PrintSystUncertainty(xjz_syst_ZEF_pb, h_pb, startbin, nBins, namecat, namepb);
                      
  //sprintf(namecat,"JES");  PrintSystUncertainty(xjz_syst_JES_bg, h_bg, startbin, nBins, namecat, namebg);
  sprintf(namecat,"JES1");   PrintSystUncertainty2(xjz_syst_JES1_bg_lo, xjz_syst_JES1_bg_hi, h_bg, startbin, nBins, namecat, namebg);
  sprintf(namecat,"JES2a");  PrintSystUncertainty2(xjz_syst_JES2_bg_lo, xjz_syst_JES2_bg_hi, h_bg, startbin, nBins, namecat, namebg);
  sprintf(namecat,"JES2b");  PrintSystUncertainty2(xjz_syst_JES2_bg_lo, xjz_syst_JES2_bg_hi, h_bg, startbin, nBins, namecat, namebg);
  sprintf(namecat,"EES");  PrintSystUncertainty(xjz_syst_EES_bg, h_bg, startbin, nBins, namecat, namebg);
  sprintf(namecat,"JER");  PrintSystUncertainty(xjz_syst_JER_bg, h_bg, startbin, nBins, namecat, namebg);
  sprintf(namecat,"ZER");  PrintSystUncertainty(xjz_syst_ZER_bg, h_bg, startbin, nBins, namecat, namebg);
  sprintf(namecat,"ZEF");  PrintSystUncertainty(xjz_syst_ZEF_bg, h_bg, startbin, nBins, namecat, namebg);

  //############################### xjz SYSTEMATICS #####################################


  //for(int i=0; i<10; i++){
  //  cout << TMath::Sqrt(xjz_syst_JES1_pb_lo[i]*xjz_syst_JES1_pb_lo[i] + xjz_syst_JES2_pb_lo[i]*xjz_syst_JES2_pb_lo[i] + xjz_syst_JES2_pb_lo[i]*xjz_syst_JES2_pb_lo[i]) << endl;
  //  cout << TMath::Sqrt(xjz_syst_JES1_pb_hi[i]*xjz_syst_JES1_pb_hi[i] + xjz_syst_JES2_pb_hi[i]*xjz_syst_JES2_pb_hi[i] + xjz_syst_JES2_pb_hi[i]*xjz_syst_JES2_pb_hi[i]) << endl;
  //  cout << xjz_syst_JES_pb[i] << endl;
  //}




  // Now make a README.... 

  cout << endl << endl;

  cout << "### make the .root workspace from the datacard: " << endl;
  cout << "# text2workspace.py combinedcard_jj.txt --X-allow-no-signal" << endl << endl;
  cout << "### to run the whole thing 'unfixed' (should return us the input values back)" << endl;
  cout << "# combine combinedcard.root -M MultiDimFit  -v 3 --freezeNuisances r --setPhysicsModelParameters r=0" << endl << endl;;
  cout << "### to fix the dx to be 1, this is like finding the best xhat (this is going to be close to the 'average' of pp and PbPb-bkg)" << endl;
  cout << "# combine combinedcard.root -M MultiDimFit  -v 3 --freezeNuisances r --setPhysicsModelParameters r=0,";
    for(int i=startbin; i<nBins+1; i++){
      if(i<nBins)       sprintf(saythis,"dx%d=1,",i);
      else if(i==nBins) sprintf(saythis,"dx%d=1",i);
      cout << saythis;
    }
    cout << " --redefineSignalPOI ";
    for(int i=startbin; i<nBins+1; i++){
      if(i<nBins)       sprintf(saythis,"dx%d,",i);
      else if(i==nBins) sprintf(saythis,"dx%d",i);
      cout << saythis;
    }
    cout << " --algo fixed" << endl << endl;
    cout << "### make the frozen workspace... similar to fixing the dx's" << endl;
    cout << "# combine combinedcard_jj.root -M MultiDimFit  -v 3  --setPhysicsModelParameters r=0,";
    for(int i=startbin; i<nBins+1; i++){
      if(i<nBins)       sprintf(saythis,"dx%d=1,",i);
      else if(i==nBins) sprintf(saythis,"dx%d=1",i);
      cout << saythis;
    }
    cout << " --freezeNuisances r,";
    for(int i=startbin; i<nBins+1; i++){
      if(i<nBins)       sprintf(saythis,"dx%d,",i);
      else if(i==nBins) sprintf(saythis,"dx%d",i);
      cout << saythis;
    }
    cout << " --saveWorkspace" << endl << endl;
    cout << "### run 10 toys" << endl;
    cout << "# combine newWS_jj.root  -M MultiDimFit  -v 3 --freezeNuisances r --setPhysicsModelParameters r=0,";
    for(int i=startbin; i<nBins+1; i++){
      if(i<nBins)       sprintf(saythis,"dx%d=1,",i);
      else if(i==nBins) sprintf(saythis,"dx%d=1",i);
      cout << saythis;
    }
    cout << " --redefineSignalPOI ";
    for(int i=startbin; i<nBins+1; i++){
      if(i<nBins)       sprintf(saythis,"dx%d,",i);
      else if(i==nBins) sprintf(saythis,"dx%d",i);
      cout << saythis;
    }
    cout << " --algo fixed --snapshotName MultiDimFit --forceRecreateNLL  -t 10  (include '--toysFrequentist'  if running 2dNLL w/ systematics. otherwise don't.) moreover, if you want to save the toys, include '--saveToys' " << endl << endl;

    cout << "### to set xhat to the pp values for a new saved workspace change the x1,x2,etc ranges to be the pp values. then they're locked." << endl;
    cout << "### run 10 toys after setting xhat to the pp values" << endl;
    cout << "# combine newWS_jj.root  -M MultiDimFit  -v 3 --freezeNuisances r --setPhysicsModelParameters r=0,";
    for(int i=startbin; i<nBins+1; i++){
      if(i<nBins)       sprintf(saythis,"dx%d=1,",i);
      else if(i==nBins) sprintf(saythis,"dx%d=1",i);
      cout << saythis;
    }
    cout << " --redefineSignalPOI ";
    for(int i=startbin; i<nBins+1; i++){
      if(i<nBins)       sprintf(saythis,"dx%d,",i);
      else if(i==nBins) sprintf(saythis,"dx%d",i);
      cout << saythis;
    }
    for(int i=startbin; i<nBins+1; i++){
      if(i<nBins)       sprintf(saythis,"x%d,",i);
      else if(i==nBins) sprintf(saythis,"x%d",i);
      cout << saythis;
    }
    cout << " --algo fixed --snapshotName MultiDimFit --forceRecreateNLL  -t 10";
    cout << " --setPhysicsModelParameterRanges ";
    for(int i=startbin; i<nBins+1; i++){
      if(i<nBins)       sprintf(saythis,"x%d=0,%2.2f:",i,10.0*h_pp->GetBinContent(h_pp->GetMaximumBin()));
      else if(i==nBins) sprintf(saythis,"x%d=0,%2.2f" ,i,10.0*h_pp->GetBinContent(h_pp->GetMaximumBin()));
      cout << saythis;
    }
    cout << endl;


/*

  for(int i=startbin; i<nBins+1; i++){
    //cout << int(floor( h1->GetBinContent(i)*h1->GetBinContent(i)/(h1->GetBinError(i)*h1->GetBinError(i))+0.5 ))/(h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i)) << endl;
    cout << (h1->GetBinContent(i)*h1->GetBinContent(i)/(h1->GetBinError(i)*h1->GetBinError(i)))/(h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i)) << endl;
  }

*/

  return;
}



void PrintSystUncertainty(double *percent_uncert, TH1 *h_data, int FirstBin, int LastBin, char *category, char *systemtype)
{

  char saythis[500];
  double adhoc_extra = 0.00;//0.01;

  sprintf(saythis,"syst%s_%s     lnN    ",category,systemtype);
  cout << saythis;
  for(int i=FirstBin; i<LastBin+1; i++){
    if(h_data->GetBinContent(i)>0 && strcmp(systemtype,"pp")==0)
      sprintf(saythis,"  %4.4f ",1.0+adhoc_extra+0.01*percent_uncert[i-1]);
    else
      sprintf(saythis,"  -     ");
    cout << saythis;
  }
  for(int i=FirstBin; i<LastBin+1; i++){
    if(h_data->GetBinContent(i)>0 && strcmp(systemtype,"Pb")==0)
      sprintf(saythis,"  %4.4f ",1.0+adhoc_extra+0.01*percent_uncert[i-1]);
    else
      sprintf(saythis,"  -     ");
    cout << saythis;
  }
  for(int i=FirstBin; i<LastBin+1; i++){
    if(h_data->GetBinContent(i)>0 && strcmp(systemtype,"Bg")==0)
      sprintf(saythis,"  %4.4f ",1.0+adhoc_extra+0.01*percent_uncert[i-1]);
    else
      sprintf(saythis,"  -     ");
    cout << saythis;
  }
  cout << endl;


  return;
}




void PrintSystUncertainty2(double *percent_uncert_lo, double *percent_uncert_hi, TH1 *h_data, int FirstBin, int LastBin, char *category, char *systemtype)
{

  char saythis[500];
  double adhoc_extra = 0.00;//0.01;

  sprintf(saythis,"syst%s_%s     lnN    ",category,systemtype);
  cout << saythis;
  for(int i=FirstBin; i<LastBin+1; i++){
    if(h_data->GetBinContent(i)>0 && strcmp(systemtype,"pp")==0)
      sprintf(saythis,"  %4.4f/%4.4f ",1.0+copysign(adhoc_extra,percent_uncert_lo[i-1])+0.01*percent_uncert_lo[i-1],1.0+copysign(adhoc_extra,percent_uncert_hi[i-1])+0.01*percent_uncert_hi[i-1]);
    else
      sprintf(saythis,"  -     ");
    cout << saythis;
  }
  for(int i=FirstBin; i<LastBin+1; i++){
    if(h_data->GetBinContent(i)>0 && strcmp(systemtype,"Pb")==0)
      sprintf(saythis,"  %4.4f/%4.4f ",1.0+adhoc_extra+0.01*percent_uncert_lo[i-1],1.0+adhoc_extra+0.01*percent_uncert_hi[i-1]);
    else
      sprintf(saythis,"  -     ");
    cout << saythis;
  }
  for(int i=FirstBin; i<LastBin+1; i++){
    if(h_data->GetBinContent(i)>0 && strcmp(systemtype,"Bg")==0)
      sprintf(saythis,"  %4.4f/%4.4f ",1.0+adhoc_extra+0.01*percent_uncert_lo[i-1],1.0+adhoc_extra+0.01*percent_uncert_hi[i-1]);
    else
      sprintf(saythis,"  -     ");
    cout << saythis;
  }
  cout << endl;


  return;
}

