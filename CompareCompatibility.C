// RooFit
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooAbsData.h"
#include "RooHist.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "Roo1DTable.h"
#include "RooConstVar.h"
#include "RooProduct.h"
#include "RooRandom.h"
#include "TStopwatch.h"
#include "RooNLLVar.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"

// RooStat
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileInspector.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/RooStatsUtils.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;



TH1D *deltaNLL_xjz  ;
TH1D *deltaNLL_dphi ;
TH1D *h_Chisq_dphi  ;
TH1D *h_Chisq_xjz   ;
TH1D *h_KolSm_dphi  ;
TH1D *h_KolSm_xjz   ;
TH1D *h_AndDr_dphi  ;
TH1D *h_AndDr_xjz   ;
TH1D *h_xjz_chi2_KS_AD  ;
TH1D *h_dphi_chi2_KS_AD ;

TH1D *xjz_raw_stat_pp ;
TH1D *xjz_raw_stat_pb ;
TH1D *xjz_bkg_stat_pb ;

TH1D *dphi_raw_stat_pp;
TH1D *dphi_raw_stat_pb;
TH1D *dphi_bkg_stat_pb;

TH1D *dphi_sum_stat_pp;
TH1D *dphi_sub_stat_pb;
TH1D *xjz_sum_stat_pp ;
TH1D *xjz_sub_stat_pb ;

TH1D *eff_pb_xjz;
TH1D *eff_pb_dphi;
TH1D *eff_pp_xjz;
TH1D *eff_pp_dphi;

TH1D *dphi_xhat;
TH1D *xjz_xhat;

TH1D *combinetoys_dphi[20];
TH1D *combinetoys_xjz[20];

double q_CS_xjz  = 20.7958;
double q_CS_dphi = 14.2788;
double q_KS_xjz  = 1.0;
double q_KS_dphi = 1.0;
double q_AD_xjz  = 1.0;
double q_AD_dphi = 1.0;
double q_HC_xjz  = 17.23;// these get reset later on !
double q_HC_dphi = 11.71;// these get reset later on !

double p_CS_xjz  = 0.0;
double p_CS_dphi = 0.0;
double p_KS_xjz  = 0.0;
double p_KS_dphi = 0.0;
double p_AD_xjz  = 0.0;
double p_AD_dphi = 0.0;
double p_HC_xjz  = 0.0;
double p_HC_dphi = 0.0;

double xhat_dphi[8] = {
  0.0342393,
  0.0533977,
  0.0694632,
  0.142984 ,
  0.259036 ,
  0.543947 ,
  1.15056	 ,
  2.33174	 
};
double xhat_xjz[16] = {
  1.80288e-12,
  0.108666	 ,
  0.238436	 ,
  0.461928	 ,
  0.69946	   ,
  0.684536	 ,
  0.74338	   ,
  0.660996	 ,
  0.496388	 ,
  0.270276	 ,
  0.155466	 ,
  0.10244	   ,
  0.0410323	 ,
  0.0356612	 ,
  0.0237741	 ,
  0.0237741	 
};

double exhat_dphi[8] = {
  0.008756,
  0.010994,
  0.014591,
  0.020697,
  0.032494,
  0.046624,
  0.093693,
  0.134107
};
double exhat_xjz[16] = {
  0.0059419,
  0.0347425,
  0.0525188,
  0.0689161,
  0.0821144,
  0.0793953,
  0.0829732,
  0.0779541,
  0.0675646,
  0.0470781,
  0.0357688,
  0.0324704,
  0.0183484,
  0.0205836,
  0.0168064,
  0.0168064
};

void PutToyInHisto(RooDataSet *rooToy, TH1D *htoy, TH1D *eff, int ch);

void CompareCompatibility()
{

  char saythis[500];

  char infilename[500];
  //sprintf(infilename,"JasonToyMC_norm1_ppflucts1_nExp1000_Error1.00_RenormalizeToy0_normpp2pb0_xjzcut_2_12.root");
  //sprintf(infilename,"Fall2016/JasonToyMC_31Aug_norm1_ppflucts1_nExp1000_Error1.00_RenormalizeToy0_normpp2pb0_xjzcut_2_12.root");
  //sprintf(infilename,"Jan2017/JasonToyMC_31Jan_norm1_ppflucts1_nExp1000_Error1.00_RenormalizeToy0_normpp2pb0_xjzcut_1_9.root");
  sprintf(infilename,"Jan2017/JasonToyMC_31Jan_norm1_ppflucts1_nExp1000_Error1.00_RenormalizeToy0_normpp2pb0_xjzcut_1_8.root");
  TFile *f_JasonToyMC   = TFile::Open(infilename);
  //TFile *f_JasonToyMC   = TFile::Open("JasonToyMC_norm1_ppflucts1_nExp1000_Error1.00_RenormalizeToy1_normpp2pb0_xjzcut_2_12.root");
  //TFile *f_JasonToyMC   = TFile::Open("JasonToyMC_norm1_ppflucts1_nExp1000_Error1.00_RenormalizeToy0_normpp2pb1_xjzcut_2_12.root");
  

  //TFile *f_combine_toy_xjz  = TFile::Open("higgsCombineTest.MultiDimFit.fitxhat.xjz_cut_2_12.limits.123456.root");
  //  TFile *f_combine_toy_xjz  = TFile::Open("higgsCombineTest.MultiDimFit.ppxhat.xjz.limitsdx.seeds.root");
  //TFile *f_combine_toy_xjz  = TFile::Open("higgsCombineTest.MultiDimFit.ppxhat2.xjz_cut_2_12.6051981.root");
  //TFile *f_combine_val_xjz  = TFile::Open("higgsCombineTest.MultiDimFit.fitxhat.xjz_cut_2_12.limits.root");
  //TFile *f_combine_val_xjz  = TFile::Open("higgsCombineTest.MultiDimFit.fitxhat.xjz.limitsdx.root");
  //  TFile *f_combine_val_xjz  = TFile::Open("higgsCombineTest.MultiDimFit.ppxhat.xjz_cut_2_12.root");
  //TFile *f_combine_toy_xjz  = TFile::Open("Fall2016/higgsCombineTest.MultiDimFit.fitxhat.xjz_cut_2_12.1000toys.root");
  //TFile *f_combine_val_xjz  = TFile::Open("Fall2016/higgsCombineTest.MultiDimFit.fitxhat.xjz_cut_2_12.root");
  //TFile *f_combine_toy_xjz  = TFile::Open("Jan2017/higgsCombineTest.MultiDimFit.fitxhat.xjz_cut_1_8.1000toys.root");
  //TFile *f_combine_val_xjz  = TFile::Open("Jan2017/higgsCombineTest.MultiDimFit.fitxhat.xjz_cut_1_8.root");
  TFile *f_combine_toy_xjz  = TFile::Open("Jan2017/higgsCombineTest.MultiDimFit.fitxhat.xjz.1000toys.root");
  TFile *f_combine_val_xjz  = TFile::Open("Jan2017/higgsCombineTest.MultiDimFit.fitxhat.xjz.root");
  //TFile *f_combine_toy_xjz  = TFile::Open("Jan2017/higgsCombineTest.MultiDimFit.fitxhat.xjz_frankenstein.1000toys.root");
  //TFile *f_combine_val_xjz  = TFile::Open("Jan2017/higgsCombineTest.MultiDimFit.fitxhat.xjz_frankenstein.root");

  //TFile *f_combine_toy_xjz  = TFile::Open("higgsCombineTest.MultiDimFit.fitxhat.xjz.limitsdx.123456.root");
  //TFile *f_combine_val_xjz  = TFile::Open("higgsCombineTest.MultiDimFit.fitxhat.xjz.limitsdx.root");

  //TFile *f_combine_toy_dphi = TFile::Open("higgsCombineTest.MultiDimFit.fitxhat.dphi.limitsdx.123456.root");
  //  TFile *f_combine_toy_dphi = TFile::Open("higgsCombineTest.MultiDimFit.ppxhat.dphi.limitsdx.123456.root");
  //TFile *f_combine_val_dphi = TFile::Open("higgsCombineTest.MultiDimFit.fitxhat.dphi.limitsdx.root");
  //  TFile *f_combine_val_dphi = TFile::Open("higgsCombineTest.MultiDimFit.ppxhat.dphi.root");
  //TFile *f_combine_toy_dphi = TFile::Open("Fall2016/higgsCombineTest.MultiDimFit.fitxhat.dphi.1000toys.root");
  //TFile *f_combine_val_dphi = TFile::Open("Fall2016/higgsCombineTest.MultiDimFit.fitxhat.dphi.root");
  TFile *f_combine_toy_dphi = TFile::Open("Jan2017/higgsCombineTest.MultiDimFit.fitxhat.dphi.1000toys.root");
  TFile *f_combine_val_dphi = TFile::Open("Jan2017/higgsCombineTest.MultiDimFit.fitxhat.dphi.root");


  TF1 *jg  = new TF1("jg", "gaus(0)",0,5);
  jg->SetLineWidth(1);
  jg->SetParameters(10,10,5);
  TF1 *f_chisq = new TF1("f_chisq","[1]*(pow(x,[0]/2-1)*TMath::Exp(-x/2))/(pow(2,[0]/2)*TMath::Gamma([0]/2))",0,50);
  f_chisq->SetParameter(1,1000);
  f_chisq->SetParameter(0,8);
  f_chisq->SetLineColor(1);
  f_chisq->SetLineWidth(2);
  f_chisq->SetLineStyle(2);

  ((TTree*)f_combine_toy_xjz->Get("limit"))->Draw("2*deltaNLL>>htemp(200,0,100)","quantileExpected>-0.5");
  deltaNLL_xjz = (TH1D*)gPad->GetPrimitive("htemp")->Clone("deltaNLL_xjz");
  ((TTree*)f_combine_toy_dphi->Get("limit"))->Draw("2*deltaNLL>>htemp(100,0,50)","quantileExpected>-0.5");
  deltaNLL_dphi = (TH1D*)gPad->GetPrimitive("htemp")->Clone("deltaNLL_dphi");
  ((TTree*)f_combine_val_xjz ->Get("limit"))->Draw("2*deltaNLL>>htemp(200,0,100)","quantileExpected>-0.5");
  q_HC_xjz  = ((TH1D*)gPad->GetPrimitive("htemp"))->GetMean();
  ((TTree*)f_combine_val_dphi->Get("limit"))->Draw("2*deltaNLL>>htemp(100,0,50)","quantileExpected>-0.5");
  q_HC_dphi = ((TH1D*)gPad->GetPrimitive("htemp"))->GetMean();


  //TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_5Jul2016.root");
  //TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSumPP_5Jul2016.root");
  //TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_HI_PP_31Aug2016.root");
  //TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_HI_PP_31Aug2016.root");
  //TFile *f_pb = TFile::Open("$CODE/ZJet/Frankenstein.root");
  //TFile *f_pp = TFile::Open("$CODE/ZJet/Frankenstein.root");
  TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/Jan23_2017/zJetHistogramSum_HI_PP_20170118.root");
  TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/Jan23_2017/zJetHistogramSum_HI_PP_20170118.root");  

  xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_binJER_ptBin0_hiBin1_jetRAW_final_norm");
  xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_binJER_ptBin0_hiBin1_jetBKG_final_norm");
  xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get("h1D_xjz_binJER_ptBin0_hiBin0_jetRAW_final_norm");
  //xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("xjz_raw_stat_pb");
  //xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("xjz_bkg_stat_pb");
  //xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get("xjz_raw_stat_pp");

  dphi_raw_stat_pb = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_dphi_rebin_ptBin0_hiBin1_jetRAW_final_norm");
  dphi_bkg_stat_pb = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_dphi_rebin_ptBin0_hiBin1_jetBKG_final_norm");
  dphi_raw_stat_pp = (TH1D*)f_pp->GetDirectory("PP")->Get("h1D_dphi_rebin_ptBin0_hiBin0_jetRAW_final_norm");

  //xjz_raw_stat_pp  = (TH1D*)f_JasonToyMC  ->Get("xjz_raw_stat_pp");
  //xjz_raw_stat_pb  = (TH1D*)f_JasonToyMC  ->Get("xjz_raw_stat_pb");
  //xjz_bkg_stat_pb  = (TH1D*)f_JasonToyMC  ->Get("xjz_bkg_stat_pb");

  //dphi_raw_stat_pp = (TH1D*)f_JasonToyMC  ->Get("dphi_raw_stat_pp");
  //dphi_raw_stat_pb = (TH1D*)f_JasonToyMC  ->Get("dphi_raw_stat_pb");
  //dphi_bkg_stat_pb = (TH1D*)f_JasonToyMC  ->Get("dphi_bkg_stat_pb");

  dphi_sum_stat_pp  = (TH1D*)f_JasonToyMC  ->Get("dphi_sum_stat_pp");
  xjz_sum_stat_pp   = (TH1D*)f_JasonToyMC  ->Get("xjz_sum_stat_pp");
  h_Chisq_dphi      = (TH1D*)f_JasonToyMC  ->Get("h_Chisq_dphi");
  h_Chisq_xjz       = (TH1D*)f_JasonToyMC  ->Get("h_Chisq_xjz");
  h_KolSm_dphi      = (TH1D*)f_JasonToyMC  ->Get("h_KolSm_dphi");
  h_KolSm_xjz       = (TH1D*)f_JasonToyMC  ->Get("h_KolSm_xjz");
  h_AndDr_dphi      = (TH1D*)f_JasonToyMC  ->Get("h_AndDr_dphi");
  h_AndDr_xjz       = (TH1D*)f_JasonToyMC  ->Get("h_AndDr_xjz");
  h_xjz_chi2_KS_AD  = (TH1D*)f_JasonToyMC  ->Get("h_xjz_chi2_KS_AD");
  h_dphi_chi2_KS_AD = (TH1D*)f_JasonToyMC  ->Get("h_dphi_chi2_KS_AD");

  dphi_sub_stat_pb  = (TH1D*)dphi_raw_stat_pb->Clone("dphi_sub_stat_pb");
  xjz_sub_stat_pb   = (TH1D*)xjz_raw_stat_pb ->Clone("xjz_sub_stat_pb");
  dphi_sub_stat_pb  ->Add(dphi_bkg_stat_pb,-1);
  xjz_sub_stat_pb   ->Add(xjz_bkg_stat_pb, -1);


  TFile *f_combine_toys_ex_xjz  = TFile::Open("Summer2016/higgsCombineTest.MultiDimFit.fitxhat.xjz.limitsdx.50saveToys.123456.root");
  TFile *f_combine_toys_ex_dphi = TFile::Open("Summer2016/higgsCombineTest.MultiDimFit.fitxhat.dphi.limitsdx.50saveToys.123456.root");
  TDirectory *mytoysdir_xjz  = (TDirectory*)f_combine_toys_ex_xjz->GetDirectory("toys");
  TDirectory *mytoysdir_dphi = (TDirectory*)f_combine_toys_ex_dphi->GetDirectory("toys");

  eff_pb_xjz  = (TH1D*)xjz_raw_stat_pb ->Clone("eff_pb_xjz");
  for(int i=1; i<eff_pb_xjz->GetNbinsX()+1; i++){
    eff_pb_xjz->SetBinContent(i,eff_pb_xjz->GetBinError(i)*eff_pb_xjz->GetBinError(i)/eff_pb_xjz->GetBinContent(i));
    eff_pb_xjz->SetBinError  (i,0.0);
  }
  eff_pb_dphi = (TH1D*)dphi_raw_stat_pb->Clone("eff_pb_dphi");
  for(int i=1; i<eff_pb_dphi->GetNbinsX()+1; i++){
    eff_pb_dphi->SetBinContent(i,eff_pb_dphi->GetBinError(i)*eff_pb_dphi->GetBinError(i)/eff_pb_dphi->GetBinContent(i));
    eff_pb_dphi->SetBinError  (i,0.0);
  }
  eff_pp_xjz  = (TH1D*)xjz_sum_stat_pp ->Clone("eff_pp_xjz");
  for(int i=1; i<eff_pp_xjz->GetNbinsX()+1; i++){
    eff_pp_xjz->SetBinContent(i,eff_pp_xjz->GetBinError(i)*eff_pp_xjz->GetBinError(i)/eff_pp_xjz->GetBinContent(i));
    eff_pp_xjz->SetBinError  (i,0.0);
  }
  eff_pp_dphi = (TH1D*)dphi_sum_stat_pp->Clone("eff_pp_dphi");
  for(int i=1; i<eff_pp_dphi->GetNbinsX()+1; i++){
    eff_pp_dphi->SetBinContent(i,eff_pp_dphi->GetBinError(i)*eff_pp_dphi->GetBinError(i)/eff_pp_dphi->GetBinContent(i));
    eff_pp_dphi->SetBinError  (i,0.0);
  }

  xjz_raw_stat_pb  ->SetMarkerColor(1);
  xjz_bkg_stat_pb  ->SetMarkerColor(1);
  xjz_raw_stat_pp  ->SetMarkerColor(4);

  dphi_raw_stat_pb ->SetMarkerColor(1);
  dphi_bkg_stat_pb ->SetMarkerColor(1);
  dphi_raw_stat_pp ->SetMarkerColor(4);

  xjz_raw_stat_pb  ->SetLineColor  (1);
  xjz_bkg_stat_pb  ->SetLineColor  (1);
  xjz_raw_stat_pp  ->SetLineColor  (4);

  dphi_raw_stat_pb ->SetLineColor  (1);
  dphi_bkg_stat_pb ->SetLineColor  (1);
  dphi_raw_stat_pp ->SetLineColor  (4);


  TCanvas *c3 = new TCanvas("c3","c3");
  c3->Divide(5,4);
  TCanvas *c4 = new TCanvas("c4","c4");
  c4->Divide(5,4);

  for(int i=0; i<20; i++){
    sprintf(saythis,"toy_%d",i+1);
    RooDataSet *rooToy_dphi = (RooDataSet*)mytoysdir_dphi->Get(saythis);
    sprintf(saythis,"combinetoys_dphi_%d",i);
    combinetoys_dphi[i] = (TH1D*)dphi_raw_stat_pb->Clone(saythis);
    combinetoys_dphi[i]->Reset();
    PutToyInHisto(rooToy_dphi, combinetoys_dphi[i], eff_pp_dphi, 1);

    sprintf(saythis,"toy_%d",i+1);
    RooDataSet *rooToy_xjz = (RooDataSet*)mytoysdir_xjz->Get(saythis);
    sprintf(saythis,"combinetoys_xjz_%d",i);
    combinetoys_xjz[i] = (TH1D*)xjz_raw_stat_pb->Clone(saythis);
    combinetoys_xjz[i]->Reset();
    PutToyInHisto(rooToy_xjz, combinetoys_xjz[i], eff_pp_xjz, 1);

    c3->cd(i+1);
    xjz_raw_stat_pb->SetMarkerColor(19);
    xjz_raw_stat_pb->SetLineColor  (19);
    xjz_raw_stat_pb->DrawCopy();
    xjz_raw_stat_pb->SetMarkerColor(1);
    xjz_raw_stat_pb->SetLineColor  (1);
    xjz_sum_stat_pp->DrawCopy("same");
    combinetoys_xjz[i]->DrawCopy("same");

    c4->cd(i+1);
    dphi_raw_stat_pb->SetMarkerColor(19);
    dphi_raw_stat_pb->SetLineColor  (19);
    dphi_raw_stat_pb->DrawCopy();
    dphi_raw_stat_pb->SetMarkerColor(1);
    dphi_raw_stat_pb->SetLineColor  (1);
    dphi_sum_stat_pp->DrawCopy("same");
    combinetoys_dphi[i]->DrawCopy("same");

  }





  dphi_xhat  = (TH1D*)dphi_raw_stat_pb->Clone("dphi_xhat");
  xjz_xhat   = (TH1D*)xjz_raw_stat_pb ->Clone("xjz_xhat");
  dphi_xhat  ->Reset();
  xjz_xhat   ->Reset();
  dphi_xhat  ->SetMarkerColor(2);
  xjz_xhat   ->SetMarkerColor(2);
  dphi_xhat  ->SetLineColor  (2);
  xjz_xhat   ->SetLineColor  (2);
  dphi_xhat  ->SetMarkerStyle(24);
  xjz_xhat   ->SetMarkerStyle(24);
  for(int i=0; i<dphi_xhat->GetNbinsX(); i++){
    dphi_xhat->SetBinContent(i+1,xhat_dphi[i]);
    dphi_xhat->SetBinError  (i+1,exhat_dphi[i]);
  }
  for(int i=0; i<xjz_xhat->GetNbinsX(); i++){
    xjz_xhat->SetBinContent(i+1,xhat_xjz[i]);
    xjz_xhat->SetBinError  (i+1,exhat_xjz[i]);
  }


  q_CS_xjz  = h_xjz_chi2_KS_AD ->GetBinContent(1);
  q_CS_dphi = h_dphi_chi2_KS_AD->GetBinContent(1);
  q_KS_xjz  = h_xjz_chi2_KS_AD ->GetBinContent(2);
  q_KS_dphi = h_dphi_chi2_KS_AD->GetBinContent(2);
  q_AD_xjz  = h_xjz_chi2_KS_AD ->GetBinContent(3);
  q_AD_dphi = h_dphi_chi2_KS_AD->GetBinContent(3);



  p_CS_xjz  = h_Chisq_xjz  ->Integral(h_Chisq_xjz  ->GetXaxis()->FindBin(q_CS_xjz), h_Chisq_xjz  ->GetNbinsX())/h_Chisq_xjz  ->Integral(1,h_Chisq_xjz  ->GetNbinsX());
  p_CS_dphi = h_Chisq_dphi ->Integral(h_Chisq_dphi ->GetXaxis()->FindBin(q_CS_dphi),h_Chisq_dphi ->GetNbinsX())/h_Chisq_dphi ->Integral(1,h_Chisq_dphi ->GetNbinsX());
  p_KS_xjz  = h_KolSm_xjz  ->Integral(h_KolSm_xjz  ->GetXaxis()->FindBin(q_KS_xjz), h_KolSm_xjz  ->GetNbinsX())/h_KolSm_xjz  ->Integral(1,h_KolSm_xjz  ->GetNbinsX());
  p_KS_dphi = h_KolSm_dphi ->Integral(h_KolSm_dphi ->GetXaxis()->FindBin(q_KS_dphi),h_KolSm_dphi ->GetNbinsX())/h_KolSm_dphi ->Integral(1,h_KolSm_dphi ->GetNbinsX());
  p_AD_xjz  = h_AndDr_xjz  ->Integral(h_AndDr_xjz  ->GetXaxis()->FindBin(q_AD_xjz), h_AndDr_xjz  ->GetNbinsX())/h_AndDr_xjz  ->Integral(1,h_AndDr_xjz  ->GetNbinsX());
  p_AD_dphi = h_AndDr_dphi ->Integral(h_AndDr_dphi ->GetXaxis()->FindBin(q_AD_dphi),h_AndDr_dphi ->GetNbinsX())/h_AndDr_dphi ->Integral(1,h_AndDr_dphi ->GetNbinsX());
  p_HC_xjz  = deltaNLL_xjz ->Integral(deltaNLL_xjz ->GetXaxis()->FindBin(q_HC_xjz), deltaNLL_xjz ->GetNbinsX())/deltaNLL_xjz ->Integral(1,deltaNLL_xjz ->GetNbinsX());
  p_HC_dphi = deltaNLL_dphi->Integral(deltaNLL_dphi->GetXaxis()->FindBin(q_HC_dphi),deltaNLL_dphi->GetNbinsX())/deltaNLL_dphi->Integral(1,deltaNLL_dphi->GetNbinsX());


  deltaNLL_xjz  ->Rebin(2);
  deltaNLL_dphi ->Rebin(2);
  h_Chisq_dphi  ->Rebin(2);
  h_Chisq_xjz   ->Rebin(2);

  deltaNLL_xjz  ->SetLineColor(2);
  deltaNLL_dphi ->SetLineColor(2);
  h_Chisq_dphi  ->SetLineColor(4);
  h_Chisq_xjz   ->SetLineColor(4);

  TLine *line_q_CS_xjz  = new TLine(q_CS_xjz ,0,q_CS_xjz ,h_Chisq_xjz ->GetMaximum());
  TLine *line_q_CS_dphi = new TLine(q_CS_dphi,0,q_CS_dphi,h_Chisq_dphi->GetMaximum());
  TLine *line_q_HC_xjz  = new TLine(q_HC_xjz ,0,q_HC_xjz ,h_Chisq_xjz ->GetMaximum());
  TLine *line_q_HC_dphi = new TLine(q_HC_dphi,0,q_HC_dphi,h_Chisq_dphi->GetMaximum());

  line_q_CS_xjz  ->SetLineWidth(2);
  line_q_CS_dphi ->SetLineWidth(2);
  line_q_HC_xjz  ->SetLineWidth(2);
  line_q_HC_dphi ->SetLineWidth(2);

  line_q_CS_xjz  ->SetLineStyle(2);
  line_q_CS_dphi ->SetLineStyle(2);
  line_q_HC_xjz  ->SetLineStyle(2);
  line_q_HC_dphi ->SetLineStyle(2);

  line_q_CS_xjz  ->SetLineColor(4);
  line_q_CS_dphi ->SetLineColor(4);
  line_q_HC_xjz  ->SetLineColor(2);
  line_q_HC_dphi ->SetLineColor(2);


  gStyle->SetOptStat(0);
  TH1D *h_frame = new TH1D("h_frame"," ; q (compatibility)",200,0,100);
  h_frame->GetXaxis()->SetRangeUser(0,50);

  TLegend *leg1a = new TLegend(0.52,0.3,0.89,0.89);
  //sprintf(saythis,"X_{JZ}  %d toys",int(deltaNLL_xjz->GetEntries()));
  //leg1a->AddEntry(deltaNLL_xjz,saythis,"");
  leg1a->AddEntry(deltaNLL_xjz,"2#DeltaNLL (combine)","L");
  sprintf(saythis,"toy mean: %2.2f",deltaNLL_xjz->GetMean());
  leg1a->AddEntry(deltaNLL_xjz,saythis,"");
  jg->SetLineColor(46);
  deltaNLL_xjz->Fit("jg","QN","",0.25*deltaNLL_xjz->GetMean(),1.5*deltaNLL_xjz->GetMean());
  sprintf(saythis,"toy MPV: %2.2f",jg->GetParameter(1));
  leg1a->AddEntry(deltaNLL_xjz,saythis,"");
  sprintf(saythis,"data 2#DeltaNLL: %2.2f",q_HC_xjz);
  leg1a->AddEntry(deltaNLL_xjz,saythis,"");
  sprintf(saythis,"P-Val: %2.3f",p_HC_xjz);
  leg1a->AddEntry(deltaNLL_xjz,saythis,"");
  leg1a->AddEntry(line_q_HC_xjz,"data","L");

  leg1a->AddEntry(deltaNLL_xjz," ","");
  leg1a->AddEntry(h_Chisq_xjz,"#chi^{2} (standalone)","L");
  sprintf(saythis,"toy mean: %2.2f",h_Chisq_xjz->GetMean());
  leg1a->AddEntry(h_Chisq_xjz,saythis,"");
  jg->SetLineColor(38);
  h_Chisq_xjz->Fit("jg","QN","",0.25*h_Chisq_xjz->GetMean(),1.5*h_Chisq_xjz->GetMean());
  sprintf(saythis,"toy MPV: %2.2f",jg->GetParameter(1));
  leg1a->AddEntry(h_Chisq_xjz,saythis,"");
  sprintf(saythis,"data #chi^{2}: %2.2f",q_CS_xjz);
  leg1a->AddEntry(h_Chisq_xjz,saythis,"");
  sprintf(saythis,"P-Val: %2.2f",p_CS_xjz);
  leg1a->AddEntry(h_Chisq_xjz,saythis,"");
  leg1a->AddEntry(line_q_CS_xjz,"data","L");

  TLegend *leg1b = new TLegend(0.52,0.3,0.89,0.89);
  //sprintf(saythis,"#Delta#varphi_{JZ}  %d toys",int(deltaNLL_dphi->GetEntries()));
  //leg1b->AddEntry(deltaNLL_dphi,saythis,"");
  leg1b->AddEntry(deltaNLL_dphi,"2#DeltaNLL (combine)","L");
  sprintf(saythis,"toy mean: %2.2f",deltaNLL_dphi->GetMean());
  leg1b->AddEntry(deltaNLL_dphi,saythis,"");
  jg->SetLineColor(46);
  deltaNLL_dphi->Fit("jg","QN","",0.25*deltaNLL_dphi->GetMean(),1.5*deltaNLL_dphi->GetMean());
  sprintf(saythis,"toy MPV: %2.2f",jg->GetParameter(1));
  leg1b->AddEntry(deltaNLL_dphi,saythis,"");
  sprintf(saythis,"data 2#DeltaNLL: %2.2f",q_HC_dphi);
  leg1b->AddEntry(h_Chisq_dphi,saythis,"");
  sprintf(saythis,"P-Val: %2.2f",p_HC_dphi);
  leg1b->AddEntry(deltaNLL_dphi,saythis,"");
  leg1b->AddEntry(line_q_HC_dphi,"data","L");

  leg1b->AddEntry(deltaNLL_dphi," ","");
  leg1b->AddEntry(h_Chisq_dphi,"#chi^{2} (standalone)","L");
  sprintf(saythis,"toy mean: %2.2f",h_Chisq_dphi->GetMean());
  leg1b->AddEntry(h_Chisq_dphi,saythis,"");
  jg->SetLineColor(38);
  h_Chisq_dphi->Fit("jg","QN","",0.25*h_Chisq_dphi->GetMean(),1.5*h_Chisq_dphi->GetMean());
  sprintf(saythis,"toy MPV: %2.2f",jg->GetParameter(1));
  leg1b->AddEntry(h_Chisq_dphi,saythis,"");
  sprintf(saythis,"data #chi^{2}: %2.2f",q_CS_dphi);
  leg1b->AddEntry(h_Chisq_dphi,saythis,"");
  sprintf(saythis,"P-Val: %2.2f",p_CS_dphi);
  leg1b->AddEntry(h_Chisq_dphi,saythis,"");
  leg1b->AddEntry(line_q_CS_dphi,"data","L");

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1",400,800);
  c1->Divide(1,2);
  c1->cd(1);
  //h_frame->GetYaxis()->SetRangeUser(0,1.2*deltaNLL_xjz->GetMaximum());
  h_frame->GetYaxis()->SetRangeUser(0,1.15*deltaNLL_xjz->GetMaximum());
  //h_frame->GetXaxis()->SetRangeUser(0,100);
  sprintf(saythis,"X_{JZ}  %d toys",int(deltaNLL_xjz->GetEntries()));
  h_frame->SetTitle(saythis);
  h_frame->DrawCopy();
  h_frame->GetXaxis()->SetRangeUser(0,50);
  //deltaNLL_xjz->Scale(1.0/0.75);
  deltaNLL_xjz->DrawCopy("same");
  h_Chisq_xjz ->DrawCopy("same");
  line_q_CS_xjz  ->Draw("same");
  line_q_HC_xjz  ->Draw("same");
  //f_chisq->SetParameter(0,26);//Frankenstein
  f_chisq->SetParameter(0,9);//Kaya reduced the number of bins.
  //f_chisq->SetParameter(0,11);
  f_chisq->DrawCopy("same");
  sprintf(saythis,"ideal #chi^{2},  %d dof",int(f_chisq->GetParameter(0)+0.1));
  leg1a->AddEntry(f_chisq," ","");
  leg1a->AddEntry(f_chisq,saythis,"L");
  leg1a->Draw();

  c1->cd(2);
  //h_frame->GetYaxis()->SetRangeUser(0,1.2*deltaNLL_dphi->GetMaximum());
  h_frame->GetYaxis()->SetRangeUser(0,140);
  sprintf(saythis,"#Delta#varphi_{JZ}  %d toys",int(deltaNLL_xjz->GetEntries()));
  h_frame->SetTitle(saythis);
  h_frame->DrawCopy();
  deltaNLL_dphi->DrawCopy("same");
  h_Chisq_dphi ->DrawCopy("same");
  line_q_CS_dphi ->Draw("same");
  line_q_HC_dphi ->Draw("same");
  f_chisq->SetParameter(0,8);
  f_chisq->DrawCopy("same");
  sprintf(saythis,"ideal #chi^{2},  %d dof",int(f_chisq->GetParameter(0)+0.1));
  leg1b->AddEntry(f_chisq," ","");
  leg1b->AddEntry(f_chisq,saythis,"L");
  leg1b->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",800,400);
  c2->Divide(2,1);
  c2->cd(1);
  xjz_sub_stat_pb->Draw();
  xjz_raw_stat_pp->Draw("same");
  xjz_xhat       ->Draw("same");
  c2->cd(2);
  dphi_sub_stat_pb->Draw();
  dphi_raw_stat_pp->Draw("same");
  dphi_xhat       ->Draw("same");

  TLegend *leg2 = new TLegend(0.14,0.65,0.61,0.88);
  leg2->AddEntry(xjz_sub_stat_pb,"PbPb-bkg","P");
  leg2->AddEntry(xjz_raw_stat_pp,"pp","P");
  leg2->AddEntry(xjz_xhat       ,"best fit #hat{x}","P");
  c2->cd(2);
  leg2->Draw();


  cout << endl;
  cout << "using : " << infilename << endl;
  cout << endl;
  cout << "___ DPHI ___" << endl;
  sprintf(saythis,"2dNLL -->  q= %2.2f  pVal= %2.2f",q_HC_dphi,p_HC_dphi);
  cout << saythis << endl;
  sprintf(saythis,"ChiSq -->  q= %2.2f  pVal= %2.2f",q_CS_dphi,p_CS_dphi);
  cout << saythis << endl;
  sprintf(saythis,"KolSm -->  q= %2.2f  pVal= %2.2f",q_KS_dphi,p_KS_dphi);
  cout << saythis << endl << endl;

  cout << "___ XJZ  ___" << endl;
  sprintf(saythis,"2dNLL -->  q= %2.2f  pVal= %2.2f",q_HC_xjz,p_HC_xjz);
  cout << saythis << endl;
  sprintf(saythis,"ChiSq -->  q= %2.2f  pVal= %2.2f",q_CS_xjz,p_CS_xjz);
  cout << saythis << endl;
  sprintf(saythis,"KolSm -->  q= %2.2f  pVal= %2.2f",q_KS_xjz,p_KS_xjz);
  cout << saythis << endl << endl;


}



void PutToyInHisto(RooDataSet *rooToy, TH1D *htoy, TH1D *eff, int ch)
{
  char saythis[500];
  double ValueCounts = 0.0;

  for(int i=0; i<eff->GetNbinsX(); i++){
    sprintf(saythis,"n_obs_binch%d_b%d",ch,i+1);
    ValueCounts = ((RooRealVar*)rooToy->get(0)->find(saythis))->getVal();
    //cout << saythis << "  ~~  " << ValueCounts << endl;
    htoy->SetBinContent(i+1,ValueCounts);
    htoy->SetBinError  (i+1,TMath::Sqrt(ValueCounts));
  }
  htoy->Multiply(eff);

  return;
}
