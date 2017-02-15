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

TH1D *h_Chisq_dphi;
TH1D *h_KolSm_dphi;
TH1D *h_AndDr_dphi;
TH1D *h_Chisq_xjz;
TH1D *h_KolSm_xjz;
TH1D *h_AndDr_xjz;
TH1D* h_KoloHist[2];

TH1D *h_xjz_chi2_KS_AD;
TH1D *h_dphi_chi2_KS_AD;

#include "CompatibilityHelperFunctions.C"


void UsePPAsComparison()
{

  int writeoutputfile = 1;

  int normalize_pp_or_pbpb    = 1;
  int include_ppflucts_in_toy = 1;
  int renormalize_toy         = 0;

  int normalize_pp2pb         = 0;
  int nPseudoExp = 1000;
  double BlowUpErrors = 1.0;
  int FirstBin_xjz = 1;//2
  int LastBin_xjz  = 8;//12
  char saythis[500];
  float yaxis_dphi[2] = {0.0,0.0};
  float yaxis_xjz [2] = {0.0,0.0};

  if(normalize_pp_or_pbpb==0){
    yaxis_dphi[0] = 0.0;
    yaxis_xjz [0] = 0.0;
    yaxis_dphi[1] = 1200;
    yaxis_xjz [1] = 250;
  }
  else if(normalize_pp_or_pbpb==1){
    yaxis_dphi[0] = 0.0;
    yaxis_xjz [0] = 0.0;
    yaxis_dphi[1] = 100;
    yaxis_xjz [1] = 40;
  }

  h_xjz_chi2_KS_AD  = new TH1D("h_xjz_chi2_KS_AD" ,"h_xjz_chi2_KS_AD" ,3,0.5,3.5);
  h_dphi_chi2_KS_AD = new TH1D("h_dphi_chi2_KS_AD","h_dphi_chi2_KS_AD",3,0.5,3.5);

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

  //xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("xjz_raw_stat_pb");
  //xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("xjz_bkg_stat_pb");
  //xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get("xjz_raw_stat_pp");
  xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_binJER_ptBin0_hiBin1_jetRAW_final_norm");
  xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_binJER_ptBin0_hiBin1_jetBKG_final_norm");
  xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get("h1D_xjz_binJER_ptBin0_hiBin0_jetRAW_final_norm");
  //xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_ptBin0_hiBin1_jetRAW_final_norm");
  //xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_ptBin0_hiBin1_jetBKG_final_norm");
  //xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get("h1D_xjz_ptBin0_hiBin0_jetRAW_final_norm");

  dphi_raw_stat_pb = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_dphi_rebin_ptBin0_hiBin1_jetRAW_final_norm");
  dphi_bkg_stat_pb = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_dphi_rebin_ptBin0_hiBin1_jetBKG_final_norm");
  dphi_raw_stat_pp = (TH1D*)f_pp->GetDirectory("PP")->Get("h1D_dphi_rebin_ptBin0_hiBin0_jetRAW_final_norm");

  for(int i=1; i<xjz_raw_stat_pp->GetNbinsX()+1; i++){
    xjz_raw_stat_pp->SetBinError(i,BlowUpErrors*xjz_raw_stat_pp->GetBinError(i));
    xjz_bkg_stat_pb->SetBinError(i,BlowUpErrors*xjz_bkg_stat_pb->GetBinError(i));

    //if(i==1 || i>13){
    //xjz_raw_stat_pp->SetBinContent(i,0);
    //xjz_raw_stat_pb->SetBinContent(i,0);
    //xjz_bkg_stat_pb->SetBinContent(i,0);
    //xjz_raw_stat_pp->SetBinError  (i,0);
    //xjz_raw_stat_pb->SetBinError  (i,0);
    //xjz_bkg_stat_pb->SetBinError  (i,0);
    //}
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
  double RatioTotalIntegrals_pb_over_pp_dphi = QuickIntegral_pb/(QuickIntegral_pp+QuickIntegral_bg);
  cout << "RatioTotalIntegrals_pb_over_pp_dphi: " << RatioTotalIntegrals_pb_over_pp_dphi << endl;

  QuickIntegral_pp = xjz_raw_stat_pp->IntegralAndError(-1,-1,QuickError,"width");
  cout << " integral rawpp: " << QuickIntegral_pp << " +/- " << QuickError << endl;
  QuickIntegral_pb = xjz_raw_stat_pb->IntegralAndError(-1,-1,QuickError,"width");
  TotalEqvStatsPbPb = QuickIntegral_pb*QuickIntegral_pb/(QuickError*QuickError);
  cout << " integral rawpb: " << QuickIntegral_pb << " +/- " << QuickError << endl;
  QuickIntegral_bg = xjz_bkg_stat_pb->IntegralAndError(-1,-1,QuickError,"width");
  cout << " integral bkgpb: " << QuickIntegral_bg << " +/- " << QuickError << endl;
  double RatioTotalIntegrals_pb_over_pp_xjz = QuickIntegral_pb/(QuickIntegral_pp+QuickIntegral_bg);
  cout << "RatioTotalIntegrals_pb_over_pp_xjz: " << RatioTotalIntegrals_pb_over_pp_xjz << endl;

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
  h_frame  ->GetXaxis()->SetTitle("X_{JZ}");
  h_frame  ->GetYaxis()->SetTitle("#frac{1}{N_{Z}} #frac{dN_{JZ}}{dx_{JZ}}");
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
  h_frame  ->GetXaxis()->SetTitle("#Delta#varphi");
  h_frame  ->GetYaxis()->SetTitle("#frac{1}{N_{Z}} #frac{dN_{JZ}}{d#Delta#varphi_{JK}}");
  h_frame->DrawCopy();
  dphi_raw_stat_pb ->DrawCopy("same");
  dphi_bkg_stat_pb ->DrawCopy("same");
  dphi_raw_stat_pp ->DrawCopy("same");


  dphi_sum_stat_pp = (TH1D*)dphi_raw_stat_pp->Clone("dphi_sum_stat_pp");
  dphi_sum_stat_pp ->Add(dphi_bkg_stat_pb);
  //dphi_sum_stat_pp->Scale(RatioTotalIntegrals_pb_over_pp_dphi);
  if(normalize_pp2pb)  dphi_sum_stat_pp ->Scale(dphi_raw_stat_pb->Integral(-1,-1)/dphi_sum_stat_pp->Integral(-1,-1));

  xjz_sum_stat_pp = (TH1D*)xjz_raw_stat_pp->Clone("xjz_sum_stat_pp");
  xjz_sum_stat_pp ->Add(xjz_bkg_stat_pb);
  //xjz_sum_stat_pp->Scale(RatioTotalIntegrals_pb_over_pp_xjz);
  if(normalize_pp2pb)  xjz_sum_stat_pp ->Scale(xjz_raw_stat_pb->Integral(-1,-1)/xjz_sum_stat_pp->Integral(-1,-1));


  for(int i=0; i<xjz_raw_stat_pb->GetNbinsX()+1; i++){
    xjz_raw_stat_pp->SetBinContent(i,xjz_raw_stat_pp->GetBinContent(i)*xjz_raw_stat_pp->GetXaxis()->GetBinWidth(i));
    xjz_raw_stat_pp->SetBinError  (i,xjz_raw_stat_pp->GetBinError  (i)*xjz_raw_stat_pp->GetXaxis()->GetBinWidth(i));
    xjz_sum_stat_pp->SetBinContent(i,xjz_sum_stat_pp->GetBinContent(i)*xjz_sum_stat_pp->GetXaxis()->GetBinWidth(i));
    xjz_sum_stat_pp->SetBinError  (i,xjz_sum_stat_pp->GetBinError  (i)*xjz_sum_stat_pp->GetXaxis()->GetBinWidth(i));
    xjz_raw_stat_pb->SetBinContent(i,xjz_raw_stat_pb->GetBinContent(i)*xjz_raw_stat_pb->GetXaxis()->GetBinWidth(i));
    xjz_raw_stat_pb->SetBinError  (i,xjz_raw_stat_pb->GetBinError  (i)*xjz_raw_stat_pb->GetXaxis()->GetBinWidth(i));
    xjz_bkg_stat_pb->SetBinContent(i,xjz_bkg_stat_pb->GetBinContent(i)*xjz_bkg_stat_pb->GetXaxis()->GetBinWidth(i));
    xjz_bkg_stat_pb->SetBinError  (i,xjz_bkg_stat_pb->GetBinError  (i)*xjz_bkg_stat_pb->GetXaxis()->GetBinWidth(i));

    float additional_correction = 1.0;
    if(normalize_pp_or_pbpb==0)
      additional_correction = xjz_sum_stat_pp->GetBinContent(i)/(xjz_sum_stat_pp ->GetBinError(i)*xjz_sum_stat_pp ->GetBinError(i));
    else if(normalize_pp_or_pbpb==1)
      additional_correction = xjz_raw_stat_pb->GetBinContent(i)/(xjz_raw_stat_pb ->GetBinError(i)*xjz_raw_stat_pb ->GetBinError(i));

    //if(additional_correction>0.0001 && additional_correction<1e8){
    if(xjz_raw_stat_pb->GetBinContent(i)){
      xjz_raw_stat_pp->SetBinContent(i,xjz_raw_stat_pp->GetBinContent(i)*additional_correction);
      xjz_raw_stat_pp->SetBinError  (i,xjz_raw_stat_pp->GetBinError  (i)*additional_correction);
      xjz_sum_stat_pp->SetBinContent(i,xjz_sum_stat_pp->GetBinContent(i)*additional_correction);
      xjz_sum_stat_pp->SetBinError  (i,xjz_sum_stat_pp->GetBinError  (i)*additional_correction);
      xjz_raw_stat_pb->SetBinContent(i,xjz_raw_stat_pb->GetBinContent(i)*additional_correction);
      xjz_raw_stat_pb->SetBinError  (i,xjz_raw_stat_pb->GetBinError  (i)*additional_correction);
      xjz_bkg_stat_pb->SetBinContent(i,xjz_bkg_stat_pb->GetBinContent(i)*additional_correction);
      xjz_bkg_stat_pb->SetBinError  (i,xjz_bkg_stat_pb->GetBinError  (i)*additional_correction);
    }
    else
      cout << "XJZ :: careful, the correction for bin " << i << " is messed up !" << endl;
  }
  for(int i=0; i<dphi_raw_stat_pb->GetNbinsX()+1; i++){
    dphi_raw_stat_pp->SetBinContent(i,dphi_raw_stat_pp->GetBinContent(i)*dphi_raw_stat_pp->GetXaxis()->GetBinWidth(i));
    dphi_raw_stat_pp->SetBinError  (i,dphi_raw_stat_pp->GetBinError  (i)*dphi_raw_stat_pp->GetXaxis()->GetBinWidth(i));
    dphi_sum_stat_pp->SetBinContent(i,dphi_sum_stat_pp->GetBinContent(i)*dphi_sum_stat_pp->GetXaxis()->GetBinWidth(i));
    dphi_sum_stat_pp->SetBinError  (i,dphi_sum_stat_pp->GetBinError  (i)*dphi_sum_stat_pp->GetXaxis()->GetBinWidth(i));
    dphi_raw_stat_pb->SetBinContent(i,dphi_raw_stat_pb->GetBinContent(i)*dphi_raw_stat_pb->GetXaxis()->GetBinWidth(i));
    dphi_raw_stat_pb->SetBinError  (i,dphi_raw_stat_pb->GetBinError  (i)*dphi_raw_stat_pb->GetXaxis()->GetBinWidth(i));
    dphi_bkg_stat_pb->SetBinContent(i,dphi_bkg_stat_pb->GetBinContent(i)*dphi_bkg_stat_pb->GetXaxis()->GetBinWidth(i));
    dphi_bkg_stat_pb->SetBinError  (i,dphi_bkg_stat_pb->GetBinError  (i)*dphi_bkg_stat_pb->GetXaxis()->GetBinWidth(i));

    float additional_correction = 1.0;
    if(normalize_pp_or_pbpb==0)
      additional_correction = dphi_sum_stat_pp->GetBinContent(i)/(dphi_sum_stat_pp ->GetBinError(i)*dphi_sum_stat_pp ->GetBinError(i));
    else if(normalize_pp_or_pbpb==1)
      additional_correction = dphi_raw_stat_pb->GetBinContent(i)/(dphi_raw_stat_pb ->GetBinError(i)*dphi_raw_stat_pb ->GetBinError(i));


    //if(additional_correction>0.0001 && additional_correction<1e8){
    if(dphi_raw_stat_pb->GetBinContent(i)){
      dphi_raw_stat_pp->SetBinContent(i,dphi_raw_stat_pp->GetBinContent(i)*additional_correction);
      dphi_raw_stat_pp->SetBinError  (i,dphi_raw_stat_pp->GetBinError  (i)*additional_correction);
      dphi_sum_stat_pp->SetBinContent(i,dphi_sum_stat_pp->GetBinContent(i)*additional_correction);
      dphi_sum_stat_pp->SetBinError  (i,dphi_sum_stat_pp->GetBinError  (i)*additional_correction);
      dphi_raw_stat_pb->SetBinContent(i,dphi_raw_stat_pb->GetBinContent(i)*additional_correction);
      dphi_raw_stat_pb->SetBinError  (i,dphi_raw_stat_pb->GetBinError  (i)*additional_correction);
      dphi_bkg_stat_pb->SetBinContent(i,dphi_bkg_stat_pb->GetBinContent(i)*additional_correction);
      dphi_bkg_stat_pb->SetBinError  (i,dphi_bkg_stat_pb->GetBinError  (i)*additional_correction);
    }
    else
      cout << "DPHI :: careful, the correction for bin " << i << " is messed up !" << endl;
  }

  c1->cd(2);
  h_frame  ->GetYaxis()->SetRangeUser(yaxis_xjz[0],yaxis_xjz[1]);
  h_frame  ->GetXaxis()->SetRangeUser(0,2.2);
  h_frame  ->GetXaxis()->SetTitle("X_{JZ}");
  if(normalize_pp_or_pbpb==0)      h_frame  ->GetYaxis()->SetTitle("effective counts (pp)");
  else if(normalize_pp_or_pbpb==1) h_frame  ->GetYaxis()->SetTitle("effective counts (PbPb)");
  h_frame->DrawCopy();
  xjz_raw_stat_pb  ->DrawCopy("same");
  xjz_bkg_stat_pb  ->DrawCopy("same");
  xjz_raw_stat_pp  ->DrawCopy("same");

  c1->cd(5);
  h_frame  ->GetYaxis()->SetRangeUser(yaxis_dphi[0],yaxis_dphi[1]);
  h_frame  ->GetXaxis()->SetRangeUser(0,3.2);
  h_frame  ->GetXaxis()->SetTitle("#Delta#varphi");
  if(normalize_pp_or_pbpb==0)      h_frame  ->GetYaxis()->SetTitle("effective counts (pp)");
  else if(normalize_pp_or_pbpb==1) h_frame  ->GetYaxis()->SetTitle("effective counts (PbPb)");
  h_frame->DrawCopy();
  dphi_raw_stat_pb ->DrawCopy("same");
  dphi_bkg_stat_pb ->DrawCopy("same");
  dphi_raw_stat_pp ->DrawCopy("same");

  c1->cd(3);
  h_frame  ->GetYaxis()->SetRangeUser(yaxis_xjz[0],yaxis_xjz[1]);
  h_frame  ->GetXaxis()->SetRangeUser(0,2.2);
  h_frame  ->GetXaxis()->SetTitle("X_{JZ}");
  if(normalize_pp_or_pbpb==0)      h_frame  ->GetYaxis()->SetTitle("effective counts (pp)");
  else if(normalize_pp_or_pbpb==1) h_frame  ->GetYaxis()->SetTitle("effective counts (PbPb)");
  h_frame->DrawCopy();
  xjz_sum_stat_pp  ->SetMarkerStyle(24);
  xjz_sum_stat_pp  ->SetMarkerColor(8);
  xjz_sum_stat_pp  ->SetLineColor  (8);
  xjz_raw_stat_pb  ->DrawCopy("same");
  xjz_sum_stat_pp  ->DrawCopy("same");

  TLegend *c1c = new TLegend(0.4,0.65,0.89,0.89);
  c1c->AddEntry(xjz_raw_stat_pb,"PbPb non-sub","P");
  c1c->AddEntry(xjz_sum_stat_pp,"pp + mix evts","P");
  c1c->Draw();

  c1->cd(6);
  h_frame  ->GetYaxis()->SetRangeUser(yaxis_dphi[0],yaxis_dphi[1]);
  h_frame  ->GetXaxis()->SetRangeUser(0,3.2);
  h_frame  ->GetXaxis()->SetTitle("#Delta#varphi");
  if(normalize_pp_or_pbpb==0)      h_frame  ->GetYaxis()->SetTitle("effective counts (pp)");
  else if(normalize_pp_or_pbpb==1) h_frame  ->GetYaxis()->SetTitle("effective counts (PbPb)");
  h_frame->DrawCopy();
  dphi_sum_stat_pp ->SetMarkerStyle(24);
  dphi_sum_stat_pp ->SetMarkerColor(8);
  dphi_sum_stat_pp ->SetLineColor  (8);
  dphi_raw_stat_pb ->DrawCopy("same");
  dphi_sum_stat_pp ->DrawCopy("same");



  cout << endl << endl;

  cout << "imax " << dphi_raw_stat_pb->GetNbinsX() << "  number of bins" << endl;
  cout << "jmax " << 1 << "  number of backgrounds" << endl;
  cout << "kmax " << 16 << "  number of nuisance parameters (sources of systematical uncertainties)" << endl;
  cout << "------------" << endl;


  cout << "bin           ";
  for(int i=1; i<dphi_raw_stat_pb->GetNbinsX()+1; i++)
    cout << i << " \t";
  cout << endl << "observation   ";
  for(int i=1; i<dphi_raw_stat_pb->GetNbinsX()+1; i++)
    cout << floor(0.5+dphi_raw_stat_pb->GetBinContent(i)) << " \t";
  cout << endl << "------------" << endl;

  cout << "bin               ";
  for(int i=1; i<dphi_raw_stat_pb->GetNbinsX()+1; i++)
    cout << i << " \t" << i << " \t";
  cout << endl << "process           ";
  for(int i=1; i<dphi_raw_stat_pb->GetNbinsX()+1; i++)
    cout << "pp " << " \t" << "bkg" << " \t";
  cout << endl << "process           ";
  for(int i=1; i<dphi_raw_stat_pb->GetNbinsX()+1; i++)
    cout << "0" << " \t" << "1" << " \t";
  cout << endl << "rate              ";
  for(int i=1; i<dphi_raw_stat_pb->GetNbinsX()+1; i++){
    sprintf(saythis,"%1.2f \t%1.2f",dphi_raw_stat_pp->GetBinContent(i),dphi_bkg_stat_pb->GetBinContent(i));
    cout << saythis << " \t";
  }
  cout << endl << "------------" << endl;

  for(int i=1; i<dphi_raw_stat_pp->GetNbinsX()+1; i++){
    sprintf(saythis,"statHIpp%d lnN\t",i);
    cout << saythis;
    for(int j=1; j<dphi_raw_stat_pp->GetNbinsX()+1; j++){
      if(i==j){
        sprintf(saythis,"%1.2f\t -  \t",dphi_raw_stat_pp->GetBinError(j)/dphi_raw_stat_pp->GetBinContent(j)+1);
        cout << saythis;
      }
      else
        cout << " - " << "\t" << " - " << " \t";
    }
    cout << endl;
  }
  for(int i=1; i<dphi_bkg_stat_pb->GetNbinsX()+1; i++){
    sprintf(saythis,"statHIbkg%d lnN\t",i);
    cout << saythis;
    for(int j=1; j<dphi_bkg_stat_pb->GetNbinsX()+1; j++){
      if(i==j){
        sprintf(saythis," - \t%1.2f \t",dphi_bkg_stat_pb->GetBinError(j)/dphi_bkg_stat_pb->GetBinContent(j)+1);
        cout << saythis;
      }
      else
        cout << " - " << "\t" << " - " << " \t";
    }
    cout << endl;
  }
  cout << "------------" << endl;







  TF1 *jpI = new TF1("jpI", "TMath::PoissonI(x,[0])",0,20);
  TF1 *jph = new TF1("jph", "TMath::Poisson(x,[0])",0,2000);
  TF1 *jpl = new TF1("jpl", "TMath::Poisson(x,[0])",0,100);
  TF1 *jg  = new TF1("jg", "gaus(0)",0,5);
  jg->SetParameters(1,1,0.05);

  h_Chisq_dphi = new TH1D("h_Chisq_dphi","#chi^{2};#chi^{2} of #Delta#varphi_{JZ}",100,0,50);
  h_KolSm_dphi = new TH1D("h_KolSm_dphi","Kolmogorov-Smirnov;KS of #Delta#varphi_{JZ}"      ,550,0,1.1);
  h_AndDr_dphi = new TH1D("h_AndDr_dphi","Anderson-Darling;AD of #Delta#varphi_{JZ}"      ,100,0,50);

  h_Chisq_xjz = new TH1D("h_Chisq_xjz","#chi^{2};#chi^{2} of x_{JZ}",100,0,50);
  h_KolSm_xjz = new TH1D("h_KolSm_xjz","Kolmogorov-Smirnov;KS of x_{JZ}"      ,550,0,1.1);
  //h_KolSm_xjz = new TH1D("h_KolSm_xjz","Kolmogorov-Smirnov;KS of x_{JZ}"      ,110,0,55);
  h_AndDr_xjz = new TH1D("h_AndDr_xjz","Anderson-Darling;AD of x_{JZ}"      ,100,0,50);

  TCanvas *c3 = new TCanvas("c3","c3");
  c3->Divide(5,4);

  TCanvas *c4 = new TCanvas("c4","c4");
  c4->Divide(5,4);


  //_________________________________________xjz toys_________________________________________________________
  TH1D *h_toy_xjz = (TH1D*)xjz_raw_stat_pb->Clone("h_toy_xjz");
  TH1D* h_toyEx_xjz[10];
  h_toy_xjz->SetMarkerStyle(24);
  h_toy_xjz->SetMarkerColor(14);
  h_toy_xjz->SetLineColor  (14);

  double qbinsx_xjz[17] = {0.0};
  for(int i=0; i<16; i++)
    qbinsx_xjz[i+1] = h_toy_xjz->GetXaxis()->GetBinUpEdge(i+1);
  //TH2D *h_TrueVals_xjz = new TH2D("h_TrueVals_xjz","h_TrueVals_xjz",16, qbinsx_xjz, 1000,0,100);
  TH2D *h_TrueVals_xjz = new TH2D("h_TrueVals_xjz","h_TrueVals_xjz",100,0,2.0, 1000,0,100);

  for(int i=0; i<nPseudoExp; i++){

    h_toy_xjz->Reset();
    double ThisTrueVal = 0.0;

    //for(int j=1; j<h_toy_xjz->GetNbinsX()+1; j++){
    for(int j=FirstBin_xjz; j<LastBin_xjz+1; j++){
      if(include_ppflucts_in_toy==0)
        ThisTrueVal  = xjz_sum_stat_pp->GetBinContent(j);
      else if(include_ppflucts_in_toy==1){
        jpl ->SetParameter(0,xjz_sum_stat_pp->GetBinContent(j)*xjz_sum_stat_pp->GetBinContent(j)/(xjz_sum_stat_pp ->GetBinError(j)*xjz_sum_stat_pp ->GetBinError(j)));
        jph ->SetParameter(0,xjz_sum_stat_pp->GetBinContent(j)*xjz_sum_stat_pp->GetBinContent(j)/(xjz_sum_stat_pp ->GetBinError(j)*xjz_sum_stat_pp ->GetBinError(j)));
        if(xjz_sum_stat_pp->GetBinContent(j)*xjz_sum_stat_pp->GetBinContent(j)/(xjz_sum_stat_pp ->GetBinError(j)*xjz_sum_stat_pp ->GetBinError(j)) < 60)
          ThisTrueVal  = jpl->GetRandom()*xjz_sum_stat_pp ->GetBinError(j)*xjz_sum_stat_pp ->GetBinError(j)/xjz_sum_stat_pp->GetBinContent(j);
        else
          ThisTrueVal  = jph->GetRandom()*xjz_sum_stat_pp ->GetBinError(j)*xjz_sum_stat_pp ->GetBinError(j)/xjz_sum_stat_pp->GetBinContent(j);

        //jg->SetParameter(2,xjz_sum_stat_pp->GetBinError(j)/xjz_sum_stat_pp->GetBinContent(j));
        ////jg->SetParameter(2, (1.0-TMath::Sqrt( xjz_raw_stat_pb->GetBinError(j)*xjz_raw_stat_pb->GetBinError(j) + xjz_sum_stat_pp->GetBinError(j)*xjz_sum_stat_pp->GetBinError(j) ) / xjz_raw_stat_pb->GetBinError(j) )*xjz_sum_stat_pp->GetBinError(j)/xjz_sum_stat_pp->GetBinContent(j) );
        ////ThisTrueVal = xjz_sum_stat_pp->GetBinContent(j)*jg->GetRandom();
      }
      h_TrueVals_xjz->Fill(xjz_sum_stat_pp->GetBinCenter(j)+1.0*h_TrueVals_xjz->GetXaxis()->GetBinWidth(3),ThisTrueVal);

      if(xjz_raw_stat_pb ->GetBinError(j)==0 || xjz_raw_stat_pb ->GetBinError(j)==0)
        cout << " #####  !! We have a zero-counts bin (" << j << ") !!  #####" << endl;

      double ThisPoisStat = xjz_raw_stat_pb->GetBinContent(j)*xjz_raw_stat_pb->GetBinContent(j) / 
                            (xjz_raw_stat_pb ->GetBinError(j)*xjz_raw_stat_pb ->GetBinError(j));
      jpl ->SetParameter(0,ThisPoisStat);
      jph ->SetParameter(0,ThisPoisStat);
      jpI->SetParameter(0,0.25);
      double MyValue = 0.0;
      if(xjz_raw_stat_pb->GetBinContent(j)>0){
        if(ThisPoisStat<60)
          MyValue = jpl->GetRandom();
        else
          MyValue = jph->GetRandom();
        h_toy_xjz->SetBinContent(j,MyValue);
        h_toy_xjz->SetBinError  (j,sqrt(MyValue));
      }
      else{
        MyValue = floor(jpI->GetRandom());
        h_toy_xjz->SetBinContent(j,MyValue);
        h_toy_xjz->SetBinError  (j,sqrt(MyValue));
      }
      h_toy_xjz->SetBinContent(j,MyValue*ThisTrueVal/ThisPoisStat);
      h_toy_xjz->SetBinError  (j,sqrt(MyValue)*ThisTrueVal/ThisPoisStat);
    }
    //jp ->SetParameter(0,TotalEqvStatsPbPb);
    //h_toy_xjz->Scale(RatioTotalIntegrals_pb_over_pp*jp->GetRandom()/TotalEqvStatsPbPb);
    if(renormalize_toy==1)  h_toy_xjz->Scale(RatioTotalIntegrals_pb_over_pp_xjz);
    //h_toy_xjz->Scale(10);

    double chi2;
    int ndf,igood;
    //double Pvalue = xjz_sum_stat_pp->Chi2TestX(h_toy_xjz,chi2,ndf,igood,"WW");
    double Pvalue = 1.0;
    chi2 = CalcChiSq(xjz_sum_stat_pp,h_toy_xjz, FirstBin_xjz,LastBin_xjz);
    double kolo = xjz_sum_stat_pp->KolmogorovTest(h_toy_xjz,"M");//"MN");
    kolo = CalcKoloSm(xjz_sum_stat_pp, h_toy_xjz, FirstBin_xjz,LastBin_xjz, 1, -1);
    //double andr = xjz_sum_stat_pp->AndersonDarlingTest(h_toy_xjz,"T");
    double andr = 1;

    h_Chisq_xjz->Fill(chi2);
    h_KolSm_xjz->Fill(kolo);
    h_AndDr_xjz->Fill(andr);

    if(i<20) c3->cd(i+1);
    else     c3->cd(20);
    xjz_raw_stat_pb->SetMarkerColor(19);
    xjz_raw_stat_pb->SetLineColor  (19);
    xjz_raw_stat_pb->DrawCopy();
    xjz_raw_stat_pb->SetMarkerColor(1);
    xjz_raw_stat_pb->SetLineColor  (1);
    xjz_sum_stat_pp->DrawCopy("same");
    h_toy_xjz->DrawCopy("same");

    if(i<10){
      cout << "XJZ:   " << "pval: " << Pvalue << "   chi2: " << chi2 << "   kolo: " << kolo << "   andr: " << andr << endl;
      sprintf(saythis,"h_toyEx_xjz_%d",i);
      h_toyEx_xjz[i] = (TH1D*)h_toy_xjz->Clone(saythis);
    }

  }

  double chi2_xjz;
  int ndf,igood;
  //double Pvalue = xjz_raw_stat_pb->Chi2TestX(xjz_sum_stat_pp,chi2_xjz,ndf,igood,"WW");
  double Pvalue = 1.0;
  chi2_xjz = CalcChiSq(xjz_raw_stat_pb,xjz_sum_stat_pp, FirstBin_xjz,LastBin_xjz);
  //h_pb[i]->Chi2Test(h_pp[i],"WW,P");
  double kolo_xjz = xjz_raw_stat_pb->KolmogorovTest(xjz_sum_stat_pp,"M");//"MN");
  kolo_xjz = CalcKoloSm(xjz_raw_stat_pb, xjz_sum_stat_pp, FirstBin_xjz,LastBin_xjz, 1, 0);
  //double andr = xjz_raw_stat_pb->AndersonDarlingTest(xjz_sum_stat_pp,"T");
  double andr_xjz = 1;

  cout << endl << "Real Values! " << endl << "chi2: " << chi2_xjz << "   kolo: " << kolo_xjz << "   andr: " << andr_xjz << endl;
  cout << endl << "ch " << 1.0-h_Chisq_xjz->Integral(1,h_Chisq_xjz->GetXaxis()->FindBin(chi2_xjz))/h_Chisq_xjz->Integral(1,h_Chisq_xjz->GetNbinsX()) << "   "
               << "ks " << 1.0-h_KolSm_xjz->Integral(1,h_KolSm_xjz->GetXaxis()->FindBin(kolo_xjz))/h_KolSm_xjz->Integral(1,h_KolSm_xjz->GetNbinsX()) << "   "
               << "ad " << 1.0-h_AndDr_xjz->Integral(1,h_AndDr_xjz->GetXaxis()->FindBin(andr_xjz))/h_AndDr_xjz->Integral(1,h_AndDr_xjz->GetNbinsX()) << endl;

  cout << "compared to built in P-Values: " << endl;
  cout << "ch " << xjz_raw_stat_pb->Chi2Test(xjz_sum_stat_pp,"WW") << "   "
       << "ks " << xjz_raw_stat_pb->KolmogorovTest(xjz_sum_stat_pp) << "   "
       << "ad " << xjz_raw_stat_pb->AndersonDarlingTest(xjz_sum_stat_pp) << endl;
  h_xjz_chi2_KS_AD  ->SetBinContent(1,chi2_xjz);
  h_xjz_chi2_KS_AD  ->SetBinContent(2,kolo_xjz);
  h_xjz_chi2_KS_AD  ->SetBinContent(3,andr_xjz);
  //_________________________________________xjz toys_________________________________________________________



  //________________________________________dphi toys_________________________________________________________
  TH1D *h_toy_dphi = (TH1D*)dphi_raw_stat_pb->Clone("h_toy_dphi");
  TH1D* h_toyEx_dphi[10];
  h_toy_dphi->SetMarkerStyle(24);
  h_toy_dphi->SetMarkerColor(14);
  h_toy_dphi->SetLineColor  (14);

  double qbinsx_dphi[9] = {0.0};
  for(int i=0; i<8; i++)
    qbinsx_dphi[i+1] = h_toy_dphi->GetXaxis()->GetBinUpEdge(i+1);
  //TH2D *h_TrueVals_dphi = new TH2D("h_TrueVals_dphi","h_TrueVals_dphi",8, qbinsx_dphi, 1000,0,100);
  TH2D *h_TrueVals_dphi = new TH2D("h_TrueVals_dphi","h_TrueVals_dphi",100,0,3.14159, 1000,0,100);

  for(int i=0; i<nPseudoExp; i++){

    h_toy_dphi->Reset();
    double ThisTrueVal = 0.0;

    for(int j=1; j<h_toy_dphi->GetNbinsX()+1; j++){
      double randompoisson = 0.0;
      if(include_ppflucts_in_toy==0)
        ThisTrueVal  = dphi_sum_stat_pp->GetBinContent(j);
      else if(include_ppflucts_in_toy==1){
        jpl ->SetParameter(0,dphi_sum_stat_pp->GetBinContent(j)*dphi_sum_stat_pp->GetBinContent(j)/(dphi_sum_stat_pp ->GetBinError(j)*dphi_sum_stat_pp ->GetBinError(j)));
        jph ->SetParameter(0,dphi_sum_stat_pp->GetBinContent(j)*dphi_sum_stat_pp->GetBinContent(j)/(dphi_sum_stat_pp ->GetBinError(j)*dphi_sum_stat_pp ->GetBinError(j)));
        if(dphi_sum_stat_pp->GetBinContent(j)*dphi_sum_stat_pp->GetBinContent(j)/(dphi_sum_stat_pp ->GetBinError(j)*dphi_sum_stat_pp ->GetBinError(j)) < 60)
          randompoisson = jpl->GetRandom();
        else
          randompoisson = jph->GetRandom();
        ThisTrueVal  = randompoisson*dphi_sum_stat_pp ->GetBinError(j)*dphi_sum_stat_pp ->GetBinError(j)/dphi_sum_stat_pp->GetBinContent(j);
        //jp ->SetParameter(0,dphi_sum_stat_pp->GetBinContent(j));
        //ThisTrueVal  = jp->GetRandom();

        //jg->SetParameter(2,dphi_sum_stat_pp->GetBinError(j)/dphi_sum_stat_pp->GetBinContent(j));
        ////jg->SetParameter(2, (1.0-TMath::Sqrt( dphi_raw_stat_pb->GetBinError(j)*dphi_raw_stat_pb->GetBinError(j) + dphi_sum_stat_pp->GetBinError(j)*dphi_sum_stat_pp->GetBinError(j) ) / dphi_raw_stat_pb->GetBinError(j) )*dphi_sum_stat_pp->GetBinError(j)/dphi_sum_stat_pp->GetBinContent(j) );
        ////ThisTrueVal = dphi_sum_stat_pp->GetBinContent(j)*jg->GetRandom();
      }
      //cout << "dphi true value: " << ThisTrueVal << endl;
      h_TrueVals_dphi->Fill(dphi_sum_stat_pp->GetBinCenter(j)+1.0*h_TrueVals_dphi->GetXaxis()->GetBinWidth(3),ThisTrueVal);

      double ThisPoisStat = dphi_raw_stat_pb->GetBinContent(j)*dphi_raw_stat_pb->GetBinContent(j) / 
                            (dphi_raw_stat_pb ->GetBinError(j)*dphi_raw_stat_pb ->GetBinError(j));
      jpl ->SetParameter(0,ThisPoisStat);
      jph ->SetParameter(0,ThisPoisStat);
      jpI->SetParameter(0,0.25);
      double MyValue = 0.0;
      if(dphi_raw_stat_pb->GetBinContent(j)>0){
        if(ThisPoisStat<60)
          MyValue = jpl->GetRandom();
        else
          MyValue = jph->GetRandom();
        h_toy_dphi->SetBinContent(j,MyValue);
        h_toy_dphi->SetBinError  (j,sqrt(MyValue));
      }
      else{
        MyValue = floor(jpI->GetRandom());
        h_toy_dphi->SetBinContent(j,MyValue);
        h_toy_dphi->SetBinError  (j,sqrt(MyValue));
      }
      h_toy_dphi->SetBinContent(j,MyValue*ThisTrueVal/ThisPoisStat);
      h_toy_dphi->SetBinError  (j,sqrt(MyValue)*ThisTrueVal/ThisPoisStat);
    }
    //jp ->SetParameter(0,TotalEqvStatsPbPb);
    //h_toy_dphi->Scale(RatioTotalIntegrals_pb_over_pp*jp->GetRandom()/TotalEqvStatsPbPb);
    if(renormalize_toy==1)  h_toy_dphi->Scale(RatioTotalIntegrals_pb_over_pp_dphi);
    //h_toy_dphi->Scale(10);

    double chi2;
    int ndf,igood;
    double Pvalue = 1.0;
    //double Pvalue = dphi_sum_stat_pp->Chi2TestX(h_toy_dphi,chi2,ndf,igood,"WW");
    chi2 = CalcChiSq(dphi_sum_stat_pp,h_toy_dphi, -1,-1);
    //h_pb[i]->Chi2Test(h_pp[i],"WW,P");
    double kolo = dphi_sum_stat_pp->KolmogorovTest(h_toy_dphi,"M");
    kolo = CalcKoloSm(dphi_sum_stat_pp, h_toy_dphi, -1,-1, 1, -1);
    double andr = dphi_sum_stat_pp->AndersonDarlingTest(h_toy_dphi,"T");

    h_Chisq_dphi->Fill(chi2);
    h_KolSm_dphi->Fill(kolo);
    h_AndDr_dphi->Fill(andr);

    if(i<20) c4->cd(i+1);
    else     c4->cd(20);
    dphi_raw_stat_pb->SetMarkerColor(19);
    dphi_raw_stat_pb->SetLineColor  (19);
    dphi_raw_stat_pb->DrawCopy();
    dphi_raw_stat_pb->SetMarkerColor(1);
    dphi_raw_stat_pb->SetLineColor  (1);
    dphi_sum_stat_pp->DrawCopy("same");
    h_toy_dphi->DrawCopy("same");

    if(i<10) {
      cout << "DPHI:  " << "pval: " << Pvalue << "   chi2: " << chi2 << "   kolo: " << kolo << "   andr: " << andr << endl;
      sprintf(saythis,"h_toyEx_dphi_%d",i);
      h_toyEx_dphi[i] = (TH1D*)h_toy_dphi->Clone(saythis);
    }

  }

  double chi2_dphi;
  //Pvalue = dphi_raw_stat_pb->Chi2TestX(dphi_sum_stat_pp,chi2_dphi,ndf,igood,"WW");
  Pvalue = 1.0;
  chi2_dphi = CalcChiSq(dphi_raw_stat_pb,dphi_sum_stat_pp, -1,-1);
  //h_pb[i]->Chi2Test(h_pp[i],"WW,P");
  double kolo_dphi = dphi_raw_stat_pb->KolmogorovTest(dphi_sum_stat_pp,"M");
  kolo_dphi = CalcKoloSm(dphi_raw_stat_pb, dphi_sum_stat_pp, -1,-1, 1, 1);
  double andr_dphi = dphi_raw_stat_pb->AndersonDarlingTest(dphi_sum_stat_pp,"T");

  cout << endl << "Real Values! " << endl << "chi2: " << chi2_dphi << "   kolo: " << kolo_dphi << "   andr: " << andr_dphi << endl;
  cout << endl << "ch " << 1.0-h_Chisq_dphi->Integral(1,h_Chisq_dphi->GetXaxis()->FindBin(chi2_dphi))/h_Chisq_dphi->Integral(1,h_Chisq_dphi->GetNbinsX()) << "   "
               << "ks " << 1.0-h_KolSm_dphi->Integral(1,h_KolSm_dphi->GetXaxis()->FindBin(kolo_dphi))/h_KolSm_dphi->Integral(1,h_KolSm_dphi->GetNbinsX()) << "   "
               << "ad " << 1.0-h_AndDr_dphi->Integral(1,h_AndDr_dphi->GetXaxis()->FindBin(andr_dphi))/h_AndDr_dphi->Integral(1,h_AndDr_dphi->GetNbinsX()) << endl;

  cout << "compared to built in P-Values: " << endl;
  cout << "ch " << dphi_raw_stat_pb->Chi2Test(dphi_sum_stat_pp,"WW") << "   "
       << "ks " << dphi_raw_stat_pb->KolmogorovTest(dphi_sum_stat_pp) << "   "
       << "ad " << dphi_raw_stat_pb->AndersonDarlingTest(dphi_sum_stat_pp) << endl;
  h_dphi_chi2_KS_AD  ->SetBinContent(1,chi2_dphi);
  h_dphi_chi2_KS_AD  ->SetBinContent(2,kolo_dphi);
  h_dphi_chi2_KS_AD  ->SetBinContent(3,andr_dphi);
  //________________________________________dphi toys_________________________________________________________


  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->Divide(3,2);
  c2->cd(1);
  h_Chisq_xjz->Draw();
  TLine *line_chi2_xjz = new TLine(chi2_xjz,0,chi2_xjz,h_Chisq_xjz->GetMaximum());
  line_chi2_xjz->SetLineColor(2);
  line_chi2_xjz->SetLineStyle(2);
  line_chi2_xjz->Draw("same");

  c2->cd(2);
  //h_KolSm_xjz->GetXaxis()->SetRangeUser(0,0.2);
  h_KolSm_xjz->GetXaxis()->SetRangeUser(0,50.0);
  h_KolSm_xjz->Draw();
  TLine *line_kolo_xjz = new TLine(kolo_xjz,0,kolo_xjz,h_KolSm_xjz->GetMaximum());
  line_kolo_xjz->SetLineColor(2);
  line_kolo_xjz->SetLineStyle(2);
  line_kolo_xjz->Draw("same");

  c2->cd(3);
  h_AndDr_xjz->Draw();

  c2->cd(4);
  h_Chisq_dphi->Draw();
  TLine *line_chi2_dphi = new TLine(chi2_dphi,0,chi2_dphi,h_Chisq_dphi->GetMaximum());
  line_chi2_dphi->SetLineColor(2);
  line_chi2_dphi->SetLineStyle(2);
  line_chi2_dphi->Draw("same");

  c2->cd(5);
  h_KolSm_dphi->GetXaxis()->SetRangeUser(0,0.2);
  h_KolSm_dphi->Draw();
  TLine *line_kolo_dphi = new TLine(kolo_dphi,0,kolo_dphi,h_KolSm_dphi->GetMaximum());
  line_kolo_dphi->SetLineColor(2);
  line_kolo_dphi->SetLineStyle(2);
  line_kolo_dphi->Draw("same");

  c2->cd(6);
  h_AndDr_dphi->Draw();

  TCanvas *c5 = new TCanvas("c5","c5",400,800);
  c5->Divide(1,2);
  c5->cd(1);
  h_KoloHist[0]->SetTitle("KS distribution X_{JZ}");
  h_KoloHist[0]->Draw();
  c5->cd(2);
  h_KoloHist[1]->SetTitle("KS distribution #Delta#varphi");
  h_KoloHist[1]->Draw();

  TCanvas *c6 = new TCanvas("c6","c6",800,400);
  c6->Divide(2,1);
  h_TrueVals_dphi->SetMarkerSize (0.2);
  h_TrueVals_xjz ->SetMarkerSize (0.2);
  h_TrueVals_dphi->SetMarkerStyle(20);
  h_TrueVals_xjz ->SetMarkerStyle(20);
  h_TrueVals_dphi->SetMarkerColor(46);
  h_TrueVals_xjz ->SetMarkerColor(46);
  dphi_sum_stat_pp->SetLineWidth(2);
  xjz_sum_stat_pp ->SetLineWidth(2);
  c6->cd(1);
  h_TrueVals_dphi ->DrawCopy("");
  dphi_sum_stat_pp->DrawCopy("same,pe");
  c6->cd(2);
  h_TrueVals_xjz  ->DrawCopy("");
  xjz_sum_stat_pp ->DrawCopy("same,pe");

  if(writeoutputfile==1){

    //sprintf(saythis,"JasonToyMC_norm%d_ppflucts%d_nExp%d_Error%2.2f.root",normalize_pp_or_pbpb,include_ppflucts_in_toy,nPseudoExp,BlowUpErrors);
    //sprintf(saythis,"Fall2016/JasonToyMC_31Aug_norm%d_ppflucts%d_nExp%d_Error%2.2f_RenormalizeToy%d_normpp2pb%d_xjzcut_%d_%d.root",normalize_pp_or_pbpb,include_ppflucts_in_toy,nPseudoExp,BlowUpErrors,renormalize_toy,normalize_pp2pb,FirstBin_xjz,LastBin_xjz);
    //sprintf(saythis,"Jan2017/JasonToyMC_31JanFrankenstein_norm%d_ppflucts%d_nExp%d_Error%2.2f_RenormalizeToy%d_normpp2pb%d_xjzcut_%d_%d.root",normalize_pp_or_pbpb,include_ppflucts_in_toy,nPseudoExp,BlowUpErrors,renormalize_toy,normalize_pp2pb,FirstBin_xjz,LastBin_xjz);
    sprintf(saythis,"../Jan2017/JasonToyMC_31Jan_norm%d_ppflucts%d_nExp%d_Error%2.2f_RenormalizeToy%d_normpp2pb%d_xjzcut_%d_%d.root",normalize_pp_or_pbpb,include_ppflucts_in_toy,nPseudoExp,BlowUpErrors,renormalize_toy,normalize_pp2pb,FirstBin_xjz,LastBin_xjz);
    TFile *fout = new TFile(saythis,"NEW");
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    c6->Write();

    xjz_raw_stat_pp ->SetName("xjz_raw_stat_pp");
    xjz_raw_stat_pb ->SetName("xjz_raw_stat_pb");
    xjz_bkg_stat_pb ->SetName("xjz_bkg_stat_pb");
                                               
    dphi_raw_stat_pp->SetName("dphi_raw_stat_pp");
    dphi_raw_stat_pb->SetName("dphi_raw_stat_pb");
    dphi_bkg_stat_pb->SetName("dphi_bkg_stat_pb");

    xjz_raw_stat_pp ->Write();
    xjz_raw_stat_pb ->Write();
    xjz_bkg_stat_pb ->Write();

    dphi_raw_stat_pp->Write();
    dphi_raw_stat_pb->Write();
    dphi_bkg_stat_pb->Write();

    dphi_sum_stat_pp->Write();
    xjz_sum_stat_pp ->Write();

    //xjz_raw_syst_pp ->Write();
    //xjz_raw_syst_pb ->Write();
    //xjz_bkg_syst_pb ->Write();

    //dphi_raw_syst_pp->Write();
    //dphi_raw_syst_pb->Write();
    //dphi_bkg_syst_pb->Write();

    h_TrueVals_xjz ->Write();
    h_TrueVals_dphi->Write();

    for(int i=0; i<10; i++)
      h_toyEx_dphi[i] ->Write();
    for(int i=0; i<10; i++)
      h_toyEx_xjz[i]  ->Write();

    h_KoloHist[0]   ->Write();
    h_KoloHist[1]   ->Write();

    h_Chisq_dphi->Write();
    h_KolSm_dphi->Write();
    h_AndDr_dphi->Write();

    h_Chisq_xjz ->Write();
    h_KolSm_xjz ->Write();
    h_AndDr_xjz ->Write();

    h_xjz_chi2_KS_AD  ->Write();
    h_dphi_chi2_KS_AD ->Write();

    cout << endl << "Wrote " << saythis << endl;
    cout << "___________________________________________________________________________" << endl << endl;

  }

}

