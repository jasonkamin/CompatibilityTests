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

// 20 July 2016
// this should make the fully combined datacard !
void MakeCombineDatacard()
{

  double BlowUpErrors = 1.0;
  char saythis[500];

  //TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_5Jul2016.root");
  //TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSumPP_5Jul2016.root");
  //TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_HI_PP_28Jul2016.root");
  //TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_HI_PP_28Jul2016.root");

  //TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_HI_PP_31Aug2016.root");
  //TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/NEW/zJetHistogramSum_HI_PP_31Aug2016.root");
  //TFile *f_pb = TFile::Open("$OUTPUT/ZJet_Kaya/Jan23_2017/zJetHistogramSum_HI_PP_20170118.root");
  //TFile *f_pp = TFile::Open("$OUTPUT/ZJet_Kaya/Jan23_2017/zJetHistogramSum_HI_PP_20170118.root");
  TFile *f_pb = TFile::Open("$CODE/ZJet/Frankenstein.root");
  TFile *f_pp = TFile::Open("$CODE/ZJet/Frankenstein.root");

  //xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_ptBin0_hiBin1_jetRAW_final_norm");
  //xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_ptBin0_hiBin1_jetBKG_final_norm");
  //xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get("h1D_xjz_ptBin0_hiBin0_jetRAW_final_norm");
  //xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_binJER_ptBin0_hiBin1_jetRAW_final_norm");
  //xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_xjz_binJER_ptBin0_hiBin1_jetBKG_final_norm");
  //xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get("h1D_xjz_binJER_ptBin0_hiBin0_jetRAW_final_norm");
  xjz_raw_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("xjz_raw_stat_pb");
  xjz_bkg_stat_pb  = (TH1D*)f_pb->GetDirectory("HI")->Get("xjz_bkg_stat_pb");
  xjz_raw_stat_pp  = (TH1D*)f_pp->GetDirectory("PP")->Get("xjz_raw_stat_pp");

  dphi_raw_stat_pb = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_dphi_rebin_ptBin0_hiBin1_jetRAW_final_norm");
  dphi_bkg_stat_pb = (TH1D*)f_pb->GetDirectory("HI")->Get("h1D_dphi_rebin_ptBin0_hiBin1_jetBKG_final_norm");
  dphi_raw_stat_pp = (TH1D*)f_pp->GetDirectory("PP")->Get("h1D_dphi_rebin_ptBin0_hiBin0_jetRAW_final_norm");

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
  MakeDataCardFromHistos(dphi_raw_stat_pp, dphi_raw_stat_pb, dphi_bkg_stat_pb, -1,-1);
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
    sprintf(saythis,"%8d",int(floor( h_pp->GetBinContent(i)*h_pp->GetBinContent(i)/(h_pp->GetBinError(i)*h_pp->GetBinError(i))+0.5 )));
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
    sprintf(saythis,"%9.4f", h_pp->GetBinContent(i)/(h_pp->GetBinError(i)*h_pp->GetBinError(i)) );
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
    if(i<10)
      sprintf(saythis,"dx%d           param  %4.2f  100  [0,100]",i,h_pb->GetBinContent(i)/h_pp->GetBinContent(i));
    else
      sprintf(saythis,"dx%d          param  %4.2f  100  [0,100]",i,h_pb->GetBinContent(i)/h_pp->GetBinContent(i));
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
    sprintf(saythis,"statbkg%d lnN    ",j);
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
    cout << " --algo fixed --snapshotName MultiDimFit --forceRecreateNLL  -t 10" << endl << endl;

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
