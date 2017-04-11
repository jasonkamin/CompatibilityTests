TH1D *hist1;
TH1D *hist2;
TH1D* h_KoloHist[2];
TH1D *h_chi2;

#include "CompatibilityHelperFunctions.C"

void NormalizationTest(int ptBinKaya = 0)
{
  double rjz_pp   = 0.0;
  double rjz_pp_e = 0.0;
  double rjz_hi   = 0.0;
  double rjz_hi_e = 0.0;

  int verbosity = 0;

  // KS shape stat only
  //double KS_SHAPE_ARRAY[7] = {
  //  0.023,
  //  0.001,
  //  0.001,
  //  0.411,
  //  0.007,
  //  0.024,
  //  0.527
  //};

  // KS shape stat+syst, baysian toys
  double KS_SHAPE_ARRAY[7] = {
    0.35,
    0.08,
    0.07,
    0.82,
    0.07,
    0.14,
    0.85
  };

  // KS shape stat+syst, frequenist toys
  // this one is not the one to use ! 
  //double KS_SHAPE_ARRAY[7] = {
  //  0.63,
  //  0.37,
  //  0.42,
  //  0.8 ,
  //  0.14,
  //  0.3 ,
  //  0.87
  //};

  // ZpT > 60
  if(ptBinKaya==0){
    rjz_pp   = 0.628201;
    rjz_pp_e = 0.030439/rjz_pp;
    rjz_hi   = 0.530373;
    rjz_hi_e = 0.055281/rjz_hi;
  }

  // ZpT > 40
  if(ptBinKaya==1){
    rjz_pp   = 0.489709;
    rjz_pp_e = 0.0181964/rjz_pp;
    rjz_hi   = 0.366976;
    rjz_hi_e = 0.0310400/rjz_hi;
  }

  // ZpT > 50
  if(ptBinKaya==2){
    rjz_pp   = 2.85481;
    rjz_pp_e = 0.12081/rjz_pp;
    rjz_hi   = 2.36623;
    rjz_hi_e = 0.210768/rjz_hi;
  }

  // 40 < ZpT < 50
  if(ptBinKaya==3){
    rjz_pp   = 1.57844;
    rjz_pp_e = 0.125511/rjz_pp;
    rjz_hi   = 0.897416;
    rjz_hi_e = 0.216023/rjz_hi;
  }

  // 50 < ZpT < 60
  if(ptBinKaya==4){
    rjz_pp   = 2.25133;
    rjz_pp_e = 0.193707/rjz_pp;
    rjz_hi   = 1.82659;
    rjz_hi_e = 0.344131/rjz_hi;
  }

  // 60 < ZpT < 80
  if(ptBinKaya==5){
    rjz_pp   = 2.53911;
    rjz_pp_e = 0.188315/rjz_pp;
    rjz_hi   = 2.0002;
    rjz_hi_e = 0.362397/rjz_hi;
  }

  // 80 < ZpT
  if(ptBinKaya==6){
    rjz_pp   = 3.77375;
    rjz_pp_e = 0.242827/rjz_pp;
    rjz_hi   = 3.17039;
    rjz_hi_e = 0.404576/rjz_hi;
  }


  // in %
  double rjz_pp_syst[7] = {
    0.021,
    0.026,
    0.022,
    0.036,
    0.027,
    0.026,
    0.021
  };
  double rjz_hi_syst[7] = {
    0.055,
    0.081,
    0.058,
    0.156,
    0.155,
    0.089,
    0.083
  };


//##### rjz SYSTEMATICS #####
//
//Bins are zpt>40, zpt>50, zpt>60
//Values are in percentage
//rjz in PP, SIG correlation :
//2.5  2.1  2.0  jet energy scale
//0.4  0.5  0.3  electron energy scale
//0.3  0.3  0.2  jet energy resolution
//0.2  0.2  0.2  jet angular resolution
//2.6  2.2  2.1  total
//
//rjz in PbPb, SIG correlation :
//4.8  4.0  3.0   jet energy scale
//0.9  0.7  0.6   electron energy scale
//4.0  3.0  2.6   jet energy resolution
//0.9  0.6  0.6   Z energy resolution
//5.0  2.9  3.7   Z reconstruction efficiency correction
//8.1  5.8  5.5   total

//I have the rjz systematics only for paper plots, the bins are 40-50, 50-60, 60-80, 80+.
//For rjz systematics, all variations are done for each bin separately (yes, the analysis is redone for each bin).
//I do not have the values for pt>40, pt>50, pt>60 in hand. I can try to calculate them, but I am not sure they will be stable as the other bins.
//
//##### rjz SYSTEMATICS #####
//Values are in percentage
//xjz in PP, SIG correlation :
//3.3  2.5  2.4  1.9   jet energy scale
//1.2  0.9  0.8  0.6   electron energy scale
//0.5  0.4  0.5  0.4   jet energy resolution
//0.3  0.2  0.2  0.1   jet angular resolution
//3.6  2.7  2.6  2.1   total  hline
//
//xjz in PbPb, SIG correlation :
//11.6  10.0  4.6  3.1   jet energy scale
//2.3   3.8   3.0  3.5   electron energy scale
//7.7   7.3   3.7  2.9   jet energy resolution
//0.8   1.1   0.8  0.8   Z energy resolution
//6.5   8.5   5.8  6.1   Z reconstruction efficiency correction
//15.6  15.5  8.9  8.3   total



  // ZpT > 60, w/ syst
  //rjz_pp_e = TMath::Sqrt(rjz_pp_e*rjz_pp_e + 0.021*0.021*rjz_pp*rjz_pp);
  //rjz_hi_e = TMath::Sqrt(rjz_hi_e*rjz_hi_e + 0.055*0.055*rjz_hi*rjz_hi);



  cout << "_____________pp____________" << endl;
  cout << " error stat(%) = " << rjz_pp_e << "   syst(%)" << rjz_pp_syst[ptBinKaya] << "   syst(%)/sqrt(2) = " << rjz_pp_syst[ptBinKaya]/TMath::Sqrt(2) << endl;
  rjz_pp_e = TMath::Sqrt(rjz_pp_e*rjz_pp_e + rjz_pp_syst[ptBinKaya]*rjz_pp_syst[ptBinKaya]);
  cout << " error total(%) = " << rjz_pp_e << endl;

  cout << "____________PbPb___________" << endl;
  cout << " error stat(%) = " << rjz_hi_e << "   syst(%)" << rjz_hi_syst[ptBinKaya] << endl;
  rjz_hi_e = TMath::Sqrt(rjz_hi_e*rjz_hi_e + rjz_hi_syst[ptBinKaya]*rjz_hi_syst[ptBinKaya]);
  cout << " error total(%) = " << rjz_hi_e << endl;
  cout << endl;




  //double prob_shape = 0.025;// ZpT > 60
  //double prob_shape = 0.001;// ZpT > 40
  double prob_shape = KS_SHAPE_ARRAY[ptBinKaya];

  double prob_norm  = 0.0;//this gets calculated later.

  hist1 = new TH1D("hist1","hist1",1,0.4,1.4);
  hist2 = new TH1D("hist2","hist2",1,0.6,1.6);

  hist1->SetMarkerStyle(20);
  hist1->SetMarkerColor(2);
  hist1->SetLineColor  (2);

  hist2->SetMarkerStyle(20);
  hist2->SetMarkerColor(4);
  hist2->SetLineColor  (4);

  hist1->SetBinContent(1,rjz_pp);
  hist1->SetBinError  (1,rjz_pp_e*rjz_pp);

  hist2->SetBinContent(1,rjz_hi);
  hist2->SetBinError  (1,rjz_hi_e*rjz_hi);

  // seems unnecessary... 
  //hist1->Sumw2();
  //hist2->Sumw2();

  // this was just a trivial cross-check
  //double delta = (rjz_pp-rjz_hi);
  //double errsq = rjz_pp_e*rjz_pp_e + rjz_hi_e*rjz_hi_e;
  //double chisq = delta*delta/errsq;
  double mychi2 = CalcChiSq(hist1,hist2,1,1);
  //cout << chisq  << "  " << mychi2 << endl;
  if(verbosity==1)  cout << mychi2 << "  with 1 dof, has a p-value = " << TMath::Prob(mychi2,1) << endl;

  TF1 *f_chisq = new TF1("f_chisq","[1]*(pow(x,[0]/2-1)*TMath::Exp(-x/2))/(pow(2,[0]/2)*TMath::Gamma([0]/2))",0,20);
  f_chisq->SetParameter(1,1);
  f_chisq->SetParameter(0,1);//NDF
  f_chisq->SetLineColor(1);
  f_chisq->SetLineWidth(2);
  f_chisq->SetLineStyle(2);
  f_chisq->GetXaxis()->SetTitle("#chi^{2}, ndf=1");

  prob_norm  = TMath::Prob(mychi2,1);

  TH1D *h_frame = new TH1D("h_frame"," ",20,0,2);
  h_frame  ->GetYaxis()->SetRangeUser(0,1);
  h_frame  ->GetXaxis()->SetRangeUser(0,2);
  h_frame  ->GetXaxis()->SetTitle(" ");
  h_frame  ->GetYaxis()->SetTitle(" ");

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1",500,1000);
  c1->Divide(1,2);
  c1->cd(1);
  h_frame->DrawCopy();
  hist1->Draw("PE,same");
  hist2->Draw("PE,same");
  c1->cd(2);
  f_chisq->Draw();

  TLine *line_data = new TLine(mychi2,0,mychi2,1.2);
  line_data->SetLineColor(2);
  line_data->SetLineStyle(2);
  line_data->SetLineWidth(2);
  line_data->Draw("same");

  //cross-check the built-in ROOT function, TMath::Prob...
  if(verbosity==1)  cout << "partial/full integrals: " << f_chisq->Integral(mychi2,50) << "/" << f_chisq->Integral(0,50) << " = " << f_chisq->Integral(mychi2,50)/f_chisq->Integral(0,50) << endl;

  double prob_total = prob_norm*prob_shape*(1-TMath::Log(prob_norm*prob_shape));
  TF2 *fcor = new TF2("fcor","x*y*(1-TMath::Log(x*y))",0.0001,1.0,0.0001,1.0);

  if(verbosity==1)  cout << "With a shape probability (from KS) of " << prob_shape << endl;
  if(verbosity==1)  cout << " and this normalization probability of " << prob_norm << endl;
  if(verbosity==1)  cout << " we follow Eadie et al., section 11.6.2, a la ROOT " << endl;
  cout << " the combined p-value = " << prob_total << endl;


}
