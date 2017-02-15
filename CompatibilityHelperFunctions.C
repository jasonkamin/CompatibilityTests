double CalcChiSq(TH1D *h1, TH1D *h2, int firstBin, int lastBin);
double CalcKoloSm(TH1D *h1, TH1D *h2, int firstBin, int lastBin, int option, int ikolo);



double CalcChiSq(TH1D *h1, TH1D *h2, int firstBin, int lastBin)
{

  double chisq = 0.0;
  double delta = 0.0;
  double errsq = 0.0;

  if(h1->GetNbinsX() != h2->GetNbinsX()){
    cout << "JJ::CalcChiSq() --> histograms don't have same number of bins !!" << endl;
    return 0;
  }

  if(firstBin<0)
    firstBin=1;
  if(lastBin<0)
    lastBin=h1->GetNbinsX();

  for(int i=firstBin; i<lastBin+1; i++){
    double h1cont = h1->GetBinContent(i);
    double h2cont = h2->GetBinContent(i);
    double h1err  = h1->GetBinError  (i);
    double h2err  = h2->GetBinError  (i);

    if(h1cont !=0 && h2cont !=0){
      delta = (h1cont-h2cont);
      errsq = h1err*h1err + h2err*h2err;
      chisq += delta*delta/errsq;
    }

  }

  return chisq;
}

double CalcKoloSm(TH1D *h1, TH1D *h2, int firstBin, int lastBin, int option, int ikolo)
{

  //
  // options: 
  // 1: normalized (this is the standard way!)
  // 2: un-normalized
  //

  char saykolo[500];
  double dfmax  = 0.0;
  double dfmaxN = 0.0;
  double sum1   = 0.0;
  double sum2   = 0.0;

  if(ikolo>=0){
    sprintf(saykolo,"h_KoloHist_%d",ikolo);
    h_KoloHist[ikolo] = (TH1D*)h1->Clone(saykolo);
    h_KoloHist[ikolo]->Reset();
  }

  if(h1->GetNbinsX() != h2->GetNbinsX()){
    cout << "JJ::CalcKoloSm() --> histograms don't have same number of bins !!" << endl;
    return 0;
  }

  if(firstBin<0)
    firstBin=1;
  if(lastBin<0)
    lastBin=h1->GetNbinsX();

  for(int i=firstBin; i<lastBin+1; i++)
    sum1 += h1->GetBinContent(i);
  for(int i=firstBin; i<lastBin+1; i++)
    sum2 += h2->GetBinContent(i);

  double h1cont  = 0.0;
  double h2cont  = 0.0;
  double h1err   = 0.0;
  double h2err   = 0.0;
  for(int i=firstBin; i<lastBin+1; i++){
    h1cont  += h1->GetBinContent(i);
    h2cont  += h2->GetBinContent(i);
    h1err   = h1->GetBinError  (i);
    h2err   = h2->GetBinError  (i);
    double mydiff  = fabs(h1cont-h2cont);
    double mydiffN = fabs(h1cont/sum1-h2cont/sum2);

    if(ikolo>=0){
      h_KoloHist[ikolo]->SetBinContent(i,mydiffN);
      h_KoloHist[ikolo]->SetBinError  (i,0.0);
    }

    if(h1cont !=0 && h2cont !=0){
      if(mydiff>dfmax)
        dfmax = mydiff;
      if(mydiffN>dfmaxN)
        dfmaxN = mydiffN;
    }

  }

  if(option==1)      return dfmaxN;
  else if(option==2) return dfmax;
  else               return -1;
}
