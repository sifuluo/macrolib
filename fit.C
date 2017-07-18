double myfunction(double *xx, double *par)
{
  Float_t x =xx[0];
  Double_t res = par[0] * (
    par[5]    * TMath::Gaus(x,par[1], par[2]) +
    (1.0-par[5]) * TMath::Gaus(x,par[1], par[2]+par[6])    )
    * ROOT::Math::normal_cdf(x,par[4],par[3]);
    return res;
}
void myfit(const int ieta, const int pt)
{
  // std::vector<float> etabins = {0., 0.087, 0.174, 0.261, 0.348, 0.435 ,0.522, 0.609, 0.696,
  // 0.783,  0.879, 0.957, 1.044, 1.131, 1.218, 1.305 ,1.392, 1.479, 1.566,
  // 1.653,  1.74,  1.83,  1.93,  2.043, 2.172, 2.322 ,2.5,   2.65,  2.853,
  // 2.964,  3.139, 3.314, 3.489, 3.664, 3.839, 4.013 ,4.191, 4.363, 4.538,
  // 4.716,4.889,5.191};
  const string eta_boundaries[] =
  {"-5.191","-4.889","-4.716","-4.538","-4.363","-4.191","-4.013","-3.839","-3.664","-3.489",
  "-3.314", "-3.139","-2.964","-2.853","-2.65", "-2.5",  "-2.322","-2.172","-2.043","-1.93",
  "-1.83",  "-1.74", "-1.653","-1.566","-1.479","-1.392","-1.305","-1.218","-1.131","-1.044",
  "-0.957", "-0.879","-0.783","-0.696","-0.609","-0.522","-0.435","-0.348","-0.261","-0.174",
  "-0.087", "0",     "0.087", "0.174", "0.261", "0.348", "0.435" ,"0.522", "0.609", "0.696",
  "0.783",  "0.879", "0.957", "1.044", "1.131", "1.218", "1.305" ,"1.392", "1.479", "1.566",
  "1.653",  "1.74",  "1.83",  "1.93",  "2.043", "2.172", "2.322" ,"2.5",   "2.65",  "2.853",
  "2.964",  "3.139", "3.314", "3.489", "3.664", "3.839", "4.013" ,"4.191", "4.363", "4.538",
  "4.716","4.889","5.191"};
  std::vector<const int> ptbins = {3, 4, 5, 6, 7, 8, 9, 10, 11};

  // string eta1 = eta_boundaries[ieta+41];
  // string eta2 = eta_boundaries[ieta+42];


  int pt1 = ptbins[pt-3];
  int pt2 = ptbins[pt-2];

  string fullname = Form("ak4pfchsl1/RelRsp_JetEta%sto%s_RefPt%dto%d",eta_boundaries[ieta+41].c_str(),eta_boundaries[ieta+42].c_str(),pt1,pt2);
  // string fullname = "ak4pfchsl1/RelRsp_JetEta" + eta1 + "to" + eta2 + "_RefPt" + pt1 + "to" + pt2;
  cout << fullname <<endl;
  TFile *fin = new TFile("jra.root");
  TH1F *h1 = (TH1F*)fin->Get(fullname.c_str());
  h1->Draw();
  TF1* fitfnc(0);
  double peak = h1->GetMean();
  double sigma = h1->GetRMS();
  // double norm = h1->GetMaximumStored();
  double norm = h1->GetEntries();
  double ptmin = pt1;
  cout << "ptmin = " <<ptmin <<endl;
  double xmin = h1->GetXaxis()->GetXmin();
  double xmax = h1->GetXaxis()->GetXmax();

  double cdfmu = 3.0/ptmin;
  double cdfsig = cdfmu/3.0;
  double tfrac = 0.3;
  double tsigma = 10*sigma;
  cout << "norm = " << norm <<endl;
  cout << "peak = " << peak <<endl;
  cout << "sigma = " << sigma <<endl;
  cout << "cdfmu = " << cdfmu <<endl;
  cout << "cdfsig = " << cdfsig <<endl;
  cout << "tfrac = " << tfrac <<endl;
  cout << "tsigma = " << tsigma <<endl;

  for (int iiter = 0; iiter < 10; iiter++){
    const float integr=h1->Integral();
    const float mean=h1->GetMean();
    const float rms=h1->GetRMS();
    fitfnc = new TF1("multigaus",myfunction,xmin,xmax,7);
    fitfnc->SetParNames("N","Core #mu","Core #sigma","CDF #mu","CDF #sigma","Tail Frac","Tail #sigma");
    fitfnc->SetParameter(0,norm);
    fitfnc->SetParLimits(0,0,2*integr);
    fitfnc->SetParameter(1,peak);
    fitfnc->SetParLimits(1,0.8*mean,1.2*mean);
    fitfnc->SetParameter(2,sigma);
    fitfnc->SetParLimits(2,0.0,1.5*rms);
    fitfnc->SetParameter(3,cdfmu);
    //fitfnc->SetParLimits(3,2./ptmin,4./ptmin);
    fitfnc->SetParameter(4,cdfsig);
    //fitfnc->SetParLimits(4,);
    fitfnc->SetParameter(5,tfrac);
    fitfnc->SetParLimits(5,0.5,1);
    fitfnc->SetParameter(6,rms);
    fitfnc->SetParLimits(6,0,rms*4);
    h1->Fit(fitfnc,"RQ0");
    // delete fitfnc;
    // fitfnc = h1->GetFunction("multigaus");
    // if (fitfnc) cout << iiter << " time fit successful" <<endl;
    norm  = fitfnc->GetParameter(0);
    peak  = fitfnc->GetParameter(1);
    sigma = fitfnc->GetParameter(2);
    cdfmu = fitfnc->GetParameter(3);
    cdfsig= fitfnc->GetParameter(4);
    tfrac = fitfnc->GetParameter(5);
    tsigma= fitfnc->GetParameter(6);

  }
  fitfnc->Draw("same");
  cout << "norm = " << norm <<endl;
  cout << "peak = " << peak <<endl;
  cout << "sigma = " << sigma <<endl;
  cout << "cdfmu = " << cdfmu <<endl;
  cout << "cdfsig = " << cdfsig <<endl;
  cout << "tfrac = " << tfrac <<endl;
  cout << "tsigma = " << tsigma <<endl;

}

void fitit(const int ieta)
{
  TCanvas *c1 = new TCanvas("c1","c1",2400,1200);
  c1->Divide(4,2);
  for (int ipad = 1; ipad <= 8; ++ipad ){
    c1->cd(ipad);
    double padpt = ipad + 2;
    myfit(ieta, padpt);
  }
}
