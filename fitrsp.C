double myfunction(double *xx, double *par)
{
  Float_t x =xx[0];
  Double_t res = par[0] * (
    par[5]    * TMath::Gaus(x,par[1], par[2]) +
    (1.0-par[5]) * TMath::Gaus(x,par[1], par[2]+par[6])    )
    * ROOT::Math::normal_cdf(x,par[4],par[3]);
    return res;
}
double myfit(const int ieta, const int pt)
{
  // std::vector<double> etabins = {0., 0.087, 0.174, 0.261, 0.348, 0.435 ,0.522, 0.609, 0.696,
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


  const int pt1 = ptbins[pt-3];
  const int pt2 = ptbins[pt-2];

  string fullname = Form("ak4pfchsl1/RelRsp_JetEta%sto%s_RefPt%dto%d",eta_boundaries[ieta+41].c_str(),eta_boundaries[ieta+42].c_str(),pt1,pt2);
  // string fullname = "ak4pfchsl1/RelRsp_JetEta" + eta1 + "to" + eta2 + "_RefPt" + pt1 + "to" + pt2;
  cout << fullname <<endl;
  TFile *fin = new TFile("jra.root");
  TH1F *h1 = (TH1F*)fin->Get(fullname.c_str());
  h1->Draw();
  TF1* fitfnc(0);
  const double ptmin = pt1;

  const double integr=h1->Integral(); //0
  const double mean=h1->GetMean(); //1
  const double rms=h1->GetRMS();  //2

  const double norm0 =integr; //0
  const double peak0 = mean;  //1
  const double sigma0 = rms;  //2
  const double cdfmu0 = 2./ ptmin; //3
  const double cdfsig0 = 0.2/ptmin; //4
  const double tfrac0 = 0.3;  //5
  const double tsigma0 = rms;  //6



  double norm[3]   = {norm0   , 0.          , 2*integr};
  double peak[3]   = {peak0   , 0.8*mean    , 1.2*mean};
  double sigma[3]  = {sigma0  , 0.2         , 1.5*rms};
  double cdfmu[3]  = {cdfmu0  , 0.2 / ptmin , 5. / ptmin};
  double cdfsig[3] = {cdfsig0 , 0.0001/ptmin, 1./ptmin};
  double tfrac[3]  = {tfrac0  , 0.1         , 1.};
  double tsigma[3] = {tsigma0 , 0.          , 4*rms};
  // double norm = h1->GetMaximumStored();

  // double xmin = h1->GetXaxis()->GetXmin();
  double xmin = 0.8/(ptmin+1.);
  double xmax = h1->GetXaxis()->GetXmax();

  for (int iiter = 0; iiter < 10; iiter++){

    fitfnc = new TF1("multigaus",myfunction,xmin,xmax,7);
    fitfnc->SetParNames("N","Core #mu","Core #sigma","CDF #mu","CDF #sigma","Tail Frac","Tail #sigma");
    fitfnc->SetParameter(0,norm[0]);
    fitfnc->SetParLimits(0,norm[1],norm[2]);

    fitfnc->SetParameter(1,peak[0]);
    fitfnc->SetParLimits(1,peak[1],peak[2]);

    fitfnc->SetParameter(2,sigma[0]);
    fitfnc->SetParLimits(2,sigma[1],sigma[2]);

    fitfnc->SetParameter(3,cdfmu[0]);
    fitfnc->SetParLimits(3,cdfmu[1],cdfmu[2]);

    fitfnc->SetParameter(4,cdfsig[0]);
    fitfnc->SetParLimits(4,cdfsig[1],cdfsig[2]);

    fitfnc->SetParameter(5,tfrac[0]);
    fitfnc->SetParLimits(5,tfrac[1],tfrac[2]);

    fitfnc->SetParameter(6,tsigma[0]);
    fitfnc->SetParLimits(6,tsigma[1],tsigma[2]);
    h1->Fit(fitfnc,"RQ0");
    // delete fitfnc;
    // fitfnc = h1->GetFunction("multigaus");
    // if (fitfnc) cout << iiter << " time fit successful" <<endl;
    norm[0]  = fitfnc->GetParameter(0);
    peak[0]  = fitfnc->GetParameter(1);
    sigma[0] = fitfnc->GetParameter(2);
    cdfmu[0] = fitfnc->GetParameter(3);
    cdfsig[0]= fitfnc->GetParameter(4);
    tfrac[0] = fitfnc->GetParameter(5);
    tsigma[0]= fitfnc->GetParameter(6);

  }
  fitfnc->Draw("same");
  cout << "norm = "   << norm0   <<" => "<< norm[0]   << "  With Limits  [  "<< norm[1]   <<"  ,  "<< norm[2]   <<"  ]"<<endl;
  cout << "peak = "   << peak0   <<" => "<< peak[0]   << "  With Limits  [  "<< peak[1]   <<"  ,  "<< peak[2]   <<"  ]"<<endl;
  cout << "sigma = "  << sigma0  <<" => "<< sigma[0]  << "  With Limits  [  "<< sigma[1]  <<"  ,  "<< sigma[2]  <<"  ]"<<endl;
  cout << "cdfmu = "  << cdfmu0  <<" => "<< cdfmu[0]  << "  With Limits  [  "<< cdfmu[1]  <<"  ,  "<< cdfmu[2]  <<"  ]"<<endl;
  cout << "cdfsig = " << cdfsig0 <<" => "<< cdfsig[0] << "  With Limits  [  "<< cdfsig[1] <<"  ,  "<< cdfsig[2] <<"  ]"<<endl;
  cout << "tfrac = "  << tfrac0  <<" => "<< tfrac[0]  << "  With Limits  [  "<< tfrac[1]  <<"  ,  "<< tfrac[2]  <<"  ]"<<endl;
  cout << "tsigma = " << tsigma0 <<" => "<< tsigma[0] << "  With Limits  [  "<< tsigma[1] <<"  ,  "<< tsigma[2] <<"  ]"<<endl;
  cout << "Chi Square / NDF = "  << fitfnc->GetChisquare()<<"/"<<fitfnc->GetNDF()<<endl;
  cout << "P = "<<fitfnc->GetProb()<<endl;
  cout << endl;

  return peak[0];

}

void fitit(const int ieta)
{
  TCanvas *c1 = new TCanvas("c1","c1",2400,1200);
  c1->Divide(4,2);
  TGraph *gr = new TGraph();
  for (int ipad = 1; ipad <= 8; ++ipad ){
    c1->cd(ipad);
    double padpt = ipad + 2;
    gr->SetPoint(gr->GetN(),ipad+2, 1 / myfit(ieta, ipad+2));

  }
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c2->cd();
  gr->Draw();
}
