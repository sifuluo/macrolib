double myfunction(double *xx, double *par)
{
//add
  Float_t x =xx[0];
  Double_t res = par[0] * (
    par[5]    * TMath::Gaus(x,par[1], par[2]) +
    (1.0-par[5]) * TMath::Gaus(x,par[1], par[2]*(1+par[6]))    )
    * ROOT::Math::normal_cdf(x,par[4],par[3]);
    return res;
}
std::pair<double, double> myfit(const int ieta, const int ipt, double nsigma)
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

  std::vector<const int> ptbins =
  {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 20, 23, 27, 30, 35, 40, 45, 57, 72, 90,
  120, 150, 200, 300, 400, 550, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000};

  const int pt1 = ptbins[ipt];
  const int pt2 = ptbins[ipt+1];

  string relrspname   = Form("ak4pfchsl1/RelRsp_JetEta%sto%s_RefPt%dto%d",eta_boundaries[ieta+41].c_str(),eta_boundaries[ieta+42].c_str(),pt1,pt2);
  string jetptname    = Form("ak4pfchsl1/JetPt_JetEta%sto%s_RefPt%dto%d",eta_boundaries[ieta+41].c_str(),eta_boundaries[ieta+42].c_str(),pt1,pt2);
  cout << relrspname <<endl;
  TFile *fin = new TFile("jra.root");
  TH1F *h1 = (TH1F*)fin->Get(relrspname.c_str());
  TH1F *h2 = (TH1F*)fin->Get(jetptname.c_str());
  double jetpt = h2->GetMean();


  TF1* fitfnc(0);
  const double ptmin = pt1;

  double integr=h1->Integral(); //0
  double mean=h1->GetMean(); //1
  const double rms=h1->GetRMS(); //2

  TSpectrum *spec = new TSpectrum(4);
  spec->Search(h1,6,"nobackground nodraw goff");
  Double_t* xpos = spec->GetPositionX();
  Double_t p = xpos[0];
  mean = p;
  double fitsigma = nsigma * rms;

  // double norm = h1->GetMaximumStored();
  struct para {
    const string name; const double ival; const double minval; const double maxval;
  };
  std::vector<para> paras = {
    {"norm      ", integr   , 0.                         , 2*integr       },
    {"core_mu   ", mean     , 0.85*mean                  , 1.2*mean       },
    {"core_sigma", rms      , 0.4*rms                    , 1.5*rms        },
    {"cdf_mu    ", 2./ptmin , 0.2/ptmin                  , 5./ptmin       },
    {"cdf_sigma ", 0.2/ptmin, 0.0001/ptmin               , 1./ptmin       },
    {"core_frac ", 0.9      , min(0.9,(ptmin-1)*0.05+0.4), 1.0            },
    {"sigma_inc ", 0.5*rms  , 0.02                       , 8.             },
  };

  std::vector<double> paravals;
  for (int ipush = 0; ipush < paras.size(); ++ipush){
    paravals.push_back(paras[ipush].ival);
  }


  // double xmin = h1->GetXaxis()->GetXmin();
  double xmin = 0.8/(ptmin+1.);

  if (pt1 > 1) xmin = max(paravals[1]-1.0*rms,xmin);
  double xmax = min(h1->GetXaxis()->GetXmax(),paravals[1] + 2.0*rms);

  for (int iiter = 0; iiter < 15; iiter++){

    fitfnc = new TF1("multigaus",myfunction,xmin,xmax,7);
    fitfnc->SetParNames("N","Core #mu","Core #sigma","CDF #mu","CDF #sigma","Tail Frac","Tail #sigma");

    for(int iset = 0; iset < paras.size(); ++iset){
      fitfnc->SetParameter(iset, paravals[iset]);
      fitfnc->SetParLimits(iset, paras[iset].minval, paras[iset].maxval);
    }
    h1->Fit(fitfnc,"RQ0");
    // delete fitfnc;
    // fitfnc = h1->GetFunction("multigaus");

    for (int iget = 0; iget < paras.size(); ++iget) {
      paravals[iget] = fitfnc->GetParameter(iget);
    }
  }
  h1->GetXaxis()->SetRangeUser(0,xmax);
  h1->Draw();
  fitfnc->Draw("same");

  //median
  double medianx, medianq;
  medianq = 0.5;
  h1->GetQuantiles(1, &medianx, &medianq);
  TLine *l1 = new TLine(medianx,0,medianx,h1->GetMaximum());
  l1->Draw("same");

  jetpt *= paravals[1]/h1->GetMean();

  for (int iprint = 0; iprint < paras.size(); ++iprint){
    cout << Form("%s = %10.5g => %10.5g With Limits [ %10.5g , %10.5g ]",paras[iprint].name.c_str(), paras[iprint].ival, paravals[iprint], paras[iprint].minval, paras[iprint].maxval) << endl;
    if(paravals[iprint] < paras[iprint].minval * 1.01) cout << "*********Hitting the LOWER limit*********" <<endl;
    if(paravals[iprint] > paras[iprint].maxval * 0.99) cout << "*********Hitting the UPPER limit*********" <<endl;
  }
  double chi2 = fitfnc->GetChisquare();
  double ndf  = fitfnc->GetNDF();
  cout << "median = " << medianx <<endl;
  cout << "Chi Square / NDF = "  << chi2 <<"/"<< ndf <<" = "<<chi2 / ndf<< endl;
  cout << "P = "<<fitfnc->GetProb()<<endl;
  cout << "JetPt = "<< h2->GetMean() << " * " << paravals[1]/h1->GetMean() << " = " << jetpt <<endl;
  cout << endl;

  std::pair<double, double> corrdots = std::make_pair(jetpt,1. / paravals[1]);

  return corrdots;

}

void fitit(const int ieta = 0,const double nsigma = 1.5)
{
  TCanvas *c1 = new TCanvas("c1","c1",3000,1500);
  int xpad = 5, ypad = 4;
  c1->Divide(xpad, ypad);

  TGraph *gr = new TGraph();
  for (int ipad = 1; ipad <= xpad* ypad; ++ipad ){
    c1->cd(ipad);
    pair<double, double> corrdots = myfit(ieta, ipad - 1, nsigma);
    // gr->SetPoint(gr->GetN(),ipad+2, 1 / myfit(ieta, ipad, nsigma);
    gr->SetPoint(gr->GetN(),corrdots.first,corrdots.second);

  }
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c2->cd();
  gr->Draw();
}
