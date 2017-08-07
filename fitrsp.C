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
double dashedfunc(double *xx, double *par)
{
  Float_t x = xx[0];
  Double_t res = par[0] * (
    par[5]    * TMath::Gaus(x,par[1], par[2]) +
    (1.0-par[5]) * TMath::Gaus(x,par[1], par[2]*(1+par[6]))    )+par[4]+par[3]-par[4]-par[3];
    return res;
}

std::vector<double> myfit(string eta1, string eta2, int pt1, int pt2, double nsigma, bool moreinfo)
{
  string relrspname   = Form("ak4pfchsl1/RelRsp_JetEta%sto%s_RefPt%dto%d",eta1.c_str(),eta2.c_str(),pt1,pt2);
  string jetptname    = Form("ak4pfchsl1/JetPt_JetEta%sto%s_RefPt%dto%d",eta1.c_str(),eta2.c_str(),pt1,pt2);
  if (moreinfo) cout << relrspname <<endl;
  TFile *fin = new TFile("jra.root");
  TH1F *h1 = (TH1F*)fin->Get(relrspname.c_str());
  TH1F *h2 = (TH1F*)fin->Get(jetptname.c_str());
  double jetpt = h2->GetMean();

  // h1->Rebin(10);

  TF1* fitfnc(0);
  const double ptmin = pt1;

  double integr=h1->Integral() * h1->GetBinWidth(0); //0
  integr = h1->GetMaximum();
  double mean=h1->GetMean(); //1
  const double rms=h1->GetRMS(); //2

  TSpectrum *spec = new TSpectrum(4);
  spec->Search(h1,6,"nobackground nodraw goff");
  Double_t* xpos = spec->GetPositionX();
  Double_t p = xpos[0];
  mean = p;
  double fitsigma = nsigma * rms;

  struct para {
    const string name; const double ival; const double minval; const double maxval;
  };
  std::vector<para> paras = {
    {"norm      ", integr   , 0.                         , 2*integr       },
    {"core_mu   ", mean     , 0.7*mean                  , 1.2*mean       },
    {"core_sigma", rms      , 0.4*rms                    , 1.5*rms        },
    {"cdf_mu    ", 3.5/ptmin , 2.0/ptmin                  , 5./ptmin       },
    {"cdf_sigma ", 0.2/ptmin, 0.1/ptmin               ,  1./ptmin       },
    {"core_frac ", 0.9      , min(0.4,(ptmin-2)*0.05+0.4), 1.0            },
    {"sigma_inc ", 0.5  , 0.02                       , 2.             },
  };

  std::vector<double> paravals;
  for (int ipush = 0; ipush < paras.size(); ++ipush){
    paravals.push_back(paras[ipush].ival);
  }
  std::vector<double> paralims;
  paralims.resize(paravals.size(),0);

  // double xmin = h1->GetXaxis()->GetXmin();
  double xmin = 0.8/(ptmin+1.);

  if (pt1 > 0) xmin = max(paravals[1]-0.5*rms,xmin);
  double xmax = min(h1->GetXaxis()->GetXmax(),paravals[1] + 2.0*rms);

  for (int iiter = 0; iiter < 15; iiter++){

    fitfnc = new TF1("multigaus",myfunction,xmin,xmax,7);
    fitfnc->SetParNames("N","Core #mu","Core #sigma","CDF #mu","CDF #sigma","Tail Frac","Tail #sigma");

    for(int iset = 0; iset < paras.size(); ++iset){
      fitfnc->SetParameter(iset, paravals[iset]);
      fitfnc->SetParLimits(iset, paras[iset].minval, paras[iset].maxval);
    }
    h1->Fit(fitfnc,"RQ0");

    for (int iget = 0; iget < paras.size(); ++iget) {
      paravals[iget] = fitfnc->GetParameter(iget);
    }

  }

  // h1->GetXaxis()->SetRangeUser(0,xmax);
  h1->Draw();
  fitfnc->Draw("same");
  TString namec="h1c_";
  namec.Append(h1->GetName());
  TH1F *h1c=(TH1F*)h1->Clone(namec);
  for(unsigned i=0; i<h1c->GetNbinsX(); ++i) h1c->SetBinContent(i+1,h1c->GetBinContent(i+1)-fitfnc->Eval(h1c->GetBinCenter(i+1)));
  h1c->SetLineColor(8);
  h1c->SetLineStyle(2);
  h1c->Draw("same,hist");
  TF1 *gaustest=new TF1("gaustest","gausn",0.,1.);
  gaustest->SetLineColor(1);
  gaustest->SetLineStyle(3);
  h1c->Fit(gaustest,"R");
  gaustest->Draw("same");

  jetpt *= paravals[1]/h1->GetMean();
  double chi2 = fitfnc->GetChisquare();
  double ndf  = fitfnc->GetNDF();
  for (int iprint = 0; iprint < paras.size(); ++iprint){
    if(paravals[iprint] < paras[iprint].minval * 1.01) {
      paralims[iprint] = 1;
    }
    if(paravals[iprint] > paras[iprint].maxval * 0.99) {
      paralims[iprint] = 2;
    }
  }
  //compose the output values
  std::vector<double> output;
  for (int iout = 0; iout < paravals.size(); ++iout) {
    output.push_back(paravals[iout]);// (ipara-1)*5
    output.push_back(paras[iout].ival); // (ipara -1)*5 +1
    output.push_back(paras[iout].minval);// (ipara -1)*5 +2
    output.push_back(paras[iout].maxval); //(ipara -1)*5 +3
    output.push_back(paralims[iout]);      //(ipara -1)*5 +4
  }
  output.push_back(chi2);   // npara * 5
  output.push_back(ndf);    // npara * 5 + 1
  output.push_back(jetpt);  // npara * 5 + 2
  output.push_back((double)pt1);    // npara * 5 + 3
  output.push_back(pt2);            // npara * 5 + 4

  //needed output all set

  if(!moreinfo) return output;

  for (int iprint = 0; iprint < paras.size(); ++iprint){
    cout << Form("%s = %10.5g => %10.5g With Limits [ %10.5g , %10.5g ]",paras[iprint].name.c_str(), paras[iprint].ival, paravals[iprint], paras[iprint].minval, paras[iprint].maxval) << endl;
    if(paravals[iprint] < paras[iprint].minval * 1.01) cout << "*********Hitting the LOWER limit*********" <<endl;
    if(paravals[iprint] > paras[iprint].maxval * 0.99) cout << "*********Hitting the UPPER limit*********" <<endl;
  }
  cout << "Threshold = " << paravals[3]*pt1 <<endl;

  // cout << "median = " << medianx <<endl;
  // cout << "median without Threshold Effect = " << medianx2 <<endl;
  cout << "Chi Square / NDF = "  << chi2 <<"/"<< ndf <<" = "<<chi2 / ndf<< endl;
  cout << "P = "<<fitfnc->GetProb()<<endl;
  cout << "JetPt = "<< h2->GetMean() << " * " << paravals[1]/h1->GetMean() << " = " << jetpt <<endl;
  cout << endl;

  // std::pair<double, double> corrdots = std::make_pair(jetpt,1. / paravals[1]);
  // std::pair<double, double> corrdots = std::make_pair(jetpt,1. / paravals[1]);



  return output;

}

void fitrsp(const int ieta = 0,const double nsigma = 1.5)
{
  cout << "============================================================================================================================================"<< endl;
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

  TCanvas *c1 = new TCanvas("c1","c1",3000,1500);
  int xpad = 5, ypad = 4;
  c1->Divide(xpad, ypad);
  bool moreinfo = 0;

  TGraph *gr = new TGraph();
  std::vector<std::vector<double> > outputs;
  string eta1 = eta_boundaries[ieta + 41];
  string eta2 = eta_boundaries[ieta + 42];
  int npara = 7;

  for (int ipad = 1; ipad <= xpad* ypad; ++ipad ){
    int ipt = ipad - 1;
    int pt1 = ptbins[ipt];
    int pt2 = ptbins[ipt+1];

    c1->cd(ipad);
    outputs.push_back(myfit(eta1,eta2,pt1,pt2,nsigma,moreinfo));
    // gr->SetPoint(gr->GetN(),ipad+2, 1 / myfit(ieta, ipad, nsigma);
  }
  if (moreinfo){
    cout <<endl;
    cout << "---------------------------------------------------------------------------------------------------------------------" <<endl<<endl;
  }
  cout << "Summary on Response of eta Range = [ " << eta1 << " , " << eta2 << " ]" <<endl;
  cout << "  RefPt"<<"       Norm "<<"     core_mu "<<"  core_sigma "<<"      cdf_mu "<<"   cdf_sigma "<<"   core_frac "<<"   sigma_inc " <<"    Chi2/NDF " <<" Threshold"<<endl;


  for (int iout = 0; iout < outputs.size(); ++iout){
    std::vector<double> printvec = outputs.at(iout);
    cout << Form("%3d-%3d ",int(printvec.at(npara*5+3)),int(printvec.at(npara*5+4)));
    for (int iout1 = 0; iout1 < npara; ++iout1){
      string limflag = " ";
      if (printvec.at(iout1*5+4) == 1) limflag = "L";
      if (printvec.at(iout1*5+4) == 2) limflag = "H";
      cout << Form("%10.3f %s ",printvec.at(iout1*5),limflag.c_str());
    }
    cout << Form("%10.3f %10.3f",printvec.at(npara*5) / printvec.at(npara*5+1), printvec.at(3*5) * printvec.at(npara*5+3));
    cout << endl;
    if (moreinfo) cout <<endl;

    double xaxis = printvec.at(npara*5+2);
    double yaxis = 1/printvec.at(1*5);
    gr->SetPoint(gr->GetN(),xaxis,yaxis);
  }

  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c2->cd();
  gr->Draw("*al");
}
