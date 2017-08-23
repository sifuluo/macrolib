// {
//   // TH1F* hmean = new TH1F("mean","mean",200,-1,1);
//   TH1F* hmedian = new TH1F("median","median",200,-1,1);
//   TH1F* hmederror = new TH1F("median_error","median error",10000,0,1);
//   TH1F* h1 = new TH1F("h1","h1",600,-5,5);
//   for (int i = 0; i < 200 ; ++i){
//     h1->Reset();
//     h1->FillRandom("gaus",1000);
//     Double_t x, q;
//     q = 0.5; // 0.5 for "median"
//     // h1->ComputeIntegral(); // just a precaution
//     h1->GetQuantiles(1, &x, &q);
//     hmedian->Fill(x);
//     // hmean->Fill(h1->GetMean());
//     // double error = TMath::Sqrt(h1->Integral()/(4 *pow(h1->GetBinContent(h1->FindBin(x)) / h1->GetBinWidth(1),2)));
//     // double error = 1.0 / TMath::Sqrt(h1->Integral("width") * 4 * pow(h1->GetBinContent(h1->FindBin(x)),2));
//     // double error = 1.0 / TMath::Sqrt(h1->Integral("width") * 4 * pow(h1->GetBinContent(h1->FindBin(x))* h1->GetBinWidth(1),2));
//     double error = 1.0 / TMath::Sqrt(4 * h1->Integral() * pow(h1->GetBinContent(h1->FindBin(x)) *   h1->GetBinWidth(1)/ h1->Integral(),2) );
//     hmederror->Fill(error);
//   }
//   // hmean->Draw();
//   // new TCanvas();
//   hmedian->Draw();
//   new TCanvas();
//   hmederror->Draw();
//
// }

{
  // TH1F* hmean = new TH1F("mean","mean",200,-1,1);
  TH1F* hmedian = new TH1F("median","median",6000,0,6);
  TH1F* hmederror = new TH1F("median_error","median error",1000,0,0.1);
  TH1F* h1 = new TH1F("h1","h1",6000,0,6);

  // Response is between 0 and 6 with a mean around 1
  TF1 * f1 = new TF1("f1","TMath::Gaus(x,1,0.5)",0,6);

  for (int i = 0; i < 200 ; ++i){
    h1->Reset();
    h1->FillRandom("f1",2000);
    Double_t x, q;
    q = 0.5; // 0.5 for "median"
    // h1->ComputeIntegral(); // just a precaution
    h1->GetQuantiles(1, &x, &q);
    hmedian->Fill(x);
    // hmean->Fill(h1->GetMean());
    // double error = TMath::Sqrt(h1->Integral()/(4 *pow(h1->GetBinContent(h1->FindBin(x)) / h1->GetBinWidth(1),2)));
    // double error = 1.0 / TMath::Sqrt(h1->Integral("width") * 4 * pow(h1->GetBinContent(h1->FindBin(x)),2));
    double n = h1->Integral();
    int bin_median = h1->FindBin(x);
    double f_median = h1->GetBinContent(bin_median) /( h1->GetBinWidth(bin_median) * n);
    double median_error = 1.0 / TMath::Sqrt(n * 4.0 * pow(f_median,2));
    cout<<"median error="<<median_error<<endl;
    hmederror->Fill(median_error);
  }

  // hmean->Draw();
  TCanvas * c = new TCanvas("Median Test","MT", 2400,800);
  c->Divide(3,1);
  c->cd(1);
  h1->Draw();

  c->cd(2);
  hmedian->Draw();

  c->cd(3);
  hmederror->Draw();

}
