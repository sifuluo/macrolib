{
  // TH1F* hmean = new TH1F("mean","mean",200,-1,1);
  TH1F* hmedian = new TH1F("median","median",200,-1,1);
  TH1F* hmederror = new TH1F("median_error","median error",10000,0,0.1);
  for (int i = 0; i < 200 ; ++i){
    TH1F* h1 = new TH1F("h1","h1",60,-0.3,0.3);
    h1->FillRandom("gaus",1000);
    Double_t x, q;
    q = 0.5; // 0.5 for "median"
    // h1->ComputeIntegral(); // just a precaution
    h1->GetQuantiles(1, &x, &q);
    hmedian->Fill(x);
    // hmean->Fill(h1->GetMean());
    double error = TMath::Sqrt(h1->Integral()/(4 *pow(h1->GetBinContent(h1->FindBin(x)) / h1->GetBinWidth(1),2)));
    // double error = 1.0 / TMath::Sqrt(h1->Integral("width") * 4 * pow(h1->GetBinContent(h1->FindBin(x)),2));
    hmederror->Fill(error);
    delete h1;
  }
  // hmean->Draw();
  // new TCanvas();
  hmedian->Draw();
  new TCanvas();
  hmederror->Draw();

}
