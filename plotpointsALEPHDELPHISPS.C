 {

  gROOT->SetStyle("Plain");


    TFile * finH=TFile::Open("ADistribfuncMBmill_eenfjnewzbin.root","READ");
  //  TFile * finH = TFile::Open("ADistribfuncMBmill_eenfjrespnoSO.root","READ");
    TH1F* hzscaffid[10];
    TH1F* hzscaffid2[10];
    TH1D* hevents=(TH1D*)finH->Get("hevents");
 
    TH1D* heventsid=(TH1D*)finH->Get("heventsid");

    Int_t neventswpi=heventsid->GetBinContent(2);
    Int_t neventswk=heventsid->GetBinContent(3);


  for(Int_t i=0; i<10;i++){
    hzscaffid[i]=(TH1F*)finH->Get(Form("hzscaffid%d",i));
    hzscaffid2[i]=(TH1F*)hzscaffid[i]->Clone(Form("hzscaffid%d",i));
  }
   Int_t nevents=hevents->GetBinContent(10);
 //  Int_t neventsonevalid=hevents->GetBinContent(2);

   cout<<"neventos="<<nevents<<endl;
  Double_t kHx[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHy[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHe[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHex[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHx[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHy[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHe[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHex[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

/*
  Double_t kHx[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHy[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHe[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHx[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHy[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHe[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
*/ 
   for(Int_t piden=1; piden<3; piden++){
    Int_t nbinzscaid=hzscaffid[1]->GetXaxis()->GetNbins();
    cout<<"nbins= "<<nbinzscaid<<endl;
    Int_t nbinzscaid2=hzscaffid2[1]->GetXaxis()->GetNbins();
    Double_t contzscaid=0;
    Double_t errzscaid=0;
    for(Int_t i=1; i<nbinzscaid+3; i++){
      cout<<"----i="<<i<<endl;
      Double_t zscaid=hzscaffid[piden]->GetBinCenter(i);
      if(piden==1)piHx[i]=zscaid;
      if(piden==2)kHx[i]=zscaid;
      Double_t zscawidthid=hzscaffid[piden]->GetBinWidth(i);
      for(Int_t j=1; j<nbinzscaid2+3; j++){
       Double_t zscaid2=hzscaffid2[piden]->GetBinCenter(j);
       Double_t zscawidthid2=hzscaffid2[piden]->GetBinWidth(j);
       if(i==j){
         cout<<"----j="<<j<<endl;
       Double_t contzscaid=hzscaffid[piden]->GetBinContent(i);
       Double_t errzscaid=hzscaffid[piden]->GetBinError(i);
   cout<<"cont= "<<contzscaid<<", zscawidthid2="<<zscawidthid2<<endl;
       hzscaffid2[piden]->SetBinContent(j,contzscaid/zscawidthid2);
       hzscaffid2[piden]->SetBinError(j,errzscaid/zscawidthid2);
       Double_t ratioc=contzscaid/(zscawidthid2*nevents); //*nevents); // /(zscawidthid2*80000); ///(zscawidthid2*nevents); //zscawidthid2
       Double_t ratioe=errzscaid/(zscawidthid2*nevents); //*nevents);  // /(zscawidthid2*80000); /// (zscawidthid2*nevents);
       cout<<"ratioc="<<ratioc<<", ratioe="<<ratioe<<endl;
       if(piden==2){
       kHy[i]=ratioc; ///neventswk; // /80000; // /(80000*0.05); // /55960;
       kHe[i]=ratioe; ///neventswk; // /80000; // /(80000*0.05); // /55960;
       kHex[i]=zscawidthid2/2;
       cout<<"kHx="<<kHx[i]<<", kHy="<<kHy[i]<<", kHe="<<kHe[i]<<endl;
       }
       if(piden==1){
       piHy[i]=ratioc; ///neventswpi; // /80000; // /(80000*0.05); // /35380;
       piHe[i]=ratioe; ///neventswpi; // /80000; // /(80000*0.05); // /35380;
       piHex[i]=zscawidthid2/2;
       cout<<"piHx="<<piHx[i]<<", piHy="<<piHy[i]<<", piHe="<<piHe[i]<<endl;
       }

       }
      }

     }
//     hzscaffid2[piden]->Scale(1./(nevents));
    } //end loop of pid

 /*
  for(Int_t i=0; i<nbinzscaid+1; i++){
    cout<<"x="<<kHx[i]<<"y="<<kHy[i]<<endl;
  }
 */
  int kudsdata_np = 21;
  double kudsdataset1z[]={0.0190,0.0245,0.0300,0.0355,0.0410,0.0465,0.0520,
               0.0605,0.0715,0.0825,0.0935,0.1045,0.1210,0.1480,
               0.2190,0.2520,0.2960,0.3510,0.4275,0.5370,0.6855};
  double kudsdataset1sig[]={22.600,19.200,18.600,17.000,14.600,15.300,14.500,
                  10.290, 9.000, 7.380, 6.120, 6.000, 4.780, 3.300,
                   2.290, 1.498, 1.272, 0.925, 0.548, 0.266, 0.101};
  double kudsdataset1stat[]={1.200,1.100,1.100,1.000,1.100,1.200,1.300,0.850,
                   0.730,0.700,0.700,0.750,0.570,0.610,0.170,0.089,
                   0.068,0.046,0.032,0.020,0.015};
  // dataset1sys/21*0.0/
 /*
   int kdata_np = 21;
  double ktotdataset1z[]={
   };
  double ktotdataset1sig[]={
   };
  double ktotdataset1stat[]={
   };
  double ktotdataset1sys[]={
  };
 */
  
     int kudstotdata_np = 32;
  double kudstotdataset4z[]={0.0190,0.0245,0.0300,0.0355,0.0410,0.0465,0.0520,
                0.0575,0.0630,0.0685,0.0740,0.0795,0.0850,0.0935,
                0.1045,0.1155,0.1320,0.1535,0.2190,0.2410,0.2630,
                0.2850,0.3070,0.3345,0.3675,0.4005,0.4335,0.4660,
                0.5040,0.5480,0.6140,0.7130};
  double kudstotdataset4sig[]={22.280,21.220,19.100,18.020,16.450,16.130,13.980,
                  11.960,11.030,11.370, 9.380, 9.510, 8.680, 7.960,
                   7.370, 6.740, 5.240, 4.080, 2.340, 1.947, 1.603,
                   1.281, 1.127, 0.933, 0.669, 0.541, 0.493, 0.357,
                   0.273, 0.180, 0.122, 0.048};
  double kudstotdataset4stat[]={0.470,0.450,0.430,0.430,0.450,0.490,0.530,0.540,
                   0.490,0.460,0.440,0.440,0.440,0.310,0.340,0.370,
                   0.290,0.400,0.080,0.065,0.057,0.050,0.045,0.034,
                   0.029,0.026,0.026,0.023,0.019,0.016,0.011,0.009};
  double kudstotdataset4sys[]={0.530,0.620,0.640,0.800,0.940,1.030,1.140,2.110,
                  1.270,0.870,0.790,0.760,0.720,0.680,0.630,0.600,
                  0.610,0.800,0.300,0.057,0.042,0.034,0.030,0.027,
                  0.023,0.023,0.023,0.020,0.018,0.014,0.011,0.006};

     int kalephtotdata_np = 29;
  double kalephtotdataset1z[]={0.0057,0.0063,0.0067,0.0072,0.0077,0.0083,0.0088,
    		 0.0092,0.0135,0.0150,0.0170,0.0725,0.0775,0.0825,
     		 0.0875,0.0950,0.1050,0.1150,0.1250,0.1350,0.1450,
     		 0.1550,0.1700,0.1900,0.2250,0.2750,0.3500,0.5000,
     		 0.7000};
  double kalephtotdataset1sig[]={12.400,13.270,15.330,17.430,18.330,19.620,20.020,
    		 21.660,25.840,27.460,27.630,10.600,9.530,9.150,
     		 8.410,7.960,7.260,6.340,5.630,4.940,4.390,4.220,
     		 3.630,3.100,2.245,1.538,0.841,0.294,0.060};
  double kalephtotdataset1stat[]={1.120,0.910,0.900,0.920,0.880,0.900,0.860,0.880,
     		0.660,0.470,0.530,0.300,0.260,0.230,0.210,0.140,
    	        0.130,0.110,0.110,0.100,0.090,0.090,0.060,0.050,
     		0.029,0.025,0.013,0.005,0.002};
  double kaplehtotdataset1sys[]={0.010,0.010,0.010,0.020,0.020,0.020,0.050,0.120,
     		0.500,0.680,2.200,1.280,0.980,0.830,0.710,0.560,
     		0.470,0.370,0.320,0.280,0.240,0.220,0.180,0.150,
     		0.109,0.076,0.043,0.015,0.003 };


  int kdelphitotdata_np = 23;
  double kdelphitotdataset1z[]={0.0175,0.0225,0.0275,0.0325,0.0375,0.0450,0.0550,
               0.0650,0.0750,0.0850,0.0950,0.1100,0.1300,0.1500,
               0.1700,0.1900,0.2300,0.2800,0.3500,0.4500,0.5500,
               0.7000,0.9000};

  double kdelphitotdataset1sig[]={25.460,24.900,21.270,20.540,18.260,16.300,13.310,
               11.200,9.840,8.910,8.000,6.580,5.120,4.130,3.290,
               2.850,1.980,1.331,0.775,0.361,0.166,0.059,0.007
   };
  double kdelphitotdataset1stat[]={3.770,3.460,3.130,2.050,2.760,2.870,2.940,0.790,
               0.760,0.600,0.600,0.470,0.330,0.450,0.490,0.420,
       	       0.100,0.078,0.056,0.039,0.020,0.015,0.006};
  //double kdelphitotdataset1sys[]={
  //};

//---------- PIONS --------------

     int piudstotdata_np = 38;
  double piudstotdataset4z[]={0.0090,0.0110,0.0130,0.0150,0.0190,0.0245,0.0300,
               0.0355,0.0410,0.0465,0.0520,0.0575,0.0630,0.0685,
               0.0740,0.0795,0.0850,0.0935,0.1045,0.1155,0.1320,
               0.1535,0.1750,0.1970,0.2190,0.2410,0.2630,0.2850,
               0.3070,0.3345,0.3675,0.4005,0.4335,0.4660,0.5040,
               0.5480,0.6140,0.7130 }; 
  double piudstotdataset4sig[]={482.300,439.000,400.500,356.100,292.800,228.500,
                  176.600,144.400,121.700,102.500, 89.200, 75.300,
                  66.000, 57.810, 51.630, 45.950, 41.350, 35.240,
                  28.120, 23.570, 18.320, 13.220,  9.840,  7.470,
                   5.711,  4.414,  3.612,  2.886,  2.206,  1.739,
                   1.350,  0.874,  0.670,  0.520,  0.355,  0.248,
                   0.146,  0.071};
  double piudstotdataset4stat[]={2.300,2.300,2.000,1.900,1.000,1.000,0.900,0.800,
                  0.800,0.900,0.900,0.800,0.700,0.600,0.560,0.520,
                  0.490,0.320,0.290,0.270,0.170,0.140,0.110,0.090,
                  0.083,0.074,0.068,0.061,0.054,0.040,0.036,0.031,
                  0.029,0.026,0.021,0.018,0.012,0.009};
  double piudstotdataset4sys[]={7.200,3.700,3.300,3.000,2.400,1.900,1.400,1.200,
                  1.600,1.900,2.200,2.800,1.400,0.810,0.600,0.540,
                  0.490,0.420,0.350,0.300,0.240,0.190,0.150,0.110,
                  0.080,0.063,0.072,0.060,0.049,0.044,0.040,0.032,
                  0.029,0.025,0.020,0.016,0.012,0.007};

     int pialephtotdata_np = 39;
  double pialephtotdataset1z[]={0.0052,0.0057,0.0063,0.0067,0.0072,0.0077,0.0083,
                0.0088,0.0092,0.0097,0.0105,0.0115,0.0125,0.0135,
                0.0150,0.0170,0.0475,0.0525,0.0575,0.0625,0.0675,
                0.0725,0.0775,0.0825,0.0875,0.0950,0.1050,0.1150,
                0.1250,0.1350,0.1450,0.1550,0.1700,0.1900,0.2250,
                0.2750,0.3500,0.5000,0.7000};
  double pialephtotdataset1sig[]={482.900,462.600,496.500,511.200,507.700,538.500,
                  484.200,499.700,494.600,473.900,460.900,425.600,
                  420.700,380.500,360.800,324.000,103.960, 89.950,
                   78.960, 69.360, 61.350, 55.270, 49.910, 44.330,
                   40.240, 35.380, 29.510, 24.910, 21.060, 18.160,
                   15.460, 13.640, 11.000,  8.484,  5.621,  3.181,
                   1.563,  0.450,  0.077};
  double pialephtotdataset1stat[]={5.900,4.800,4.600,4.400,4.200,4.400,3.900,3.900,
                  3.800,3.600,2.500,2.300,2.300,2.200,1.500,1.400,
                  0.610,0.530,0.500,0.350,0.330,0.320,0.300,0.290,
                  0.270,0.180,0.170,0.160,0.140,0.130,0.120,0.120,
                  0.070,0.066,0.035,0.026,0.013,0.005,0.002 };
  double piaplehtotdataset1sys[]={1.300,0.900,0.800,0.800,0.700,0.700,0.600,0.700,
                 0.600,0.500,0.500,0.500,0.400,0.400,0.600,1.800,
                 2.090,1.020,0.900,0.720,0.600,0.490,0.440,0.380,
                 0.340,0.300,0.250,0.220,0.180,0.160,0.150,0.130,
                 0.110,0.094,0.071,0.047,0.028,0.010,0.002};
  int pidelphitotdata_np = 23;
  double pidelphitotdataset1z[]={0.0175,0.0225,0.0275,0.0325,0.0375,0.0450,0.0550,
                0.0650,0.0750,0.0850,0.0950,0.1100,0.1300,0.1500,
                0.1700,0.1900,0.2300,0.2800,0.3500,0.4500,0.5500,
                0.7000,0.9000};
  double pidelphitotdataset1sig[]={315.100,250.660,204.610,168.230,140.850,110.650,
                   83.880, 65.830, 52.530, 42.710, 35.190, 27.190,
                  19.640, 14.650, 11.240,  8.560,  5.450,  3.160,
                   1.670,  0.680,  0.307,  0.097,  0.0136}; //0.014}; Data was wrong accordind paper
  double pidelphitotdataset1stat[]={11.100,4.830,4.540,2.790,2.610,2.230,2.100,
                    1.330,1.230,1.040,0.890,0.740,0.590,0.480,
                    0.400,0.340,0.210,0.130,0.077,0.038,0.022,
                    0.014,0.0046}; // 0.005};


/*

     int kdata_np = 21;
  double ktotdataset1z[]={
   };
  double ktotdataset1sig[]={
   };
  double ktotdataset1stat[]={
   };
  double ktotdataset1sys[]={
  };
*/

 // TGraph* pkuds = new TGraph(kudsdata_np, kudsdataset1z, kudsdataset1sig);      
   TGraphErrors* pkuds = new TGraphErrors(kudsdata_np, kudsdataset1z, kudsdataset1sig,0,kudsdataset1stat);
   TGraphErrors* pkudstot = new TGraphErrors(kudstotdata_np, kudstotdataset4z, kudstotdataset4sig,0,kudstotdataset4stat);
   TGraphErrors* pkalephtot = new TGraphErrors(kalephtotdata_np, kalephtotdataset1z, kalephtotdataset1sig,0,kalephtotdataset1stat);
   TGraphErrors* pkdelphitot = new TGraphErrors(kdelphitotdata_np, kdelphitotdataset1z, kdelphitotdataset1sig,0,kdelphitotdataset1stat);
   TGraphErrors* pkpy6 = new TGraphErrors(26,kHx,kHy,kHex,kHe);

    TGraphErrors* ppiudstot = new TGraphErrors(piudstotdata_np, piudstotdataset4z, piudstotdataset4sig,0,piudstotdataset4stat);
   TGraphErrors* ppialephtot = new TGraphErrors(pialephtotdata_np, pialephtotdataset1z, pialephtotdataset1sig,0,pialephtotdataset1stat);
   TGraphErrors* ppidelphitot = new TGraphErrors(pidelphitotdata_np, pidelphitotdataset1z, pidelphitotdataset1sig,0,pidelphitotdataset1stat);
   TGraphErrors* ppipy6 = new TGraphErrors(26,piHx,piHy,piHex,piHe);
 
  TCanvas* canvaszk = new TCanvas("canvaszk","canvas zscal");
  canvaszk->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  pkuds->SetMarkerStyle(20);
  pkuds->SetMarkerColor(3);
  pkudstot->SetMarkerStyle(21);
  pkudstot->SetMarkerColor(1);
  pkalephtot->SetMarkerStyle(23);
  pkalephtot->SetMarkerColor(kBlue);
  pkdelphitot->SetMarkerStyle(24);
  pkdelphitot->SetMarkerColor(kGreen);
  pkpy6->SetMarkerStyle(20);
  pkpy6->SetMarkerColor(kRed);
     TLegend* lff = new TLegend(0.1,0.7,0.48,0.9);
   lff->SetHeader("e^{+}e^{-} #sqrt{s}=91 GeV");// Form("Pythia6, %s at #sqrt{s}=%d GeV",system,ECM),"C"); //"C" to center header
   lff->SetTextFont(42);
   lff->SetTextSize(0.042);
   lff->SetLineColorAlpha(kWhite,0);
   lff->SetLineWidth(0);
   lff->SetLineStyle(3);
   lff->SetShadowColor(kYellow);
   lff->SetFillColorAlpha(kWhite,0.05);
    lff->AddEntry(pkudstot,"UDS total","lep");
    lff->AddEntry(pkalephtot,"ALEPH total","lep"); 
    lff->AddEntry(pkdelphitot,"DELPHI total","lep");
    lff->AddEntry( pkpy6,"Pythia6","lep");
  pkudstot->SetName("K_UDS");
  pkdelphitot->SetTitle("fragmentation functions fo K");
  pkdelphitot->Draw("A P");
  pkudstot->Draw("same P");
  pkalephtot->Draw("same P");
 // pkdelphitot->Draw("same P");
  pkpy6->Draw("same P");
 // hzscaffid2[2]->Draw("same P");
  lff->Draw();

 
  TCanvas* canvaszpi = new TCanvas("canvaszpi","canvas zscal");
  canvaszpi->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  ppiudstot->SetMarkerStyle(21);
  ppiudstot->SetMarkerColor(1);
  ppialephtot->SetMarkerStyle(23);
  ppialephtot->SetMarkerColor(kBlue);
  ppidelphitot->SetMarkerStyle(24);
  ppidelphitot->SetMarkerColor(kGreen);
  ppipy6->SetMarkerStyle(20);
  ppipy6->SetMarkerColor(kRed);
     TLegend* lffpi = new TLegend(0.1,0.7,0.48,0.9);
   lffpi->SetHeader("e^{+}e^{-} #sqrt{s}=91 GeV");// Form("Pythia6, %s at #sqrt{s}=%d GeV",system,ECM),"C"); //"C" to center header
   lffpi->SetTextFont(42);
   lffpi->SetTextSize(0.042);
   lffpi->SetLineColorAlpha(kWhite,0);
   lffpi->SetLineWidth(0);
   lffpi->SetLineStyle(3);
   lffpi->SetShadowColor(kYellow);
   lffpi->SetFillColorAlpha(kWhite,0.05);
    lffpi->AddEntry(ppiudstot,"UDS total","lep");
    lffpi->AddEntry(ppialephtot,"ALEPH total","lep");
    lffpi->AddEntry(ppidelphitot,"DELPHI total","lep");
    lffpi->AddEntry( ppipy6,"Pythia6","lep");
  ppiudstot->SetName("pi_UDS");
  ppidelphitot->SetTitle("fragmentation functions fo #pi");
  ppidelphitot->Draw("AP");
    ppipy6->Draw("same P");
  ppiudstot->Draw("same P");
  ppialephtot->Draw("same P");
//  ppipy6->Draw("same P");
  lffpi->Draw();



};
