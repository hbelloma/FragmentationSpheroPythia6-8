void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05); 


void plotpointsnew(){

/*
  #ifndef __CINT__

#include "TFile.h"
#include "TError.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TDecayChannel.h"
#include <cstdlib>

using namespace std;
#endif
*/


  gROOT->SetStyle("Plain");


   
 //  TFile * finH=TFile::Open("ADistribfuncMBmill_tocompare_zdef2pqsobq2_new100m.root","READ");
//  TFile * finH=TFile::Open("ADistribfuncMBmill_tocompare_zdefEsobre91_new100m.root","READ");

// TFile * finH=TFile::Open("ADistribfuncMBmill_tocompare_zdefpsobre91_new100m.root","READ");
// TFile * finH=TFile::Open("ADistribfuncMBmill_eenewzbin_3MillTdef_GOOD.root","READ");       
 //TFile * finH=TFile::Open("ADistribfuncMBmill_eenewzbin_80mil_so_fordelphi.root","READ");  

// TFile * finH=TFile::Open("ADistribfuncMBmill_eenewzbin_test.root","READ");
// TFile * finH=TFile::Open("ADistribfuncMBmill_eenewzbin.root","READ");
  TFile * finFJ=TFile::Open("ADistribfuncMBmill_eenewzbin_80mil_so_fordelphi.root","READ"); 
  TFile * finSLD=TFile::Open("ADistribfuncMBmill_eenewzbin_SLD_SLAC.root","READ");
  TFile * finAL=TFile::Open("ADistribfuncMBmill_eenewzbin_ALEPH_LEP.root","READ");
  TFile * finH=TFile::Open("ADistribfuncMBmill_eenewzbin_DELPHI_LEP.root","READ");
  TFile *finSo=TFile::Open("ADistribfuncMBmill_eenewzbin_30MILL_goodSO.root","READ");

// TFile * finH=TFile::Open("ADistribfuncMBmill_eenewzbin_30MILL_goodSO.root","READ");

 // TFile * finH=TFile::Open("ADistribfuncMBmill_eenewzbin.root","READ");
 //  TFile * finH=TFile::Open("ADistribfuncMBmill_eenewzbin_3MillTdaf_good.root","READ"); 
 //   TFile * finH=TFile::Open("ADistribfuncMBmill_eenewzbin_1Mill_t350.root","READ"); 
 //  TFile * finH=TFile::Open("ADistribfuncMBmill_eenewzbin_py8.root","READ"); 
   //ADistribfuncMBmill_eenfjnewzbin.root","READ");
  //  TFile * finH = TFile::Open("ADistribfuncMBmill_eenfjrespnoSO.root","READ");
 // Not sense to compare but just to look difference in e+e- and pp
 // TFile * finH=TFile::Open("ADistribfuncMBmill_ppnewzbin.root","READ");

    const char* systema="e^{+}e^{-}";
    #define ECM  91 
    const Int_t nSobins = 10;
Double_t Sobins[nSobins+1] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
const Int_t nMultbins=13;
Double_t Multbins[nMultbins+1]={ 0.0, 1.0, 4.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 140.0};
const Int_t nPtBins      = 59;
Double_t xBins[nPtBins+1] = {
        0.01,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
        0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,
        1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,
        4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20
};
 const char *idt[10]={"Charged","#pi","K","p","K^{0}_{s}","#phi","#Lambda","#Xi","#Omega","rest"};
 const char *ptint[7]={"p^{ch jet}_{T} 5-50 GeV/c (inclusive)","p^{ch jet}_{T} 5-10 GeV/c","p^{ch jet}_{T} 10-15 GeV/c","p^{ch jet}_{T} 15-20 GeV/c","p^{ch jet}_{T} 20-30 GeV/c","p^{ch jet}_{T} 30-50 GeV/c","p^{ch jet}_{T} 40-100 GeV/c"};
 const char *Sopcin[nSobins+1]={"0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"};


    TH1F* hzscaffid[10];
    TH1F* hzscaffid2[10];
    TH1F* hzscaffidSLD[10];
    TH1F* hzscaffid2SLD[10];
    TH1F* hzscaffidAL[10];
    TH1F* hzscaffid2AL[10];
    TH1D* hevents=(TH1D*)finH->Get("hevents");
    TH1D* heventsAL=(TH1D*)finAL->Get("hevents");
    TH1D* heventsSLD=(TH1D*)finSLD->Get("hevents");

    TH1D* hjets=(TH1D*)finH->Get("hjets");



     TH1F* hzchjet[7];
     TH1F* hzscaffSo[nSobins];
     TH1F* hzscaffSoNbin[nSobins][nMultbins];
     TH1F* hzchjet2[7];
     TH1F* hzchjet2reb0;
     TH1F* hzchjet2reb02;
     TH1F* hzscaffSo2[nSobins];
     TH1F* hzscaffSoNbin2[nSobins][nMultbins];
     TH1F* hzjetRat[7];
     TH1F* hzSoRat[nSobins];
     TH1F* hzSoNbinRat[nSobins][nMultbins];


   for(Int_t i=0; i<7; i++){
    hzchjet[i]=(TH1F*)finFJ->Get(Form("hzchjet%d",i));
    hzchjet2[i]=(TH1F*)hzchjet[i]->Clone(Form("hzchjet2%d",i));
//    hzjetRat[i]=(TH1F*)hzchjet[i]->Clone(Form("hzjetRat%d",i));
   }

   hzchjet2reb0=(TH1F*)hzchjet2[0]->Rebin(2,"hzchjet2reb0");
   hzchjet2reb02=(TH1F*)hzchjet2reb0->Clone("hzchjet2reb02");

   for(Int_t i=0; i<nSobins; i++){
     hzscaffSo[i]=(TH1F*)finSo->Get(Form("hzscaffSo%d",i));
     hzscaffSo2[i]=(TH1F*)hzscaffSo[i]->Clone(Form("hzscaffSo2%d",i));
//     hzSoRat[i]=(TH1F*)hzscaffSo[i]->Clone(Form("hzSoRat%d",i));
     for(Int_t j=0; j<nMultbins; j++){
        hzscaffSoNbin[i][j]=(TH1F*)finSo->Get(Form("hzscaffSoNbin%d_%d",i,j));
        hzscaffSoNbin2[i][j]=(TH1F*)hzscaffSoNbin[i][j]->Clone(Form("hzscaffSoNbin2%d_%d",i,j));
        hzSoNbinRat[i][j]=(TH1F*)hzscaffSoNbin[i][j]->Clone(Form("hzSoNbinRat%d_%d",i,j));
     }
    }


 
    TH1D* heventsid=(TH1D*)finH->Get("heventsid");
    TH1D* hjetsptcut=(TH1D*)finFJ->Get("hjetsptcut");
    TH1D* heventsSo=(TH1D*)finSo->Get("heventsSo");
    TH1D* heventsSoN[nSobins];
    for(Int_t i=0; i<nSobins; i++){
      heventsSoN[i]=(TH1D*)finSo->Get(Form("heventsSoN%d",i));
    }



  for(Int_t i=0; i<10;i++){
    hzscaffid[i]=(TH1F*)finH->Get(Form("hzscaffid%d",i));
    hzscaffid2[i]=(TH1F*)hzscaffid[i]->Clone(Form("hzscaffid%d",i));
    hzscaffidSLD[i]=(TH1F*)finSLD->Get(Form("hzscaffid%d",i));
    hzscaffid2SLD[i]=(TH1F*)hzscaffidSLD[i]->Clone(Form("hzscaffidSLD%d",i));
    hzscaffidAL[i]=(TH1F*)finAL->Get(Form("hzscaffid%d",i));
    hzscaffid2AL[i]=(TH1F*)hzscaffidAL[i]->Clone(Form("hzscaffidAL%d",i));

  }
    Int_t nevents=hevents->GetBinContent(10);
    Int_t neventsAL=heventsAL->GetBinContent(10);
    Int_t neventsSLD=heventsSLD->GetBinContent(10);
 //  Int_t neventsonevalid=hevents->GetBinContent(2);
    Int_t neventswpi=nevents;//heventsid->GetBinContent(2);
    Int_t neventswk=nevents;//heventsid->GetBinContent(3);
    Int_t njets=hjets->GetBinContent(1);
    Int_t njetspt[7];
    for(Int_t i=1; i<8; i++){
          njetspt[i]=hjetsptcut->GetBinContent(i);
    }
    //cout<<"number of events="<<nevents<<endl; 
     Int_t neventsSo[nSobins];
     Int_t neventsSoNch[nSobins][nMultbins];
    for( Int_t i=0; i<nSobins; i++){
      neventsSo[i]=heventsSo->GetBinContent(i+1);
      for(Int_t j=0; j<nMultbins; j++){
       neventsSoNch[i][j]=heventsSoN[i]->GetBinContent(j+1);
      }
    }
    

    Int_t nbinzjr0= hzchjet2reb0->GetXaxis()->GetNbins();
    Double_t contzjr0=0;    
    Double_t errzjr0=0;
    for(Int_t i=0; i<nbinzjr0+1; i++){
      Double_t zjr0= hzchjet2reb0->GetBinCenter(i);
      Double_t zjwidthr0= hzchjet2reb0->GetBinWidth(i);
      Double_t contzjr0= hzchjet2reb0->GetBinContent(i);
      Double_t errzjr0= hzchjet2reb0->GetBinError(i);
       hzchjet2reb02->SetBinContent(i,contzjr0/zjwidthr0);
       hzchjet2reb02->SetBinError(i, errzjr0/zjwidthr0);
    } 
      hzchjet2reb02->Scale(1./(njetspt[1]));


    for(Int_t piden=0; piden<7; piden++){
    Int_t nbinzj=hzchjet[0]->GetXaxis()->GetNbins();
    Int_t nbinzj2=hzchjet2[0]->GetXaxis()->GetNbins();
    Double_t contzj=0;
    Double_t errzj=0;
    for(Int_t i=0; i<nbinzj+1; i++){
      Double_t zj=hzchjet[piden]->GetBinCenter(i);
      Double_t zjwidth=hzchjet[piden]->GetBinWidth(i);
      for(Int_t j=0; j<nbinzj2+1; j++){
      Double_t zj2=hzchjet2[piden]->GetBinCenter(j);
      Double_t zjwidth2=hzchjet2[piden]->GetBinWidth(j);
       if(i==j){
       Double_t contzj=hzchjet[piden]->GetBinContent(i);
       Double_t errzj=hzchjet[piden]->GetBinError(i);
      // cout<<"cont"<<contzjid<<endl;
       hzchjet2[piden]->SetBinContent(j,contzj/zjwidth2);
       hzchjet2[piden]->SetBinError(j,errzj/zjwidth2);
       }
      }
     }

     hzchjet2[piden]->Scale(1./(njetspt[piden+1]));
     hzjetRat[piden]=(TH1F*)hzchjet2[piden]->Clone(Form("hzjetRat%d",piden));
     //hzjetRat[piden]->Scale(1./(njetspt[piden+1]));
     
     // TO CHECK IN THE OTHER CODE, REBINING FOR JETS IN LARGE Z INTERVALS
     if(piden<4)
     hzjetRat[piden]->Divide(hzchjet2[0]);
     else hzjetRat[piden]->Divide(hzchjet2reb02); 
   } //end loop of pt range for zjet

 

    for(Int_t so=0; so<nSobins; so++){
    Int_t nbinzscaSo=hzscaffSo[0]->GetXaxis()->GetNbins();
    Int_t nbinzscaSo2=hzscaffSo2[0]->GetXaxis()->GetNbins();
    Double_t contzscaSo=0;
    Double_t errzscaSo=0;
    for(Int_t i=0; i<nbinzscaSo+1; i++){
      Double_t zscaSo=hzscaffSo[so]->GetBinCenter(i);
      Double_t zscawidthSo=hzscaffSo[so]->GetBinWidth(i);
      for(Int_t j=0; j<nbinzscaSo2+1; j++){
      Double_t zscaSo2=hzscaffSo2[so]->GetBinCenter(j);
      Double_t zscawidthSo2=hzscaffSo2[so]->GetBinWidth(j);
       if(i==j){
       Double_t contzscaSo=hzscaffSo[so]->GetBinContent(i);
       Double_t errzscaSo=hzscaffSo[so]->GetBinError(i);
      // cout<<"cont"<<contzscaid<<endl;
       hzscaffSo2[so]->SetBinContent(j,contzscaSo/zscawidthSo2);
       hzscaffSo2[so]->SetBinError(j,errzscaSo/zscawidthSo2);
       }
      }
     }
        hzscaffSo2[so]->Scale(1./(neventsSo[so]));
        hzSoRat[so]=(TH1F*)hzscaffSo2[so]->Clone(Form("hzSoRat%d",so));
        //cout<<"so="<<so<<endl; 
        hzSoRat[so]->Divide(hzscaffSo2[0]);
    } //end loop of So


      for(Int_t so=0; so<nSobins; so++){
    for(Int_t n=0; n<nMultbins; n++){
    Int_t nbinzscaSoNbin=hzscaffSoNbin[0][0]->GetXaxis()->GetNbins();
    Int_t nbinzscaSoNbin2=hzscaffSoNbin2[0][0]->GetXaxis()->GetNbins();
    Double_t contzscaSoNbin=0;
    Double_t errzscaSoNbin=0;
    for(Int_t i=0; i<nbinzscaSoNbin+1; i++){
      Double_t zscaSoNbin=hzscaffSoNbin[so][n]->GetBinCenter(i);
      Double_t zscawidthSoNbin=hzscaffSoNbin[so][n]->GetBinWidth(i);
      for(Int_t j=0; j<nbinzscaSoNbin2+1; j++){
      Double_t zscaSoNbin2=hzscaffSoNbin2[so][n]->GetBinCenter(j);
      Double_t zscawidthSoNbin2=hzscaffSoNbin2[so][n]->GetBinWidth(j);
       if(i==j){
       Double_t contzscaSoNbin=hzscaffSoNbin[so][n]->GetBinContent(i);
       Double_t errzscaSoNbin=hzscaffSoNbin[so][n]->GetBinError(i);
      // cout<<"cont"<<contzscaid<<endl;
       hzscaffSoNbin2[so][n]->SetBinContent(j,contzscaSoNbin/zscawidthSoNbin2);
       hzscaffSoNbin2[so][n]->SetBinError(j,errzscaSoNbin/zscawidthSoNbin2);
       }
      }
     }
    } //end loop of mult
     // hzSoNbinfor[so]=;
     // hzSoNbinRat[so][0]->Divide() ->Divide(hzscaffSo2[so][12]);
    } // end loop of So




   cout<<"neventos="<<nevents<<endl;
  Double_t kHx[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHy[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHe[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHex[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHx[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHy[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHe[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHex[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

 // KSLD=32, PISLD=38, KAL=29, PIAL=39
  Double_t kHxSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHySLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHeSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHexSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHxSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHySLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHeSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHexSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  Double_t kHxAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHyAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHeAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHexAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHxAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHyAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHeAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHexAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  Double_t kHyr[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHer[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHexr[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHyr[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHer[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHexr[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  Double_t kHyrSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHerSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHexrSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHyrSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHerSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHexrSLD[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  Double_t kHyrAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHerAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHexrAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHyrAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHerAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHexrAL[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

/*
*/
/*
  Double_t kHx[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHy[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t kHe[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHx[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHy[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t piHe[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
*/ 
 

  for(Int_t piden=1; piden<3; piden++){
    Int_t nbinzscaid=hzscaffid[1]->GetXaxis()->GetNbins();
    cout<<"nbins= "<<nbinzscaid<<endl;
    Int_t nbinzscaid2=hzscaffid2[1]->GetXaxis()->GetNbins();
    Double_t contzscaid=0;
    Double_t errzscaid=0;
    for(Int_t i=1; i<nbinzscaid+1; i++){   //1- +3
     // cout<<"----i="<<i<<endl;
      Double_t zscaid=hzscaffid[piden]->GetBinCenter(i);
      if(piden==1)piHx[i]=zscaid;
      if(piden==2)kHx[i]=zscaid;
      Double_t zscawidthid=hzscaffid[piden]->GetBinWidth(i);
      for(Int_t j=1; j<nbinzscaid2+1; j++){
       Double_t zscaid2=hzscaffid2[piden]->GetBinCenter(j);
       Double_t zscawidthid2=hzscaffid2[piden]->GetBinWidth(j);
       if(i==j){
       //  cout<<"----j="<<j<<endl;
       Double_t contzscaid=hzscaffid[piden]->GetBinContent(i);
       Double_t errzscaid=hzscaffid[piden]->GetBinError(i);
      //cout<<"cont= "<<contzscaid<<", zscawidthid2="<<zscawidthid2<<endl;
       hzscaffid2[piden]->SetBinContent(j,contzscaid/zscawidthid2);
       hzscaffid2[piden]->SetBinError(j,errzscaid/zscawidthid2);
       Double_t ratioc=contzscaid/(zscawidthid2); //*nevents); // /(zscawidthid2*80000); ///(zscawidthid2*nevents); //zscawidthid2
       Double_t ratioe=errzscaid/(zscawidthid2); //*nevents);  // /(zscawidthid2*80000); /// (zscawidthid2*nevents);
       // cout<<"-------------------------err="<<errzscaid<<", Perc%="<<errzscaid/contzscaid<<endl;
       // cout<<"ratioc="<<ratioc<<", ratioe="<<ratioe<<endl;
       if(piden==2){
       kHy[i]=ratioc/neventswk; // /80000; // /(80000*0.05); // /55960;
       kHe[i]=ratioe/neventswk; // /80000; // /(80000*0.05); // /55960;
       kHex[i]=zscawidthid2/2;
       //cout<<"kHx="<<kHx[i]<<", kHy="<<kHy[i]<<", kHe="<<kHe[i]<<endl;
       }
       if(piden==1){
       piHy[i]=ratioc/neventswpi; // /80000; // /(80000*0.05); // /35380;
       piHe[i]=ratioe/neventswpi; // /80000; // /(80000*0.05); // /35380;
       piHex[i]=zscawidthid2/2;
       //cout<<"piHx="<<piHx[i]<<", piHy="<<piHy[i]<<", piHe="<<piHe[i]<<endl;
       }

       }
      }

     }
      

      Int_t nbinzscaidAL=hzscaffidAL[piden]->GetXaxis()->GetNbins();
     cout<<"nbinsAL= "<<nbinzscaidAL<<endl;
    Double_t contzscaidAL=0;
    Double_t errzscaidAL=0;
    for(Int_t i=1; i<nbinzscaidAL+1; i++){ 
      cout<<"----i="<<i<<endl;
      Double_t zscaidAL=hzscaffidAL[piden]->GetBinCenter(i);
      Double_t zscawidthidAL=hzscaffidAL[piden]->GetBinWidth(i);
       Double_t contzscaidAL=hzscaffidAL[piden]->GetBinContent(i);
         cout<<"contAL= "<<contzscaidAL<<", zscawidthidAL="<<zscawidthidAL<<endl;
       Double_t errzscaidAL=hzscaffidAL[piden]->GetBinError(i);
       Double_t ratiocAL=contzscaidAL/(zscawidthidAL);
       Double_t ratioeAL=errzscaidAL/(zscawidthidAL);
        cout<<"ratiocAL="<<ratiocAL<<", ratioeAL="<<ratioeAL<<endl;
       if(piden==2){
        if(i==9 || i==13){ 
         kHxAL[i]=zscaidAL;
         kHyAL[i]=0.0;
         kHeAL[i]=0.0;
         kHexAL[i]=0.0;
	}else{     
        kHxAL[i]=zscaidAL;
        kHyAL[i]=100*ratiocAL/neventsAL; 
        kHeAL[i]=ratioeAL/neventsAL;
        kHexAL[i]=zscawidthidAL/2;
       }
       cout<<"kHxAL="<<kHxAL[i]<<", kHyAL="<<kHyAL[i]<<", kHeAL="<<kHeAL[i]<<endl;
       }
       if(piden==1){
       if(i==17){
        piHxAL[i]=zscaidAL;
        piHyAL[i]=0.0;
        piHeAL[i]=0.0;
        piHexAL[i]=0.0;
       }else{
       piHxAL[i]=zscaidAL;
       piHyAL[i]=100*ratiocAL/neventsAL; 
       piHeAL[i]=ratioeAL/neventsAL; 
       piHexAL[i]=zscawidthidAL/2;
       }
       cout<<"piHxAL="<<piHxAL[i]<<", piHyAL="<<piHyAL[i]<<", piHeAL="<<piHeAL[i] <<", piHexAL="<<piHexAL[i]<<endl;
       }
     }


   
    Int_t nbinzscaidSLD=hzscaffidSLD[piden]->GetXaxis()->GetNbins();
    Double_t contzscaidSLD=0;
    Double_t errzscaidSLD=0;
     cout<<"nbinsSLD= "<<nbinzscaidSLD<<endl;
    for(Int_t i=1; i<nbinzscaidSLD+1; i++){
       cout<<"----i="<<i<<endl;
      Double_t zscaidSLD=hzscaffidSLD[piden]->GetBinCenter(i);
      Double_t zscawidthidSLD=hzscaffidSLD[piden]->GetBinWidth(i);
       Double_t contzscaidSLD=hzscaffidSLD[piden]->GetBinContent(i);
        cout<<"contSLD= "<<contzscaidSLD<<", zscawidthidSLD="<<zscawidthidSLD<<endl;
       Double_t errzscaidSLD=hzscaffidSLD[piden]->GetBinError(i);
       Double_t ratiocSLD=contzscaidSLD/(zscawidthidSLD);
       Double_t ratioeSLD=errzscaidSLD/(zscawidthidSLD);
         cout<<"ratiocSLD="<<ratiocSLD<<", ratioeSLD="<<ratioeSLD<<endl;
       if(piden==2){
       if(i==19){ 
         kHxSLD[i]=zscaidSLD;
         kHySLD[i]=0.0;
         kHeSLD[i]=0.0;
         kHexSLD[i]=0.0;
        }else{
        kHxSLD[i]=zscaidSLD;
        kHySLD[i]=1000*ratiocSLD/neventsSLD;
        kHeSLD[i]=ratioeSLD/neventsSLD;
        kHexSLD[i]=zscawidthidSLD/2;
        cout<<"kHxSLD="<<kHxSLD[i]<<", kHySLD="<<kHySLD[i]<<", kHeSLD="<<kHeSLD[i]<<endl;
        }
       }
       if(piden==1){
       piHxSLD[i]=zscaidSLD;
       piHySLD[i]=1000*ratiocSLD/neventsSLD;
       piHeSLD[i]=ratioeSLD/neventsSLD;
       piHexSLD[i]=zscawidthidSLD/2;
         cout<<"piHxSLD="<<piHxSLD[i]<<", piHySLD="<<piHySLD[i]<<", piHeSLD="<<piHeSLD[i]<<endl;
       }
     }


    } //end loop of pid

 /*
  for(Int_t i=0; i<nbinzscaid+1; i++){
    cout<<"x="<<kHx[i]<<"y="<<kHy[i]<<endl;
  }
 */


 // SLD (Sqrt(s)=91.2), PRD59(1999)052001, 1/SIGMAH DS/DXP
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

   // SLD (Sqrt(s)=91.2), PRD59(1999)052001, 1/SIGMAH DS/DXP
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
    for(Int_t i=0; i<kudstotdata_np; i++){
      kudstotdataset4sig[i]=1000*kudstotdataset4sig[i];
   }

  double kudstotdataset4stat[]={0.470,0.450,0.430,0.430,0.450,0.490,0.530,0.540,
                   0.490,0.460,0.440,0.440,0.440,0.310,0.340,0.370,
                   0.290,0.400,0.080,0.065,0.057,0.050,0.045,0.034,
                   0.029,0.026,0.026,0.023,0.019,0.016,0.011,0.009};
  double kudstotdataset4sys[]={0.530,0.620,0.640,0.800,0.940,1.030,1.140,2.110,
                  1.270,0.870,0.790,0.760,0.720,0.680,0.630,0.600,
                  0.610,0.800,0.300,0.057,0.042,0.034,0.030,0.027,
                  0.023,0.023,0.023,0.020,0.018,0.014,0.011,0.006};

 //ALEPH (Sqrt(s)=91.2), ZPhys C66(1995) 355, 1/SIGMAH DS/DXP
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
    for(Int_t i=0; i<kalephtotdata_np; i++){
      kalephtotdataset1sig[i]=100*kalephtotdataset1sig[i];
    }

  double kalephtotdataset1stat[]={1.120,0.910,0.900,0.920,0.880,0.900,0.860,0.880,
     		0.660,0.470,0.530,0.300,0.260,0.230,0.210,0.140,
    	        0.130,0.110,0.110,0.100,0.090,0.090,0.060,0.050,
     		0.029,0.025,0.013,0.005,0.002};
  double kaplehtotdataset1sys[]={0.010,0.010,0.010,0.020,0.020,0.020,0.050,0.120,
     		0.500,0.680,2.200,1.280,0.980,0.830,0.710,0.560,
     		0.470,0.370,0.320,0.280,0.240,0.220,0.180,0.150,
     		0.109,0.076,0.043,0.015,0.003 };

//DELPHI (Sqrt(s)=91.2), EPJ C5,585 (1998) 1/SIGMAH DS/DXP 
  int kdelphitotdata_np = 23;
  double kdelphitotdataset1z[]={0.0175,0.0225,0.0275,0.0325,0.0375,0.0450,0.0550,
               0.0650,0.0750,0.0850,0.0950,0.1100,0.1300,0.1500,
               0.1700,0.1900,0.2300,0.2800,0.3500,0.4500,0.5500,
               0.7000,0.9000};

  double kdelphitotdataset1sig[]={25.460,24.900,21.270,20.540,18.260,16.300,13.310,
               11.200,9.840,8.910,8.000,6.580,5.120,4.130,3.290,
               2.850,1.980,1.331,0.775,0.361,0.166,0.059,0.007
   };

  for(Int_t i=0; i<kdelphitotdata_np; i++){
      kdelphitotdataset1sig[i]=10*kdelphitotdataset1sig[i];
   }

  double kdelphitotdataset1stat[]={3.770,3.460,3.130,2.050,2.760,2.870,2.940,0.790,
               0.760,0.600,0.600,0.470,0.330,0.450,0.490,0.420,
       	       0.100,0.078,0.056,0.039,0.020,0.015,0.006};
  //double kdelphitotdataset1sys[]={
  //};

//---------- PIONS --------------
  //base for all plots
//   int pibase_np=40;
//   double pibasez[]={0.004, 0.0090,0.0110,0.0130,0.0150,0.0190,0.0245,0.0300,
   int pibase_np=38;
   double pibasez[]={0.0110,0.0130,0.0150,0.0190,0.0245,0.0300,
               0.0355,0.0410,0.0465,0.0520,0.0575,0.0630,0.0685,
               0.0740,0.0795,0.0850,0.0935,0.1045,0.1155,0.1320,
               0.1535,0.1750,0.1970,0.2190,0.2410,0.2630,0.2850,
               0.3070,0.3345,0.3675,0.4005,0.4335,0.4660,0.5040,
               0.5480,0.6140,0.7130,1.01 };

 //  double pibasezr[]={0.004, 0.0090,0.0110,0.0130,0.0150,0.0190,0.0245,0.0300,
   double pibasezr[]={0.0110,0.0130,0.0150,0.0190,0.0245,0.0300,
               0.0355,0.0410,0.0465,0.0520,0.0575,0.0630,0.0685,
               0.0740,0.0795,0.0850,0.0935,0.1045,0.1155,0.1320,
               0.1535,0.1750,0.1970,0.2190,0.2410,0.2630,0.2850,
               0.3070,0.3345,0.3675,0.4005,0.4335,0.4660,0.5040,
               0.5480,0.6140,0.7130,1.01 };
 
//  double pibasesig[]={1000000, 482.300,439.000,400.500,356.100,292.800,228.500,
  double pibasesig[]={1000000,400.500,356.100,292.800,228.500,
                  176.600,144.400,121.700,102.500, 89.200, 75.300,
                  66.000, 57.810, 51.630, 45.950, 41.350, 35.240,
                  28.120, 23.570, 18.320, 13.220,  9.840,  7.470,
                   5.711,  4.414,  3.612,  2.886,  2.206,  1.739,
                   1.350,  0.874,  0.670,  0.520,  0.355,  0.248,
                   0.146,  0.071, 0.0001};

//  double pibasesigr[]={2,2,2,2,2,2,2,
  double pibasesigr[]={ 2,2,2,2,2,
                        2,2,2,2,2,2,
			2,2,2,2,2,2,
			2,2,2,2,2,2,
			2,2,2,2,2,2,
                        2,2,2,2,2,2,
    			2,2,0};

//SLD (Sqrt(s)=91.2), PRD59(1999)052001, 1/SIGMAH DS/DXP norm 1.7%
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
  
  for(Int_t i=0; i<piudstotdata_np; i++){
      piudstotdataset4sig[i]=1000*piudstotdataset4sig[i];
   }

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

//ALEPH (Sqrt(s)=91.2), ZPhys C66(1995) 355, 1/SIGMAH DS/DXP norm 3%
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

   for(Int_t i=0; i<pialephtotdata_np; i++){
      pialephtotdataset1sig[i]=100*pialephtotdataset1sig[i];
    }

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


//DELPHI (Sqrt(s)=91.2), EPJ C5,585 (1998) 1/SIGMAH DS/DXP  no norm
  int pidelphitotdata_np = 23;
  double pidelphitotdataset1z[]={0.0175,0.0225,0.0275,0.0325,0.0375,0.0450,0.0550,
                0.0650,0.0750,0.0850,0.0950,0.1100,0.1300,0.1500,
                0.1700,0.1900,0.2300,0.2800,0.3500,0.4500,0.5500,
                0.7000,0.9000};
  double pidelphitotdataset1sig[]={315.100,250.660,204.610,168.230,140.850,110.650,
                   83.880, 65.830, 52.530, 42.710, 35.190, 27.190,
                  19.640, 14.650, 11.240,  8.560,  5.450,  3.160,
                   1.670,  0.680,  0.307,  0.097,  0.0136}; //0.014}; Data was wrong accordind paper
 
   for(Int_t i=0; i<pidelphitotdata_np; i++){
      pidelphitotdataset1sig[i]=10*pidelphitotdataset1sig[i];
   }


  double pidelphitotdataset1stat[]={11.100,4.830,4.540,2.790,2.610,2.230,2.100,
                    1.330,1.230,1.040,0.890,0.740,0.590,0.480,
                    0.400,0.340,0.210,0.130,0.077,0.038,0.022,
                    0.014,0.0046}; // 0.005};
  // ratio to DELPHI
   for(Int_t i=0; i<23; i++){
     cout<<"****** DELPHI ratio ******"<<endl;
     cout<<"bin i= "<< i <<endl;
     cout<<"piHy[i]="<<piHy[i]<<", pidelphitotdataset1sig[i]="<<pidelphitotdataset1sig[i] <<endl;
     cout<<"piHex[i]="<<piHex[i]<<", pidelphitotwidth[i]="<<hzscaffid[1]->GetBinWidth(i) <<endl;
     cout<<"piHx[i]"<< piHx[i] <<"GetBinCenter(i)=" << hzscaffid[1]->GetBinCenter(i) <<endl;

     Double_t percerrdel=pidelphitotdataset1stat[i]/pidelphitotdataset1z[i];
     cout<<"piHy[i+2]="<<piHy[i+2]<<", pidelphitotdataset1sig[i]="<<pidelphitotdataset1sig[i]*0.1 <<endl;
     cout<<"piHex[i+2]="<<piHex[i+2]<<", pidelphitotwidth[i]="<<hzscaffid[1]->GetBinWidth(i) <<endl;
     cout<<"piHx[i+2]"<< piHx[i+2] <<", GetBinCenter(i)=" << hzscaffid[1]->GetBinCenter(i) <<endl;

     kHyr[i]=kHy[i+2]/(0.1*kdelphitotdataset1sig[i]);
     kHer[i]=kHyr[i]*TMath::Sqrt( ((kHe[i+2]/ kHy[i+2])*(kHe[i+2]/ kHy[i+2]))+((kdelphitotdataset1stat[i]/kdelphitotdataset1sig[i]))*(kdelphitotdataset1stat[i]/kdelphitotdataset1sig[i]));            
     kHexr[i]=kHex[i+2];
     
     piHyr[i]=piHy[i+2]/(0.1*pidelphitotdataset1sig[i]);
     piHer[i]=piHyr[i]*TMath::Sqrt( ((piHe[i+2]/ piHy[i+2])*(piHe[i+2]/ piHy[i+2]))+((pidelphitotdataset1stat[i]/pidelphitotdataset1sig[i]))*(pidelphitotdataset1stat[i]/pidelphitotdataset1sig[i]));        
     piHexr[i]=piHex[i+2];
    
    // cout<<"piHyr="<<piHyr[i]<<", piHer="<<piHer[i]<<", kHyr[i]="<<kHyr[i]<<", kHer=" <<kHer[i]<<endl;
    //  cout<<"perceerrdel="<<percerrdel<<endl;

   }

  // RATIO TO ALEPH
   for(Int_t i=1; i<pialephtotdata_np; i++){
      if(i>15){
         piHyrAL[i]=piHyAL[i+2]/(pialephtotdataset1sig[i]);
         piHerAL[i]=piHyrAL[i]*TMath::Sqrt( ((piHeAL[i+2]/ piHyAL[i+2])*(piHeAL[i+2]/ piHyAL[i+2]))+((pialephtotdataset1stat[i]/pialephtotdataset1sig[i]))*(pialephtotdataset1stat[i]/pialephtotdataset1sig[i]));
         piHexrAL[i]=piHexAL[i+2];
      }
      else{      
         piHyrAL[i]=piHyAL[i+1]/(pialephtotdataset1sig[i]);
         piHerAL[i]=piHyrAL[i]*TMath::Sqrt( ((piHeAL[i+1]/ piHyAL[i+1])*(piHeAL[i+1]/ piHyAL[i+1]))+((pialephtotdataset1stat[i]/pialephtotdataset1sig[i]))*(pialephtotdataset1stat[i]/pialephtotdataset1sig[i]));
         piHexrAL[i]=piHexAL[i+1];
      }

     cout<<"***** ALEPH ratio *******"<<endl;
     cout<<"bin i= "<< i <<endl;
     cout<<"piHyAL[i+2]="<<piHyAL[i+2]<<", pialephtotdataset1sig[i]="<<pialephtotdataset1sig[i] <<endl;
     cout<<"piHexAL[i+2]="<<piHexAL[i+2]<<", pialephtotwidth[i+2]="<<hzscaffidAL[1]->GetBinWidth(i+2) <<endl;
     cout<<"piHxAL[i]"<< piHxAL[i] <<"GetBinCenter(i)=" << hzscaffidAL[1]->GetBinCenter(i) <<endl;
      cout<<"i=" << i<< "piHeAL[i]="<< piHeAL[i] <<", piHerAL[i]=" << piHerAL[i] <<endl;
    }
    for(Int_t i=1; i<kalephtotdata_np; i++){
      if(i<8){
        kHyrAL[i]=kHyAL[i+1]/(kalephtotdataset1sig[i]);
        kHerAL[i]=kHyrAL[i]*TMath::Sqrt( ((kHeAL[i+1]/ kHyAL[i+1])*(kHeAL[i+1]/ kHyAL[i+1]))+((kalephtotdataset1stat[i]/kalephtotdataset1sig[i]))*(kalephtotdataset1stat[i]/kalephtotdataset1sig[i]));
        kHexrAL[i]=kHexAL[i+1];
      }
      else if(i>7 && i<11){
        kHyrAL[i]=kHyAL[i+2]/(kalephtotdataset1sig[i]);
        kHerAL[i]=kHyrAL[i]*TMath::Sqrt( ((kHeAL[i+2]/ kHyAL[i+2])*(kHeAL[i+2]/ kHyAL[i+2]))+((kalephtotdataset1stat[i]/kalephtotdataset1sig[i]))*(kalephtotdataset1stat[i]/kalephtotdataset1sig[i]));
        kHexrAL[i]=kHexAL[i+2];
      }
      if(i>10){
        kHyrAL[i]=kHyAL[i+3]/(kalephtotdataset1sig[i]);
        kHerAL[i]=kHyrAL[i]*TMath::Sqrt( ((kHeAL[i+3]/ kHyAL[i+3])*(kHeAL[i+3]/ kHyAL[i+3]))+((kalephtotdataset1stat[i]/kalephtotdataset1sig[i]))*(kalephtotdataset1stat[i]/kalephtotdataset1sig[i]));
        kHexrAL[i]=kHexAL[i+3];
      }
     
     cout<<"bin i= "<< i <<endl;
     cout<<"kHyAL[i+3]="<<kHyAL[i+3]<<", kalephtotdataset1sig[i]="<<kalephtotdataset1sig[i] <<endl;
     cout<<"kHexAL[i+3]="<<kHexAL[i+3]<<", kalephtotwidth[i+3]="<<hzscaffidAL[2]->GetBinWidth(i+2) <<endl;
     cout<<"kHxAL[i]"<< kHxAL[i] <<"GetBinCenter(i)=" << hzscaffidAL[2]->GetBinCenter(i) <<endl;
     cout<<"i=" << i<< "kHeAL[i]="<< kHeAL[i] <<", kHerAL[i]=" << kHerAL[i] <<endl;
    }


  // RATIO TO SLD 
   for(Int_t i=1; i<piudstotdata_np; i++){
     piHyrSLD[i]=piHySLD[i+1]/(piudstotdataset4sig[i]);
     piHerSLD[i]=piHyrSLD[i]*TMath::Sqrt( ((piHeSLD[i+1]/ piHySLD[i+1])*(piHeSLD[i+1]/ piHySLD[i+1]))+((piudstotdataset4stat[i]/piudstotdataset4sig[i]))*(piudstotdataset4stat[i]/piudstotdataset4sig[i])); 
     piHexrSLD[i]=piHexSLD[i+1]; 
     cout<<"***** SLD ratio *******"<<endl;
     cout<<"bin i= "<< i <<endl;
     cout<<"piHySLD[i+1]="<<piHySLD[i+1]<<", piudstotdataset4sig[i]="<<piudstotdataset4sig[i] <<endl;
     cout<<"piHexSLD[i+1]="<<piHexSLD[i+1]<<", piudstotwidth[i+1]="<<hzscaffid[1]->GetBinWidth(i+1) <<endl;
     cout<<"piHxSLD[i]"<< piHxSLD[i] <<"GetBinCenter(i)=" << hzscaffid[1]->GetBinCenter(i) <<endl;
     cout<<"i=" << i<< "piHeSLD[i]="<< piHeSLD[i] <<", piHerSLD[i]=" << piHerSLD[i] <<endl;

    }
    for(Int_t i=1; i<kudstotdata_np; i++){
      if(i>17){
        kHyrSLD[i]=kHySLD[i+2]/(kudstotdataset4sig[i]);
        kHerSLD[i]=kHyrSLD[i]*TMath::Sqrt( ((kHeSLD[i+2]/ kHySLD[i+2])*(kHeSLD[i+2]/ kHySLD[i+2]))+((kudstotdataset4stat[i]/kudstotdataset4sig[i]))*(kudstotdataset4stat[i]/kudstotdataset4sig[i]));
        kHexrSLD[i]=kHexSLD[i+2];
      }
      else{
        kHyrSLD[i]=kHySLD[i+1]/(kudstotdataset4sig[i]);
        kHerSLD[i]=kHyrSLD[i]*TMath::Sqrt( ((kHeSLD[i+1]/ kHySLD[i+1])*(kHeSLD[i+1]/ kHySLD[i+1]))+((kudstotdataset4stat[i]/kudstotdataset4sig[i]))*(kudstotdataset4stat[i]/kudstotdataset4sig[i]));  
        kHexrSLD[i]=kHexSLD[i+1];
      }
     cout<<"***** ALEPH ratio *******"<<endl;
     cout<<"bin i= "<< i <<endl;
     cout<<"kHySLD[i+1]="<<kHySLD[i+1]<<", kudstotdataset4sig[i]="<<kudstotdataset4sig[i] <<endl;
     cout<<"kHexSLD[i+1]="<<kHexSLD[i+1]<<", kudstotwidth[i+1]="<<hzscaffidSLD[1]->GetBinWidth(i+1) <<endl;
     cout<<"kHxSLD[i]"<< kHxSLD[i] <<"GetBinCenter(i)=" << hzscaffid[1]->GetBinCenter(i) <<endl;
      cout<<"i=" << i<< "kHeSLD[i]="<< kHeSLD[i] <<", kHerSLD[i]=" << kHerSLD[i] <<endl;


    }
 


 // TGraph* pkuds = new TGraph(kudsdata_np, kudsdataset1z, kudsdataset1sig);      
   TGraphErrors* pkuds = new TGraphErrors(kudsdata_np, kudsdataset1z, kudsdataset1sig,0,kudsdataset1stat);
   TGraphErrors* pkudstot = new TGraphErrors(kudstotdata_np, kudstotdataset4z, kudstotdataset4sig,0,kudstotdataset4stat);
   TGraphErrors* pkalephtot = new TGraphErrors(kalephtotdata_np, kalephtotdataset1z, kalephtotdataset1sig,0,kalephtotdataset1stat);
   TGraphErrors* pkdelphitot = new TGraphErrors(kdelphitotdata_np, kdelphitotdataset1z, kdelphitotdataset1sig,0,kdelphitotdataset1stat);
   TGraphErrors* pkpy6 = new TGraphErrors(26,kHx,kHy,kHex,kHe);
   TGraphErrors* pkAL = new TGraphErrors(kalephtotdata_np+3,kHxAL,kHyAL,kHexAL,kHeAL);
   TGraphErrors* pkSLD = new TGraphErrors(kudstotdata_np+3,kHxSLD,kHySLD,kHexSLD,kHeSLD);
   TGraphErrors* pk = new TGraphErrors(26,kHx,kHy,kHex,kHe);

// KSLD=32, PISLD=38, KAL=29, PIAL=39

   TGraphErrors* ppiudstot = new TGraphErrors(piudstotdata_np, piudstotdataset4z, piudstotdataset4sig,0,piudstotdataset4stat);
 TGraphErrors* ppialephtot = new TGraphErrors(pialephtotdata_np, pialephtotdataset1z, pialephtotdataset1sig,0,pialephtotdataset1stat);
   TGraphErrors* ppidelphitot = new TGraphErrors(pidelphitotdata_np, pidelphitotdataset1z, pidelphitotdataset1sig,0,pidelphitotdataset1stat);

   TGraphErrors* ppipy6 = new TGraphErrors(26,piHx,piHy,piHex,piHe);
   TGraphErrors* ppiAL = new TGraphErrors(pialephtotdata_np+3,piHxAL,piHyAL,piHexAL,piHeAL);
   TGraphErrors* ppiSLD = new TGraphErrors(piudstotdata_np+3,piHxSLD,piHySLD,piHexSLD,piHeSLD);
   TGraphErrors* ppi = new TGraphErrors(26,piHx,piHy,piHex,piHe);
   
   TGraphErrors* ppibase = new TGraphErrors(pibase_np, pibasez, pibasesig,0,0);   
   TGraphErrors* ppibaser= new TGraphErrors(pibase_np, pibasezr, pibasesigr,0,0);

   TGraphErrors* ppiratio = new TGraphErrors(23,pidelphitotdataset1z,piHyr,piHexr,piHer);
   TGraphErrors* pkratio = new TGraphErrors(23,kdelphitotdataset1z,kHyr,kHexr,kHer);
  //39, 29
   TGraphErrors* ppiratioAL = new TGraphErrors(pialephtotdata_np,pialephtotdataset1z,piHyrAL,piHexrAL,piHerAL);
   TGraphErrors* pkratioAL = new TGraphErrors(kalephtotdata_np,kalephtotdataset1z,kHyrAL,kHexrAL,kHerAL);
   //38, 32
   TGraphErrors* ppiratioSLD = new TGraphErrors(piudstotdata_np,piudstotdataset4z,piHyrSLD,piHexrSLD,piHerSLD);
   TGraphErrors* pkratioSLD = new TGraphErrors(kudstotdata_np,kudstotdataset4z,kHyrSLD,kHexrSLD,kHerSLD);

  
  TH1D *hframezpi;
  for(Int_t j=0;j<2;++j)
    hframezpi = new TH1D("hframezpi","hframezpi",20,0.01,1.01);

  TCanvas* canvaszpi = new TCanvas("canvaszpi","canvas zscal",550, 800);
  canvaszpi->SetLogy();
   canvaszpi->SetLogx();
  gPad->SetTickx();
  gPad->SetTicky();
 // ppidelphitot->SetTitle("None");
    TPad* pczpi;
    Float_t shiftzpi = 0.05;
  const Int_t nPadXzpi = 1;
  const Int_t nPadYzpi = 1;
  Float_t noMarginzpi    = 0.0;
  Float_t topMarginzpi   = 0.005;
  Float_t botMarginzpi   = 0.09;
  Float_t leftMarginzpi  = 0.06;
  Float_t rightMarginzpi = 0.01;
  Float_t widthzpi       = (1-rightMarginzpi-leftMarginzpi)/nPadXzpi;
  Float_t heightzpi      = (1-botMarginzpi-topMarginzpi)/nPadYzpi;
  const char* yTitzpi = "  1/#sigma d#sigma/dz ";
  const char* xTitzpi = "z";
  canvaszpi->cd();

  TLatex* latexzpi = new TLatex();
  latexzpi->SetNDC();
  latexzpi->SetTextAlign(22);
  latexzpi->SetTextAngle(0);
  latexzpi->SetTextFont(42);
  latexzpi->DrawLatex(0.8,0.04,xTitzpi);
  latexzpi->SetTextAngle(90);
  latexzpi->SetTextSize(0.055);
  latexzpi->DrawLatex(0.05,0.52,yTitzpi);
  latexzpi->SetTextSize(0.07);
  latexzpi->SetTextAngle(270);
  latexzpi->SetTextSize(0.04);

  Printf("iniciando for");
  for (Int_t i = 1; i < 3; i++) {
    Int_t iX = (i-1)%nPadXzpi;
    Int_t iY = (i-1)/nPadXzpi;
    Float_t x1 = leftMarginzpi + iX*widthzpi;
    if(iX==2) x1-=0.01;//it was 0.01
    else if(iX==1) x1+=0.01;//it was 0.01
    Float_t x2 = leftMarginzpi + (iX + 1)*widthzpi;
    if(iX==0) x2+=0.01;
    else if(iX==1) x2-=0.01;
     Float_t y1 = 1 - topMarginzpi - (iY +1)*heightzpi;
     Float_t y2 = 1 - topMarginzpi - iY*heightzpi;
     pczpi = new TPad("pcMCD","", x1, y1, x2, y2, 0, 0, 0);
     Float_t mBot = noMarginzpi;
     Float_t mTop = noMarginzpi;
     Float_t mLeft = noMarginzpi;
     Float_t mRight = noMarginzpi;
     if(iY==0)       mTop   = shiftzpi;
     if(iY==nPadYzpi-1) mBot   = shiftzpi;
     if(iX==0)       mLeft  = 0.08;
     if(iX==nPadXzpi-1) mRight = 0.08;
     pczpi->SetBottomMargin(mBot);
     pczpi->SetTopMargin(mTop);
     pczpi->SetLeftMargin(mLeft);
     pczpi->SetRightMargin(mRight);
     pczpi->SetFillStyle(0);
     pczpi->SetLeftMargin(0.13);
     if(i==2)
       pczpi->SetBottomMargin(0.08);
     canvaszpi->cd();
     pczpi->Draw();
  }
    pczpi->cd();
    pczpi->SetTickx(1);
    pczpi->SetTicky(1);
    pczpi->SetLogy(1);
    pczpi->SetLogx(1);
    hframezpi->GetYaxis()->SetLabelSize(0.045);
    hframezpi->GetXaxis()->SetTitleSize(0.08);
    hframezpi->GetXaxis()->SetTitle("");
    hframezpi->GetYaxis()->SetTitle("");
    hframezpi->SetLineColor(0);
      hframezpi->Draw();
      hframezpi->GetYaxis()->SetRangeUser(0.0001,1000000);
//  ppibase->GetYaxis()->SetTitle("1/#sigma d#sigma/dz");
//  ppibase->GetXaxis()->SetTitle("z");
//   ppibase->GetXaxis()->SetTitleSize(0.08);
//   ppibase->GetYaxis()->SetTitleSize(0.076);
//  ppibase->SetMarkerStyle(29);
//  ppibase->SetMarkerColor(kWhite);
//  ppidelphitot->GetYaxis()->SetTitle("1/#sigma d#sigma/dz");
//  ppidelphitot->GetXaxis()->SetTitle("z");
  ppiAL->SetMarkerStyle(22);
  ppiSLD->SetMarkerStyle(29);
  pkAL->SetMarkerStyle(26);
  pkSLD->SetMarkerStyle(30);

  ppiudstot->SetMarkerStyle(29);
  ppiudstot->SetMarkerColor(kMagenta);
  ppialephtot->SetMarkerStyle(22);
  ppialephtot->SetMarkerColor(kBlue);
  ppidelphitot->SetMarkerStyle(20);
  ppidelphitot->SetMarkerColor(kGreen);
  ppipy6->SetMarkerStyle(21);
  ppipy6->SetMarkerColor(kRed);
  pkudstot->SetMarkerStyle(30);
  pkudstot->SetMarkerColor(kMagenta);
  pkalephtot->SetMarkerStyle(26);
  pkalephtot->SetMarkerColor(kBlue);
  pkdelphitot->SetMarkerStyle(24);
  pkdelphitot->SetMarkerColor(kGreen);
  pkpy6->SetMarkerStyle(25);
  pkpy6->SetMarkerColor(kRed);
     TLegend* lff2 = new TLegend(0.1,0.7,0.48,0.9);
    lff2->SetHeader("e^{+}e^{-} at #sqrt{s}=91 GeV");// Form("Pythia6, %s at #sqrt{s}=%d GeV",systema,ECM),"C"); //"C" to center header
   TLegend* lfffull = new TLegend(0.1,0.7,0.48,0.9); 
   TLegend* lffempty = new TLegend(0.1,0.7,0.48,0.9);
   lfffull->SetHeader("#pi^{+}+#pi^{-}");
   lffempty->SetHeader("K^{+}+K^{-}");
   lff2->SetTextFont(42);
   lff2->SetTextSize(0.042);
   lff2->SetLineColorAlpha(kWhite,0);
   lff2->SetLineWidth(0);
   lff2->SetLineStyle(3);
   lff2->SetShadowColor(kYellow);
   lff2->SetFillColorAlpha(kWhite,0.05);
   lfffull->SetTextFont(42);
   lfffull->SetTextSize(0.042);
   lfffull->SetLineColorAlpha(kWhite,0);
   lfffull->SetLineWidth(0);
   lfffull->SetLineStyle(3);
   lfffull->SetShadowColor(kYellow);
   lfffull->SetFillColorAlpha(kWhite,0.05);
   lffempty->SetTextFont(42);
   lffempty->SetTextSize(0.042);
   lffempty->SetLineColorAlpha(kWhite,0);
   lffempty->SetLineWidth(0);
   lffempty->SetLineStyle(3);
   lffempty->SetShadowColor(kYellow);
   lffempty->SetFillColorAlpha(kWhite,0.05);
//   lff2->AddEntry(ppi,"Full #pi^{+}+#pi^{-}","p");
//   lff2->AddEntry(pk,"Empty K^{+}+K^{-}","p");
    lfffull->AddEntry(ppiudstot,"SLD total (X1000)","lep");
    lfffull->AddEntry(ppialephtot,"ALEPH total (X100)","lep");
    lfffull->AddEntry(ppidelphitot,"DELPHI total (X10)","lep");
    lfffull->AddEntry(ppipy6,"Pythia6","lep");
    lffempty->AddEntry(pkudstot,"SLD total (X1000)","lep");
    lffempty->AddEntry(pkalephtot,"ALEPH total(X100)","lep");
    lffempty->AddEntry(pkdelphitot,"DELPHI total (X10) ","lep");
    lffempty->AddEntry(pkpy6,"Pythia6","lep");

//  ppibase->Draw("AP");
  ppidelphitot->Draw("same P");
  ppipy6->Draw("same P");
  ppiudstot->Draw("same P");
  ppialephtot->Draw("same P");
  
  pkdelphitot->Draw("P same");
  pkpy6->Draw("same P");
  pkudstot->Draw("same P");
  pkalephtot->Draw("same P");

//  ppiAL->Draw("same P");
//  pkAL->Draw("same P");
//  ppiSLD->Draw("same P");
//  pkSLD->Draw("same P");

  lff2->Draw();
  lfffull->Draw();
  lffempty->Draw();


//inicia plot para razones MC/data

 TH1D *hframeMCD[2];
  TF1 *funo = new TF1("funo","1.0+pol0",0.01,1.0);
  funo->SetLineColor(1);
  funo->SetLineStyle(2);
  funo->SetLineWidth(2);
  for(Int_t j=0;j<2;++j)
    hframeMCD[j] = new TH1D(Form("hframeMCD_%d",j),Form("hframeMCD_%d",j),20,0.01,1.01);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  const Int_t nPadX = 1;
  const Int_t nPadY = 2;
  Float_t noMargin    = 0.0;
  Float_t topMargin   = 0.005;
  Float_t botMargin   = 0.09;
  Float_t leftMargin  = 0.06;
  Float_t rightMargin = 0.01;
  Float_t width       = (1-rightMargin-leftMargin)/nPadX;
  Float_t height      = (1-botMargin-topMargin)/nPadY;

  TCanvas* canvaszMCD = new TCanvas("canvaszMCD","canvas z for MC/DATA",550, 800);
  canvaszMCD->SetFillStyle(0);
  TPad* pcMCD[nPadX*nPadY];
  Float_t shift = 0.05;
  const char* yTitMCD = " Monte Carlo/Data (K^{#pm})  Monte Carlo/Data (#pi^{#pm} )";
  const char* xTitMCD = "z";
  canvaszMCD->cd();

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(22);
  latex->SetTextAngle(0);
  latex->SetTextFont(42);
  latex->DrawLatex(0.8,0.04,xTitMCD);
  latex->SetTextAngle(90);
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.05,0.52,yTitMCD);
  latex->SetTextSize(0.07);
  latex->SetTextAngle(270);
  latex->SetTextSize(0.04);
  Printf("iniciando for");
  for (Int_t i = 1; i < 3; i++) {
    Int_t iX = (i-1)%nPadX;
    Int_t iY = (i-1)/nPadX;
    Float_t x1 = leftMargin + iX*width;
    if(iX==2) x1-=0.01;//it was 0.01
    else if(iX==1) x1+=0.01;//it was 0.01
    Float_t x2 = leftMargin + (iX + 1)*width;
    if(iX==0) x2+=0.01;
    else if(iX==1) x2-=0.01;
     Float_t y1 = 1 - topMargin - (iY +1)*height;
     Float_t y2 = 1 - topMargin - iY*height;
     pcMCD[i-1] = new TPad(Form("pcMCD%d",i),"", x1, y1, x2, y2, 0, 0, 0);
     Float_t mBot = noMargin;
     Float_t mTop = noMargin;
     Float_t mLeft = noMargin;
     Float_t mRight = noMargin;
     if(iY==0)       mTop   = shift;
     if(iY==nPadY-1) mBot   = shift;
     if(iX==0)       mLeft  = 0.08;
     if(iX==nPadX-1) mRight = 0.08;
     pcMCD[i-1]->SetBottomMargin(mBot);
     pcMCD[i-1]->SetTopMargin(mTop);
     pcMCD[i-1]->SetLeftMargin(mLeft);
     pcMCD[i-1]->SetRightMargin(mRight);
     pcMCD[i-1]->SetFillStyle(0);
     pcMCD[i-1]->SetLeftMargin(0.13);
     if(i==2)
       pcMCD[i-1]->SetBottomMargin(0.08);
     canvaszMCD->cd();
     pcMCD[i-1]->Draw();
  }
  for (Int_t index = 1; index < 3; index++) {
    pcMCD[index-1]->cd();
    pcMCD[index-1]->SetTickx(1);
    pcMCD[index-1]->SetTicky(1);
    pcMCD[index-1]->SetLogy(0);
    pcMCD[index-1]->SetLogx(1);
    pcMCD[index-1]->SetGridy();
    hframeMCD[index-1]->GetYaxis()->SetLabelSize(0.076);
    hframeMCD[index-1]->GetXaxis()->SetTitleSize(0.08);
    hframeMCD[index-1]->GetXaxis()->SetTitle("");
    hframeMCD[index-1]->GetYaxis()->SetTitle("");
    hframeMCD[index-1]->SetLineColor(0);
   if(index == 1){
      Printf("if");
      hframeMCD[index-1]->Draw();
      hframeMCD[index-1]->GetYaxis()->SetRangeUser(0.01,2);
      ppiratio->SetMarkerColor(kGreen); //(kBlack);
      ppiratio->SetLineColor(kGreen);
      ppiratio->SetMarkerStyle(20);
      ppiratioAL->SetMarkerStyle(22);
      ppiratioAL->SetMarkerColor(kBlue);
      ppiratioAL->SetLineColor(kBlue);
      ppiratioSLD->SetMarkerStyle(29);
      ppiratioSLD->SetMarkerColor(kMagenta);
      ppiratioSLD->SetLineColor(kMagenta);
      ppiratio->Draw("P same");
      ppiratioAL->Draw("P same");
      ppiratioSLD->Draw("P same");
 
     TLegend* lzMCD = new TLegend(0.1,0.7,0.48,0.9);
     lzMCD->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",systema,ECM),"C"); //"C" to center header
     lzMCD->SetTextFont(42);
     lzMCD->SetTextSize(0.052);
     lzMCD->SetLineColorAlpha(kWhite,0);
     lzMCD->SetLineWidth(0);
     lzMCD->SetLineStyle(3);
     lzMCD->SetShadowColor(kWhite);
     lzMCD->SetFillColorAlpha(kWhite,0.05);
      lzMCD->AddEntry(ppiratio ,"#pi^{#pm}","");
      lzMCD->AddEntry(pkratio ,"K^{#pm}","");
     lzMCD->Draw();
   }else{
     Printf("else");
      hframeMCD[index-1]->GetYaxis()->SetRangeUser(0.01,2.05);
      hframeMCD[index-1]->GetXaxis()->SetLabelSize(0.08);
      hframeMCD[index-1]->Draw();
      funo->Draw("same");
      pkratio->SetMarkerColor(kGreen);
      pkratio->SetLineColor(kGreen);
      pkratio->SetMarkerStyle(24);  //4
      pkratioAL->SetMarkerStyle(26);
      pkratioAL->SetMarkerColor(kBlue);
      pkratioAL->SetLineColor(kBlue);
      pkratioSLD->SetMarkerStyle(30);
      pkratioSLD->SetMarkerColor(kMagenta);
      pkratioSLD->SetLineColor(kMagenta);
      pkratio->Draw("P same");
      pkratioAL->Draw("P same");
      pkratioSLD->Draw("P same");
    }
  }






/*
 TCanvas* canvaszratio = new TCanvas("canvaszratio","canvas zscalratio");
   canvaszratio->SetLogx();
  gPad->SetTickx();
  gPad->SetTicky();
  ppibaser->GetYaxis()->SetTitle("MC/Data");
  ppibaser->GetXaxis()->SetTitle("z");
  ppibaser->SetMarkerStyle(29);
  ppibaser->SetMarkerColor(kWhite);
  ppibaser->SetLineColorAlpha(kWhite,0);
  ppibaser->SetFillColorAlpha(kWhite,0.05);
  ppiratio->SetMarkerColor(kGreen); //(kBlack);
  ppiratio->SetLineColor(kGreen); 
  ppiratio->SetMarkerStyle(20);
  ppiratioAL->SetMarkerStyle(22);
  ppiratioAL->SetMarkerColor(kBlue);
  ppiratioAL->SetLineColor(kBlue);
  ppiratioSLD->SetMarkerStyle(29);
  ppiratioSLD->SetMarkerColor(kMagenta);
  ppiratioSLD->SetLineColor(kMagenta);
  ppibaser->Draw("AP");
  ppiratio->Draw("P same");
  ppiratioAL->Draw("P same");
  ppiratioSLD->Draw("P same");


 TCanvas* canvaszratioK = new TCanvas("canvaszratioK","canvas zscalratioK");
   canvaszratioK->SetLogx();
  gPad->SetTickx();
  gPad->SetTicky();
  ppibaser->GetYaxis()->SetTitle("MC/Data");
  ppibaser->GetXaxis()->SetTitle("z");
  ppibaser->SetMarkerStyle(29);
  ppibaser->SetMarkerColor(kWhite);
  ppibaser->SetLineColorAlpha(kWhite,0);
  ppibaser->SetFillColorAlpha(kWhite,0.05);
  pkratio->SetMarkerColor(kGreen);
  pkratio->SetLineColor(kGreen);
  pkratio->SetMarkerStyle(24);  //4
  pkratioAL->SetMarkerStyle(26);
  pkratioAL->SetMarkerColor(kBlue);
  pkratioAL->SetLineColor(kBlue);  
  pkratioSLD->SetMarkerStyle(30);
  pkratioSLD->SetMarkerColor(kMagenta);
  pkratioSLD->SetLineColor(kMagenta);
  ppibaser->Draw("AP");
  pkratio->Draw("P same");
  pkratioAL->Draw("P same");
  pkratioSLD->Draw("P same");


*/




//inicia plot para comparacion 1,2  with Jets

 TH1D *hframej[2];
 // TF1 *funo = new TF1("funo","1.0+pol0",0.01,1.0);
  funo->SetLineColor(1);
  funo->SetLineStyle(2);
  funo->SetLineWidth(2);
  for(Int_t j=0;j<2;++j)
    hframej[j] = new TH1D(Form("hframej_%d",j),Form("hframej_%d",j),20,0.01,1.01);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
 /*
  const Int_t nPadX = 1;
  const Int_t nPadY = 2;
  Float_t noMargin    = 0.0;
  Float_t topMargin   = 0.005;
  Float_t botMargin   = 0.09;
  Float_t leftMargin  = 0.06;
  Float_t rightMargin = 0.01;
  Float_t width       = (1-rightMargin-leftMargin)/nPadX;
  Float_t height      = (1-botMargin-topMargin)/nPadY;
*/

  TCanvas* canvaszj = new TCanvas("canvaszj","canvas z for jets",550, 800);
  canvaszj->SetFillStyle(0);
  TPad* pcjet[nPadX*nPadY];
 // Float_t shift = 0.05;
  const char* yTitj = " Ratio to inclusive  #kern[0.5]{ (1/N_{jets})dN_{ch}/dz^{ch} } ";
  const char* xTitj = "z^{ch}";
  canvaszj->cd();

 // TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(22);
  latex->SetTextAngle(0);
  latex->SetTextFont(42);
  latex->DrawLatex(0.8,0.04,xTitj);
  latex->SetTextAngle(90);
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.05,0.52,yTitj);
  latex->SetTextSize(0.07);
  latex->SetTextAngle(270);
  latex->SetTextSize(0.04);

  Printf("iniciando for");
  for (Int_t i = 1; i < 3; i++) {
    Int_t iX = (i-1)%nPadX;
    Int_t iY = (i-1)/nPadX;
    Float_t x1 = leftMargin + iX*width;
    if(iX==2) x1-=0.01;//it was 0.01
    else if(iX==1) x1+=0.01;//it was 0.01
    Float_t x2 = leftMargin + (iX + 1)*width;
    if(iX==0) x2+=0.01;
    else if(iX==1) x2-=0.01;
     Float_t y1 = 1 - topMargin - (iY +1)*height;
     Float_t y2 = 1 - topMargin - iY*height;
     pcjet[i-1] = new TPad(Form("pcjet%d",i),"", x1, y1, x2, y2, 0, 0, 0);
     Float_t mBot = noMargin;
     Float_t mTop = noMargin;
     Float_t mLeft = noMargin;
     Float_t mRight = noMargin;
     if(iY==0)       mTop   = shift;
     if(iY==nPadY-1) mBot   = shift;
     if(iX==0)       mLeft  = 0.08;
     if(iX==nPadX-1) mRight = 0.08;
     pcjet[i-1]->SetBottomMargin(mBot);
     pcjet[i-1]->SetTopMargin(mTop);
     pcjet[i-1]->SetLeftMargin(mLeft);
     pcjet[i-1]->SetRightMargin(mRight);
     pcjet[i-1]->SetFillStyle(0);
     pcjet[i-1]->SetLeftMargin(0.13);
     if(i==2)
       pcjet[i-1]->SetBottomMargin(0.08);
     canvaszj->cd();
     pcjet[i-1]->Draw();
  }
  for (Int_t index = 1; index < 3; index++) {
    pcjet[index-1]->cd();
    pcjet[index-1]->SetTickx(1);
    pcjet[index-1]->SetTicky(1);
    if(index == 1) pcjet[index-1]->SetLogy(1);
    if(index==2) pcjet[index-1]->SetGridy();
    hframej[index-1]->GetYaxis()->SetLabelSize(0.076);
    hframej[index-1]->GetXaxis()->SetTitleSize(0.08);
    hframej[index-1]->GetXaxis()->SetTitle("");
    hframej[index-1]->GetYaxis()->SetTitle("");
    hframej[index-1]->SetLineColor(0);
   if(index == 1){
      Printf("if");
      hframej[index-1]->Draw();
      hframej[index-1]->GetYaxis()->SetRangeUser(0.005,100);
     hzchjet2[0]->SetMarkerColor(kBlack);
     hzchjet2[1]->SetMarkerColor(kRed);
     hzchjet2[2]->SetMarkerColor(kGreen+2);
     hzchjet2[3]->SetMarkerColor(kBlue);
     hzchjet2[4]->SetMarkerColor(kCyan);
     hzchjet2[5]->SetMarkerColor(kPink+10);
     hzchjet2[6]->SetMarkerColor(kYellow+2);
     for(Int_t i=0; i<7; i++){
      hzchjet2[i]->SetMarkerStyle(21);
      hzchjet2[i]->Draw("psame");
      }
 
     TLegend* lzj = new TLegend(0.1,0.7,0.48,0.9);
     lzj->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",systema,ECM),"C"); //"C" to center header
     lzj->SetTextFont(42);
     lzj->SetTextSize(0.052);
     lzj->SetLineColorAlpha(kWhite,0);
     lzj->SetLineWidth(0);
     lzj->SetLineStyle(3);
     lzj->SetShadowColor(kWhite);
     lzj->SetFillColorAlpha(kWhite,0.05);
     for(Int_t i=0; i<6; i++){
            lzj->AddEntry(hzchjet2[i],Form("%s",ptint[i]),"lep");
      }
     lzj->Draw();
   }else{
     Printf("else");
      hframej[index-1]->GetYaxis()->SetRangeUser(0.005,2.05);
      hframej[index-1]->GetXaxis()->SetLabelSize(0.08);
      hframej[index-1]->GetXaxis()->SetTitleSize(0.08);
      hframej[index-1]->Draw();
      funo->Draw("same");
   hzjetRat[0]->SetMarkerColor(kBlack);
   hzjetRat[1]->SetMarkerColor(kRed);
   hzjetRat[2]->SetMarkerColor(kGreen+2);
   hzjetRat[3]->SetMarkerColor(kBlue);
   hzjetRat[4]->SetMarkerColor(kCyan);
   hzjetRat[5]->SetMarkerColor(kPink+10);
   hzjetRat[6]->SetMarkerColor(kYellow+2);
     for(Int_t i=0; i<7; i++){
      hzjetRat[i]->SetMarkerStyle(21);
      hzjetRat[i]->Draw("samep");
     }
    }
  }




/*
  TCanvas* canvaszj = new TCanvas("canvaszj","canvas z for jets");
  canvaszj->SetLogy();
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  hzchjet2[0]->SetTitle("Jet Fragmentation function");
  hzchjet2[0]->GetXaxis()->SetTitle("z^{ch}");
  hzchjet2[0]->GetYaxis()->SetTitle("(1/N_{jets})dN_{ch}/dz^{ch}");
  hzchjet2[0]->GetYaxis()->SetRangeUser(0.001,100);
  hzchjet2[0]->GetXaxis()->SetRangeUser(0.01,1.0);
   hzchjet2[0]->GetYaxis()->SetLabelSize(0.04);
   hzchjet2[0]->GetXaxis()->SetLabelSize(0.04);
   hzchjet2[0]->GetXaxis()->SetTitleSize(0.045);
   hzchjet2[0]->GetYaxis()->SetTitleSize(0.045);
  for(Int_t i=0; i<7; i++){
 //  cout<<"i= "<<i+1<<" njetspt= "<<njetspt[i+1]<<endl;
//   hzchjet2[i]->Scale(1./(njetspt[i+1]));
   hzchjet2[i]->SetMarkerStyle(21);
  }
   hzchjet2[0]->SetMarkerColor(kBlack);
   hzchjet2[1]->SetMarkerColor(kRed);
   hzchjet2[2]->SetMarkerColor(kGreen+2);
   hzchjet2[3]->SetMarkerColor(kBlue);
   hzchjet2[4]->SetMarkerColor(kCyan);
   hzchjet2[5]->SetMarkerColor(kPink+10);
   hzchjet2[6]->SetMarkerColor(kYellow+2);
   TLegend* lzj = new TLegend(0.1,0.7,0.48,0.9);
   lzj->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",systema,ECM),"C");
      lzj->SetTextFont(42);
   lzj->SetTextSize(0.03);
   lzj->SetLineColorAlpha(kWhite,0);
   lzj->SetLineWidth(0);
   lzj->SetLineStyle(3);
   lzj->SetShadowColor(kYellow);
   lzj->SetFillColorAlpha(kWhite,0.05);
    for(Int_t i=0; i<6; i++){
            lzj->AddEntry(hzchjet2[i],Form("%s",ptint[i]),"lep");
    }
   hzchjet2[0]->Draw("p");
  for(Int_t i=1; i<6; i++){
   hzchjet2[i]->Draw("psame");
  }
  lzj->Draw();



 TCanvas* crat1 = new TCanvas("crat1","canvas zjrat");
  gPad->SetTickx();
  gPad->SetTicky();
  hzjetRat[0]->GetYaxis()->SetTitle("Ratio to inclusive");
  hzjetRat[0]->GetXaxis()->SetTitle("z");
 for(Int_t i=0; i<7; i++){
   hzjetRat[i]->SetMarkerStyle(21);
  }
   hzjetRat[0]->SetMarkerColor(kBlack);
   hzjetRat[1]->SetMarkerColor(kRed);
   hzjetRat[2]->SetMarkerColor(kGreen+2);
   hzjetRat[3]->SetMarkerColor(kBlue);
   hzjetRat[4]->SetMarkerColor(kCyan);
   hzjetRat[5]->SetMarkerColor(kPink+10);
   hzjetRat[6]->SetMarkerColor(kYellow+2);
   TLegend* lzjr = new TLegend(0.1,0.7,0.48,0.9);
   lzjr->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",systema,ECM),"C");
   lzjr->SetTextFont(42);
   lzjr->SetTextSize(0.03);
   lzjr->SetLineColorAlpha(kWhite,0);
   lzjr->SetLineWidth(0);
   lzjr->SetLineStyle(3);
   lzjr->SetShadowColor(kYellow);
   lzjr->SetFillColorAlpha(kWhite,0.05);
    for(Int_t i=0; i<6; i++){
      lzjr->AddEntry(hzjetRat[i],Form("%s",ptint[i]),"lep");
    }
   hzjetRat[0]->Draw("p");
  for(Int_t i=1; i<6; i++){
   hzjetRat[i]->Draw("psame");
  }
  lzjr->Draw();
*/

//inicia plot para comparacion 1,2  with So

 TH1D *hframe[2];
  //TF1 *funo = new TF1("funo","1.0+pol0",0.001,1.0);
  funo->SetLineColor(1);
  funo->SetLineStyle(2);
  funo->SetLineWidth(2);
  for(Int_t j=0;j<2;++j)
    hframe[j] = new TH1D(Form("hframe_%d",j),Form("hframe_%d",j),20,0.01,1.01);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
/*
  const Int_t nPadX = 1;
  const Int_t nPadY = 2;
  Float_t noMargin    = 0.0;
  Float_t topMargin   = 0.005;
  Float_t botMargin   = 0.09;
  Float_t leftMargin  = 0.06;
  Float_t rightMargin = 0.01;
  Float_t width       = (1-rightMargin-leftMargin)/nPadX;
  Float_t height      = (1-botMargin-topMargin)/nPadY;
*/
  TCanvas* canvaszSo = new TCanvas("canvaszSo","canvas z for So",550, 800);
  canvaszSo->SetFillStyle(0);
  TPad* pcSo[nPadX*nPadY];
 // Float_t shift = 0.05;
  const char* yTitSo = " Ratio to jetty events  #kern[0.5]{ (1/N_{ev})dN_{ch}/dz } ";
  const char* xTitSo = "z";
  canvaszSo->cd();
 
 // TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(22);
  latex->SetTextAngle(0);
  latex->SetTextFont(42);
  latex->DrawLatex(0.8,0.04,xTitSo);
  latex->SetTextAngle(90);
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.05,0.52,yTitSo);
  latex->SetTextSize(0.07);
  latex->SetTextAngle(270);
  latex->SetTextSize(0.04);
 
  Printf("iniciando for");
  for (Int_t i = 1; i < 3; i++) {
    Int_t iX = (i-1)%nPadX;
    Int_t iY = (i-1)/nPadX;
    Float_t x1 = leftMargin + iX*width;
    if(iX==2) x1-=0.01;//it was 0.01
    else if(iX==1) x1+=0.01;//it was 0.01
    Float_t x2 = leftMargin + (iX + 1)*width;
    if(iX==0) x2+=0.01;
    else if(iX==1) x2-=0.01;
     Float_t y1 = 1 - topMargin - (iY +1)*height;
     Float_t y2 = 1 - topMargin - iY*height;
     pcSo[i-1] = new TPad(Form("pcSo%d",i),"", x1, y1, x2, y2, 0, 0, 0);
     Float_t mBot = noMargin;
     Float_t mTop = noMargin;
     Float_t mLeft = noMargin;
     Float_t mRight = noMargin;
     if(iY==0)       mTop   = shift;
     if(iY==nPadY-1) mBot   = shift;
     if(iX==0)       mLeft  = 0.08;
     if(iX==nPadX-1) mRight = 0.08;
     pcSo[i-1]->SetBottomMargin(mBot);
     pcSo[i-1]->SetTopMargin(mTop);
     pcSo[i-1]->SetLeftMargin(mLeft);
     pcSo[i-1]->SetRightMargin(mRight);
     pcSo[i-1]->SetFillStyle(0);
     pcSo[i-1]->SetLeftMargin(0.13);
     if(i==2)
       pcSo[i-1]->SetBottomMargin(0.08);
     canvaszSo->cd();
     pcSo[i-1]->Draw();
  }
  for (Int_t index = 1; index < 3; index++) {
    pcSo[index-1]->cd();
    pcSo[index-1]->SetTickx(1);
    pcSo[index-1]->SetTicky(1);
    if(index == 1) pcSo[index-1]->SetLogy(1);
    if(index==2) pcSo[index-1]->SetGridy();
    hframe[index-1]->GetYaxis()->SetLabelSize(0.076);
    hframe[index-1]->GetXaxis()->SetTitleSize(0.08);
    hframe[index-1]->GetXaxis()->SetTitle("");
    hframe[index-1]->GetYaxis()->SetTitle("");
    hframe[index-1]->SetLineColor(0);
   if(index == 1){
      Printf("if");
      hframe[index-1]->Draw();
      hframe[index-1]->GetYaxis()->SetRangeUser(0.005,1000);
     
     for(Int_t i=0; i<nSobins; i++){
       hzscaffSo2[i]->SetMarkerStyle(20+i);
       hzscaffSo2[i]->SetMarkerColor(99-5*(i));
       hzscaffSo2[i]->Draw("psame");
     }
     TLegend* lzSo = new TLegend(0.1,0.7,0.48,0.9);
     lzSo->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",systema,ECM),"C"); //"C" to center header
     lzSo->SetTextFont(42);
     lzSo->SetTextSize(0.052);
     lzSo->SetLineColorAlpha(kWhite,0);
     lzSo->SetLineWidth(0);
     lzSo->SetLineStyle(3);
     lzSo->SetShadowColor(kWhite);
     lzSo->SetFillColorAlpha(kWhite,0.05);
     for(Int_t i=0; i<nSobins; i++){
            lzSo->AddEntry(hzscaffSo2[i],Form("%s<S_{0pc}<%s",Sopcin[i],Sopcin[i+1]),"lep");
     }
     lzSo->Draw(); 
   }else{
      Printf("else");
      hframe[index-1]->GetYaxis()->SetRangeUser(0.005,2.05);
      hframe[index-1]->GetXaxis()->SetLabelSize(0.08);
      hframe[index-1]->GetXaxis()->SetTitleSize(0.08);
      hframe[index-1]->Draw();
      funo->Draw("same");

      for(Int_t i=0; i<nSobins; i++){ //7
        hzSoRat[i]->SetMarkerStyle(20+i);
        hzSoRat[i]->SetMarkerColor(99-5*(i));
        hzSoRat[i]->Draw("samep");
      }
    }
  }
 













/*
  TCanvas* canvaszSo = new TCanvas("canvaszSo","canvas z for So");
  canvaszSo->Divide(1,2);
//  CanvasPartition(canvaszSo,1,2,lMar,rMar,bMar,tMar);
   canvaszSo->cd(1);
   canvaszSo->SetLogy();
   gPad->SetTickx();
   gPad->SetTicky();
   gStyle->SetOptStat(0);
  hzscaffSo2[0]->SetTitle("Fragmentation functions in S_{0pc} bins");
//  hzscaffSo2[0]->GetXaxis()->SetTitle("z");
  hzscaffSo2[0]->GetYaxis()->SetTitle("(1/N_{ev})dN_{ch}/dz");
  hzscaffSo2[0]->GetYaxis()->SetRangeUser(0.001,100);
  hzscaffSo2[0]->GetXaxis()->SetRangeUser(0.001,1.0);

   hzscaffSo2[0]->GetXaxis()->SetTitleOffset(0.3); 
   hzscaffSo2[0]->GetXaxis()->SetTitleSize(0.09);
   hzscaffSo2[0]->GetXaxis()->SetLabelSize(0.05);
   hzscaffSo2[0]->GetYaxis()->SetLabelSize(0.05);
   hzscaffSo2[0]->GetYaxis()->SetTitleSize(0.09);

  for(Int_t i=0; i<nSobins; i++){
  // hzscaffSo2[i]->Scale(1./(neventsSo[i]));
   hzscaffSo2[i]->SetMarkerStyle(20+i);
   hzscaffSo2[i]->SetMarkerColor(99-5*(i));
  }
   TLegend* lzSo = new TLegend(0.1,0.7,0.48,0.9);
   lzSo->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",systema,ECM),"C"); //"C" to center header
   lzSo->SetTextFont(42);
   lzSo->SetTextSize(0.042);
   lzSo->SetLineColorAlpha(kWhite,0);
   lzSo->SetLineWidth(0);
   lzSo->SetLineStyle(3);
   lzSo->SetShadowColor(kYellow);
   lzSo->SetFillColorAlpha(kWhite,0.05);
    for(Int_t i=0; i<nSobins; i++){
            lzSo->AddEntry(hzscaffSo2[i],Form("%s<S_{0pc}<%s",Sopcin[i],Sopcin[i+1]),"lep");
    }
   hzscaffSo2[0]->Draw("p");
  for(Int_t i=1; i<nSobins; i++){
   hzscaffSo2[i]->Draw("psame");
  }
  lzSo->Draw();
  canvaszSo->Update();

  canvaszSo->cd(2);

// TCanvas* crat2 = new TCanvas("crat2","canvas zjrat2");
  gPad->SetTickx();
  gPad->SetTicky();
  hzSoRat[0]->GetYaxis()->SetTitle("ratio to jetty");
  hzSoRat[0]->GetXaxis()->SetTitle("z");

  hzSoRat[0]->GetYaxis()->SetRangeUser(0.01,2.0);
  hzSoRat[0]->GetXaxis()->SetRangeUser(0.001,1.0);

  hzSoRat[0]->GetXaxis()->SetTitleOffset(0.3);
  hzSoRat[0]->GetXaxis()->SetTitleSize(0.09);
  hzSoRat[0]->GetXaxis()->SetLabelSize(0.05);
  hzSoRat[0]->GetYaxis()->SetTitleSize(0.09);
  hzSoRat[0]->GetYaxis()->SetLabelSize(0.05);

 for(Int_t i=0; i<nSobins; i++){ //7
   hzSoRat[i]->SetMarkerStyle(20+i);
   hzSoRat[i]->SetMarkerColor(99-5*(i));
  }
   TLegend* lzjr2 = new TLegend(0.1,0.7,0.48,0.9);
   lzjr2->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",systema,ECM),"C");
   lzjr2->SetTextFont(42);
   lzjr2->SetTextSize(0.03);
   lzjr2->SetLineColorAlpha(kWhite,0);
   lzjr2->SetLineWidth(0);
   lzjr2->SetLineStyle(3);
   lzjr2->SetShadowColor(kYellow);
   lzjr2->SetFillColorAlpha(kWhite,0.05);
    for(Int_t i=0; i<nSobins; i++){ //6
      lzjr2->AddEntry(hzSoRat[i],Form("%s-%s",Sopcin[i],Sopcin[i+1]),"lep");
    }
   hzSoRat[0]->Draw("p");
  for(Int_t i=1; i<nSobins; i++){
   hzSoRat[i]->Draw("psame");
  }
  lzjr2->Draw();
   canvaszSo->Update();
*/

/*
//---------------------------------
  TCanvas* canvaszSoNbin[10];
  TLegend* lzSoNbins[10];
  for(Int_t j=0; j<10; j++){
    canvaszSoNbin[j] = new TCanvas(Form("canvaszSoNbin%d",j),Form("canvas z for Sobin=%d vs Nch",j));
    canvaszSoNbin[j]->SetLogy();
    gPad->SetTickx();
    gPad->SetTicky();
    gStyle->SetOptStat(0);
  hzscaffSoNbin2[j][2]->SetTitle(Form("fragmentation function for %s<S_{0pc}<%s",Sopcin[j],Sopcin[j+1]));
  hzscaffSoNbin2[j][2]->GetXaxis()->SetTitle("z");
  hzscaffSoNbin2[j][2]->GetYaxis()->SetTitle("(1/N_{ev})dN_{ch}/dz");
  hzscaffSoNbin2[j][2]->GetYaxis()->SetRangeUser(0.001,1000);
  hzscaffSoNbin2[j][2]->GetXaxis()->SetRangeUser(0.001,1.0);
  hzscaffSoNbin2[j][2]->GetYaxis()->SetLabelSize(0.04);
  hzscaffSoNbin2[j][2]->GetXaxis()->SetLabelSize(0.04);
  hzscaffSoNbin2[j][2]->GetXaxis()->SetTitleSize(0.045);
  hzscaffSoNbin2[j][2]->GetYaxis()->SetTitleSize(0.045);
  for(Int_t i=2; i<nMultbins; i++){
   if(neventsSoNch[j][i]>0)hzscaffSoNbin2[j][i]->Scale(1./(neventsSoNch[j][i]));
   hzscaffSoNbin2[j][i]->SetMarkerStyle(18+i);
  }
   hzscaffSoNbin2[j][2]->SetMarkerColor(kMagenta-1);
   hzscaffSoNbin2[j][3]->SetMarkerColor(kViolet+3);
   hzscaffSoNbin2[j][4]->SetMarkerColor(kViolet+2);
   hzscaffSoNbin2[j][5]->SetMarkerColor(kBlue);
   hzscaffSoNbin2[j][6]->SetMarkerColor(kGreen+2);
   hzscaffSoNbin2[j][7]->SetMarkerColor(kGreen);
   hzscaffSoNbin2[j][8]->SetMarkerColor(kSpring+7);
   hzscaffSoNbin2[j][9]->SetMarkerColor(kYellow+1);
   hzscaffSoNbin2[j][10]->SetMarkerColor(kOrange);
   hzscaffSoNbin2[j][11]->SetMarkerColor(kOrange+10);
   hzscaffSoNbin2[j][12]->SetMarkerColor(kRed);
   lzSoNbins[j] = new TLegend(0.1,0.7,0.48,0.9);
   lzSoNbins[j]->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",systema,ECM),"C"); //"C" to center header
   lzSoNbins[j]->SetTextFont(42);
   lzSoNbins[j]->SetTextSize(0.042);
   lzSoNbins[j]->SetLineColorAlpha(kWhite,0);
   lzSoNbins[j]->SetLineWidth(0);
   lzSoNbins[j]->SetLineStyle(3);
   lzSoNbins[j]->SetShadowColor(kYellow);
   lzSoNbins[j]->SetFillColorAlpha(kWhite,0.05);
    for(Int_t i=2; i<nMultbins; i++){
            lzSoNbins[j]->AddEntry(hzscaffSoNbin2[j][i],Form("%1.0f#leq N_{ch}<%1.0f",Multbins[i],Multbins[i+1]),"lep");
    }
   hzscaffSoNbin2[j][2]->Draw("p");
  for(Int_t i=2; i<nMultbins; i++){
   hzscaffSoNbin2[j][i]->Draw("psame");
  }
  lzSoNbins[j]->Draw();
 }


*/

}







void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;
   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;
   for (Int_t i=0;i<Nx;i++) {
      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }
      for (Int_t j=0;j<Ny;j++) {
         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }
         C->cd(0);
         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);
         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);
         pad->Draw();
      }
   }
}
