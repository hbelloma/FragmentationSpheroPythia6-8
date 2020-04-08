#include <TFile.h>
#include <TList.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH2D.h>
#include <TProfile.h>
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <string>


void MakefiletogetPerC(const Char_t * inFileName1="ADistribfuncMBmill_eenfj.root"){
  TFile *fin1 = TFile::Open(inFileName1);
  if(!fin1)
   return;
  TFile* outFile = new TFile("HistPercPP_eenfj30Mill.root", "RECREATE");
 
  const Int_t nMultbins=13;
  Double_t Multbins[nMultbins+1]={ 0.0, 1.0, 4.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 140.0};

  const Int_t npc= nMultbins-1; // last was 12
  TH1D* histSOT[npc+1];
  
  Printf("loop npc");
  for(Int_t i=0; i<npc+1; i++){
   histSOT[i]=(TH1D*)fin1->Get(Form("hSOT%d",i));
  }
  TH1F* hSOT[npc+1];
  for(Int_t i=0; i<npc+1; i++){
   hSOT[i]= new TH1F(Form("hSOTPerc%d",i),";;", 1000,0,1);
  }
  printf("out of loop npc ...\n");
  Int_t binstot =1000;
  printf("binstot=%d\n",binstot);

  Double_t intehistSOT[npc+1]; 
  Double_t conthistSOT[npc+1];
  Double_t contSOT[npc+1];
  Int_t intehistSOTt[npc+1];
  printf("entrando  loop of normalization\n");
 for(Int_t i=0; i<npc+1; i++){
  cout<<"i= "<<i<<endl;
  intehistSOTt[i]=histSOT[i]->Integral(1,binstot+1);
  cout<<"intehistSOTt[i]= "<<intehistSOTt[i]<< endl;
  for(Int_t bin=1; bin<binstot+1; bin++){
   intehistSOT[i]=histSOT[i]->Integral(1,bin+1);
   conthistSOT[i]=histSOT[i]->GetBinContent(bin);
   if(i>0){
        if( intehistSOTt[i]>1){
	contSOT[i]=intehistSOT[i]/intehistSOTt[i];
        hSOT[i]->SetBinContent(bin,contSOT[i]);
	}
       	else{ hSOT[i]->SetBinContent(bin,0); 
		cout<<"this was set to cont 0"<<endl;}

   }
   else{ hSOT[i]->SetBinContent(bin,0); cout<<"i= 0"<<endl;}
  }
 }
  
 cout<<"now printing canvas"<<endl; 
 // */
 TCanvas* csot = new TCanvas("csot","SoT");
 hSOT[2]->Draw();

 outFile->cd();
 for(Int_t i=0; i<npc+1; i++){
  hSOT[i]->Write();
 }
 outFile->Close();


}
