//____________________________________________________________________
//
// Generation process: MB looking for fragmentation functions
// Using Pythia6 with AliROOT
//
// Author: Héctor Bello Martínez <hbelloma@cern.ch>
// Update: 2018-10-12
//
// Modification by: 
// Date:
//
//-------------------------------------------------------------------
// To make an event sample (of size 1000) do
//
//    shell> aliroot
//    root [0] .L pythia6Distrfunc_eenofjnewz.C
//    root [1]  gSystem->Load("libpythia6");
//    root [2]  gSystem->Load("libAliPythia6");
//    root [3] makeEventSample(1000)
//
// To start the tree view on the generated tree, do
//  Works only if the file has been produced
//    shell> aliroot
//    root [0] .L pythia6Distrfunc_eenofjnewz.C
//    root [1]  gSystem->Load("libpythia6");
//    root [2]  gSystem->Load("libAliPythia6");
//    root [3] showEventSample()
//
// The following session:
//    shell> aliroot
//    root [0] .L pythia6Distrfunc_eenofjnewz.C
//    root [1]  gSystem->Load("libpythia6");
//    root [2]  gSystem->Load("libAliPythia6");
//    root [3] .x pythiaExample(500)
// will execute makeEventSample(500) and showEventSample()
//
//____________________________________________________________________


#ifndef __CINT__
#include "TApplication.h"
#include "TPythia6.h"
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TDecayChannel.h"
#include "TParticlePDG.h"
#include <cstdlib>

//#include "fastjet/ClusterSequence.hh"
//#include <fastjet/PseudoJet.hh>
#include "TMCParticle.h"
//using namespace fastjet;
using namespace std;
#endif

//------------------------ define parametros ---
#define FILENAME    "ADistribfuncMBmill_eenfjnewz.root"
// ADistribfuncMBmill AApy6H2ZZ24MU_onlyggh_onlymufromz.root
#define OUTDIREC     "MBfragmentationFPY6_eenfjnewz"
#define OUTFILEPOSTPRO "ADistribfuncMB_postprocmill_eenfjnewz.root"
#define IMAGEF     "pdf"
#define ECM        91 
#define TREENAME   "tree"
#define BRANCHNAME "particles"
#define INTEVT 10051



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
 const char *ptint[7]={"p^{ch jet}_{T} 5-100 GeV/c (inclusive)","p^{ch jet}_{T} 5-10 GeV/c","p^{ch jet}_{T} 10-15 GeV/c","p^{ch jet}_{T} 15-20 GeV/c","p^{ch jet}_{T} 20-30 GeV/c","p^{ch jet}_{T} 30-100 GeV/c","p^{ch jet}_{T} 40-100 GeV/c"};
 const char *Sopcin[nSobins+1]={"0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"};


TH1D* hSOTPerc[nMultbins];
Double_t GetESPerc( Double_t valES, Int_t multbin); 

#define HISTNAME   "ptSpectraformu"
#define HISTNAMEZ  "ptSpectraforZ"
#define PDGNUMBER 13 


// This function just load the needed libraries if we're executing from
// an interactive session.
void loadLibraries()
{
#ifdef __CINT__
  // Load the Event Generator abstraction library, Pythia 6
  // library, and the Pythia 8 interface library in AliRoot.
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  gSystem->Load("libJETAN");
  gSystem->Load("libFASTJETAN");

  gSystem->Load("libEG");
 // gSystem->Load("libqpythia");
  
  gSystem->Load("libEGPythia6"); 
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
//  gSystem->Load("libpythia8");
//  gSystem->Load("libAliPythia8");
//  gSystem->Load("libqpythia.so");


#endif
}
//---------------------------------------------------------------------------
void CreateDir(const Char_t* dirName)
{
   TString pwd(gSystem->pwd());
   gSystem->cd(pwd.Data());

   if(gSystem->cd(dirName)) {
   gSystem->cd(pwd.Data());
   } else {
   gSystem->mkdir(dirName, kTRUE); // kTRUE means recursive
   }
}
//----------------------------------------------------------------------------
/*
void DeclareHistsEqual(){
 TH1D* hcontadoreventos = new TH1D("hcontadoreventos","",10,0,10);
 TH1D* hNch = new TH1D("hNch","",300,0,300);
 TH1D* hSOT[nMultbins];
 TH1D* hSTT[nMultbins];

 for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){
    hSOT[ibin] = new TH1D(Form("hSOT%d",ibin),"true spherocity", 1000,0,1.0);
    hSTT[ibin] = new TH1D(Form("hSTT%d",ibin),"true sphericity", 1000,0,1.0);
 }
}
*/
//______________________________
TString CompleteDirName( TString nameOut3 ){
   TString nameOut = nameOut3;
   nameOut += "histos";
   return nameOut;
}
//---------------------------------------------------------------------------
Double_t GetSpherocity(Double_t pxA[], Double_t pyA[], Double_t sumaptev, Int_t Mult, Double_t size_stepe){
  
  Double_t spherocity=-10;
  if(Mult<3) return -0.5;

  Double_t pFull = 0;
  Double_t Spherocitytr = 2;
  //Getting spherocity
  for(Int_t i = 0; i < 360/(size_stepe); ++i){
	Double_t numerador = 0;
	Double_t phiparam  = 0;
	Double_t nx = 0;
	Double_t ny = 0;
	phiparam=((TMath::Pi()) * i * size_stepe ) / 180; // parametrization of the angle
	nx = TMath::Cos(phiparam);            // x component of an unitary vector n
	ny = TMath::Sin(phiparam);            // y component of an unitary vector n
	for(Int_t i1 = 0; i1 < Mult; ++i1){
            numerador += TMath::Abs(ny * pxA[i1] - nx * pyA[i1]);
        }
        pFull=TMath::Power( (numerador / sumaptev),2 );
	 if(pFull < Spherocitytr) Spherocitytr = pFull;//minimization of pFull
   }
   spherocity=((Spherocitytr)*TMath::Pi()*TMath::Pi())/4.0;
   return spherocity;

}

//____________________________________________________________
Int_t GetMeasMultiplicityBin(Double_t multIn){
	Int_t binOut = -1;
		for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){
			if( (multIn >= Multbins[ibin]) && (multIn < Multbins[ibin+1]) )
				binOut = ibin;
		}
	return binOut;
}
//_____________________________________________________________
Double_t GetESPerc( Double_t valES, Int_t multbin){
	if(multbin<0)
		return -1.0;
	Double_t porcentaje = -1.0;
   for(Int_t bin = 1; bin <= hSOTPerc[multbin]->GetNbinsX(); ++bin){
	Double_t minEs = 0.0;
	Double_t maxEs = 0.0;
	minEs = hSOTPerc[multbin]->GetBinLowEdge(bin);
	maxEs = minEs + hSOTPerc[multbin]->GetBinWidth(bin);
	if( (valES >= minEs) && (valES < maxEs) )
	porcentaje = hSOTPerc[multbin]->GetBinContent(bin);
    }
	return porcentaje;
}


//--------------------------------------------------------------------------------



// nEvents is how many events we want.
int makeEventSample(Int_t nEvents)
{

  Bool_t isESPerc=kFALSE;
  Bool_t isnotequalized=kTRUE;
  Bool_t isfast=kFALSE; //kTRUE;// kFALSE; 
  if(nEvents<101) isfast=kTRUE;

// Load needed libraries
  loadLibraries();
  
  Double_t size_step = 0.1;   //this is for the computation of spherocity


  //if(IsESPerc){
  TFile * finPercent = TFile::Open("HistPercPP.root");
  //}
  TDatabasePDG *pdg = TDatabasePDG::Instance(); // TDataBasePDG contains info on all particle properties
  TPythia6* pythia = new TPythia6;

//-------Realistic process during the generation (1=on,0=off) ---------------
///*  //Plot 1 W/o MPI  //see Pythia62 pag320
    pythia->SetMSTP(61,1);  // switch off ISR
    pythia->SetMSTP(71,1);  // switch off FSR
    pythia->SetMSTP(81,1);  // switch off multiple interactions
    pythia->SetMSTP(111,1); // Switch off fragmentation
//*/

 //------- controled subprocess ----------------------------------------

  //for e+e- according Aliroot Pythia6
   pythia->SetMSEL(1);
    for (Int_t i=1; i<= 200; i++) {
     pythia->SetMSUB(i,0);
    }

pythia->SetMSTJ(101, 2); //qqbar,qqbarg,qqbargg,qqbarq'qbar'events are generated according to second-order QCD
pythia->SetMSTJ(104, 6); //to allow ttbar production
pythia->SetMSTJ(107, 1); //to include initial-state photon radiation

//pythia->SetPARJ(123,91.187); // to change the Z0 mass
//pythia->SetPARJ(124,2.489); //width of Z0 
//pythia->SetPARJ(81,0.3); //to change the parton-shower Λ value
//pythia->SetPARJ(82,1.5); //to change the parton-shower cut-off.


/*pythia->SetMSUB(11, 1); //qiqj->qiqj
pythia->SetMSUB(11, 1); //qiqj->qiqj
   pythia->SetMSUB(12, 1); //qibarqi->qkbarqk
   pythia->SetMSUB(13, 1); //qibarqi->gg
   pythia->SetMSUB(28, 1); //qig->qig
   pythia->SetMSUB(53, 1); //gg->qibarqi
   pythia->SetMSUB(68, 1); //gg->gg
*/

// No weak decays
  pythia->SetMDCY(130 ,1,0); // k0L
  pythia->SetMDCY(321 ,1,0); // k+
  pythia->SetMDCY(313 ,1,0); // k*0
  pythia->SetMDCY(323 ,1,0); // k*+
  pythia->SetMDCY(311 ,1,0); // k0
  pythia->SetMDCY(310 ,1,0); // K0S
  pythia->SetMDCY(3122,1,0); // Lambda
  pythia->SetMDCY(3112,1,0); // sigma -
  pythia->SetMDCY(3222,1,0); // sigma +
  pythia->SetMDCY(3212,1,0); // sigma 0
  pythia->SetMDCY(3312,1,0); // xi - 
  pythia->SetMDCY(3322,1,0); // xi 0
  pythia->SetMDCY(3334,1,0); // omega-

  /*
   SetMSTP(41,1);   // all resonance decays switched on
    // Select structure function 
  pythia->SetMSTP(52,2);
  pythia->SetMSTP(51,kCTEQ5L);//kCTEQ5L; AliStructFuncType::PDFsetIndex(strucfunc));
    // Particles produced in string fragmentation point directly to either of the two endpoints of the string (depending in the side they were generated from).
  pythia->SetMSTU(16,2);
// */
//  CR color reconnection strength
   //pythia->SetMSTP(95,0);

  //-----------  ATLAS TUNE

// for lund frag mod see py82_arx0603175 p436

 // B) Hard QCD subprocess 
 // pythia->SetMSEL(2); //use MSEL=2 insted 0/1 (agregue processos MSUB 91-96)
/*
   pythia->SetMSUB(11, 1); //qiqj->qiqj
   pythia->SetMSUB(12, 1); //qibarqi->qkbarqk
   pythia->SetMSUB(13, 1); //qibarqi->gg
   pythia->SetMSUB(28, 1); //qig->qig
   pythia->SetMSUB(53, 1); //gg->qibarqi
   pythia->SetMSUB(68, 1); //gg->gg
 */


//----- setting general switches via SetMSTP()

// pythia->SetMSTP(1,3);//max num of gen D=3 set<4
// pythia->SetMSTP(2,1);//calculation alpha_s D=1(0=fixed,1=1stOrder,2=2ndOrder)
// pythia->SetMSTP(3,2);//Lambda in alpha_s D=2 (acording pdf)
// pythia->SetMSTP(4,0);//neutral higgs sector D=0 StandarMod couplings
// pythia->SetMSTP(5,0);//anomalous couplings see pythia6 ec 8.6.5
// pythia->SetMSTP(7,0);// HF subprocess D=0
// pythia->SetMSTP(8,0);// EW param D=0 (a_EM(q2) and fixed sin2theta_W)
// pythia->SetMSTP(9,0);//4th generation inclusion D=0 (no)
// pythia->SetMSTP(11,1);//use electron pdf in e+e- and ep interact.
// pythia->SetMSTP(12,0);// use e´,e+,quark and gluon dist func inside e- D=0 
// pythia->SetMSTP(13,1);//choice of Q2 range for e- to radiate gammas
// pythia->SetMSTP(14,30);//structure of inc. phot beam (D=30)
// pythia->SetMSTP(15,0);//anomalous photon component
// pythia->SetMSTP(16,1);//fract of moment tak by phot radiated of a lept
// pythia->SetMSTP(17,4);//extract factor of virtualphotons
// pythia->SetMSTP(18,3);//choice of pTmin for dir procc
// pythia->SetMSTP(19,4);//partonic cross section in DIS process 99
// pythia->SetMSTP(20,3);//supress of reslvd cross sect to comp overlap
// pythia->SetMSTP(21,1);//nature of fermion fermion scatt
// pythia->SetMSTP(22,0);//spetial overide of norm Q2 se sec 10.4 pythia6
// pythia->SetMSTP(23,1);// DIS (10y83)
// pythia->SetMSTP(25,0);//ang decay corr in Higgs dec to 4fermion(0=scalar dec)
// pythia->SetMSTP(31,1);//param of tot elast and diffrac cross sect 
// pythia->SetMSTP(32,8);//Q2 def in hard scatt
// pythia->SetMSTP(33,0);//K factors for hard cross section
// pythia->SetMSTP(34,1);//interference term in matrix elemnts of QCD procc D=1 inc
// pythia->SetMSTP(35,0);//threshold of HF production
// pythia->SetMSTP(36,2);//reg of alpha_s in Q2->0
// pythia->SetMSTP(37,1);//incl of quark masses in h0 prod
// pythia->SetMSTP(38,5);//q-loop masses in box (5-efective masses)
// pythia->SetMSTP(39,2);//q2 for pd and IS part show
// pythia->SetMSTP(40,0);//coloumb correction in W+W-
// pythia->SetMSTP(41,2);//master switch of all resonance decays
// pythia->SetMSTP(42,1);//mass treatment for finite width FS resonance(D=1 Breit-Wigner)

  // ----------- try uncoment this two --------------------
//   pythia->SetMSTP(43,2);//Z0/gamma* interference matrix D=3 (D=2 only Z0) 
//   pythia->SetMSTP(44,2);//Z'/Z0/gamm* interf matrix elements D=7


// ----------------------------
// ######    Important  ####### 
// check switches MSTU and MSTJ
// ----------------------------


//  pythia.SetMSTP(51,partondist);  // Use structure fcn (11 for GRV LO)
//
// Parton distribution functions PDF
// LO  EHLQ, DO, GRV 92L, CTEQ 3L, GRV94L and CTEQ5L 
// NLO CTEQ 3D and GRV94D
// MSTP(56) = 2, MSTP(54) = 2 and  a pion distribution code in MSTP(53),
// The Q0 scale can be set freely in the PARP(15) parameter

  //--- Setting decay modes by SetMDME()

// set (and use) random gen. seed to fixed value for reproducible results
//  pythia->SetMRPY(1,12345678);
//  pythia->SetMRPY(2,0);
  //--------------------------------------------

  pythia->Initialize("cms", "e+", "e-", ECM);

 
  // Open an output file
  TFile* file = TFile::Open(FILENAME, "RECREATE");
  if (!file || !file->IsOpen()) {
   // Error("makeEventSample", "Couldn't open file %s", FILENAME);
    return 1;
  }
 // if(isnotequalized){
   /*TString nameOut2 = "Equal_";
   nameOut2 = CompleteDirName( nameOut2 );
   CreateDir(nameOut2.Data());
   TFile* dataOutFile = new TFile(Form("%s/histosforEqual.root", nameOut2.Data()), "RECREATE");
  //declare histograms
//  */
   TH1D* hcontadoreventos = new TH1D("hcontadoreventos","",10,0,10);
   TH1D* hNch = new TH1D("hNch","",300,0,300);
   TH1D* hSOT[nMultbins];
  // TH1D* hSTT[nMultbins];


//   TH1D* hSoTInNch[nMultbins][nSobins+1];
   for( Int_t ibin = 0; ibin < nMultbins; ++ibin ){
     hSOT[ibin] = new TH1D(Form("hSOT%d",ibin),"true spherocity", 1000,0,1.0);
   //  hSTT[ibin] = new TH1D(Form("hSTT%d",ibin),"true sphericity", 1000,0,1.0);

/*
 * for(Int_t i_So = 0; i_So < nSobins; ++i_So ){ 
			hSoTInNch[ibin][i_So] = 0;
			hSoTInNch[ibin][i_So]  = new TH1D(Form("hESTInNch%dES%d",ibin,i_So),"",1000,0.0,1.0);
      }
      */
   }
 //  if(isESPerc){
    for( Int_t i_mult = 0; i_mult < nMultbins; ++i_mult ){
	hSOTPerc[i_mult] = 0;
	hSOTPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSOTPerc%d",i_mult));
    }
 //  }

  // DeclareHistEqual();
  // Make a tree in that file ...
  // read https://root.cern.ch/trees-five-steps
  TTree* tree = new TTree(TREENAME, "Pythia 6 tree");

  int NbinIVplot = 100;
  int Nbinsz = 100;
  int Nbins4mu = 500;  
  int NbinsH =500;
  int Nbinsy =50;
  // Histograms
 TH1D* hevents=new TH1D("hevents","Eventos",10,0,10);
 TH1D* heventsid=new TH1D("heventsid","Eventos por clase de part",10,0,10);
 TH1D* heventsSo=new TH1D("heventsSo","Eventos with So",10,0,10);
 TH1D* hjets=new TH1D("hjets","Jets producidos",10,0,10);
 TH1D* hjetsptcut=new TH1D("hjetsptcut","jets in a pt cut",10,0,10);
 TH1F* invptch = new TH1F("invptch","Invariant yield", nPtBins,xBins); 
 TH1F* invptchid[10];
 TH1F* hzscaffid[10];
 TH1F* hzchjetid[10];
 TH1F* hzchjet[7];
 TH1F* hzscaffSo[nSobins];
 TH1D* heventsSoN[nSobins];
 TH1F* hzscaffSoNbin[nSobins][nMultbins];


 for(Int_t i=0; i<10; i++){
       invptchid[i] = new TH1F(Form("invptchid%d",i),Form("%s Invariant yield",idt[i]), nPtBins,xBins);
       hzscaffid[i] = new TH1F(Form("hzscaffid%d",i),Form("%s fragmentation function",idt[i]),22,0.0,1.1001); //200,0,2
       hzchjetid[i] = new TH1F(Form("hzchjetid%d",i),Form("%s jet fragmentation function",idt[i]),22,0.0,1.1001);
 }
 for(Int_t i=0; i<4; i++){
       hzchjet[i]=new TH1F(Form("hzchjet%d",i), Form("%s jet fragmentation",ptint[i]),22,0.0,1.1001);
 }
 for(Int_t i=4; i<7; i++){//diferents limmits
       hzchjet[i]=new TH1F(Form("hzchjet%d",i), Form("%s jet fragmentation",ptint[i]),11,0.0,1.1001);
 }
 for(Int_t i=0; i<nSobins; i++){
    if(isESPerc)
       hzscaffSo[i]=new TH1F(Form("hzscaffSo%d",i),Form("fragmentation function for %s<S_{0pc}<%s",Sopcin[i],Sopcin[i+1]),22,0.0,1.1001);
    else
       hzscaffSo[i]=new TH1F(Form("hzscaffSo%d",i),Form("fragmentation function for %f<S_{0}<%f",Sobins[i],Sobins[i+1]),22,0.0,1.1001);
    for(Int_t j=0; j<nMultbins; j++){
      if(isESPerc)
      hzscaffSoNbin[i][j]=new TH1F(Form("hzscaffSoNbin%d_%d",i,j),Form("fragmentation function for %s<S_{0pc}<%s and %f#leq N_{ch}<%f",Sopcin[i],Sopcin[i+1],Multbins[j],Multbins[j+1]),22,0.0,1.1001);
      else
      hzscaffSoNbin[i][j]=new TH1F(Form("hzscaffSoNbin%d_%d",i,j),Form("fragmentation function for %f<S_{0}<%f and %f#leq N_{ch}<%f",Sobins[i],Sobins[i+1],Multbins[j],Multbins[j+1]),22,0.0,1.1001);
    }
    heventsSoN[i]=new TH1D(Form("heventsSoN%d",i),Form("Eventos with %s<So<%s",Sopcin[i],Sopcin[i+1]),14,0,14);
 }


 TH1F* hetach = new TH1F("hetach", "Pseudorapidity", 240, -12., 12.);
 TH1F* hych = new TH1F("hych", "rapidity", 240, -12., 12.);
 TH1F* hphich = new TH1F("hphich", "phiangle", 700, 0.0, 7.0);
 TH1F* hzscaff=new TH1F("hzscaff","dN/dz", 22, 0.0,1.1001);
 TH1F* hzscaffcop=new TH1F("hzscaffcop","dN/dz", 22,0.0,1.1001);

 TH2F* hyphiNfj= new TH2F("hyphiNfj"," y vs #phi vs N_{ch} Fastjet",100,-7.0,7.0,100,0.0,7.0);
 TH2F* hyphiwptfj= new TH2F("hyphiwptfj","y vs #phi weigthed by p_{T} Fastjet",100,-7.0,7.0,100,0.0,7.0);
 TH2F* hyphiN= new TH2F("hyphiN","y vs #phi vs N_{ch}",100,-7.0,7.0,100,0.0,7.0);
 TH2F* hyphiwpt= new TH2F("hyphiwpt","y vs #phi weigthed by p_{T}",100,-7.0,7.0,100,0.0,7.0);

//not usefull
 TH3F* hyphip= new TH3F("hyphip","jets momentum",50,-5.0,5.0,50,-1.0,9.0,20,0.0,150.0);
 TH3F* hyphipt= new TH3F("hyphipt","pt jets",50,-5.0,5.0,50,-1.0,9.0,20,0.0,80.0);
  

  invptch->Sumw2();
  invptch->SetYTitle("1/(2#pip_{T}N_{ev}) d^{2}N_{ch}/dp_{T}d#eta");
  invptch->SetXTitle("p_{T}(GeV/c)");
for(Int_t i=0; i<10; i++){
    invptchid[i]->Sumw2();
    invptchid[i]->SetYTitle("dN_{ch}/dp_{T}");
    invptchid[i]->SetXTitle("p_{T}");
}


  //Array for TMCParticle objects (only pythia6). 
  TClonesArray* particles = (TClonesArray*)pythia->GetListOfParticles();
  tree->Branch(BRANCHNAME, &particles);
  //rama de rapidity
   Float_t fy=-999999; 
   // Need to fix-----------------------ME
   // TLeaf *fyl=particles("fy",&fy);
//   particles->Leaf("fy",&fy);
 //  TLeaf(TBranch* BRANCHNAME, "fy", "f/y")
   TBranch *newBranch = tree->Branch("fy", &fy, "fy/F");

//  Float_t dphi;
//  tree->Branch("dphi", &dphi, "dphi/F");


  pythia->Pystat(4);  // prints a table of kinem cuts CKIN(I) p.135 and p.251

  Int_t nHiggs=0; // contador de Higgs totales
  Int_t multhmax=0;
  Int_t maxjetpt=0;
  Int_t evtmulthmax=0;
  Int_t evtmaxjetpt=0;
  Int_t evtspheromax=0;
  Int_t evtspheromin=0;
  Int_t evtmercedes=0;
  Double_t spheromax=0;
  Double_t spheromin=1;


 Int_t eventidcount[10];
 for(Int_t i=0; i<11; i++){eventidcount[i]=0;}

  // Now we make some events
  for (Int_t ievt = 0; ievt < nEvents; ievt++) {
    // Show how many events we have.
    if(isfast) cout << "Event # " << ievt << endl;
    else{if(ievt % 100 == 0) cout << "Event # " << ievt << endl;}

    // Make event.
    pythia->GenerateEvent();

    //----For 1 event --------
    if (ievt < 1) {
      // print out particle listing information if desired (this shows decay modes and code numbers)
      pythia->Pylist(1);   // list full pythia generated event
     // pythia->Pylist(2);   // list of particle created during the history of the first event
//      pythia->Pylist(3);  // V matrix contains decay vertice
//      pythia->Pylist(12);  // list of all decay modes for all defined particles
    }
    hevents->Fill(0);
   
//    vector <PseudoJet> jparticles;
   // inline void reset(double px, double py, double pz, double E);
//    vector <PseudoJet> jpions;
/*  vector <PseudoJet> jkaons;
    vector <PseudoJet> jprotons;
    vector <PseudoJet> jlambdas;
    vector <PseudoJet> jxis;
    vector <PseudoJet> jomegas;
    vector <PseudoJet> jK0s;
    vector <PseudoJet> jphis;
*/
    //    TLorentzVector vechiggs,vecza,veczb,vecmum1,vecmum2,vecmup1,vecmup2,veczatomu,veczbtomu,vechtoz,vecto4mu;
    
    Int_t np = particles->GetEntriesFast();
    Int_t nmultch=0;
    Int_t nh0=0;

    Double_t etach=-999.;
    Double_t ych=-999.;
    Double_t phichmv=-999.;
    Double_t phich=-999.;
    // for So eq
    Int_t counter=0;
    Double_t sumapt=0;
    Double_t *pxA=new Double_t[np];
    Double_t *pyA=new Double_t[np];

//    cout<<"starting loop of particles np= "<<np<<endl;
    // loop over particles
    for (Int_t ip=0; ip<np; ip++){
      
      TMCParticle* part = (TMCParticle*) particles->At(ip);
      if(part->GetParent()==0) continue;
      TMCParticle* parent=dynamic_cast<TMCParticle*>((*particles)[part->GetParent()-1]);     
      Int_t idp = part->GetKF();//  GetPdgCode(); //(Pythia8) pdg number
      if(TMath::Abs(idp)>4000012) continue; // protection to not known pdg codes in root
    //  cout<<"idp="<<idp<<endl;
    //  cout<<"parent="<<part->GetParent()<<",parent pdgcode="<<parent->GetKF()<<endl;
      Int_t pstat = part->GetKS();// GetStatusCode();//(Pythia8) particle status 
  //    newBranch->Fill();


       Double_t pxch = part->GetPx();
       Double_t pych = part->GetPy();
       Double_t pzch = part->GetPz();
       Double_t ech = part->GetEnergy();
       Double_t ptch=TMath::Sqrt(pxch*pxch+pych*pych);
       Double_t costhetach=pzch/(TMath::Sqrt(ptch*ptch+pzch*pzch));
       ych= TMath::ATan(pzch/ech);
       etach=TMath::ATanH(costhetach);
       if(ptch!=0){
	       /*phichmv=TMath::ATan2(pych,pxch);
                  if (phichmv < 0.0) {phich =phichmv+2*TMath::Pi();}
                  if (phichmv >= 2*TMath::Pi()){phich =phichmv-2*TMath::Pi();}
		  */
	          phichmv=TMath::ACos(pxch/ptch);
                  if(pych<0) phich=2*TMath::Pi()-phichmv;
                  else phich=phichmv;
                  }
       // Filling particles info for jet

       Double_t papx= parent->GetPx();
       Double_t papy= parent->GetPy();
       Double_t papz= parent->GetPz();
       Double_t pae = parent->GetEnergy();
       Double_t zscaling=(ech*ech-(pxch*pxch+pych*pych+pzch*pzch))/(pae*pae-(papx*papx+papy*papy+papz*papz));
      // Double_t zscaling=2*(ech*pae-(pxch*papx+pych*papy+pzch*papz))/(pae*pae-(papx*papx+papy*papy+papz*papz));
       //Double_t zscaling=2*(ech/TMath::Sqrt(pae*pae-(papx*papx+papy*papy+papz*papz)));
       //Double_t zscaling=(ech*pae-(pxch*papx+pych*papy+pzch*papz))/(pae*pae-(papx*papx+papy*papy+papz*papz));

  //  cout<<"starting cuts for particle= "<<ip<<endl;
    if(pdg->GetParticle(idp)->Charge()!=0){
       if(ptch<0.15) continue; //vary cut to 5 GeV see dep on So //newfix
       if(pstat>10)continue;
       hetach->Fill(etach);
       if(TMath::Abs(etach)>0.8) continue;//0.8 //ALICE paper eta 0.9
      // if(pstat!=1) continue; //only status 1-10 (unfragmented or not decayed particle, i.e final state)	    
       nmultch+=1;
       pxA[counter]=pxch;
       pyA[counter]=pych;
       sumapt+=ptch;
       counter++;

  //     jparticles.push_back(PseudoJet(pxch, pych, pzch, ech) );
       invptch->Fill(ptch);
       hych->Fill(ych);
       hphich->Fill(phich);
       if(parent->GetKF()==92){ 
        hzscaff->Fill(zscaling);
	hzscaffcop->Fill(zscaling);
       }
       if(isfast){
	       cout<<"Particle= "<<idp<<"pt="<<ptch<<endl;
               cout<<"parent= "<<part->GetParent()<<",parent pdgcode="<<parent->GetKF()<<endl;
               cout<<"zscaling= "<<zscaling<<endl;
       }
       if(isfast){if(zscaling>1.0) cout<<"greater than 1 zscaling= "<<zscaling<<endl;} //to garanty no z>1


      
       if(ievt==INTEVT){
	      // cout<<"pt="<<ptch<<"pid="<<idp<<"Status: "<<pstat<< endl;
        hyphiN->Fill(etach,phich);
        hyphiwpt->Fill(etach,phich,ptch);
       }

       invptchid[0]->Fill(ptch); //All charged  and stable
       if(parent->GetKF()==92) hzscaffid[0]->Fill(zscaling); 
       if(TMath::Abs(idp)==211 ){ invptchid[1]->Fill(ptch); //pions 211
	   //eventidcount[1]+=1;    
	   ////if(eventidcount[1]=ievt) heventsid->Fill(1);	   

//	   jpions.push_back(PseudoJet(pxch, pych, pzch, ech) );
        if(parent->GetKF()==92) hzscaffid[1]->Fill(zscaling);}
       else if(TMath::Abs(idp)==321){ invptchid[2]->Fill(ptch); //kaons 321
//	   jkaons.push_back(PseudoJet(pxch, pych, pzch, ech) );    
	if(parent->GetKF()==92) hzscaffid[2]->Fill(zscaling);}
       else if(TMath::Abs(idp)==2212 ){ invptchid[3]->Fill(ptch); //protons 2212
//	   jprotons.push_back(PseudoJet(pxch, pych, pzch, ech) );   
	if(parent->GetKF()==92) hzscaffid[3]->Fill(zscaling);}
       else if(TMath::Abs(idp)==3122 ){ invptchid[6]->Fill(ptch);//lambdas 3122
//           jlambdas.push_back(PseudoJet(pxch, pych, pzch, ech) );
        if(parent->GetKF()==92) hzscaffid[6]->Fill(zscaling);}
       else if(TMath::Abs(idp)==3312 ){ invptchid[7]->Fill(ptch); //Xi 3312 
//	   jxis.push_back(PseudoJet(pxch, pych, pzch, ech) );
	if(parent->GetKF()==92) hzscaffid[7]->Fill(zscaling);}
       else if(TMath::Abs(idp)==3334 ){ invptchid[8]->Fill(ptch); //omega 3334
//	   jomegas.push_back(PseudoJet(pxch, pych, pzch, ech) );
	if(parent->GetKF()==92) hzscaffid[8]->Fill(zscaling);}
       else {invptchid[9]->Fill(ptch);
        if(parent->GetKF()==92) hzscaffid[9]->Fill(zscaling);}
     }  //end loop of charged within kin cuts

      if(TMath::Abs(idp)==310 ){ invptchid[4]->Fill(ptch); //K0S 310
//	  jK0s.push_back(PseudoJet(pxch, pych, pzch, ech) );
        if(parent->GetKF()==92) hzscaffid[4]->Fill(zscaling);}
      if(TMath::Abs(idp)==333 ){ invptchid[5]->Fill(ptch); //phi(1020) 333 
//          jphis.push_back(PseudoJet(pxch, pych, pzch, ech) );
        if(parent->GetKF()==92) hzscaffid[5]->Fill(zscaling);}
      if (idp==21 && pstat>0) //detect produced gluons
      if(isfast) cout<<"gluon is produced"<<endl; 
      

    }// end particleloop
    if(isfast) cout<<"out particle loop"<<endl;
    if(multhmax<nmultch){ multhmax=nmultch; evtmulthmax=ievt; }

    //if(nmultch>2)continue;
    if(nmultch>2){
     //cout<<"mult="<<nmultch<<endl;
     Int_t binNch  = -1;
     Int_t BinSot = -1;
     binNch = GetMeasMultiplicityBin(nmultch);
     hNch->Fill(nmultch);

     Double_t Sot = GetSpherocity(pxA, pyA, sumapt, nmultch, size_step); 
     if(ievt==INTEVT)cout<<"interest event with Sot="<<Sot<<endl;
     if(isfast) cout<<"binNch="<<binNch<<", Sot="<<Sot<<endl;
     for( Int_t i_so = 0; i_so < nSobins; ++i_so ){
       if(isESPerc){
	        //cout<<"isESPerc"<<endl;
	        Double_t SotPerc = GetESPerc(Sot,binNch);//,hSOTPerc[binNch]);
		if(isfast) cout<<"SotPerc="<<SotPerc<<endl;
	        if( ( SotPerc >= Sobins[i_so] && SotPerc <= Sobins[i_so+1] ) )  //changeH
                {
		 BinSot = i_so;
	        }
       }else{
	       if( ( Sot >= Sobins[i_so] && Sot <= Sobins[i_so+1] ) )  //changeH
               {
		BinSot = i_so;
	       }
       }
     }

      heventsSo->Fill(BinSot);
      heventsSoN[BinSot]->Fill(binNch);
     if(isfast) cout<<"BinSot="<<BinSot<<endl;
      //separation of zscaling for fragmentation functions
      hzscaffSo[BinSot]->Add(hzscaffcop);
      hzscaffSoNbin[BinSot][binNch]->Add(hzscaffcop);
      for(Int_t i=0; i<hzscaffcop->GetNbinsX(); i++){
         hzscaffcop->SetBinContent(i,0);
      }

     if(isfast) cout<<"zscaffSo y Sonbin added"<<endl;
     if(nmultch==3){if(Sot>0.6) evtmercedes=ievt;}
     if(spheromax<Sot){spheromax=Sot; evtspheromax=ievt;}
     if(spheromin>Sot){spheromin=Sot; evtspheromin=ievt;}
    //  if(isfast)
    //  if(Sot<0.01)
    //  if(Sot>0.90)	
   //   if(Sot>0.80 && nmultch==4)   
//	    cout<<"event # "<<ievt <<", Spherocity="<<Sot<<", Mult="<<nmultch<<endl;    
       if(Sot>0.6 && nmultch==3) 
	    cout<<"Mercedesevent # "<<ievt <<", Spherocity="<<Sot<<", Mult="<<nmultch<<endl; 
      
        if(pxA){// clean up array memory used for TMath::Sort
             delete[] pxA;
             pxA=0;
        }
        if(pyA){// clean up array memory used for TMath::Sort
             delete[] pyA;
             pyA=0;
        }

        if(isnotequalized){
            if(binNch>=0){
               hSOT[binNch]->Fill(Sot);
         //       hSTT[binNch]->Fill(STt);
	       if( BinSot >= 0 && BinSot < nSobins ){
		//		hSoTInNch[binNch][BinSot]->Fill(Sot);
	       }
            }
        }

	
    }//close if nch>2 for So
   if(isfast) cout<<"closed if for Nch>2"<<endl;

    //to find max mult event (3712,274)
    //72
   // if(multhmax==72)cout<<"event #" << ievt <<" with mult max= "<<multhmax<<endl;

/*
    double R= 0.4;
    double ptmin=5.0;//5 GeV ALICE pap
    JetDefinition jet_def(antikt_algorithm,R,pt_scheme,Best);
    if(isfast) cout<<"-----jets--->> "<<endl;
    ClusterSequence cs(jparticles,jet_def);//,pt_scheme);
		    //pt_scheme);
//--- $$ ClusterSequence cspi(jpions,jet_def);
    vector <PseudoJet> jets=sorted_by_pt(cs.inclusive_jets(ptmin));
//--- $$     vector <PseudoJet> jetspi=sorted_by_pt(cspi.inclusive_jets());
    // print out some info
    if(isfast) cout << "Clustered with " <<jet_def.description()<< endl;
    // print the jets
    if(isfast) cout<< " pt, y, phi" <<endl;
    
    for (unsigned ji= 0; ji<jets.size(); ji++) {
      if(isfast) cout<< "jet " <<ji<< ": "<<jets[ji].perp() << " "<<jets[ji].rap() << " " <<jets[ji].phi() <<endl;
      Double_t jrap=jets[ji].rap();
      Double_t jeta=jets[ji].pseudorapidity();
      Double_t jphi=jets[ji].phi();
      Double_t jperp=jets[ji].perp();
      Double_t jpmag=jets[ji].modp2();
      if(maxjetpt<jperp){ maxjetpt=jperp;
          evtmaxjetpt=ievt;
      }
     if(ievt==INTEVT){ //3712 35194){
      //hyphipt->Fill(jrap,jphi,jperp);
      hyphiNfj->Fill(jeta,jphi);
      hyphiwptfj->Fill(jeta,jphi,jperp);
      }
      //hyphip->Fill(jrap,jphi,jpmag);
       hjets->Fill(0);
        if(jperp>=5.0&&jperp<100.0)hjetsptcut->Fill(0);
	if(jperp>=5.0&&jperp<10.)hjetsptcut->Fill(1); 
        if(jperp>=10.0&&jperp<15.)hjetsptcut->Fill(2);
        if(jperp>=15.0&&jperp<20.)hjetsptcut->Fill(3);
        if(jperp>=20.0&&jperp<30.)hjetsptcut->Fill(4);
        if(jperp>=30.0&&jperp<100.)hjetsptcut->Fill(5);
      //  if(jperp>=40.0&&jperp<100.)hjetsptcut->Fill(6);

      if(isfast)cout<<"Filling histo for jet="<<ji<<endl;
      vector <PseudoJet> constituents=jets[ji].constituents();
      for (unsigned jc= 0; jc<constituents.size(); jc++) {
         if(isfast) cout<< "    constituent " <<jc<< "’s pt: "<<constituents[jc].perp()<<" R= "<<constituents[jc].delta_R(jets[ji])<<endl;
      if(ievt==INTEVT){   
 //     hyphiNfj->Fill(constituents[jc].pseudorapidity(),constituents[jc].phi());
 //     hyphiwptfj->Fill(constituents[jc].pseudorapidity(),constituents[jc].phi(),constituents[jc].perp());
      }
       //hyphipt->Fill(constituents[jc].rap(),constituents[jc].phi(),constituents[jc].perp());
       //hyphip->Fill(constituents[jc].rap(),constituents[jc].phi(),constituents[jc].modp2());
       Double_t zchrgd=constituents[jc].perp()/jperp;
   
       //if(zchrgd>=1) cout<<"zchrgd= "<<zchrgd<<endl;
       hzchjetid[0]->Fill(zchrgd);
       if(jperp>=5.0&&jperp<100.0){hzchjet[0]->Fill(zchrgd);}
       if(jperp>=5.0&&jperp<10.){hzchjet[1]->Fill(zchrgd);}
       if(jperp>=10.0&&jperp<15.){hzchjet[2]->Fill(zchrgd);}
       if(jperp>=15.0&&jperp<20.){hzchjet[3]->Fill(zchrgd);}
       if(jperp>=20.0&&jperp<30.){hzchjet[4]->Fill(zchrgd);}
       if(jperp>=30.0&&jperp<40.){hzchjet[5]->Fill(zchrgd);}
      // if(jperp>=40.0&&jperp<100.){hzchjet[6]->Fill(zchrgd);}
       //}
      }
     }
     // to find the max jet event (pt 173) 35194
     //if(maxjetpt==173)cout<<"Event for maxjetpt ="<<ievt <<endl; 

     for (unsigned i= 0; i<jetspi.size(); i++) {
      if(isfast) cout<< "jet " <<i<< ": "<<jetspi[i].perp() << " "<<jetspi[i].rap() << " " <<jetspi[i].phi() <<endl;
      Double_t jrappi=jetspi[i].rap();
      Double_t jphipi=jetspi[i].phi();
      Double_t jperppi=jetspi[i].perp();
      hjets->Fill(1);
      vector <PseudoJet> constpi=jets[i].constituents();
      for (unsigned j= 0; j<constpi.size(); j++) {
       Double_t zchrgdpi=constpi[j].perp()/jperppi;
       hzchjetid[1]->Fill(zchrgdpi);
      }
     }
*/

//*/
    
    // Invariant mass 


    // filling the tree, and the event is over.
    tree->Fill();
  }//end loop of events


   if(isfast){
   cout<<"---------------------------------------"<<endl;
   cout<<"--  no of generated Higgs="<<nHiggs<<endl;
   cout<<"---------------------------------------"<<endl;
   }

   cout<<"event #"<< evtmulthmax <<" with max mult= "<<multhmax<<endl;
   cout<<"event #"<<evtmaxjetpt<<" with max pt jet= " << maxjetpt <<endl;
   cout<<"Isot event #"<<evtspheromax<<" with spheromax= "<<spheromax <<endl;
   cout<<"Jetty event #"<<evtspheromin<<" with spheromin= "<<spheromin <<endl;

   // Show tree structure
  tree->Print();
  

  // some summary transverse momentum  plots:
  
// /* 
  TH1D* hist = new TH1D(HISTNAME, "p_{T}  spectrum for  #mu^{+}",
                        1000, 0, 150);
  hist->SetXTitle("p_{T}[GeV/c]");
  hist->SetYTitle("dN/dp_{T}(GeV^{-1}c)");
  char expression[64];
  sprintf(expression,"sqrt(pow(%s.fPx,2)+pow(%s.fPy,2))>>%s",
          BRANCHNAME, BRANCHNAME, HISTNAME);
  char selection[64];
  sprintf(selection,"%s.fKF==%d", BRANCHNAME, PDGNUMBER); //selection of mu-
  tree->Draw(expression,selection);

  hist->Sumw2();
//*/

//  hist->Fit("expo", "QO+", "", .25, 1.75);
//  TF1* func = hist->GetFunction("expo");
//  func->SetParNames("A", "- 1 / T");
/*
   TF1 *f2 = new TF1("f2","[0]*([1]*(3/4)*(1-x**2)+(3/8)*(1-[1])*(1+x**2))",-1,1);
   TF1 *f3 = new TF1("f3","[0]*( (1/(1+[1]))*(3/4)*(1-x**2) + (3/8)*([1]/(1+[1]))*(1+x**2) )",-1,1);
   TF1 *f4 = new TF1("f4","[0]*( ((-1-[1])/([1]-3)) * (3/4)*(1-x**2) + (3/8)*((2*[1]-2)/([1]-3))*(1+x**2) )",-1,1);
   f2->SetParNames("Normalization","alpha");
   f3->SetParNames("Normalization","R");
   f4->SetParNames("Normalization","tullyR");
   f2->SetParameter(1,0.9);
   f3->SetParameter(1,0.1);
   f4->SetParameter(1,0.8);
   f2->FixParameter(0,Cosine1plot->Integral("bin width"));
   f3->FixParameter(0,Cosine1plot->Integral("bin width"));
   f4->FixParameter(0,Cosine1plot->Integral("bin width"));
   Cosine1plot->Fit("f2","L");
   Cosine1plot->Fit("f3","L+");
   Cosine1plot->Fit("f4","L+");
 */
/*
*/



  // and now we flush and close the file
  file->Write();
  file->Close();
/*  
    if(isnotequalized){
        dataOutFile->cd();
        dataOutFile->Write();
        dataOutFile->Close();
        delete dataOutFile;
        dataOutFile = 0;
    }	
*/
  return 0;
}

// ------ second function to show canvas 
// Show the Pt spectra, and start the tree viewer.
int showEventSample()
{
  // Load needed libraries
  loadLibraries();

   Double_t pi =TMath::Pi();
   Double_t deltaeta=16; //1.6;
   const char* system="e^{+}e^{-}";

  // Open the file
  TFile* file = TFile::Open(FILENAME, "READ");
  if (!file || !file->IsOpen()) {
    //Error("showEventSample", "Couldn't open file %s", FILENAME);
    return 1;
  }

  // Get the tree
  TTree* tree = (TTree*)file->Get(TREENAME);
  if (!tree) {
    //Error("showEventSample", "couldn't get TTree %s", TREENAME);
    return 2;
  }

  // Start the viewer.
  tree->StartViewer();


  const char *outDir="MBfragmentationFPY6";
  TString suffix=IMAGEF;	
  const char *outf=OUTFILEPOSTPRO;
  CreateDir(OUTDIREC);

  TFile* outFile = new TFile(  Form("%s/%s",outDir,outf), "RECREATE");

 //declare postprocessing histograms
  Double_t eta2max=2.;  //2
  Double_t eta2min=-2.; //-2
  Int_t nbineseta2=20; //20
  Double_t deta2=(eta2max-eta2min)/nbineseta2;
  TH1F* hetach2 = new TH1F("hetach2", "Pseudorapidity", nbineseta2, eta2min,eta2max);

//  TH1F* hinvch2=(TH1F*)hinvch->Clone("hinvch2");
	  //new TH1F("hinvch2","Invariant yield",nPtBins,xBins);
		  //nbinespt2,ptch2min,ptch2max);

  // Get the histogram
  TH1D* hist = (TH1D*)file->Get(HISTNAME);
  
  TH1D* hevents=(TH1D*)file->Get("hevents");
  TH1D* heventsSo=(TH1D*)file->Get("heventsSo");
  TH1D* heventsSoN[nSobins];
  for(Int_t i=0; i<nSobins; i++){
    heventsSoN[i]=(TH1D*)file->Get(Form("heventsSoN%d",i));
  }
  TH1D* hjets=(TH1D*)file->Get("hjets");
  TH1D* hjetsptcut=(TH1D*)file->Get("hjetsptcut");
  TH1F* hinvch=(TH1F*)file->Get("invptch");
  TH1F* hinvchid[10];
  TH1F* hzscaffid[10];
  TH1F* hzchjetid[10];
  TH1F* hzchjet[7];
  TH1F* hzscaffSo[nSobins];
  TH1F* hzscaffSoNbin[nSobins][nMultbins];

  TH1F* hSOT[nMultbins];

  for(Int_t i=0; i<10; i++){
    hinvchid[i]=(TH1F*)file->Get(Form("invptchid%d",i));
    hzscaffid[i]=(TH1F*)file->Get(Form("hzscaffid%d",i));
    hzchjetid[i]=(TH1F*)file->Get(Form("hzchjetid%d",i));
  }
  for(Int_t i=0; i<7; i++){
    hzchjet[i]=(TH1F*)file->Get(Form("hzchjet%d",i));
  }
  for(Int_t i=0; i<nSobins; i++){
     hzscaffSo[i]=(TH1F*)file->Get(Form("hzscaffSo%d",i));
     for(Int_t j=0; j<nMultbins; j++){
        hzscaffSoNbin[i][j]=(TH1F*)file->Get(Form("hzscaffSoNbin%d_%d",i,j));
     }
  }

  TH1F* hetach=(TH1F*)file->Get("hetach");
  TH1F* hych=(TH1F*)file->Get("hych");
  
  TH1F* hinvch2=(TH1F*)hinvch->Clone("hinvch2"); 
  TH1F* hinvchid2[10];
  TH1F* hzscaffid2[10];
  TH1F* hzchjetid2[10];
  TH1F* hzchjet2[7];
  TH1F* hzscaffSo2[nSobins];
  TH1F* hzscaffSoNbin2[nSobins][nMultbins];

  for(Int_t i=0; i<10;i++){
    hinvchid2[i]=(TH1F*)hinvchid[i]->Clone(Form("hinvchid%d",i));
    hzscaffid2[i]=(TH1F*)hzscaffid[i]->Clone(Form("hzscaffid%d",i));
    hzchjetid2[i]=(TH1F*)hzchjetid[i]->Clone(Form("hzchjetid%d",i));
  }
  for(Int_t i=0; i<7; i++){
    hzchjet2[i]=(TH1F*)hzchjet[i]->Clone(Form("hzchjet%d",i));
  }
  for(Int_t i=0; i<nSobins; i++){
     hzscaffSo2[i]=(TH1F*)hzscaffSo[i]->Clone(Form("hzscaffSo%d",i));
     for(Int_t j=0; j<nMultbins; j++){
        hzscaffSoNbin2[i][j]=(TH1F*)hzscaffSoNbin[i][j]->Clone(Form("hzscaffSoNbin%d_%d",i,j));
     }
  }


    //TH3F* hyphipt2=(TH3F*)file->Get("hyphipt");
    //TH3F* hyphip2=(TH3F*)file->Get("hyphip");
    TH2F* hyphiN2=(TH2F*)file->Get("hyphiN");
    TH2F* hyphiwpt2=(TH2F*)file->Get("hyphiwpt");
    TH2F* hyphiNfj2=(TH2F*)file->Get("hyphiNfj");
    TH2F* hyphiwptfj2=(TH2F*)file->Get("hyphiwptfj");


  Int_t nevents=hevents->GetBinContent(1);
  Int_t neventsSo[nSobins];
  Int_t neventsSoNch[nSobins][nMultbins];

  for( Int_t i=0; i<nSobins; i++){
    neventsSo[i]=heventsSo->GetBinContent(i+1);
    for(Int_t j=0; j<nMultbins; j++){
       neventsSoNch[i][j]=heventsSoN[i]->GetBinContent(j+1);
     }
  }
  


  Int_t njets=hjets->GetBinContent(1);

  Int_t njetspt[7];
  for(Int_t i=1; i<8; i++){
	  njetspt[i]=hjetsptcut->GetBinContent(i);
  }
  cout<<"number of events="<<nevents<<endl;

  if (!hist) {
    //Error("showEventSample", "couldn't get TH1D %s", HISTNAME);
    return 4;
  }

  //---- Hand rebinning pT -----  
    Int_t nbinvch=hinvch->GetXaxis()->GetNbins();
    Int_t nbinvch2=hinvch2->GetXaxis()->GetNbins();
    Double_t contacpt=0;
    Double_t erracpt=0;
    Int_t binApt=1;
    for(Int_t i=0; i<nbinvch+1; i++){
      Double_t invchcent=hinvch->GetBinCenter(i);
      Double_t invchwidth=hinvch->GetBinWidth(i);
      //simple (same histogram)
      for(Int_t j=0; j<nbinvch2+1; j++){
      Double_t invchcent2=hinvch2->GetBinCenter(j);
      Double_t invchwidth2=hinvch2->GetBinWidth(j);
       if(i==j){
       Double_t cont=hinvch->GetBinContent(i);
       Double_t err=hinvch->GetBinError(i);
    //   cout<<"cont"<<cont<<endl;
       hinvch2->SetBinContent(j,(1/(deltaeta*2*pi*invchcent))*cont/invchwidth2);
       hinvch2->SetBinError(j,(1/(deltaeta*2*pi*invchcent))*err/invchwidth2);
       }
       /*
      if(invchcent>invchcent2-0.5*invchwidth2 && invchcent<invchcent2+0.5*invchwidth2){
     //cout<<"centpt1="<<etacent<<",cent2="<<etacent2<<",width="<<etawidth2<<endl;
      Double_t cont=hinvch->GetBinContent(i);
      //cout<<"j="<<j<<",bingo!....cont="<<cont<<endl;
      Double_t err=hinvch->GetBinError(i);
      if(j==binApt){
        contacpt+=cont;
        erracpt+=err;
        //cout<<"contac"<<contac<<endl;
      hinvch2->SetBinContent(j,contacpt);//invchwidth2);
      hinvch2->SetBinError(j,erracpt);//invchwidth2);
      }else{
        contacpt=cont;
        erracpt=err;
      }
      binApt=j;
    }
*/
    }
   }    

    //Hand rebinning pT (hinvchid[i])
//   /* 
    for(Int_t piden=0; piden<10; piden++){
    Int_t nbinvchid=hinvchid[0]->GetXaxis()->GetNbins();
    Int_t nbinvchid2=hinvchid2[0]->GetXaxis()->GetNbins();
    Double_t contacptid=0;
    Double_t erracptid=0;
    Int_t binAptid=1;
    for(Int_t i=0; i<nbinvchid+1; i++){
      Double_t invchcentid=hinvchid[piden]->GetBinCenter(i);
      Double_t invchwidthid=hinvchid[piden]->GetBinWidth(i);
      for(Int_t j=0; j<nbinvchid2+1; j++){
      Double_t invchidcent2=hinvchid2[piden]->GetBinCenter(j);
      Double_t invchidwidth2=hinvchid2[piden]->GetBinWidth(j);
       if(i==j){
       Double_t contid=hinvchid[piden]->GetBinContent(i);
       Double_t errid=hinvchid[piden]->GetBinError(i);
      // cout<<"cont"<<contid<<endl;
       hinvchid2[piden]->SetBinContent(j,(1.0/(deltaeta*2*pi*invchcentid))*contid/invchidwidth2);
       hinvchid2[piden]->SetBinError(j,(1.0/(deltaeta*2*pi*invchcentid))*errid/invchidwidth2);
       }
      }
     }
    } //end loop of pid


 //Hand rebinning hzscaffid2
//   /* 
    for(Int_t piden=0; piden<10; piden++){
    Int_t nbinzscaid=hzscaffid[0]->GetXaxis()->GetNbins();
    Int_t nbinzscaid2=hzscaffid2[0]->GetXaxis()->GetNbins();
    Double_t contzscaid=0;
    Double_t errzscaid=0;
    for(Int_t i=0; i<nbinzscaid+1; i++){
      Double_t zscaid=hzscaffid[piden]->GetBinCenter(i);
      Double_t zscawidthid=hzscaffid[piden]->GetBinWidth(i);
      for(Int_t j=0; j<nbinzscaid2+1; j++){
      Double_t zscaid2=hzscaffid2[piden]->GetBinCenter(j);
      Double_t zscawidthid2=hzscaffid2[piden]->GetBinWidth(j);
       if(i==j){
       Double_t contzscaid=hzscaffid[piden]->GetBinContent(i);
       Double_t errzscaid=hzscaffid[piden]->GetBinError(i);
      // cout<<"cont"<<contzscaid<<endl;
       hzscaffid2[piden]->SetBinContent(j,contzscaid/zscawidthid2);
       hzscaffid2[piden]->SetBinError(j,errzscaid/zscawidthid2);
       }
      }
     }
    } //end loop of pid


    for(Int_t piden=0; piden<10; piden++){
    Int_t nbinzjid=hzchjetid[0]->GetXaxis()->GetNbins();
    Int_t nbinzjid2=hzchjetid2[0]->GetXaxis()->GetNbins();
    Double_t contzjid=0;
    Double_t errzjid=0;
    for(Int_t i=0; i<nbinzjid+1; i++){
      Double_t zjid=hzchjetid[piden]->GetBinCenter(i);
      Double_t zjwidthid=hzchjetid[piden]->GetBinWidth(i);
      for(Int_t j=0; j<nbinzjid2+1; j++){
      Double_t zjid2=hzchjetid2[piden]->GetBinCenter(j);
      Double_t zjwidthid2=hzchjetid2[piden]->GetBinWidth(j);
       if(i==j){
       Double_t contzjid=hzchjetid[piden]->GetBinContent(i);
       Double_t errzjid=hzchjetid[piden]->GetBinError(i);
      // cout<<"cont"<<contzjid<<endl;
       hzchjetid2[piden]->SetBinContent(j,contzjid/zjwidthid2);
       hzchjetid2[piden]->SetBinError(j,errzjid/zjwidthid2);
       }
      }
     }
    } //end loop of pid

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
    } //end loop of pid
// */

// /* 
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
    } // end loop of So
//*/


  //---- Hand Rebinning eta density------
   Int_t nbinseta=hetach->GetXaxis()->GetNbins();
   Int_t nbinseta2=hetach2->GetXaxis()->GetNbins();
   //cout<<"nbineseta="<<nbinseta<<"nbineseta2="<<nbinseta2<<endl;
   Double_t contac=0;
   Double_t errac=0;
   Int_t binA=1;
  for(Int_t i=0; i<nbinseta+1; i++){
    Double_t etacent=hetach->GetBinCenter(i);
    //cout<<"----i="<<i<<endl;
    if(TMath::Abs(etacent)>8)continue; //heta limit take account
    for(Int_t j=1; j<nbinseta2+1; j++){
     Double_t etacent2=hetach2->GetBinCenter(j);
     Double_t etawidth2=hetach2->GetBinWidth(j);
     if(etacent>etacent2-0.5*etawidth2 && etacent<etacent2+0.5*etawidth2){
     //cout<<"cent1="<<etacent<<",cent2="<<etacent2<<",width="<<etawidth2<<endl;
      Double_t cont=hetach->GetBinContent(i);
      //cout<<"j="<<j<<",bingo!....cont="<<cont<<endl;
      Double_t err=hetach->GetBinError(i);
      if(j==binA){
	contac+=cont;
        errac+=err;
        //cout<<"contac"<<contac<<endl;
      hetach2->SetBinContent(j,contac);
      hetach2->SetBinError(j,errac);
      }else{
        contac=cont;
	errac=err;
      }
      binA=j;
     }
    }
  }


/*
    outFile->cd();
  hetach2->Write();
  hinvch2->Write();
  outFile->Close();
*/

  //--- Draw all the histogram in canvas ------
  gStyle->SetOptStat(1);


 const char * id[10]={"charged |#Delta#eta|<0.8","#pi^{+}+#pi^{-}","K^{+}+K^{-}","p+#bar{p}","K^{0}_{s}","#phi","#Lambda+#bar{#Lambda}","#Xi^{+}+#Xi^{-}","#Omega^{+}+#Omega^{-}","rest"};

   const char *ptint[7]={"p^{ch jet}_{T} 5-100 GeV/c (inclusive)","p^{ch jet}_{T} 5-10 GeV/c","p^{ch jet}_{T} 10-15 GeV/c","p^{ch jet}_{T} 15-20 GeV/c","p^{ch jet}_{T} 20-30 GeV/c","p^{ch jet}_{T} 30-100 GeV/c","p^{ch jet}_{T} 40-100 GeV/c"};


   TCanvas* canvasz = new TCanvas("canvasz","canvas zscal");
  canvasz->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  hzscaffid2[1]->SetTitle("Fragmentation functions for identified particles");
  hzscaffid2[1]->GetXaxis()->SetTitle("z");
  hzscaffid2[1]->GetYaxis()->SetTitle("(1/N_{ev})dN_{ch}/dz");
  hzscaffid2[1]->GetYaxis()->SetRangeUser(0.001,100);
  hzscaffid2[1]->GetXaxis()->SetRangeUser(0.01,1.0);
   hzscaffid2[1]->GetYaxis()->SetLabelSize(0.04);
   hzscaffid2[1]->GetXaxis()->SetLabelSize(0.04);
   hzscaffid2[1]->GetXaxis()->SetTitleSize(0.045);
   hzscaffid2[1]->GetYaxis()->SetTitleSize(0.045);
  for(Int_t i=0; i<10; i++){
   hzscaffid2[i]->Scale(1./(nevents));
   hzscaffid2[i]->SetMarkerStyle(20+i);
  }
  hzscaffid2[0]->SetMarkerColor(kBlack);
  hzscaffid2[1]->SetMarkerColor(kRed);
  hzscaffid2[2]->SetMarkerColor(kGreen);
  hzscaffid2[3]->SetMarkerColor(kBlue);
  hzscaffid2[4]->SetMarkerColor(kViolet);
  hzscaffid2[5]->SetMarkerColor(kPink);
  hzscaffid2[6]->SetMarkerColor(kGray);
  hzscaffid2[7]->SetMarkerColor(kOrange);
  hzscaffid2[8]->SetMarkerColor(kYellow);
  hzscaffid2[9]->SetMarkerColor(kBlue-2);
   TLegend* lzsca = new TLegend(0.1,0.7,0.48,0.9);
   lzsca->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",system,ECM),"C"); //"C" to center header
   lzsca->SetTextFont(42);
   lzsca->SetTextSize(0.042);
   lzsca->SetLineColorAlpha(kWhite,0);
   lzsca->SetLineWidth(0);
   lzsca->SetLineStyle(3);
   lzsca->SetShadowColor(kYellow);
   lzsca->SetFillColorAlpha(kWhite,0.05);
    for(Int_t i=0; i<9; i++){ //9
            lzsca->AddEntry(hzscaffid2[i],Form("%s",id[i]),"lep");
    }
  hzscaffid2[1]->Draw("p");
   hzscaffid2[0]->Draw("psame");
  // hzscaffid2[2]->Draw("psame");
  // hzscaffid2[3]->Draw("psame");
  // hzscaffid2[8]->Draw("psame");
// commetn for only pi k p
 for(Int_t i=2; i<9; i++){
  hzscaffid2[i]->Draw("psame");
  }
  
  lzsca->Draw();


  TCanvas* canvaszjid = new TCanvas("canvaszjid","canvas z for jets");
  canvaszjid->SetLogy();
   gPad->SetTickx();
   gPad->SetTicky();
  hzchjetid2[0]->GetXaxis()->SetTitle("z^{ch}");
  hzchjetid2[0]->GetYaxis()->SetTitle("(1/N_{jets})dN_{ch}/dz^{ch}");
  hzchjetid2[0]->GetYaxis()->SetRangeUser(0.01,100);
  hzchjetid2[0]->GetXaxis()->SetRangeUser(0.01,1.0);
   hzchjetid2[0]->GetYaxis()->SetLabelSize(0.04);
   hzchjetid2[0]->GetXaxis()->SetLabelSize(0.04);
   hzchjetid2[0]->GetXaxis()->SetTitleSize(0.045);
   hzchjetid2[0]->GetYaxis()->SetTitleSize(0.045);
  for(Int_t i=0; i<10; i++){
   hzchjetid2[i]->Scale(1./(njets));
   hzchjetid2[i]->SetMarkerStyle(20+i);
  }
   hzchjetid2[0]->SetMarkerColor(kBlack);
   hzchjetid2[1]->SetMarkerColor(kRed);
   hzchjetid2[2]->SetMarkerColor(kGreen);
   hzchjetid2[3]->SetMarkerColor(kBlue);
   hzchjetid2[4]->SetMarkerColor(kViolet);
   hzchjetid2[5]->SetMarkerColor(kPink);
   hzchjetid2[6]->SetMarkerColor(kGray);
   hzchjetid2[7]->SetMarkerColor(kOrange);
   hzchjetid2[8]->SetMarkerColor(kYellow);
   hzchjetid2[9]->SetMarkerColor(kBlue-2);
   TLegend* lzjid = new TLegend(0.1,0.7,0.48,0.9);
   lzjid->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",system,ECM),"C"); //"C" to center header
   lzjid->SetTextFont(42);
   lzjid->SetTextSize(0.042);
   lzjid->SetLineColorAlpha(kWhite,0);
   lzjid->SetLineWidth(0);
   lzjid->SetLineStyle(3);
   lzjid->SetShadowColor(kYellow);
   lzjid->SetFillColorAlpha(kWhite,0.05);
    for(Int_t i=0; i<10; i++){
            lzjid->AddEntry(hzchjetid2[i],Form("%s",id[i]),"lep");
    }
   hzchjetid2[0]->Draw("p");
  for(Int_t i=1; i<10; i++){
   hzchjetid2[i]->Draw("psame");
  }
  lzjid->Draw();


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
   cout<<"i= "<<i+1<<" njetspt= "<<njetspt[i+1]<<endl;	  
   hzchjet2[i]->Scale(1./(njetspt[i+1]));
   hzchjet2[i]->SetMarkerStyle(21);
  }
 // for(Int_t j=0; j<nSobins; i++){
 //  hzchjet2[j]->SetMarkerColor(100-5*(j-1));
// }
   hzchjet2[0]->SetMarkerColor(kBlack);
   hzchjet2[1]->SetMarkerColor(kRed);
   hzchjet2[2]->SetMarkerColor(kGreen+2);
   hzchjet2[3]->SetMarkerColor(kBlue);
   hzchjet2[4]->SetMarkerColor(kCyan);
   hzchjet2[5]->SetMarkerColor(kPink+10);
   hzchjet2[6]->SetMarkerColor(kYellow+2);
   TLegend* lzj = new TLegend(0.1,0.7,0.48,0.9);
   lzj->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",system,ECM),"C"); //"C" to center header
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

for(Int_t i=0; i<nSobins; i++){
 cout<<"N Spherocity events= "<< neventsSo[i]<<endl;
}

  TCanvas* canvaszSo = new TCanvas("canvaszSo","canvas z for So");
   canvaszSo->SetLogy();
   gPad->SetTickx();
   gPad->SetTicky();
   gStyle->SetOptStat(0);
  hzscaffSo2[0]->SetTitle("Fragmentation functions in S_{0pc} bins");
  hzscaffSo2[0]->GetXaxis()->SetTitle("z");
  hzscaffSo2[0]->GetYaxis()->SetTitle("(1/N_{ev})dN_{ch}/dz");
  hzscaffSo2[0]->GetYaxis()->SetRangeUser(0.001,100);
  hzscaffSo2[0]->GetXaxis()->SetRangeUser(0.001,1.0);
   hzscaffSo2[0]->GetYaxis()->SetLabelSize(0.04);
   hzscaffSo2[0]->GetXaxis()->SetLabelSize(0.04);
   hzscaffSo2[0]->GetXaxis()->SetTitleSize(0.045);
   hzscaffSo2[0]->GetYaxis()->SetTitleSize(0.045);
  for(Int_t i=0; i<nSobins; i++){
   hzscaffSo2[i]->Scale(1./(neventsSo[i])); 
   hzscaffSo2[i]->SetMarkerStyle(20+i);
   hzscaffSo2[i]->SetMarkerColor(99-5*(i));
  }
  /*
   hzscaffSo2[0]->SetMarkerColor(kBlack);
   hzscaffSo2[1]->SetMarkerColor(kRed);
   hzscaffSo2[2]->SetMarkerColor(kGreen);
   hzscaffSo2[3]->SetMarkerColor(kBlue);
   hzscaffSo2[4]->SetMarkerColor(kViolet);
   hzscaffSo2[5]->SetMarkerColor(kPink);
   hzscaffSo2[6]->SetMarkerColor(kGray);
   hzscaffSo2[7]->SetMarkerColor(kOrange);
   hzscaffSo2[8]->SetMarkerColor(kYellow);
   hzscaffSo2[9]->SetMarkerColor(kBlue-2);
   */
   TLegend* lzSo = new TLegend(0.1,0.7,0.48,0.9);
   lzSo->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",system,ECM),"C"); //"C" to center header
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
//----- sovsn  jetty
  
  for(Int_t i=1; i<nMultbins; i++){
   cout<<"events for jetty="<<neventsSoNch[0][i]<<", with"<<Multbins[i]<<"<N<"<<Multbins[i+1]<<endl;
  }

  TCanvas* c3part= new TCanvas("c3part","canvas z for 3 particles structure");
  c3part->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
   gStyle->SetOptStat(0);
  hzscaffSoNbin2[0][1]->SetTitle("3-particle fragmentation function in S_{0pc} bins");
  hzscaffSoNbin2[0][1]->GetXaxis()->SetTitle("z");
  hzscaffSoNbin2[0][1]->GetYaxis()->SetTitle("(1/N_{ev})dN_{ch}/dz");
  hzscaffSoNbin2[0][1]->GetYaxis()->SetRangeUser(0.001,100);
  hzscaffSoNbin2[0][1]->GetXaxis()->SetRangeUser(0.001,1.0);
  hzscaffSoNbin2[0][1]->GetYaxis()->SetLabelSize(0.04);
  hzscaffSoNbin2[0][1]->GetXaxis()->SetLabelSize(0.04);
  hzscaffSoNbin2[0][1]->GetXaxis()->SetTitleSize(0.045);
  hzscaffSoNbin2[0][1]->GetYaxis()->SetTitleSize(0.045);
  for(Int_t j=0; j<nSobins; j++){
   if(neventsSoNch[j][1]>0)hzscaffSoNbin2[j][1]->Scale(1./(neventsSoNch[j][1]));
   hzscaffSoNbin2[j][1]->SetMarkerStyle(20+j);
   hzscaffSoNbin2[j][1]->SetMarkerColor(99-5*j);
  }
  TLegend* lz3pSo = new TLegend(0.1,0.7,0.48,0.9);
   lz3pSo->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",system,ECM),"C"); //"C" to center header
   lz3pSo->SetTextFont(42);
   lz3pSo->SetTextSize(0.042);
   lz3pSo->SetLineColorAlpha(kWhite,0);
   lz3pSo->SetLineWidth(0);
   lz3pSo->SetLineStyle(3);
   lz3pSo->SetShadowColor(kYellow);
   lz3pSo->SetFillColorAlpha(kWhite,0.05);
  for(Int_t i=0; i<nSobins; i++){
            lz3pSo->AddEntry(hzscaffSoNbin2[i][1],Form("%s#leq S_{0}<%s",Sopcin[i],Sopcin[i+1]),"lep");
   }
   hzscaffSoNbin2[0][1]->Draw("p");
  for(Int_t i=1; i<nSobins; i++){
   hzscaffSoNbin2[i][1]->Draw("psame");
  }
  lz3pSo->Draw();


/*
  TCanvas* canvaszSoNbin0 = new TCanvas("canvaszSoNbin0","canvas z for jetty vs Nch");
  canvaszSoNbin0->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  hzscaffSoNbin2[0][2]->GetXaxis()->SetTitle("z");
  hzscaffSoNbin2[0][2]->GetYaxis()->SetTitle("(1/N_{ev})dN_{ch}/dz^{ch}");
  hzscaffSoNbin2[0][2]->GetYaxis()->SetRangeUser(0.000001,100000000);
  hzscaffSoNbin2[0][2]->GetXaxis()->SetRangeUser(0.001,1.0);
   hzscaffSoNbin2[0][2]->GetYaxis()->SetLabelSize(0.04);
   hzscaffSoNbin2[0][2]->GetXaxis()->SetLabelSize(0.04);
   hzscaffSoNbin2[0][2]->GetXaxis()->SetTitleSize(0.045);
   hzscaffSoNbin2[0][2]->GetYaxis()->SetTitleSize(0.045);
  for(Int_t i=2; i<nMultbins; i++){
   if(neventsSoNch[0][i]>0)hzscaffSoNbin2[0][i]->Scale(1./(neventsSoNch[0][i]));
   hzscaffSoNbin2[0][i]->SetMarkerStyle(18+i);
  }
//   hzscaffSoNbin2[0][0]->SetMarkerColor(kBlack);
   hzscaffSoNbin2[0][2]->SetMarkerColor(kRed);
   hzscaffSoNbin2[0][3]->SetMarkerColor(kGreen);
   hzscaffSoNbin2[0][4]->SetMarkerColor(kBlue);
   hzscaffSoNbin2[0][5]->SetMarkerColor(kViolet);
   hzscaffSoNbin2[0][6]->SetMarkerColor(kPink);
   hzscaffSoNbin2[0][7]->SetMarkerColor(kGray);
   hzscaffSoNbin2[0][8]->SetMarkerColor(kOrange);
   hzscaffSoNbin2[0][9]->SetMarkerColor(kYellow);
   hzscaffSoNbin2[0][10]->SetMarkerColor(kBlue-2);
   TLegend* lzSoNbin = new TLegend(0.1,0.7,0.48,0.9);
   lzSoNbin->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",system,ECM),"C"); //"C" to center header
   lzSoNbin->SetTextFont(42);
   lzSoNbin->SetTextSize(0.042);
   lzSoNbin->SetLineColorAlpha(kWhite,0);
   lzSoNbin->SetLineWidth(0);
   lzSoNbin->SetLineStyle(3);
   lzSoNbin->SetShadowColor(kYellow);
   lzSoNbin->SetFillColorAlpha(kWhite,0.05);
    for(Int_t i=2; i<nMultbins; i++){
            lzSoNbin->AddEntry(hzscaffSoNbin2[0][i],Form("%1.0f#leq N_{ch}<%1.0f",Multbins[i],Multbins[i+1]),"lep");
    }
   hzscaffSoNbin2[0][2]->Draw("p");
  for(Int_t i=1; i<nMultbins; i++){
   hzscaffSoNbin2[0][i]->Draw("psame");
  }
  lzSoNbin->Draw();
*/

//---------------------------------
  TCanvas* canvaszSoNbin[10];
  TLegend* lzSoNbins[10];
  for(Int_t j=0; j<10; j++){
	  canvaszSoNbin[j] = new TCanvas(Form("canvaszSoNbin%d",j),Form("canvas z for Sobin=%d vs Nch",j));
	  canvaszSoNbin[j]->SetLogy(); 
	  gPad->SetTickx();
          gPad->SetTicky();
	  gStyle->SetOptStat(0);
 // TCanvas* canvaszSoNbin1 = new TCanvas("canvaszSoNbin1","canvas z for s0bin1 vs Nch");
 // canvaszSoNbin1->SetLogy();
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
//   hzscaffSoNbin2[0][0]->SetMarkerColor(kBlack);
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
   //TLegend* lzSoNbin1 = new TLegend(0.1,0.7,0.48,0.9);
   lzSoNbins[j]->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",system,ECM),"C"); //"C" to center header
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




//-------
///*
  TCanvas* canvas4 = new TCanvas("canvas4","canvas charge");
  canvas4->SetLogy();
  hinvchid2[1]->GetYaxis()->SetTitle("1/(2#pip_{T}N_{ev}) d^{2}N_{ch}/dp_{T}d#eta");
  hinvchid2[1]->GetXaxis()->SetTitle("p_{T} [GeV/c] |#eta<0.8|");
  hinvchid2[1]->GetYaxis()->SetRangeUser(0.000001,100000000);
  for(Int_t i=0; i<10; i++){
   hinvchid2[i]->Scale(1./(nevents));
   hinvchid2[i]->SetMarkerStyle(20+i);
  }
  hinvchid2[0]->SetMarkerColor(kBlack);
  hinvchid2[1]->SetMarkerColor(kRed);
  hinvchid2[2]->SetMarkerColor(kGreen);
  hinvchid2[3]->SetMarkerColor(kBlue);
  hinvchid2[4]->SetMarkerColor(kViolet);
  hinvchid2[5]->SetMarkerColor(kPink);
  hinvchid2[6]->SetMarkerColor(kGray);
  hinvchid2[7]->SetMarkerColor(kOrange);
  hinvchid2[8]->SetMarkerColor(kYellow);
  hinvchid2[9]->SetMarkerColor(kBlue-2);
   TLegend* linvch = new TLegend(0.1,0.7,0.48,0.9);
   linvch->SetHeader(Form("Pythia6, %s at #sqrt{s}=%d GeV",system,ECM),"C"); //"C" to center header
  // linvch->AddEntry(hinvch2,"charged |#Delta#eta| not restricted","lep");
    for(Int_t i=0; i<10; i++){
	    linvch->AddEntry(hinvchid2[i],Form("%s",id[i]),"lep");
    }
  hinvchid2[1]->Draw();
  for(Int_t i=2; i<10; i++){
    hinvchid2[i]->Draw("same");
  }
  linvch->Draw();


//*/

  TCanvas* ceta=new TCanvas("ceta","canvas eta");
  gStyle->SetOptStat(0);
  //ceta->SetLogy();
  hetach2->GetYaxis()->SetTitle("#LT dN/d#eta #GT");
  hetach2->GetXaxis()->SetTitle("#eta");
  hetach2->GetYaxis()->SetRangeUser(4.6,9.0);
  hetach2->SetMarkerStyle(23);
  hetach2->SetMarkerColor(kBlue);
   TLegend* legeta = new TLegend(0.1,0.7,0.48,0.9);
   legeta->SetHeader("pseudorapidity density","C"); //"C" to center header
   legeta->AddEntry(hetach2,Form("Pythia6 %s at #sqrt{s}=%d GeV",system,ECM),"lep");
  hetach2->Scale(1./(deta2*nevents));
  hetach2->Draw();
  legeta->Draw();
//  hetach2->SaveAs(Form("%s/etach.%s",outDir,suffix.Data()));
///*
//
 gStyle->SetPalette(1,0);
/*
 TCanvas* cetaphipt=new TCanvas("cetaphipt","canvas for jets");
 hyphipt2->GetXaxis()->SetTitle("y");
 hyphipt2->GetYaxis()->SetTitle("phi");
 hyphipt2->GetZaxis()->SetTitle("p_{T} (GeV/c)");
// hyphipt2->Draw("lego2");
 hyphipt2->Draw("BOX2Z");
*/

 TCanvas* cetaphiNfj=new TCanvas("cetaphiNfj","canvas for N y phi fj");
 hyphiNfj2->GetXaxis()->SetTitle("#eta");
 hyphiNfj2->GetYaxis()->SetTitle("#phi (rads)");
 hyphiNfj2->GetZaxis()->SetTitle("N_{ch}");
 hyphiNfj2->Draw("lego2z");

 TCanvas* cetaphiwptfj=new TCanvas("cetaphiwptfj","canvas for y phi and pt fj");
 hyphiwptfj2->GetXaxis()->SetTitle("#eta");
 hyphiwptfj2->GetYaxis()->SetTitle("#phi (rads)");
 hyphiwptfj2->GetZaxis()->SetTitle("p_{T} (GeV/c)");
 hyphiwptfj2->Draw("lego2z");


 TCanvas* cetaphiN=new TCanvas("cetaphiN","canvas for N y phi");
 hyphiN2->GetXaxis()->SetTitle("#eta");
 hyphiN2->GetYaxis()->SetTitle("#phi (rads)");
 hyphiN2->GetZaxis()->SetTitle("N_{ch}");
 hyphiN2->Draw("lego2z");

 TCanvas* cetaphiwpt=new TCanvas("cetaphiwpt","canvas for y phi and pt");
 hyphiwpt2->GetXaxis()->SetTitle("#eta");
 hyphiwpt2->GetYaxis()->SetTitle("#phi (rads)");
 hyphiwpt2->GetZaxis()->SetTitle("p_{T} (GeV/c)");
 hyphiwpt2->Draw("lego2z");

/*
 TCanvas* cetaphip=new TCanvas("cetaphip","canvas for jets momentum mag");
 hyphip2->GetXaxis()->SetTitle("y");
 hyphip2->GetYaxis()->SetTitle("phi");
 hyphip2->GetZaxis()->SetTitle("|p| (GeV/c)");
 hyphip2->Draw("BOX2Z");
*/

//*/
//  hist->Draw("e1");
//  TF1* func = hist->GetFunction("expo");

//  char expression[64];
//  sprintf(expression,"T #approx %5.1f", -1000 / func->GetParameter(1));
//  TLatex* latex = new TLatex(1.5, 1e-4, expression);
//  latex->SetTextSize(.1);
//  latex->SetTextColor(4);
//  latex->Draw();

/*
  outFile->cd();
  hetach2->Write();
  hinvch2->Write(); 
  outFile->Close();
*/
  return 0;
}

void pythiaExample(Int_t n=10000) {
   makeEventSample(n);
   showEventSample();
}



//______________________________




#ifndef __CINT__
int main(int argc, char** argv)
{
  TApplication app("app", &argc, argv);

  Int_t n = 100;
  if (argc > 1)
    n = strtol(argv[1], NULL, 0);

  int retVal = 0;
  if (n > 0)
    retVal = makeEventSample(n);
  else {
    retVal = showEventSample();
    app.Run();
  }

  return retVal;
}
#endif

//____________________________________________________________________
//
// EOF
//

