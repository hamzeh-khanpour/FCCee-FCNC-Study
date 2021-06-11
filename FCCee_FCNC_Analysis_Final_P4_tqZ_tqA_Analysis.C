
// Single top quark production as a probe of anomalous tqγ and tqZ couplings at the FCC-ee
// Hamzeh Khanpour, Sara Khatibi, Morteza Khatiri Yanehsari and Mojtaba Mohammadi Najafabadi,
// arXiv:1408.2090 [hep-ph]

/*
To run:
root -l examples/FCCee_FCNC_Analysis_Final_P4_tqZ_tqA_Analysis.C\(\"delphes_file.root\"\)

or Simply execute with: root -l FCCee_FCNC_Analysis_Final_P4_tqZ_tqA_Analysis.C
*/

// ------------------------------

// Stdlib header file for input and output.

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

// include statements for all needed dependencies

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector2.h"
#include "TMath.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TImage.h"
#include "vector"
#include "iomanip"
#include "Math/LorentzVector.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include <TStopwatch.h>
#include <TDatime.h>
#include <TPaveText.h>;


// include statements for Delphes

//#include "external/ExRootAnalysis/ExRootResult.h"        // delphes
//#include "external/ExRootAnalysis/ExRootTreeReader.h"    // delphes
//#include "classes/DelphesClasses.h"                      // delphes
//#include "modules/Delphes.h"                             // delphes

//------------------------------


Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

Float_t deltaPhi( const Float_t phi1, const Float_t phi2 );


void FCCee_FCNC_Analysis_Final_P4_tqZ_tqA_Analysis()
{


  gStyle->SetOptStat(0);   // Dellet Statistical Output in Plots

  TStopwatch Watch;
  Watch.Start();


// Load shared library

    gSystem->Load("/root/ExRootAnalysis/libExRootAnalysis.so");
    gSystem->Load("libPhysics");

    gSystem->Load("libDelphes");

for ( Int_t ifile = 1; ifile < 5; ++ifile )
{

if( ifile == 1 ) {
  const TString intputFile="Delphes_Output-Signal-tqphoton-350GeV.root";

  // Create chain of ROOT trees and read the input file

  TChain chain("Delphes");
  chain.Add(intputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader -> GetEntries();

}

else if( ifile == 2 ) {
  const TString intputFile="Delphes_Output-Background-wwjj-350GeV.root";

  // Create chain of ROOT trees and read the input file

  TChain chain("Delphes");
  chain.Add(intputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader -> GetEntries();

}

else if( ifile == 3 ) {
  const TString intputFile="Delphes_Output-Background-ttbar-350GeV.root";

  // Create chain of ROOT trees and read the input file

  TChain chain("Delphes");
  chain.Add(intputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader -> GetEntries();

}

else if( ifile == 4 ) {
  const TString intputFile="Delphes_Output-Background-zll-350GeV.root";

  // Create chain of ROOT trees and read the input file

  TChain chain("Delphes");
  chain.Add(intputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader -> GetEntries();

}

  // set up branches to read in from file

  TClonesArray *branchMuon = treeReader -> UseBranch("Muon");
  TClonesArray *branchElectron = treeReader -> UseBranch("Electron");
  TClonesArray *branchJet = treeReader -> UseBranch("Jet");
  TClonesArray *branchMissingET = treeReader -> UseBranch("MissingET");


  if (!(branchJet)) {
    cout << "File broken" << endl;
    return;
  }


// Declaring histos and set up output histogram

  TH1D *histbJetPt[4];
  TH1D *histbJetEta[4];

  TH1D *histLeptonPt[4];
  TH1D *histLeptonEta[4];
  TH1D *histLeptonE[4];

  TH1D *histWbosonMass[4];

  TH1D *histTopMass[4];
  TH1D *histTopPt[4];

  TH1D *histDeltaRLeptonbJet[4];

  TH1D *histlightJetE[4];

  TH1D *histbJetMulti[4];


  histbJetPt[ifile]      =  new TH1D( "bJet_PT" , "bJet P_{T}", 80, 0.0, 200.0 );
  histbJetEta[ifile]     =  new TH1D( "bJet_Eta", "bJet Eta"  , 80, -3.50, 3.50 );

  histLeptonPt[ifile]    =  new TH1D( "Lepton_PT" , "Lepton P_{T}", 80, 0.0, 200.0 );
  histLeptonEta[ifile]   =  new TH1D( "Lepton_Eta", "Lepton Eta"  , 80, -3.50, 3.50 );
  histLeptonE[ifile]     =  new TH1D( "Energy of lepton", "Energy of lepton"  , 80, 0.0, 200.0 );

  histTopMass[ifile]     =  new TH1D( "Top_Mass", "Top_Mass"  , 80, 10.0, 400.0 );
  histTopPt[ifile]       =  new TH1D( "Top_PT", "Top_PT", 80, 0.0, 200.0 );

  histDeltaRLeptonbJet[ifile]  =  new TH1D( "DeltaR_Lepton b-Jet", "DeltaR_Lepton b-Jet"  , 80, 0.0, 5.0 );

  histWbosonMass[ifile]   =  new TH1D( "Wboson_Mass", "Wboson_Mass"  , 80, 10.0, 160.0 );

  histlightJetE[ifile]    =  new TH1D( "light Jet E", "light Jet E"  , 80, 0.0, 200.0 );

  histbJetMulti[ifile]    =  new TH1D( "Jet Multiplicity", "Jet Multiplicity"  , 8, 0.0, 8.0 );



// Input Variables to TMVA

  TFile *MyFile;
  TTree *MyTree;

  Float_t  weight, Masstop, Pttop, Etalepton, Elepton, dRleptonbJet, ElightJet;
  Float_t  PtbJet, EtabJet;

// ------------------------------

if( ifile == 1 ) {
    MyFile = new TFile("Branch-Signal-tqphoton-350GeV.root","RECREATE");
    MyTree = new TTree("Signaltqphoton","variables");
// MyTree = new TTree("SignaltqZSigma","variables");
// MyTree = new TTree("SignaltqZLRgamma","variables");

    MyTree -> Branch("weight"  , &weight, "weight/F");
    MyTree -> Branch("Masstop" ,&Masstop, "Masstop/F");
    MyTree -> Branch("Pttop"   ,&Pttop, "Pttop/F");
    MyTree -> Branch("Etalepton",&Etalepton, "Etalepton/F");
    MyTree -> Branch("Elepton"  ,&Elepton, "Elepton/F");
    MyTree -> Branch("dRleptonbJet" ,&dRleptonbJet, "dRleptonbJet/F");
    MyTree -> Branch("PtbJet"  ,&PtbJet, "PtbJet/F");
    MyTree -> Branch("EtabJet" ,&EtabJet, "EtabJet/F");
    MyTree -> Branch("ElightJet" ,&ElightJet, "ElightJet/F");
}

else if( ifile == 2 ) {
    MyFile = new TFile("Branch-Background-wwjj-350GeV.root","RECREATE");
    MyTree = new TTree("Backgroundwwjj","variables");

    MyTree -> Branch("weight"  , &weight, "weight/F");
    MyTree -> Branch("Masstop" ,&Masstop, "Masstop/F");
    MyTree -> Branch("Pttop"   ,&Pttop, "Pttop/F");
    MyTree -> Branch("Etalepton",&Etalepton, "Etalepton/F");
    MyTree -> Branch("Elepton"  ,&Elepton, "Elepton/F");
    MyTree -> Branch("dRleptonbJet" ,&dRleptonbJet, "dRleptonbJet/F");
    MyTree -> Branch("PtbJet"  ,&PtbJet, "PtbJet/F");
    MyTree -> Branch("EtabJet" ,&EtabJet, "EtabJet/F");
    MyTree -> Branch("ElightJet" ,&ElightJet, "ElightJet/F");
}

else if( ifile == 3 ) {
    MyFile = new TFile("Branch-Background-ttbar-350GeV.root","RECREATE");
    MyTree = new TTree("Backgroundttbar","variables");

    MyTree -> Branch("weight"  , &weight, "weight/F");
    MyTree -> Branch("Masstop" ,&Masstop, "Masstop/F");
    MyTree -> Branch("Pttop"   ,&Pttop, "Pttop/F");
    MyTree -> Branch("Etalepton",&Etalepton, "Etalepton/F");
    MyTree -> Branch("Elepton"  ,&Elepton, "Elepton/F");
    MyTree -> Branch("dRleptonbJet" ,&dRleptonbJet, "dRleptonbJet/F");
    MyTree -> Branch("PtbJet"  ,&PtbJet, "PtbJet/F");
    MyTree -> Branch("EtabJet" ,&EtabJet, "EtabJet/F");
    MyTree -> Branch("ElightJet" ,&ElightJet, "ElightJet/F");
}

else if( ifile == 4 ) {
    MyFile = new TFile("Branch-Background-zll-350GeV.root","RECREATE");
    MyTree = new TTree("Backgroundzll","variables");

    MyTree -> Branch("weight"  , &weight, "weight/F");
    MyTree -> Branch("Masstop" ,&Masstop, "Masstop/F");
    MyTree -> Branch("Pttop"   ,&Pttop, "Pttop/F");
    MyTree -> Branch("Etalepton",&Etalepton, "Etalepton/F");
    MyTree -> Branch("Elepton"  ,&Elepton, "Elepton/F");
    MyTree -> Branch("dRleptonbJet" ,&dRleptonbJet, "dRleptonbJet/F");
    MyTree -> Branch("PtbJet"  ,&PtbJet, "PtbJet/F");
    MyTree -> Branch("EtabJet" ,&EtabJet, "EtabJet/F");
    MyTree -> Branch("ElightJet" ,&ElightJet, "ElightJet/F");
}

// ------------------------------

  // Pt and Eta cut parameters

     double ElectronMinPt  = 10.0;
     double ElectronMaxEta = 2.5;
     double MuonMinPt  = 10.0;
     double MuonMaxEta = 2.5;
     double JetMinPt   = 10.0;
     double JetMaxEta  = 2.5;
     double MissingETMinMET = 10.0;

     double deltaRCut = 0.4;

  // counters for b-quark, light Jet and lepton in the events

     Int_t nlepton = 0;
     Int_t nJets = 0;
     Int_t nbtagged = 0;

     
     weight = 1.0;

	// My Good Electrons, Muons and Jets (bJets) which pass the Pt & Eta cuts

    std::vector<Muon*> MyGoodMuon;
    std::vector<Electron*> MyGoodElectron;
    std::vector<MissingET*> MyGoodMissingET;

    std::vector<Jet>  MyPerfectbJet;
    std::vector<Jet>  MyPerfectlightJet;
    std::vector<Jet>  MyPerfectAllJet;


    TLorentzVector MyGoodLepton;
    TLorentzVector MyGoodNeutrino;
    TLorentzVector MyGoodWboson;

    TLorentzVector MyFinalbJet;
    TLorentzVector MyFinallightJet;

    TLorentzVector tmp_top;
    TLorentzVector leptonicTop;

  // set up storage variables
     Int_t iMaxPtJet = -1;
     Int_t bJetIndex = -1;
     Int_t lightJetIndex = -1;

     
   // iEntry is per event
      Long64_t iEntry;
    
  std::cout << "*** Chain Contains " <<  numberOfEntries << "  Events" << std::endl;
      

  for ( iEntry = 0;   iEntry  <  numberOfEntries;  ++iEntry ) // numberOfEntries start event loop
    {
    treeReader -> ReadEntry(iEntry);

      if (iEntry % 50000 == 0)

      std::cout << " --- iEntry =  " << iEntry << std::endl;


    MyGoodMuon.clear();
    MyGoodElectron.clear();
    MyGoodMissingET.clear();

    MyPerfectbJet.clear();
    MyPerfectlightJet.clear();
    MyPerfectAllJet.clear();

	// To check both positive and negative leptons (including electron and muon) pass the Pt and Eta cuts
	bool havePositiveLepton = false;
	bool haveNegativeLepton = false;

// ==================================================================================


     for  (Int_t iMuon = 0; iMuon < branchMuon -> GetEntriesFast(); ++iMuon) {    // reco Muon loop
          Muon* MyMuon = (Muon*) branchMuon -> At(iMuon);

if ( MyMuon->PT  <  MuonMinPt ) { continue; }
if ( fabs(MyMuon->Eta) > MuonMaxEta ) { continue; }

			MyGoodMuon.push_back(  (Muon*) branchMuon -> At(iMuon)  );

				if      (MyMuon -> Charge == 1)
				{
					havePositiveLepton = true;
				}
				else if (MyMuon -> Charge == -1)
				{
					haveNegativeLepton = true;
				}

    } // end Muon loop


// ==================================================================================

    for  (Int_t iElectron = 0; iElectron < branchElectron -> GetEntriesFast(); ++iElectron) {    // reco Electron loop
           Electron * MyElectron = (Electron*) branchElectron -> At(iElectron);

if ( MyElectron->PT  <  ElectronMinPt ) { continue; }
if ( fabs(MyElectron->Eta) > ElectronMaxEta ) { continue; }

			MyGoodElectron.push_back(  (Electron*) branchElectron -> At(iElectron)  );

				if      (MyElectron -> Charge == 1)
				{
					havePositiveLepton = true;
				}
				else if (MyElectron -> Charge == -1)
				{
					haveNegativeLepton = true;
				}

    } // end Electron loop


// ==================================================================================

    // Select the Lepton
    if ( MyGoodMuon.size() <  1  &&  MyGoodElectron.size() <  1  ) { continue; }
    if ( MyGoodMuon.size() >= 1  &&  MyGoodElectron.size() >= 1  ) { continue; }
    if ( MyGoodElectron.size() > 0 ) {

  MyGoodLepton.SetPtEtaPhiE( MyGoodElectron.at(0)->PT, MyGoodElectron.at(0)->Eta, MyGoodElectron.at(0)->Phi, MyGoodElectron.at(0)->P4().E() );


    } else {

  MyGoodLepton.SetPtEtaPhiE( MyGoodMuon.at(0)->PT, MyGoodMuon.at(0)->Eta, MyGoodMuon.at(0)->Phi, MyGoodMuon.at(0)->P4().E() );


    }

// ==================================================================================

    // Get Missing ET

    for  (Int_t iMissingET = 0; iMissingET < branchMissingET -> GetEntriesFast(); ++iMissingET) {    // reco MissingET loop
            MissingET * MyMissingET = (MissingET*) branchMissingET -> At(iMissingET);

         // "MET" Selection
      if ( MyMissingET->MET < MissingETMinMET ) { continue; }

	MyGoodMissingET.push_back(  (MissingET*) branchMissingET -> At(iMissingET)  );

    }  // end MissingET loop

    // Get Missing ET

    if ( MyGoodMissingET.size() < 1 )  { continue; }

    // Neutrino Solution
    MyGoodNeutrino.SetPtEtaPhiE( MyGoodMissingET.at(0)->MET, MyGoodMissingET.at(0)->Eta, MyGoodMissingET.at(0)->Phi, MyGoodMissingET.at(0)->P4().E()  );


         MyGoodWboson =  MyGoodLepton + MyGoodNeutrino;

   nlepton++;  //  lepton selection (counters for lepton candidate)

// ==================================================================================

                //  Reset counter variable
                    iMaxPtJet = -1;

    for  (Int_t iJet = 0; iJet < branchJet -> GetEntriesFast(); ++iJet) {   //   Start reconstructed hadronic jets
	  Jet * MyJet = (Jet*) branchJet -> At(iJet);

   	if ( MyJet->PT < JetMinPt ) { continue; }
   	if ( fabs(MyJet->Eta) > JetMaxEta ) { continue; }

/*
      // for first jet in event
      if (iMaxPtJet == -1) {
	// assume it has the highest pT, store its # from list
	  iMaxPtJet =  iJet;
	// and its object
	  Jet* MyMaxPtJet = (Jet*) branchJet->At(iMaxPtJet);
      }
      // for all other jets
      else {
	// check if this jet's pT is greater than the current max jet pT
	if ( (MyPerfectAllJet[iJet].P4()).Pt() > MyMaxPtJet->PT ) {
	  // if so, store its # from list
	  iMaxPtJet =  iJet;
	  // and its object
	  Jet* MyMaxPtJet = (Jet*) branchJet->At(iMaxPtJet);
	}
      }
*/

		MyPerfectAllJet.push_back(*MyJet);  // MyMaxPtJet

// ----------------------------------------------------

		//  Apply light-tag cut to Jets
		//  Set up storage variables
		    lightJetIndex = -1;

    for (Int_t ilightJet = 0; ilightJet < MyPerfectAllJet.size(); ++ilightJet)  {    //    Start to find light-jet


		if (MyPerfectAllJet[ilightJet].BTag == 0)  // light tag switch
		    {

			MyPerfectlightJet.push_back(MyPerfectAllJet[ilightJet]);

			lightJetIndex = ilightJet;

		    }
		}    //    End light-jet loop

// ----------------------------------------------------

		//  Apply b-tag cut to Jets
		//  Set up storage variables
		    bJetIndex = -1;

    for (Int_t ibJet = 0; ibJet < MyPerfectAllJet.size(); ++ibJet)  {     //    Start to find b-jet


		if (MyPerfectAllJet[ibJet].BTag != 0)  // b tag switch
		    {

		      MyPerfectbJet.push_back(MyPerfectAllJet[ibJet]);

			bJetIndex = ibJet;

		    }
		}    //    End b-jet loop

}   //  End reconstructed jet loop

// ==================================================================================

    if ( MyPerfectAllJet.size()  >= 3 ) { continue; }
    if ( MyPerfectlightJet.size() < 1 ) { continue; }


MyFinallightJet.SetPtEtaPhiE( (MyPerfectlightJet[lightJetIndex].P4()).Pt(), (MyPerfectlightJet[lightJetIndex].P4()).Eta(), (MyPerfectlightJet[lightJetIndex].P4()).Phi(),
 (MyPerfectlightJet[lightJetIndex].P4()).E()  );

   if (  deltaR( (MyPerfectlightJet[lightJetIndex].P4()).Eta(), MyGoodLepton.Eta(), (MyPerfectlightJet[lightJetIndex].P4()).Phi(), MyGoodLepton.Phi() ) < deltaRCut  ) continue;


 nJets++;   // Jets selection (counters for Jets candidate)


  if ( MyPerfectbJet.size()  < 1 ) { continue; }

 // Find the best and highest Pt b-Jet to make the top mass closest to 173.34

double DeltaM  = 1000000.0;
double Ptb = 1.0;

      for ( Int_t i = 0; i <= bJetIndex; ++i ){

      MyFinalbJet.SetPtEtaPhiE( (MyPerfectbJet[i].P4()).Pt(), (MyPerfectbJet[i].P4()).Eta(), (MyPerfectbJet[i].P4()).Phi(), (MyPerfectbJet[i].P4()).E()  );

      if (  deltaR( (MyPerfectbJet[i].P4()).Eta(), MyGoodLepton.Eta(), (MyPerfectbJet[i].P4()).Phi(), MyGoodLepton.Phi() ) < deltaRCut  ) continue;

      tmp_top =  MyGoodWboson + MyFinalbJet;

	  if ( (MyPerfectbJet[i].P4()).Pt() > Ptb ){

	if ( fabs( tmp_top.M() - 173.34 )  <  DeltaM ){

	  DeltaM  =  fabs ( tmp_top.M() - 173.34 );

	  Ptb =  (MyPerfectbJet[i].P4()).Pt();

	  leptonicTop = tmp_top;

	  }
        }
      }


  MyFinalbJet.SetPtEtaPhiE( (MyPerfectbJet[bJetIndex].P4()).Pt(), (MyPerfectbJet[bJetIndex].P4()).Eta(), (MyPerfectbJet[bJetIndex].P4()).Phi(), (MyPerfectbJet[bJetIndex].P4()).E()  );

    if (  deltaR( (MyPerfectbJet[bJetIndex].P4()).Eta(), MyGoodLepton.Eta(), (MyPerfectbJet[bJetIndex].P4()).Phi(), MyGoodLepton.Phi() ) < deltaRCut  ) continue;


	nbtagged++;  // b-jet selection (counters for b-quark candidate)


double luminosity = 100.0;
double WeightS  = 0;
double WeightBwwjj  = 0;
double WeightBttbar = 0;
double WeightBzll  = 0;


double SigmaS = ( 0.017116 * 1000 )*2.0;      // top q A  // We need to Sigma*2.0 for u and c quark
// double SigmaS = ( 0.0096372 * 1000 )*2.0;  // top q Z Sigma mu nu
// double SigmaS = ( 0.0035164 * 1000 )*2.0;  // top q Z Gamma mu
double SigmaBwwjj  = 3.2212  * 1000;
double SigmaBttbar = 0.06255 * 1000;
double SigmaBzll   = 4.0849519 * 1000;

WeightS = 1.0;      // luminosity*SigmaS/numberOfEntries;
WeightBwwjj = 1.0;  // luminosity*SigmaBwwjj/numberOfEntries;
WeightBttbar = 1.0; // luminosity*SigmaBttbar/numberOfEntries;
WeightBzll = 1.0; // luminosity*SigmaBttbar/numberOfEntries;

// Input Variables for Plot


            double TopMass = leptonicTop.M();

if ( 0 < TopMass  &&  TopMass < 500) {

if( ifile == 1 ) {
  histTopMass[ifile] -> Fill( TopMass, WeightS );
                 }
else if( ifile == 2 ) {
  histTopMass[ifile] -> Fill( TopMass, WeightBwwjj );
                      }
else if( ifile == 3 ) {
  histTopMass[ifile] -> Fill( TopMass, WeightBttbar );
                      }
else if( ifile == 4 ) {
  histTopMass[ifile] -> Fill( TopMass, WeightBzll );
                      }
}


	    double TopPt   = leptonicTop.Pt();
	    
                histTopPt[ifile]   -> Fill( TopPt );

	        histWbosonMass[ifile]   -> Fill( MyGoodWboson.M() );

	    double bJetPt  = MyFinalbJet.Pt();
	    double bJetEta = MyFinalbJet.Eta();

	        histbJetPt[ifile]  -> Fill( bJetPt );
	        histbJetEta[ifile] -> Fill( bJetEta );

	    double LeptonPt  = MyGoodLepton.Pt();
	    double LeptonEta = MyGoodLepton.Eta();
	    double LeptonE   = MyGoodLepton.E();
	    double LeptonPhi = MyGoodLepton.Phi();

	        histLeptonPt[ifile]  -> Fill( LeptonPt );
                histLeptonEta[ifile] -> Fill( LeptonEta );
	        histLeptonE[ifile]   -> Fill( LeptonE );

	    double DeltaRLeptonbJet = deltaR( MyFinalbJet.Eta(), MyGoodLepton.Eta(), MyFinalbJet.Phi(), MyGoodLepton.Phi() );

	        histDeltaRLeptonbJet[ifile] -> Fill( DeltaRLeptonbJet );

            double lightJetE  = MyFinallightJet.E();

                histlightJetE[ifile] -> Fill ( lightJetE );


                histbJetMulti[ifile] -> Fill ( MyPerfectbJet.size() );


// ==================================================================================


// Input Variables to TMVA

/*

        Masstop   = leptonicTop.M();
        Pttop     = leptonicTop.Pt();

	    Etalepton = MyGoodLepton.Eta();
	    Elepton   = MyGoodLepton.E();

        dRleptonbJet = deltaR( MyFinalbJet.Eta(), MyGoodLepton.Eta(), MyFinalbJet.Phi(), MyGoodLepton.Phi() );

	    PtbJet  =  MyFinalbJet.Pt();
	    EtabJet =  MyFinalbJet.Eta();

	    ElightJet = MyFinallightJet.E();
*/

// ==================================================================================

    MyTree -> Fill();


    }    // end event loop


    MyFile  -> cd();
    MyTree  -> Write();
    MyFile  -> Close();


cout << " ****************************************************" << endl;

    cout << " Number of lepton candidate (lepton selection) = " <<  nlepton << endl;
        cout << " lepton selection efficiency = " <<  nlepton*1.0/numberOfEntries << endl;
	

    cout << " Number of Jets candidate (Jets selection)  = "  << nJets << endl;
        cout << " Jets selection efficiency  = "  << nJets*1.0/numberOfEntries << endl;
	

    cout << " Number of b-tagged Jets (b-quark candidate) = " <<  nbtagged << endl;
        cout << " b-quark selection efficiency = " <<  nbtagged*1.0/numberOfEntries << endl;


cout << " ****************************************************" << endl;

// back efficiency = all passed bkg events / all input bkg

 }     // End ifile loop


// ==================================================================================


  // initialize a canvas to draw on  // draw histogram on canvas
  
  
    TCanvas *c1 = new TCanvas("c1", "bJet_PT",10, 10, 900, 700);
    c1->cd();
    histbJetPt[1] -> Scale(1/histbJetPt[1] -> Integral());
    histbJetPt[1] -> Draw("HIST");
c1 -> SetTickx(1);
c1 -> SetTicky(1);
    histbJetPt[1] -> SetLineColor(2);
    histbJetPt[1] -> SetLineWidth(2);
    histbJetPt[1] -> SetLineStyle(1);
    histbJetPt[1] -> SetFillColor(3);
    histbJetPt[1] -> SetFillStyle(3001);

    histbJetPt[2] -> Scale(1/histbJetPt[2] -> Integral());
    histbJetPt[2] -> Draw("HISTsames");
    histbJetPt[2] -> SetLineColor(4);
    histbJetPt[2] -> SetLineWidth(2);
    histbJetPt[2] -> SetLineStyle(1);
    histbJetPt[2] -> SetFillColor(5);
    histbJetPt[2] -> SetFillStyle(3002);

    histbJetPt[3] -> Scale(1/histbJetPt[3] -> Integral());
    histbJetPt[3] -> Draw("HISTsames");
    histbJetPt[3] -> SetLineColor(6);
    histbJetPt[3] -> SetLineWidth(2);
    histbJetPt[3] -> SetLineStyle(1);
    histbJetPt[3] -> SetFillColor(7);
    histbJetPt[3] -> SetFillStyle(3017);

    histbJetPt[4] -> Scale(1/histbJetPt[4] -> Integral());
    histbJetPt[4] -> Draw("HISTsames");
    histbJetPt[4] -> SetLineColor(kMagenta+4);
    histbJetPt[4] -> SetLineWidth(2);
    histbJetPt[4] -> SetLineStyle(1);
    histbJetPt[4] -> SetFillColor(9);
    histbJetPt[4] -> SetFillStyle(3018);

      // add axis labels
  histbJetPt[1] -> GetXaxis() -> SetTitle("P_{T}^{b-jet} [GeV]");
  histbJetPt[1] -> GetXaxis() -> SetTitleFont(22);
  histbJetPt[1] -> GetXaxis() -> SetTitleSize(0.04);
  histbJetPt[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histbJetPt[1] -> GetXaxis() -> SetLabelFont(22);
  histbJetPt[1] -> GetXaxis() -> SetLabelSize(0.035);
  histbJetPt[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histbJetPt[1] -> GetYaxis() -> SetTitleFont(22);
  histbJetPt[1] -> GetYaxis() -> SetTitleSize(0.04);
  histbJetPt[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histbJetPt[1] -> GetYaxis() -> SetLabelFont(22);
  histbJetPt[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histbJetPt[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend1=new TLegend(0.65,0.65,0.85,0.85);
   legend1 -> SetFillColor(0);
   legend1 -> SetFillStyle(0);
   legend1 -> SetLineStyle(0);
   legend1 -> SetLineColor(0);
   legend1 -> SetTextFont(42);
   legend1 -> SetTextSize(0.04);
//   legend1->SetHeader("The Legend Title");
   legend1->AddEntry(histbJetPt[1],"tq#gamma","f")->SetTextColor(2);
   legend1->AddEntry(histbJetPt[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend1->AddEntry(histbJetPt[3],"t tbar","f")->SetTextColor(6);
   legend1->AddEntry(histbJetPt[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend1->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c1 -> Print("PTbjet.pdf","pdf");
  c1 -> Print("PTbjet.eps","eps");

// _________________


    TCanvas *c2 = new TCanvas("c2", "bJet_Eta",10, 10, 900, 700);
    c2->cd();
    histbJetEta[1] -> Scale(1/histbJetEta[1] -> Integral());
    histbJetEta[1] -> Draw("HIST");
c2 -> SetTickx(1);
c2 -> SetTicky(1);
    histbJetEta[1] -> SetLineColor(2);
    histbJetEta[1] -> SetLineWidth(2);
    histbJetEta[1] -> SetLineStyle(1);
    histbJetEta[1] -> SetFillColor(3);
    histbJetEta[1] -> SetFillStyle(3001);

    histbJetEta[2] -> Scale(1/histbJetEta[2] -> Integral());
    histbJetEta[2] -> Draw("HISTsames");
    histbJetEta[2] -> SetLineColor(4);
    histbJetEta[2] -> SetLineWidth(2);
    histbJetEta[2] -> SetLineStyle(1);
    histbJetEta[2] -> SetFillColor(5);
    histbJetEta[2] -> SetFillStyle(3002);

    histbJetEta[3] -> Scale(1/histbJetEta[3] -> Integral());
    histbJetEta[3] -> Draw("HISTsames");
    histbJetEta[3] -> SetLineColor(6);
    histbJetEta[3] -> SetLineWidth(2);
    histbJetEta[3] -> SetLineStyle(1);
    histbJetEta[3] -> SetFillColor(7);
    histbJetEta[3] -> SetFillStyle(3017);

    histbJetEta[4] -> Scale(1/histbJetEta[4] -> Integral());
    histbJetEta[4] -> Draw("HISTsames");
    histbJetEta[4] -> SetLineColor(kMagenta+4);
    histbJetEta[4] -> SetLineWidth(2);
    histbJetEta[4] -> SetLineStyle(1);
    histbJetEta[4] -> SetFillColor(9);
    histbJetEta[4] -> SetFillStyle(3018);

      // add axis labels
  histbJetEta[1] -> GetXaxis() -> SetTitle("#eta^{b-jet}");
  histbJetEta[1] -> GetXaxis() -> SetTitleFont(22);
  histbJetEta[1] -> GetXaxis() -> SetTitleSize(0.04);
  histbJetEta[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histbJetEta[1] -> GetXaxis() -> SetLabelFont(22);
  histbJetEta[1] -> GetXaxis() -> SetLabelSize(0.035);
  histbJetEta[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histbJetEta[1] -> GetYaxis() -> SetTitleFont(22);
  histbJetEta[1] -> GetYaxis() -> SetTitleSize(0.04);
  histbJetEta[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histbJetEta[1] -> GetYaxis() -> SetLabelFont(22);
  histbJetEta[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histbJetEta[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend2=new TLegend(0.65,0.65,0.85,0.85);
   legend2 -> SetFillColor(0);
   legend2 -> SetFillStyle(0);
   legend2 -> SetLineStyle(0);
   legend2 -> SetLineColor(0);
   legend2 -> SetTextFont(42);
   legend2 -> SetTextSize(0.04);
//   legend2->SetHeader("The Legend Title");
   legend2->AddEntry(histbJetEta[1],"tq#gamma","f")->SetTextColor(2);
   legend2->AddEntry(histbJetEta[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend2->AddEntry(histbJetEta[3],"t tbar","f")->SetTextColor(6);
   legend2->AddEntry(histbJetEta[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend2->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c2 -> Print("Etabjet.pdf","pdf");
  c2 -> Print("Etabjet.eps","eps");


// _________________

    TCanvas *c3 = new TCanvas("c3", "p_{T} of Lepton",10, 10, 900, 700);
    c3->cd();
    histLeptonPt[1] -> Scale(1/histLeptonPt[1] -> Integral());
    histLeptonPt[1] -> Draw("HIST");
c3 -> SetTickx(1);
c3 -> SetTicky(1);
    histLeptonPt[1] -> SetLineColor(2);
    histLeptonPt[1] -> SetLineWidth(2);
    histLeptonPt[1] -> SetLineStyle(1);
    histLeptonPt[1] -> SetFillColor(3);
    histLeptonPt[1] -> SetFillStyle(3001);

    histLeptonPt[2] -> Scale(1/histLeptonPt[2] -> Integral());
    histLeptonPt[2] -> Draw("HISTsames");
    histLeptonPt[2] -> SetLineColor(4);
    histLeptonPt[2] -> SetLineWidth(2);
    histLeptonPt[2] -> SetLineStyle(1);
    histLeptonPt[2] -> SetFillColor(5);
    histLeptonPt[2] -> SetFillStyle(3002);

    histLeptonPt[3] -> Scale(1/histLeptonPt[3] -> Integral());
    histLeptonPt[3] -> Draw("HISTsames");
    histLeptonPt[3] -> SetLineColor(6);
    histLeptonPt[3] -> SetLineWidth(2);
    histLeptonPt[3] -> SetLineStyle(1);
    histLeptonPt[3] -> SetFillColor(7);
    histLeptonPt[3] -> SetFillStyle(3017);

    histLeptonPt[4] -> Scale(1/histLeptonPt[4] -> Integral());
    histLeptonPt[4] -> Draw("HISTsames");
    histLeptonPt[4] -> SetLineColor(kMagenta+4);
    histLeptonPt[4] -> SetLineWidth(2);
    histLeptonPt[4] -> SetLineStyle(1);
    histLeptonPt[4] -> SetFillColor(9);
    histLeptonPt[4] -> SetFillStyle(3018);

      // add axis labels
  histLeptonPt[1] -> GetXaxis() -> SetTitle("P_{T}^{lepton} [GeV]");
  histLeptonPt[1] -> GetXaxis() -> SetTitleFont(22);
  histLeptonPt[1] -> GetXaxis() -> SetTitleSize(0.04);
  histLeptonPt[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histLeptonPt[1] -> GetXaxis() -> SetLabelFont(22);
  histLeptonPt[1] -> GetXaxis() -> SetLabelSize(0.035);
  histLeptonPt[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histLeptonPt[1] -> GetYaxis() -> SetTitleFont(22);
  histLeptonPt[1] -> GetYaxis() -> SetTitleSize(0.04);
  histLeptonPt[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histLeptonPt[1] -> GetYaxis() -> SetLabelFont(22);
  histLeptonPt[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histLeptonPt[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend3=new TLegend(0.65,0.65,0.85,0.85);
   legend3 -> SetFillColor(0);
   legend3 -> SetFillStyle(0);
   legend3 -> SetLineStyle(0);
   legend3 -> SetLineColor(0);
   legend3 -> SetTextFont(42);
   legend3 -> SetTextSize(0.04);
//   legend3->SetHeader("The Legend Title");
   legend3->AddEntry(histLeptonPt[1],"tq#gamma","f")->SetTextColor(2);
   legend3->AddEntry(histLeptonPt[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend3->AddEntry(histLeptonPt[3],"t tbar","f")->SetTextColor(6);
   legend3->AddEntry(histbJetEta[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend3->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c3 -> Print("LeptonPt.pdf","pdf");
  c3 -> Print("LeptonPt.eps","eps");

  
// _________________

    TCanvas *c4 = new TCanvas("c4", "Lepton_Eta",10, 10, 900, 700);
    c4->cd();
    histLeptonEta[1] -> Scale(1/histLeptonEta[1] -> Integral());
    histLeptonEta[1] -> Draw("HIST");
c4 -> SetTickx(1);
c4 -> SetTicky(1);
    histLeptonEta[1] -> SetLineColor(2);
    histLeptonEta[1] -> SetLineWidth(2);
    histLeptonEta[1] -> SetLineStyle(1);
    histLeptonEta[1] -> SetFillColor(3);
    histLeptonEta[1] -> SetFillStyle(3001);

    histLeptonEta[2] -> Scale(1/histLeptonEta[2] -> Integral());
    histLeptonEta[2] -> Draw("HISTsames");
    histLeptonEta[2] -> SetLineColor(4);
    histLeptonEta[2] -> SetLineWidth(2);
    histLeptonEta[2] -> SetLineStyle(1);
    histLeptonEta[2] -> SetFillColor(5);
    histLeptonEta[2] -> SetFillStyle(3002);

    histLeptonEta[3] -> Scale(1/histLeptonEta[3] -> Integral());
    histLeptonEta[3] -> Draw("HISTsames");
    histLeptonEta[3] -> SetLineColor(6);
    histLeptonEta[3] -> SetLineWidth(2);
    histLeptonEta[3] -> SetLineStyle(1);
    histLeptonEta[3] -> SetFillColor(7);
    histLeptonEta[3] -> SetFillStyle(3017);

    histLeptonEta[4] -> Scale(1/histLeptonEta[4] -> Integral());
    histLeptonEta[4] -> Draw("HISTsames");
    histLeptonEta[4] -> SetLineColor(kMagenta+4);
    histLeptonEta[4] -> SetLineWidth(2);
    histLeptonEta[4] -> SetLineStyle(1);
    histLeptonEta[4] -> SetFillColor(9);
    histLeptonEta[4] -> SetFillStyle(3018);

      // add axis labels
  histLeptonEta[1] -> GetXaxis() -> SetTitle("#eta^{lepton}");
  histLeptonEta[1] -> GetXaxis() -> SetTitleFont(22);
  histLeptonEta[1] -> GetXaxis() -> SetTitleSize(0.04);
  histLeptonEta[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histLeptonEta[1] -> GetXaxis() -> SetLabelFont(22);
  histLeptonEta[1] -> GetXaxis() -> SetLabelSize(0.035);
  histLeptonEta[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histLeptonEta[1] -> GetYaxis() -> SetTitleFont(22);
  histLeptonEta[1] -> GetYaxis() -> SetTitleSize(0.04);
  histLeptonEta[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histLeptonEta[1] -> GetYaxis() -> SetLabelFont(22);
  histLeptonEta[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histLeptonEta[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend4=new TLegend(0.65,0.65,0.85,0.85);
   legend4 -> SetFillColor(0);
   legend4 -> SetFillStyle(0);
   legend4 -> SetLineStyle(0);
   legend4 -> SetLineColor(0);
   legend4 -> SetTextFont(42);
   legend4 -> SetTextSize(0.04);
//   legend4->SetHeader("The Legend Title");
   legend4->AddEntry(histLeptonEta[1],"tq#gamma","f")->SetTextColor(2);
   legend4->AddEntry(histLeptonEta[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend4->AddEntry(histLeptonEta[3],"t tbar","f")->SetTextColor(6);
   legend4->AddEntry(histbJetEta[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend4->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c4 -> Print("LeptonEta.pdf","pdf");
  c4 -> Print("LeptonEta.eps","eps");

// _________________


    TCanvas *c5 = new TCanvas("c5", "Pt of Top",10, 10, 900, 700);
    c5->cd();
    histTopPt[1] -> Scale(1/histTopPt[1] -> Integral());
    histTopPt[1] -> Draw("HIST");
c5 -> SetTickx(1);
c5 -> SetTicky(1);
    histTopPt[1] -> SetLineColor(2);
    histTopPt[1] -> SetLineWidth(2);
    histTopPt[1] -> SetLineStyle(1);
    histTopPt[1] -> SetFillColor(3);
    histTopPt[1] -> SetFillStyle(3001);

    histTopPt[2] -> Scale(1/histTopPt[2] -> Integral());
    histTopPt[2] -> Draw("HISTsames");
    histTopPt[2] -> SetLineColor(4);
    histTopPt[2] -> SetLineWidth(2);
    histTopPt[2] -> SetLineStyle(1);
    histTopPt[2] -> SetFillColor(5);
    histTopPt[2] -> SetFillStyle(3002);

    histTopPt[3] -> Scale(1/histTopPt[3] -> Integral());
    histTopPt[3] -> Draw("HISTsames");
    histTopPt[3] -> SetLineColor(6);
    histTopPt[3] -> SetLineWidth(2);
    histTopPt[3] -> SetLineStyle(1);
    histTopPt[3] -> SetFillColor(7);
    histTopPt[3] -> SetFillStyle(3017);

    histTopPt[4] -> Scale(1/histTopPt[4] -> Integral());
    histTopPt[4] -> Draw("HISTsames");
    histTopPt[4] -> SetLineColor(kMagenta+4);
    histTopPt[4] -> SetLineWidth(2);
    histTopPt[4] -> SetLineStyle(1);
    histTopPt[4] -> SetFillColor(9);
    histTopPt[4] -> SetFillStyle(3018);

      // add axis labels
  histTopPt[1] -> GetXaxis() -> SetTitle("P_{T}^{top} [GeV]");
  histTopPt[1] -> GetXaxis() -> SetTitleFont(22);
  histTopPt[1] -> GetXaxis() -> SetTitleSize(0.04);
  histTopPt[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histTopPt[1] -> GetXaxis() -> SetLabelFont(22);
  histTopPt[1] -> GetXaxis() -> SetLabelSize(0.035);
  histTopPt[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histTopPt[1] -> GetYaxis() -> SetTitleFont(22);
  histTopPt[1] -> GetYaxis() -> SetTitleSize(0.04);
  histTopPt[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histTopPt[1] -> GetYaxis() -> SetLabelFont(22);
  histTopPt[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histTopPt[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend5=new TLegend(0.65,0.65,0.85,0.85);
   legend5 -> SetFillColor(0);
   legend5 -> SetFillStyle(0);
   legend5 -> SetLineStyle(0);
   legend5 -> SetLineColor(0);
   legend5 -> SetTextFont(42);
   legend5 -> SetTextSize(0.04);
//   legend5->SetHeader("The Legend Title");
   legend5->AddEntry(histTopPt[1],"tq#gamma","f")->SetTextColor(2);
   legend5->AddEntry(histTopPt[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend5->AddEntry(histTopPt[3],"t tbar","f")->SetTextColor(6);
   legend5->AddEntry(histbJetEta[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend5->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c5 -> Print("TopPt.pdf","pdf");
  c5 -> Print("TopPt.eps","eps");

 
// _________________


    TCanvas *c6 = new TCanvas("c6", "Energy of Lepton",10, 10, 900, 700);
    c6->cd();
    histLeptonE[1] -> Scale(1/histLeptonE[1] -> Integral());
    histLeptonE[1] -> Draw("HIST");
c6 -> SetTickx(1);
c6 -> SetTicky(1);
    histLeptonE[1] -> SetLineColor(2);
    histLeptonE[1] -> SetLineWidth(2);
    histLeptonE[1] -> SetLineStyle(1);
    histLeptonE[1] -> SetFillColor(3);
    histLeptonE[1] -> SetFillStyle(3001);

    histLeptonE[2] -> Scale(1/histLeptonE[2] -> Integral());
    histLeptonE[2] -> Draw("HISTsames");
    histLeptonE[2] -> SetLineColor(4);
    histLeptonE[2] -> SetLineWidth(2);
    histLeptonE[2] -> SetLineStyle(1);
    histLeptonE[2] -> SetFillColor(5);
    histLeptonE[2] -> SetFillStyle(3002);

    histLeptonE[3] -> Scale(1/histLeptonE[3] -> Integral());
    histLeptonE[3] -> Draw("HISTsames");
    histLeptonE[3] -> SetLineColor(6);
    histLeptonE[3] -> SetLineWidth(2);
    histLeptonE[3] -> SetLineStyle(1);
    histLeptonE[3] -> SetFillColor(7);
    histLeptonE[3] -> SetFillStyle(3017);

    histLeptonE[4] -> Scale(1/histLeptonE[4] -> Integral());
    histLeptonE[4] -> Draw("HISTsames");
    histLeptonE[4] -> SetLineColor(kMagenta+4);
    histLeptonE[4] -> SetLineWidth(2);
    histLeptonE[4] -> SetLineStyle(1);
    histLeptonE[4] -> SetFillColor(9);
    histLeptonE[4] -> SetFillStyle(3018);

      // add axis labels
  histLeptonE[1] -> GetXaxis() -> SetTitle("E^{lepton} [GeV]");
  histLeptonE[1] -> GetXaxis() -> SetTitleFont(22);
  histLeptonE[1] -> GetXaxis() -> SetTitleSize(0.04);
  histLeptonE[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histLeptonE[1] -> GetXaxis() -> SetLabelFont(22);
  histLeptonE[1] -> GetXaxis() -> SetLabelSize(0.035);
  histLeptonE[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histLeptonE[1] -> GetYaxis() -> SetTitleFont(22);
  histLeptonE[1] -> GetYaxis() -> SetTitleSize(0.04);
  histLeptonE[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histLeptonE[1] -> GetYaxis() -> SetLabelFont(22);
  histLeptonE[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histLeptonE[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend6=new TLegend(0.65,0.65,0.85,0.85);
   legend6 -> SetFillColor(0);
   legend6 -> SetFillStyle(0);
   legend6 -> SetLineStyle(0);
   legend6 -> SetLineColor(0);
   legend6 -> SetTextFont(42);
   legend6 -> SetTextSize(0.04);
//   legend6->SetHeader("The Legend Title");
   legend6->AddEntry(histLeptonE[1],"tq#gamma","f")->SetTextColor(2);
   legend6->AddEntry(histLeptonE[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend6->AddEntry(histLeptonE[3],"t tbar","f")->SetTextColor(6);
   legend6->AddEntry(histLeptonE[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend6->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c6 -> Print("LeptonE.pdf","pdf");
  c6 -> Print("LeptonE.eps","eps");

// _________________

    TCanvas *c7 = new TCanvas("c7", "Delta_R Lepton b-Jet",10, 10, 900, 700);
    c7->cd();
    histDeltaRLeptonbJet[1] -> Scale(1/histDeltaRLeptonbJet[1] -> Integral());
    histDeltaRLeptonbJet[1] -> Draw("HIST");
c7 -> SetTickx(1);
c7 -> SetTicky(1);
    histDeltaRLeptonbJet[1] -> SetLineColor(2);
    histDeltaRLeptonbJet[1] -> SetLineWidth(2);
    histDeltaRLeptonbJet[1] -> SetLineStyle(1);
    histDeltaRLeptonbJet[1] -> SetFillColor(3);
    histDeltaRLeptonbJet[1] -> SetFillStyle(3001);

    histDeltaRLeptonbJet[2] -> Scale(1/histDeltaRLeptonbJet[2] -> Integral());
    histDeltaRLeptonbJet[2] -> Draw("HISTsames");
    histDeltaRLeptonbJet[2] -> SetLineColor(4);
    histDeltaRLeptonbJet[2] -> SetLineWidth(2);
    histDeltaRLeptonbJet[2] -> SetLineStyle(1);
    histDeltaRLeptonbJet[2] -> SetFillColor(5);
    histDeltaRLeptonbJet[2] -> SetFillStyle(3002);

    histDeltaRLeptonbJet[3] -> Scale(1/histDeltaRLeptonbJet[3] -> Integral());
    histDeltaRLeptonbJet[3] -> Draw("HISTsames");
    histDeltaRLeptonbJet[3] -> SetLineColor(6);
    histDeltaRLeptonbJet[3] -> SetLineWidth(2);
    histDeltaRLeptonbJet[3] -> SetLineStyle(1);
    histDeltaRLeptonbJet[3] -> SetFillColor(7);
    histDeltaRLeptonbJet[3] -> SetFillStyle(3017);

    histDeltaRLeptonbJet[4] -> Scale(1/histDeltaRLeptonbJet[4] -> Integral());
    histDeltaRLeptonbJet[4] -> Draw("HISTsames");
    histDeltaRLeptonbJet[4] -> SetLineColor(kMagenta+4);
    histDeltaRLeptonbJet[4] -> SetLineWidth(2);
    histDeltaRLeptonbJet[4] -> SetLineStyle(1);
    histDeltaRLeptonbJet[4] -> SetFillColor(9);
    histDeltaRLeptonbJet[4] -> SetFillStyle(3018);

      // add axis labels
  histDeltaRLeptonbJet[1] -> GetXaxis() -> SetTitle("#DeltaR(lepton, b-jet)");
  histDeltaRLeptonbJet[1] -> GetXaxis() -> SetTitleFont(22);
  histDeltaRLeptonbJet[1] -> GetXaxis() -> SetTitleSize(0.04);
  histDeltaRLeptonbJet[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histDeltaRLeptonbJet[1] -> GetXaxis() -> SetLabelFont(22);
  histDeltaRLeptonbJet[1] -> GetXaxis() -> SetLabelSize(0.035);
  histDeltaRLeptonbJet[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histDeltaRLeptonbJet[1] -> GetYaxis() -> SetTitleFont(22);
  histDeltaRLeptonbJet[1] -> GetYaxis() -> SetTitleSize(0.04);
  histDeltaRLeptonbJet[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histDeltaRLeptonbJet[1] -> GetYaxis() -> SetLabelFont(22);
  histDeltaRLeptonbJet[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histDeltaRLeptonbJet[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend7=new TLegend(0.65,0.65,0.85,0.85);
   legend7 -> SetFillColor(0);
   legend7 -> SetFillStyle(0);
   legend7 -> SetLineStyle(0);
   legend7 -> SetLineColor(0);
   legend7 -> SetTextFont(42);
   legend7 -> SetTextSize(0.04);
//   legend7->SetHeader("The Legend Title");
   legend7->AddEntry(histDeltaRLeptonbJet[1],"tq#gamma","f")->SetTextColor(2);
   legend7->AddEntry(histDeltaRLeptonbJet[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend7->AddEntry(histDeltaRLeptonbJet[3],"t tbar","f")->SetTextColor(6);
   legend7->AddEntry(histDeltaRLeptonbJet[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend7->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c7 -> Print("DeltaRLeptonbJet.pdf","pdf");
  c7 -> Print("DeltaRLeptonbJet.eps","eps");

 
// _________________


    TCanvas *c8 = new TCanvas("c8", "Top quark Mass",10, 10, 900, 700);
    c8->cd();
    histTopMass[1] -> Scale(1/histTopMass[1] -> Integral());
    histTopMass[1] -> Draw("HIST");
c8 -> SetTickx(1);
c8 -> SetTicky(1);
    histTopMass[1] -> SetLineColor(2);
    histTopMass[1] -> SetLineWidth(2);
    histTopMass[1] -> SetLineStyle(1);
    histTopMass[1] -> SetFillColor(3);
    histTopMass[1] -> SetFillStyle(3001);

    histTopMass[2] -> Scale(1/histTopMass[2] -> Integral());
    histTopMass[2] -> Draw("HISTsames");
    histTopMass[2] -> SetLineColor(4);
    histTopMass[2] -> SetLineWidth(2);
    histTopMass[2] -> SetLineStyle(1);
    histTopMass[2] -> SetFillColor(5);
    histTopMass[2] -> SetFillStyle(3002);

    histTopMass[3] -> Scale(1/histTopMass[3] -> Integral());
    histTopMass[3] -> Draw("HISTsames");
    histTopMass[3] -> SetLineColor(6);
    histTopMass[3] -> SetLineWidth(2);
    histTopMass[3] -> SetLineStyle(1);
    histTopMass[3] -> SetFillColor(7);
    histTopMass[3] -> SetFillStyle(3017);

    histTopMass[4] -> Scale(1/histTopMass[4] -> Integral());
    histTopMass[4] -> Draw("HISTsames");
    histTopMass[4] -> SetLineColor(kMagenta+4);
    histTopMass[4] -> SetLineWidth(2);
    histTopMass[4] -> SetLineStyle(1);
    histTopMass[4] -> SetFillColor(9);
    histTopMass[4] -> SetFillStyle(3018);

      // add axis labels
  histTopMass[1] -> GetXaxis() -> SetTitle("M_{top}^{reco} [GeV]");
  histTopMass[1] -> GetXaxis() -> SetTitleFont(22);
  histTopMass[1] -> GetXaxis() -> SetTitleSize(0.04);
  histTopMass[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histTopMass[1] -> GetXaxis() -> SetLabelFont(22);
  histTopMass[1] -> GetXaxis() -> SetLabelSize(0.035);
  histTopMass[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histTopMass[1] -> GetYaxis() -> SetTitleFont(22);
  histTopMass[1] -> GetYaxis() -> SetTitleSize(0.04);
  histTopMass[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histTopMass[1] -> GetYaxis() -> SetLabelFont(22);
  histTopMass[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histTopMass[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend8=new TLegend(0.65,0.65,0.85,0.85);
   legend8 -> SetFillColor(0);
   legend8 -> SetFillStyle(0);
   legend8 -> SetLineStyle(0);
   legend8 -> SetLineColor(0);
   legend8 -> SetTextFont(42);
   legend8 -> SetTextSize(0.04);
//   legend8->SetHeader("The Legend Title");
   legend8->AddEntry(histTopMass[1],"tq#gamma","f")->SetTextColor(2);
   legend8->AddEntry(histTopMass[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend8->AddEntry(histTopMass[3],"t tbar","f")->SetTextColor(6);
   legend8->AddEntry(histTopMass[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend8->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c8 -> Print("Mtop.pdf","pdf");
  c8 -> Print("Mtop.eps","eps");

// __________________



    TCanvas *c9 = new TCanvas("c9", "W boson Mass",10, 10, 900, 700);
    c9->cd();
    histWbosonMass[1] -> Scale(1/histWbosonMass[1] -> Integral());
    histWbosonMass[1] -> Draw("HIST");
c9 -> SetTickx(1);
c9 -> SetTicky(1);
    histWbosonMass[1] -> SetLineColor(2);
    histWbosonMass[1] -> SetLineWidth(2);
    histWbosonMass[1] -> SetLineStyle(1);
    histWbosonMass[1] -> SetFillColor(3);
    histWbosonMass[1] -> SetFillStyle(3001);

    histWbosonMass[2] -> Scale(1/histWbosonMass[2] -> Integral());
    histWbosonMass[2] -> Draw("HISTsames");
    histWbosonMass[2] -> SetLineColor(4);
    histWbosonMass[2] -> SetLineWidth(2);
    histWbosonMass[2] -> SetLineStyle(1);
    histWbosonMass[2] -> SetFillColor(5);
    histWbosonMass[2] -> SetFillStyle(3002);

    histWbosonMass[3] -> Scale(1/histWbosonMass[3] -> Integral());
    histWbosonMass[3] -> Draw("HISTsames");
    histWbosonMass[3] -> SetLineColor(6);
    histWbosonMass[3] -> SetLineWidth(2);
    histWbosonMass[3] -> SetLineStyle(1);
    histWbosonMass[3] -> SetFillColor(7);
    histWbosonMass[3] -> SetFillStyle(3017);

    histWbosonMass[4] -> Scale(1/histWbosonMass[4] -> Integral());
    histWbosonMass[4] -> Draw("HISTsames");
    histWbosonMass[4] -> SetLineColor(kMagenta+4);
    histWbosonMass[4] -> SetLineWidth(2);
    histWbosonMass[4] -> SetLineStyle(1);
    histWbosonMass[4] -> SetFillColor(9);
    histWbosonMass[4] -> SetFillStyle(3018);

      // add axis labels
  histWbosonMass[1] -> GetXaxis() -> SetTitle("M_{W}^{reco} [GeV]");
  histWbosonMass[1] -> GetXaxis() -> SetTitleFont(22);
  histWbosonMass[1] -> GetXaxis() -> SetTitleSize(0.04);
  histWbosonMass[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histWbosonMass[1] -> GetXaxis() -> SetLabelFont(22);
  histWbosonMass[1] -> GetXaxis() -> SetLabelSize(0.035);
  histWbosonMass[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histWbosonMass[1] -> GetYaxis() -> SetTitleFont(22);
  histWbosonMass[1] -> GetYaxis() -> SetTitleSize(0.04);
  histWbosonMass[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histWbosonMass[1] -> GetYaxis() -> SetLabelFont(22);
  histWbosonMass[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histWbosonMass[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend9=new TLegend(0.65,0.65,0.85,0.85);
   legend9 -> SetFillColor(0);
   legend9 -> SetFillStyle(0);
   legend9 -> SetLineStyle(0);
   legend9 -> SetLineColor(0);
   legend9 -> SetTextFont(42);
   legend9 -> SetTextSize(0.04);
//   legend->SetHeader("The Legend Title");
   legend9->AddEntry(histWbosonMass[1],"tq#gamma","f")->SetTextColor(2);
   legend9->AddEntry(histWbosonMass[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend9->AddEntry(histWbosonMass[3],"t tbar","f")->SetTextColor(6);
   legend9->AddEntry(histbJetEta[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend9->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c9 -> Print("WbosonMass.pdf","pdf");
  c9 -> Print("WbosonMass.eps","eps");

// __________________


    TCanvas *c10 = new TCanvas("c10", "Eenrgy of light Jet",10, 10, 900, 700);
    c10->cd();
    histlightJetE[1] -> Scale(1/histlightJetE[1] -> Integral());
    histlightJetE[1] -> Draw("HIST");
c10 -> SetTickx(1);
c10 -> SetTicky(1);
    histlightJetE[1] -> SetLineColor(2);
    histlightJetE[1] -> SetLineWidth(2);
    histlightJetE[1] -> SetLineStyle(1);
    histlightJetE[1] -> SetFillColor(3);
    histlightJetE[1] -> SetFillStyle(3001);

    histlightJetE[2] -> Scale(1/histlightJetE[2] -> Integral());
    histlightJetE[2] -> Draw("HISTsames");
    histlightJetE[2] -> SetLineColor(4);
    histlightJetE[2] -> SetLineWidth(2);
    histlightJetE[2] -> SetLineStyle(1);
    histlightJetE[2] -> SetFillColor(5);
    histlightJetE[2] -> SetFillStyle(3002);

    histlightJetE[3] -> Scale(1/histlightJetE[3] -> Integral());
    histlightJetE[3] -> Draw("HISTsames");
    histlightJetE[3] -> SetLineColor(6);
    histlightJetE[3] -> SetLineWidth(2);
    histlightJetE[3] -> SetLineStyle(1);
    histlightJetE[3] -> SetFillColor(7);
    histlightJetE[3] -> SetFillStyle(3017);

    histlightJetE[4] -> Scale(1/histlightJetE[4] -> Integral());
    histlightJetE[4] -> Draw("HISTsames");
    histlightJetE[4] -> SetLineColor(kMagenta+4);
    histlightJetE[4] -> SetLineWidth(2);
    histlightJetE[4] -> SetLineStyle(1);
    histlightJetE[4] -> SetFillColor(9);
    histlightJetE[4] -> SetFillStyle(3018);

      // add axis labels
  histlightJetE[1] -> GetXaxis() -> SetTitle("E^{light Jet} [GeV]");
  histlightJetE[1] -> GetXaxis() -> SetTitleFont(22);
  histlightJetE[1] -> GetXaxis() -> SetTitleSize(0.04);
  histlightJetE[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histlightJetE[1] -> GetXaxis() -> SetLabelFont(22);
  histlightJetE[1] -> GetXaxis() -> SetLabelSize(0.035);
  histlightJetE[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histlightJetE[1] -> GetYaxis() -> SetTitleFont(22);
  histlightJetE[1] -> GetYaxis() -> SetTitleSize(0.04);
  histlightJetE[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histlightJetE[1] -> GetYaxis() -> SetLabelFont(22);
  histlightJetE[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histlightJetE[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend10=new TLegend(0.65,0.65,0.85,0.85);
   legend10 -> SetFillColor(0);
   legend10 -> SetFillStyle(0);
   legend10 -> SetLineStyle(0);
   legend10 -> SetLineColor(0);
   legend10 -> SetTextFont(42);
   legend10 -> SetTextSize(0.04);
//   legend->SetHeader("The Legend Title");
   legend10->AddEntry(histlightJetE[1],"tq#gamma","f")->SetTextColor(2);
   legend10->AddEntry(histlightJetE[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend10->AddEntry(histlightJetE[3],"t tbar","f")->SetTextColor(6);
   legend10->AddEntry(histbJetEta[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend10->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c10 -> Print("lightJetE.pdf","pdf");
  c10 -> Print("lightJetE.eps","eps");


 
// __________________


    TCanvas *c11 = new TCanvas("c11", "b-Jet Multiplicity",10, 10, 900, 700);
    c11->cd();
    histbJetMulti[1] -> Scale(1/histbJetMulti[1] -> Integral());
    histbJetMulti[1] -> Draw("HIST");
c11 -> SetTickx(1);
c11 -> SetTicky(1);
    histbJetMulti[1] -> SetLineColor(2);
    histbJetMulti[1] -> SetLineWidth(2);
    histbJetMulti[1] -> SetLineStyle(1);
    histbJetMulti[1] -> SetFillColor(3);
    histbJetMulti[1] -> SetFillStyle(3001);

    histbJetMulti[2] -> Scale(1/histbJetMulti[2] -> Integral());
    histbJetMulti[2] -> Draw("HISTsames");
    histbJetMulti[2] -> SetLineColor(4);
    histbJetMulti[2] -> SetLineWidth(2);
    histbJetMulti[2] -> SetLineStyle(1);
    histbJetMulti[2] -> SetFillColor(5);
    histbJetMulti[2] -> SetFillStyle(3002);

    histbJetMulti[3] -> Scale(1/histbJetMulti[3] -> Integral());
    histbJetMulti[3] -> Draw("HISTsames");
    histbJetMulti[3] -> SetLineColor(6);
    histbJetMulti[3] -> SetLineWidth(2);
    histbJetMulti[3] -> SetLineStyle(1);
    histbJetMulti[3] -> SetFillColor(7);
    histbJetMulti[3] -> SetFillStyle(3017);

    histbJetMulti[4] -> Scale(1/histbJetMulti[4] -> Integral());
    histbJetMulti[4] -> Draw("HISTsames");
    histbJetMulti[4] -> SetLineColor(kMagenta+4);
    histbJetMulti[4] -> SetLineWidth(2);
    histbJetMulti[4] -> SetLineStyle(1);
    histbJetMulti[4] -> SetFillColor(9);
    histbJetMulti[4] -> SetFillStyle(3018);

      // add axis labels
  histbJetMulti[1] -> GetXaxis() -> SetTitle("b-jet Multiplicity");
  histbJetMulti[1] -> GetXaxis() -> SetTitleFont(22);
  histbJetMulti[1] -> GetXaxis() -> SetTitleSize(0.04);
  histbJetMulti[1] -> GetXaxis() -> SetTitleOffset(1.05);
  histbJetMulti[1] -> GetXaxis() -> SetLabelFont(22);
  histbJetMulti[1] -> GetXaxis() -> SetLabelSize(0.035);
  histbJetMulti[1] -> GetYaxis() -> SetTitle("Normalized distribution");
  histbJetMulti[1] -> GetYaxis() -> SetTitleFont(22);
  histbJetMulti[1] -> GetYaxis() -> SetTitleSize(0.04);
  histbJetMulti[1] -> GetYaxis() -> SetTitleOffset(1.20);
  histbJetMulti[1] -> GetYaxis() -> SetLabelFont(22);
  histbJetMulti[1] -> GetYaxis() -> SetLabelSize(0.035);
//  histbJetMulti[1] -> SetTitle( "FCC-ee/TLEP," "" "#sqrt{S} = 350 GeV," "" "tq#gamma" ); // title on top FCCee FCNC Analysis


   // Draw the legend
   TLegend *legend11=new TLegend(0.65,0.65,0.85,0.85);
   legend11 -> SetFillColor(0);
   legend11 -> SetFillStyle(0);
   legend11 -> SetLineStyle(0);
   legend11 -> SetLineColor(0);
   legend11 -> SetTextFont(42);
   legend11 -> SetTextSize(0.04);
//   legend->SetHeader("The Legend Title");
   legend11->AddEntry(histbJetMulti[1],"tq#gamma","f")->SetTextColor(2);
   legend11->AddEntry(histbJetMulti[2],"W^{#pm}jj","f")->SetTextColor(4);
   legend11->AddEntry(histbJetMulti[3],"t tbar","f")->SetTextColor(6);
   legend11->AddEntry(histbJetMulti[4],"z ll","f")->SetTextColor(kMagenta+4);
   legend11->Draw();


   TPaveText * pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
   TLatex *   tex = new TLatex(50.0,0.17,"FCC-ee");
   tex->SetTextSize(0.04);
   tex->Draw();

   pt = new TPaveText(0.2,0.9,0.8,0.9,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   pt->Draw();
              tex = new TLatex(35.0,0.15,"#sqrt{S} = 350 GeV");
   tex->SetTextSize(0.04);
   tex->Draw();

  // outputs to PDF and postscript file

  c10 -> Print("bJetMulti.pdf","pdf");
  c10 -> Print("bJetMulti.eps","eps");

    
// __________________


  cout << endl;
  cout << "The End of the FCC-ee FCNC Analysis" << endl;



   	Int_t Min = Watch.RealTime()/60.0;

   	double Sec = Watch.RealTime() - 60.0*Min;

   	std::cout << "The Running Time Was: " << Min << "--- Minute and " << Sec << "--- Second" << std::endl;



//  delete [] histbJetPt;

}     // End void FCCee_FCNC_Analysis_Final()


// __________________


Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t Pi = 3.14159265358979;

  Float_t etaDiff = ( eta1 - eta2 );
  Float_t phiDiff = fabs( phi1 - phi2 );
  while ( phiDiff > Pi ) phiDiff = fabs( phiDiff - 2.0 * Pi );

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}

Float_t deltaPhi( const Float_t phi1, const Float_t phi2 ) {

  const Float_t Pi = 3.14159265358979;

  Float_t phiDiff = fabs( phi1 - phi2 );
  while ( phiDiff > Pi ) phiDiff = fabs(phiDiff - 2.0 * Pi);

  return phiDiff;

}
