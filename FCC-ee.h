//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 21 01:50:39 2022 by ROOT version 6.24/02
// from TTree events/events
// found on file: mgp8_ee_tbw_FCNC_tua_ecm240.root
//////////////////////////////////////////////////////////

#ifndef FCC-ee_h
#define FCC-ee_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class FCC-ee {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           n_electrons;
   Int_t           n_muons;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selmuons_p;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selmuons_theta;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selmuons_px;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selmuons_py;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selmuons_pz;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selmuons_e;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selmuons_charge;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selectrons_p;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selectrons_theta;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selectrons_px;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selectrons_py;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selectrons_pz;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selectrons_e;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *selectrons_charge;
   vector<vector<int> > *tw_jetconstituents_ee_genkt_ES;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_ES_e;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_ES_px;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_ES_py;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_ES_pz;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_ES_theta;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_ES_flavour;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_ES_btag;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_ES_btag_true;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_ES_ctag;
   vector<vector<int> > *th_jetconstituents_ee_genkt_ES;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_ES_e;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_ES_px;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_ES_py;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_ES_pz;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_ES_theta;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_ES_flavour;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_ES_btag;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_ES_btag_true;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_ES_ctag;
   vector<vector<int> > *tw_jetconstituents_ee_genkt_E0S;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_E0S_e;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_E0S_px;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_E0S_py;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_E0S_pz;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_E0S_theta;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_E0S_flavour;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_E0S_btag;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_E0S_btag_true;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_E0S_ctag;
   vector<vector<int> > *th_jetconstituents_ee_genkt_E0S;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_E0S_e;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_E0S_px;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_E0S_py;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_E0S_pz;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_E0S_theta;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_E0S_flavour;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_E0S_btag;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_E0S_btag_true;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_E0S_ctag;
   vector<vector<int> > *tw_jetconstituents_ee_genkt_PS;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_PS_e;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_PS_px;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_PS_py;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_PS_pz;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *tw_jets_ee_genkt_PS_theta;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_PS_flavour;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_PS_btag;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_PS_btag_true;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *tw_jets_ee_genkt_PS_ctag;
   vector<vector<int> > *th_jetconstituents_ee_genkt_PS;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_PS_e;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_PS_px;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_PS_py;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_PS_pz;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *th_jets_ee_genkt_PS_theta;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_PS_flavour;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_PS_btag;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_PS_btag_true;
   vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int> > *th_jets_ee_genkt_PS_ctag;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *MET;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *MET_x;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *MET_y;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *MET_z;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *MC_PDG;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *MC_status;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *MC_p;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *MC_theta;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *MC_mass;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *MC_e;
   vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float> > *RP_thrustangle;

   // List of branches
   TBranch        *b_n_electrons;   //!
   TBranch        *b_n_muons;   //!
   TBranch        *b_selmuons_p;   //!
   TBranch        *b_selmuons_theta;   //!
   TBranch        *b_selmuons_px;   //!
   TBranch        *b_selmuons_py;   //!
   TBranch        *b_selmuons_pz;   //!
   TBranch        *b_selmuons_e;   //!
   TBranch        *b_selmuons_charge;   //!
   TBranch        *b_selectrons_p;   //!
   TBranch        *b_selectrons_theta;   //!
   TBranch        *b_selectrons_px;   //!
   TBranch        *b_selectrons_py;   //!
   TBranch        *b_selectrons_pz;   //!
   TBranch        *b_selectrons_e;   //!
   TBranch        *b_selectrons_charge;   //!
   TBranch        *b_tw_jetconstituents_ee_genkt_ES;   //!
   TBranch        *b_tw_jets_ee_genkt_ES_e;   //!
   TBranch        *b_tw_jets_ee_genkt_ES_px;   //!
   TBranch        *b_tw_jets_ee_genkt_ES_py;   //!
   TBranch        *b_tw_jets_ee_genkt_ES_pz;   //!
   TBranch        *b_tw_jets_ee_genkt_ES_theta;   //!
   TBranch        *b_tw_jets_ee_genkt_ES_flavour;   //!
   TBranch        *b_tw_jets_ee_genkt_ES_btag;   //!
   TBranch        *b_tw_jets_ee_genkt_ES_btag_true;   //!
   TBranch        *b_tw_jets_ee_genkt_ES_ctag;   //!
   TBranch        *b_th_jetconstituents_ee_genkt_ES;   //!
   TBranch        *b_th_jets_ee_genkt_ES_e;   //!
   TBranch        *b_th_jets_ee_genkt_ES_px;   //!
   TBranch        *b_th_jets_ee_genkt_ES_py;   //!
   TBranch        *b_th_jets_ee_genkt_ES_pz;   //!
   TBranch        *b_th_jets_ee_genkt_ES_theta;   //!
   TBranch        *b_th_jets_ee_genkt_ES_flavour;   //!
   TBranch        *b_th_jets_ee_genkt_ES_btag;   //!
   TBranch        *b_th_jets_ee_genkt_ES_btag_true;   //!
   TBranch        *b_th_jets_ee_genkt_ES_ctag;   //!
   TBranch        *b_tw_jetconstituents_ee_genkt_E0S;   //!
   TBranch        *b_tw_jets_ee_genkt_E0S_e;   //!
   TBranch        *b_tw_jets_ee_genkt_E0S_px;   //!
   TBranch        *b_tw_jets_ee_genkt_E0S_py;   //!
   TBranch        *b_tw_jets_ee_genkt_E0S_pz;   //!
   TBranch        *b_tw_jets_ee_genkt_E0S_theta;   //!
   TBranch        *b_tw_jets_ee_genkt_E0S_flavour;   //!
   TBranch        *b_tw_jets_ee_genkt_E0S_btag;   //!
   TBranch        *b_tw_jets_ee_genkt_E0S_btag_true;   //!
   TBranch        *b_tw_jets_ee_genkt_E0S_ctag;   //!
   TBranch        *b_th_jetconstituents_ee_genkt_E0S;   //!
   TBranch        *b_th_jets_ee_genkt_E0S_e;   //!
   TBranch        *b_th_jets_ee_genkt_E0S_px;   //!
   TBranch        *b_th_jets_ee_genkt_E0S_py;   //!
   TBranch        *b_th_jets_ee_genkt_E0S_pz;   //!
   TBranch        *b_th_jets_ee_genkt_E0S_theta;   //!
   TBranch        *b_th_jets_ee_genkt_E0S_flavour;   //!
   TBranch        *b_th_jets_ee_genkt_E0S_btag;   //!
   TBranch        *b_th_jets_ee_genkt_E0S_btag_true;   //!
   TBranch        *b_th_jets_ee_genkt_E0S_ctag;   //!
   TBranch        *b_tw_jetconstituents_ee_genkt_PS;   //!
   TBranch        *b_tw_jets_ee_genkt_PS_e;   //!
   TBranch        *b_tw_jets_ee_genkt_PS_px;   //!
   TBranch        *b_tw_jets_ee_genkt_PS_py;   //!
   TBranch        *b_tw_jets_ee_genkt_PS_pz;   //!
   TBranch        *b_tw_jets_ee_genkt_PS_theta;   //!
   TBranch        *b_tw_jets_ee_genkt_PS_flavour;   //!
   TBranch        *b_tw_jets_ee_genkt_PS_btag;   //!
   TBranch        *b_tw_jets_ee_genkt_PS_btag_true;   //!
   TBranch        *b_tw_jets_ee_genkt_PS_ctag;   //!
   TBranch        *b_th_jetconstituents_ee_genkt_PS;   //!
   TBranch        *b_th_jets_ee_genkt_PS_e;   //!
   TBranch        *b_th_jets_ee_genkt_PS_px;   //!
   TBranch        *b_th_jets_ee_genkt_PS_py;   //!
   TBranch        *b_th_jets_ee_genkt_PS_pz;   //!
   TBranch        *b_th_jets_ee_genkt_PS_theta;   //!
   TBranch        *b_th_jets_ee_genkt_PS_flavour;   //!
   TBranch        *b_th_jets_ee_genkt_PS_btag;   //!
   TBranch        *b_th_jets_ee_genkt_PS_btag_true;   //!
   TBranch        *b_th_jets_ee_genkt_PS_ctag;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_MET_x;   //!
   TBranch        *b_MET_y;   //!
   TBranch        *b_MET_z;   //!
   TBranch        *b_MC_PDG;   //!
   TBranch        *b_MC_status;   //!
   TBranch        *b_MC_p;   //!
   TBranch        *b_MC_theta;   //!
   TBranch        *b_MC_mass;   //!
   TBranch        *b_MC_e;   //!
   TBranch        *b_RP_thrustangle;   //!

   FCC-ee(TTree *tree=0);
   virtual ~FCC-ee();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef FCC-ee_cxx
FCC-ee::FCC-ee(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mgp8_ee_tbw_FCNC_tua_ecm240.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("mgp8_ee_tbw_FCNC_tua_ecm240.root");
      }
      f->GetObject("events",tree);

   }
   Init(tree);
}

FCC-ee::~FCC-ee()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FCC-ee::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FCC-ee::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void FCC-ee::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   selmuons_p = 0;
   selmuons_theta = 0;
   selmuons_px = 0;
   selmuons_py = 0;
   selmuons_pz = 0;
   selmuons_e = 0;
   selmuons_charge = 0;
   selectrons_p = 0;
   selectrons_theta = 0;
   selectrons_px = 0;
   selectrons_py = 0;
   selectrons_pz = 0;
   selectrons_e = 0;
   selectrons_charge = 0;
   tw_jetconstituents_ee_genkt_ES = 0;
   tw_jets_ee_genkt_ES_e = 0;
   tw_jets_ee_genkt_ES_px = 0;
   tw_jets_ee_genkt_ES_py = 0;
   tw_jets_ee_genkt_ES_pz = 0;
   tw_jets_ee_genkt_ES_theta = 0;
   tw_jets_ee_genkt_ES_flavour = 0;
   tw_jets_ee_genkt_ES_btag = 0;
   tw_jets_ee_genkt_ES_btag_true = 0;
   tw_jets_ee_genkt_ES_ctag = 0;
   th_jetconstituents_ee_genkt_ES = 0;
   th_jets_ee_genkt_ES_e = 0;
   th_jets_ee_genkt_ES_px = 0;
   th_jets_ee_genkt_ES_py = 0;
   th_jets_ee_genkt_ES_pz = 0;
   th_jets_ee_genkt_ES_theta = 0;
   th_jets_ee_genkt_ES_flavour = 0;
   th_jets_ee_genkt_ES_btag = 0;
   th_jets_ee_genkt_ES_btag_true = 0;
   th_jets_ee_genkt_ES_ctag = 0;
   tw_jetconstituents_ee_genkt_E0S = 0;
   tw_jets_ee_genkt_E0S_e = 0;
   tw_jets_ee_genkt_E0S_px = 0;
   tw_jets_ee_genkt_E0S_py = 0;
   tw_jets_ee_genkt_E0S_pz = 0;
   tw_jets_ee_genkt_E0S_theta = 0;
   tw_jets_ee_genkt_E0S_flavour = 0;
   tw_jets_ee_genkt_E0S_btag = 0;
   tw_jets_ee_genkt_E0S_btag_true = 0;
   tw_jets_ee_genkt_E0S_ctag = 0;
   th_jetconstituents_ee_genkt_E0S = 0;
   th_jets_ee_genkt_E0S_e = 0;
   th_jets_ee_genkt_E0S_px = 0;
   th_jets_ee_genkt_E0S_py = 0;
   th_jets_ee_genkt_E0S_pz = 0;
   th_jets_ee_genkt_E0S_theta = 0;
   th_jets_ee_genkt_E0S_flavour = 0;
   th_jets_ee_genkt_E0S_btag = 0;
   th_jets_ee_genkt_E0S_btag_true = 0;
   th_jets_ee_genkt_E0S_ctag = 0;
   tw_jetconstituents_ee_genkt_PS = 0;
   tw_jets_ee_genkt_PS_e = 0;
   tw_jets_ee_genkt_PS_px = 0;
   tw_jets_ee_genkt_PS_py = 0;
   tw_jets_ee_genkt_PS_pz = 0;
   tw_jets_ee_genkt_PS_theta = 0;
   tw_jets_ee_genkt_PS_flavour = 0;
   tw_jets_ee_genkt_PS_btag = 0;
   tw_jets_ee_genkt_PS_btag_true = 0;
   tw_jets_ee_genkt_PS_ctag = 0;
   th_jetconstituents_ee_genkt_PS = 0;
   th_jets_ee_genkt_PS_e = 0;
   th_jets_ee_genkt_PS_px = 0;
   th_jets_ee_genkt_PS_py = 0;
   th_jets_ee_genkt_PS_pz = 0;
   th_jets_ee_genkt_PS_theta = 0;
   th_jets_ee_genkt_PS_flavour = 0;
   th_jets_ee_genkt_PS_btag = 0;
   th_jets_ee_genkt_PS_btag_true = 0;
   th_jets_ee_genkt_PS_ctag = 0;
   MET = 0;
   MET_x = 0;
   MET_y = 0;
   MET_z = 0;
   MC_PDG = 0;
   MC_status = 0;
   MC_p = 0;
   MC_theta = 0;
   MC_mass = 0;
   MC_e = 0;
   RP_thrustangle = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("n_electrons", &n_electrons, &b_n_electrons);
   fChain->SetBranchAddress("n_muons", &n_muons, &b_n_muons);
   fChain->SetBranchAddress("selmuons_p", &selmuons_p, &b_selmuons_p);
   fChain->SetBranchAddress("selmuons_theta", &selmuons_theta, &b_selmuons_theta);
   fChain->SetBranchAddress("selmuons_px", &selmuons_px, &b_selmuons_px);
   fChain->SetBranchAddress("selmuons_py", &selmuons_py, &b_selmuons_py);
   fChain->SetBranchAddress("selmuons_pz", &selmuons_pz, &b_selmuons_pz);
   fChain->SetBranchAddress("selmuons_e", &selmuons_e, &b_selmuons_e);
   fChain->SetBranchAddress("selmuons_charge", &selmuons_charge, &b_selmuons_charge);
   fChain->SetBranchAddress("selectrons_p", &selectrons_p, &b_selectrons_p);
   fChain->SetBranchAddress("selectrons_theta", &selectrons_theta, &b_selectrons_theta);
   fChain->SetBranchAddress("selectrons_px", &selectrons_px, &b_selectrons_px);
   fChain->SetBranchAddress("selectrons_py", &selectrons_py, &b_selectrons_py);
   fChain->SetBranchAddress("selectrons_pz", &selectrons_pz, &b_selectrons_pz);
   fChain->SetBranchAddress("selectrons_e", &selectrons_e, &b_selectrons_e);
   fChain->SetBranchAddress("selectrons_charge", &selectrons_charge, &b_selectrons_charge);
   fChain->SetBranchAddress("tw_jetconstituents_ee_genkt_ES", &tw_jetconstituents_ee_genkt_ES, &b_tw_jetconstituents_ee_genkt_ES);
   fChain->SetBranchAddress("tw_jets_ee_genkt_ES_e", &tw_jets_ee_genkt_ES_e, &b_tw_jets_ee_genkt_ES_e);
   fChain->SetBranchAddress("tw_jets_ee_genkt_ES_px", &tw_jets_ee_genkt_ES_px, &b_tw_jets_ee_genkt_ES_px);
   fChain->SetBranchAddress("tw_jets_ee_genkt_ES_py", &tw_jets_ee_genkt_ES_py, &b_tw_jets_ee_genkt_ES_py);
   fChain->SetBranchAddress("tw_jets_ee_genkt_ES_pz", &tw_jets_ee_genkt_ES_pz, &b_tw_jets_ee_genkt_ES_pz);
   fChain->SetBranchAddress("tw_jets_ee_genkt_ES_theta", &tw_jets_ee_genkt_ES_theta, &b_tw_jets_ee_genkt_ES_theta);
   fChain->SetBranchAddress("tw_jets_ee_genkt_ES_flavour", &tw_jets_ee_genkt_ES_flavour, &b_tw_jets_ee_genkt_ES_flavour);
   fChain->SetBranchAddress("tw_jets_ee_genkt_ES_btag", &tw_jets_ee_genkt_ES_btag, &b_tw_jets_ee_genkt_ES_btag);
   fChain->SetBranchAddress("tw_jets_ee_genkt_ES_btag_true", &tw_jets_ee_genkt_ES_btag_true, &b_tw_jets_ee_genkt_ES_btag_true);
   fChain->SetBranchAddress("tw_jets_ee_genkt_ES_ctag", &tw_jets_ee_genkt_ES_ctag, &b_tw_jets_ee_genkt_ES_ctag);
   fChain->SetBranchAddress("th_jetconstituents_ee_genkt_ES", &th_jetconstituents_ee_genkt_ES, &b_th_jetconstituents_ee_genkt_ES);
   fChain->SetBranchAddress("th_jets_ee_genkt_ES_e", &th_jets_ee_genkt_ES_e, &b_th_jets_ee_genkt_ES_e);
   fChain->SetBranchAddress("th_jets_ee_genkt_ES_px", &th_jets_ee_genkt_ES_px, &b_th_jets_ee_genkt_ES_px);
   fChain->SetBranchAddress("th_jets_ee_genkt_ES_py", &th_jets_ee_genkt_ES_py, &b_th_jets_ee_genkt_ES_py);
   fChain->SetBranchAddress("th_jets_ee_genkt_ES_pz", &th_jets_ee_genkt_ES_pz, &b_th_jets_ee_genkt_ES_pz);
   fChain->SetBranchAddress("th_jets_ee_genkt_ES_theta", &th_jets_ee_genkt_ES_theta, &b_th_jets_ee_genkt_ES_theta);
   fChain->SetBranchAddress("th_jets_ee_genkt_ES_flavour", &th_jets_ee_genkt_ES_flavour, &b_th_jets_ee_genkt_ES_flavour);
   fChain->SetBranchAddress("th_jets_ee_genkt_ES_btag", &th_jets_ee_genkt_ES_btag, &b_th_jets_ee_genkt_ES_btag);
   fChain->SetBranchAddress("th_jets_ee_genkt_ES_btag_true", &th_jets_ee_genkt_ES_btag_true, &b_th_jets_ee_genkt_ES_btag_true);
   fChain->SetBranchAddress("th_jets_ee_genkt_ES_ctag", &th_jets_ee_genkt_ES_ctag, &b_th_jets_ee_genkt_ES_ctag);
   fChain->SetBranchAddress("tw_jetconstituents_ee_genkt_E0S", &tw_jetconstituents_ee_genkt_E0S, &b_tw_jetconstituents_ee_genkt_E0S);
   fChain->SetBranchAddress("tw_jets_ee_genkt_E0S_e", &tw_jets_ee_genkt_E0S_e, &b_tw_jets_ee_genkt_E0S_e);
   fChain->SetBranchAddress("tw_jets_ee_genkt_E0S_px", &tw_jets_ee_genkt_E0S_px, &b_tw_jets_ee_genkt_E0S_px);
   fChain->SetBranchAddress("tw_jets_ee_genkt_E0S_py", &tw_jets_ee_genkt_E0S_py, &b_tw_jets_ee_genkt_E0S_py);
   fChain->SetBranchAddress("tw_jets_ee_genkt_E0S_pz", &tw_jets_ee_genkt_E0S_pz, &b_tw_jets_ee_genkt_E0S_pz);
   fChain->SetBranchAddress("tw_jets_ee_genkt_E0S_theta", &tw_jets_ee_genkt_E0S_theta, &b_tw_jets_ee_genkt_E0S_theta);
   fChain->SetBranchAddress("tw_jets_ee_genkt_E0S_flavour", &tw_jets_ee_genkt_E0S_flavour, &b_tw_jets_ee_genkt_E0S_flavour);
   fChain->SetBranchAddress("tw_jets_ee_genkt_E0S_btag", &tw_jets_ee_genkt_E0S_btag, &b_tw_jets_ee_genkt_E0S_btag);
   fChain->SetBranchAddress("tw_jets_ee_genkt_E0S_btag_true", &tw_jets_ee_genkt_E0S_btag_true, &b_tw_jets_ee_genkt_E0S_btag_true);
   fChain->SetBranchAddress("tw_jets_ee_genkt_E0S_ctag", &tw_jets_ee_genkt_E0S_ctag, &b_tw_jets_ee_genkt_E0S_ctag);
   fChain->SetBranchAddress("th_jetconstituents_ee_genkt_E0S", &th_jetconstituents_ee_genkt_E0S, &b_th_jetconstituents_ee_genkt_E0S);
   fChain->SetBranchAddress("th_jets_ee_genkt_E0S_e", &th_jets_ee_genkt_E0S_e, &b_th_jets_ee_genkt_E0S_e);
   fChain->SetBranchAddress("th_jets_ee_genkt_E0S_px", &th_jets_ee_genkt_E0S_px, &b_th_jets_ee_genkt_E0S_px);
   fChain->SetBranchAddress("th_jets_ee_genkt_E0S_py", &th_jets_ee_genkt_E0S_py, &b_th_jets_ee_genkt_E0S_py);
   fChain->SetBranchAddress("th_jets_ee_genkt_E0S_pz", &th_jets_ee_genkt_E0S_pz, &b_th_jets_ee_genkt_E0S_pz);
   fChain->SetBranchAddress("th_jets_ee_genkt_E0S_theta", &th_jets_ee_genkt_E0S_theta, &b_th_jets_ee_genkt_E0S_theta);
   fChain->SetBranchAddress("th_jets_ee_genkt_E0S_flavour", &th_jets_ee_genkt_E0S_flavour, &b_th_jets_ee_genkt_E0S_flavour);
   fChain->SetBranchAddress("th_jets_ee_genkt_E0S_btag", &th_jets_ee_genkt_E0S_btag, &b_th_jets_ee_genkt_E0S_btag);
   fChain->SetBranchAddress("th_jets_ee_genkt_E0S_btag_true", &th_jets_ee_genkt_E0S_btag_true, &b_th_jets_ee_genkt_E0S_btag_true);
   fChain->SetBranchAddress("th_jets_ee_genkt_E0S_ctag", &th_jets_ee_genkt_E0S_ctag, &b_th_jets_ee_genkt_E0S_ctag);
   fChain->SetBranchAddress("tw_jetconstituents_ee_genkt_PS", &tw_jetconstituents_ee_genkt_PS, &b_tw_jetconstituents_ee_genkt_PS);
   fChain->SetBranchAddress("tw_jets_ee_genkt_PS_e", &tw_jets_ee_genkt_PS_e, &b_tw_jets_ee_genkt_PS_e);
   fChain->SetBranchAddress("tw_jets_ee_genkt_PS_px", &tw_jets_ee_genkt_PS_px, &b_tw_jets_ee_genkt_PS_px);
   fChain->SetBranchAddress("tw_jets_ee_genkt_PS_py", &tw_jets_ee_genkt_PS_py, &b_tw_jets_ee_genkt_PS_py);
   fChain->SetBranchAddress("tw_jets_ee_genkt_PS_pz", &tw_jets_ee_genkt_PS_pz, &b_tw_jets_ee_genkt_PS_pz);
   fChain->SetBranchAddress("tw_jets_ee_genkt_PS_theta", &tw_jets_ee_genkt_PS_theta, &b_tw_jets_ee_genkt_PS_theta);
   fChain->SetBranchAddress("tw_jets_ee_genkt_PS_flavour", &tw_jets_ee_genkt_PS_flavour, &b_tw_jets_ee_genkt_PS_flavour);
   fChain->SetBranchAddress("tw_jets_ee_genkt_PS_btag", &tw_jets_ee_genkt_PS_btag, &b_tw_jets_ee_genkt_PS_btag);
   fChain->SetBranchAddress("tw_jets_ee_genkt_PS_btag_true", &tw_jets_ee_genkt_PS_btag_true, &b_tw_jets_ee_genkt_PS_btag_true);
   fChain->SetBranchAddress("tw_jets_ee_genkt_PS_ctag", &tw_jets_ee_genkt_PS_ctag, &b_tw_jets_ee_genkt_PS_ctag);
   fChain->SetBranchAddress("th_jetconstituents_ee_genkt_PS", &th_jetconstituents_ee_genkt_PS, &b_th_jetconstituents_ee_genkt_PS);
   fChain->SetBranchAddress("th_jets_ee_genkt_PS_e", &th_jets_ee_genkt_PS_e, &b_th_jets_ee_genkt_PS_e);
   fChain->SetBranchAddress("th_jets_ee_genkt_PS_px", &th_jets_ee_genkt_PS_px, &b_th_jets_ee_genkt_PS_px);
   fChain->SetBranchAddress("th_jets_ee_genkt_PS_py", &th_jets_ee_genkt_PS_py, &b_th_jets_ee_genkt_PS_py);
   fChain->SetBranchAddress("th_jets_ee_genkt_PS_pz", &th_jets_ee_genkt_PS_pz, &b_th_jets_ee_genkt_PS_pz);
   fChain->SetBranchAddress("th_jets_ee_genkt_PS_theta", &th_jets_ee_genkt_PS_theta, &b_th_jets_ee_genkt_PS_theta);
   fChain->SetBranchAddress("th_jets_ee_genkt_PS_flavour", &th_jets_ee_genkt_PS_flavour, &b_th_jets_ee_genkt_PS_flavour);
   fChain->SetBranchAddress("th_jets_ee_genkt_PS_btag", &th_jets_ee_genkt_PS_btag, &b_th_jets_ee_genkt_PS_btag);
   fChain->SetBranchAddress("th_jets_ee_genkt_PS_btag_true", &th_jets_ee_genkt_PS_btag_true, &b_th_jets_ee_genkt_PS_btag_true);
   fChain->SetBranchAddress("th_jets_ee_genkt_PS_ctag", &th_jets_ee_genkt_PS_ctag, &b_th_jets_ee_genkt_PS_ctag);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MET_x", &MET_x, &b_MET_x);
   fChain->SetBranchAddress("MET_y", &MET_y, &b_MET_y);
   fChain->SetBranchAddress("MET_z", &MET_z, &b_MET_z);
   fChain->SetBranchAddress("MC_PDG", &MC_PDG, &b_MC_PDG);
   fChain->SetBranchAddress("MC_status", &MC_status, &b_MC_status);
   fChain->SetBranchAddress("MC_p", &MC_p, &b_MC_p);
   fChain->SetBranchAddress("MC_theta", &MC_theta, &b_MC_theta);
   fChain->SetBranchAddress("MC_mass", &MC_mass, &b_MC_mass);
   fChain->SetBranchAddress("MC_e", &MC_e, &b_MC_e);
   fChain->SetBranchAddress("RP_thrustangle", &RP_thrustangle, &b_RP_thrustangle);
   Notify();
}

Bool_t FCC-ee::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FCC-ee::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FCC-ee::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FCC-ee_cxx
