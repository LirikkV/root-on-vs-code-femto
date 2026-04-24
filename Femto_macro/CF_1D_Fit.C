#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"

#include <iostream>
using namespace::std;

const Int_t N_hist_types_1D = 2; //A or B
const Int_t N_Charge = 2; // 0 <-> Pi_Plus; 1 <-> Pi_Minus
const Int_t N_Bins_Kt = 4;
const Int_t N_Bins_Centr = 9; 

const Double_t KtBins[N_Bins_Kt+1] = {0.15, 0.25, 0.35, 0.45, 0.60};

const TString Input_File = "/home/kirill/root-on-vs-code/Femto_output/23_04_FMRW/out_Au_Au200_23_04_FMR_proseed_1D_integral_and_cuts.root";
const TString Output_File = "/home/kirill/root-on-vs-code/Femto_output/23_04_FMRW/out_Au_Au200_23_04_FMR_proseed_1D_integral_and_cuts_Fit.root";
const TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_output/23_04_FMRW/";

void CF_1D_Fit()
{
TFile *f = TFile::Open(Input_File, "READ");
//1D total, no cuts:

TH1D *hCF_Pi_Plus = (TH1D*)f->Get("CF_Pi_Plus");
TH1D *hCF_Pi_Minus = (TH1D*)f->Get("CF_Pi_Minus");
//
hCF_Pi_Plus->Sumw2();
hCF_Pi_Plus->Sumw2();

//Let's fit CF:
const Double_t x1_fit_range = 0.00;
const Double_t x2_fit_range = 0.20;
const Double_t N_start = 5.0;
const Double_t Lambda_start = 0.6;
const Double_t R_inv_start = 5.0;


hCF_Pi_Plus->GetXaxis()->SetRangeUser(x1_fit_range,x2_fit_range);

TF1 *f1 = new TF1("f1","[0]*(1+[1]*exp(-[2]*[2]*x*x))",x1_fit_range,x2_fit_range);
f1->SetParameter(0,N_start);
f1->SetParameter(1,Lambda_start);
f1->SetParameter(2,R_inv_start);
hCF_Pi_Plus->Fit("f1");

//Write all:
TFile *f_out = new TFile(Output_File, "RECREATE");
hCF_Pi_Plus->Write();
hCF_Pi_Minus->Write();


//pictures:
TCanvas *c1 = new TCanvas("c1", "Canvas",1920,1080);
hCF_Pi_Plus->Draw();
c1->SaveAs(Output_Folder + "hCF_Pi_Plus.pdf");

hCF_Pi_Minus->Draw();
c1->SaveAs(Output_Folder + "hCF_Pi_Minus.pdf");

f_out->Close();
f->Close();
}