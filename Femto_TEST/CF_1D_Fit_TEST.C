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

//here output file from CF_1D
const TString Input_File = "/home/kirill/root-on-vs-code/Femto_TEST/Output_TEST/Output_file_set_1_pr_1D_integral_and_cuts_TEST.root";
const TString Output_File = "/home/kirill/root-on-vs-code/Femto_TEST/Output_TEST/Output_file_set_1_pr_1D_integral_and_cuts_Fit_TEST.root";
const TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_TEST/Output_TEST/";

//Fit function:
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = par[0]*(1+par[1]*TMath::Exp(-(par[2]*par[2]*x[0]*x[0])/(0.197327*0.197327)));
    return fitval;
}

void CF_1D_Fit_TEST()
{
TFile *f = TFile::Open(Input_File, "READ");
//1D total, no cuts:

TH1D *hCF_Pi_Plus = (TH1D*)f->Get("CF_Pi_Plus");
TH1D *hCF_Pi_Minus = (TH1D*)f->Get("CF_Pi_Minus");
//
hCF_Pi_Plus->Sumw2();
hCF_Pi_Minus->Sumw2();

//Let's fit CF:
const Double_t x1_fit_range = 0.00;
const Double_t x2_fit_range = 0.20;
const Double_t N_start = 5.0;
const Double_t Lambda_start = 0.6;
const Double_t R_inv_start = 5.0;
const Double_t hbar_c = 197.327/1000.; //GeV*Fm


hCF_Pi_Plus->GetXaxis()->SetRangeUser(x1_fit_range,x2_fit_range);

TF1 *f1 = new TF1("f1","[0]*(1+[1]*exp(-([2]*[2]*x*x)/(0.197327*0.197327)))",x1_fit_range,x2_fit_range);
TF1 *fMy = new TF1("f2",fitf,x1_fit_range,x2_fit_range,3);

fMy->SetParameter(1,Lambda_start);
fMy->SetParameter(0,N_start);
fMy->SetParameter(2,R_inv_start);
hCF_Pi_Plus->Fit("f2");
hCF_Pi_Minus->Fit("f2");
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