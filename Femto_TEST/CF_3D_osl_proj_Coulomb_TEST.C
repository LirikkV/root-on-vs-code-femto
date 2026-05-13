#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"

#include <iostream>
using namespace ::std;

const Int_t N_hist_types_1D = 2; // A or B
const Int_t N_Charge = 2;        // 0 <-> Pi_Plus; 1 <-> Pi_Minus
const Int_t N_Bins_Kt = 4;
const Int_t N_Bins_Centr = 9;

const Double_t KtBins[N_Bins_Kt + 1] = {0.15, 0.25, 0.35, 0.45, 0.60};

const TString Input_File =
    "/home/kirill/root-on-vs-code/Femto_input/out_Au_Au200_07_05_FMRW_article.root";
const TString Output_File =
    "/home/kirill/root-on-vs-code/Femto_output/07_05_FMRW_article/out_Au_Au200_07_05_FMR_3D_CF_Fit_CentrMost_KtLast.root";
const TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_output/07_05_FMRW_article/";

const TString Input_K_file = "/home/kirill/root-on-vs-code/Femto_macro/Kqinv_R_5fm.root";
const TString Input_CF_3D_file = 
"/home/kirill/root-on-vs-code/Femto_output/07_05_FMRW_article/out_Au_Au200_07_05_FMR_proseed_1D_3D_proj.root";

// Fit function:
const Int_t N_bins_K = 40;
const Double_t Bin_K_size = 0.01; // GeV/c
std::vector<Double_t> K_values_arr {0.52244,0.848399,0.926191,0.958098,0.973856,0.982901,0.988093,0.991415,0.993299,0.994838,0.995766,0.996424,0.997075,0.997387,0.997735,0.998068,0.998318,0.998518,0.998644,0.998784,0.998946,0.99896,0.99911,0.999183,0.999261,0.999322,0.999329,0.999418,0.999432,0.99948,0.999504,0.999489,0.999568,0.999598,0.999641,0.999632,0.999642,0.999662,0.999694,0.999708};


TFile *f_K = TFile::Open(Input_K_file, "READ");
TH1D *hQCoul = (TH1D *)f_K->Get("hQCoul");

Double_t fitf(Double_t *x, Double_t *par)
{
    Double_t K = hQCoul->GetBinContent(hQCoul->FindBin(x[0]));
    //std::cout<<K<<": "<<hQCoul->FindBin(x[0])<<std::endl;
    Double_t fitval = par[0] * ((1. - par[1]) + par[1] *K* TMath::Exp(-(par[2] * x[0] * x[0]) / (0.197327 * 0.197327)));
    return fitval;
    // Double_t fitval = par[0] * ((1.) + par[1]* TMath::Exp(-(par[2] * x[0] * x[0]) / (0.197327 * 0.197327)));
    // return fitval;
}

void CF_1D_Coulomb_Fit()
{
    

    for (Int_t i = 1; i <= N_bins_K; i++)
    {
        K_values_arr[i - 1] = hQCoul->GetBinContent(i);
    }

    TFile *f = TFile::Open(Input_File, "READ");
    // 1D total, no cuts:

    TH1D *hCF_Pi_Plus = (TH1D *)f->Get("CF_Pi_Plus");
    TH1D *hCF_Pi_Minus = (TH1D *)f->Get("CF_Pi_Minus");
    //
    hCF_Pi_Plus->Sumw2();
    hCF_Pi_Minus->Sumw2();

    // Let's fit CF:
    const Double_t x1_fit_range = 0.00;
    const Double_t x2_fit_range = 0.20;
    const Double_t N_start = 1.2;
    const Double_t Lambda_start = 0.2;
    const Double_t R_inv_2_start = 25.0;

    const Double_t hbar_c = 197.327 / 1000.; // GeV*Fm

    hCF_Pi_Plus->GetXaxis()->SetRangeUser(x1_fit_range, x2_fit_range);

    // TF1 *f1 = new TF1("f1","[0]*(1+[1]*exp(-([2]*[2]*x*x)/(0.197327*0.197327)))",x1_fit_range,x2_fit_range);
    TF1 *fMy = new TF1("f2", fitf, x1_fit_range, x2_fit_range, 3);

    fMy->SetParameter(0, N_start);
    fMy->SetParameter(1, Lambda_start);
    fMy->SetParameter(2, R_inv_2_start);
    fMy->SetParLimits(0, 0.0, 10.);
    fMy->SetParLimits(1, 0.0, 1.);
    fMy->SetParLimits(2, 0.0, 150.);
    fMy->SetNpx(20);
    hCF_Pi_Plus->Fit("f2");

    // fMy->SetParameter(0,N_start);
    // fMy->SetParameter(1,Lambda_start);
    // fMy->SetParameter(2,R_inv_2_start);
    // hCF_Pi_Minus->Fit("f2");

    // Write all:
    TFile *f_out = new TFile(Output_File, "RECREATE");
    hCF_Pi_Plus->Write();
    // hCF_Pi_Minus->Write();

    // pictures:
    TCanvas *c1 = new TCanvas("c1", "Canvas", 1920, 1080);
    hCF_Pi_Plus->Draw();
    c1->SaveAs(Output_Folder + "hCF_Pi_Plus.png");

    hCF_Pi_Minus->Draw();
    c1->SaveAs(Output_Folder + "hCF_Pi_Minus.png");

    f_out->Close();
    f->Close();
    f_K->Close();
}