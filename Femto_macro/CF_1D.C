#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>
using namespace::std;

const Int_t N_hist_types_1D = 2; //A or B
const Int_t N_Charge = 2; // 0 <-> Pi_Plus; 1 <-> Pi_Minus
const Int_t N_Bins_Kt = 4;
const Int_t N_Bins_Centr = 9; 

const Double_t KtBins[N_Bins_Kt+1] = {0.15, 0.25, 0.35, 0.45, 0.60};

const TString Input_File = "/home/kirill/root-on-vs-code/Femto_input/out_Au_Au200_07_05_FMRW_article.root";
const TString Output_File = "/home/kirill/root-on-vs-code/Femto_output/07_05_FMRW_article/out_Au_Au200_07_05_FMRW_article_pr_1D_integral_and_cuts.root";
const TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_output/07_05_FMRW_article/";

void CF_1D()
{
TFile *f = TFile::Open(Input_File, "READ");
//1D total, no cuts:

//for Pi+Pi+:
TH1D *hA_Plus = (TH1D*)f->Get("hA_Pi_Plus_q_inv_ALL");
TH1D *hB_Plus = (TH1D*)f->Get("hB_Pi_Plus_q_inv_ALL");

//for Pi-Pi-:
TH1D *hA_Minus = (TH1D*)f->Get("hA_Pi_Minus_q_inv_ALL");
TH1D *hB_Minus = (TH1D*)f->Get("hB_Pi_Minus_q_inv_ALL");


//1D cuts hists:
std::array<std::array<std::array<std::array<TH1D*,N_Bins_Kt>,N_Bins_Centr>,N_Charge>, N_hist_types_1D> h_Arr_1D = {};

TString hist_1D_Name = "h_1D_";

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            //1D:
            h_Arr_1D[0][iCh][iCent][iKt] = (TH1D *)f->Get(hist_1D_Name + Form("%i_%i_%i_%i",0,iCh,iCent,iKt));
            h_Arr_1D[1][iCh][iCent][iKt] = (TH1D *)f->Get(hist_1D_Name + Form("%i_%i_%i_%i",1,iCh,iCent,iKt));
        }
    }
}

//!!! WARNING NEXT STRINGSs are NECECERY:
hA_Plus->Sumw2();
hB_Plus->Sumw2();
hA_Minus->Sumw2();
hB_Minus->Sumw2();

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            //1D:
            h_Arr_1D[0][iCh][iCent][iKt]->Sumw2();
            h_Arr_1D[1][iCh][iCent][iKt]->Sumw2();
        }
    }
}


//let's normalize CF:
Double_t q_min = 0.2;
Double_t q_max = 0.9;

//1D no cuts:
Double_t A_Pi_Plus_integral = hA_Plus->Integral(hA_Plus->FindBin(q_min), hA_Plus->FindBin(q_max));
Double_t B_Pi_Plus_integral = hB_Plus->Integral(hB_Plus->FindBin(q_min),hB_Plus->FindBin(q_max));
Double_t Scale_factor_Pi_Plus = B_Pi_Plus_integral/A_Pi_Plus_integral;

Double_t A_Pi_Minus_integral = hA_Minus->Integral(hA_Minus->FindBin(q_min), hA_Minus->FindBin(q_max));
Double_t B_Pi_Minus_integral = hB_Minus->Integral(hB_Minus->FindBin(q_min),hB_Minus->FindBin(q_max));
Double_t Scale_factor_Pi_Minus = B_Pi_Minus_integral/A_Pi_Minus_integral;

//1D with cuts:
Double_t A_integrals_arr[N_Charge][N_Bins_Centr][N_Bins_Kt] = {0.0};
Double_t B_integrals_arr[N_Charge][N_Bins_Centr][N_Bins_Kt] = {0.0};
Double_t Norm_coeff_arr[N_Charge][N_Bins_Centr][N_Bins_Kt] = {0.0}; // Norm_coeff = B/A

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            A_integrals_arr[iCh][iCent][iKt]=h_Arr_1D[0][iCh][iCent][iKt]->Integral(h_Arr_1D[0][iCh][iCent][iKt]->FindBin(q_min),
                                                                                    h_Arr_1D[0][iCh][iCent][iKt]->FindBin(q_max));
            B_integrals_arr[iCh][iCent][iKt]=h_Arr_1D[1][iCh][iCent][iKt]->Integral(h_Arr_1D[1][iCh][iCent][iKt]->FindBin(q_min),
                                                                                    h_Arr_1D[1][iCh][iCent][iKt]->FindBin(q_max));

            Norm_coeff_arr[iCh][iCent][iKt] = B_integrals_arr[iCh][iCent][iKt]/A_integrals_arr[iCh][iCent][iKt];
            std::cout<< A_integrals_arr[iCh][iCent][iKt]<< " "<< B_integrals_arr[iCh][iCent][iKt]<< " " <<iCent<<" "<< Norm_coeff_arr[iCh][iCent][iKt]<<std::endl;
        }
    }
}
//For hists scales:
Double_t x1_range = 0.0;
Double_t x2_range = 0.2;
Double_t y1_range = 0.4;
Double_t y2_range = 1.6;
/*
Method: 
1)Get relation of A/B in range of q_inv where no correlations 
2)Scale A by relation from 1)
3)CF=A_normilized/B
*/

//CF 1D no cuts:
TH1D *hA_Plus_normalized = (TH1D*)hA_Plus->Clone("CF_Pi_Plus");
hA_Plus_normalized->SetTitle("Numerator of Corr.Funct Pi+ Pi+ normalized");
hA_Plus_normalized->Scale(Scale_factor_Pi_Plus);
TH1D *CF_Pi_Plus = (TH1D*)hA_Plus_normalized->Clone("CF_Pi_Plus");
CF_Pi_Plus->Divide(hB_Plus);
CF_Pi_Plus->SetTitle("Corr.Funct Pi+ Pi+");
CF_Pi_Plus->GetXaxis()->SetTitle("q_inv");
CF_Pi_Plus->GetYaxis()->SetTitle("CF");
CF_Pi_Plus->GetYaxis()->SetRangeUser(y1_range,y2_range);
CF_Pi_Plus->GetXaxis()->SetRangeUser(x1_range,x2_range);

TH1D *hA_Minus_normalized = (TH1D*)hA_Minus->Clone("CF_Pi_Minus");
hA_Minus_normalized->SetTitle("Numerator of Corr.Funct Pi- Pi- normalized");
hA_Minus_normalized->Scale(Scale_factor_Pi_Minus);
TH1D *CF_Pi_Minus = (TH1D*)hA_Minus_normalized->Clone("CF_Pi_Minus");
CF_Pi_Minus->Divide(hB_Minus);
CF_Pi_Minus->SetTitle("Corr.Funct Pi- Pi-");
CF_Pi_Minus->GetXaxis()->SetTitle("q_inv");
CF_Pi_Minus->GetYaxis()->SetTitle("CF");
CF_Pi_Minus->GetYaxis()->SetRangeUser(y1_range,y2_range);
CF_Pi_Minus->GetXaxis()->SetRangeUser(x1_range,x2_range);

//CF 1D with cuts:
std::array<std::array<std::array<TH1D*,N_Bins_Kt>,N_Bins_Centr>,N_Charge> h_CF_1D = {};
for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            h_Arr_1D[0][iCh][iCent][iKt]->Scale(Norm_coeff_arr[iCh][iCent][iKt]);
            h_CF_1D[iCh][iCent][iKt] = (TH1D*)h_Arr_1D[0][iCh][iCent][iKt]->Clone(Form("CF_1D_%i_%i_%i",iCh,iCent,iKt) );
            h_CF_1D[iCh][iCent][iKt]->Divide(h_Arr_1D[1][iCh][iCent][iKt]);
            h_CF_1D[iCh][iCent][iKt]->SetTitle(Form("CF_1D Charge: %i Centrality: %i K_t: %i",iCh,iCent,iKt));
            h_CF_1D[iCh][iCent][iKt]->GetXaxis()->SetTitle("q_inv");
            h_CF_1D[iCh][iCent][iKt]->GetYaxis()->SetTitle("CF");
        }
    }
}

//Write all:
TFile *f_out = new TFile(Output_File, "RECREATE");
CF_Pi_Plus->Write();
CF_Pi_Minus->Write();

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            h_CF_1D[iCh][iCent][iKt]->Write();
        }
    }
}

//pictures:
TCanvas *c1 = new TCanvas("c1", "Canvas",1920,1080);
{
// hA_Plus->Draw();
// c1->SaveAs(Output_Folder + "hA_Plus.pdf");

// hA_Plus_normalized->Draw();
// c1->SaveAs(Output_Folder + "hA_Norm_Pi_Plus.pdf");

// hB_Plus->Draw();
// c1->SaveAs(Output_Folder + "hB_Plus.pdf");

// gStyle->SetOptStat(0);

// hB_Plus->SetLineColor(kRed);
// hA_Plus_normalized->SetLineColor(kGreen);

// hB_Plus->GetXaxis()->SetRangeUser(0.0,0.1);
// hA_Plus_normalized->GetXaxis()->SetRangeUser(0.0,0.1);

// hB_Plus->SetMarkerStyle(21);
// hA_Plus_normalized->SetMarkerStyle(20);

// hB_Plus->SetMarkerSize(0.5);
// hA_Plus_normalized->SetMarkerSize(0.5);

// hB_Plus->SetTitle("Normalized numerator(A) and denumerator (B)");
// hB_Plus->Draw();
// hA_Plus_normalized->Draw("SAME");

// TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
// leg->AddEntry(hB_Plus, "B", "lp");
// leg->AddEntry(hA_Plus_normalized, "A normalized" , "lp");
// leg->Draw();
// c1->SaveAs(Output_Folder + "hA_B_difference.pdf");
}
CF_Pi_Plus->Draw();
c1->SaveAs(Output_Folder + "hCF_Pi_Plus.pdf");

CF_Pi_Minus->Draw();
c1->SaveAs(Output_Folder + "hCF_Pi_Minus.pdf");

// TCanvas *c_1D_Plus = new TCanvas("c_1D_Plus", "Canvas",9000,4000);
// c_1D_Plus->Divide(N_Bins_Centr,N_Bins_Kt,0,0);

// for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
// {
//     for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
//     {
//         c_1D_Plus->cd((iCent + 1) * (iKt + 1));
//         h_CF_1D[0][iCent][iKt]->GetXaxis()->SetRangeUser(0.0, 0.3);
//         h_CF_1D[0][iCent][iKt]->SetTitle(Form("CF_1D Charge %i, Centrality %i, K_t [%.2f,%.2f]", 0, iCent,
//                                                        KtBins[iKt], KtBins[iKt + 1]));
//         h_CF_1D[0][iCent][iKt]->Draw();
//     }
// }

// c_1D_Plus->SaveAs(Output_Folder + "c_1D_Plus.pdf");

// c_1D_Plus->Write();
f_out->Close();
f->Close();
}