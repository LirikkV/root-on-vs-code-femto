#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>
using namespace::std;

const Int_t N_hist_types_3D = 3; // A or B or B_weighted
const Int_t N_hist_types_1D = 2; //A or B
const Int_t N_Charge = 2; // 0 <-> Pi_Plus; 1 <-> Pi_Minus
const Int_t N_Bins_Kt = 4;
const Int_t N_Bins_Centr = 9; 

const Double_t KtBins[N_Bins_Kt+1] = {0.15, 0.25, 0.35, 0.45, 0.60};

TString Input_File = "/home/kirill/root-on-vs-code/Femto_output/out_Au_Au200_15_04_FMR.root";
TString Output_File = "/home/kirill/root-on-vs-code/Femto_output/out_Au_Au200_15_04_FMR_proseed_1D_3D.root";
TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_output/";
TString OSL_xyz[3] {"x","y","z"};
TString OSL_names[3] {"out","side","long"};

void CF_3D_osl_proj()
{
TFile *f = TFile::Open(Input_File, "READ");
if (!f || f->IsZombie()) {std::cout << "Ohhhhh" << std::endl; exit(0); }
//osl 3D:
// std::array<std::array<std::array<std::array<TH3F*,N_Bins_Kt>,N_Bins_Centr>,N_Charge>, N_hist_types_3D> h_Arr_pionters_3D = {}; 
// std::array<std::array<std::array<std::array<TH3F,N_Bins_Kt>,N_Bins_Centr>,N_Charge>, N_hist_types_3D> h_Arr_3D = {}; 
TString hist_3D_Name = "h_3D";

TH3F* h_Arr_pionters_3D[N_hist_types_3D][N_Charge][N_Bins_Centr][N_Bins_Kt];
TH3F* h_Arr_3D[N_hist_types_3D][N_Charge][N_Bins_Centr][N_Bins_Kt];


for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            for (Int_t ihType=0; ihType < N_hist_types_3D; ihType++)
                {
            //3D:
            h_Arr_pionters_3D[ihType][iCh][iCent][iKt] = (TH3F* )f->Get(Form("%s_%i_%i_%i_%i",hist_3D_Name.Data(),ihType,iCh,iCent,iKt));
            // cout<<h_Arr_pionters_3D[ihType][iCh][iCent][iKt]->GetEntries()<<endl;
            if (!h_Arr_pionters_3D[ihType][iCh][iCent][iKt]) {std::cout << "Ohhhhhist" << std::endl; exit(0); }
            h_Arr_3D[ihType][iCh][iCent][iKt] = (TH3F*)h_Arr_pionters_3D[ihType][iCh][iCent][iKt]->Clone(Form("CF_3D_non_norm_%i_%i_%i_%i",ihType,iCh,iCent,iKt));
        }
        }
    }
}

//!!! WARNING NEXT STRINGSs are NECECERY:
for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            //3D:
            h_Arr_3D[0][iCh][iCent][iKt]->Sumw2();
            h_Arr_3D[1][iCh][iCent][iKt]->Sumw2();
            h_Arr_3D[2][iCh][iCent][iKt]->Sumw2();
        }
    }
}
/*
Method: 
1)Get relation of A/B in range of q_inv where no correlations 
2)Scale A by relation from 1)
3)CF=A_normilized/B
*/

//Method 2 for out side long: first project A and B then divide:
//out

//3D with bins:                         
TH1F* h_Arr_3D_Projects_OSL[N_hist_types_3D][N_Charge][N_Bins_Centr][N_Bins_Kt][3];
//Getting projections:
for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            for (Int_t iOSL = 0; iOSL < 3; iOSL++)
            {
                for (Int_t ihType=0; ihType < N_hist_types_3D; ihType++)
                {
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][iOSL] = (TH1F*)h_Arr_3D[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[iOSL]);
                }
                h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->Divide(h_Arr_3D_Projects_OSL[1][iCh][iCent][iKt][iOSL]);
            }
        }
    }
}
TFile *f_out = new TFile(Output_File, "RECREATE");
TCanvas *c_3D= new TCanvas("c_3D_Plus", "Canvas",1920,1080);
c_3D->Divide(3);

for (Int_t iOSL = 0; iOSL < 3; iOSL++)
{
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->GetXaxis()->SetRangeUser(-0.1, 0.1);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->GetYaxis()->SetRangeUser(0.096, 0.11);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->SetTitle(Form("CF_3D_Projections_Non_norm Charge +, Centrality 0-5%, K_t [0.45,0.60]" + OSL_names[iOSL] ));
    c_3D->cd(iOSL+1);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->Draw();
    
}
c_3D->SaveAs(Output_Folder + "c_3D_1.pdf");
for (Int_t iOSL = 0; iOSL < 3; iOSL++)
{
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->GetXaxis()->SetRangeUser(-0.4, 0.4);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->GetYaxis()->SetRangeUser(0.096, 0.11);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->SetTitle(Form("CF_3D_Projections_Non_norm Charge +, Centrality 0-5%, K_t [0.45,0.60]" + OSL_names[iOSL] ));
    c_3D->cd(iOSL+1);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->Draw();
    
}
c_3D->SaveAs(Output_Folder + "c_3D_2.pdf");

c_3D->Write();
f_out->Close();
f->Close();
}