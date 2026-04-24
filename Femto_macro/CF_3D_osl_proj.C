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

const TString Input_File = "/home/kirill/root-on-vs-code/Femto_input/out_Au_Au200_23_04_FMRW.root";
const TString Output_File = "/home/kirill/root-on-vs-code/Femto_output/23_04_FMRW/out_Au_Au200_23_04_FMR_proseed_1D_3D_proj.root";
const TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_output/23_04_FMRW/";
const TString OSL_xyz[3] {"x","y","z"};
const TString OSL_names[3] {"out","side","long"};

void CF_3D_osl_proj()
{
TFile *f = TFile::Open(Input_File, "READ");
if (!f || f->IsZombie()) {std::cout << "Ohhhhh" << std::endl; exit(0); }
//osl 3D: 
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
            if (!h_Arr_pionters_3D[ihType][iCh][iCent][iKt]) {std::cout << "Ohhhh hist" << std::endl; exit(0); }
            h_Arr_3D[ihType][iCh][iCent][iKt] = (TH3F*)h_Arr_pionters_3D[ihType][iCh][iCent][iKt]->Clone(Form("CF_3D_non_norm_%i_%i_%i_%i",ihType,iCh,iCent,iKt));
        }
        }
    }
}

// //!!! WARNING NEXT STRINGSs are NECECERY:
// for (Int_t iCh = 0; iCh < N_Charge; iCh++)
// {
//     for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
//     {
//         for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
//         {
//             //3D:
//             h_Arr_3D[0][iCh][iCent][iKt]->Sumw2();
//             h_Arr_3D[1][iCh][iCent][iKt]->Sumw2();
//             h_Arr_3D[2][iCh][iCent][iKt]->Sumw2();
//         }
//     }
// }
// 1) Calculating Norm. coeff:
const Double_t q_i_min = 0.1;
const Double_t q_i_max = 0.2;

//A and B hists must be with same binning via all directions - out,side,long
Double_t q_i_min_max_bins[4][N_Charge][N_Bins_Centr][N_Bins_Kt][3] = {0.0};
// ----[q_min_max[0], q_min_max[1]]---0---[q_min_max[2], q_min_max[3]]----

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {

            q_i_min_max_bins[0][iCh][iCent][iKt][0] = h_Arr_3D[0][iCh][iCent][iKt]->GetXaxis()->FindBin(-q_i_max);
            q_i_min_max_bins[0][iCh][iCent][iKt][1] = h_Arr_3D[0][iCh][iCent][iKt]->GetYaxis()->FindBin(-q_i_max);
            q_i_min_max_bins[0][iCh][iCent][iKt][2] = h_Arr_3D[0][iCh][iCent][iKt]->GetZaxis()->FindBin(-q_i_max);


            q_i_min_max_bins[1][iCh][iCent][iKt][0] = h_Arr_3D[0][iCh][iCent][iKt]->GetXaxis()->FindBin(-q_i_min);
            q_i_min_max_bins[1][iCh][iCent][iKt][1] = h_Arr_3D[0][iCh][iCent][iKt]->GetYaxis()->FindBin(-q_i_min);
            q_i_min_max_bins[1][iCh][iCent][iKt][2] = h_Arr_3D[0][iCh][iCent][iKt]->GetZaxis()->FindBin(-q_i_min);

            q_i_min_max_bins[2][iCh][iCent][iKt][0] = h_Arr_3D[0][iCh][iCent][iKt]->GetXaxis()->FindBin(q_i_min);
            q_i_min_max_bins[2][iCh][iCent][iKt][1] = h_Arr_3D[0][iCh][iCent][iKt]->GetYaxis()->FindBin(q_i_min);
            q_i_min_max_bins[2][iCh][iCent][iKt][2] = h_Arr_3D[0][iCh][iCent][iKt]->GetZaxis()->FindBin(q_i_min);

            q_i_min_max_bins[3][iCh][iCent][iKt][0] = h_Arr_3D[0][iCh][iCent][iKt]->GetXaxis()->FindBin(q_i_max);
            q_i_min_max_bins[3][iCh][iCent][iKt][1] = h_Arr_3D[0][iCh][iCent][iKt]->GetYaxis()->FindBin(q_i_max);
            q_i_min_max_bins[3][iCh][iCent][iKt][2] = h_Arr_3D[0][iCh][iCent][iKt]->GetZaxis()->FindBin(q_i_max);

            // for(int i=0;i<4;i++)
            // {
            // std::cout<<q_i_min_max_bins[i][iCh][iCent][iKt][0]<<std::endl;
            // }
        }
    }
}

Double_t A_integrals_arr[N_Charge][N_Bins_Centr][N_Bins_Kt] = {0.0};
Double_t B_integrals_arr[N_Charge][N_Bins_Centr][N_Bins_Kt] = {0.0};
Double_t Norm_coeff_arr[N_Charge][N_Bins_Centr][N_Bins_Kt] = {0.0}; // Norm_coeff = B/A


for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {

            A_integrals_arr[iCh][iCent][iKt]= h_Arr_3D[0][iCh][iCent][iKt]->Integral(
                q_i_min_max_bins[0][iCh][iCent][iKt][0], q_i_min_max_bins[1][iCh][iCent][iKt][0],
                q_i_min_max_bins[0][iCh][iCent][iKt][1], q_i_min_max_bins[1][iCh][iCent][iKt][1],
                q_i_min_max_bins[0][iCh][iCent][iKt][2], q_i_min_max_bins[1][iCh][iCent][iKt][2]) 

                + h_Arr_3D[0][iCh][iCent][iKt]->Integral(
                q_i_min_max_bins[2][iCh][iCent][iKt][0], q_i_min_max_bins[3][iCh][iCent][iKt][0],
                q_i_min_max_bins[2][iCh][iCent][iKt][1], q_i_min_max_bins[3][iCh][iCent][iKt][1],
                q_i_min_max_bins[2][iCh][iCent][iKt][2], q_i_min_max_bins[3][iCh][iCent][iKt][2]
                );
            B_integrals_arr[iCh][iCent][iKt]= h_Arr_3D[1][iCh][iCent][iKt]->Integral(
                q_i_min_max_bins[0][iCh][iCent][iKt][0], q_i_min_max_bins[1][iCh][iCent][iKt][0],
                q_i_min_max_bins[0][iCh][iCent][iKt][1], q_i_min_max_bins[1][iCh][iCent][iKt][1],
                q_i_min_max_bins[0][iCh][iCent][iKt][2], q_i_min_max_bins[1][iCh][iCent][iKt][2]) 

                + h_Arr_3D[1][iCh][iCent][iKt]->Integral(
                q_i_min_max_bins[2][iCh][iCent][iKt][0], q_i_min_max_bins[3][iCh][iCent][iKt][0],
                q_i_min_max_bins[2][iCh][iCent][iKt][1], q_i_min_max_bins[3][iCh][iCent][iKt][1],
                q_i_min_max_bins[2][iCh][iCent][iKt][2], q_i_min_max_bins[3][iCh][iCent][iKt][2]
                );

            Norm_coeff_arr[iCh][iCent][iKt] = B_integrals_arr[iCh][iCent][iKt]/A_integrals_arr[iCh][iCent][iKt];
            std::cout<< A_integrals_arr[iCh][iCent][iKt]<< " "<< B_integrals_arr[iCh][iCent][iKt]<< " " <<iCent<<" "<< Norm_coeff_arr[iCh][iCent][iKt]<<std::endl;
        }
    }
}
//Scaling 3D A in CF:
for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            h_Arr_3D[0][iCh][iCent][iKt]->Scale(Norm_coeff_arr[iCh][iCent][iKt]);
        }
    }
}


// 2) Getting projections:
const Double_t range = 0.05;//Gev/c              
TH1F* h_Arr_3D_Projects_OSL[N_hist_types_3D][N_Charge][N_Bins_Centr][N_Bins_Kt][3];

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
                for (Int_t ihType=0; ihType < N_hist_types_3D; ihType++)
                {
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(-range,range);
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][0] = (TH1F*)h_Arr_3D[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[0]);
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][0]->SetName(Form("CF_3D_non_norm_%i_%i_%i_%i_%s",ihType,iCh,iCent,iKt,OSL_names[0].Data()));
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetYaxis()->SetRange();
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetZaxis()->SetRange();

                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(-range,range);
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][1] = (TH1F*)h_Arr_3D[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[1]);
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][1]->SetName(Form("CF_3D_non_norm_%i_%i_%i_%i_%s",ihType,iCh,iCent,iKt,OSL_names[1].Data()));
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetXaxis()->SetRange();
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetZaxis()->SetRange();

                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(-range,range);
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][2] = (TH1F*)h_Arr_3D[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[2]);
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][2]->SetName(Form("CF_3D_non_norm_%i_%i_%i_%i_%s",ihType,iCh,iCent,iKt,OSL_names[2].Data()));
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetXaxis()->SetRange();
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetYaxis()->SetRange();
                }
        }
    }
}
// //Deviding A/B:
for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
                    for(Int_t iOSL =0;iOSL<3;iOSL++)
                    {
                        // h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->Sumw2();
                        // h_Arr_3D_Projects_OSL[1][iCh][iCent][iKt][iOSL]->Sumw2();
                        h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->Divide(h_Arr_3D_Projects_OSL[1][iCh][iCent][iKt][iOSL]);
                    }
        }
    }
}


TFile *f_out = new TFile(Output_File, "RECREATE");
//Writing CF projections:
for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
                    for(Int_t iOSL =0;iOSL<3;iOSL++)
                    {
                        h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->Write();
                    }
        }
    }
}

TCanvas *c_3D= new TCanvas("c_3D_Plus", "Canvas",1920,1080);
c_3D->Divide(3);

for (Int_t iOSL = 0; iOSL < 3; iOSL++)
{
    
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->GetXaxis()->SetRangeUser(-0.1, 0.1);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->GetYaxis()->SetRangeUser(0.8, 1.8);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->SetTitle(Form("CF_3D_Projections_Non_norm Charge +, Centrality 0-5%%,  K_t [0.45,0.60] %s", OSL_names[iOSL].Data()));
    c_3D->cd(iOSL+1);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->Draw();
    
}
c_3D->SaveAs(Output_Folder + "c_3D_1.pdf");
for (Int_t iOSL = 0; iOSL < 3; iOSL++)
{
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->GetXaxis()->SetRangeUser(-0.4, 0.4);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->GetYaxis()->SetRangeUser(0.8, 1.8);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->SetTitle(Form("CF_3D_Projections_Non_norm Charge +, Centrality 0-5%%, K_t [0.45,0.60] %s", OSL_names[iOSL].Data()));
    c_3D->cd(iOSL+1);
    h_Arr_3D_Projects_OSL[0][0][N_Bins_Centr-1][3][iOSL]->Draw();
    
}
c_3D->SaveAs(Output_Folder + "c_3D_2.pdf");

c_3D->Write();
f_out->Close();
f->Close();
}