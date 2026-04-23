#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>
using namespace::std;

const Int_t N_hist_types_3D_integral = 2; // A or B
const Int_t N_Charge = 2; // 0 <-> Pi_Plus; 1 <-> Pi_Minus


const TString Input_File = "/home/kirill/root-on-vs-code/Femto_input/out_Au_Au200_23_04_FMRW.root";
const TString Output_File = "/home/kirill/root-on-vs-code/Femto_output/23_04_FMRW/out_Au_Au200_23_04_FMR_proseed_1D_3D.root";
const TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_output/23_04_FMRW/";
const TString OSL_xyz[3] {"x","y","z"};
const TString OSL_names[3] {"out","side","long"};

void CF_3D_osl_integral()
{
TFile *f = TFile::Open(Input_File, "READ");
if (!f || f->IsZombie()) {std::cout << "Ohhhhh" << std::endl; exit(0); }
//osl 3D: 
TString hist_3D_Name_AB[2] = {"hA_Pi_Plus_q_osl","hB_Pi_Plus_q_osl"};


TH3F* h_Arr_pionters_3D_integral[N_hist_types_3D_integral];
TH3F* h_Arr_3D_integral[N_hist_types_3D_integral];



for (Int_t ihType = 0; ihType < N_hist_types_3D_integral; ihType++)
{
    h_Arr_pionters_3D_integral[ihType] = (TH3F *)f->Get(hist_3D_Name_AB[ihType]);
    // 3D:
    if (!h_Arr_pionters_3D_integral[ihType])
    {
        std::cout << "Ohhhh hist" << std::endl;
        exit(0);
    }
    h_Arr_3D_integral[ihType] = (TH3F *)h_Arr_pionters_3D_integral[ihType]->Clone(hist_3D_Name_AB[ihType]);
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
const Double_t q_i_min = 0.05;
const Double_t q_i_max = 0.15;

//A and B hists must be with same binning via all directions - out,side,long
Double_t q_i_min_max_bins[4][N_Charge][3] = {0.0};
// ----[q_min_max[0], q_min_max[1]]---0---[q_min_max[2], q_min_max[3]]----

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{

            q_i_min_max_bins[0][iCh][0] = h_Arr_3D_integral[0]->GetXaxis()->FindBin(-q_i_max);
            q_i_min_max_bins[0][iCh][1] = h_Arr_3D_integral[0]->GetYaxis()->FindBin(-q_i_max);
            q_i_min_max_bins[0][iCh][2] = h_Arr_3D_integral[0]->GetZaxis()->FindBin(-q_i_max);


            q_i_min_max_bins[1][iCh][0] = h_Arr_3D_integral[0]->GetXaxis()->FindBin(-q_i_min);
            q_i_min_max_bins[1][iCh][1] = h_Arr_3D_integral[0]->GetYaxis()->FindBin(-q_i_min);
            q_i_min_max_bins[1][iCh][2] = h_Arr_3D_integral[0]->GetZaxis()->FindBin(-q_i_min);

            q_i_min_max_bins[2][iCh][0] = h_Arr_3D_integral[0]->GetXaxis()->FindBin(q_i_min);
            q_i_min_max_bins[2][iCh][1] = h_Arr_3D_integral[0]->GetYaxis()->FindBin(q_i_min);
            q_i_min_max_bins[2][iCh][2] = h_Arr_3D_integral[0]->GetZaxis()->FindBin(q_i_min);

            q_i_min_max_bins[3][iCh][0] = h_Arr_3D_integral[0]->GetXaxis()->FindBin(q_i_max);
            q_i_min_max_bins[3][iCh][1] = h_Arr_3D_integral[0]->GetYaxis()->FindBin(q_i_max);
            q_i_min_max_bins[3][iCh][2] = h_Arr_3D_integral[0]->GetZaxis()->FindBin(q_i_max);

            // for(int i=0;i<4;i++)
            // {
            // std::cout<<q_i_min_max_bins[i][0][0]<<std::endl;
            // }
}

Double_t A_integrals_arr[N_Charge] = {0.0};
Double_t B_integrals_arr[N_Charge] = {0.0};
Double_t Norm_coeff_arr[N_Charge] = {0.0}; // Norm_coeff = B/A


for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
            A_integrals_arr[iCh]= h_Arr_3D_integral[0]->Integral(
                q_i_min_max_bins[0][iCh][0], q_i_min_max_bins[1][iCh][0],
                q_i_min_max_bins[0][iCh][1], q_i_min_max_bins[1][iCh][1],
                q_i_min_max_bins[0][iCh][2], q_i_min_max_bins[1][iCh][2]) 

                + h_Arr_3D_integral[0]->Integral(
                q_i_min_max_bins[2][iCh][0], q_i_min_max_bins[3][iCh][0],
                q_i_min_max_bins[2][iCh][1], q_i_min_max_bins[3][iCh][1],
                q_i_min_max_bins[2][iCh][2], q_i_min_max_bins[3][iCh][2]
                );
            B_integrals_arr[iCh]= h_Arr_3D_integral[1]->Integral(
                q_i_min_max_bins[0][iCh][0], q_i_min_max_bins[1][iCh][0],
                q_i_min_max_bins[0][iCh][1], q_i_min_max_bins[1][iCh][1],
                q_i_min_max_bins[0][iCh][2], q_i_min_max_bins[1][iCh][2]) 

                + h_Arr_3D_integral[1]->Integral(
                q_i_min_max_bins[2][iCh][0], q_i_min_max_bins[3][iCh][0],
                q_i_min_max_bins[2][iCh][1], q_i_min_max_bins[3][iCh][1],
                q_i_min_max_bins[2][iCh][2], q_i_min_max_bins[3][iCh][2]
                );

            Norm_coeff_arr[iCh] = B_integrals_arr[iCh]/A_integrals_arr[iCh];
            std::cout<< A_integrals_arr[iCh]<< " "<< B_integrals_arr[iCh]<<" "<< Norm_coeff_arr[iCh]<<std::endl;
}
//Scaling 3D A in CF:
// for (Int_t iCh = 0; iCh < N_Charge; iCh++)
// {
            h_Arr_3D_integral[0]->Scale(Norm_coeff_arr[0]);
// }


// 2) Getting projections:
const Double_t range = 0.05;//Gev/c              
TH1F* h_Arr_3D_Projects_integral_OSL[N_hist_types_3D_integral][N_Charge][3];

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
                for (Int_t ihType=0; ihType < N_hist_types_3D_integral; ihType++)
                {
                    h_Arr_3D_integral[ihType]->GetYaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_integral[ihType]->GetZaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_Projects_integral_OSL[ihType][iCh][0] = (TH1F*)h_Arr_3D_integral[ihType]->Project3D(OSL_xyz[0]);
                    h_Arr_3D_Projects_integral_OSL[ihType][iCh][0]->SetName(hist_3D_Name_AB[ihType]+OSL_names[0]);
                    h_Arr_3D_integral[ihType]->GetYaxis()->SetRange();
                    h_Arr_3D_integral[ihType]->GetZaxis()->SetRange();

                    h_Arr_3D_integral[ihType]->GetXaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_integral[ihType]->GetZaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_Projects_integral_OSL[ihType][iCh][1] = (TH1F*)h_Arr_3D_integral[ihType]->Project3D(OSL_xyz[1]);
                    h_Arr_3D_Projects_integral_OSL[ihType][iCh][1]->SetName(hist_3D_Name_AB[ihType]+OSL_names[1]);
                    h_Arr_3D_integral[ihType]->GetXaxis()->SetRange();
                    h_Arr_3D_integral[ihType]->GetZaxis()->SetRange();

                    h_Arr_3D_integral[ihType]->GetXaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_integral[ihType]->GetYaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_Projects_integral_OSL[ihType][iCh][2] = (TH1F*)h_Arr_3D_integral[ihType]->Project3D(OSL_xyz[2]);
                    h_Arr_3D_Projects_integral_OSL[ihType][iCh][2]->SetName(hist_3D_Name_AB[ihType]+OSL_names[2]);
                    h_Arr_3D_integral[ihType]->GetXaxis()->SetRange();
                    h_Arr_3D_integral[ihType]->GetYaxis()->SetRange();
                }
}
// //Deviding A/B:
for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
                    for(Int_t iOSL =0;iOSL<3;iOSL++)
                    {
                        // h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->Sumw2();
                        // h_Arr_3D_Projects_OSL[1][iCh][iCent][iKt][iOSL]->Sumw2();
                        h_Arr_3D_Projects_integral_OSL[0][iCh][iOSL]->Divide(h_Arr_3D_Projects_integral_OSL[1][iCh][iOSL]);
                    }
}


TFile *f_out = new TFile(Output_File, "RECREATE");
TCanvas *c_3D= new TCanvas("c_3D_Plus", "Canvas",1920,1080);
c_3D->Divide(3);

for (Int_t iOSL = 0; iOSL < 3; iOSL++)
{
    
    h_Arr_3D_Projects_integral_OSL[0][0][iOSL]->GetXaxis()->SetRangeUser(-0.1, 0.1);
    //h_Arr_3D_Projects_integral_OSL[0][0][iOSL]->GetYaxis()->SetRangeUser(0.8, 1.8);
    h_Arr_3D_Projects_integral_OSL[0][0][iOSL]->SetTitle(Form("CF_3D_Projections Charge + %s", OSL_names[iOSL].Data()));
    c_3D->cd(iOSL+1);
    h_Arr_3D_Projects_integral_OSL[0][0][iOSL]->Draw();
    
}
c_3D->SaveAs(Output_Folder + "c_3D_integral_1.pdf");
for (Int_t iOSL = 0; iOSL < 3; iOSL++)
{
    h_Arr_3D_Projects_integral_OSL[0][0][iOSL]->GetXaxis()->SetRangeUser(-0.4, 0.4);
    //h_Arr_3D_Projects_integral_OSL[0][0][iOSL]->GetYaxis()->SetRangeUser(0.8, 1.8);
    h_Arr_3D_Projects_integral_OSL[0][0][iOSL]->SetTitle(Form("CF_3D_Projections Charge +, %s", OSL_names[iOSL].Data()));
    c_3D->cd(iOSL+1);
    h_Arr_3D_Projects_integral_OSL[0][0][iOSL]->Draw();
    
}
c_3D->SaveAs(Output_Folder + "c_3D_integral_2.pdf");
c_3D->SaveAs(Output_Folder + "c_3D_integral_2.png");

c_3D->Write();
f_out->Close();
f->Close();
}