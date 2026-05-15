//this macro creates .root file with projections of 3D_CF for each [charge][Centrality][K_t][projection]
//also this .root file will contain 3D_CF & B_weighted for each [charge][Centrality][K_t] 

//this macro Fits with &without Coulomb 3D_CF 
//this macro projects Fit surface on the OSL projections
//this macro does this for each [charge][centrality][k_t]
#include "TFile.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include <iostream>
using namespace::std;

const Int_t N_hist_types_3D = 3; // A or B or B_weighted
//const Int_t N_hist_types_1D = 2; //A or B
const Int_t N_Charge = 2; // 0 <-> Pi_Plus; 1 <-> Pi_Minus
const Int_t N_Bins_Kt = 4;
const Int_t N_Bins_Centr = 9; 

const Int_t N_of_Pars_for_fit = 5;

//const Double_t KtBins[N_Bins_Kt+1] = {0.15, 0.25, 0.35, 0.45, 0.60};

const TString Input_File = "/home/kirill/root-on-vs-code/Femto_input/out_Au_Au200_07_05_FMRW_article.root";
const TString Output_File = "/home/kirill/root-on-vs-code/Femto_output/15_05_FMRW_article/out_Au_Au200_15_05_FMR_pr_3D_proj.root";
const TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_output/15_05_FMRW_article/";
const TString Input_K_file = "/home/kirill/root-on-vs-code/Femto_macro/Kqinv_R_5fm.root";
const TString OSL_xyz[3] {"x","y","z"};
const TString OSL_names[3] {"out","side","long"};

//Info about K coulomb hist:
const Int_t N_bins_K_coul = 40;
const Double_t Bin_K_coul_size = 0.01; // GeV/c
std::vector<Double_t> K_values_arr(N_bins_K_coul, 0.);

TFile *f_K = TFile::Open(Input_K_file, "READ");
TH1D *hQCoul = (TH1D *)f_K->Get("hQCoul");
// Fit function:
// parameters: N,lambda; R_o,R_s,R_l; R_os,R_ol <-> 0,1; 2,3,4; 5,6
Double_t fitf_No_Coulomb(Double_t *x, Double_t *par)
{
    Double_t fitval = par[0]  * ((1.) + par[1] * 
    TMath::Exp(-(par[2] * x[0] * x[0] + par[3]*x[1]*x[1] + par[4] * x[2]*x[2]
    ) / (0.197327 * 0.197327)));

    return fitval;
}


//3D map of q_inv !!!fill it in code!!!
TH3F* h_Map_3D = nullptr;

Double_t fitf_Coulomb(Double_t *x, Double_t *par)
{
    Double_t q_inv_from_B_w = h_Map_3D->GetBinContent(h_Map_3D->FindBin(x[0],x[1],x[2]));
    Double_t K = hQCoul->GetBinContent(hQCoul->FindBin(q_inv_from_B_w));
    Double_t fitval = par[0]  * ((1. - par[1]) + par[1] * K *(1.+
    TMath::Exp(-(par[2] * x[0] * x[0] + par[3]*x[1]*x[1] + par[4] * x[2]*x[2]
    ) / (0.197327 * 0.197327))));

    return fitval;
}

void CF_3D_osl_proj()
{
TFile *f = TFile::Open(Input_File, "READ");
if (!f || f->IsZombie()) {std::cout << "Ohhhhh" << std::endl; exit(0); }
//osl 3D: 
TString hist_3D_Name = "h_3D";

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
            h_Arr_3D[ihType][iCh][iCent][iKt] = (TH3F* )f->Get(Form("%s_%i_%i_%i_%i",hist_3D_Name.Data(),ihType,iCh,iCent,iKt));
        }
        }
    }
}

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
            //std::cout<< A_integrals_arr[iCh][iCent][iKt]<< " "<< B_integrals_arr[iCh][iCent][iKt]<< " " <<iCent<<" "<< Norm_coeff_arr[iCh][iCent][iKt]<<std::endl;
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
                        h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->Divide(h_Arr_3D_Projects_OSL[1][iCh][iCent][iKt][iOSL]);
                    }
                    h_Arr_3D[0][iCh][iCent][iKt]->Divide(h_Arr_3D[1][iCh][iCent][iKt]);
        }
    }
}

//Now let's fit:

//Fill array of K(q_inv):
for (Int_t i = 1; i <= N_bins_K_coul; i++)
{
    K_values_arr[i - 1] = hQCoul->GetBinContent(i);
}


//Fit params:
Double_t par[N_Charge][N_Bins_Centr][N_Bins_Kt][N_of_Pars_for_fit];
TF3* f_No_Coul_Arr[N_Charge][N_Bins_Centr][N_Bins_Kt];

const Double_t x1_fit_range = -0.20;
const Double_t x2_fit_range = 0.20;
const Double_t N_start = 0.997251;
const Double_t Lambda_start = 0.31636;
const Double_t R_o_2_start = 43.7445;
const Double_t R_s_2_start = 19.1359;
const Double_t R_l_2_start = 25.6675;
// const Double_t R_os_2_start = 4.0;
// const Double_t R_ol_2_start = 4.0;

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {

            h_Arr_3D[0][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(x1_fit_range, x2_fit_range);
            h_Arr_3D[0][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(x1_fit_range, x2_fit_range);
            h_Arr_3D[0][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(x1_fit_range, x2_fit_range);

            TString fitf_name = TString::Format("f_No_Coul_%d_%d_%d", iCh, iCent, iKt);
            f_No_Coul_Arr[iCh][iCent][iKt]= new TF3(fitf_name, fitf_No_Coulomb,  x1_fit_range, x2_fit_range, 
                                                                            x1_fit_range, x2_fit_range, 
                                                                            x1_fit_range, x2_fit_range, 5);

            // f_No_Coul->FixParameter(0,N_start);
            f_No_Coul_Arr[iCh][iCent][iKt]->SetParameter(0, N_start);
            f_No_Coul_Arr[iCh][iCent][iKt]->SetParameter(1, Lambda_start);

            f_No_Coul_Arr[iCh][iCent][iKt]->SetParameter(2, R_o_2_start);
            f_No_Coul_Arr[iCh][iCent][iKt]->SetParameter(3, R_s_2_start);
            f_No_Coul_Arr[iCh][iCent][iKt]->SetParameter(4, R_l_2_start);

            // f_No_Coul->SetParameter(5, R_os_2_start);
            // f_No_Coul->SetParameter(6, R_ol_2_start);
            f_No_Coul_Arr[iCh][iCent][iKt]->SetParLimits(0, 0., 2.);

            f_No_Coul_Arr[iCh][iCent][iKt]->SetParLimits(1, 0., 1.);

            f_No_Coul_Arr[iCh][iCent][iKt]->SetParLimits(2, 0.01, 200.);
            f_No_Coul_Arr[iCh][iCent][iKt]->SetParLimits(3, 0.01, 200.);
            f_No_Coul_Arr[iCh][iCent][iKt]->SetParLimits(4, 0.01, 200.);

            f_No_Coul_Arr[iCh][iCent][iKt]->SetNpx(20);
            h_Arr_3D[0][iCh][iCent][iKt]->Fit(fitf_name);

            f_No_Coul_Arr[iCh][iCent][iKt]->GetParameters(par[iCh][iCent][iKt]);
}}}


TF3 *f_Coul = new TF3("f_Coul", fitf_Coulomb, x1_fit_range, x2_fit_range, x1_fit_range, x2_fit_range, x1_fit_range, x2_fit_range, 5);
f_Coul->SetParameter(0, N_start);
f_Coul->SetParameter(1, Lambda_start);

f_Coul->SetParameter(2, R_o_2_start);
f_Coul->SetParameter(3, R_s_2_start);
f_Coul->SetParameter(4, R_l_2_start);

// f_Coul->SetParameter(5, R_os_2_start);
// f_Coul->SetParameter(6, R_ol_2_start);
f_Coul->SetParLimits(0,0.,2.0);

f_Coul->SetParLimits(1,0.,1.);

f_Coul->SetParLimits(2,0.01,200.);
f_Coul->SetParLimits(3,0.01,200.);
f_Coul->SetParLimits(4,0.01,200.);

//!!! fill Map hist !!!
h_Map_3D = (TH3F*)h_Arr_3D[2][0][8][1]->Clone("h_Map_3D");
h_Map_3D->Divide(h_Arr_3D[1][0][8][1]);

f_Coul->SetNpx(20);
h_Arr_3D[0][0][8][1]->Fit("f_Coul");

// TFitResultPtr fitres = h_Arr_3D[0][0][8][1]->Fit("f_Coul");
// fitres->PrintCovMatrix(std::cout);




TFile *f_out = new TFile(Output_File, "RECREATE");
//Writing CF projections & CF_3D B_weighted:
for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
                    for(Int_t iOSL =0;iOSL<3;iOSL++)
                    {
                        //CF_projections:
                        h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->Write();
                    }
                    //CF:
                    h_Arr_3D[0][iCh][iCent][iKt]->Write();
                    //B_weighted:
                    h_Arr_3D[2][iCh][iCent][iKt]->Write();
        }
    }
}

//Hists for drawing projections of fit:
TH3F* h_Fit_No_Coul[2][N_Charge][N_Bins_Centr][N_Bins_Kt];
TH1F* h_Fit_No_Coul_Proj[2][N_Charge][N_Bins_Centr][N_Bins_Kt][3];

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            h_Fit_No_Coul[0][iCh][iCent][iKt] = (TH3F*)h_Arr_3D[0][iCh][iCent][iKt]->Clone("hA_Fit_No_Coul");
            h_Fit_No_Coul[1][iCh][iCent][iKt] = (TH3F*)h_Arr_3D[0][iCh][iCent][iKt]->Clone("hB_Fit_No_Coul");

                    
//Fill A hist:
double x_arr[3];
for(Int_t i = 1;i<=h_Fit_No_Coul[0][iCh][iCent][iKt]->GetNbinsX();i++)
{
    x_arr[0]=h_Fit_No_Coul[0][iCh][iCent][iKt]->GetXaxis()->GetBinCenter(i);
    for(Int_t j = 1;j<=h_Fit_No_Coul[0][iCh][iCent][iKt]->GetNbinsY();j++)
    {
        x_arr[1]=h_Fit_No_Coul[0][iCh][iCent][iKt]->GetYaxis()->GetBinCenter(j);
        for(Int_t k = 1;k<=h_Fit_No_Coul[0][iCh][iCent][iKt]->GetNbinsZ();k++)
        {
            x_arr[2]=h_Fit_No_Coul[0][iCh][iCent][iKt]->GetZaxis()->GetBinCenter(k);

            Int_t bin = h_Fit_No_Coul[0][iCh][iCent][iKt]->GetBin(i,j,k);
            h_Fit_No_Coul[0][iCh][iCent][iKt]->SetBinContent(bin, f_No_Coul_Arr[iCh][iCent][iKt]->Eval(x_arr[0],x_arr[1],x_arr[2]));
        }
    }
}

//Fill B hist:
for(Int_t i = 1;i<=h_Fit_No_Coul[1][iCh][iCent][iKt]->GetNbinsX();i++)
{
    for(Int_t j = 1;j<=h_Fit_No_Coul[1][iCh][iCent][iKt]->GetNbinsY();j++)
    {
        for(Int_t k = 1;k<=h_Fit_No_Coul[1][iCh][iCent][iKt]->GetNbinsZ();k++)
        {
            h_Fit_No_Coul[1][iCh][iCent][iKt]->SetBinContent(i,j,k,1.0);
        }
    }
}

//For Divide fit function A/B:
for (Int_t ihType = 0; ihType < 2; ihType++)
{
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][0] = (TH1F *)h_Fit_No_Coul[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[0]);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][0]->SetName(Form("CF_3D_non_norm_%i_%i_%i_%i_%s", ihType, iCh, iCent, iKt, OSL_names[0].Data()));
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRange();
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRange();

    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][1] = (TH1F *)h_Fit_No_Coul[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[1]);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][1]->SetName(Form("CF_3D_non_norm_%i_%i_%i_%i_%s", ihType, iCh, iCent, iKt, OSL_names[1].Data()));
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRange();
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRange();

    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][2] = (TH1F *)h_Fit_No_Coul[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[2]);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][2]->SetName(Form("CF_3D_non_norm_%i_%i_%i_%i_%s", ihType, iCh, iCent, iKt, OSL_names[2].Data()));
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRange();
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRange();
}
// //Deviding A/B:

for (Int_t iOSL = 0; iOSL < 3; iOSL++)
{
    h_Fit_No_Coul_Proj[0][iCh][iCent][iKt][iOSL]->Divide(h_Fit_No_Coul_Proj[1][iCh][iCent][iKt][iOSL]);
    h_Fit_No_Coul_Proj[0][iCh][iCent][iKt][iOSL]->GetXaxis()->SetRangeUser(x1_fit_range,x2_fit_range);  
}
}}}
//----------------------------------------------------------------------------------------------------//

//Now Same for Coulomb:

//Hists for drawing projections of fit:
TH3F* h_Fit_Coul[2];
TH1F* h_Fit_Coul_Proj[2][3];
h_Fit_Coul[0] = (TH3F*)h_Arr_3D[0][0][8][1]->Clone("hA_Fit_Coul");
h_Fit_Coul[1] = (TH3F*)h_Arr_3D[0][0][8][1]->Clone("hB_Fit_Coul");

//Fill A hist:
double x_arr[3];
for(Int_t i = 1;i<=h_Fit_Coul[0]->GetNbinsX();i++)
{
    x_arr[0]=h_Fit_Coul[0]->GetXaxis()->GetBinCenter(i);
    for(Int_t j = 1;j<=h_Fit_Coul[0]->GetNbinsY();j++)
    {
        x_arr[1]=h_Fit_Coul[0]->GetYaxis()->GetBinCenter(j);
        for(Int_t k = 1;k<=h_Fit_Coul[0]->GetNbinsZ();k++)
        {
            x_arr[2]=h_Fit_Coul[0]->GetZaxis()->GetBinCenter(k);

            Int_t bin = h_Fit_Coul[0]->GetBin(i,j,k);
            h_Fit_Coul[0]->SetBinContent(bin, f_Coul->Eval(x_arr[0],x_arr[1],x_arr[2]));
        }
    }
}

//Fill B hist:
for(Int_t i = 1;i<=h_Fit_Coul[1]->GetNbinsX();i++)
{
    for(Int_t j = 1;j<=h_Fit_Coul[1]->GetNbinsY();j++)
    {
        for(Int_t k = 1;k<=h_Fit_Coul[1]->GetNbinsZ();k++)
        {
            h_Fit_Coul[1]->SetBinContent(i,j,k,1.0);
        }
    }
}

//For Divide fit function A/B:
for (Int_t ihType = 0; ihType < 2; ihType++)
{
    h_Fit_Coul[ihType]->GetYaxis()->SetRangeUser(-range, range);
    h_Fit_Coul[ihType]->GetZaxis()->SetRangeUser(-range, range);
    h_Fit_Coul_Proj[ihType][0] = (TH1F *)h_Fit_Coul[ihType]->Project3D(OSL_xyz[0]);
    h_Fit_Coul_Proj[ihType][0]->SetName(Form("CF_3D_non_norm_%i_%i_%i_%i_%s", ihType, 0, 8, 1, OSL_names[0].Data()));
    h_Fit_Coul[ihType]->GetYaxis()->SetRange();
    h_Fit_Coul[ihType]->GetZaxis()->SetRange();

    h_Fit_Coul[ihType]->GetXaxis()->SetRangeUser(-range, range);
    h_Fit_Coul[ihType]->GetZaxis()->SetRangeUser(-range, range);
    h_Fit_Coul_Proj[ihType][1] = (TH1F *)h_Fit_Coul[ihType]->Project3D(OSL_xyz[1]);
    h_Fit_Coul_Proj[ihType][1]->SetName(Form("CF_3D_non_norm_%i_%i_%i_%i_%s", ihType, 0, 8, 1, OSL_names[1].Data()));
    h_Fit_Coul[ihType]->GetXaxis()->SetRange();
    h_Fit_Coul[ihType]->GetZaxis()->SetRange();

    h_Fit_Coul[ihType]->GetXaxis()->SetRangeUser(-range, range);
    h_Fit_Coul[ihType]->GetYaxis()->SetRangeUser(-range, range);
    h_Fit_Coul_Proj[ihType][2] = (TH1F *)h_Fit_Coul[ihType]->Project3D(OSL_xyz[2]);
    h_Fit_Coul_Proj[ihType][2]->SetName(Form("CF_3D_non_norm_%i_%i_%i_%i_%s", ihType, 0, 8, 1, OSL_names[2].Data()));
    h_Fit_Coul[ihType]->GetXaxis()->SetRange();
    h_Fit_Coul[ihType]->GetYaxis()->SetRange();
}
// //Deviding A/B:

for (Int_t iOSL = 0; iOSL < 3; iOSL++)
{
    h_Fit_Coul_Proj[0][iOSL]->Divide(h_Fit_Coul_Proj[1][iOSL]);
    h_Fit_Coul_Proj[0][iOSL]->GetXaxis()->SetRangeUser(x1_fit_range,x2_fit_range);  
}

//Drawing projections:
TCanvas *c_3D_8_3_X= new TCanvas("c_3D_0_0_8_3_X", "Canvas",1920,1080);
TCanvas *c_3D_8_3_Y= new TCanvas("c_3D_0_0_8_3_Y", "Canvas",1920,1080);
TCanvas *c_3D_8_3_Z= new TCanvas("c_3D_0_0_8_3_Z", "Canvas",1920,1080);

c_3D_8_3_X->cd();
h_Arr_3D_Projects_OSL[0][0][8][1][0]->Draw();

h_Fit_Coul_Proj[0][0]->SetLineColor(kRed);
h_Fit_Coul_Proj[0][0]->Draw("L SAME");

h_Fit_No_Coul_Proj[0][0][8][1][0]->SetLineColor(kGreen);
h_Fit_No_Coul_Proj[0][0][8][1][0]->Draw("L SAME");

c_3D_8_3_X->SaveAs(Output_Folder + "c_3D_8_3_X.pdf");
c_3D_8_3_X->Write();

c_3D_8_3_Y->cd();
h_Arr_3D_Projects_OSL[0][0][8][1][1]->Draw();

h_Fit_Coul_Proj[0][1]->SetLineColor(kRed);
h_Fit_Coul_Proj[0][1]->Draw("L SAME");

h_Fit_No_Coul_Proj[0][0][8][1][1]->SetLineColor(kGreen);
h_Fit_No_Coul_Proj[0][0][8][1][1]->Draw("L SAME");

c_3D_8_3_Y->SaveAs(Output_Folder + "c_3D_8_3_Y.pdf");
c_3D_8_3_Y->Write();

c_3D_8_3_Z->cd();
h_Arr_3D_Projects_OSL[0][0][8][1][2]->Draw();

h_Fit_Coul_Proj[0][2]->SetLineColor(kRed);
h_Fit_Coul_Proj[0][2]->Draw("SAME L");

h_Fit_No_Coul_Proj[0][0][8][1][2]->SetLineColor(kGreen);
h_Fit_No_Coul_Proj[0][0][8][1][2]->Draw("L SAME");

c_3D_8_3_Z->SaveAs(Output_Folder + "c_3D_8_3_Z.pdf");
c_3D_8_3_Z->Write();



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


delete f_Coul;
f_Coul = nullptr;
for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            delete f_No_Coul_Arr[iCh][iCent][iKt];
            f_No_Coul_Arr[iCh][iCent][iKt] = nullptr;
}}}


f->Close();
}