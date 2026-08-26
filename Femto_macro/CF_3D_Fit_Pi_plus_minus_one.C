//This macro for plus & minus pions data in one CF
//this macro creates .root file with fit data
//also this .root file will contain 3D_CF & B_weighted for each [Centrality][K_t] 

//this macro Fits with &without Coulomb 3D_CF 
//this macro projects Fit surface on the OSL projections
//this macro does this for each [centrality][k_t]

//after this macro you can execute CF_3D_osl_proj_paint_No_Fit
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
#include "TPaveText.h"
#include "TTree.h"

#include <iostream>
using namespace std;

const Int_t N_hist_types_3D = 3; // A or B or B_weighted
//const Int_t N_hist_types_1D = 2; //A or B
const Int_t N_Charge = 2; // 0 <-> Pi_Plus; 1 <-> Pi_Minus
const Int_t N_Bins_Kt = 4;
const Int_t N_Bins_Centr = 9; 

const Int_t N_of_Pars_for_fit = 5;

const TString Charged_particles_pairs_titles[N_Charge] = {"Pi+Pi+","Pi-Pi-"};
const TString Centrality_region_titles[N_Bins_Centr] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
const TString Kt_bins_titles[N_Bins_Kt] = {"[0.15, 0.25]","[0.25, 0.35]","[0.35, 0.45]","[0.45, 0.60]"};

struct FitResult {
    Int_t    centrality;
    Int_t    Kt;
    Int_t    fitType;     // 0 = NoCoul, 1 = Coul
    Double_t par[N_of_Pars_for_fit];
    Double_t parErr[N_of_Pars_for_fit];
    Double_t chi2ndf;
};

const TString Input_File = "/home/kirill/root-on-vs-code/Femto_input/out_7p7gev_bes_1_FULL.root";
const TString Output_params_File = "/home/kirill/root-on-vs-code/Femto_output/26_08_Plus_Minus_One_article/out_7p7gev_bes_1_FULL_3D_Fit_params_Plus_Min_One.root";
const TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_output/26_08_Plus_Minus_One_article/";
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
    TMath::Exp(-(par[2]*par[2] * x[0] * x[0] + par[3]*par[3]*x[1]*x[1] + par[4]*par[4] * x[2]*x[2]
    ) / (0.197327 * 0.197327)));

    return fitval;
}


//3D map of q_inv !!!fill it in code!!!
TH3F* h_Map_3D[N_Bins_Centr][N_Bins_Kt];
Int_t cur_iCent = -1;
Int_t cur_iKt = -1;

Double_t fitf_Coulomb(Double_t *x, Double_t *par)
{
    Double_t q_inv_from_B_w = h_Map_3D[cur_iCent][cur_iKt]->GetBinContent(h_Map_3D[cur_iCent][cur_iKt]->FindBin(x[0],x[1],x[2]));
    Double_t K = hQCoul->GetBinContent(hQCoul->FindBin(q_inv_from_B_w));
    Double_t fitval = par[0]  * ((1. - par[1]) + par[1] * K *(1.+
    TMath::Exp(-(par[2]*par[2] * x[0] * x[0] + par[3]*par[3]*x[1]*x[1] + par[4]*par[4] * x[2]*x[2]
    ) / (0.197327 * 0.197327))));

    return fitval;
}

void CF_3D_Fit_Pi_plus_minus_one()
{
TFile *f = TFile::Open(Input_File, "READ");
if (!f || f->IsZombie()) {std::cout << "Ohhhhh" << std::endl; exit(0); }
//osl 3D: 
TString hist_3D_Name = "h_3D";
TString hist_3D_Plus_Minus_One_Name = "h_3D_Plus_Minus_One";

TH3F* h_Arr_3D[N_hist_types_3D][N_Charge][N_Bins_Centr][N_Bins_Kt];
TH3F* h_Arr_3D_Plus_Minus_One[N_hist_types_3D][N_Bins_Centr][N_Bins_Kt];

for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
{
    for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
    {
        for (Int_t ihType = 0; ihType < N_hist_types_3D; ihType++)
        {
            for (Int_t iCh = 0; iCh < N_Charge; iCh++)
            {
                // 3D:
                h_Arr_3D[ihType][iCh][iCent][iKt] = (TH3F *)f->Get(Form("%s_%i_%i_%i_%i", hist_3D_Name.Data(), ihType, iCh, iCent, iKt));
            }
            //Plus Minus in One hist:
            h_Arr_3D_Plus_Minus_One[ihType][iCent][iKt] = (TH3F *)h_Arr_3D[ihType][0][iCent][iKt]->Clone(Form("%s_%i_%i_%i", hist_3D_Plus_Minus_One_Name.Data(), ihType, iCent, iKt));
            h_Arr_3D_Plus_Minus_One[ihType][iCent][iKt]->Add(h_Arr_3D[ihType][1][iCent][iKt]);
        }
    }
}

// 1) Calculating Norm. coeff:
const Double_t q_i_0 = -0.195;
const Double_t q_i_1 = -0.095;
const Double_t q_i_2 =  0.095;
const Double_t q_i_3 =  0.195;

//A and B hists must be with same binning via all directions - out,side,long
Double_t q_i_min_max_bins[4][N_Bins_Centr][N_Bins_Kt][3] = {0.0};
// ----[q_min_max[0], q_min_max[1]]---0---[q_min_max[2], q_min_max[3]]----

    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {

            q_i_min_max_bins[0][iCent][iKt][0] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetXaxis()->FindBin(q_i_0);
            q_i_min_max_bins[0][iCent][iKt][1] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetYaxis()->FindBin(q_i_0);
            q_i_min_max_bins[0][iCent][iKt][2] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetZaxis()->FindBin(q_i_0);

            q_i_min_max_bins[1][iCent][iKt][0] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetXaxis()->FindBin(q_i_1);
            q_i_min_max_bins[1][iCent][iKt][1] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetYaxis()->FindBin(q_i_1);
            q_i_min_max_bins[1][iCent][iKt][2] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetZaxis()->FindBin(q_i_1);

            q_i_min_max_bins[2][iCent][iKt][0] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetXaxis()->FindBin(q_i_2);
            q_i_min_max_bins[2][iCent][iKt][1] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetYaxis()->FindBin(q_i_2);
            q_i_min_max_bins[2][iCent][iKt][2] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetZaxis()->FindBin(q_i_2);

            q_i_min_max_bins[3][iCent][iKt][0] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetXaxis()->FindBin(q_i_3);
            q_i_min_max_bins[3][iCent][iKt][1] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetYaxis()->FindBin(q_i_3);
            q_i_min_max_bins[3][iCent][iKt][2] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetZaxis()->FindBin(q_i_3);

        }
    }

Double_t A_integrals_Plus_Minus_One_arr[N_Bins_Centr][N_Bins_Kt] = {0.0};
Double_t B_integrals_Plus_Minus_One_arr[N_Bins_Centr][N_Bins_Kt] = {0.0};
Double_t Norm_coeff_Plus_Minus_One_arr[N_Bins_Centr][N_Bins_Kt] = {0.0}; // Norm_coeff = B/A

// Double_t Integral_of_cube_NO_center = 0.0; //For denumerator of fit function

// TH3F* h_Hist_with_ones = (TH3F*)h_Arr_3D[1][0][0][0]->Clone("h_Hist_with_ones");

// //Fill hist with ONES:
// for(Int_t i = 1;i<=h_Hist_with_ones->GetNbinsX();i++)
// {
//     for(Int_t j = 1;j<=h_Hist_with_ones->GetNbinsY();j++)
//     {
//         for(Int_t k = 1;k<=h_Hist_with_ones->GetNbinsZ();k++)
//         {
//             h_Hist_with_ones->SetBinContent(i,j,k,1.0);
//         }
//     }
// }
// //Count Integral for B fit normalize:
// Integral_of_cube_NO_center= h_Hist_with_ones->Integral(
//                 q_i_min_max_bins[0][0][0][0][0], q_i_min_max_bins[3][0][0][0][0],
//                 q_i_min_max_bins[0][0][0][0][1], q_i_min_max_bins[3][0][0][0][1],
//                 q_i_min_max_bins[0][0][0][0][2], q_i_min_max_bins[3][0][0][0][2]) 

//                 - h_Hist_with_ones->Integral(
//                 q_i_min_max_bins[1][0][0][0][0], q_i_min_max_bins[2][0][0][0][0],
//                 q_i_min_max_bins[1][0][0][0][1], q_i_min_max_bins[2][0][0][0][1],
//                 q_i_min_max_bins[1][0][0][0][2], q_i_min_max_bins[2][0][0][0][2]
//                 );          
// std::cout<< Integral_of_cube_NO_center <<std::endl;

for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
{
    for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
    {
        //For Plus Minus One:
        A_integrals_Plus_Minus_One_arr[iCent][iKt] = h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->Integral(
                q_i_min_max_bins[0][iCent][iKt][0], q_i_min_max_bins[3][iCent][iKt][0],
                q_i_min_max_bins[0][iCent][iKt][1], q_i_min_max_bins[3][iCent][iKt][1],
                q_i_min_max_bins[0][iCent][iKt][2], q_i_min_max_bins[3][iCent][iKt][2])

            - h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->Integral(
                q_i_min_max_bins[1][iCent][iKt][0], q_i_min_max_bins[2][iCent][iKt][0],
                q_i_min_max_bins[1][iCent][iKt][1], q_i_min_max_bins[2][iCent][iKt][1],
                q_i_min_max_bins[1][iCent][iKt][2], q_i_min_max_bins[2][iCent][iKt][2]);

            B_integrals_Plus_Minus_One_arr[iCent][iKt] = h_Arr_3D_Plus_Minus_One[1][iCent][iKt]->Integral(
                q_i_min_max_bins[0][iCent][iKt][0], q_i_min_max_bins[3][iCent][iKt][0],
                q_i_min_max_bins[0][iCent][iKt][1], q_i_min_max_bins[3][iCent][iKt][1],
                q_i_min_max_bins[0][iCent][iKt][2], q_i_min_max_bins[3][iCent][iKt][2])

            - h_Arr_3D_Plus_Minus_One[1][iCent][iKt]->Integral(
                q_i_min_max_bins[1][iCent][iKt][0], q_i_min_max_bins[2][iCent][iKt][0],
                q_i_min_max_bins[1][iCent][iKt][1], q_i_min_max_bins[2][iCent][iKt][1],
                q_i_min_max_bins[1][iCent][iKt][2], q_i_min_max_bins[2][iCent][iKt][2]);

            Norm_coeff_Plus_Minus_One_arr[iCent][iKt] = B_integrals_Plus_Minus_One_arr[iCent][iKt] / A_integrals_Plus_Minus_One_arr[iCent][iKt];
    }
}
//Scaling 3D A in CF:

for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
{
    for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
    {
        h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->Scale(Norm_coeff_Plus_Minus_One_arr[iCent][iKt]);
    }
}

// //Deviding A/B:

for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
{
    for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
    {
        h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->Divide(h_Arr_3D_Plus_Minus_One[1][iCent][iKt]);
    }
}

//Now let's fit:
//Fill array of K(q_inv):
for (Int_t i = 1; i <= N_bins_K_coul; i++)
{
    K_values_arr[i - 1] = hQCoul->GetBinContent(i);
}


//Fit params for Plus Minus One:
Double_t par_No_Coul[N_Bins_Centr][N_Bins_Kt][N_of_Pars_for_fit];
Double_t par_Coul[N_Bins_Centr][N_Bins_Kt][N_of_Pars_for_fit];

Double_t par_No_Coul_Err[N_Bins_Centr][N_Bins_Kt][N_of_Pars_for_fit];
Double_t par_Coul_Err[N_Bins_Centr][N_Bins_Kt][N_of_Pars_for_fit];

Double_t par_No_Coul_Chi_div_NDF[N_Bins_Centr][N_Bins_Kt];
Double_t par_Coul_Chi_div_NDF[N_Bins_Centr][N_Bins_Kt];

TF3* f_No_Coul_Arr[N_Bins_Centr][N_Bins_Kt];
TF3* f_Coul_Arr[N_Bins_Centr][N_Bins_Kt];

const Double_t x1_fit_range = -0.20;
const Double_t x2_fit_range = 0.20;
const Double_t N_start = 1.;
const Double_t Lambda_start = 0.3;
const Double_t R_o_start = 5.;
const Double_t R_s_start = 5.;
const Double_t R_l_start = 5.;

const Double_t R_os_start = 4.0;
const Double_t R_ol_start = 4.0;

for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
{
    for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
    {

        h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetXaxis()->SetRangeUser(x1_fit_range, x2_fit_range);
        h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetYaxis()->SetRangeUser(x1_fit_range, x2_fit_range);
        h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetZaxis()->SetRangeUser(x1_fit_range, x2_fit_range);

        TString fitf_name = TString::Format("f_No_Coul_%d_%d", iCent, iKt);
        f_No_Coul_Arr[iCent][iKt] = new TF3(fitf_name, fitf_No_Coulomb, x1_fit_range, x2_fit_range,
                                                 x1_fit_range, x2_fit_range,
                                                 x1_fit_range, x2_fit_range, N_of_Pars_for_fit);

        f_No_Coul_Arr[iCent][iKt]->SetParameter(0, N_start);
        f_No_Coul_Arr[iCent][iKt]->SetParameter(1, Lambda_start);

        f_No_Coul_Arr[iCent][iKt]->SetParameter(2, R_o_start);
        f_No_Coul_Arr[iCent][iKt]->SetParameter(3, R_s_start);
        f_No_Coul_Arr[iCent][iKt]->SetParameter(4, R_l_start);

        // f_No_Coul_Arr[iCh][iCent][iKt]->SetParameter(5, R_os_start);
        // f_No_Coul_Arr[iCh][iCent][iKt]->SetParameter(6, R_ol_start);

        f_No_Coul_Arr[iCent][iKt]->SetParLimits(0, 0., 2.);
        f_No_Coul_Arr[iCent][iKt]->SetParLimits(1, 0., 2.);

        f_No_Coul_Arr[iCent][iKt]->SetParLimits(2, 0.01, 15.);
        f_No_Coul_Arr[iCent][iKt]->SetParLimits(3, 0.01, 15.);
        f_No_Coul_Arr[iCent][iKt]->SetParLimits(4, 0.01, 15.);

        // f_No_Coul_Arr[iCh][iCent][iKt]->SetParLimits(5, 0.01, 15.);
        // f_No_Coul_Arr[iCh][iCent][iKt]->SetParLimits(6, 0.01, 15.);

        f_No_Coul_Arr[iCent][iKt]->SetNpx(20);
        h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->Fit(fitf_name);

        f_No_Coul_Arr[iCent][iKt]->GetParameters(par_No_Coul[iCent][iKt]);

        const Double_t *errs = f_No_Coul_Arr[iCent][iKt]->GetParErrors();
        // Copying errors:
        for (Int_t p = 0; p < N_of_Pars_for_fit; ++p)
        {
            par_No_Coul_Err[iCent][iKt][p] = errs[p];
        }
        par_No_Coul_Chi_div_NDF[iCent][iKt] = f_No_Coul_Arr[iCent][iKt]->GetChisquare() / f_No_Coul_Arr[iCent][iKt]->GetNDF();
    }
}

for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
{
    for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
    {

        h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetXaxis()->SetRangeUser(x1_fit_range, x2_fit_range);
        h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetYaxis()->SetRangeUser(x1_fit_range, x2_fit_range);
        h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->GetZaxis()->SetRangeUser(x1_fit_range, x2_fit_range);

        TString fitf_name = TString::Format("f_Coul_%d_%d", iCent, iKt);
        f_Coul_Arr[iCent][iKt] = new TF3(fitf_name, fitf_Coulomb, x1_fit_range, x2_fit_range,
                                              x1_fit_range, x2_fit_range,
                                              x1_fit_range, x2_fit_range, N_of_Pars_for_fit);

        f_Coul_Arr[iCent][iKt]->SetParameter(0, N_start);
        f_Coul_Arr[iCent][iKt]->SetParameter(1, Lambda_start);

        f_Coul_Arr[iCent][iKt]->SetParameter(2, R_o_start);
        f_Coul_Arr[iCent][iKt]->SetParameter(3, R_s_start);
        f_Coul_Arr[iCent][iKt]->SetParameter(4, R_l_start);

        // f_Coul_Arr[iCh][iCent][iKt]->SetParameter(5, R_os_start);
        // f_Coul_Arr[iCh][iCent][iKt]->SetParameter(6, R_ol_start);

        f_Coul_Arr[iCent][iKt]->SetParLimits(0, 0., 2.);
        f_Coul_Arr[iCent][iKt]->SetParLimits(1, 0., 2.);

        f_Coul_Arr[iCent][iKt]->SetParLimits(2, 0.01, 15.);
        f_Coul_Arr[iCent][iKt]->SetParLimits(3, 0.01, 15.);
        f_Coul_Arr[iCent][iKt]->SetParLimits(4, 0.01, 15.);

        // f_Coul_Arr[iCh][iCent][iKt]->SetParLimits(5, 0.01, 15.);
        // f_Coul_Arr[iCh][iCent][iKt]->SetParLimits(6, 0.01, 15.);

        //!!! fill Map hist !!!
        h_Map_3D[iCent][iKt] = (TH3F *)h_Arr_3D_Plus_Minus_One[2][iCent][iKt]->Clone(TString::Format("h_Map_3D_%d_%d", iCent, iKt));
        h_Map_3D[iCent][iKt]->Divide(h_Arr_3D_Plus_Minus_One[1][iCent][iKt]);
        cur_iCent = iCent;
        cur_iKt = iKt;

        f_Coul_Arr[iCent][iKt]->SetNpx(20);
        h_Arr_3D_Plus_Minus_One[0][iCent][iKt]->Fit(fitf_name);

        f_Coul_Arr[iCent][iKt]->GetParameters(par_Coul[iCent][iKt]);

        const Double_t *errs = f_Coul_Arr[iCent][iKt]->GetParErrors();
        // Copying errors:
        for (Int_t p = 0; p < N_of_Pars_for_fit; ++p)
        {
            par_Coul_Err[iCent][iKt][p] = errs[p];
        }

        par_Coul_Chi_div_NDF[iCent][iKt] = f_Coul_Arr[iCent][iKt]->GetChisquare() / f_Coul_Arr[iCent][iKt]->GetNDF();
    }
}

//Let's Fill fit parameters:
TFile *f_params = new TFile(Output_params_File,"RECREATE");
TTree *tree_params = new TTree("fit_params_Tree", "Tree of fit parameters");


FitResult res;
tree_params->Branch("centrality",&res.centrality,"centrality/I");
tree_params->Branch("Kt",        &res.Kt,        "Kt/I");
tree_params->Branch("fitType",   &res.fitType,   "fitType/I");
tree_params->Branch("par",    res.par,    Form("par[%d]/D", N_of_Pars_for_fit));
tree_params->Branch("parErr", res.parErr, Form("parErr[%d]/D", N_of_Pars_for_fit));
tree_params->Branch("chi2ndf",   &res.chi2ndf,   "chi2ndf/D");

for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
{
    for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
    {

        res.centrality = iCent;
        res.Kt = iKt;
        res.fitType = 0; // No coulomb
        for (int p = 0; p < N_of_Pars_for_fit; p++)
        {
            res.par[p] = par_No_Coul[iCent][iKt][p];
            res.parErr[p] = par_No_Coul_Err[iCent][iKt][p];
        }
        res.chi2ndf = par_No_Coul_Chi_div_NDF[iCent][iKt];
        tree_params->Fill();

        res.centrality = iCent;
        res.Kt = iKt;
        res.fitType = 1;
        for (int p = 0; p < N_of_Pars_for_fit; p++)
        {
            res.par[p] = par_Coul[iCent][iKt][p];
            res.parErr[p] = par_Coul_Err[iCent][iKt][p];
        }
        res.chi2ndf = par_Coul_Chi_div_NDF[iCent][iKt];
        tree_params->Fill();
    }
}

tree_params->Write();
f_params->Close();

for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
{
    for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
    {
        delete f_No_Coul_Arr[iCent][iKt];
        f_No_Coul_Arr[iCent][iKt] = nullptr;

        delete f_Coul_Arr[iCent][iKt];
        f_Coul_Arr[iCent][iKt] = nullptr;
    }
}

f_K->Close();
delete f_K;
f->Close();
delete f;
}