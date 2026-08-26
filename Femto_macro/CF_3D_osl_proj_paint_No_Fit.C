//run this macro after CF_3D_osl_proj.C
//synchronize IO files
//this macro creates .root file with projections of 3D_CF for each [charge][Centrality][K_t][projection]
//this macro require .root file with fit data from CF_3D_osl_proj
//also this .root file will contain 3D_CF & B_weighted for each [charge][Centrality][K_t] 

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
    Int_t    charge;
    Int_t    centrality;
    Int_t    Kt;
    Int_t    fitType;     // 0 = NoCoul, 1 = Coul
    Double_t par[N_of_Pars_for_fit];
    Double_t parErr[N_of_Pars_for_fit];
    Double_t chi2ndf;
};

const TString Input_File = "/home/kirill/root-on-vs-code/Femto_input/out_7p7gev_bes_1_FULL.root";
const TString Output_File = "/home/kirill/root-on-vs-code/Femto_output/26_08_FMR_article/out_7p7gev_bes_1_FULL_FMR_pr_3D_proj_From_paint.root";
const TString Input_params_File = "/home/kirill/root-on-vs-code/Femto_output/26_08_FMR_article/out_7p7gev_bes_1_FULL_FMR_3D_Fit_params.root";
const TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_output/26_08_FMR_article/";
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
TH3F* h_Map_3D[N_Charge][N_Bins_Centr][N_Bins_Kt];
Int_t cur_iCh = -1;
Int_t cur_iCent = -1;
Int_t cur_iKt = -1;

Double_t fitf_Coulomb(Double_t *x, Double_t *par)
{
    Double_t q_inv_from_B_w = h_Map_3D[cur_iCh][cur_iCent][cur_iKt]->GetBinContent(h_Map_3D[cur_iCh][cur_iCent][cur_iKt]->FindBin(x[0],x[1],x[2]));
    Double_t K = hQCoul->GetBinContent(hQCoul->FindBin(q_inv_from_B_w));
    Double_t fitval = par[0]  * ((1. - par[1]) + par[1] * K *(1.+
    TMath::Exp(-(par[2]*par[2] * x[0] * x[0] + par[3]*par[3]*x[1]*x[1] + par[4]*par[4] * x[2]*x[2]
    ) / (0.197327 * 0.197327))));

    return fitval;
}

void CF_3D_osl_proj_paint_No_Fit()
{
TFile *f = TFile::Open(Input_File, "READ");
if (!f || f->IsZombie()) {std::cout << "Ohhhhh" << std::endl; exit(0); }
//osl 3D: 
TString hist_3D_Name = "h_3D";

TH3F* h_Arr_3D[N_hist_types_3D][N_Charge][N_Bins_Centr][N_Bins_Kt];
Double_t Entries_in_A_arr[N_Charge][N_Bins_Centr][N_Bins_Kt];

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
            Entries_in_A_arr[iCh][iCent][iKt] = h_Arr_3D[0][iCh][iCent][iKt]->GetEntries();
        }
    }
}

// 1) Calculating Norm. coeff:
const Double_t q_i_0 = -0.195;
const Double_t q_i_1 = -0.095;
const Double_t q_i_2 =  0.095;
const Double_t q_i_3 =  0.195;

//A and B hists must be with same binning via all directions - out,side,long
Double_t q_i_min_max_bins[4][N_Charge][N_Bins_Centr][N_Bins_Kt][3] = {0.0};
// ----[q_min_max[0], q_min_max[1]]---0---[q_min_max[2], q_min_max[3]]----

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {

            q_i_min_max_bins[0][iCh][iCent][iKt][0] = h_Arr_3D[0][iCh][iCent][iKt]->GetXaxis()->FindBin(q_i_0);
            q_i_min_max_bins[0][iCh][iCent][iKt][1] = h_Arr_3D[0][iCh][iCent][iKt]->GetYaxis()->FindBin(q_i_0);
            q_i_min_max_bins[0][iCh][iCent][iKt][2] = h_Arr_3D[0][iCh][iCent][iKt]->GetZaxis()->FindBin(q_i_0);


            q_i_min_max_bins[1][iCh][iCent][iKt][0] = h_Arr_3D[0][iCh][iCent][iKt]->GetXaxis()->FindBin(q_i_1);
            q_i_min_max_bins[1][iCh][iCent][iKt][1] = h_Arr_3D[0][iCh][iCent][iKt]->GetYaxis()->FindBin(q_i_1);
            q_i_min_max_bins[1][iCh][iCent][iKt][2] = h_Arr_3D[0][iCh][iCent][iKt]->GetZaxis()->FindBin(q_i_1);

            q_i_min_max_bins[2][iCh][iCent][iKt][0] = h_Arr_3D[0][iCh][iCent][iKt]->GetXaxis()->FindBin(q_i_2);
            q_i_min_max_bins[2][iCh][iCent][iKt][1] = h_Arr_3D[0][iCh][iCent][iKt]->GetYaxis()->FindBin(q_i_2);
            q_i_min_max_bins[2][iCh][iCent][iKt][2] = h_Arr_3D[0][iCh][iCent][iKt]->GetZaxis()->FindBin(q_i_2);

            q_i_min_max_bins[3][iCh][iCent][iKt][0] = h_Arr_3D[0][iCh][iCent][iKt]->GetXaxis()->FindBin(q_i_3);
            q_i_min_max_bins[3][iCh][iCent][iKt][1] = h_Arr_3D[0][iCh][iCent][iKt]->GetYaxis()->FindBin(q_i_3);
            q_i_min_max_bins[3][iCh][iCent][iKt][2] = h_Arr_3D[0][iCh][iCent][iKt]->GetZaxis()->FindBin(q_i_3);

        }
    }
}

Double_t A_integrals_arr[N_Charge][N_Bins_Centr][N_Bins_Kt] = {0.0};
Double_t B_integrals_arr[N_Charge][N_Bins_Centr][N_Bins_Kt] = {0.0};
Double_t Norm_coeff_arr[N_Charge][N_Bins_Centr][N_Bins_Kt] = {0.0}; // Norm_coeff = B/A

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

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {

            A_integrals_arr[iCh][iCent][iKt]= h_Arr_3D[0][iCh][iCent][iKt]->Integral(
                q_i_min_max_bins[0][iCh][iCent][iKt][0], q_i_min_max_bins[3][iCh][iCent][iKt][0],
                q_i_min_max_bins[0][iCh][iCent][iKt][1], q_i_min_max_bins[3][iCh][iCent][iKt][1],
                q_i_min_max_bins[0][iCh][iCent][iKt][2], q_i_min_max_bins[3][iCh][iCent][iKt][2]) 

                - h_Arr_3D[0][iCh][iCent][iKt]->Integral(
                q_i_min_max_bins[1][iCh][iCent][iKt][0], q_i_min_max_bins[2][iCh][iCent][iKt][0],
                q_i_min_max_bins[1][iCh][iCent][iKt][1], q_i_min_max_bins[2][iCh][iCent][iKt][1],
                q_i_min_max_bins[1][iCh][iCent][iKt][2], q_i_min_max_bins[2][iCh][iCent][iKt][2]
                );
            B_integrals_arr[iCh][iCent][iKt]= h_Arr_3D[1][iCh][iCent][iKt]->Integral(
                q_i_min_max_bins[0][iCh][iCent][iKt][0], q_i_min_max_bins[3][iCh][iCent][iKt][0],
                q_i_min_max_bins[0][iCh][iCent][iKt][1], q_i_min_max_bins[3][iCh][iCent][iKt][1],
                q_i_min_max_bins[0][iCh][iCent][iKt][2], q_i_min_max_bins[3][iCh][iCent][iKt][2]) 

                - h_Arr_3D[1][iCh][iCent][iKt]->Integral(
                q_i_min_max_bins[1][iCh][iCent][iKt][0], q_i_min_max_bins[2][iCh][iCent][iKt][0],
                q_i_min_max_bins[1][iCh][iCent][iKt][1], q_i_min_max_bins[2][iCh][iCent][iKt][1],
                q_i_min_max_bins[1][iCh][iCent][iKt][2], q_i_min_max_bins[2][iCh][iCent][iKt][2]
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
const Double_t range = 0.045;//Gev/c              
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
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][0]->SetName(Form("CF_3D_proj_%i_%i_%i_%i_%s",ihType,iCh,iCent,iKt,OSL_names[0].Data()));
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetYaxis()->SetRange();
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetZaxis()->SetRange();

                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(-range,range);
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][1] = (TH1F*)h_Arr_3D[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[1]);
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][1]->SetName(Form("CF_3D_proj_%i_%i_%i_%i_%s",ihType,iCh,iCent,iKt,OSL_names[1].Data()));
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetXaxis()->SetRange();
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetZaxis()->SetRange();

                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(-range,range);
                    h_Arr_3D[ihType][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(-range,range);
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][2] = (TH1F*)h_Arr_3D[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[2]);
                    h_Arr_3D_Projects_OSL[ihType][iCh][iCent][iKt][2]->SetName(Form("CF_3D_proj_%i_%i_%i_%i_%s",ihType,iCh,iCent,iKt,OSL_names[2].Data()));
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


//Fit params from params file:
Double_t par_No_Coul[N_Charge][N_Bins_Centr][N_Bins_Kt][N_of_Pars_for_fit];
Double_t par_Coul[N_Charge][N_Bins_Centr][N_Bins_Kt][N_of_Pars_for_fit];

Double_t par_No_Coul_Err[N_Charge][N_Bins_Centr][N_Bins_Kt][N_of_Pars_for_fit];
Double_t par_Coul_Err[N_Charge][N_Bins_Centr][N_Bins_Kt][N_of_Pars_for_fit];

Double_t par_No_Coul_Chi_div_NDF[N_Charge][N_Bins_Centr][N_Bins_Kt];
Double_t par_Coul_Chi_div_NDF[N_Charge][N_Bins_Centr][N_Bins_Kt];

TF3* f_No_Coul_Arr[N_Charge][N_Bins_Centr][N_Bins_Kt];
TF3* f_Coul_Arr[N_Charge][N_Bins_Centr][N_Bins_Kt];

// Open input file and get tree
TFile *fin_par = TFile::Open(Input_params_File, "READ");
if (!fin_par || fin_par->IsZombie()) {
    std::cerr << "Cannot open input params file " << Input_params_File << std::endl;
    return;
}
TTree *tree_of_params = (TTree*)fin_par->Get("fit_params_Tree");
if (!tree_of_params) {
    std::cerr << "Tree 'fit_params_Tree' not found!" << std::endl;
    fin_par->Close();
    return;
}
// Read tree and fill arrays for chosen fitType

FitResult res;
tree_of_params->SetBranchAddress("charge",     &res.charge);
tree_of_params->SetBranchAddress("centrality", &res.centrality);
tree_of_params->SetBranchAddress("Kt",         &res.Kt);
tree_of_params->SetBranchAddress("fitType",    &res.fitType);
tree_of_params->SetBranchAddress("par",        res.par);
tree_of_params->SetBranchAddress("parErr",     res.parErr);
tree_of_params->SetBranchAddress("chi2ndf",    &res.chi2ndf);

Long64_t nEntries = tree_of_params->GetEntries();
for (Long64_t i = 0; i < nEntries; i++)
{
    tree_of_params->GetEntry(i);
    if (res.fitType == 0)
    {   
        for(int i=0;i<N_of_Pars_for_fit;i++)
        {
        par_No_Coul[res.charge][res.centrality][res.Kt][i] = res.par[i]; 
        par_No_Coul_Err[res.charge][res.centrality][res.Kt][i] = res.parErr[i];
        }
        par_No_Coul_Chi_div_NDF[res.charge][res.centrality][res.Kt] = res.chi2ndf; 
    }
    else if (res.fitType == 1)
    {
        for(int i=0;i<N_of_Pars_for_fit;i++)
        {
        par_Coul[res.charge][res.centrality][res.Kt][i] = res.par[i]; 
        par_Coul_Err[res.charge][res.centrality][res.Kt][i] = res.parErr[i];
        }
        par_Coul_Chi_div_NDF[res.charge][res.centrality][res.Kt] = res.chi2ndf; 
    }
    
}

fin_par->Close();


const Double_t x1_fit_range = -0.20;
const Double_t x2_fit_range = 0.20;
// const Double_t N_start = 1.;
// const Double_t Lambda_start = 0.3;
// const Double_t R_o_start = 5.;
// const Double_t R_s_start = 5.;
// const Double_t R_l_start = 5.;

// const Double_t R_os_start = 4.0;
// const Double_t R_ol_start = 4.0;

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
                                                                            x1_fit_range, x2_fit_range, N_of_Pars_for_fit);

            f_No_Coul_Arr[iCh][iCent][iKt]->SetParameters(par_No_Coul[iCh][iCent][iKt]);
            f_No_Coul_Arr[iCh][iCent][iKt]->SetParErrors(par_No_Coul_Err[iCh][iCent][iKt]);
            //f_No_Coul_Arr[iCh][iCent][iKt]->SetChisquare

            f_No_Coul_Arr[iCh][iCent][iKt]->SetNpx(20);

}}}

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {

            h_Arr_3D[0][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(x1_fit_range, x2_fit_range);
            h_Arr_3D[0][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(x1_fit_range, x2_fit_range);
            h_Arr_3D[0][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(x1_fit_range, x2_fit_range);

            TString fitf_name = TString::Format("f_Coul_%d_%d_%d", iCh, iCent, iKt);
            f_Coul_Arr[iCh][iCent][iKt]= new TF3(fitf_name, fitf_Coulomb,  x1_fit_range, x2_fit_range, 
                                                                            x1_fit_range, x2_fit_range, 
                                                                            x1_fit_range, x2_fit_range, N_of_Pars_for_fit);
            
            //!!! fill Map hist !!!
            h_Map_3D[iCh][iCent][iKt] = (TH3F*)h_Arr_3D[2][iCh][iCent][iKt]->Clone(TString::Format("h_Map_3D_%d_%d_%d",iCh,iCent,iKt));
            h_Map_3D[iCh][iCent][iKt]->Divide(h_Arr_3D[1][iCh][iCent][iKt]);
            cur_iCh=iCh;
            cur_iCent=iCent;
            cur_iKt=iKt;

            f_Coul_Arr[iCh][iCent][iKt]->SetParameters(par_Coul[iCh][iCent][iKt]);
            f_Coul_Arr[iCh][iCent][iKt]->SetParErrors(par_Coul_Err[iCh][iCent][iKt]);

            f_Coul_Arr[iCh][iCent][iKt]->SetNpx(20);

}}}


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
            TString name_hA_Fit_No_Coul = TString::Format("hA_Fit_No_Coul_%d_%d_%d",iCh,iCent,iKt);
            TString name_hB_Fit_No_Coul = TString::Format("hB_Fit_No_Coul_%d_%d_%d",iCh,iCent,iKt);
            h_Fit_No_Coul[0][iCh][iCent][iKt] = (TH3F*)h_Arr_3D[0][iCh][iCent][iKt]->Clone(name_hA_Fit_No_Coul);
            h_Fit_No_Coul[1][iCh][iCent][iKt] = (TH3F*)h_Arr_3D[0][iCh][iCent][iKt]->Clone(name_hB_Fit_No_Coul);

                    
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
            //Here we must du Cf = B*Cf/B
            //A <-> B*Cf:
            Double_t B_value_for_FitFunc = h_Arr_3D[1][iCh][iCent][iKt]->GetBinContent(bin);
            h_Fit_No_Coul[0][iCh][iCent][iKt]->SetBinContent(bin, 
                                              B_value_for_FitFunc*f_No_Coul_Arr[iCh][iCent][iKt]->Eval(x_arr[0],x_arr[1],x_arr[2]));
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
            h_Fit_No_Coul[1][iCh][iCent][iKt]->SetBinContent(i,j,k,h_Arr_3D[1][iCh][iCent][iKt]->GetBinContent(i,j,k));
        }
    }
}

//For Divide fit function A/B:
for (Int_t ihType = 0; ihType < 2; ihType++)
{
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][0] = (TH1F *)h_Fit_No_Coul[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[0]);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][0]->SetName(Form("CF_3D_norm_No_Coul_proj_%i_%i_%i_%i_%s", ihType, iCh, iCent, iKt, OSL_names[0].Data()));
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRange();
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRange();

    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][1] = (TH1F *)h_Fit_No_Coul[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[1]);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][1]->SetName(Form("CF_3D_norm_No_Coul_proj_%i_%i_%i_%i_%s", ihType, iCh, iCent, iKt, OSL_names[1].Data()));
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRange();
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRange();

    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(-range, range);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][2] = (TH1F *)h_Fit_No_Coul[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[2]);
    h_Fit_No_Coul_Proj[ihType][iCh][iCent][iKt][2]->SetName(Form("CF_3D_norm_No_Coul_proj_%i_%i_%i_%i_%s", ihType, iCh, iCent, iKt, OSL_names[2].Data()));
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
TH3F* h_Fit_Coul[2][N_Charge][N_Bins_Centr][N_Bins_Kt];
TH1F* h_Fit_Coul_Proj[2][N_Charge][N_Bins_Centr][N_Bins_Kt][3];

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            TString name_hA_Fit_Coul = TString::Format("hA_Fit_Coul_%d_%d_%d",iCh,iCent,iKt);
            TString name_hB_Fit_Coul = TString::Format("hB_Fit_Coul_%d_%d_%d",iCh,iCent,iKt);
            h_Fit_Coul[0][iCh][iCent][iKt] = (TH3F*)h_Arr_3D[0][iCh][iCent][iKt] ->Clone(name_hA_Fit_Coul);
            h_Fit_Coul[1][iCh][iCent][iKt] = (TH3F*)h_Arr_3D[0][iCh][iCent][iKt] ->Clone(name_hB_Fit_Coul);

//Fill A hist:
double x_arr[3];
for(Int_t i = 1;i<=h_Fit_Coul[0][iCh][iCent][iKt]->GetNbinsX();i++)
{
    x_arr[0]=h_Fit_Coul[0][iCh][iCent][iKt]->GetXaxis()->GetBinCenter(i);
    for(Int_t j = 1;j<=h_Fit_Coul[0][iCh][iCent][iKt]->GetNbinsY();j++)
    {
        x_arr[1]=h_Fit_Coul[0][iCh][iCent][iKt]->GetYaxis()->GetBinCenter(j);
        for(Int_t k = 1;k<=h_Fit_Coul[0][iCh][iCent][iKt]->GetNbinsZ();k++)
        {
            x_arr[2]=h_Fit_Coul[0][iCh][iCent][iKt]->GetZaxis()->GetBinCenter(k);

            Int_t bin = h_Fit_Coul[0][iCh][iCent][iKt]->GetBin(i,j,k);
            //WARNING global variables for evaluation of fit function:
            cur_iCh=iCh;
            cur_iCent=iCent;
            cur_iKt=iKt;
            Double_t B_value_for_FitFunc = h_Arr_3D[1][iCh][iCent][iKt]->GetBinContent(bin);
            h_Fit_Coul[0][iCh][iCent][iKt]->SetBinContent(bin, 
                                            B_value_for_FitFunc*f_Coul_Arr[iCh][iCent][iKt]->Eval(x_arr[0],x_arr[1],x_arr[2]));
        }
    }
}

//Fill B hist:
for(Int_t i = 1;i<=h_Fit_Coul[1][iCh][iCent][iKt]->GetNbinsX();i++)
{
    for(Int_t j = 1;j<=h_Fit_Coul[1][iCh][iCent][iKt]->GetNbinsY();j++)
    {
        for(Int_t k = 1;k<=h_Fit_Coul[1][iCh][iCent][iKt]->GetNbinsZ();k++)
        {
            h_Fit_Coul[1][iCh][iCent][iKt]->SetBinContent(i,j,k,h_Arr_3D[1][iCh][iCent][iKt]->GetBinContent(i,j,k));
        }
    }
}

//For Divide fit function A/B:
for (Int_t ihType = 0; ihType < 2; ihType++)
{
    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(-range, range);
    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(-range, range);
    h_Fit_Coul_Proj[ihType][iCh][iCent][iKt][0] = (TH1F *)h_Fit_Coul[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[0]);
    h_Fit_Coul_Proj[ihType][iCh][iCent][iKt][0]->SetName(Form("CF_3D_norm_Coul_proj_%i_%i_%i_%i_%s", ihType, iCh, iCent, iKt, OSL_names[0].Data()));
    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRange();
    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRange();

    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(-range, range);
    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRangeUser(-range, range);
    h_Fit_Coul_Proj[ihType][iCh][iCent][iKt][1] = (TH1F *)h_Fit_Coul[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[1]);
    h_Fit_Coul_Proj[ihType][iCh][iCent][iKt][1]->SetName(Form("CF_3D_norm_Coul_proj_%i_%i_%i_%i_%s", ihType, iCh, iCent, iKt, OSL_names[1].Data()));
    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRange();
    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetZaxis()->SetRange();

    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRangeUser(-range, range);
    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRangeUser(-range, range);
    h_Fit_Coul_Proj[ihType][iCh][iCent][iKt][2] = (TH1F *)h_Fit_Coul[ihType][iCh][iCent][iKt]->Project3D(OSL_xyz[2]);
    h_Fit_Coul_Proj[ihType][iCh][iCent][iKt][2]->SetName(Form("CF_3D_norm_Coul_proj_%i_%i_%i_%i_%s", ihType, iCh, iCent, iKt, OSL_names[2].Data()));
    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetXaxis()->SetRange();
    h_Fit_Coul[ihType][iCh][iCent][iKt]->GetYaxis()->SetRange();
}
// //Deviding A/B:

for (Int_t iOSL = 0; iOSL < 3; iOSL++)
{
    h_Fit_Coul_Proj[0][iCh][iCent][iKt][iOSL]->Divide(h_Fit_Coul_Proj[1][iCh][iCent][iKt][iOSL]);
    h_Fit_Coul_Proj[0][iCh][iCent][iKt][iOSL]->GetXaxis()->SetRangeUser(x1_fit_range,x2_fit_range);  
}

}}}
//Drawing projections:

TCanvas* c_3D_Proj_Arr[N_Charge][N_Bins_Centr][N_Bins_Kt][3];

for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
                    for(Int_t iOSL =0;iOSL<3;iOSL++)
                    {
                        TString canv_name = TString::Format("c_3D_%d_%d_%d_%d",iCh,iCent,iKt,iOSL);
                        TString canv_title = Charged_particles_pairs_titles[iCh] +", Centr: "+ Centrality_region_titles[iCent]
                        + ", Kt: " + Kt_bins_titles[iKt] +", "+ OSL_names[iOSL]+" projection";
                        c_3D_Proj_Arr[iCh][iCent][iKt][iOSL] = new TCanvas(canv_name, canv_title,1920,1080);


                        c_3D_Proj_Arr[iCh][iCent][iKt][iOSL]->cd();
                        h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->SetTitle(c_3D_Proj_Arr[iCh][iCent][iKt][iOSL]->GetTitle());
                        h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->SetStats(0);
                        h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->GetYaxis()->SetRangeUser(0.7,1.5);
                        h_Arr_3D_Projects_OSL[0][iCh][iCent][iKt][iOSL]->Draw("E");

                        h_Fit_No_Coul_Proj[0][iCh][iCent][iKt][iOSL]->SetLineColor(kGreen);
                        h_Fit_No_Coul_Proj[0][iCh][iCent][iKt][iOSL]->Draw("P L SAME HIST");

                        h_Fit_Coul_Proj[0][iCh][iCent][iKt][iOSL]->SetLineColor(kRed);
                        h_Fit_Coul_Proj[0][iCh][iCent][iKt][iOSL]->Draw("P L SAME HIST");


                        TLegend* leg = new TLegend(0.55,0.55,0.9,0.9);
                        leg->SetTextSize(0.04);
                        leg->SetFillStyle(0);     // Transparent

                        //Entries in A:
                        leg->AddEntry((TObject*)0, Form(" Entries in A: %.3f Millions", Entries_in_A_arr[iCh][iCent][iKt]/1000000.),"");
                        // Gauss:
                        leg->AddEntry(h_Fit_No_Coul_Proj[0][iCh][iCent][iKt][iOSL], "Gauss", "F");

                        leg->AddEntry((TObject*)0, Form(" Chi2/NDF = %.3f", par_No_Coul_Chi_div_NDF[iCh][iCent][iKt]),"");

                        leg->AddEntry((TObject*)0, Form(" N=%.4f #pm %.4f", par_No_Coul[iCh][iCent][iKt][0],par_No_Coul_Err[iCh][iCent][iKt][0]),"");
                        leg->AddEntry((TObject*)0, Form(" #lambda=%.4f #pm %.4f", par_No_Coul[iCh][iCent][iKt][1],par_No_Coul_Err[iCh][iCent][iKt][1]),"");
                        leg->AddEntry((TObject*)0, Form(" R_{o}=%.3f #pm %.3f", par_No_Coul[iCh][iCent][iKt][2],par_No_Coul_Err[iCh][iCent][iKt][2]),"");
                        leg->AddEntry((TObject*)0, Form(" R_{s}=%.3f #pm %.3f", par_No_Coul[iCh][iCent][iKt][3],par_No_Coul_Err[iCh][iCent][iKt][3]),"");
                        leg->AddEntry((TObject*)0, Form(" R_{l}=%.3f #pm %.3f", par_No_Coul[iCh][iCent][iKt][4],par_No_Coul_Err[iCh][iCent][iKt][4]),"");

                            
                        // Bowler-Sinyukov:
                        leg->AddEntry(h_Fit_Coul_Proj[0][iCh][iCent][iKt][iOSL], "Bowler-Sinyukov", "F");
                            
                        leg->AddEntry((TObject*)0, Form(" Chi2/NDF = %.3f", par_Coul_Chi_div_NDF[iCh][iCent][iKt]),"");

                        leg->AddEntry((TObject*)0, Form(" N=%.4f #pm %.4f", par_Coul[iCh][iCent][iKt][0],par_Coul_Err[iCh][iCent][iKt][0]),"");
                        leg->AddEntry((TObject*)0, Form(" #lambda=%.3f #pm %.3f", par_Coul[iCh][iCent][iKt][1],par_Coul_Err[iCh][iCent][iKt][1]),"");
                        leg->AddEntry((TObject*)0, Form(" R_{o}=%.3f #pm %.3f", par_Coul[iCh][iCent][iKt][2],par_Coul_Err[iCh][iCent][iKt][2]),"");
                        leg->AddEntry((TObject*)0, Form(" R_{s}=%.3f #pm %.3f", par_Coul[iCh][iCent][iKt][3],par_Coul_Err[iCh][iCent][iKt][3]),"");
                        leg->AddEntry((TObject*)0, Form(" R_{l}=%.3f #pm %.3f", par_Coul[iCh][iCent][iKt][4],par_Coul_Err[iCh][iCent][iKt][4]),"");

                        leg->Draw();
                        
                        c_3D_Proj_Arr[iCh][iCent][iKt][iOSL]->Write();
                    }
        }
    }
}



f_out->Close();


for (Int_t iCh = 0; iCh < N_Charge; iCh++)
{
    for (Int_t iCent = 0; iCent < N_Bins_Centr; iCent++)
    {
        for (Int_t iKt = 0; iKt < N_Bins_Kt; iKt++)
        {
            delete f_No_Coul_Arr[iCh][iCent][iKt];
            f_No_Coul_Arr[iCh][iCent][iKt] = nullptr;

            delete f_Coul_Arr[iCh][iCent][iKt];
            f_Coul_Arr[iCh][iCent][iKt] = nullptr;
}}}

f_K->Close();
delete f_K;
f->Close();
delete f;
}