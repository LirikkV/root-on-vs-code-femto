// macro to read fit parameters from the tree and plot R_out, R_side, R_long vs K_t
// for Pi+Pi+ and Pi-Pi- separately, 6 canvases total, 9 centralities each.

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include <iostream>
#include <vector>

// ---------------- User settings ----------------
const TString Input_File  = "/home/kirill/root-on-vs-code/Femto_output/20_05_FMRW_article/out_Au_Au200_20_05_FMR_pr_3D_proj_params_21_05.root";
const TString Output_File = "/home/kirill/root-on-vs-code/Femto_output/20_05_FMRW_article/out_Au_Au200_20_05_FMR_pr_3D_proj_R_vs_Kt.root";
const TString Output_Folder = "/home/kirill/root-on-vs-code/Femto_output/20_05_FMRW_article/";

// Choose which fit type to plot: 0 = No Coulomb, 1 = Coulomb (Bowler-Sinyukov)
const Int_t fitType_to_plot = 1;

// ------------- Constants from original macro -----
const Int_t N_Charge  = 2; // 0 = Pi+Pi+, 1 = Pi-Pi-
const Int_t N_Centr   = 9; // 9 centrality bins
const Int_t N_Kt      = 4; // 4 K_t bins
const Int_t N_Pars    = 5; // N, lambda, R_o, R_s, R_l

// Radii names and corresponding parameter indices (0-based in par array)
const TString rad_names[3] = {"out", "side", "long"};
const Int_t   rad_idx[3]   = {2, 3, 4};   // R_out=par[2], R_side=par[3], R_long=par[4]

// K_T bin centres and half-widths
const Double_t Kt_center[4]    = {0.20, 0.30, 0.40, 0.525};
const Double_t Kt_halfwidth[4] = {0.05, 0.05, 0.05, 0.075};

// Centrality titles (for legend)
const TString cent_titles[9] = {
    "70-80%", "60-70%", "50-60%", "40-50%", "30-40%",
    "20-30%", "10-20%", "5-10%", "0-5%"
};

// Charge pair titles
const TString charge_titles[2] = {"#pi^{+}#pi^{+}", "#pi^{-}#pi^{-}"};

// Distinguishable colours for 9 centralities
const Int_t cent_colors[9] = {
    kRed, kBlue, kBlack, kMagenta, kGreen,
    kGray+2, kBlue, kRed, kBlack
};
const Int_t cent_markers[9] = {
    30, 31, 27, 29, 25,
    22, 24, 28, 20
};
// Celected centralities:
const Int_t selected_cents[] = {8,7,6,5,4};  // 0 <-> 70-80% ... 7 <-> 5-10%, 8 <-> 0-5%
const Int_t N_selected_cents = sizeof(selected_cents)/sizeof(selected_cents[0]);

// --------------- Data storage -------------------
Double_t R_val [N_Charge][N_Centr][N_Kt][3];  // [charge][cent][kt][rad]
Double_t R_err [N_Charge][N_Centr][N_Kt][3];
Bool_t   filled[N_Charge][N_Centr][N_Kt][3] = {false};

void CF_3D_osl_proj_R_paint() {
    // Open input file and get tree
    TFile *fin = TFile::Open(Input_File, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open input file " << Input_File << std::endl;
        return;
    }
    TTree *tree = (TTree*)fin->Get("fit_params_Tree");
    if (!tree) {
        std::cerr << "Tree 'fit_params_Tree' not found!" << std::endl;
        fin->Close();
        return;
    }

    // Tree variables
    Int_t   charge, centrality, Kt, fitType;
    Double_t par[N_Pars], parErr[N_Pars], chi2ndf;
    tree->SetBranchAddress("charge",     &charge);
    tree->SetBranchAddress("centrality", &centrality);
    tree->SetBranchAddress("Kt",         &Kt);
    tree->SetBranchAddress("fitType",    &fitType);
    tree->SetBranchAddress("par",        par);
    tree->SetBranchAddress("parErr",     parErr);
    tree->SetBranchAddress("chi2ndf",    &chi2ndf);

    // Reset arrays
    for (Int_t ic=0; ic<N_Charge; ic++)
        for (Int_t icent=0; icent<N_Centr; icent++)
            for (Int_t ikt=0; ikt<N_Kt; ikt++)
                for (Int_t ir=0; ir<3; ir++)
                    filled[ic][icent][ikt][ir] = false;

    // Read tree and fill arrays for chosen fitType
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i=0; i<nEntries; i++) {
        tree->GetEntry(i);
        if (fitType != fitType_to_plot) continue;
        if (charge<0 || charge>=N_Charge) continue;
        if (centrality<0 || centrality>=N_Centr) continue;
        if (Kt<0 || Kt>=N_Kt) continue;

        for (Int_t ir=0; ir<3; ir++) {
            Int_t idx = rad_idx[ir];
            R_val[charge][centrality][Kt][ir] = par[idx];
            R_err[charge][centrality][Kt][ir] = parErr[idx];
            filled[charge][centrality][Kt][ir] = true;
        }
    }

    fin->Close();

    // Create output file
    TFile *fout = new TFile(Output_File, "RECREATE");
    gStyle->SetOptStat(0);

    // Loop over charges (2) and radii (3) => 6 canvases
    for (Int_t ic=0; ic<N_Charge; ic++) {
        for (Int_t ir=0; ir<3; ir++) {

            // Build graph array for this (charge, radius) – one graph per centrality
            TGraphErrors* gr_cent[N_Centr] = {nullptr};
            bool has_any_data = false;

            for (Int_t icent=0; icent<N_Centr; icent++) {
                // Count valid points for this centrality
                std::vector<Double_t> x, y, ex, ey;
                for (Int_t ikt=0; ikt<N_Kt; ikt++) {
                    if (filled[ic][icent][ikt][ir]) {
                        x.push_back(Kt_center[ikt]);
                        y.push_back(R_val[ic][icent][ikt][ir]);
                        ex.push_back(Kt_halfwidth[ikt]);
                        ey.push_back(R_err[ic][icent][ikt][ir]);
                    }
                }
                if (x.empty()) continue;

                gr_cent[icent] = new TGraphErrors(x.size(), &x[0], &y[0], 0, &ey[0]);
                gr_cent[icent]->SetMarkerStyle(cent_markers[icent]);
                gr_cent[icent]->SetMarkerColor(cent_colors[icent]);
                gr_cent[icent]->SetLineColor(cent_colors[icent]);
                gr_cent[icent]->SetMarkerSize(1.9);
                has_any_data = true;
            }

            if (!has_any_data) {
                std::cout << "No data for charge " << ic << ", radius " << rad_names[ir] << std::endl;
                continue;
            }

            // Create canvas
            TString canvName = Form("c_%s_%s", charge_titles[ic].Data(), rad_names[ir].Data());
            // Make a safe canvas name (remove special chars)
            TString safeName = canvName;
            safeName.ReplaceAll("#","");
            safeName.ReplaceAll("{","");
            safeName.ReplaceAll("}","");
            safeName.ReplaceAll(" ","_");
            TCanvas *c = new TCanvas(safeName, canvName, 900, 700);
            c->SetGrid(0, 1);

            // Draw the first non-null graph to set axes
            TGraphErrors *firstGr = nullptr;
            for (Int_t idx_2=0; idx_2<N_selected_cents; idx_2++) {
            Int_t icent = selected_cents[idx_2];
                if (gr_cent[icent]) { firstGr = gr_cent[icent]; break; }
            }
            if (firstGr) {
                firstGr->SetTitle(Form("%s; K_{T} (GeV/c); R_{%s} (fm)",
                                       canvName.Data(), rad_names[ir].Data()));
                firstGr->Draw("AP");
                // Set reasonable Y range (maybe auto, but can be adjusted)
                firstGr->GetYaxis()->SetRangeUser(2.0, 7.0);
            }

            // Draw the rest
            for (Int_t idx_2=0; idx_2<N_selected_cents; idx_2++) {
            Int_t icent = selected_cents[idx_2];
                if (gr_cent[icent] && gr_cent[icent] != firstGr)
                    gr_cent[icent]->Draw("P SAME");
            }

            // Legend
            TLegend *leg = new TLegend(0.63, 0.45, 0.77, 0.89);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetBorderSize(1);
            leg->SetTextSize(0.05);
            for (Int_t idx_2=0; idx_2<N_selected_cents; idx_2++) {
            Int_t icent = selected_cents[idx_2];
                if (gr_cent[icent])
                    leg->AddEntry(gr_cent[icent], cent_titles[icent], "p");
            }
            leg->Draw();

            c->Write();
            // Optionally save as png/pdf
            // c->SaveAs(Output_Folder + safeName + ".png");
        }
    }

    fout->Close();
    std::cout << "All canvases saved to " << Output_File << std::endl;
}