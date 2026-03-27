#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>
using namespace::std;


void corrFunct_devide()
{
TFile *f = TFile::Open("/home/lirikk/root-on-vscode/Output_Data_1/Output_file_1.root", "READ");
//for Pi+Pi+:
TH1D *hA_Plus = (TH1D*)f->Get("hA_Pi_Plus_q_inv_ALL");
TH1D *hB_Plus = (TH1D*)f->Get("hB_Pi_Plus_q_inv_ALL");

TH3F *hA_osl_Plus = (TH3F*)f->Get("hA_Pi_Plus_q_osl");
TH3F *hB_osl_Plus = (TH3F*)f->Get("hB_Pi_Plus_q_osl");

//for Pi-Pi-:
TH1D *hA_Minus = (TH1D*)f->Get("hA_Pi_Minus_q_inv_ALL");
TH1D *hB_Minus = (TH1D*)f->Get("hB_Pi_Minus_q_inv_ALL");

TH3F *hA_osl_Minus = (TH3F*)f->Get("hA_Pi_Minus_q_osl");
TH3F *hB_osl_Minus = (TH3F*)f->Get("hB_Pi_Minus_q_osl");


//!!! WARNING NEXT STRINGSs are NECECERY:
hA_Plus->Sumw2();
hB_Plus->Sumw2();
hA_Minus->Sumw2();
hB_Minus->Sumw2();
hA_osl_Plus->Sumw2();
hB_osl_Plus->Sumw2();
hA_osl_Minus->Sumw2();
hB_osl_Minus->Sumw2();


TH1D *CF_Pi_Plus_non_nrom = (TH1D*)hA_Plus->Clone("CF_Pi_Plus_non_normalized");
CF_Pi_Plus_non_nrom->Divide(hB_Plus);
CF_Pi_Plus_non_nrom->SetTitle("Corr.Funct Pi+ Pi+ A/B - non normalized");
CF_Pi_Plus_non_nrom->GetXaxis()->SetTitle("q_inv");
CF_Pi_Plus_non_nrom->GetYaxis()->SetTitle("CF");


//let's normalize CF:
Double_t q_min = 0.2;
Double_t q_max = 0.9;

Double_t A_Pi_Plus_integral = hA_Plus->Integral(hA_Plus->FindBin(q_min), hA_Plus->FindBin(q_max));
Double_t B_Pi_Plus_integral = hB_Plus->Integral(hB_Plus->FindBin(q_min),hB_Plus->FindBin(q_max));
Double_t Scale_factor_Pi_Plus = B_Pi_Plus_integral/A_Pi_Plus_integral;

Double_t A_Pi_Minus_integral = hA_Minus->Integral(hA_Minus->FindBin(q_min), hA_Minus->FindBin(q_max));
Double_t B_Pi_Minus_integral = hB_Minus->Integral(hB_Minus->FindBin(q_min),hB_Minus->FindBin(q_max));
Double_t Scale_factor_Pi_Minus = B_Pi_Minus_integral/A_Pi_Minus_integral;

//For hists scales:
Double_t x1_range = 0.0;
Double_t x2_range = 0.2;
Double_t y1_range = 0.4;
Double_t y2_range = 1.6;
/*
Method #1: 
1)Get relation of A/B in range of q_inv where no correlations 
2)Scale A by relation from 1)
3)CF=A_normilized/B
*/
TH1D *hA_Plus_normalized = (TH1D*)hA_Plus->Clone("CF_Pi_Plus");
hA_Plus_normalized->SetTitle("Numerator of Corr.Funct Pi+ Pi+ with both TPC & TPC+TOF methods normalized");


hA_Plus_normalized->Scale(Scale_factor_Pi_Plus);
TH1D *CF_Pi_Plus = (TH1D*)hA_Plus_normalized->Clone("CF_Pi_Plus");
CF_Pi_Plus->Divide(hB_Plus);

CF_Pi_Plus->SetTitle("Corr.Funct Pi+ Pi+ scale A then A/B");
CF_Pi_Plus->GetXaxis()->SetTitle("q_inv");
CF_Pi_Plus->GetYaxis()->SetTitle("CF");
CF_Pi_Plus->GetYaxis()->SetRangeUser(y1_range,y2_range);
CF_Pi_Plus->GetXaxis()->SetRangeUser(x1_range,x2_range);


TH1D *hA_Minus_normalized = (TH1D*)hA_Minus->Clone("CF_Pi_Minus");

hA_Minus_normalized->Scale(Scale_factor_Pi_Minus);
TH1D *CF_Pi_Minus_Meth_1 = (TH1D*)hA_Minus_normalized->Clone("CF_Pi_Minus_Meth_1");
CF_Pi_Minus_Meth_1->Divide(hB_Minus);

CF_Pi_Minus_Meth_1->SetTitle("Corr.Funct Pi- Pi- scale A then A/B");
CF_Pi_Minus_Meth_1->GetXaxis()->SetTitle("q_inv");
CF_Pi_Minus_Meth_1->GetYaxis()->SetTitle("CF");
CF_Pi_Minus_Meth_1->GetYaxis()->SetRangeUser(y1_range,y2_range);
CF_Pi_Minus_Meth_1->GetXaxis()->SetRangeUser(x1_range,x2_range);


//out-side-long analysis:
//hA_osl_Plus->GetYaxis()->SetRangeUser()

// TH3F *CF_Pi_Plus_osl_non_nrom = (TH3F*)hA_osl_Plus->Clone("CF_Pi_Plus_osl_non_normalized");
// CF_Pi_Plus_osl_non_nrom->Divide(hB_osl_Plus);
// TH1F *CF_Pi_Plus_side_non_norm = (TH1F*)CF_Pi_Plus_osl_non_nrom->Project3D("y");

TH1F *hA_Plus_out = (TH1F*)hA_osl_Plus->Project3D("x");
hA_Plus_out->SetName("hA_Plus_out");
TH1F *hB_Plus_out = (TH1F*)hB_osl_Plus->Project3D("x");

TH1F *hA_Plus_side = (TH1F*)hA_osl_Plus->Project3D("y");
TH1F *hB_Plus_side = (TH1F*)hB_osl_Plus->Project3D("y");

TH1F *hA_Plus_long = (TH1F*)hA_osl_Plus->Project3D("z");
TH1F *hB_Plus_long = (TH1F*)hB_osl_Plus->Project3D("z");


TH1F *CF_Pi_Plus_out_non_nrom = (TH1F*)hA_Plus_out->Clone("CF_Pi_Plus_out_non_normalized");
CF_Pi_Plus_out_non_nrom->Divide(hB_Plus_out);
CF_Pi_Plus_out_non_nrom->SetTitle("CF Pi+ Pi+ A/B - q_out non normalized");
CF_Pi_Plus_out_non_nrom->GetXaxis()->SetTitle("q_out");
CF_Pi_Plus_out_non_nrom->GetYaxis()->SetTitle("CF");

TH1F *CF_Pi_Plus_side_non_nrom = (TH1F*)hA_Plus_side->Clone("CF_Pi_Plus_side_non_normalized");
CF_Pi_Plus_side_non_nrom->Divide(hB_Plus_side);
CF_Pi_Plus_side_non_nrom->SetTitle("CF Pi+ Pi+ A/B - q_side non normalized");
CF_Pi_Plus_side_non_nrom->GetXaxis()->SetTitle("q_side");
CF_Pi_Plus_side_non_nrom->GetYaxis()->SetTitle("CF");

TH1F *CF_Pi_Plus_long_non_nrom = (TH1F*)hA_Plus_long->Clone("CF_Pi_Plus_long_non_normalized");
CF_Pi_Plus_long_non_nrom->Divide(hB_Plus_long);
CF_Pi_Plus_long_non_nrom->SetTitle("CF Pi+ Pi+ A/B - q_long non normalized");
CF_Pi_Plus_long_non_nrom->GetXaxis()->SetTitle("q_long");
CF_Pi_Plus_long_non_nrom->GetYaxis()->SetTitle("CF");



TFile *f_out = new TFile("/home/lirikk/root-on-vscode/Output_Data_1/Output_file_1_pr.root", "RECREATE");
CF_Pi_Plus_non_nrom->Write();
CF_Pi_Plus->Write();
CF_Pi_Minus_Meth_1->Write();

//CF_Pi_Plus_osl_non_nrom->Write();
hA_Plus_out->Write();
hB_Plus_out->Write("hB_Plus_out");
CF_Pi_Plus_out_non_nrom->Write();

hA_Plus_side->Write();
hB_Plus_side->Write();
CF_Pi_Plus_side_non_nrom->Write();
//CF_Pi_Plus_side_non_norm->Write();
hA_Plus_long->Write();
hB_Plus_long->Write();
CF_Pi_Plus_long_non_nrom->Write();

TCanvas *c1 = new TCanvas("c1", "Canvas",1920,1080);
hA_Plus->Draw();
c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hA_Plus.pdf");
hB_Plus->Draw();
c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hB_Plus.pdf");
CF_Pi_Plus_non_nrom->Draw();
c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/CF_Pi_Plus_non_norm.pdf");
CF_Pi_Plus->Draw();
c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hCF_Pi_Plus.pdf");
hA_Plus_normalized->Draw();

c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hA_Norm_Pi_Plus.pdf");
gStyle->SetOptStat(0);

hB_Plus->SetLineColor(kRed);
hA_Plus_normalized->SetLineColor(kGreen);

hB_Plus->GetXaxis()->SetRangeUser(0.0,0.1);
hA_Plus_normalized->GetXaxis()->SetRangeUser(0.0,0.1);

hB_Plus->SetMarkerStyle(21);
hA_Plus_normalized->SetMarkerStyle(20);

hB_Plus->SetMarkerSize(0.5);
hA_Plus_normalized->SetMarkerSize(0.5);

hB_Plus->SetTitle("Normalized numerator(A) and denumerator (B)");
hB_Plus->Draw();
hA_Plus_normalized->Draw("SAME");

TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
leg->AddEntry(hB_Plus, "B", "lp");
leg->AddEntry(hA_Plus_normalized, "A normalized" , "lp");
leg->Draw();

TString Output_data_folder = "/home/lirikk/root-on-vscode/Output_Data_1/";
c1->SaveAs(Output_data_folder + "hA_B_difference.pdf");
//Now series of saved pdf-s:

// TH1D* hRefMult = (TH1D*)f->Get("hRefMult"); 
// TH2D* hVtxXvsY = (TH2D*)f->Get("hVtxXvsY"); 
// TH2D* hVtxXvsY_cut = (TH2D*)f->Get("hVtxXvsY_cut"); 
// TH1D* hVtxZ = (TH1D*)f->Get("hVtxZ"); 
// TH1D* hVtxZ_cut = (TH1D*)f->Get("hVtxZ_cut"); 
// TH1D* hPrimaryPtot = (TH1D*)f->Get("hPrimaryPtot"); 
// TH1D* hPrimaryPtot_cut = (TH1D*)f->Get("hPrimaryPtot_cut"); 
// TH1D* hNFitHits = (TH1D*)f->Get("hNFitHits"); 
// TH1D* hNFitHits_cut = (TH1D*)f->Get("hNFitHits_cut"); 
// TH1D* hDCA = (TH1D*)f->Get("hDCA"); 
// TH1D* hDCA_cut = (TH1D*)f->Get("hDCA_cut"); 
// TH1D* hPrimaryPtrans = (TH1D*)f->Get("hPrimaryPtrans"); 
// TH1D* hPrimaryPttrans_cut = (TH1D*)f->Get("hPrimaryPttrans_cut"); 
// TH1D* hPrimaryPseudorap = (TH1D*)f->Get("hPrimaryPseudorap"); 
// TH1D* hPrimaryPseudorap_cut = (TH1D*)f->Get("hPrimaryPseudorap_cut"); 
// TH2D* hNSigmPion_vs_pPrimTotDevQ = (TH2D*)f->Get("hNSigmPion_vs_pPrimTotDevQ"); 
// TH2D* hNSigmPion_vs_pPrimTotDevQ_cut_PID = (TH2D*)f->Get("hNSigmPion_vs_pPrimTotDevQ_cut_PID"); 
// TH2D* hNSigmKaon_vs_pPrimTotDevQ = (TH2D*)f->Get("hNSigmKaon_vs_pPrimTotDevQ"); 
// TH2D* hNSigmKaon_vs_pPrimTotDevQ_cut_PID = (TH2D*)f->Get("hNSigmKaon_vs_pPrimTotDevQ_cut_PID"); 
// TH2D* hNSigmProton_vs_pPrimTotDevQ = (TH2D*)f->Get("hNSigmProton_vs_pPrimTotDevQ"); 
// TH2D* hNSigmProton_vs_pPrimTotDevQ_cut_PID = (TH2D*)f->Get("hNSigmProton_vs_pPrimTotDevQ_cut_PID"); 
// TH2D* hNSigmElectron_vs_pPrimTotDevQ = (TH2D*)f->Get("hNSigmElectron_vs_pPrimTotDevQ"); 
// TH2D* hNSigmElectron_vs_pPrimTotDevQ_cut_PID = (TH2D*)f->Get("hNSigmElectron_vs_pPrimTotDevQ_cut_PID"); 
// TH2D* h1_OverBeta_vs_pPrimTotDevQ = (TH2D*)f->Get("1_OverBeta_vs_pPrimTotDevQ"); 
// TH2D* h1_OverBeta_vs_pPrimTotDevQ_cut_PID = (TH2D*)f->Get("1_OverBeta_vs_pPrimTotDevQ_cut_PID"); 
// TH2D* hm2_vs_pPrimTotDevQ = (TH2D*)f->Get("hm2_vs_pPrimTotDevQ"); 
// TH2D* hm2_vs_pPrimTotDevQ_cut_PID = (TH2D*)f->Get("hm2_vs_pPrimTotDevQ_cut_PID"); 
// TH1D* hTEST_P2_pPrimTotDevQ = (TH1D*)f->Get("hTEST_P2_pPrimTotDevQ"); 
// TH2D* hTEST_sqrt_Beta_pPrimTotDevQ = (TH2D*)f->Get("hTEST_sqrt_Beta_pPrimTotDevQ"); 
// TH1D* h1_OverBetaDelta_vs_pPrimTotDevQ = (TH1D*)f->Get("1_OverBetaDelta_vs_pPrimTotDevQ"); 
// TH1D* h1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID = (TH1D*)f->Get("1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID"); 
// TH2D* hdEdx_vs_pPrimTotDevQ = (TH2D*)f->Get("hdEdx_vs_pPrimTotDevQ"); 
// TH2D* hdEdx_vs_pPrimTotDevQ_cut_PID = (TH2D*)f->Get("hdEdx_vs_pPrimTotDevQ_cut_PID"); 
// TH2D* h2DpPrimTr_vs_etaPrim = (TH2D*)f->Get("h2DpPrimTr_vs_etaPtim"); 
// TH2D* h2DpPrimTr_vs_etaPrim_cut = (TH2D*)f->Get("h2DpPrimTr_vs_etaPtim_cut"); 
// TH2D* hTEST2DpPrimTr_vs_etaPrim_equal_P = (TH2D*)f->Get("hTEST2DpPrimTr_vs_etaPtim_equal_P"); 
// TH1D* hTofBeta = (TH1D*)f->Get("hTofBeta"); 
// TH1D* hA_Pi_Plus_q_inv_ALL = (TH1D*)f->Get("hA_Pi_Plus_q_inv_ALL"); 
// TH1D* hB_Pi_Plus_q_inv_ALL = (TH1D*)f->Get("hB_Pi_Plus_q_inv_ALL"); 
// TH1D* hA_Pi_Minus_q_inv_ALL = (TH1D*)f->Get("hA_Pi_Minus_q_inv_ALL"); 
// TH1D* hB_Pi_Minus_q_inv_ALL = (TH1D*)f->Get("hB_Pi_Minus_q_inv_ALL"); 

// c1->Clear();
// TString data_of_report = "02.03.2026_Report/";
// hRefMult->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hRefMult.png");
// hVtxXvsY->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hVtxXvsY.png");
// hVtxXvsY_cut->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hVtxXvsY_cut.png");
// hVtxZ->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hVtxZ.png");
// hVtxZ_cut->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hVtxZ_cut.png");
// hPrimaryPtot->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hPrimaryPtot.png");
// hPrimaryPtot_cut->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hPrimaryPtot_cut.png");
// hNFitHits->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hNFitHits.png");
// hNFitHits_cut->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hNFitHits_cut.png");
// hDCA->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hDCA.png");
// hDCA_cut->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hDCA_cut.png");
// hPrimaryPtrans->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hPrimaryPtrans.png");
// hPrimaryPttrans_cut->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hPrimaryPttrans_cut.png");
// hPrimaryPseudorap->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hPrimaryPseudorap.png");
// hPrimaryPseudorap_cut->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hPrimaryPseudorap_cut.png");
// hNSigmPion_vs_pPrimTotDevQ->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hNSigmPion_vs_pPrimTotDevQ.png");
// hNSigmPion_vs_pPrimTotDevQ_cut_PID->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hNSigmPion_vs_pPrimTotDevQ_cut_PID.png");
// hNSigmKaon_vs_pPrimTotDevQ->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hNSigmKaon_vs_pPrimTotDevQ.png");
// hNSigmKaon_vs_pPrimTotDevQ_cut_PID->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hNSigmKaon_vs_pPrimTotDevQ_cut_PID.png");
// hNSigmProton_vs_pPrimTotDevQ->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hNSigmProton_vs_pPrimTotDevQ.png");
// hNSigmProton_vs_pPrimTotDevQ_cut_PID->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hNSigmProton_vs_pPrimTotDevQ_cut_PID.png");
// hNSigmElectron_vs_pPrimTotDevQ->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hNSigmElectron_vs_pPrimTotDevQ.png");
// hNSigmElectron_vs_pPrimTotDevQ_cut_PID->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hNSigmElectron_vs_pPrimTotDevQ_cut_PID.png");
// h1_OverBeta_vs_pPrimTotDevQ->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "h1_OverBeta_vs_pPrimTotDevQ.png");
// h1_OverBeta_vs_pPrimTotDevQ_cut_PID->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "h1_OverBeta_vs_pPrimTotDevQ_cut_PID.png");
// hm2_vs_pPrimTotDevQ->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hm2_vs_pPrimTotDevQ.png");
// hm2_vs_pPrimTotDevQ_cut_PID->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hm2_vs_pPrimTotDevQ_cut_PID.png");
// h1_OverBetaDelta_vs_pPrimTotDevQ->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "h1_OverBetaDelta_vs_pPrimTotDevQ.png");
// h1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "h1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID.png");
// hdEdx_vs_pPrimTotDevQ->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hdEdx_vs_pPrimTotDevQ.png");
// hdEdx_vs_pPrimTotDevQ_cut_PID->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hdEdx_vs_pPrimTotDevQ_cut_PID.png");
// h2DpPrimTr_vs_etaPrim->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "h2DpPrimTr_vs_etaPrim.png");
// h2DpPrimTr_vs_etaPrim_cut->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "h2DpPrimTr_vs_etaPrim_cut.png");
// hTEST2DpPrimTr_vs_etaPrim_equal_P->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hTEST2DpPrimTr_vs_etaPrim_equal_P.png");
// hTofBeta->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hTofBeta.png");
// hA_Pi_Plus_q_inv_ALL->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hA_Pi_Plus_q_inv_ALL.png");
// hB_Pi_Plus_q_inv_ALL->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hB_Pi_Plus_q_inv_ALL.png");
// hA_Pi_Minus_q_inv_ALL->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hA_Pi_Minus_q_inv_ALL.png");
// hB_Pi_Minus_q_inv_ALL->Draw("COLZ");
// c1->SaveAs(Output_data_folder + data_of_report + "hB_Pi_Minus_q_inv_ALL.png");



f_out->Close();
f->Close();
}