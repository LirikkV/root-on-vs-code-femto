#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>
using namespace::std;


void CF_3D()
{
TFile *f = TFile::Open("/home/lirikk/root-on-vscode/Output_Data_1/Output_file_4.root", "READ");
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

// //Method 1 for out side long: first divide 3D, then project:
// TH3F *CF_Pi_Plus_osl_non_nrom = (TH3F*)hA_osl_Plus->Clone("CF_Pi_Plus_osl_non_norm");
// CF_Pi_Plus_osl_non_nrom->Sumw2();
// CF_Pi_Plus_osl_non_nrom->Divide(hB_osl_Plus);
// //out
// TH1F *CF_Pi_Plus_out_3D_first_non_norm = (TH1F*)CF_Pi_Plus_osl_non_nrom->Project3D("x");
// CF_Pi_Plus_out_3D_first_non_norm->Sumw2();
// CF_Pi_Plus_out_3D_first_non_norm->SetName("CF_Pi_Plus_out_3D_first_non_norm");
// CF_Pi_Plus_out_3D_first_non_norm->SetTitle("CF Pi+ Pi+ A/B - q_out 3D first non norm");
// CF_Pi_Plus_out_3D_first_non_norm->GetXaxis()->SetTitle("q_out, GeV/c");
// CF_Pi_Plus_out_3D_first_non_norm->GetYaxis()->SetTitle("CF");
// //side
// TH1F *CF_Pi_Plus_side_3D_first_non_norm = (TH1F*)CF_Pi_Plus_osl_non_nrom->Project3D("y");
// CF_Pi_Plus_side_3D_first_non_norm->Sumw2();
// CF_Pi_Plus_side_3D_first_non_norm->SetName("CF_Pi_Plus_side_3D_first_non_norm");
// CF_Pi_Plus_side_3D_first_non_norm->SetTitle("CF Pi+ Pi+ A/B - q_side 3D first non norm");
// CF_Pi_Plus_side_3D_first_non_norm->GetXaxis()->SetTitle("q_side, GeV/c");
// CF_Pi_Plus_side_3D_first_non_norm->GetYaxis()->SetTitle("CF");
// //long
// TH1F *CF_Pi_Plus_long_3D_first_non_norm = (TH1F*)CF_Pi_Plus_osl_non_nrom->Project3D("z");
// CF_Pi_Plus_long_3D_first_non_norm->Sumw2();
// CF_Pi_Plus_long_3D_first_non_norm->SetName("CF_Pi_Plus_long_3D_first_non_norm");
// CF_Pi_Plus_long_3D_first_non_norm->SetTitle("CF Pi+ Pi+ A/B - q_long 3D first non norm");
// CF_Pi_Plus_long_3D_first_non_norm->GetXaxis()->SetTitle("q_long, GeV/c");
// CF_Pi_Plus_long_3D_first_non_norm->GetYaxis()->SetTitle("CF");

//Method 2 for out side long: first project A and B then divide:
//out
TH1F *hA_Plus_out = (TH1F*)hA_osl_Plus->Project3D("x");
TH1F *hB_Plus_out = (TH1F*)hB_osl_Plus->Project3D("x");

hA_Plus_out->Sumw2();
hB_Plus_out->Sumw2();

TH1F *CF_Pi_Plus_out_Project_first_non_norm = (TH1F*)hA_Plus_out->Clone("CF_Pi_Plus_out_Project_first_non_norm");
CF_Pi_Plus_out_Project_first_non_norm->Sumw2();
CF_Pi_Plus_out_Project_first_non_norm->Divide(hB_Plus_out);
CF_Pi_Plus_out_Project_first_non_norm->SetTitle("CF Pi+ Pi+ A/B - q_out Project first non norm");
CF_Pi_Plus_out_Project_first_non_norm->GetXaxis()->SetTitle("q_out, GeV/c");
CF_Pi_Plus_out_Project_first_non_norm->GetYaxis()->SetTitle("CF");
//side
TH1F *hA_Plus_side = (TH1F*)hA_osl_Plus->Project3D("y");
TH1F *hB_Plus_side = (TH1F*)hB_osl_Plus->Project3D("y");

hA_Plus_side->Sumw2();
hB_Plus_side->Sumw2();

TH1F *CF_Pi_Plus_side_Project_first_non_norm = (TH1F*)hA_Plus_side->Clone("CF_Pi_Plus_side_Project_first_non_norm");
CF_Pi_Plus_side_Project_first_non_norm->Sumw2();
CF_Pi_Plus_side_Project_first_non_norm->Divide(hB_Plus_side);
CF_Pi_Plus_side_Project_first_non_norm->SetTitle("CF Pi+ Pi+ A/B - q_side Project first non norm");
CF_Pi_Plus_side_Project_first_non_norm->GetXaxis()->SetTitle("q_side, GeV/c");
CF_Pi_Plus_side_Project_first_non_norm->GetYaxis()->SetTitle("CF");
//long
TH1F *hA_Plus_long = (TH1F*)hA_osl_Plus->Project3D("z");
TH1F *hB_Plus_long = (TH1F*)hB_osl_Plus->Project3D("z");

hA_Plus_long->Sumw2();
hB_Plus_long->Sumw2();

TH1F *CF_Pi_Plus_long_Project_first_non_norm = (TH1F*)hA_Plus_long->Clone("CF_Pi_Plus_long_Project_first_non_norm");
CF_Pi_Plus_long_Project_first_non_norm->Sumw2();
CF_Pi_Plus_long_Project_first_non_norm->Divide(hB_Plus_long);
CF_Pi_Plus_long_Project_first_non_norm->SetTitle("CF Pi+ Pi+ A/B - q_long Project first non norm");
CF_Pi_Plus_long_Project_first_non_norm->GetXaxis()->SetTitle("q_long, GeV/c");
CF_Pi_Plus_long_Project_first_non_norm->GetYaxis()->SetTitle("CF");




TFile *f_out = new TFile("/home/lirikk/root-on-vscode/Output_Data_1/Output_file_4_pr.root", "RECREATE");
CF_Pi_Plus_non_nrom->Write();
CF_Pi_Plus->Write();
CF_Pi_Minus_Meth_1->Write();
//saving osl analysis hists:
//out
hA_Plus_out->Write();
hB_Plus_out->Write();
CF_Pi_Plus_out_Project_first_non_norm->Write();
//CF_Pi_Plus_out_3D_first_non_norm->Write();
//side
hA_Plus_side->Write();
hB_Plus_side->Write();
CF_Pi_Plus_side_Project_first_non_norm->Write();
//CF_Pi_Plus_side_3D_first_non_norm->Write();
//long
hA_Plus_long->Write();
hB_Plus_long->Write();
CF_Pi_Plus_long_Project_first_non_norm->Write();
//CF_Pi_Plus_long_3D_first_non_norm->Write();

//normalize CF:
double q_min_osl =0.15;
double q_max_osl = 0.4;

//out normalize
Double_t A_Pi_Plus_out_integral = hA_Plus_out->Integral(hA_Plus_out->FindBin(q_min_osl), hA_Plus_out->FindBin(q_max_osl));
Double_t B_Pi_Plus_out_integral = hB_Plus_out->Integral(hB_Plus_out->FindBin(q_min_osl),hB_Plus_out->FindBin(q_max_osl));
Double_t Scale_factor_Pi_Plus_out = A_Pi_Plus_out_integral/B_Pi_Plus_out_integral;

TH1D *hB_Plus_out_normalized = (TH1D*)hB_Plus_out->Clone("B_Pi_Plus_out");
hB_Plus_out_normalized->SetTitle("B of CF out Pi+ Pi+ normalized");

hB_Plus_out_normalized->Scale(Scale_factor_Pi_Plus_out);

//side normalize
Double_t A_Pi_Plus_side_integral = hA_Plus_side->Integral(hA_Plus_side->FindBin(q_min_osl), hA_Plus_side->FindBin(q_max_osl));
Double_t B_Pi_Plus_side_integral = hB_Plus_side->Integral(hB_Plus_side->FindBin(q_min_osl),hB_Plus_side->FindBin(q_max_osl));
Double_t Scale_factor_Pi_Plus_side = A_Pi_Plus_side_integral/B_Pi_Plus_side_integral;

TH1D *hB_Plus_side_normalized = (TH1D*)hB_Plus_side->Clone("B_Pi_Plus_side");
hB_Plus_side_normalized->SetTitle("B of CF side Pi+ Pi+ normalized");

hB_Plus_side_normalized->Scale(Scale_factor_Pi_Plus_side);
//long normalize
Double_t A_Pi_Plus_long_integral = hA_Plus_long->Integral(hA_Plus_long->FindBin(q_min_osl), hA_Plus_long->FindBin(q_max_osl));
Double_t B_Pi_Plus_long_integral = hB_Plus_long->Integral(hB_Plus_long->FindBin(q_min_osl),hB_Plus_long->FindBin(q_max_osl));
Double_t Scale_factor_Pi_Plus_long = A_Pi_Plus_long_integral/B_Pi_Plus_long_integral;

TH1D *hB_Plus_long_normalized = (TH1D*)hB_Plus_long->Clone("B_Pi_Plus_long");
hB_Plus_long_normalized->SetTitle("B of CF long Pi+ Pi+ normalized");

hB_Plus_long_normalized->Scale(Scale_factor_Pi_Plus_long);
//osl draw:
TCanvas *c_osl = new TCanvas("c_osl", "Canvas_osl",1920,1080);
c_osl->Divide(3,1);
gStyle->SetOptStat(0);

//out
c_osl->cd(1);
hA_Plus_out->SetLineColor(kGreen);
hA_Plus_out->SetMarkerStyle(21);
hA_Plus_out->SetMarkerSize(0.1);

hB_Plus_out_normalized->SetLineColor(kRed);
hB_Plus_out_normalized->SetMarkerStyle(20);
hB_Plus_out_normalized->SetMarkerSize(0.1);

hA_Plus_out->Draw();
hB_Plus_out_normalized->Draw("SAME");

//side
c_osl->cd(2);
hA_Plus_side->SetLineColor(kGreen);
hA_Plus_side->SetMarkerStyle(21);
hA_Plus_side->SetMarkerSize(0.1);

hB_Plus_side_normalized->SetLineColor(kRed);
hB_Plus_side_normalized->SetMarkerStyle(20);
hB_Plus_side_normalized->SetMarkerSize(0.1);

hA_Plus_side->Draw();
hB_Plus_side_normalized->Draw("SAME");

//long
c_osl->cd(3);
hA_Plus_long->SetLineColor(kGreen);
hA_Plus_long->SetMarkerStyle(21);
hA_Plus_long->SetMarkerSize(0.1);

hB_Plus_long_normalized->SetLineColor(kRed);
hB_Plus_long_normalized->SetMarkerStyle(20);
hB_Plus_long_normalized->SetMarkerSize(0.1);

hA_Plus_long->Draw();
hB_Plus_long_normalized->Draw("SAME");


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

c_osl->Write();
f_out->Close();
f->Close();
}