
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <TStyle.h>
#include <TLine.h>

#include "ACCSEL.h"

void EffSig() {

// MC FILEs
const char * files[] = {
   // "/lstore/cms/hlegoinha/Bmeson/MC_DATA/MC_ppRef_Bmeson/Bu_phat5_Bfinder.root",     //ppRef                         
   // "/lstore/cms/hlegoinha/Bmeson/MC_DATA/MC_ppRef_Bmeson/Bd_phat5_Bfinder.root"         //ppRef old
    //"/lstore/cms/lekai/Bmeson/MC/ppRef/Bd_phat5_Bfinder.root"  //ppRef New 
   //"/lstore/cms/lekai/Bmeson/MC/ppRef/Bs_phat5_Bfinder.root"     //ppRef New
     "/lstore/cms/lekai/Bmeson/MC/ppRef/Bu_phat5_Bfinder.root" //ppRef New 
   // "/lstore/cms/hlegoinha/Bmeson/MC_DATA/MC_ppRef_Bmeson/Bs_phat5_Bfinder.root"
   // "/lstore/cms/hlegoinha/X3872/MC_DATA/prompt_PSI2S_to_Jpsi_pipi_phat5_Bfinder.root" //dados MC PSI2S
   // "/lstore/cms/hlegoinha/X3872/MC_DATA/prompt_X3872_to_Jpsi_Rho_phat5_Bfinder.root" //dados MC X
           
};

const char* files_data={
     "/lstore/cms/hlegoinha/Bmeson/MC_DATA/DATA_ppref_Bmeson/DATA_ppref_Bmeson.root" // Dados completos
    //"/user/u/u25pedrochan/HiForestMINIAOD_ppRefData_100.root" // Real data
    // "/lstore/cms/hlegoinha/X3872/MC_DATA/DATA_ppRef_X3872.root"
};

//VARIABLES
//VARIABLES
const char * variables[] = {/*"Balpha",  "BQvalueuj"    ,  "Bcos_dtheta", "BtrkPtimb","Bchi2cl", "Btrk1dR", "Btrk2dR", "Btrk1Pt",    "Btrk2Pt",*/ "Bnorm_svpvDistance_2D", "Bnorm_svpvDistance"/* , "Bnorm_trk1Dxy"  , "Bnorm_trk2Dxy"*/};
const double ranges[][2] = {/*{0,3.15},    {0,2.5}    ,     {0.99,1}  ,    {0,1},       {0,1},  {0,4.5},    {0,4.5},  {0, 10} , {0, 10},  */ {0,85},                       {0,85}   ,  /*       {-22,22}  ,          {-22,22}*/};
int SELplots = 0; //mudar para 1 com ruído e descomentar a linha acima

//const char * variables[] = {"Bmass"/*,  "Btktkmass" , "Bpt", "By", "nSelectedChargedTracks"*/};//comentar o Btkmass
//const double ranges[][2] = { {5 , 6}/*,{0,2.5} ,{5, 50}, {-2.4, 2.4}, {0,150}*/};
//VARIABLES
//VARIABLES

/////////////////////////////////  ///////////////////////////  ////////////////

TString cutlevel = ""; // "_RAW", "_ACC", "_SEL", "_TRG", "", 

/////////////////////////////////  ///////////////////////////  ///////////

TString path_to_file = "";

const int nVars = sizeof(variables)/sizeof(variables[0]);

for (int ifile = 0; ifile < sizeof(files)/sizeof(files[0]); ++ifile) {
    path_to_file = Form("%s", files[ifile]);
    //path_to_file = Form("/eos/user/h/hmarques/MC_ppRef_Bmeson/MC_ppRef_Bmeson/%s_Bfinder.root", files[ifile]);

     

    TFile *file = TFile::Open(path_to_file.Data());
    TFile *file_data = TFile::Open(files_data);
    // Get the trees from the file
    TTree *treeMix;
    TTree *treedata;
    if (path_to_file.Contains("Bs")){                             //Bs
        file->GetObject("Bfinder/ntphi", treeMix);
        file_data->GetObject("Bfinder/ntphi", treedata);
    }else if (path_to_file.Contains("Bd")){                      //Bd
        file->GetObject("Bfinder/ntKstar", treeMix);
        file_data->GetObject("Bfinder/ntKstar", treedata);
    }else if(path_to_file.Contains("Bu")){                       //Bu
        file->GetObject("Bfinder/ntKp", treeMix);
        file_data->GetObject("Bfinder/ntKp", treedata);
    }else{                                                        //X3872
         file->GetObject("Bfinder/ntmix", treeMix);//PSI2S  
        // filex->GetObject("Bfinder/ntmix", treex);//X3872                                       
         file_data->GetObject("Bfinder/ntmix", treedata);
    }

    std::cout << "\n" << "Entries in treeMix: " << treeMix->GetEntries() << std::endl;
    std::cout << "\n" << "Entries in treedata: " << treedata->GetEntries() << std::endl;

    for (int i = 0; i < nVars; ++i) {
        TString var = variables[i];

        if(path_to_file.Contains("Bu") && ((var.Contains("trk2") || var.Contains("Ptimb")))) continue; // B+ has less one track!

        // Create a canvas to draw the histograms
        /*TCanvas *canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetLeftMargin(0.15);
        canvas->SetTopMargin(0.05);
        canvas->SetRightMargin(0.05); */

        double hist_Xhigh      = ranges[i][1];
        double hist_Xlow       = ranges[i][0];
        int hist_Nbin          = 50000 ;
        if (var == "nSelectedChargedTracks") {
            hist_Nbin = (hist_Xhigh - hist_Xlow)/0.01;
        } 
        double bin_length_MEV  = (hist_Xhigh - hist_Xlow) / hist_Nbin;
        if(SELplots){ hist_Nbin = 50; }
        
        TString Xlabel ;
        if (var == "Bmass"){ 
            if (path_to_file.Contains("Bs")){
                Xlabel = "m_{J/#Psi K^{+} K^{-}} [GeV/c^{2}]";
            } else if (path_to_file.Contains("Bd")){
                Xlabel = "m_{J/#Psi K^{+} #pi^{-}} [GeV/c^{2}]";
            } else if(path_to_file.Contains("PSI2S")){
                Xlabel="m_{J/#Psi#pi^{+}#pi^{-}} [GeV/c^{2}]";
            } else if(path_to_file.Contains("Rho")){
                Xlabel = "m_{J/#Psi#rho} [GeV/c^{2}]";
            } else {
                Xlabel = "m_{J/#Psi K^{+}} [GeV/c^{2}]";
            }
        } else if (var == "Bpt"){ 
            Xlabel = "p_{T} [GeV/c]";
        } else { 
            Xlabel = var.Data();
        }

        // Create histograms
        TH1F *hist_SIG = new TH1F("hist_SIG"      , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh); 
        //TH1F *hist_sig = new TH1F("hist_sig"      , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh);
        TH1F *hist_BKG = new TH1F("hist_BKG"      , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh);
        TH1F *hist     = new TH1F("hist"          , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh);
        TH1F *hist_SIG_WT   = new TH1F("hist_SIG_WT"  , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh);        
        TH1F *hist_SIG_BOTH = new TH1F("hist_SIG_BOTH", Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , hist_Nbin, hist_Xlow ,hist_Xhigh);        

        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        //SELECT THE acc + presel CUT 

        TString dirNAME = "";
        TString Final = "1";      
        TString trgmatches = TRGmatching.Data();   //TRG matching only in ppRef
        TString ACCcuts = "" ;
        TString SELcuts = "" ;

        if (path_to_file.Contains("Bu")){
            ACCcuts    = ACCcuts_ppRef_Bu.Data(); //ppRef
            SELcuts    = SELcuts_ppRef_Bu.Data(); //ppRef
            if (path_to_file.Contains("PbPb")) { 
                ACCcuts = ACCcuts_PbPb_Bu.Data();
                SELcuts = SELcuts_PbPb_Bu.Data();
                trgmatches = "1";
            }
        }
        else {
            ACCcuts    = ACCcuts_ppRef.Data(); //ppRef
            SELcuts    = SELcuts_ppRef.Data(); //ppRef
            if (path_to_file.Contains("PbPb")) { 
                ACCcuts = ACCcuts_PbPb.Data();
                SELcuts = SELcuts_PbPb.Data();
                trgmatches = "1";
            }
        }

        TString cut = "";
        if (cutlevel == "_RAW")       {cut = Form(" %s "                   ,FIDreg.Data());}                                                              //RAW (inside fid reg only)
        else if (cutlevel == "_ACC")  {cut = Form(" %s && %s "             ,FIDreg.Data(), ACCcuts.Data());}                                              //ACC
        else if (cutlevel == "_SEL")  {cut = Form(" %s && %s && %s "       ,FIDreg.Data(), ACCcuts.Data(), SELcuts.Data());}                              //SEL
        else if (cutlevel == "_TRG")  {cut = Form(" %s && %s && %s && %s " ,FIDreg.Data(), ACCcuts.Data(), SELcuts.Data(), trgmatches.Data());}           //TRG
        else if (cutlevel == ""){
            if (!SELplots) {dirNAME  = "";}
            cut = Form(" %s && %s && %s && %s", ACCcuts.Data(), SELcuts.Data(), trgmatches.Data(), Final.Data());                   //Final
        }
        else{
            std::cerr << "Invalid cut level specified: " << cutlevel << std::endl;
            return;
        }                                                                                                 

        TString sepcCASES = "1";
        if (path_to_file.Contains("Bs")){ 
            sepcCASES = "abs(Btktkmass - 1.019455) < 0.015"; // phi meson mass cut
             treeMix->Draw(Form("%s >> hist_SIG", var.Data()), Form( "Bnorm_svpvDistance_2D>3.2232 && Bchi2cl>0.003 && Bnorm_svpvDistance>2 && %s && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));
             treedata->Draw(Form("%s >> hist_BKG", var.Data()), Form("Bnorm_svpvDistance_2D>3.2232 && Bchi2cl>0.003 && ((Bmass <5.3148002594534231054413939881266) || (Bmass > 5.4200597405465768945586060118734)) && %s && %s", cut.Data(), sepcCASES.Data()));
        } else if (path_to_file.Contains("Bd")){ 
            sepcCASES = "abs(Btktkmass - 0.89594) < 0.25"; // Kstar meson mass cut
        }
      
          
        //treeMix->Draw(Form("%s >> hist_SIG", var.Data()), Form("Balpha<0.173754 && Bnorm_svpvDistance_2D>4.473 && Bnorm_svpvDistance>2 && %s && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));  // SIG
        if (path_to_file.Contains("Bd")){ 
          //treeMix->Draw(Form("%s >> hist_SIG_WT"  , var.Data()), Form(" (Bgen == 41000) && %s && %s", cut.Data(), sepcCASES.Data()));                              // WT component
	   treedata->Draw(Form("%s >> hist_BKG" , var.Data()), Form("Bchi2cl>0.003 && ((Bmass < 5.1) || (Bmass > 5.6)) && %s && %s", cut.Data(), sepcCASES.Data())); // BKG -- (notice the *!* in the first %s)
           treeMix->Draw(Form("%s >> hist_SIG_BOTH", var.Data()), Form("Bchi2cl>0.003 && Bnorm_svpvDistance>2 && (%s || (Bgen == 41000)) && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));  // SIG + WT
        } else if(path_to_file.Contains("Bu")) {
            treeMix->Draw(Form("%s >> hist_SIG", var.Data()), Form( "Bchi2cl>0.003 && Bnorm_svpvDistance>2 && %s && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));
            treedata->Draw(Form("%s >> hist_BKG", var.Data()), Form("Bchi2cl>0.003  && (Bmass>5.380091232) && %s && %s", cut.Data(), sepcCASES.Data())); // BKG -- (notice the *!* in the first %s)4        
        } else if(path_to_file.Contains("Rho")) {
            treeMix->Draw(Form("%s >> hist_SIG", var.Data()), Form("%s && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));
            //treex->Draw(Form("%s >> hist_sig", var.Data()), Form("%s && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));
            //treedata->Draw(Form("%s >> hist_BKG", var.Data()), Form("((Bmass < 3.65826) || (Bmass > 3.7156 && Bmass < 3.83176) || (Bmass > 3.91336)) && %s && %s", cut.Data(),sepcCASES.Data()));
              treedata->Draw(Form("%s >> hist_BKG", var.Data()), Form("Bmass > 3.77369 && (( Bmass < 3.83176) || (Bmass > 3.91336)) && %s", cut.Data()));
        } else if(path_to_file.Contains("PSI2S")){  
            treeMix->Draw(Form("%s >> hist_SIG", var.Data()), Form("%s && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));
            treedata->Draw(Form("%s >> hist_BKG", var.Data()), Form(" Bmass < 3.77369 && (( Bmass < 3.65826) || (Bmass > 3.71562))  && %s", cut.Data()));

              }
        //treeMix->Draw(Form("%s >> hist"    , var.Data()), Form(" %s && %s", cut.Data(), sepcCASES.Data()) );                          // ALL 

        //SELECT THE acc + presel CUT 
        ///////////////////////////////////////////////////////////////////////////////////////////////////////// 

    
 
        // Set display titles (for legend and stat box)
        hist_SIG->SetName("MC_SIG");  // <--- This affects the stat box label

        hist_BKG->SetName("DATA_BKG");  // <--- Also affects the stat box

        hist_SIG_BOTH->SetName("MC_SIG");
        if(path_to_file.Contains("Rho")){
            hist_SIG->SetName("MC_SIG_X3872");
            hist_SIG->SetName("MC_SIG_PSI2S");}

         // <--- This affects the stat box label

         // <--- Also affects the stat box

        


       if (SELplots == 1) { // NORMALIZE
              double int_sig     = hist_SIG->Integral();
              double int_bkg     = hist_BKG->Integral();
              double int_sig_wt  = hist_SIG_BOTH->Integral();

         if (int_sig > 0 || int_sig_wt > 0){     
                hist_SIG->Scale(1.0 / int_sig);
                hist_BKG->Scale(1.0 / int_bkg);
                hist_SIG_BOTH->Scale(1.0 / int_sig_wt); 
          }
        }

        if(1){// set the y-axis maximum if needed
            Double_t     max_val = TMath::Max(hist->GetMaximum(), TMath::Max(hist_BKG->GetMaximum(), hist_SIG->GetMaximum()));
            if(SELplots) {
                if (path_to_file.Contains("Bd")) {
                    max_val = TMath::Max( hist_BKG->GetMaximum(), hist_SIG_BOTH->GetMaximum()) ;
                    hist_SIG_BOTH->SetMaximum(max_val * 1.1);  // Increase the max range to give some space
                } else {
                    max_val = TMath::Max( hist_SIG->GetMaximum(), hist_BKG->GetMaximum());
                    hist_SIG->SetMaximum(max_val * 1.1);  // Increase the max range to give some space
                }
                hist_BKG->SetMaximum(max_val * 1.1); 
            } else {
                hist_SIG->SetMaximum(max_val * 1.1);  // Increase the max range to give some space
                hist_BKG->SetMaximum(max_val * 1.1);
            }
        }

        // Draw the histograms
        hist->SetStats(0);
        if (SELplots && path_to_file.Contains("Bd")){
            hist_SIG_BOTH->Draw("HIST");
        } else if(SELplots && path_to_file.Contains("Rho")){
            hist_SIG->Draw("HIST"); 
            hist_BKG->Draw("HIST SAMES");
        }else if(SELplots && path_to_file.Contains("PSI2S")){
            hist_SIG->Draw("HIST"); 
            hist_BKG->Draw("HIST SAMES");
        }else {
            hist_SIG->Draw("HIST");
            if (path_to_file.Contains("Bd")) {
                hist_SIG_WT->Draw("HIST SAMES");
            }
        }
        hist_BKG->Draw("HIST SAMES");

        if(!SELplots) hist->Draw("HIST SAME");

        


        /*double fs;
        double fb ;
        if (path_to_file.Contains("Bu")){
            fs = 9.022;
            fb = 0.349;  
        }else if (path_to_file.Contains("Bd")){
            fs=3.769;
            fb=0.376;                
        } else if (path_to_file.Contains("Bs")){
            fs=;
            fb=;
        } else if (path_to_file.Contains("Rho")){
        
        } else if (path_to_file.Contains("PSI2S")){

        }*/
           
        if (0) { // NORMALIZE
              double int_sig     = hist_SIG->Integral();
              double int_bkg     = hist_BKG->Integral();
              double int_sig_wt  = hist_SIG_BOTH->Integral();

          if (int_sig > 0 || int_sig_wt > 0){     
                hist_SIG->Scale(1.0 / int_sig);
                hist_BKG->Scale(1.0 / int_bkg);
                hist_SIG_BOTH->Scale(1.0 / int_sig_wt); 
          }
        }


        TString particleNAME = "Bu";
        TString systemNAME = "ROC_ppRef_";
        if (path_to_file.Contains("Bs")){
            particleNAME = "Bs";
        } else if (path_to_file.Contains("Bd")){
            particleNAME = "Bd";
        } else if (path_to_file.Contains("Rho")){
            particleNAME = "X3872";
        } else if (path_to_file.Contains("PSI2S")){
            particleNAME = "PSI2S";
        }
        if (path_to_file.Contains("PbPb"))  { systemNAME = "ROC_PbPb_";}

        std::vector<double> effBkg_right, effSig_right;
        std::vector<double> effBkg_left, effSig_left;

        /*double maxSig1 = -1;
        double bestCut1;
        double B_Cut1;
        double S_Cut1;
        TString bestDirection1 = "";

        double maxSig2 = -1;
        double bestCut2;
        double B_Cut2;
        double S_Cut2;
        TString bestDirection2="";

        double maxSig_both;
        double bestCut_both;
        TString bestDirection="";*/

        TH1F *histogram = nullptr;
        int nbins;
        // Save the canvas as an image
        if(path_to_file.Contains("Bd")){
            nbins = hist_SIG_BOTH->GetNbinsX();
            std::cout << " nbins = " << nbins << std::endl;
            histogram = hist_SIG_BOTH;
        }else{
            nbins = hist_SIG->GetNbinsX();
            std::cout << " nbins = " << nbins << std::endl;
            histogram = hist_SIG;
        }

        double S_total = histogram->Integral(1, nbins);  
        double B_total = hist_BKG->Integral(1, nbins);                                

        // === Left-side cut: cut < x  ===
        for (int ibin = 1; ibin <= nbins; ++ibin) {                                           //std::cout << "  B = " << B << std::endl;
            double S_pass = histogram->Integral(ibin, nbins);                                  //std::cout << "  Significance1 = " << significance << std::endl;
            double B_pass = hist_BKG->Integral(ibin, nbins);                                       //std::cout << "  Cut Direction1    = " << "cut < x" << std::endl;
            double effSig =  S_pass/S_total;        
            double effBkg = (B_pass/B_total);

            effBkg_right.push_back(effBkg);
            effSig_right.push_back(effSig);
        }   

double minDistRight = 999.0;
int bestBinRight = -1;

for (int i = 0; i < nbins; ++i) {
    double fpr = effBkg_right[i];
    double tpr = effSig_right[i];
    double dist = std::sqrt(std::pow(fpr, 2) + std::pow(1 - tpr, 2));

    if (dist < minDistRight) {
        minDistRight = dist;
        bestBinRight = i;
    }
}
double cutValueRight = histogram->GetBinLowEdge(bestBinRight + 1); // +1 because bins start at 1 in ROOT
if (bestBinRight >= 0) {
    double cutValueRight = histogram->GetBinLowEdge(bestBinRight + 1); // +1 because bins start at 1 in ROOT
    std::cout << "\n[RIGHT CUT] Best ROC Point (Closest to (0,1)):" << std::endl;
    std::cout << "  Cut Direction: x > " << cutValueRight << std::endl;
    std::cout << "  FPR: " << effBkg_right[bestBinRight] << ", TPR: " << effSig_right[bestBinRight] << std::endl;
    std::cout << "  Distance to (0,1): " << minDistRight << std::endl;
}


             TCanvas* cRight = new TCanvas("c", "", 800, 600);
             cRight->SetGrid();
             TGraph* rocRight = new TGraph(effSig_right.size(), &effBkg_right[0], &effSig_right[0]);

            rocRight->SetTitle(";False Postive Rate;True Positive Rate");
            rocRight->SetLineColor(kBlue);
            rocRight->SetLineWidth(2);
            rocRight->GetXaxis()->SetLimits(0.0, 1.0);
            rocRight->GetYaxis()->SetRangeUser(0.0, 1.0);

            rocRight->Draw("AL");

// --- Marker + label at chosen point ---
if (bestBinRight >= 0) {
    double xR = effBkg_right[bestBinRight];
    double yR = effSig_right[bestBinRight];

    TMarker* bestPtRight = new TMarker(xR, yR, 20);
    bestPtRight->SetMarkerColor(kBlue+2);
    bestPtRight->SetMarkerSize(1.4);
    bestPtRight->Draw("SAME");

    TLatex labR;
    labR.SetTextSize(0.028);
    labR.SetTextColor(kBlue+2);
    double xr = std::min(0.95, std::max(0.05, xR + 0.02));
    double yr = std::min(0.95, std::max(0.05, yR - 0.04));
    labR.DrawLatex(xr, yr, Form("cut = %.5g", cutValueRight));
}

std::vector<std::pair<double, double>> rocPoints_right;
for (size_t i = 0; i < effBkg_right.size(); ++i)
    rocPoints_right.emplace_back(effBkg_right[i], effSig_right[i]);

std::sort(rocPoints_right.begin(), rocPoints_right.end());  // sort by FPR (x)

// Reassign sorted values
for (size_t i = 0; i < rocPoints_right.size(); ++i) {
    effBkg_right[i] = rocPoints_right[i].first;
    effSig_right[i] = rocPoints_right[i].second;
}


// --- Compute AUC (Left-side cut) ---
double aucRight = 0.0;
for (size_t i = 1; i < effBkg_right.size(); ++i) {
    double x1 = effBkg_right[i - 1], x2 = effBkg_right[i];
    double y1 = effSig_right[i - 1], y2 = effSig_right[i];
    aucRight += 0.5 * (x2 - x1) * (y1 + y2); // Trapezoid rule
}

// Draw AUC on canvas
/*TLatex latexRight;
latexRight.SetNDC();
latexRight.SetTextSize(0.03);
latexRight.DrawLatex(0.2, 0.88, Form("AUC = %.4f", aucRight));*/

auto legend_r = new TLegend(0.6, 0.75, 0.88, 0.88);
    legend_r->SetHeader(particleNAME, "");
    legend_r->AddEntry((TObject*)0,Form("Cut: x > %s", var.Data()));
    legend_r->AddEntry((TObject*)0, Form("AUC = %.4f", aucRight));
    legend_r->AddEntry((TObject*)0, Form("FPR = %.4f, TPR = %.4f", effBkg_right[bestBinRight], effSig_right[bestBinRight]), "");
    // make the text smaller
legend_r->SetTextSize(0.035);

    // shrink & move the legend box by tweaking its NDC corners
legend_r->SetX1NDC(0.65);  // left
legend_r->SetX2NDC(0.15);  // right
legend_r->SetY1NDC(0.80);  // bottom
legend_r->SetY2NDC(0.30);  // top

// optionally remove the fill to make it less obtrusive
//legend_r->SetFillStyle(0);
//legend_r->SetBorderSize(0);

   legend_r->Draw();


            cRight->SaveAs(Form("./%s%s%s_%s%s_right.png", dirNAME.Data(), systemNAME.Data() , var.Data(), particleNAME.Data(), cutlevel.Data()));


        
        //=== Right-side cut: x < cut ===
        for (int ibin = 1; ibin <= nbins; ++ibin) {
           double S_pass = histogram->Integral(1, ibin);                                          //std::cout << "  S    = " << S << std::endl;
           double B_pass = hist_BKG->Integral(1, ibin);                                               //std::cout << "  Significance2 = " << significance << std::endl;       
           double effSig =  S_pass/S_total;   // True Positive Rate     
           double effBkg = (B_pass/B_total); //False Postive Rate
                                                                                                //std::cout << "  Cut Value2        = " << cut << std::endl;      
          effBkg_left.push_back(effBkg);
          effSig_left.push_back(effSig);

    }

double minDistLeft = 999.0;
int bestBinLeft = -1;

for (int i = 0; i < nbins; ++i) {
    double fpr = effBkg_left[i];
    double tpr = effSig_left[i];
    double dist = std::sqrt(std::pow(fpr, 2) + std::pow(1 - tpr, 2));

    if (dist < minDistLeft) {
        minDistLeft = dist;
        bestBinLeft = i;
    }
}
double cutValueLeft = histogram->GetBinLowEdge(bestBinLeft + 1);
if (bestBinLeft >= 0) {
    double cutValueLeft = histogram->GetBinLowEdge(bestBinLeft + 1);
    std::cout << "\n[LEFT CUT] Best ROC Point (Closest to (0,1)):" << std::endl;
    std::cout << "  Cut Direction: x < " << cutValueLeft << std::endl;
    std::cout << "  FPR: " << effBkg_left[bestBinLeft] << ", TPR: " << effSig_left[bestBinLeft] << std::endl;
    std::cout << "  Distance to (0,1): " << minDistLeft << std::endl;
}


// Plot left-side ROC
TCanvas* cLeft = new TCanvas("cLeft", "", 800, 600);
cLeft->SetGrid();

TGraph* rocLeft = new TGraph(effSig_left.size(), &effBkg_left[0], &effSig_left[0]);
rocLeft->SetTitle(";False Postive Rate;True Positive Rate");
rocLeft->SetLineColor(kRed);
rocLeft->SetLineWidth(2);

rocLeft->GetXaxis()->SetLimits(0.0, 1.0);
rocLeft->GetYaxis()->SetRangeUser(0.0, 1.0);
rocLeft->Draw("AL");

// --- Marker + label at chosen point ---
if (bestBinLeft >= 0) {
    double xL = effBkg_left[bestBinLeft];
    double yL = effSig_left[bestBinLeft];

    TMarker* bestPtLeft = new TMarker(xL, yL, 20);
    bestPtLeft->SetMarkerColor(kRed+1);
    bestPtLeft->SetMarkerSize(1.4);
    bestPtLeft->Draw("SAME");

    TLatex labL;
    labL.SetTextSize(0.028);
    labL.SetTextColor(kRed+1);
    double xl = std::min(0.95, std::max(0.05, xL + 0.02));
    double yl = std::min(0.95, std::max(0.05, yL - 0.04));
    labL.DrawLatex(xl, yl, Form("cut = %.5g", cutValueLeft));
}



std::vector<std::pair<double, double>> rocPoints_left;
for (size_t i = 0; i < effBkg_left.size(); ++i)
    rocPoints_left.emplace_back(effBkg_left[i], effSig_left[i]);

std::sort(rocPoints_left.begin(), rocPoints_left.end());  // sort by FPR (x)

// Reassign sorted values
for (size_t i = 0; i < rocPoints_left.size(); ++i) {
    effBkg_left[i] = rocPoints_left[i].first;
    effSig_left[i] = rocPoints_left[i].second;
}


// --- Compute AUC (Left-side cut) ---
double aucLeft = 0.0;
for (size_t i = 1; i < effBkg_left.size(); ++i) {
    double x1 = effBkg_left[i - 1], x2 = effBkg_left[i];
    double y1 = effSig_left[i - 1], y2 = effSig_left[i];
    aucLeft += 0.5 * (x2 - x1) * (y1 + y2); // Trapezoid rule
}

// Draw AUC on canvas
/*TLatex latexLeft;
latexLeft.SetNDC();
latexLeft.SetTextSize(0.03);
latexLeft.DrawLatex(0.2, 0.88, Form("AUC = %.4f", aucLeft));*/

auto legend = new TLegend(0.6, 0.75, 0.88, 0.88);
    legend->SetHeader(particleNAME, "");
    legend->AddEntry((TObject*)0,Form("Cut: x < %s", var.Data()));
    legend->AddEntry((TObject*)0, Form("AUC = %.4f", aucLeft));
    legend->AddEntry((TObject*)0, Form("FPR = %.4f, TPR = %.4f", effBkg_left[bestBinLeft], effSig_left[bestBinLeft]));
    // make the text smaller
legend->SetTextSize(0.035);  

    // shrink & move the legend box by tweaking its NDC corners
legend->SetX1NDC(0.65);  // left
legend->SetX2NDC(0.15);  // right
legend->SetY1NDC(0.80);  // bottom
legend->SetY2NDC(0.30);  // top

// optionally remove the fill to make it less obtrusive
//legend->SetFillStyle(0);  
//legend->SetBorderSize(0);

   legend->Draw();

cLeft->SaveAs(Form("./%s%s%s_%s%s_left.png", dirNAME.Data(), systemNAME.Data() , var.Data(), particleNAME.Data(), cutlevel.Data()));


       
        




        /*std::cout << "  Variable  = " << var.Data() << std::endl;
        std::cout << "  B = " << B_Cut1 << std::endl;
        std::cout << "  S    = " << S_Cut1 << std::endl;
        std::cout << "  Max Significance1 = " << maxSig1 << std::endl;
        std::cout << "  Cut Direction1    = " << bestDirection1 << std::endl;
        std::cout << "  Cut Value1        = " << bestCut1 << std::endl;


        std::cout << "  B = " << B_Cut2 << std::endl;
        std::cout << "  S    = " << S_Cut2 << std::endl;
        std::cout << "  Max Significance2 = " << maxSig2 << std::endl;
        std::cout << "  Cut Direction    = " << bestDirection2 << std::endl;
        std::cout << "  Cut Value        = " << bestCut2 << std::endl;

        std::cout << "\n=== Best Significance ===" << std::endl;
        std::cout << "  Max Significance = " << maxSig_both << std::endl;
        std::cout << "  Cut Direction    = " << bestDirection << std::endl;
        std::cout << "  Cut Value        = " << bestCut_both << std::endl;*/
     
    


        /*double x_cut = bestCut_both;        
        double y_min = 0;
        double y_max = hist_BKG->GetMaximum();*/

       /* TLine* line = new TLine(x_cut, y_min, x_cut, y_max);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->SetLineStyle(2);  // dashed line
        line->Draw("SAME");*/

        /*TString particleNAME = "Bu";
        TString systemNAME = "SG_ppRef_";
        if (path_to_file.Contains("Bs")){
            particleNAME = "Bs";
        } else if (path_to_file.Contains("Bd")){
            particleNAME = "Bd";
        } else if (path_to_file.Contains("Rho")){
            particleNAME = "X3872";
        } else if (path_to_file.Contains("PSI2S")){
            particleNAME = "PSI2S";
        }
        if (path_to_file.Contains("PbPb"))  { systemNAME = "SG_PbPb_";}*/


          // === Plotting ===
        /*TCanvas* c = new TCanvas("c", "", 900, 700);
        c->SetGrid();
    
        double xMin_R = cuts_right.front();
        double xMax_R = cuts_right.back();
        double xMin_L = cuts_left.front();
        double xMax_L = cuts_left.back();

        double yMax_R = *std::max_element(sig_right.begin(), sig_right.end());

        TGraph* gRight = new TGraph(cuts_right.size(), &cuts_right[0], &sig_right[0]);
        gRight->SetLineColor(kBlue);
        gRight->SetLineWidth(2);
    // gRight->SetTitle("Significance Vs Cut Value");
        gRight->GetXaxis()->SetTitle(Form("Cut Value, %s", var.Data()));
        gRight->GetYaxis()->SetTitle("Significance");

        gRight->GetXaxis()->SetTitleSize(0.045);
        gRight->GetXaxis()->SetTitleOffset(1.1);
        gRight->GetYaxis()->SetTitleSize(0.045);
        gRight->GetYaxis()->SetTitleOffset(1.2);


        TGraph* gLeft = new TGraph(cuts_left.size(), &cuts_left[0], &sig_left[0]);
        gLeft->SetLineColor(kRed);
        gLeft->SetLineWidth(2);

        gRight->GetXaxis()->SetLimits(xMin_R, xMax_R);
        gRight->GetYaxis()->SetRangeUser(0, maxSig_both*1.1);
        gLeft->GetXaxis()->SetLimits(xMin_L, xMax_L);
        gLeft->GetYaxis()->SetRangeUser(0, maxSig_both*1.1);

        gRight->Draw("AL");
        gLeft->Draw("L SAME");

    // draw a vertical dashed line at the best cut
TLine* bestLine = new TLine(bestCut_both, 0, bestCut_both, maxSig_both*1.1);
bestLine->SetLineColor(kOrange+7);
bestLine->SetLineStyle(2);
bestLine->SetLineWidth(2);
bestLine->Draw("SAME");

// draw a marker at the maximum-significance point
TMarker* bestMark = new TMarker(bestCut_both, maxSig_both, 29);
bestMark->SetMarkerColor(kMagenta);
bestMark->SetMarkerSize(1.2);
bestMark->Draw("SAME");

    auto legend = new TLegend(0.6, 0.75, 0.88, 0.88);
    legend->SetHeader(particleNAME, "");
    legend->AddEntry(gRight, "Cut: x > value", "l");
    legend->AddEntry(gLeft,  "Cut: x < value", "l");

    // make the text smaller
legend->SetTextSize(0.025);  

    // shrink & move the legend box by tweaking its NDC corners
legend->SetX1NDC(0.65);  // left
legend->SetX2NDC(0.85);  // right
legend->SetY1NDC(0.80);  // bottom
legend->SetY2NDC(0.90);  // top

// optionally remove the fill to make it less obtrusive
legend->SetFillStyle(0);  
legend->SetBorderSize(0);

        legend->Draw();*/
      
     
    //c->SaveAs(Form("./%s%s%s_%s%s.pdf", dirNAME.Data(), systemNAME.Data() , var.Data(), particleNAME.Data(), cutlevel.Data()));
    

       // canvas->SaveAs(Form("./%s%s%s_%s%ssig.pdf", dirNAME.Data(), systemNAME.Data() , var.Data(), particleNAME.Data(), cutlevel.Data()));
        
      
        
      
      //  c->SaveAs("significance_scan_both.pdf");
        delete hist_SIG;
        delete hist_SIG_WT;
        delete hist_SIG_BOTH;
        delete hist_BKG;
        delete hist;
        delete cLeft;
        delete cRight;
        
    }
    }

 }

        // Clean up
       
       
  



int main() {
    EffSig();
    return 0;
}
