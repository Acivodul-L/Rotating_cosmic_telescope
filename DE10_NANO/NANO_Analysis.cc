#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath> // Necessario per sqrt e pow
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"

// =============================================================
// PARAMETRI GLOBALI
// =============================================================
const double NS_PER_CLOCK = 5.0; 

// =============================================================
// FUNZIONE DI LETTURA
// =============================================================
void process_file(const char* filename, TH1F* histogram) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cout << "Errore critico: Impossibile aprire il file " << filename << std::endl;
        return;
    }

    Int_t epoch = 0; 
    const Long64_t MARKER_KEY = 2147483648LL; 
    const Long64_t CLOCK_WRAP = 1073741824LL; 

    Long64_t last_time_pmt1 = -1;
    Long64_t last_time_pmt4 = -1;

    std::string line;
    Long64_t col1, col2;

    while (std::getline(infile, line)) {
        if (line.empty() || line.find("//") == 0) continue;
        
        std::stringstream ss(line);
        if (!(ss >> col1 >> col2)) continue;

        if (col1 == MARKER_KEY) {
            epoch = col2 - MARKER_KEY;
        } else {
            Long64_t current_global_time = ((Long64_t)epoch * CLOCK_WRAP) + (Long64_t)col2;
            Int_t mask = (Int_t)col1;

            if (mask == 3) {
                last_time_pmt1 = current_global_time;
                last_time_pmt4 = current_global_time;
            } else if (mask == 1) {
                last_time_pmt1 = current_global_time;
            } else if (mask == 2) {
                last_time_pmt4 = current_global_time;
            }
            
            if (mask == 4) {
                if (last_time_pmt1 != -1 && last_time_pmt4 != -1) {
                    Long64_t diff_cycles = last_time_pmt4 - last_time_pmt1;
                    double diff_ns = (double)diff_cycles * NS_PER_CLOCK;
                    histogram->Fill(diff_ns);
                }
            }
        }
    }
    infile.close();
}

// =============================================================
// FUNZIONE PRINCIPALE
// =============================================================
void make_hist(const char* fileUp = "NANO_0_UP.txt", const char* fileDown = "NANO_0_DOWN.txt") {
    
    gStyle->SetOptStat(0); 

    // --- CONFIGURAZIONE ASSI ---
    const double BIN_WIDTH = 5.0;   
    const double MAX_RANGE = 25.0; 

    int n_bins = (int)((2 * MAX_RANGE) / BIN_WIDTH) + 1;
    double min_val = -(n_bins * BIN_WIDTH) / 2.0;
    double max_val =  (n_bins * BIN_WIDTH) / 2.0;

    std::cout << "--- Configurazione Binning ---" << std::endl;
    std::cout << "Bin totali: " << n_bins << " (Range: " << min_val << " to " << max_val << " ns)" << std::endl;

    TH1F *hUp = new TH1F("hUp", "Confronto Coincidenze (in ns)", n_bins, min_val, max_val);
    TH1F *hDown = new TH1F("hDown", "Confronto Coincidenze (in ns)", n_bins, min_val, max_val);

    hUp->SetLineColor(kBlue);
    hUp->SetLineWidth(2);

    hDown->SetLineColor(kRed);
    hDown->SetLineWidth(2);
    hDown->SetLineStyle(2); 

    // --- ELABORAZIONE ---
    process_file(fileUp, hUp);
    process_file(fileDown, hDown);

    // --- CALCOLO TEMPO DI VOLO E ERRORI ---
    
    // Recupero valori base
    double meanUp = hUp->GetMean();
    double rmsUp  = hUp->GetStdDev();
    double nUp    = hUp->GetEntries();

    double meanDown = hDown->GetMean();
    double rmsDown  = hDown->GetStdDev();
    double nDown    = hDown->GetEntries();

    // Calcolo Errore sulla Media (Standard Error of Mean = RMS / sqrt(N))
    double errMeanUp = 0;
    if (nUp > 0) errMeanUp = rmsUp / std::sqrt(nUp);

    double errMeanDown = 0;
    if (nDown > 0) errMeanDown = rmsDown / std::sqrt(nDown);

    // Calcolo Tempo di Volo (Differenza delle medie)
    double timeOfFlight = (meanUp - meanDown)/2;

    // Calcolo Errore sul ToF (Somma in quadratura degli errori sulle medie)
    double errTimeOfFlight = std::sqrt(std::pow(errMeanUp, 2) + std::pow(errMeanDown, 2))/2;


    // --- STAMPA RISULTATI DETTAGLIATI ---
    std::cout << "\n========================================" << std::endl;
    std::cout << "           RISULTATI ANALISI            " << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::cout << "FILE UP:" << std::endl;
    std::cout << "  Eventi (N): " << nUp << std::endl;
    std::cout << "  Media:      " << meanUp << " +/- " << errMeanUp << " ns" << std::endl;
    std::cout << "  RMS:        " << rmsUp << " ns" << std::endl;
    
    std::cout << "\nFILE DOWN:" << std::endl;
    std::cout << "  Eventi (N): " << nDown << std::endl;
    std::cout << "  Media:      " << meanDown << " +/- " << errMeanDown << " ns" << std::endl;
    std::cout << "  RMS:        " << rmsDown << " ns" << std::endl;

    std::cout << "----------------------------------------" << std::endl;
    std::cout << " TEMPO DI VOLO (Media UP - Media DOWN): " << std::endl;
    std::cout << "   " << timeOfFlight << " +/- " << errTimeOfFlight << " ns" << std::endl;
    std::cout << "========================================\n" << std::endl;


    // --- DISEGNO ---
    TCanvas *c1 = new TCanvas("c1", "Confronto Nanosecondi", 1000, 600);
    
    if (hUp->GetMaximum() > hDown->GetMaximum()) {
        hUp->Draw("HIST");       
        hDown->Draw("HIST SAME"); 
    } else {
        hDown->Draw("HIST");
        hUp->Draw("HIST SAME");
    }

    hUp->GetXaxis()->SetTitle("Differenza Temporale (ns)");
    hUp->GetYaxis()->SetTitle("Conteggi");

    // Legenda
    TLegend *leg = new TLegend(0.60, 0.70, 0.89, 0.89); 
    
    TString labelUp = TString::Format("UP (#mu=%.2f, #sigma=%.2f)", meanUp, rmsUp);
    TString labelDown = TString::Format("DOWN (#mu=%.2f, #sigma=%.2f)", meanDown, rmsDown);

    leg->AddEntry(hUp, labelUp, "l");    
    leg->AddEntry(hDown, labelDown, "l");
    leg->Draw();
    
    std::cout << "Grafico generato." << std::endl;
}
