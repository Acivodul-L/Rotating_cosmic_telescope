#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath> // Necessario per sqrt e pow
#include <algorithm> // Necessario per std::max
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <cmath>
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TString.h"
#include "TPaveStats.h"


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
void make_hist(const char* fileUp = "NANO_90_UP.txt", const char* fileDown = "NANO_90_DOWN.txt") {
    
    gStyle->SetOptStat(0); 

    // --- CONFIGURAZIONE ASSI ---
    const double BIN_WIDTH = 5.0;   
    const double MAX_RANGE = 100.0; 

    int n_bins = (int)((2 * MAX_RANGE) / BIN_WIDTH) + 1;
    double min_val = -(n_bins * BIN_WIDTH) / 2.0;
    double max_val =  (n_bins * BIN_WIDTH) / 2.0;

    std::cout << "--- Configurazione Binning ---" << std::endl;
    std::cout << "Bin totali: " << n_bins << " (Range: " << min_val << " to " << max_val << " ns)" << std::endl;

    TH1F *hUp = new TH1F("hUp", "Confronto Coincidenze (Normalizzato)", n_bins, min_val, max_val);
    TH1F *hDown = new TH1F("hDown", "Confronto Coincidenze (Normalizzato)", n_bins, min_val, max_val);

    // Stile hUp
    hUp->SetLineColor(kBlue);
    hUp->SetLineWidth(2);

    // Stile hDown
    hDown->SetLineColor(kRed);
    hDown->SetLineWidth(2);
    hDown->SetLineStyle(2); 

    // --- ELABORAZIONE ---
    process_file(fileUp, hUp);
    process_file(fileDown, hDown);

    // --- CALCOLO STATISTICHE (PRIMA DELLA NORMALIZZAZIONE) ---
    // È importante calcolare mean, rms ed errori sui dati grezzi
    
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

    // Calcolo Tempo di Volo (Differenza delle medie / 2)
    double timeOfFlight = (meanUp - meanDown)/2.0;

    // Calcolo Errore sul ToF (Somma in quadratura degli errori sulle medie)
    double errTimeOfFlight = std::sqrt(std::pow(errMeanUp, 2) + std::pow(errMeanDown, 2))/2.0;


    // --- STAMPA RISULTATI DETTAGLIATI ---
    std::cout << "\n========================================" << std::endl;
    std::cout << "            RISULTATI ANALISI             " << std::endl;
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
    std::cout << " TEMPO DI VOLO (Media UP - Media DOWN)/2: " << std::endl;
    std::cout << "   " << timeOfFlight << " +/- " << errTimeOfFlight << " ns" << std::endl;
    std::cout << "========================================\n" << std::endl;


    // --- NORMALIZZAZIONE DEGLI ISTOGRAMMI ---
    // Scaliamo l'area degli istogrammi a 1 per confrontare la "forma"
    // indipendentemente dal numero di eventi.
    if (hUp->Integral() > 0)   hUp->Scale(1.0 / hUp->Integral());
    if (hDown->Integral() > 0) hDown->Scale(1.0 / hDown->Integral());


    // --- DISEGNO ---
    TCanvas *c1 = new TCanvas("c1", "Confronto Nanosecondi", 1000, 600);
    
    // Calcoliamo il massimo valore Y tra i due istogrammi normalizzati
    // e aggiungiamo un margine (es. 20%) per l'estetica.
    double maxVal = std::max(hUp->GetMaximum(), hDown->GetMaximum());
    hUp->GetYaxis()->SetRangeUser(0, maxVal * 1.2);
    
    // Titoli assi aggiornati
    hUp->GetXaxis()->SetTitle("Differenza Temporale (ns)");
    hUp->GetYaxis()->SetTitle("Frazione di Eventi (Normalizzati)");

    // Disegno
    hUp->Draw("HIST");        
    hDown->Draw("HIST SAME"); 

    // --- LEGENDA ---
    TLegend *leg = new TLegend(0.55, 0.70, 0.89, 0.89); // Leggermente più largo
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // OPZIONE B: Se vuoi anche il numero di eventi (N), usa questa versione:
    TString labelUp = TString::Format("UP (#mu=%.2f, #sigma=%.2f, N=%.0f)", meanUp, rmsUp, nUp);
    TString labelDown = TString::Format("DOWN (#mu=%.2f, #sigma=%.2f, N=%.0f)", meanDown, rmsDown, nDown);

    leg->AddEntry(hUp, labelUp, "l");     
    leg->AddEntry(hDown, labelDown, "l");
    leg->Draw();
    
    std::cout << "Grafico normalizzato generato." << std::endl;
}

void make_plot(){

  double times_fly[] = {4.58, 4.35, 3.92, 0.54};
  double times_fly_err[] = {0.14, 0.06, 0.11, 0.16};
  double angles[] = {0.0, 1.0/6.0, 1.0/3.0, 1.0/2.0};
  double angles_err[] = {0.0, 0.0, 0.0, 0.0};
  
  const int n = 4;

  auto *fly_times = new TGraphErrors(n, angles, times_fly, angles_err, times_fly_err);

    fly_times->SetTitle("Distribuzione TOF;Angolo [rad/pi];TOF [ns]");
    fly_times->SetMarkerStyle(21); 
    fly_times->SetMarkerColor(kBlue);
    fly_times->SetLineColor(kBlue);
   
    fly_times->GetXaxis()->SetLimits(-0.1, (1.0/2.0)*1.1); 
    fly_times->GetYaxis()->SetRangeUser(0.0, 5.0);
    
   fly_times->Draw("ALP");
}



