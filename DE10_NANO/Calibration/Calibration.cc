//.L Calibration.cc
//AnalyzeCalibration()

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"

void AnalyzeCalibration() {
    // Vettori per memorizzare i dati per il grafico
    std::vector<double> x_delay_ns; // Asse X: Ritardo indotto in ns
    std::vector<double> y_cycles;   // Asse Y: Ritardo letto in cicli
    std::vector<double> y_errors;   // Errore sulla media (SEM)

    int fileCounter = 1;

    // Loop per leggere i file Calibration_1, Calibration_2, ecc.
    while (true) {
        TString filename = Form("Calibration_%d.txt", fileCounter);
        std::ifstream infile(filename.Data());

        if (!infile.is_open()) {
            if (fileCounter == 1) {
                std::cout << "Errore: Nessun file 'Calibration_1.txt' trovato." << std::endl;
                return;
            }
            break; // Fine dei file
        }

        // Variabili per la lettura
        int id;
        long long time_val;
        int prev_id = -1;
        long long prev_time = 0;
        
        std::vector<double> diffs;

        // Legge il file e calcola le differenze (T_tag1 - T_tag2)
        while (infile >> id >> time_val) {
            // Se troviamo un 1 ed è preceduto da un 2
            if (id == 1 && prev_id == 2) {
                double diff = (double)(time_val - prev_time);
                diffs.push_back(diff);
            }
            prev_id = id;
            prev_time = time_val;
        }
        infile.close();

        // Analisi statistica
        double n = diffs.size();
        if (n > 0) {
            double sum = 0;
            for (double val : diffs) sum += val;
            double mean = sum / n;

            double sum_sq_diff = 0;
            for (double val : diffs) sum_sq_diff += (val - mean) * (val - mean);
            double std_dev = (n > 1) ? std::sqrt(sum_sq_diff / (n - 1)) : 0.0;
            double std_err = (n > 0) ? std_dev / std::sqrt(n) : 0.0;

            // Calcolo X: Calibration_1 -> 5ns, Calibration_2 -> 10ns, ecc.
            double delay_val = 5.0 + fileCounter * 5.0;

            x_delay_ns.push_back(delay_val);
            y_cycles.push_back(mean);
            y_errors.push_back(std_err);
            
            std::cout << "File " << filename << " (" << delay_val << " ns): " 
                      << "Media = " << mean << " cicli +/- " << std_err << std::endl;
        }
        fileCounter++;
    }

    if (x_delay_ns.empty()) return;

    // --- Creazione Grafico ---
    TCanvas *c1 = new TCanvas("c1", "Clock Calibration", 800, 600);
    c1->SetGrid();

    TGraphErrors *gr = new TGraphErrors(x_delay_ns.size(), 
                                        &x_delay_ns[0], 
                                        &y_cycles[0], 
                                        0, // Nessun errore su X (assunto esatto)
                                        &y_errors[0]);

    gr->SetTitle("Verifica Durata Clock;Ritardo indotto (ns);Ritardo letto (in cicli di clock)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);

    // --- Fit Lineare ---
    // Usiamo pol1: p0 + p1*x
    TF1 *fitFunc = new TF1("fitFunc", "pol1", 0, x_delay_ns.back() + 5);
    fitFunc->SetLineColor(kRed);
    fitFunc->SetParName(0, "Offset (cicli)");
    fitFunc->SetParName(1, "Slope (cicli/ns)");
    
    std::cout << "\n--- Esecuzione Fit Lineare ---" << std::endl;
    gr->Fit(fitFunc);
    gr->Draw("AP");

    // --- Calcolo del Periodo del Clock ---
    // Slope (m) = cicli / ns
    // Periodo (T) = ns / ciclo = 1 / m
    double slope = fitFunc->GetParameter(1);
    double slope_err = fitFunc->GetParError(1);

    if (slope != 0) {
        double clock_period = 1.0 / slope;
        // Propagazione errore: sigma_T = T^2 * sigma_m
        double clock_period_err = (clock_period * clock_period) * slope_err;

        std::cout << "\n========================================" << std::endl;
        std::cout << " RISULTATI VERIFICA CLOCK" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Pendenza (m)  : " << slope << " cicli/ns" << std::endl;
        std::cout << "Periodo stimato: " << clock_period << " +/- " << clock_period_err << " ns" << std::endl;
        
        double expected = 5.0;
        double sigma_diff = std::abs(clock_period - expected) / clock_period_err;
        
        std::cout << "Valore atteso  : " << expected << " ns" << std::endl;
        std::cout << "Discrepanza    : " << sigma_diff << " sigma" << std::endl;
        std::cout << "========================================" << std::endl;

        // Aggiunge una legenda con il risultato
        TLegend *leg = new TLegend(0.15, 0.75, 0.55, 0.88);
        leg->AddEntry(gr, "Dati misurati", "ep");
        leg->AddEntry(fitFunc, Form("Fit: T = %.3f #pm %.3f ns", clock_period, clock_period_err), "l");
        leg->Draw();
    } else {
        std::cout << "Attenzione: Pendenza zero, impossibile calcolare il periodo." << std::endl;
    }
}