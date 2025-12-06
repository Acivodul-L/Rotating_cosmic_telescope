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

// --- DEFINIZIONE FUNZIONE DI FIT ---
Double_t fitFunction(Double_t *x, Double_t *par) {
    Double_t angle = x[0];
    Double_t A = par[0];
    Double_t b = par[1];
    Double_t c = par[2];
    Double_t d = par[3];

    return A * TMath::Power(TMath::Abs(TMath::Cos(angle + b)), c) + d;
}

// --- FUNZIONE PRINCIPALE ---
void fit_flux() {
    // 1. Impostazioni Grafiche
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);

    // 2. DEFINIZIONE DATI (Tutte le variabili devono essere dichiarate qui)
    const int n_points = 10;
    
    // Dati grezzi
    double pairs_raw[] = {3464, 6025, 10708, 18842, 3272, 2677, 25397, 11643, 5266, 3278};
    double triples_raw[] = {1513, 4519, 8404, 15843, 2785, 2267, 21100, 9321, 3819, 1442};
    
    // Conteggi singoli per il calcolo del buio
    double singles_1[] = {1046157, 591662, 673483, 0, 108863, 0, 0, 765830, 0, 1043397}; //MANCANO DEI DATI
    double singles_2[] = {977130, 554116, 641943, 0, 104272, 0, 0, 774706, 0, 928795};
    double singles_4[] = {5494318, 3050526, 3189050, 0, 448241, 0, 0, 3345766, 0, 5114312};
    
    // Tempi in minuti
    double times_min[] = {1680, 855, 870, 1065, 144, 122, 1678, 1020, 840, 1620}; 
    
    // Angoli in radianti
    double angles[] = {
        -TMath::Pi()/2, -TMath::Pi()/3, -TMath::Pi()/4, -TMath::Pi()/6, 
        0.0, 
        TMath::Pi()/12, TMath::Pi()/6, TMath::Pi()/4, TMath::Pi()/3, TMath::Pi()/2
    };

    // 3. COSTANTI FISICHE E GEOMETRICHE
    double sigma = 0.043;   // Steradianti
    double S = 0.2 * 0.4;   // Area m^2

    // Efficienze (Stima)
    double eff_1 = 0.95; 
    double eff_4 = 0.95; 
    double eff_coinc = eff_1 * eff_4; // Efficienza totale coincidenza

    // Finestre temporali per accidentali (in secondi)
    double width_disc_1 = 51e-9; 
    double width_disc_4 = 78e-9; 
    double coincidence_window = width_disc_1 + width_disc_4; 

    // 4. VETTORI DI SUPPORTO (Devono essere definiti prima del loop)
    std::vector<double> x_val(n_points);
    std::vector<double> x_err(n_points, 0.0);
    
    std::vector<double> y_pairs(n_points);
    std::vector<double> y_pairs_err(n_points);
    
    std::vector<double> y_triples(n_points);
    std::vector<double> y_triples_err(n_points);

    std::vector<double> pairs_dark(n_points);
    std::vector<double> pairs_dark_err(n_points);

    // 5. CALCOLO FLUSSI
    std::cout << "--- Inizio Analisi ---" << std::endl;
    std::cout << "Efficienza Coincidenza usata: " << eff_coinc << std::endl;

    for(int i = 0; i < n_points; i++) {
        double time_sec = times_min[i] * 60.0;
        x_val[i] = angles[i];

        // --- A. Calcolo Buio (Accidentali) ---
        if (singles_1[i] > 0 && singles_4[i] > 0 && time_sec > 0) {
            // Rate accidentale = (R1 * R4 * deltaT) -> N_acc = (N1 * N4 * deltaT) / T
            pairs_dark[i] = (singles_1[i] * singles_4[i] * coincidence_window) / time_sec;
            
            // Errore relativo accidentali
            double rel_err = sqrt(1.0/singles_1[i] + 1.0/singles_4[i]);
            pairs_dark_err[i] = pairs_dark[i] * rel_err;
        } else {
            pairs_dark[i] = 0.0;
            pairs_dark_err[i] = 0.0;
        }

        // --- B. Calcolo Flusso Pairs ---
        // Denominatore comune che include Tempo, Accettanza ed Efficienza
        double denom = time_sec * sigma * S * eff_coinc;

        if (denom > 0) {
            y_pairs[i] = (pairs_raw[i] - pairs_dark[i]) / denom;
            
            // Errore: somma in quadratura di Poisson (raw) e errore accidentali
            double err_raw_sq = pairs_raw[i]; 
            double err_dark_sq = pow(pairs_dark_err[i], 2.0); // Nota: 2.0 per evitare errori pow
            
            y_pairs_err[i] = sqrt(err_raw_sq + err_dark_sq) / denom;
        } else {
            y_pairs[i] = 0;
            y_pairs_err[i] = 0;
        }

        // --- C. Calcolo Flusso Triples ---
        // Nota: Qui non applichiamo l'efficienza coincidenza doppia, 
        // ma idealmente servirebbe l'efficienza tripla. Lascio base per ora.
        double denom_trip = time_sec * sigma * S; // * eff_triple (se nota)
        
        if (denom_trip > 0) {
            y_triples[i] = triples_raw[i] / denom_trip;
            
            double term_inside = triples_raw[i] * (1.0 - triples_raw[i] / (pairs_raw[i] > 0 ? pairs_raw[i] : 1.0));
            if(term_inside < 0) term_inside = 0; 
            y_triples_err[i] = sqrt(term_inside) / denom_trip;
        }
    }

    // 6. CREAZIONE GRAFICI
    TGraphErrors *gr_pairs = new TGraphErrors(n_points, x_val.data(), y_pairs.data(), x_err.data(), y_pairs_err.data());
    TGraphErrors *gr_triples = new TGraphErrors(n_points, x_val.data(), y_triples.data(), x_err.data(), y_triples_err.data());

    // 7. FIT
    TF1 *f_pairs = new TF1("f_pairs", fitFunction, -TMath::Pi()/2, TMath::Pi()/2, 4);
    f_pairs->SetParNames("A", "b", "c", "d");
    f_pairs->SetLineColor(kBlue);
    f_pairs->SetLineStyle(2); 
    f_pairs->SetParameters(100, 0.0, 3.5, 20);

    TF1 *f_triples = new TF1("f_triples", fitFunction, -TMath::Pi()/2, TMath::Pi()/2, 4);
    f_triples->SetParNames("A", "b", "c", "d");
    f_triples->SetLineColor(kRed);
    f_triples->SetLineStyle(2);
    f_triples->SetParameters(100, 0.0, 2.5, 2.0);

    std::cout << "\n--- Fit Pairs ---" << std::endl;
    gr_pairs->Fit(f_pairs, "SMER");

    std::cout << "\n--- Fit Triples ---" << std::endl;
    gr_triples->Fit(f_triples, "SMER");

    // 8. DISEGNO
    TCanvas *c1 = new TCanvas("c1", "Fit Flusso vs Angolo", 1000, 600);
    c1->SetGrid();

    gr_pairs->SetTitle("Fit del Flusso;Angolo [rad];Flusso [eventi/s/str/m^{2}]");
    gr_pairs->SetMarkerStyle(21); 
    gr_pairs->SetMarkerColor(kBlue);
    gr_pairs->SetLineColor(kBlue);
    
    // Pairs
    gr_pairs->Draw("AP");

    // Triples
    gr_triples->SetMarkerStyle(22); 
    gr_triples->SetMarkerColor(kRed);
    gr_triples->SetLineColor(kRed);
    gr_triples->Draw("P SAME");

    // Legenda
    TLegend *leg = new TLegend(0.15, 0.7, 0.45, 0.85);
    leg->AddEntry(gr_pairs, "Pairs (Corrected)", "ep");
    leg->AddEntry(f_pairs, "Fit Pairs", "l");
    leg->AddEntry(gr_triples, "Triples", "ep");
    leg->AddEntry(f_triples, "Fit Triples", "l");
    leg->Draw();
    
    // Sistemazione Box Statistiche
    gPad->Update();
    TPaveStats *st1 = (TPaveStats*)gr_pairs->FindObject("stats");
    if(st1) {
        st1->SetY1NDC(0.7); st1->SetY2NDC(0.9);
        st1->SetTextColor(kBlue);
    }
    TPaveStats *st2 = (TPaveStats*)gr_triples->FindObject("stats");
    if(st2) {
        st2->SetY1NDC(0.5); st2->SetY2NDC(0.7);
        st2->SetTextColor(kRed);
    }
    gPad->Modified();
    gPad->Update();
}