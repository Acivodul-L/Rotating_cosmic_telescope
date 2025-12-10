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
    // Dati grezzi
    double pairs_raw[] = {3464, 6025, 10708, 18842, 3574, 3272, 2677, 25397, 11643, 5266, 3278};
    double triples_raw[] = {1513, 4519, 8404, 15843, 3042, 2785, 2267, 21100, 9321, 3819, 1442};
    
    // Conteggi singoli per il calcolo del buio
    double singles_1[] = {1046157, 591662, 673483, 893636, 144500, 108863, 108863, 1294018, 765830, 652622, 1043397};
    double singles_2[] = {977130, 554116, 641943, 849766, 138248, 104272, 104272, 1240390, 774706, 567167, 928795};
    double singles_4[] = {5494318, 3050526, 3189050, 3918389, 603573, 448241, 448241, 5737493, 3345766, 2861250, 5114312};
    
    // Tempi in minuti
    double times_min[] = {1680, 855, 870, 1065, 170, 144, 122, 1678, 1020, 840, 1620}; 
    
    // Angoli in radianti
    double angles[] = {
        -TMath::Pi()/2, -TMath::Pi()/3, -TMath::Pi()/4, -TMath::Pi()/6, -TMath::Pi()/12,
        0.0, 
        TMath::Pi()/12, TMath::Pi()/6, TMath::Pi()/4, TMath::Pi()/3, TMath::Pi()/2
    };

    //sizeof(pairs_raw) è la grandezza in byte di tutto l'array
    // sizeof(pairs_raw[0]) è la grandezza in byte di un solo double
    const int n_points = sizeof(times_min) / sizeof(times_min[0]);
    
    // 3. COSTANTI FISICHE E GEOMETRICHE
    double sigma = 0.043;   // Steradianti
    double S = 0.2 * 0.4;   // Area m^2

    // Efficienze (Stima)
    double eff_1 = 0.95; 
    double eff_2 = 0.95;
    double eff_4 = 0.95; 
    double eff_err_1 = 0.07;
    double eff_err_2 = 0.07;
    double eff_err_4 = 0.07;
    
    double eff_coinc = eff_1 * eff_4; // Efficienza totale coincidenza doppia
    double eff_coinc_err = (eff_coinc)*pow((pow(eff_err_1/eff_1, 2.0) + pow((eff_err_4/eff_4),2.0)), 0.5); //Errore sull'efficienza di doppia
    
    double eff_coinc_trip = eff_1 * eff_2 * eff_4;
    double eff_coinc_trip_err = (eff_coinc_trip)*pow(pow(eff_err_1/eff_1, 2.0) + (pow(eff_err_2/eff_2, 2.0) + pow((eff_err_4/eff_4),2.0)), 0.5);
  

    // Finestre temporali per accidentali (in secondi)
    double width_disc_1 = 51e-9; 
    double width_disc_4 = 78e-9; 
    double coincidence_window = width_disc_1 + width_disc_4 - pow(4, -9); 

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
        
        y_pairs[i] = (pairs_raw[i] - pairs_dark[i])/ denom;    
        y_pairs_err[i] =  y_pairs[i] * sqrt((1/pairs_raw[i]) + pow(pairs_dark_err[i]/pairs_dark[i] , 2.0) + pow((eff_coinc_err/eff_coinc), 2.0));
  

        // --- C. Calcolo Flusso Triples ---
        // Nota: Qui ignoriamo le coincidenze triple di buio, se erano trascurabili prima lo saranno anche ora 
        double denom_trip = time_sec * sigma * S * eff_coinc_trip;
        
        y_triples[i] = (triples_raw[i])/ denom;    
        y_triples_err[i] =  y_triples[i] * sqrt((1/triples_raw[i])+ pow((eff_coinc_trip_err/eff_coinc_trip), 2.0));
        
    }

    // 6. CREAZIONE GRAFICI
    TGraphErrors *gr_pairs = new TGraphErrors(n_points, x_val.data(), y_pairs.data(), x_err.data(), y_pairs_err.data());
    TGraphErrors *gr_triples = new TGraphErrors(n_points, x_val.data(), y_triples.data(), x_err.data(), y_triples_err.data());

    // 7. FIT
    TF1 *f_pairs = new TF1("f_pairs", fitFunction, -TMath::Pi()/2-0.5, TMath::Pi()/2 + 0.5, 4);
    f_pairs->SetParNames("A", "b", "c", "d");
    f_pairs->SetLineColor(kBlue);
    f_pairs->SetLineStyle(2); 
    f_pairs->SetParameters(100, 0.0, 3.5, 20);

    TF1 *f_triples = new TF1("f_triples", fitFunction, -TMath::Pi()/2-1, TMath::Pi()/2 + 1, 4);
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

    // --- MODIFICA QUI: Impostazione Limiti ---
    
    // 1. Allarghiamo l'asse X per vedere bene i punti a +/- Pi/2
    // Impostiamo da -2.0 a +2.0 (Pi/2 è circa 1.57, quindi ci stiamo larghi)
    gr_pairs->GetXaxis()->SetLimits(-2.0, 2.0); 

    // 2. Impostiamo l'asse Y (Minimo e Massimo)
    // Esempio: da 0 a 160 (o un valore poco sopra il tuo massimo flusso)
    gr_pairs->GetYaxis()->SetRangeUser(0.0, 160.0);

    // Titoli assi (opzionale ma consigliato per chiarezza)
    gr_pairs->GetXaxis()->SetTitle("Angolo [rad]");
    gr_pairs->GetYaxis()->SetTitle("Flusso [Hz/sr/m^2]");

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
