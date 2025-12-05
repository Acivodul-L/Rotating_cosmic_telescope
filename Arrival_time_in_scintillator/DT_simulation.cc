//.L DT_simulation.cc
// Per eseguire senza TIR (default): ScintillatorSim(100000)
// Per eseguire CON TIR:             ScintillatorSim(100000, kTRUE)

#include "TCanvas.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TVector3.h"
#include "TPaveText.h"
#include <iostream>

/**
 * Simulazione Monte Carlo di uno scintillatore plastico.
 * * @param nPhotons Numero di fotoni da simulare.
 * @param useTIR   Se kTRUE, abilita la perdita di fotoni per mancata Riflessione Totale Interna.
 */
void ScintillatorSim(Int_t nPhotons, Bool_t useTIR = kFALSE) {

    // ==========================================
    // 1. DEFINIZIONE GEOMETRIA E FISICA
    // ==========================================
    
    // Dimensioni dello scintillatore (cm)
    Double_t L = 40.0; // Lunghezza (Asse X)
    Double_t W = 20.0; // Larghezza (Asse Y)
    Double_t H = 1.0;  // Spessore (Asse Z)
    
    // Raggio del fotomoltiplicatore (PMT) posizionato sulla faccia X+
    Double_t pmtRadius = 2.5; 
    
    // Proprietà ottiche
    Double_t n_ref = 1.58;            // Indice di rifrazione scintillatore (BC-408 circa)
    Double_t n_air = 1.0;             // Indice aria
    Double_t c_vac = 29.979;          // Velocità luce nel vuoto (cm/ns)
    Double_t v_scint = c_vac / n_ref; // Velocità luce nel mezzo
    Double_t attLength = 380.0;       // Lunghezza di attenuazione (cm)

    // --- Calcolo parametri per la Riflessione Totale Interna (TIR) ---
    // Legge di Snell: n1 * sin(theta_1) = n2 * sin(theta_2)
    // Angolo critico: sin(theta_c) = n2 / n1
    Double_t theta_crit = TMath::ASin(n_air / n_ref);
    Double_t cos_theta_crit = TMath::Cos(theta_crit);

    std::cout << "--- Configurazione ---" << std::endl;
    std::cout << "Simulazione: " << nPhotons << " fotoni." << std::endl;
    std::cout << "Modello TIR: " << (useTIR ? "ATTIVO" : "DISATTIVATO") << std::endl;
    std::cout << "Angolo Critico: " << theta_crit * TMath::RadToDeg() << " deg" << std::endl;
    std::cout << "Soglia Coseno:  " << cos_theta_crit << std::endl;

    Int_t maxBounces = 1000; // Limite massimo rimbalzi per evitare loop infiniti
    
    // ==========================================
    // 2. SETUP ISTOGRAMMI E CANVAS
    // ==========================================
    
    TCanvas *c1 = new TCanvas("c1", "Simulazione Scintillatore", 800, 600);
    
    // Istogramma temporale: Range 0-20ns
    TH1F *hTime = new TH1F("hTime", "Distribuzione Tempi di Arrivo;Tempo (ns);Conteggi", 200, 0, 20);
    
    // Generatore numeri casuali
    TRandom3 *rnd = new TRandom3(0); 

    // ==========================================
    // 3. LOOP MONTE CARLO PRINCIPALE
    // ==========================================
    
    for (Int_t i = 0; i < nPhotons; i++) {
        
        // --- A. Generazione Posizione Iniziale (Isotropa nel volume) ---
        Double_t x = rnd->Uniform(0, L);
        Double_t y = rnd->Uniform(-W/2.0, W/2.0);
        Double_t z = rnd->Uniform(-H/2.0, H/2.0);
        TVector3 pos(x, y, z);

        // --- B. Generazione Direzione (Isotropa 4pi) ---
        Double_t costheta = rnd->Uniform(-1.0, 1.0);
        Double_t theta = TMath::ACos(costheta);
        Double_t phi = rnd->Uniform(0, 2 * TMath::Pi());
        
        TVector3 dir;
        dir.SetMagThetaPhi(1.0, theta, phi); // Vettore unitario direzione

        Double_t time = 0.0;
        Bool_t detected = kFALSE;
        
        // --- C. Ray Tracing (Rimbalzi) ---
        for (Int_t b = 0; b < maxBounces; b++) {
            
            // 1. Calcolo distanza verso tutte le 6 pareti
            Double_t d_min = 1e9;
            Int_t wall_hit = -1; // ID parete colpita: 0=X+, 1=X-, 2=Y+, 3=Y-, 4=Z+, 5=Z-

            // Lambda function per trovare l'intersezione più vicina
            auto checkWall = [&](Double_t distToWall, Double_t dirComp, Int_t id) {
                if (dirComp != 0) { // Evita divisione per zero
                    Double_t d = distToWall / dirComp;
                    // Se d > 0 (è davanti a noi) e d < d_min (è la più vicina trovata finora)
                    if (d > 1e-9 && d < d_min) { d_min = d; wall_hit = id; }
                }
            };

            // Controllo intersezioni assi X, Y, Z
            if (dir.X() > 0) checkWall(L - pos.X(), dir.X(), 0);
            else             checkWall(0 - pos.X(), dir.X(), 1);
            
            if (dir.Y() > 0) checkWall(W/2.0 - pos.Y(), dir.Y(), 2);
            else             checkWall(-W/2.0 - pos.Y(), dir.Y(), 3);
            
            if (dir.Z() > 0) checkWall(H/2.0 - pos.Z(), dir.Z(), 4);
            else             checkWall(-H/2.0 - pos.Z(), dir.Z(), 5);

            // 2. Avanzamento del fotone alla parete
            pos = pos + (dir * d_min);
            time += d_min / v_scint;

            // 3. Simulazione Assorbimento (Attenuazione esponenziale)
            // Probabilità di sopravvivere = exp(-distanza / lunghezza_attenuazione)
            if (rnd->Uniform(0, 1) > TMath::Exp(-d_min / attLength)) {
                break; // Fotone morto (assorbito dal materiale)
            }

            // 4. Controllo Rivelazione (Solo sulla parete X+ -> wall_hit == 0)
            if (wall_hit == 0) { 
                // Distanza dal centro della faccia (0,0 locale su YZ)
                Double_t distFromCenter = TMath::Sqrt(pos.Y()*pos.Y() + pos.Z()*pos.Z());
                if (distFromCenter <= pmtRadius) {
                    detected = kTRUE;
                    break; // Fotone rivelato, usciamo dal loop rimbalzi
                }
            }

            // 5. Gestione Riflessione Totale Interna (TIR) - SE ATTIVA
            if (useTIR) {
                Double_t cosIncidence = 0.0;
                
                // Calcoliamo il coseno dell'angolo rispetto alla normale della parete colpita.
                // Essendo un parallelepipedo allineato agli assi, è semplicemente
                // il valore assoluto della componente del vettore direzione.
                if (wall_hit <= 1)      cosIncidence = TMath::Abs(dir.X()); // Pareti X
                else if (wall_hit <= 3) cosIncidence = TMath::Abs(dir.Y()); // Pareti Y
                else                    cosIncidence = TMath::Abs(dir.Z()); // Pareti Z

                // CONDIZIONE DI FUGA:
                // Se theta_incidenza < theta_critico, il fotone scappa.
                // Matematicamente equivale a: cos(theta_incidenza) > cos(theta_critico)
                // (perché il coseno decresce tra 0 e 90 gradi).
                if (cosIncidence > cos_theta_crit) {
                    break; // Fotone perso (rifratto fuori dallo scintillatore)
                }
            }

            // 6. Calcolo Riflessione Speculare
            // Invertiamo solo la componente della velocità perpendicolare alla parete
            if (wall_hit <= 1)      dir.SetX(-dir.X());
            else if (wall_hit <= 3) dir.SetY(-dir.Y());
            else                    dir.SetZ(-dir.Z());
        }

        // Se il fotone è stato rivelato, riempiamo l'istogramma
        if (detected) {
            hTime->Fill(time);
        }
    }

    // ==========================================
    // 4. VISUALIZZAZIONE E STATISTICHE
    // ==========================================
    
    // Disattiviamo il box statistiche standard di ROOT
    hTime->SetStats(0);

    hTime->SetFillColor(kAzure-4);
    hTime->SetLineColor(kBlack);
    hTime->Draw();

    // --- Calcolo Mediana ---
    Double_t quantiles[1];
    Double_t probs[1] = {0.5}; // 0.5 indica il 50% (Mediana)
    // Se l'istogramma è vuoto (es. TIR molto aggressiva), evitiamo crash
    Double_t median = 0.0;
    if (hTime->GetEntries() > 0) {
        hTime->GetQuantiles(1, quantiles, probs);
        median = quantiles[0];
    }

    // --- Creazione Box Statistiche Personalizzato ---
    // Coordinate (0-1) relative alla finestra del canvas
    TPaveText *pt = new TPaveText(0.65, 0.65, 0.89, 0.89, "NDC"); 
    
    pt->SetBorderSize(1);
    pt->SetFillColor(kWhite); 
    pt->SetTextAlign(12); // Allinea a sx
    pt->SetTextFont(42);  // Font standard

    // Aggiunta righe di testo formattate
    pt->AddText(Form("Config: %s", useTIR ? "TIR ON" : "TIR OFF")); // Info utile
    pt->AddText(Form("Entries = %.0f", hTime->GetEntries()));
    pt->AddText(Form("Mean    = %.4f ns", hTime->GetMean()));
    pt->AddText(Form("Std Dev = %.4f ns", hTime->GetStdDev()));
    pt->AddText(Form("Median  = %.4f ns", median));

    pt->Draw();

    c1->Modified();
    c1->Update();
}
