#include "PhotonFlux.h"
#include "DipoleCrossSection.h"

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"

#include <cmath>
#include <vector>

using namespace std;


void plotOnCanvas(vector<TGraph*> graphs, TString xAxis, TString yAxis,  TString outName,bool setLogY = true){
    if(graphs.size() == 0) return;

    TCanvas *c1 = new TCanvas("c1","Photon Flux vs Rapidity",800,600);
    auto legend = new TLegend(0.6,0.7,0.88,0.88);
    legend->SetBorderSize(0);
    
    int i = 2;
    for( auto g : graphs ){
        g->SetLineWidth(2);
        g->SetLineColor(i);
        i++;
        legend->AddEntry(g, g->GetName(), "l");
    }
    
    graphs[0]->GetYaxis()->SetTitle(yAxis);
    graphs[0]->GetXaxis()->SetTitle(xAxis);
    graphs[0]->Draw("AL");

    for(unsigned int i = 1; i < graphs.size(); ++i){
        graphs[i]->Draw("same L");
    }

    legend->Draw("same");

    if(setLogY)  c1->SetLogy();
    
    c1->Update();
    c1->SaveAs(outName + ".pdf");
}

void plotRapidity(){
    // initialize PhotonFlux
    PhotonFlux* photonFlux = new PhotonFlux();

    int nPoints = 1000;
    double yMin = -4.0, yMax = 4.0;
    TGraph *grDZ = new TGraph();
    grDZ->SetName("Drees-Zeppenfeld");
    TGraph *grDZM = new TGraph();
    grDZM->SetName("Modified Drees-Zeppenfeld");
    TGraph *grK = new TGraph();
    grK->SetName("Kniehl");
    TGraph *grJ = new TGraph();
    grJ->SetName("Jackson");
    for(int i = 0; i < nPoints; ++i){
        double y = yMin + (yMax - yMin) * i / (nPoints - 1);
        grDZ->SetPoint(i, y, photonFlux->getDreesZeppenfeldFromRapidity(y));
        grDZM->SetPoint(i, y, photonFlux->getModifiedDZFromRapidity(y));
        grK->SetPoint(i, y, photonFlux->getKniehlFromRapidity(y));
        grJ->SetPoint(i, y, photonFlux->getJacksonFromRapidity(y));
    }

    plotOnCanvas({grDZ, grDZM, grK, grJ},"y","dN/dy", "PhotonFlux_vs_Rapidity", true);

}


void plotFractionOfMomentum(){
    // initialize PhotonFlux
    PhotonFlux* photonFlux = new PhotonFlux();
    TGraph *grDZ = new TGraph();
    grDZ->SetName("Drees-Zeppenfeld");
    TGraph *grDZM = new TGraph();
    grDZM->SetName("Modified Drees-Zeppenfeld");
    TGraph *grK = new TGraph();
    grK->SetName("Kniehl");
    TGraph *grJ = new TGraph();
    grJ->SetName("Jackson");
    int nPoints = 1000;
    double xMin = 0.001, xMax = 0.2;
    for(int i = 0; i < nPoints; ++i){
        double x = xMin + (xMax - xMin) * i / (nPoints - 1);
        grDZ->SetPoint(i, x, x*photonFlux->getDreesZeppenfeld(x));
        grDZM->SetPoint(i, x, x*photonFlux->getModifiedDZ(x));
        grK->SetPoint(i, x, x*photonFlux->getKniehl(x));
        grJ->SetPoint(i, x, x*photonFlux->getJackson(x));
    }

    plotOnCanvas({grDZ, grDZM, grK, grJ}, "x", "x * f(x)", "PhotonFlux_vs_FractionOfMomentum", true);

}

void plotPhotonEnergy(){
    // initialize PhotonFlux
    PhotonFlux* photonFlux = new PhotonFlux();

    int nPoints = 1000;
    double omegaMin = 0.01, omegaMax = 250;  //GeV
    TGraph *grDZ = new TGraph();
    TGraph *grDZM = new TGraph();
    TGraph *grK = new TGraph();
    TGraph *grJ = new TGraph();
    grDZ->SetName("Drees-Zeppenfeld");
    grDZM->SetName("Modified Drees-Zeppenfeld");
    grK->SetName("Kniehl");
    grJ->SetName("Jackson");
    for(int i = 0; i < nPoints; ++i){
        double omega = omegaMin + (omegaMax - omegaMin) * i / (nPoints - 1);
        grDZ->SetPoint(i, omega, photonFlux->getDreesZeppenfeldFromPhotonEnergy(omega));
        grDZM->SetPoint(i, omega, photonFlux->getModifiedDZFromPhotonEnergy(omega));
        grK->SetPoint(i, omega, photonFlux->getKniehlFromPhotonEnergy(omega));
        grJ->SetPoint(i, omega, photonFlux->getJacksonFromPhotonEnergy(omega));
    }

    plotOnCanvas({grDZ, grDZM, grK, grJ},"#omega" , "dN/d#omega", "PhotonFlux_vs_PhotonEnergy", true);

}


void plotDipoleCrossSection(double x, double rMin, double rMax, bool setLogY) {

    DipoleCrossSection *mDipoleCS = new DipoleCrossSection();

    if(!mDipoleCS){
        cout << "ERROR: DipoleCrossSection not initialized" << endl;
        return;
    }

    string gridName = "../grids_gbw_dglap/dipcs_dglap_2.dat";
    mDipoleCS->loadGrid(gridName);

    TGraph *grGBWOld = new TGraph();
    grGBWOld->SetName("GBW Old");
    TGraph *grGBWNew = new TGraph();
    grGBWNew->SetName("GBW New");
    TGraph *grKST = new TGraph();
    grKST->SetName("KST");
    TGraph *grGBWNewNew = new TGraph();
    grGBWNewNew->SetName("GBW New New");
    TGraph *grGBWDglap = new TGraph();
    grGBWDglap->SetName("GBW Dglap");
    int nPoints = 1000;
    for(int i = 0; i < nPoints; ++i){
        double rVal = rMin + (rMax - rMin) * i / (nPoints - 1);
        grGBWOld->SetPoint(i, rVal, mDipoleCS->getGBWOld(x, rVal));
        grGBWNew->SetPoint(i, rVal, mDipoleCS->getGBWNew(x, rVal));
        grKST->SetPoint(i, rVal, mDipoleCS->getKST(x, rVal));
        grGBWNewNew->SetPoint(i, rVal, mDipoleCS->getGBWNewNew(x, rVal));
        grGBWDglap->SetPoint(i, rVal, mDipoleCS->getGBWDglap(x, rVal));
    }

    plotOnCanvas({grGBWOld, grGBWNew, grKST, grGBWNewNew, grGBWDglap}, "r", "dN/d#sigma" , "DipoleCrossSection_x_" + to_string(x), setLogY);

    delete mDipoleCS;
}
