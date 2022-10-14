#include "Style.h"

void DrawScore00(){
	
	gStyle->SetOptStat(0);

	ifstream fdata;
	fdata.open("score002.txt");

	float score;

	TH1D *hscore = new TH1D("hscore","",24,-5,115);

	while ( fdata >> score ){
		hscore->Fill(score/2);
	}

	TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)hscore;
	SetHistoStyle("Homework score","","",22,20);

	hscore->SetLineWidth(2);
	hscore->SetLineColor(1);
	hscore->Draw("same");

	TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.045);
	leg->AddEntry("","Average score, March 19","");
	leg->Draw();


}
