#include "Style.h"

void Draw02(){

	gStyle->SetOptStat(0);
	gRandom = new TRandom3(1);

	int nevent = 100;

	TFile *infile = new TFile("outfile_Read02.root","read");

	TH1D *h1 = (TH1D*)infile->Get(Form("h1_%d",nevent));

	TCanvas *c1 = new TCanvas("c1","c1",1.2*3*500,500);
	c1->Divide(3,1);

	c1->cd(1);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0,100,1.3*h1->GetMaximum());
	SetHistoStyle("N_{ch}^{|#eta|<1}","N_{event}","",24,20);

	h1->SetLineColor(1);
	h1->Draw("same");

	float nch_int = 0;
	for (int ii=0; ii<h1->GetNbinsX(); ii++){

		float xx = h1->GetBinCenter(ii+1);
		float yy = h1->GetBinContent(ii+1);

		nch_int += xx*yy;
	}

	c1->cd(2);
	SetPadStyle();

	TLegend *leg = new TLegend(0.25,0.25,0.85,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(20);
	leg->SetTextFont(43);
	leg->AddEntry("",Form("Mean=%4.2f",h1->GetMean()),"h");
	leg->AddEntry("",Form("StdDev=%4.2f",h1->GetStdDev()),"h");
	leg->AddEntry("","","");
	leg->AddEntry("",Form("N_{ch}^{%d evt}=%d",nevent,int(nch_int)),"h");
	leg->AddEntry("",Form("#sqrt{N_{ch}^{%d evt}}=%4.2f",nevent,sqrt(nch_int)),"h");
	leg->AddEntry("",Form("N_{ch}^{%d evt}/%d=%4.2f",nevent,nevent,nch_int/nevent),"h");
	leg->AddEntry("",Form("#sqrt{N_{ch}^{%d evt}}/%d=%4.2f",nevent,nevent,sqrt(nch_int)/nevent),"h");
	leg->AddEntry("",Form("#LTN_{ch}#GT in 0-1%% V0M =40.14#pm0.63"),"h");
	leg->Draw();

	//return;

	//bootstrapping
	
	TH1D *h2 = new TH1D("h2","",100,35,45);

	for (int ii=0; ii<5000; ii++){

		float boot_nch_int = 0.0;

		for (int jj=0; jj<nevent; jj++){

			boot_nch_int += h1->GetRandom();

		}

		h2->Fill(boot_nch_int/nevent);

	}

	c1->cd(3);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(35,0,45,1.3*h2->GetMaximum());
	SetHistoStyle("Bootstrapped #LTN_{ch}#GT in 0-1%% V0M","Counts","",24,20);

	h2->SetLineColor(1);
	h2->Draw("same");

	{
		TLegend *leg = new TLegend(0.23,0.75,0.85,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(20);
		leg->SetTextFont(43);
		leg->AddEntry("",Form("Mean=%4.2f",h2->GetMean()),"h");
		leg->AddEntry("",Form("StdDev=%4.2f",h2->GetStdDev()),"h");
		leg->Draw();
	}


	return;


}
