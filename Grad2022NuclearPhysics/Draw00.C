#include "Style.h"

void Draw00(){

	gStyle->SetOptStat(0);

	const int ncent = 3;
	const int centbin_lo[ncent] = {0, 0, 60};
	const int centbin_hi[ncent] = {100, 10, 100};

	const int nMarker[ncent] = {1, 2, 4};

	TFile *infile = new TFile("outfile_Read00.root","read");

	TH1D *hcent = (TH1D*)infile->Get("hcent");
	TH2D *hcent_pt = (TH2D*)infile->Get("hcent_pt");

	TH1D *hpt[ncent];
	TH1D *hratio[ncent];

	for (int ii=0; ii<ncent; ii++){

		int binmin = hcent_pt->GetXaxis()->FindBin(centbin_lo[ii]+0.001);
		int binmax = hcent_pt->GetXaxis()->FindBin(centbin_hi[ii]-0.001);
		
		hpt[ii] = (TH1D*)hcent_pt->ProjectionY(Form("hpt_cent%d",ii),binmin,binmax); 
		hpt[ii]->Sumw2();

		float dpt = hpt[ii]->GetBinWidth(1);
		float deta = 1.0;
		float nevent = hcent->Integral(binmin, binmax);

		hpt[ii]->Scale(1./dpt/deta/nevent);
		hpt[ii]->SetMarkerStyle(24);
		hpt[ii]->SetMarkerColor(nMarker[ii]);
		hpt[ii]->SetLineColor(nMarker[ii]);
	}//

	for (int ii=0; ii<ncent; ii++){
		hratio[ii] = (TH1D*)hpt[ii]->Clone(Form("hratio_%d",ii));
		hratio[ii]->Divide(hpt[0]);
	}

	TCanvas *c1 = new TCanvas("c1","c1",1.2*2*500,500);
	c1->Divide(2,1);

	c1->cd(1);
	SetPadStyle();
	gPad->SetLogy();

	htmp = (TH1D*)gPad->DrawFrame(0,1e-4,5,200);
	SetHistoStyle("p_{T} (GeV/c)","#frac{1}{N_{event}} #frac{d^{2}N_{ch}}{dp_{T}d#eta} (GeV/c)^{-1}","",22,20);
	htmp->GetYaxis()->SetTitleOffset(1.3);

	hpt[0]->Draw("p same");
	hpt[1]->Draw("p same");
	hpt[2]->Draw("p same");

	{
		TLegend *leg = new TLegend(0.6,0.7,0.93,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("","#it{pp} 13 TeV","h");
		leg->AddEntry("","|#eta|<0.5","h");
		leg->AddEntry(hpt[0],"0-100%","P");
		leg->AddEntry(hpt[1],"0-10%","P");
		leg->AddEntry(hpt[2],"60-100%","P");
		leg->Draw();
	}

	c1->cd(2);
	SetPadStyle();
	//gPad->SetLogy();

	htmp = (TH1D*)gPad->DrawFrame(0,0,5,5);
	SetHistoStyle("p_{T} (GeV/c)","Ratio","",22,20);
	htmp->GetYaxis()->SetTitleOffset(1.3);

	//hratio[0]->Draw("p same");
	hratio[1]->Draw("p same");
	hratio[2]->Draw("p same");

	{
		TLegend *leg = new TLegend(0.6,0.5,0.93,0.6);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry(hpt[1],"0-10%/0-100%","P");
		leg->AddEntry(hpt[2],"60-100%/0-100%","P");
		leg->Draw();
	}

}
