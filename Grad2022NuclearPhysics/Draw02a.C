#include "Style.h"

void Draw02a(){

	gStyle->SetOptStat(0);
	gRandom = new TRandom3(1);

	const int nset = 7;
	const int nevent[nset] = {100, 200, 500, 1000, 2000, 5000, 10000};

	TFile *infile = new TFile("outfile_Read02.root","read");

	TH1D *h1[nset];

	for (int iset=0; iset<nset; iset++){
		h1[iset] = (TH1D*)infile->Get(Form("h1_%d",nevent[iset]));
	}


	float nch_int[nset] = {0};

	TCanvas *c1 = new TCanvas("c1","c1",1.2*3*500,2*500);
	c1->Divide(3,2);

	for (int iset=0; iset<6; iset++){
		c1->cd(iset+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,100,1.2*h1[iset]->GetMaximum());
		SetHistoStyle("N_{ch}^{|#eta|<1}","N_{event}","",24,20);
		htmp->GetYaxis()->SetTitleOffset(2.5);
		htmp->GetXaxis()->SetTitleOffset(2.5);

		h1[iset]->SetLineColor(1);
		h1[iset]->Draw("same");

		for (int ii=0; ii<h1[iset]->GetNbinsX(); ii++){
			float xx = h1[iset]->GetBinCenter(ii+1);
			float yy = h1[iset]->GetBinContent(ii+1);

			nch_int[iset] += xx*yy;
		}

		TLegend *leg = new TLegend(0.55,0.8,0.95,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N_{event}=%d",nevent[iset]),"h");
		leg->AddEntry("",Form("Mean=%4.2f",h1[iset]->GetMean()),"h");
		leg->AddEntry("",Form("StdDev=%4.2f",h1[iset]->GetStdDev()),"h");
		leg->Draw();
	}

	//bootstrapping
	
	TH1D *h2[nset];

	for (int iset=0; iset<nset; iset++){
		h2[iset] = new TH1D(Form("h2_%d",iset),"",100,35,45);

		for (int ii=0; ii<1000; ii++){

			float boot_nch_int = 0.0;

			for (int jj=0; jj<nevent[iset]; jj++){

				boot_nch_int += h1[iset]->GetRandom();

			}

			h2[iset]->Fill(boot_nch_int/nevent[iset]);

		}//ii
	}//iset

	TCanvas *c2 = new TCanvas("c2","c2",1.2*3*500,2*500);
	c2->Divide(3,2);

	for (int iset=0; iset<6; iset++){
		c2->cd(iset+1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(35,0,45,1.2*h2[iset]->GetMaximum());
		SetHistoStyle("Bootstrap mean","","",24,20);
		htmp->GetYaxis()->SetTitleOffset(2.5);
		htmp->GetXaxis()->SetTitleOffset(2.5);

		h2[iset]->SetLineColor(1);
		h2[iset]->Draw("same");

		TLegend *leg = new TLegend(0.20,0.95-0.05*6,0.5,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N_{event}=%d",nevent[iset]),"h");
		leg->AddEntry("",Form("Mean=%4.2f",h2[iset]->GetMean()),"h");
		leg->AddEntry("",Form("StdDev=%4.2f",h2[iset]->GetStdDev()),"h");
		leg->AddEntry("","","h");
		leg->AddEntry("",Form("#LTN#GT=%d/%d=%4.2f",int(nch_int[iset]),nevent[iset],nch_int[iset]/nevent[iset]),"h");
		leg->AddEntry("",Form("#delta#LTN#GT=#sqrt{%d}/%d=%4.2f",int(nch_int[iset]),nevent[iset],sqrt(nch_int[iset])/nevent[iset]),"h");
		leg->Draw();

	}

	TProfile *hprof1 = (TProfile*)infile->Get("hprof_cent_ntrk1");
	TProfile *hprof2 = (TProfile*)infile->Get("hprof_cent_ntrk2");
	TProfile *hprof3 = (TProfile*)infile->Get("hprof_cent_ntrk3");

	TGraphErrors *gprof1 = new TGraphErrors;
	TGraphErrors *gprof2 = new TGraphErrors;
	TGraphErrors *gprof3 = new TGraphErrors;

	for (int ii=0; ii<hprof1->GetNbinsX(); ii++){
		float xx = hprof1->GetBinCenter(ii+1);
		float yy = hprof1->GetBinContent(ii+1);

		float yy_err1 = hprof1->GetBinError(ii+1);
		float yy_err2 = hprof2->GetBinError(ii+1);
		float yy_err3 = hprof3->GetBinError(ii+1);

		gprof1->SetPoint(ii+1, xx, yy);
		gprof2->SetPoint(ii+1, xx, yy);
		gprof3->SetPoint(ii+1, xx, yy);

		gprof1->SetPointError(ii+1, 0, yy_err1);
		gprof2->SetPointError(ii+1, 0, yy_err2);
		gprof3->SetPointError(ii+1, 0, yy_err3);
	}

	TCanvas *c3 = new TCanvas("c3","c3",1.2*3*500,500);
	c3->Divide(3,1);

	c3->cd(1);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0,40,50);
	SetHistoStyle("V0M percentile","N_{ch}^{|#eta|<1}","",24,20);

	gprof1->SetMarkerStyle(24);
	gprof1->SetMarkerSize(0.2);
	gprof1->SetLineColor(1);
	gprof1->SetMarkerColor(1);
	gprof1->Draw("p");

	{
		TLegend *leg = new TLegend(0.50,0.95-0.05*5,0.95,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("Option=#doublequote #doublequote"),"h");
		leg->AddEntry("","0-1% centrality","h");
		leg->AddEntry("",Form("N_{event}=%d",int(hprof1->GetBinEntries(1))),"h");
		leg->AddEntry("",Form("BinContent=%4.2f",hprof1->GetBinContent(1)),"h");
		leg->AddEntry("",Form("BinError=%4.2f",hprof1->GetBinError(1)),"h");
		leg->Draw();
	}

	c3->cd(2);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0,40,50);
	SetHistoStyle("V0M percentile","N_{ch}^{|#eta|<1}","",24,20);

	gprof2->SetMarkerStyle(24);
	gprof2->SetMarkerSize(0.2);
	gprof2->SetLineColor(1);
	gprof2->SetMarkerColor(1);
	gprof2->Draw("p");

	{
		TLegend *leg = new TLegend(0.50,0.95-0.05*5,0.95,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("Option=#doublequotes#doublequote"),"h");
		leg->AddEntry("","0-1% centrality","h");
		leg->AddEntry("",Form("N_{event}=%d",int(hprof2->GetBinEntries(1))),"h");
		leg->AddEntry("",Form("BinContent=%4.2f",hprof2->GetBinContent(1)),"h");
		leg->AddEntry("",Form("BinError=%4.2f",hprof2->GetBinError(1)),"h");
		leg->Draw();
	}

	c3->cd(3);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0,40,50);
	SetHistoStyle("V0M percentile","N_{ch}^{|#eta|<1}","",24,20);

	gprof3->SetMarkerStyle(24);
	gprof3->SetMarkerSize(0.2);
	gprof3->SetLineColor(1);
	gprof3->SetMarkerColor(1);
	gprof3->Draw("p");

	{
		TLegend *leg = new TLegend(0.50,0.95-0.05*5,0.95,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("Option=#doublequoteg#doublequote"),"h");
		leg->AddEntry("","0-1% centrality","h");
		leg->AddEntry("",Form("N_{event}=%d",int(hprof3->GetBinEntries(1))),"h");
		leg->AddEntry("",Form("BinContent=%4.2f",hprof3->GetBinContent(1)),"h");
		leg->AddEntry("",Form("BinError=%4.2f",hprof3->GetBinError(1)),"h");
		leg->Draw();
	}

	return;


}
