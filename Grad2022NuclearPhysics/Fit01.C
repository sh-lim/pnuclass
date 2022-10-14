#include "Style.h"

#include "TMinuit.h"

void Fit01(){

	//TVirtualFitter::SetDefaultFitter("Minuit");

	gStyle->SetOptStat(0);

	gRandom = new TRandom3(1);

	TF1 *fpol1 = new TF1("fpo1","pol1",0,10);
	fpol1->SetParameters(2,1);

	TGraphErrors *g1 = new TGraphErrors;

	float data_yy[10] = {0.};

	for (int ii=0; ii<10; ii++){

		float xx = ii + 0.5;
		float yy = fpol1->Eval(xx) + gRandom->Gaus(0,1); 

		data_yy[ii] = yy;

		g1->SetPoint(ii, xx, yy);
		g1->SetPointError(ii, 0, 1);

	}//ii

	TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0,10,15);
	SetHistoStyle("x","y","",24,20);

	g1->SetMarkerStyle(24);
	g1->SetMarkerSize(0.8);
	g1->SetMarkerColor(1);
	g1->SetLineColor(1);
	g1->Draw("p");

	TF1 *f1 = new TF1("f1","pol1",0,10);
	g1->Fit(f1,"R0E");
	f1->Draw("same");

	gMinuit->mnmatu(1);

	gMinuit->SetErrorDef(4);
	TGraph *gcor2sig = (TGraph*)gMinuit->Contour(80,0,1);

	gMinuit->SetErrorDef(1);
	TGraph *gcor1sig = (TGraph*)gMinuit->Contour(80,0,1);

	TH1D *h2sig = new TH1D("h2sig","",100,0,10);
	TH1D *h1sig = new TH1D("h1sig","",100,0,10);

	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h2sig);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h1sig, 0.67);

	TLegend *leg = new TLegend(0.2,0.9-4*0.06,0.65,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.045);
	leg->AddEntry(f1,"y=p_{0} + p_{1}x","L");
	leg->AddEntry("",Form("p_{0}=%4.3f#pm%4.3f",f1->GetParameter(0),f1->GetParError(0)),"");
	leg->AddEntry("",Form("p_{1}=%4.3f#pm%4.3f",f1->GetParameter(1),f1->GetParError(1)),"");
	leg->AddEntry("",Form("#chi^{2}/NDF=%4.3f/%d",f1->GetChisquare(),f1->GetNDF()),"");
	leg->Draw();

	TH2D *hchi2 = new TH2D("hchi2","",250,0,5,250,0,2);

	for (int ii=0; ii<250; ii++){
		float p0 = hchi2->GetXaxis()->GetBinCenter(ii+1);
		for (int jj=0; jj<250; jj++){
			float p1 = hchi2->GetYaxis()->GetBinCenter(jj+1);

			float chi2 = 0.0;
			for (int kk=0; kk<10; kk++){
				float xx = kk + 0.5;
				float _yy = p0 + p1*xx;

				chi2 += (data_yy[kk] - _yy)*(data_yy[kk] - _yy); 
			}

			hchi2->SetBinContent(ii+1, jj+1, chi2); 
		}//jj
	}//ii

	TCanvas *c2 = new TCanvas("c2","c2",1.2*500,500);
	SetPadStyle(1);

	htmp = (TH1D*)hchi2;
	SetHistoStyle("p_{0}","p_{1}","",24,20);

	hchi2->Draw("colz");

	TCanvas *c3 = new TCanvas("c3","c3",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0,5,2);
	SetHistoStyle("p_{0}","p_{1}","",24,20);

	gcor2sig->SetFillColor(42);
	gcor2sig->SetLineColor(42);
	gcor2sig->Draw("f");

	gcor1sig->SetFillColor(38);
	gcor1sig->SetLineColor(38);
	gcor1sig->Draw("f");

	{
		TLegend *leg = new TLegend(0.25,0.9-0.06*2,0.65,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry(gcor2sig,"#Delta#chi^{2}<2","F");
		leg->AddEntry(gcor1sig,"#Delta#chi^{2}<1","F");
		leg->Draw();
	}

	TCanvas *c4 = new TCanvas("c4","c4",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0,10,15);
	SetHistoStyle("x","y","",24,20);

	g1->Draw("p");

	h2sig->SetFillColorAlpha(42,0.3);
	h2sig->SetLineColorAlpha(42,0.3);
	h2sig->Draw("e3 same");

	h1sig->SetFillColorAlpha(38,0.3);
	h1sig->SetLineColorAlpha(38,0.3);
	h1sig->Draw("e3 same");

	f1->Draw("same");

	{
		TLegend *leg = new TLegend(0.25,0.9-0.06*2,0.65,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry(h2sig,"95% Confidence Interval","F");
		leg->AddEntry(h1sig,"67% Confidence Interval","F");
		leg->Draw();
	}

}
