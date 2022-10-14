#include "Style.h"

void Fit00(){

	gRandom = new TRandom3(1);

	TF1 *fpol0 = new TF1("fpo0","pol0",0,10);
	fpol0->SetParameter(0, 5);

	TGraphErrors *g1 = new TGraphErrors;

	float data_yy[10] = {0.};

	for (int ii=0; ii<10; ii++){

		float xx = ii + 0.5;
		float yy = fpol0->Eval(xx) + gRandom->Gaus(0,1); 

		data_yy[ii] = yy;

		g1->SetPoint(ii, xx, yy);
		g1->SetPointError(ii, 0, 1);

	}//ii

	TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0,10,10);
	SetHistoStyle("x","y","",24,20);

	g1->SetMarkerStyle(24);
	g1->SetMarkerSize(0.8);
	g1->SetMarkerColor(1);
	g1->SetLineColor(1);
	g1->Draw("p");

	for (int ii=0; ii<11; ii++){
		TF1 *f1 = new TF1("f1","pol0",0,10);
		f1->SetLineWidth(1);
		f1->SetLineStyle(2);
		f1->SetParameter(0, 2.5+0.5*ii);
		f1->Draw("same");
	}

	TGraphErrors *gchi2 = new TGraphErrors;

	for (int ii=0; ii<110; ii++){
		float _yy = 2.5 + 0.05*ii;

		float chi2 = 0.0;
		for (int jj=0; jj<10; jj++){

			chi2 += (data_yy[jj] - _yy)*(data_yy[jj] - _yy);

		}//jj

		gchi2->SetPoint(ii, _yy, chi2);
	}//ii

	TCanvas *c2 = new TCanvas("c2","c2",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(2.5,0,7.5,100);
	SetHistoStyle("p_{0}","#chi^{2}","",24,20);

	gchi2->SetMarkerStyle(24);
	gchi2->SetMarkerSize(0.8);
	gchi2->SetMarkerColor(1);
	gchi2->Draw("p");

	TF1 *fchi2 = new TF1("fchi2","pol2",2.5,7.5);
	gchi2->Fit(fchi2,"R0Q");
	fchi2->SetLineWidth(1);
	fchi2->Draw("same");

	float fchi2_min = fchi2->GetMinimum();
	float fchi2_minx = fchi2->GetMinimumX();
	float fchi2_errx1 = fchi2->GetX(fchi2_min+1,0,fchi2_minx);
	float fchi2_errx2 = fchi2->GetX(fchi2_min+1,fchi2_minx,7.5);

	{
		TLegend *leg = new TLegend(0.25,0.9-0.06*3,0.65,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("minimum #chi^{2}=%4.2f",fchi2->GetMinimum()),"h");
		leg->AddEntry("",Form("p_{0} for the minimum #chi^{2}=%4.2f",fchi2->GetMinimumX()),"h");
		leg->AddEntry("",Form("p_{0} for the minimum #chi^{2}+1=%4.2f, %4.2f",fchi2_errx1,fchi2_errx2),"h");
		leg->Draw();
	}

	TCanvas *c3 = new TCanvas("c3","c3",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0,10,10);
	SetHistoStyle("x","y","",24,20);

	g1->SetMarkerStyle(24);
	g1->SetMarkerSize(0.8);
	g1->SetMarkerColor(1);
	g1->SetLineColor(1);
	g1->Draw("p");

	TF1 *ffit = new TF1("ffit","pol0",0,10);
	g1->Fit(ffit,"R0Q");

	ffit->SetLineStyle(2);
	ffit->Draw("same");

	TLegend *leg = new TLegend(0.2,0.93-3*0.06,0.65,0.93);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.045);
	leg->AddEntry(ffit,"y=p_{0}","L");
	leg->AddEntry("",Form("p_{0}=%4.3f#pm%4.3f",ffit->GetParameter(0),ffit->GetParError(0)),"");
	leg->AddEntry("",Form("#chi^{2}/NDF=%4.3f/%d",ffit->GetChisquare(),ffit->GetNDF()),"");
	leg->Draw();

	return;

	TF1 *f1 = new TF1("f1","pol1",0,10);
	g1->Fit(f1,"R0");
	//f1->Draw("same");

}
