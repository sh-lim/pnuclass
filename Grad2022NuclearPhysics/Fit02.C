#include "Style.h"

void Fit02(){

	const int ndata = 100;
	const float const_pi = TMath::Pi();

	gStyle->SetOptStat(0);

	gRandom = new TRandom3(1);

	TH1D *h1 = new TH1D("h1","Gaus",100,-5,5);
	h1->Sumw2();

	TH1D *h2 = new TH1D("h2","Gaus",100,-5,5);
	h2->Sumw2();

	TH1D *h3 = new TH1D("h3","Gaus",100,-5,5);
	h3->Sumw2();

	TH1D *h4 = new TH1D("h4","Gaus",100,-5,5);
	h4->Sumw2();

	for (int ii=0; ii<ndata; ii++){
		h1->Fill(gRandom->Gaus(0,1.25));
	}

	for (int ii=0; ii<10*ndata; ii++){
		h2->Fill(gRandom->Gaus(0,1.25));
	}

	for (int ii=0; ii<ndata/10; ii++){
		h3->Fill(gRandom->Gaus(0,1.25));
	}

	for (int ii=0; ii<100*ndata; ii++){
		h4->Fill(gRandom->Gaus(0,1.25));
	}

	float sum_sq = 0.0;
	for (int ibin=0; ibin<h1->GetNbinsX(); ibin++){
		float xx = h1->GetBinCenter(ibin+1);
		float yy = h1->GetBinContent(ibin+1);

		sum_sq += yy*xx*xx;
	}

	cout << "sum sq: " << sum_sq << endl;

	TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-5,0,5,10);
	SetHistoStyle("x","Counts / 0.1","",24,20);

	h1->SetLineColor(1);
	h1->Draw("same");

	TLegend *leg = new TLegend(0.2,0.7,0.6,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.045);
	leg->AddEntry("",Form("N=%d",ndata),"h");
	leg->AddEntry("","#mu=0","h");
	leg->AddEntry("","#sigma=1.25","h");
	leg->Draw();

	TF1 *fgaus = new TF1("fgaus","gaus",-5,5);
	h1->Fit(fgaus,"R0Q");
	//fgaus->Draw("same");

	TGraph *gLL = new TGraph; 
	gLL->SetMarkerStyle(24);
	gLL->SetMarkerSize(0.8);

	TGraph *gLL2 = new TGraph; 
	gLL2->SetMarkerStyle(24);
	gLL2->SetMarkerSize(0.8);

	for (int ipar=0; ipar<100; ipar++){

		float par0 = 0 + 5./100*(ipar+1);
		float logL = -ndata*log(sqrt(2*const_pi)) - ndata/2.0*log(par0*par0) - sum_sq/(2*par0*par0);

		//cout << par0 << " " << logL << endl;

		gLL->SetPoint(ipar, par0, logL);
	}

	for (int ipar=0; ipar<100; ipar++){

		float par0 = -2.5 + 5./100*(ipar+1);
		float tmp_dev = 0.0;
		for (int ibin=0; ibin<h1->GetNbinsX(); ibin++){
			float xx = h1->GetBinCenter(ibin+1);
			float yy = h1->GetBinContent(ibin+1);

			tmp_dev += yy*(xx-par0)*(xx-par0);
		}

		float logL = -ndata*log(sqrt(2*const_pi)) - ndata/2.0*log(1.25*1.25) - tmp_dev/(2*1.25*1.25);

		//cout << par0 << " " << logL << endl;

		gLL2->SetPoint(ipar, par0, logL);
	}

	TCanvas *c2 = new TCanvas("c2","c2",1.2*500,500);
	SetPadStyle(0);

	htmp = (TH1D*)gPad->DrawFrame(0,-700,2.5,0);
	SetHistoStyle("log-likelihood","","",24,20);

	gLL->Draw("P");

	TF1 *fLL = new TF1("fLL", "pol8", 0.5, 5);
	gLL->Fit(fLL, "R0Q");
	fLL->Draw("same");

	float fLL_max = fLL->GetMaximum();
	float fLL_maxx = fLL->GetMaximumX();
	float fLL_errx1 = fLL->GetX(fLL_max-0.5,0,fLL_maxx);
	float fLL_errx2 = fLL->GetX(fLL_max-0.5,fLL_maxx,2.5);

	{
		TLegend *leg = new TLegend(0.4,0.25,0.9,0.25+0.06*3);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("maximum LL=%4.2f",fLL_max),"h");
		leg->AddEntry("",Form("#sigma for the maximum LL=%4.2f",fLL_maxx),"h");
		leg->AddEntry("",Form("#sigma for the maximum LL-0.5=%4.2f, %4.2f",fLL_errx1,fLL_errx2),"h");
		leg->Draw();
	}

	TCanvas *c3 = new TCanvas("c3","c3",1.2*500,500);
	SetPadStyle(0);

	htmp = (TH1D*)gPad->DrawFrame(-2.5,-700,2.5,0);
	SetHistoStyle("log-likelihood","","",24,20);

	gLL2->Draw("P");

	TF1 *fLL2 = new TF1("fLL2", "pol4", -2.5, 2.5);
	gLL2->Fit(fLL2, "R0Q");
	fLL2->Draw("same");

	float fLL2_max = fLL2->GetMaximum();
	float fLL2_maxx = fLL2->GetMaximumX();
	float fLL2_errx1 = fLL2->GetX(fLL2_max-0.5,-2.5,fLL2_maxx);
	float fLL2_errx2 = fLL2->GetX(fLL2_max-0.5,fLL2_maxx,2.5);

	{
		TLegend *leg = new TLegend(0.4,0.25,0.9,0.25+0.06*3);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("maximum LL=%4.2f",fLL2_max),"h");
		leg->AddEntry("",Form("#mu for the maximum LL=%4.2f",fLL2_maxx),"h");
		leg->AddEntry("",Form("#mu for the maximum LL-0.5=%4.2f, %4.2f",fLL2_errx1,fLL2_errx2),"h");
		leg->Draw();
	}

	//return;


	TF1 *ftn = new TF1("ftn","gaus",-5,5);
	h1->Fit(ftn, "R0QL");

	c1->cd();
	ftn->Draw("same");

	{
		TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("","Fit result","h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",ftn->GetParameter(1),ftn->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",ftn->GetParameter(2),ftn->GetParError(2)),"h");
		leg->Draw();
	}

	TCanvas *c4 = new TCanvas("c4","c4",1.2*2*500,2*500);
	c4->Divide(2,2);

	c4->cd(1);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-5,0,5,4*h3->GetMaximum());
	SetHistoStyle("x","Counts / 0.1","",24,20);
	htmp->GetYaxis()->SetTitleOffset(2);

	h3->SetLineColor(1);
	h3->Draw("same");

	{
		TF1 *f1 = new TF1("f1","gaus",-5,5);
		h3->Fit(f1, "R0Q");
		f1->Draw("same");

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=10"),"h");
		leg->AddEntry("",Form("Least-squares method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->Draw();
	}

	c4->cd(2);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-5,0,5,2*h1->GetMaximum());
	SetHistoStyle("x","Counts / 0.1","",24,20);
	htmp->GetYaxis()->SetTitleOffset(2);

	h1->SetLineColor(1);
	h1->Draw("same");

	{
		TF1 *f1 = new TF1("f1","gaus",-5,5);
		h1->Fit(f1, "R0Q");
		f1->Draw("same");

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=100"),"h");
		leg->AddEntry("",Form("Least-squares method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->Draw();
	}

	c4->cd(3);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-5,0,5,1.5*h2->GetMaximum());
	SetHistoStyle("x","Counts / 0.1","",24,20);
	htmp->GetYaxis()->SetTitleOffset(2);

	h2->SetLineColor(1);
	h2->Draw("same");

	{
		TF1 *f1 = new TF1("f1","gaus",-5,5);
		h2->Fit(f1, "R0Q");
		f1->Draw("same");

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=1000"),"h");
		leg->AddEntry("",Form("Least-squares method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->Draw();
	}

	c4->cd(4);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-5,0,5,1.5*h4->GetMaximum());
	SetHistoStyle("x","Counts / 0.1","",24,20);
	htmp->GetYaxis()->SetTitleOffset(2);

	h4->SetLineColor(1);
	h4->Draw("same");

	{
		TF1 *f1 = new TF1("f1","gaus",-5,5);
		h4->Fit(f1, "R0Q");
		f1->Draw("same");

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=10000"),"h");
		leg->AddEntry("",Form("Least-squares method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->Draw();
	}

	TCanvas *c5 = new TCanvas("c5","c5",1.2*2*500,2*500);
	c5->Divide(2,2);

	c5->cd(1);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-5,0,5,4*h3->GetMaximum());
	SetHistoStyle("x","Counts / 0.1","",24,20);
	htmp->GetYaxis()->SetTitleOffset(2);

	h3->SetLineColor(1);
	h3->Draw("same");

	{
		TF1 *f1 = new TF1("f1","gaus",-5,5);
		h3->Fit(f1, "R0L");
		f1->Draw("same");

		cout << h3->Integral() << " " << f1->Integral(-5,5)/h3->GetBinWidth(1) << endl;

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=10"),"h");
		leg->AddEntry("",Form("Log-likelihood method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->Draw();
	}

	c5->cd(2);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-5,0,5,2*h1->GetMaximum());
	SetHistoStyle("x","Counts / 0.1","",24,20);
	htmp->GetYaxis()->SetTitleOffset(2);

	h1->SetLineColor(1);
	h1->Draw("same");

	{
		TF1 *f1 = new TF1("f1","gaus",-5,5);
		h1->Fit(f1, "R0QL");
		f1->Draw("same");

		cout << h1->Integral() << " " << f1->Integral(-5,5)/h3->GetBinWidth(1) << endl;

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=100"),"h");
		leg->AddEntry("",Form("Log-likelihood method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->Draw();
	}

	c5->cd(3);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-5,0,5,1.5*h2->GetMaximum());
	SetHistoStyle("x","Counts / 0.1","",24,20);
	htmp->GetYaxis()->SetTitleOffset(2);

	h2->SetLineColor(1);
	h2->Draw("same");

	{
		TF1 *f1 = new TF1("f1","gaus",-5,5);
		h2->Fit(f1, "R0QL");
		f1->Draw("same");

		cout << h2->Integral() << " " << f1->Integral(-5,5)/h3->GetBinWidth(1) << endl;

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=1000"),"h");
		leg->AddEntry("",Form("Log-likelihood method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->Draw();
	}

	c5->cd(4);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-5,0,5,1.5*h4->GetMaximum());
	SetHistoStyle("x","Counts / 0.1","",24,20);
	htmp->GetYaxis()->SetTitleOffset(2);

	h4->SetLineColor(1);
	h4->Draw("same");

	{
		TF1 *f1 = new TF1("f1","gaus",-5,5);
		h4->Fit(f1, "R0QL");
		f1->Draw("same");

		cout << h4->Integral() << " " << f1->Integral(-5,5)/h3->GetBinWidth(1) << endl;

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*4);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=10000"),"h");
		leg->AddEntry("",Form("Log-likelihood method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->Draw();
	}


	return;

}
