#include "Style.h"

#include "TMinuit.h"

void Fit02a(){

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
	h1->Fit(ftn, "R0L");

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

		TH1D *h2sig = new TH1D("h2sig_LS1","",100,-5,5);
		TH1D *h1sig = new TH1D("h1sig_LS1","",100,-5,5);

		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h2sig);
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h1sig, 0.67);

		h2sig->SetFillColorAlpha(42,0.3);
		h2sig->SetLineColorAlpha(42,0.3);
		h2sig->Draw("e3 same");

		h1sig->SetFillColorAlpha(38,0.3);
		h1sig->SetLineColorAlpha(38,0.3);
		h1sig->Draw("e3 same");

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*5);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=10"),"h");
		leg->AddEntry("",Form("Least-squares method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->AddEntry("",Form("N=%4.2f#pm%4.2f",f1->Integral(-5,5)/0.1,f1->IntegralError(-5,5)/0.1),"h");
		leg->Draw();
	}

	//return;

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

		TH1D *h2sig = new TH1D("h2sig_LS2","",100,-5,5);
		TH1D *h1sig = new TH1D("h1sig_LS2","",100,-5,5);

		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h2sig);
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h1sig, 0.67);

		h2sig->SetFillColorAlpha(42,0.3);
		h2sig->SetLineColorAlpha(42,0.3);
		h2sig->Draw("e3 same");

		h1sig->SetFillColorAlpha(38,0.3);
		h1sig->SetLineColorAlpha(38,0.3);
		h1sig->Draw("e3 same");

		float norm_err1 = 0, norm_err2 = 0;
		for (int ii=0; ii<100; ii++){
			norm_err1 += (h1sig->GetBinError(ii+1));
			norm_err2 += (h2sig->GetBinError(ii+1));
		}

		cout << h1sig->Integral() << " +/- " << norm_err1 << endl;
		cout << h2sig->Integral() << " +/- " << norm_err2 << endl;
		cout << f1->Integral(-5,5)/h1->GetBinWidth(1) << " +/- " << f1->IntegralError(-5,5)/h1->GetBinWidth(1) << endl;

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*5);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=100"),"h");
		leg->AddEntry("",Form("Least-squares method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->AddEntry("",Form("N=%4.2f#pm%4.2f",f1->Integral(-5,5)/0.1,f1->IntegralError(-5,5)/0.1),"h");
		leg->Draw();

		TCanvas *c4_2 = new TCanvas("c4_2","c4_2",1.2*3*500,500);
		c4_2->Divide(3,1);

		c4_2->cd(1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,-2,6,2);
		SetHistoStyle("p_{0}","p_{1}","",24,20);

		gMinuit->SetErrorDef(1);
		TGraph *gcor1sig = (TGraph*)gMinuit->Contour(80,0,1);

		gMinuit->SetErrorDef(4);
		TGraph *gcor2sig = (TGraph*)gMinuit->Contour(80,0,1);

		gcor2sig->SetFillColor(42);
		gcor2sig->SetLineColor(42);
		gcor2sig->Draw("f");

		gcor1sig->SetFillColor(38);
		gcor1sig->SetLineColor(38);
		gcor1sig->Draw("f");

		c4_2->cd(2);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,6,5);
		SetHistoStyle("p_{0}","p_{2}","",24,20);

		gMinuit->SetErrorDef(1);
		gcor1sig = (TGraph*)gMinuit->Contour(80,0,2);

		gMinuit->SetErrorDef(4);
		gcor2sig = (TGraph*)gMinuit->Contour(80,0,2);

		gcor2sig->SetFillColor(42);
		gcor2sig->SetLineColor(42);
		gcor2sig->Draw("f");

		gcor1sig->SetFillColor(38);
		gcor1sig->SetLineColor(38);
		gcor1sig->Draw("f");

		c4_2->cd(3);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(-2,0,2,5);
		SetHistoStyle("p_{1}","p_{2}","",24,20);

		gMinuit->SetErrorDef(1);
		gcor1sig = (TGraph*)gMinuit->Contour(80,1,2);

		gMinuit->SetErrorDef(4);
		gcor2sig = (TGraph*)gMinuit->Contour(80,1,2);

		gcor2sig->SetFillColor(42);
		gcor2sig->SetLineColor(42);
		gcor2sig->Draw("f");

		gcor1sig->SetFillColor(38);
		gcor1sig->SetLineColor(38);
		gcor1sig->Draw("f");
	}

	//return;

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

		TH1D *h2sig = new TH1D("h2sig_LS3","",100,-5,5);
		TH1D *h1sig = new TH1D("h1sig_LS3","",100,-5,5);

		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h2sig);
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h1sig, 0.67);

		h2sig->SetFillColorAlpha(42,0.3);
		h2sig->SetLineColorAlpha(42,0.3);
		h2sig->Draw("e3 same");

		h1sig->SetFillColorAlpha(38,0.3);
		h1sig->SetLineColorAlpha(38,0.3);
		h1sig->Draw("e3 same");

		float norm_err1 = 0, norm_err2 = 0;
		for (int ii=0; ii<100; ii++){
			norm_err1 += (h1sig->GetBinError(ii+1));
			norm_err2 += (h2sig->GetBinError(ii+1));
		}

		cout << h1sig->Integral() << " +/- " << norm_err1 << endl;
		cout << h2sig->Integral() << " +/- " << norm_err2 << endl;
		cout << f1->Integral(-5,5)/h1->GetBinWidth(1) << " +/- " << f1->IntegralError(-5,5)/h1->GetBinWidth(1) << endl;

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*5);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=1000"),"h");
		leg->AddEntry("",Form("Least-squares method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->AddEntry("",Form("N=%4.2f#pm%4.2f",f1->Integral(-5,5)/0.1,f1->IntegralError(-5,5)/0.1),"h");
		leg->Draw();

		TCanvas *c4_3 = new TCanvas("c4_3","c4_3",1.2*3*500,500);
		c4_3->Divide(3,1);

		c4_3->cd(1);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,-2,60,2);
		SetHistoStyle("p_{0}","p_{1}","",24,20);

		gMinuit->SetErrorDef(1);
		TGraph *gcor1sig = (TGraph*)gMinuit->Contour(80,0,1);

		gMinuit->SetErrorDef(4);
		TGraph *gcor2sig = (TGraph*)gMinuit->Contour(80,0,1);

		gcor2sig->SetFillColor(42);
		gcor2sig->SetLineColor(42);
		gcor2sig->Draw("f");

		gcor1sig->SetFillColor(38);
		gcor1sig->SetLineColor(38);
		gcor1sig->Draw("f");

		c4_3->cd(2);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,60,5);
		SetHistoStyle("p_{0}","p_{2}","",24,20);

		gMinuit->SetErrorDef(1);
		gcor1sig = (TGraph*)gMinuit->Contour(80,0,2);

		gMinuit->SetErrorDef(4);
		gcor2sig = (TGraph*)gMinuit->Contour(80,0,2);

		gcor2sig->SetFillColor(42);
		gcor2sig->SetLineColor(42);
		gcor2sig->Draw("f");

		gcor1sig->SetFillColor(38);
		gcor1sig->SetLineColor(38);
		gcor1sig->Draw("f");

		c4_3->cd(3);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(-2,0,2,5);
		SetHistoStyle("p_{1}","p_{2}","",24,20);

		gMinuit->SetErrorDef(1);
		gcor1sig = (TGraph*)gMinuit->Contour(80,1,2);

		gMinuit->SetErrorDef(4);
		gcor2sig = (TGraph*)gMinuit->Contour(80,1,2);

		gcor2sig->SetFillColor(42);
		gcor2sig->SetLineColor(42);
		gcor2sig->Draw("f");

		gcor1sig->SetFillColor(38);
		gcor1sig->SetLineColor(38);
		gcor1sig->Draw("f");
	}

	return;

	c4->cd(4);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-5,0,5,1.5*h4->GetMaximum());
	SetHistoStyle("x","Counts / 0.1","",24,20);
	htmp->GetYaxis()->SetTitleOffset(2);

	h4->SetLineColor(1);
	h4->Draw("same");

	{
		TF1 *f1 = new TF1("f1","gaus",-5,5);
		TFitResultPtr r = h4->Fit(f1, "R0Q");
		f1->Draw("same");

		TH1D *h2sig = new TH1D("h2sig_LS4","",100,-5,5);
		TH1D *h1sig = new TH1D("h1sig_LS4","",100,-5,5);

		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h2sig);
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h1sig, 0.67);

		h2sig->SetFillColorAlpha(42,0.3);
		h2sig->SetLineColorAlpha(42,0.3);
		h2sig->Draw("e3 same");

		h1sig->SetFillColorAlpha(38,0.3);
		h1sig->SetLineColorAlpha(38,0.3);
		h1sig->Draw("e3 same");

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*5);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=10000"),"h");
		leg->AddEntry("",Form("Least-squares method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->AddEntry("",Form("N=%4.2f#pm%4.2f",f1->Integral(-5,5)/0.1,f1->IntegralError(-5,5)/0.1),"h");
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
		h3->Fit(f1, "R0QL");
		f1->Draw("same");

		TH1D *h2sig = new TH1D("h2sig_LL1","",100,-5,5);
		TH1D *h1sig = new TH1D("h1sig_LL1","",100,-5,5);

		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h2sig);
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h1sig, 0.67);

		h2sig->SetFillColorAlpha(42,0.3);
		h2sig->SetLineColorAlpha(42,0.3);
		h2sig->Draw("e3 same");

		h1sig->SetFillColorAlpha(38,0.3);
		h1sig->SetLineColorAlpha(38,0.3);
		h1sig->Draw("e3 same");

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*5);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=10"),"h");
		leg->AddEntry("",Form("Log-likelihood method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->AddEntry("",Form("N=%4.2f#pm%4.2f",f1->Integral(-5,5)/0.1,f1->IntegralError(-5,5)/0.1),"h");
		leg->Draw();
	}

	//return;

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

		TH1D *h2sig = new TH1D("h2sig_LL2","",100,-5,5);
		TH1D *h1sig = new TH1D("h1sig_LL2","",100,-5,5);

		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h2sig);
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h1sig, 0.67);

		h2sig->SetFillColorAlpha(42,0.3);
		h2sig->SetLineColorAlpha(42,0.3);
		h2sig->Draw("e3 same");

		h1sig->SetFillColorAlpha(38,0.3);
		h1sig->SetLineColorAlpha(38,0.3);
		h1sig->Draw("e3 same");

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*5);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=100"),"h");
		leg->AddEntry("",Form("Log-likelihood method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->AddEntry("",Form("N=%4.2f#pm%4.2f",f1->Integral(-5,5)/0.1,f1->IntegralError(-5,5)/0.1),"h");
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

		TH1D *h2sig = new TH1D("h2sig_LL3","",100,-5,5);
		TH1D *h1sig = new TH1D("h1sig_LL3","",100,-5,5);

		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h2sig);
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h1sig, 0.67);

		h2sig->SetFillColorAlpha(42,0.3);
		h2sig->SetLineColorAlpha(42,0.3);
		h2sig->Draw("e3 same");

		h1sig->SetFillColorAlpha(38,0.3);
		h1sig->SetLineColorAlpha(38,0.3);
		h1sig->Draw("e3 same");

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*5);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=1000"),"h");
		leg->AddEntry("",Form("Log-likelihood method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->AddEntry("",Form("N=%4.2f#pm%4.2f",f1->Integral(-5,5)/0.1,f1->IntegralError(-5,5)/0.1),"h");
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

		TH1D *h2sig = new TH1D("h2sig_LL4","",100,-5,5);
		TH1D *h1sig = new TH1D("h1sig_LL4","",100,-5,5);

		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h2sig);
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h1sig, 0.67);

		h2sig->SetFillColorAlpha(42,0.3);
		h2sig->SetLineColorAlpha(42,0.3);
		h2sig->Draw("e3 same");

		h1sig->SetFillColorAlpha(38,0.3);
		h1sig->SetLineColorAlpha(38,0.3);
		h1sig->Draw("e3 same");

		float nval = f1->Integral(-5,5) / h4->GetBinWidth(1);
		float nvalerr = f1->IntegralError(-5,5) / h4->GetBinWidth(1);

		TLegend *leg = new TLegend(0.2,0.93,0.5,0.93-0.055*5);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("",Form("N=10000"),"h");
		leg->AddEntry("",Form("Log-likelihood method"),"h");
		leg->AddEntry("",Form("#mu=%4.2f#pm%4.2f",f1->GetParameter(1),f1->GetParError(1)),"h");
		leg->AddEntry("",Form("#sigma=%4.2f#pm%4.2f",f1->GetParameter(2),f1->GetParError(2)),"h");
		leg->AddEntry("",Form("N=%4.2f#pm%4.2f",f1->Integral(-5,5)/0.1,f1->IntegralError(-5,5)/0.1),"h");
		leg->Draw();
	}


	return;

}
