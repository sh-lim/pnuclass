#include "Style.h"

void DrawMC00(){

	gStyle->SetOptStat(0);

	const int npt = 7;
	const double ptbin[npt+1] = {0, 1, 2, 4, 6, 8, 12, 24};

	TFile *infile = new TFile("outfile_ReadMC00.root","read");

	TH2D *hmc_cut0 = (TH2D*)infile->Get("hmc_cut0");
	TH2D *hmc_cut1 = (TH2D*)infile->Get("hmc_cut1");

	TH1D *hmc_eta_cut0 = (TH1D*)hmc_cut0->ProjectionX("hmc_eta_cut0");
	TH1D *hmc_eta_cut1 = (TH1D*)hmc_cut1->ProjectionX("hmc_eta_cut1");
	hmc_eta_cut0->SetLineColor(1);
	hmc_eta_cut0->SetLineWidth(1);

	hmc_eta_cut1->SetLineColor(4);
	hmc_eta_cut1->SetLineWidth(1);

	int etamin = hmc_eta_cut0->FindBin(-0.5+0.001);
	int etamax = hmc_eta_cut0->FindBin(+0.5-0.001);

	TH1D *hmc_pt_cut0 = (TH1D*)hmc_cut0->ProjectionY("hmc_pt_cut0",etamin,etamax);

	TCanvas *c1 = new TCanvas("c1","c1",1.2*3*500,500);
	c1->Divide(3,1);

	{
		c1->cd(1);
		SetPadStyle(1);
		gPad->SetLeftMargin(0.12);

		htmp = (TH2D*)hmc_cut0;
		hmc_cut0->SetAxisRange(0,10,"Y");
		SetHistoStyle("#eta","p_{T} (GeV/c)","",24,20);
		htmp->GetYaxis()->SetTitleOffset(1.15);

		hmc_cut0->Draw("colz");
	}

	{
		c1->cd(2);
		SetPadStyle(1);
		gPad->SetLeftMargin(0.12);

		htmp = (TH2D*)hmc_cut1;
		hmc_cut1->SetAxisRange(0,10,"Y");
		SetHistoStyle("#eta","p_{T} (GeV/c)","",24,20);
		htmp->GetYaxis()->SetTitleOffset(1.15);

		hmc_cut1->Draw("colz");
	}

	{
		c1->cd(3);
		SetPadStyle();

		htmp = (TH1D*)hmc_eta_cut0;
		SetHistoStyle("#eta","","",24,20);

		hmc_eta_cut0->Draw("");
		hmc_eta_cut1->Draw("same");
	}

	TCanvas *c100 = new TCanvas;
	gPad->SetLogy();
	hmc_pt_cut0->Draw();

	//return;

	//
	TH2D *hmass_pt[2][2];
	TH1D *hmass[2][2][npt];

	TH2D *hmass_pt_truth[2][2];
	TH1D *hmass_truth[2][2][npt];

	for (int ii=0; ii<2; ii++){
		hmass_pt[0][ii] = (TH2D*)infile->Get(Form("hmass_pt_cut0_chg%d",ii));
		hmass_pt[1][ii] = (TH2D*)infile->Get(Form("hmass_pt_cut1_chg%d",ii));

		hmass_pt_truth[0][ii] = (TH2D*)infile->Get(Form("hmass_pt_truth_cut0_chg%d",ii));
		hmass_pt_truth[1][ii] = (TH2D*)infile->Get(Form("hmass_pt_truth_cut1_chg%d",ii));

		for (int ipt=0; ipt<npt; ipt++){
			hmass[0][ii][ipt] = (TH1D*)hmass_pt[0][ii]->ProjectionX(Form("hmass_cut0_%d_%d",ii,ipt),ipt+1,ipt+1);
			hmass[1][ii][ipt] = (TH1D*)hmass_pt[1][ii]->ProjectionX(Form("hmass_cut1_%d_%d",ii,ipt),ipt+1,ipt+1);

			hmass_truth[0][ii][ipt] = (TH1D*)hmass_pt_truth[0][ii]->ProjectionX(Form("hmass_truth_cut0_%d_%d",ii,ipt),ipt+1,ipt+1);
			hmass_truth[1][ii][ipt] = (TH1D*)hmass_pt_truth[1][ii]->ProjectionX(Form("hmass_truth_cut1_%d_%d",ii,ipt),ipt+1,ipt+1);

			hmass[0][ii][ipt]->Sumw2();
			hmass[1][ii][ipt]->Sumw2();

			hmass_truth[0][ii][ipt]->Sumw2();
			hmass_truth[1][ii][ipt]->Sumw2();

			hmass[0][ii][ipt]->SetMarkerStyle(24);
			hmass[0][ii][ipt]->SetMarkerColor(1+ii);
			hmass[0][ii][ipt]->SetLineColor(1+ii);

			hmass[1][ii][ipt]->SetMarkerStyle(24);
			hmass[1][ii][ipt]->SetMarkerColor(1+ii);
			hmass[1][ii][ipt]->SetLineColor(1+ii);

			hmass_truth[0][ii][ipt]->SetMarkerStyle(24);
			hmass_truth[0][ii][ipt]->SetMarkerColor(4);
			hmass_truth[0][ii][ipt]->SetLineColor(4);

			hmass_truth[1][ii][ipt]->SetMarkerStyle(24);
			hmass_truth[1][ii][ipt]->SetMarkerColor(4);
			hmass_truth[1][ii][ipt]->SetLineColor(4);
		}//ipt

	}//ii

	TH1D *hN1[2];
	TH1D *hN2[2];
	
	TCanvas *c0[2];

	for (int ii=1; ii<2; ii++){

		c0[ii] = new TCanvas(Form("c0_%d",ii),Form("c0_%d",ii),1.2*4*400,2*400);
		c0[ii]->Divide(4,2);

		hN1[ii] = new TH1D(Form("hN1_%d",ii),"",npt,ptbin);
		hN2[ii] = new TH1D(Form("hN2_%d",ii),"",npt,ptbin);

		for (int ipt=0; ipt<npt; ipt++){

			{
				c0[ii]->cd(ipt+1);
				SetPadStyle();
				gPad->SetTopMargin(0.05);

				htmp = (TH1D*)hmass[ii][0][ipt];
				SetHistoStyle("Mass (GeV/c^{2})","");
				htmp->GetXaxis()->SetTitleOffset(2.3);
				htmp->SetMinimum(1);
				htmp->SetMaximum(1.1*htmp->GetMaximum());

				htmp->SetAxisRange(1.7+0.001,2.2-0.001);
				htmp->SetLineColor(1);
				hmass[ii][0][ipt]->Draw("p e");
				hmass_truth[ii][0][ipt]->Draw("p e same");

				TF1 *f1 = new TF1("f1","gaus(0) + pol1(3)",1.7,2.2);
				f1->SetLineColor(1);
				f1->SetNpx(500);
				f1->SetParameter(1, 1.87);
				f1->SetParameter(2, 0.02);
				hmass[ii][0][ipt]->Fit(f1,"R0");
				f1->Draw("same");

				float mean = f1->GetParameter(1);
				float sigma = fabs(f1->GetParameter(2));

				TF1 *fsig = new TF1("fsig","gaus",1.7,2.2);
				fsig->SetParameter(0, f1->GetParameter(0));
				fsig->SetParameter(1, f1->GetParameter(1));
				fsig->SetParameter(2, f1->GetParameter(2));
				fsig->SetLineColor(2);
				fsig->Draw("same");

				float N2 = fsig->Integral(mean-3*sigma, mean+3*sigma) / hmass[ii][0][ipt]->GetBinWidth(1);
				float N1 = hmass_truth[ii][0][ipt]->Integral(hmass_truth[ii][0][ipt]->FindBin(mean-3*sigma),hmass_truth[ii][0][ipt]->FindBin(mean+3*sigma));

				float Nden = hmc_pt_cut0->Integral(hmc_pt_cut0->FindBin(ptbin[ipt]+0.001),hmc_pt_cut0->FindBin(ptbin[ipt+1]-0.001));

				cout << N1 << " " << N2 << " " << Nden << endl;

				hN1[ii]->SetBinContent(ipt+1, N1/Nden);
				hN1[ii]->SetBinError(ipt+1, sqrt(N1)/Nden);
				hN2[ii]->SetBinContent(ipt+1, N2/Nden);
				hN2[ii]->SetBinError(ipt+1, sqrt(N2)/Nden);

				{
					TLegend *leg;
					leg = new TLegend(0.5,0.94-0.06*8,0.94,0.94);
					leg->SetFillStyle(0);
					leg->SetBorderSize(0);
					leg->SetTextSize(0.05);
					leg->AddEntry("","Unlike-sign","");
					if ( ii==0 ){
						leg->AddEntry("",Form("w/o d_{0} cut"),"");
					}else{
						leg->AddEntry("",Form("w/ d_{0} cut"),"");
					}
					leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
					leg->AddEntry(hmass[ii][0][ipt],"w/o Mother PID tag","P");
					leg->AddEntry("",Form("N_{fit}=%d",int(N2)),"");
					leg->AddEntry(hmass_truth[ii][0][ipt],"w/ Mother PID tag","P");
					leg->AddEntry("",Form("N_{hist}=%d",int(N1)),"");
					leg->AddEntry("",Form("#sigma for signal=%4.4f",fabs(fsig->GetParameter(2))),"");
					leg->Draw();
				}
			}

		}//ipt

	}//ii

	//TCanvas *c10 = new TCanvas("c10","c10",1.2*500,500);
	c0[1]->cd(8);
	SetPadStyle();

	{
		htmp = (TH1D*)gPad->DrawFrame(0,0,24,0.3);
		SetHistoStyle("p_{T} (GeV/c)","Acceptance and efficiency");
		htmp->GetYaxis()->SetTitleOffset(2.5);
		htmp->GetXaxis()->SetTitleOffset(2.3);

		hN1[1]->SetMarkerStyle(24);
		hN1[1]->SetMarkerColor(1);
		hN1[1]->SetLineColor(1);
		hN1[1]->Draw("p same");

		hN2[1]->SetMarkerStyle(25);
		hN2[1]->SetMarkerColor(2);
		hN2[1]->SetLineColor(2);
		hN2[1]->Draw("p same");

		TLegend *leg;
		leg = new TLegend(0.25,0.94-0.05*6,0.6,0.94);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("","D^{0}#rightarrow#piK","");
		leg->AddEntry("","|y|<0.5","");
		leg->AddEntry("","w/ d_{0} cut","");
		leg->AddEntry(hN1[1],"Fit function","P");
		leg->AddEntry(hN2[1],"Histogram","P");
		leg->Draw();
	}

	return;

	TFile *outfile = new TFile("outfile_acceff.root","recreate");
	hN1[1]->Write("hacceff_1");
	hN2[1]->Write("hacceff_2");
	outfile->Close();


}

