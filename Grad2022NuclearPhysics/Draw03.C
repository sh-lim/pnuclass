#include "Style.h"

void Draw03(){

	gStyle->SetOptStat(0);

	const int npt = 7;
	const double ptbin[npt+1] = {0, 1, 2, 4, 6, 8, 12, 24};

	TFile *infile = new TFile("outfile_runhist_pp13TeV_v2_try100.root","read");

	TH2D *hmass_pt[2][2];

	TH1D *hmass[2][2][npt];
	TH1D *hmass_sub[2][npt];

	for (int ii=0; ii<2; ii++){

		hmass_pt[0][ii] = (TH2D*)infile->Get(Form("hmass_pt_cut0_chg%d",ii));
		hmass_pt[1][ii] = (TH2D*)infile->Get(Form("hmass_pt_cut1_chg%d",ii));

		for (int ipt=0; ipt<npt; ipt++){
			hmass[0][ii][ipt] = (TH1D*)hmass_pt[0][ii]->ProjectionX(Form("hmass_cut0_%d_%d",ii,ipt),ipt+1,ipt+1);
			hmass[1][ii][ipt] = (TH1D*)hmass_pt[1][ii]->ProjectionX(Form("hmass_cut1_%d_%d",ii,ipt),ipt+1,ipt+1);

			hmass[0][ii][ipt]->Sumw2();
			hmass[1][ii][ipt]->Sumw2();

			hmass[0][ii][ipt]->SetMarkerStyle(24);
			hmass[0][ii][ipt]->SetMarkerColor(1+ii);
			hmass[0][ii][ipt]->SetLineColor(1+ii);

			hmass[1][ii][ipt]->SetMarkerStyle(24);
			hmass[1][ii][ipt]->SetMarkerColor(1+ii);
			hmass[1][ii][ipt]->SetLineColor(1+ii);
		}//ipt
	}//ii

	for (int jj=0; jj<2; jj++){
		for (int ipt=0; ipt<npt; ipt++){

			hmass_sub[jj][ipt] = (TH1D*)hmass[jj][0][ipt]->Clone(Form("hmass_sub_%d_%d",jj,ipt));
			hmass_sub[jj][ipt]->Add(hmass[jj][1][ipt], -1);

		}//ipt
	}//jj

	TCanvas *c0[2];
	TCanvas *c1[2];

	for (int ii=0; ii<2; ii++){

		c0[ii] = new TCanvas(Form("c0_%d",ii),Form("c0_%d",ii),1.2*3*400,2*400);
		c0[ii]->Divide(3,2);

		c1[ii] = new TCanvas(Form("c1_%d",ii),Form("c1_%d",ii),1.2*3*400,2*400);
		c1[ii]->Divide(3,2);

		for (int ipt=1; ipt<npt; ipt++){

			c0[ii]->cd(ipt);
			SetPadStyle();

			htmp = (TH1D*)hmass[ii][0][ipt];
			SetHistoStyle("Mass (GeV/c^{2})","");

			htmp->SetAxisRange(1.7+0.001,2.2-0.001);
			htmp->SetLineColor(1);
			hmass[ii][0][ipt]->Draw("p e");
			hmass[ii][1][ipt]->Draw("p e same");

			c1[ii]->cd(ipt);
			SetPadStyle();

			htmp = (TH1D*)hmass_sub[ii][ipt];
			SetHistoStyle("Mass (GeV/c^{2})","");

			htmp->SetAxisRange(1.7+0.001,2.2-0.001);
			htmp->SetLineColor(1);
			hmass_sub[ii][ipt]->Draw("p e");

			TF1 *f1 = new TF1("f1","gaus(0) + pol2(3)",1.7,2.2);
			f1->SetParameter(0, hmass_sub[ii][ipt]->GetMaximum());
			f1->SetParameter(1, 1.87);
			f1->SetParameter(2, 0.02);

			hmass_sub[ii][ipt]->Fit(f1,"R0");
			f1->Draw("same");
		}

	}



}
