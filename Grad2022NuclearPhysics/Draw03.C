#include "Style.h"

void Draw03(){

	gStyle->SetOptStat(0);

	const int npt = 7;
	const double ptbin[npt+1] = {0, 1, 2, 4, 6, 8, 12, 24};

	TFile *infile = new TFile("outfile_Read03.root","read");

	TH2D *hmass_pt_cut0[2];
	TH2D *hmass_pt_cut1[2];

	TH1D *hmass_cut0[2][npt];
	TH1D *hmass_cut1[2][npt];

	for (int ii=0; ii<2; ii++){

		hmass_pt_cut0[ii] = (TH2D*)infile->Get(Form("hmass_pt_cut0_chg%d",ii));
		hmass_pt_cut1[ii] = (TH2D*)infile->Get(Form("hmass_pt_cut1_chg%d",ii));

		for (int ipt=0; ipt<npt; ipt++){
			hmass_cut0[ii][ipt] = (TH1D*)hmass_pt_cut0[ii]->ProjectionX(Form("hmass_cut0_%d_%d",ii,ipt),ipt+1,ipt+1);
			hmass_cut1[ii][ipt] = (TH1D*)hmass_pt_cut1[ii]->ProjectionX(Form("hmass_cut1_%d_%d",ii,ipt),ipt+1,ipt+1);
		}//ipt

	}//ii

	TCanvas *c0[2];

	for (int ii=0; ii<2; ii++){

		c0[ii] = new TCanvas(Form("c0_%d",ii),Form("c0_%d",ii),1.2*3*400,2*400);
		c0[ii]->Divide(3,2);

		for (int ipt=1; ipt<npt; ipt++){

			c0[ii]->cd(ipt);
			SetPadStyle();

			htmp = (TH1D*)hmass_cut0[ii][ipt];
			SetHistoStyle("Mass (GeV/c^{2})","");

			htmp->SetAxisRange(1.7+0.001,2.2-0.001);
			htmp->Sumw2();
			htmp->SetLineColor(1);
			hmass_cut0[ii][ipt]->Draw("");

		}

	}



}
