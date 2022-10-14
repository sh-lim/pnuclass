#include "Style.h"

void Draw03(){

	gStyle->SetOptStat(0);

	const int npt = 7;
	const double ptbin[npt+1] = {0, 1, 2, 4, 6, 8, 12, 24};

	TFile *infile = new TFile("outfile_runhist_pp13TeV_v2_try101.root","read");

	TH2D *hmass_pt[2][2];
	TH2D *hmass_pt_rot[2][2];
	TH2D *hmass_pt_mix[2][2];

	TH1D *hmass[2][2][npt];
	TH1D *hmass_rot[2][2][npt];
	TH1D *hmass_mix[2][2][npt];

	TH1D *hmass_sub[2][npt];

	for (int ii=0; ii<2; ii++){

		hmass_pt[0][ii] = (TH2D*)infile->Get(Form("hmass_pt_cut0_chg%d",ii));
		hmass_pt[1][ii] = (TH2D*)infile->Get(Form("hmass_pt_cut1_chg%d",ii));

		hmass_pt_rot[0][ii] = (TH2D*)infile->Get(Form("hmass_pt_rot_cut0_chg%d",ii));
		hmass_pt_rot[1][ii] = (TH2D*)infile->Get(Form("hmass_pt_rot_cut1_chg%d",ii));

		hmass_pt_mix[0][ii] = (TH2D*)infile->Get(Form("hmass_pt_mix_cut0_chg%d",ii));
		hmass_pt_mix[1][ii] = (TH2D*)infile->Get(Form("hmass_pt_mix_cut1_chg%d",ii));

		for (int ipt=0; ipt<npt; ipt++){
			hmass[0][ii][ipt] = (TH1D*)hmass_pt[0][ii]->ProjectionX(Form("hmass_cut0_%d_%d",ii,ipt),ipt+1,ipt+1);
			hmass[1][ii][ipt] = (TH1D*)hmass_pt[1][ii]->ProjectionX(Form("hmass_cut1_%d_%d",ii,ipt),ipt+1,ipt+1);

			hmass_rot[0][ii][ipt] = (TH1D*)hmass_pt_rot[0][ii]->ProjectionX(Form("hmass_rot_cut0_%d_%d",ii,ipt),ipt+1,ipt+1);
			hmass_rot[1][ii][ipt] = (TH1D*)hmass_pt_rot[1][ii]->ProjectionX(Form("hmass_rot_cut1_%d_%d",ii,ipt),ipt+1,ipt+1);

			hmass_mix[0][ii][ipt] = (TH1D*)hmass_pt_mix[0][ii]->ProjectionX(Form("hmass_mix_cut0_%d_%d",ii,ipt),ipt+1,ipt+1);
			hmass_mix[1][ii][ipt] = (TH1D*)hmass_pt_mix[1][ii]->ProjectionX(Form("hmass_mix_cut1_%d_%d",ii,ipt),ipt+1,ipt+1);

			hmass[0][ii][ipt]->Sumw2();
			hmass[1][ii][ipt]->Sumw2();

			hmass_rot[0][ii][ipt]->Sumw2();
			hmass_rot[1][ii][ipt]->Sumw2();

			hmass_mix[0][ii][ipt]->Sumw2();
			hmass_mix[1][ii][ipt]->Sumw2();

			hmass[0][ii][ipt]->SetMarkerStyle(24);
			hmass[0][ii][ipt]->SetMarkerColor(1+ii);
			hmass[0][ii][ipt]->SetLineColor(1+ii);
			hmass[1][ii][ipt]->SetMarkerStyle(24);
			hmass[1][ii][ipt]->SetMarkerColor(1+ii);
			hmass[1][ii][ipt]->SetLineColor(1+ii);

			hmass_rot[0][ii][ipt]->SetMarkerStyle(25);
			hmass_rot[0][ii][ipt]->SetMarkerColor(1+ii);
			hmass_rot[0][ii][ipt]->SetLineColor(1+ii);
			hmass_rot[1][ii][ipt]->SetMarkerStyle(25);
			hmass_rot[1][ii][ipt]->SetMarkerColor(1+ii);
			hmass_rot[1][ii][ipt]->SetLineColor(1+ii);

			hmass_mix[0][ii][ipt]->SetMarkerStyle(27);
			hmass_mix[0][ii][ipt]->SetMarkerColor(1+ii);
			hmass_mix[0][ii][ipt]->SetLineColor(1+ii);
			hmass_mix[1][ii][ipt]->SetMarkerStyle(27);
			hmass_mix[1][ii][ipt]->SetMarkerColor(1+ii);
			hmass_mix[1][ii][ipt]->SetLineColor(1+ii);
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
	TCanvas *c2[2];
	TCanvas *c3[2];

	for (int ii=0; ii<2; ii++){

		c0[ii] = new TCanvas(Form("c0_%d",ii),Form("c0_%d",ii),1.2*3*400,2*400);
		c0[ii]->Divide(3,2);

		c1[ii] = new TCanvas(Form("c1_%d",ii),Form("c1_%d",ii),1.2*3*400,2*400);
		c1[ii]->Divide(3,2);

		c2[ii] = new TCanvas(Form("c2_%d",ii),Form("c2_%d",ii),1.2*3*400,2*400);
		c2[ii]->Divide(3,2);

		c3[ii] = new TCanvas(Form("c3_%d",ii),Form("c3_%d",ii),1.2*3*400,2*400);
		c3[ii]->Divide(3,2);

		for (int ipt=1; ipt<npt; ipt++){

			c0[ii]->cd(ipt);
			SetPadStyle();
			gPad->SetTopMargin(0.05);

			htmp = (TH1D*)hmass[ii][0][ipt];
			SetHistoStyle("Mass (GeV/c^{2})","");
			htmp->GetXaxis()->SetTitleOffset(2.3);
			htmp->SetMinimum(0);
			htmp->SetMaximum(1.1*htmp->GetMaximum());

			htmp->SetAxisRange(1.7+0.001,2.2-0.001);
			htmp->SetLineColor(1);
			hmass[ii][0][ipt]->Draw("p e");
			hmass[ii][1][ipt]->Draw("p e same");

			{
				TLegend *leg;
				if ( ipt>3 ){
					leg = new TLegend(0.55,0.94-0.06*4,0.94,0.94);
				}else{
					leg = new TLegend(0.55,0.20,0.94,0.20+0.06*4);
				}
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.05);
				if ( ii==0 ){
					leg->AddEntry("",Form("w/o d_{0} cut"),"");
				}else{
					leg->AddEntry("",Form("w/ d_{0} cut"),"");
				}
				leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
				leg->AddEntry(hmass[ii][0][ipt],"Unlike-sign","P");
				leg->AddEntry(hmass[ii][1][ipt],"Like-sign","P");
				leg->Draw();
			}

			c1[ii]->cd(ipt);
			SetPadStyle();
			gPad->SetTopMargin(0.05);

			htmp = (TH1D*)hmass[ii][1][ipt];
			SetHistoStyle("Mass (GeV/c^{2})","");
			htmp->GetXaxis()->SetTitleOffset(2.3);

			htmp->SetAxisRange(1.7+0.001,2.2-0.001);
			htmp->SetMinimum(0);
			htmp->SetMaximum(1.1*htmp->GetMaximum());

			hmass[ii][1][ipt]->Draw("p e");
			hmass_mix[ii][0][ipt]->Draw("p e same");
			hmass_mix[ii][1][ipt]->Draw("p e same");
			hmass_rot[ii][0][ipt]->Draw("p e same");
			hmass_rot[ii][1][ipt]->Draw("p e same");

			if ( ipt<4 ){
				TLegend *leg = new TLegend(0.35,0.15,0.94,0.15+0.05*7);
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.045);
				if ( ii==0 ){
					leg->AddEntry("",Form("w/o d_{0} cut"),"");
				}else{
					leg->AddEntry("",Form("w/ d_{0} cut"),"");
				}
				leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
				leg->AddEntry(hmass[ii][1][ipt],"Same event, Like-sign","P");
				leg->AddEntry(hmass_mix[ii][0][ipt],"Mixed event, Unlike-sign","P");
				leg->AddEntry(hmass_mix[ii][1][ipt],"Mixed event, Like-sign","P");
				leg->AddEntry(hmass_rot[ii][0][ipt],"Track rotation, Unlike-sign","P");
				leg->AddEntry(hmass_rot[ii][1][ipt],"Track rotation, Like-sign","P");
				leg->Draw();
			}

			c2[ii]->cd(ipt);
			SetPadStyle();
			gPad->SetTopMargin(0.05);

			htmp = (TH1D*)hmass_sub[ii][ipt];
			SetHistoStyle("Mass (GeV/c^{2})","");
			htmp->GetXaxis()->SetTitleOffset(2.3);
			htmp->SetMaximum(1.1*htmp->GetMaximum());

			htmp->SetAxisRange(1.7+0.001,2.2-0.001);
			htmp->SetLineColor(1);
			hmass_sub[ii][ipt]->Draw("p e");

			TF1 *f1 = new TF1("f1","gaus(0) + pol1(3)",1.7,2.2);
			f1->SetParameter(0, hmass_sub[ii][ipt]->GetMaximum());
			f1->SetParameter(1, 1.87);
			f1->SetParameter(2, 0.02);

			TF1 *f2 = new TF1("f2","gaus(0) + pol2(3)",1.7,2.2);
			f2->SetParameter(0, hmass_sub[ii][ipt]->GetMaximum());
			f2->SetParameter(1, 1.87);
			f2->SetParameter(2, 0.02);

			hmass_sub[ii][ipt]->Fit(f1,"R0");
			hmass_sub[ii][ipt]->Fit(f2,"R0");
			f1->Draw("same");
			f2->SetLineColor(4);
			f2->Draw("same");

			{
				TLegend *leg = new TLegend(0.5,0.90-0.05*4,0.94,0.90);
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.045);
				if ( ii==0 ){
					leg->AddEntry("",Form("w/o d_{0} cut"),"");
				}else{
					leg->AddEntry("",Form("w/ d_{0} cut"),"");
				}
				leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
				leg->AddEntry(f1,"Gaus sig. + pol1 bkg.","L");
				leg->AddEntry(f2,"Gaus sig. + pol2 bkg.","L");
				leg->Draw();
			}

			c3[ii]->cd(ipt);
			SetPadStyle();
			gPad->SetTopMargin(0.05);

			htmp = (TH1D*)hmass[ii][0][ipt];
			SetHistoStyle("Mass (GeV/c^{2})","");
			htmp->GetXaxis()->SetTitleOffset(2.3);
			htmp->SetMaximum(1.1*htmp->GetMaximum());

			htmp->SetAxisRange(1.7+0.001,2.2-0.001);
			htmp->SetLineColor(1);
			hmass[ii][0][ipt]->Draw("p e");

			TF1 *f3 = new TF1("f3","gaus(0) + pol2(3)",1.7,2.2);
			f3->SetParameter(0, hmass[ii][0][ipt]->GetMaximum());
			f3->SetParameter(1, 1.87);
			f3->SetParameter(2, 0.02);

			TF1 *f4 = new TF1("f4","gaus(0) + pol3(3)",1.7,2.2);
			f4->SetParameter(0, hmass[ii][0][ipt]->GetMaximum());
			f4->SetParameter(1, 1.87);
			f4->SetParameter(2, 0.02);

			hmass[ii][0][ipt]->Fit(f3, "R0");
			hmass[ii][0][ipt]->Fit(f4, "R0");
			f3->Draw("same");
			//f4->SetLineColor(4);
			//f4->Draw("same");

			{
				TLegend *leg = new TLegend(0.5,0.90-0.05*4,0.94,0.90);
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.045);
				if ( ii==0 ){
					leg->AddEntry("",Form("w/o d_{0} cut"),"");
				}else{
					leg->AddEntry("",Form("w/ d_{0} cut"),"");
				}
				leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
				leg->AddEntry(f3,"Gaus sig. + pol2 bkg.","L");
				leg->AddEntry("","","");
				leg->Draw();
			}
		}

	}



}
