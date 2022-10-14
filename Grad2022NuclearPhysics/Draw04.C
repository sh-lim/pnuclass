#include "Style.h"

void Draw04(){

	gStyle->SetOptStat(0);

	const int npt = 7;
	const double ptbin[npt+1] = {0, 1, 2, 4, 6, 8, 12, 24};

	TFile *infile = new TFile("outfile_runhist_pp13TeV_v2_try104.root","read");

	TH2D *hmass_pt[2][2];
	TH2D *hmass_pt_rot[2][2];
	TH2D *hmass_pt_mix[2][2];

	TH1D *hmass[2][2][npt];
	TH1D *hmass_rot[2][2][npt];
	TH1D *hmass_mix[2][2][npt];

	TH1D *hmass_sub1[2][npt];
	TH1D *hmass_sub2[2][npt];

	for (int ii=0; ii<2; ii++){

		hmass_pt[0][ii] = (TH2D*)infile->Get(Form("hmass_pt_cut0_chg%d",ii));
		hmass_pt[1][ii] = (TH2D*)infile->Get(Form("hmass_pt_cut1_chg%d",ii));

		hmass_pt_rot[0][ii] = (TH2D*)infile->Get(Form("hmass_pt_rot_cut0_chg%d",ii));
		hmass_pt_rot[1][ii] = (TH2D*)infile->Get(Form("hmass_pt_rot_cut1_chg%d",ii));

		hmass_pt_mix[0][ii] = (TH2D*)infile->Get(Form("hmass_pt_mix_cut0_chg%d",ii));
		hmass_pt_mix[1][ii] = (TH2D*)infile->Get(Form("hmass_pt_mix_cut1_chg%d",ii));

		/*
		hmass_pt[0][ii]->RebinX(2);
		hmass_pt[1][ii]->RebinX(2);

		hmass_pt_rot[0][ii]->RebinX(2);
		hmass_pt_rot[1][ii]->RebinX(2);

		hmass_pt_mix[0][ii]->RebinX(2);
		hmass_pt_mix[1][ii]->RebinX(2);
		*/

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

			hmass_sub1[jj][ipt] = (TH1D*)hmass[jj][0][ipt]->Clone(Form("hmass_sub1_%d_%d",jj,ipt));

			hmass_sub1[jj][ipt]->Add(hmass[jj][1][ipt], -1);

			hmass_sub2[jj][ipt] = (TH1D*)hmass[jj][0][ipt]->Clone(Form("hmass_sub2_%d_%d",jj,ipt));

			TH1D *hbkg_mix = (TH1D*)hmass_mix[jj][0][ipt]->Clone(Form("hbkg_mix_%d_%d",jj,ipt));
			hbkg_mix->Add(hmass_mix[jj][1][ipt]);

			float norm1 = hmass[jj][0][ipt]->Integral(hmass[jj][0][ipt]->FindBin(1.6),hmass[jj][0][ipt]->FindBin(1.7));
			float norm2 = hbkg_mix->Integral(hbkg_mix->FindBin(1.6), hbkg_mix->FindBin(1.7));

			if ( ipt>1 ){
				norm1 = hmass[jj][0][ipt]->Integral(hmass[jj][0][ipt]->FindBin(2.2),hmass[jj][0][ipt]->FindBin(2.3));
				norm2 = hbkg_mix->Integral(hbkg_mix->FindBin(2.2), hbkg_mix->FindBin(2.3));
			}

			hbkg_mix->Scale(norm1/norm2);
			hmass_sub2[jj][ipt]->Add(hbkg_mix, -1);

		}//ipt
	}//jj

	TCanvas *c0[2];
	TCanvas *c1[2];
	TCanvas *c2[2];
	TCanvas *c3[2];
	TCanvas *c4[2];

	for (int ii=0; ii<2; ii++){

		c0[ii] = new TCanvas(Form("c0_%d",ii),Form("c0_%d",ii),1.2*3*400,2*400);
		c0[ii]->Divide(3,2);

		c1[ii] = new TCanvas(Form("c1_%d",ii),Form("c1_%d",ii),1.2*3*400,2*400);
		c1[ii]->Divide(3,2);

		c2[ii] = new TCanvas(Form("c2_%d",ii),Form("c2_%d",ii),1.2*3*400,2*400);
		c2[ii]->Divide(3,2);

		c3[ii] = new TCanvas(Form("c3_%d",ii),Form("c3_%d",ii),1.2*3*400,2*400);
		c3[ii]->Divide(3,2);

		c4[ii] = new TCanvas(Form("c4_%d",ii),Form("c4_%d",ii),1.2*3*400,2*400);
		c4[ii]->Divide(3,2);

		for (int ipt=1; ipt<npt; ipt++){

			{
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
			}

			{
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
				hmass_rot[ii][0][ipt]->Scale(1./10);
				hmass_rot[ii][1][ipt]->Scale(1./10);
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
					leg->AddEntry(hmass_rot[ii][0][ipt],"Track rotation, Unlike-sign #times0.1","P");
					leg->AddEntry(hmass_rot[ii][1][ipt],"Track rotation, Like-sign #times0.1","P");
					leg->Draw();
				}
			}

			{
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

				TF1 *f3 = new TF1("f3","gausn(0) + pol2(3)",1.7,2.2);
				f3->SetParameter(0, hmass[ii][0][ipt]->GetMaximum());
				f3->SetParameter(1, 1.87);
				f3->SetParameter(2, 0.01);

				float par0 = hmass[ii][0][ipt]->GetBinContent(hmass[ii][0][ipt]->FindBin(1.87))-hmass[ii][0][ipt]->GetBinContent(hmass[ii][0][ipt]->FindBin(1.95));
				float par1 = 1.87;
				float par2 = 0.015;

				TF1 *f4 = new TF1("f4","gaus(0) + pol2(3)",1.7,2.2);
				f4->SetParameter(0, par0);
				f4->SetParameter(1, par1);
				f4->SetParameter(2, par2);

				hmass[ii][0][ipt]->Fit(f4, "R0");
				f4->SetLineColor(2);
				f4->Draw("same");

				TF1 *f4_sig = new TF1("f4_sig","gaus(0)",1.7,2.2);
				f4_sig->SetParameter(0, f4->GetParameter(0));
				f4_sig->SetParameter(1, f4->GetParameter(1));
				f4_sig->SetParameter(2, f4->GetParameter(2));

				f4_sig->SetLineColor(1);
				f4_sig->Draw("same");

				float nsig = f4_sig->Integral(1.75,2.00)/hmass[ii][0][ipt]->GetBinWidth(1);

				{
					TLegend *leg = new TLegend(0.5,0.90-0.05*5,0.94,0.90);
					leg->SetFillStyle(0);
					leg->SetBorderSize(0);
					leg->SetTextSize(0.045);
					if ( ii==0 ){
						leg->AddEntry("",Form("w/o d_{0} cut"),"");
					}else{
						leg->AddEntry("",Form("w/ d_{0} cut"),"");
					}
					leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
					leg->AddEntry(f4,"Gaus sig. + pol2 bkg.","L");
					leg->AddEntry(f4_sig,"Gaus sig. only","L");
					leg->AddEntry("",Form("N sig.=%d",int(nsig)),"");
					leg->Draw();
				}
			}

			{
				c2[ii]->cd(ipt);
				SetPadStyle();
				gPad->SetTopMargin(0.05);

				htmp = (TH1D*)hmass_sub1[ii][ipt];
				SetHistoStyle("Mass (GeV/c^{2})","");
				htmp->GetXaxis()->SetTitleOffset(2.3);
				htmp->SetMaximum(1.1*htmp->GetMaximum());

				htmp->SetAxisRange(1.7+0.001,2.2-0.001);
				htmp->SetLineColor(1);
				hmass_sub1[ii][ipt]->Draw("p e");

				float par0 = hmass_sub1[ii][ipt]->GetBinContent(hmass_sub1[ii][ipt]->FindBin(1.87))-hmass_sub1[ii][ipt]->GetBinContent(hmass_sub1[ii][ipt]->FindBin(1.95));
				float par1 = 1.87;
				float par2 = 0.015;

				TF1 *f1 = new TF1("f1","gaus(0) + pol1(3)",1.7,2.2);
				f1->SetParameter(0, par0);
				f1->SetParameter(1, par1);
				f1->SetParameter(2, par2);

				TF1 *f2 = new TF1("f2","gaus(0) + pol2(3)",1.7,2.2);
				f2->SetParameter(0, par0);
				f2->SetParameter(1, par1);
				f2->SetParameter(2, par2);

				TF1 *f3 = new TF1("f3","gaus(0) + pol3(3)",1.7,2.2);
				f3->SetParameter(0, par0);
				f3->SetParameter(1, par1);
				f3->SetParameter(2, par2);

				hmass_sub1[ii][ipt]->Fit(f1,"R0");
				hmass_sub1[ii][ipt]->Fit(f2,"R0");
				hmass_sub1[ii][ipt]->Fit(f3,"R0");
				//f1->Draw("same");
				f2->SetLineColor(4);
				f2->Draw("same");
				f3->SetLineColor(1);
				f3->Draw("same");

				TF1 *fgaus = new TF1("fgaus","gaus",1.8,1.95);
				fgaus->SetParameter(0, f2->GetParameter(0));
				fgaus->SetParameter(1, f2->GetParameter(1));
				fgaus->SetParameter(2, f2->GetParameter(2));
				fgaus->Draw("same");

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
					leg->AddEntry(f2,"Gaus sig. + pol2 bkg.","L");
					leg->AddEntry(f3,"Gaus sig. + pol3 bkg.","L");
					//leg->AddEntry("",Form("N=%d",int(fgaus->Integral(1.8,1.95)/hmass_sub1[ii][ipt]->GetBinWidth(1))),"");
					leg->Draw();
				}
			}//

			{
				c4[ii]->cd(ipt);
				SetPadStyle();
				gPad->SetTopMargin(0.05);

				htmp = (TH1D*)hmass_sub2[ii][ipt];
				SetHistoStyle("Mass (GeV/c^{2})","");
				htmp->GetXaxis()->SetTitleOffset(2.3);
				htmp->SetMaximum(1.1*htmp->GetMaximum());

				htmp->SetAxisRange(1.7+0.001,2.2-0.001);
				htmp->SetLineColor(1);
				hmass_sub2[ii][ipt]->Draw("p e");

				float bkg = hmass_sub2[ii][ipt]->GetBinContent(hmass_sub2[ii][ipt]->FindBin(1.95)) + hmass_sub2[ii][ipt]->GetBinContent(hmass_sub2[ii][ipt]->FindBin(1.80));
				float par0 = hmass_sub2[ii][ipt]->GetBinContent(hmass_sub2[ii][ipt]->FindBin(1.87)) - 0.5*bkg;
				float par1 = 1.87;
				float par2 = 0.01;

				TF1 *f2 = new TF1("f2","gaus(0) + pol2(3)",1.7,2.2);
				f2->SetParameter(0, par0);
				f2->SetParameter(1, par1);
				f2->SetParameter(2, par2);

				TF1 *f3 = new TF1("f3","gaus(0) + pol3(3)",1.7,2.2);
				f3->SetParameter(0, par0);
				f3->SetParameter(1, par1);
				f3->SetParameter(2, par2);

				hmass_sub2[ii][ipt]->Fit(f2,"R0");
				hmass_sub2[ii][ipt]->Fit(f3,"R0");
				f2->SetLineColor(4);
				f2->Draw("same");
				f3->SetLineColor(1);
				f3->Draw("same");

				TF1 *fgaus = new TF1("fgaus","gaus",1.8,1.95);
				fgaus->SetParameter(0, f2->GetParameter(0));
				fgaus->SetParameter(1, f2->GetParameter(1));
				fgaus->SetParameter(2, f2->GetParameter(2));
				//fgaus->Draw("same");

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
					leg->AddEntry(f2,"Gaus sig. + pol2 bkg.","L");
					leg->AddEntry(f3,"Gaus sig. + pol3 bkg.","L");
					//leg->AddEntry("",Form("N=%d",int(fgaus->Integral(1.8,1.95)/hmass_sub2[ii][ipt]->GetBinWidth(1))),"");
					leg->Draw();
				}
			}//

		}

	}



}
