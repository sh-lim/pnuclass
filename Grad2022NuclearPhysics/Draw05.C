#include "Style.h"

void Draw05(){

	bool bFIX = true;

	gStyle->SetOptStat(0);

	const int npt = 7;
	const double ptbin[npt+1] = {0, 1, 2, 4, 6, 8, 12, 24};

	const float par_width[npt] = {0.0077, 0.0079, 0.0092, 0.0107, 0.0118, 0.0134, 0.0155};

	const float BR = 0.0395;

	const float visibleXection = 55*1000.; //ub

	//TFile *infile = new TFile("outfile_runhist_pp13TeV_v2_try105.root","read");
	TFile *infile = new TFile("outfile_runhist_pp13TeV_v2_try106.root","read");

	TH2D *hmass_pt[2][2];
	TH2D *hmass_pt_rot[2][2];
	TH2D *hmass_pt_mix[2][2];

	TH1D *hmass[2][2][npt];
	TH1D *hmass_rot[2][2][npt];
	TH1D *hmass_mix[2][2][npt];

	TH1D *hmass_sub1[2][npt];
	TH1D *hmass_sub2[2][npt];

	TH1D *hnsig_sub1[2];
	TH1D *hnsig_sub2[2];

	TH1D *hntrig = (TH1D*)infile->Get("htrig");

	for (int ii=0; ii<2; ii++){

		hmass_pt[0][ii] = (TH2D*)infile->Get(Form("hmass_pt_cut0_chg%d",ii));
		hmass_pt[1][ii] = (TH2D*)infile->Get(Form("hmass_pt_cut1_chg%d",ii));

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

			hmass_mix[0][ii][ipt] = (TH1D*)hmass_pt_mix[0][ii]->ProjectionX(Form("hmass_mix_cut0_%d_%d",ii,ipt),ipt+1,ipt+1);
			hmass_mix[1][ii][ipt] = (TH1D*)hmass_pt_mix[1][ii]->ProjectionX(Form("hmass_mix_cut1_%d_%d",ii,ipt),ipt+1,ipt+1);

			hmass[0][ii][ipt]->Sumw2();
			hmass[1][ii][ipt]->Sumw2();

			hmass_mix[0][ii][ipt]->Sumw2();
			hmass_mix[1][ii][ipt]->Sumw2();

			hmass[0][ii][ipt]->SetMarkerStyle(24);
			hmass[0][ii][ipt]->SetMarkerColor(1+ii);
			hmass[0][ii][ipt]->SetLineColor(1+ii);
			hmass[1][ii][ipt]->SetMarkerStyle(24);
			hmass[1][ii][ipt]->SetMarkerColor(1+ii);
			hmass[1][ii][ipt]->SetLineColor(1+ii);

			hmass_mix[0][ii][ipt]->SetMarkerStyle(27);
			hmass_mix[0][ii][ipt]->SetMarkerColor(1+ii);
			hmass_mix[0][ii][ipt]->SetLineColor(1+ii);
			hmass_mix[1][ii][ipt]->SetMarkerStyle(27);
			hmass_mix[1][ii][ipt]->SetMarkerColor(1+ii);
			hmass_mix[1][ii][ipt]->SetLineColor(1+ii);
		}//ipt
	}//ii

	for (int jj=1; jj<2; jj++){
		for (int ipt=0; ipt<npt; ipt++){

			hmass_sub1[jj][ipt] = (TH1D*)hmass[jj][0][ipt]->Clone(Form("hmass_sub1_%d_%d",jj,ipt));

			hmass_sub1[jj][ipt]->Add(hmass[jj][1][ipt], -1);

			hmass_sub2[jj][ipt] = (TH1D*)hmass[jj][0][ipt]->Clone(Form("hmass_sub2_%d_%d",jj,ipt));

			TH1D *hbkg_mix = (TH1D*)hmass_mix[jj][0][ipt]->Clone(Form("hbkg_mix_%d_%d",jj,ipt));
			hbkg_mix->Add(hmass_mix[jj][1][ipt]);

			float norm1 = hmass[jj][0][ipt]->Integral(hmass[jj][0][ipt]->FindBin(1.65),hmass[jj][0][ipt]->FindBin(1.7));
			float norm2 = hbkg_mix->Integral(hbkg_mix->FindBin(1.65), hbkg_mix->FindBin(1.7));

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

	TFile *infile_eff = new TFile("outfile_acceff.root","read");
	TH1D *heff[2];
	heff[0] = (TH1D*)infile_eff->Get("hacceff_1");
	heff[1] = (TH1D*)infile_eff->Get("hacceff_2");

	for (int ii=1; ii<2; ii++){

		c0[ii] = new TCanvas(Form("c0_%d",ii),Form("c0_%d",ii),1.2*4*400,2*400);
		c0[ii]->Divide(4,2);

		c1[ii] = new TCanvas(Form("c1_%d",ii),Form("c1_%d",ii),1.2*4*400,2*400);
		c1[ii]->Divide(4,2);

		c2[ii] = new TCanvas(Form("c2_%d",ii),Form("c2_%d",ii),1.2*4*400,2*400);
		c2[ii]->Divide(4,2);

		c3[ii] = new TCanvas(Form("c3_%d",ii),Form("c3_%d",ii),1.2*4*400,2*400);
		c3[ii]->Divide(4,2);

		c4[ii] = new TCanvas(Form("c4_%d",ii),Form("c4_%d",ii),1.2*4*400,2*400);
		c4[ii]->Divide(4,2);

		hnsig_sub1[ii] = new TH1D(Form("hnsig_sub1_cut%d",ii),"",npt,ptbin); 
		hnsig_sub2[ii] = new TH1D(Form("hnsig_sub2_cut%d",ii),"",npt,ptbin); 

		for (int ipt=0; ipt<npt; ipt++){

			{
				c0[ii]->cd(ipt+1);
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
					if ( ipt>3 || ipt==0 ){
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
				c1[ii]->cd(ipt+1);
				SetPadStyle();
				gPad->SetTopMargin(0.05);

				htmp = (TH1D*)hmass[ii][1][ipt];
				SetHistoStyle("Mass (GeV/c^{2})","");
				htmp->GetXaxis()->SetTitleOffset(2.3);

				htmp->SetAxisRange(1.7+0.001,2.2-0.001);
				htmp->SetMinimum(0);
				htmp->SetMaximum(1.1*htmp->GetMaximum());

				hmass[ii][1][ipt]->Draw("p e");
				hmass_mix[ii][0][ipt]->Scale(1./10);
				hmass_mix[ii][1][ipt]->Scale(1./10);
				hmass_mix[ii][0][ipt]->Draw("p e same");
				hmass_mix[ii][1][ipt]->Draw("p e same");

				if ( ipt<4 ){
					TLegend *leg = new TLegend(0.35,0.15,0.94,0.15+0.05*5);
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
					leg->Draw();
				}
			}

			{
				c3[ii]->cd(ipt+1);
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

				hmass[ii][0][ipt]->Fit(f4, "R0Q");
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
				c2[ii]->cd(ipt+1);
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

				if ( bFIX ){
					f2->FixParameter(2, par_width[ipt]);
				}

				TF1 *f3 = new TF1("f3","gaus(0) + pol3(3)",1.7,2.2);
				f3->SetParameter(0, par0);
				f3->SetParameter(1, par1);
				f3->SetParameter(2, par2);

				//hmass_sub1[ii][ipt]->Fit(f1,"R0Q");
				hmass_sub1[ii][ipt]->Fit(f2,"R0");
				//hmass_sub1[ii][ipt]->Fit(f3,"R0");
				//f1->Draw("same");
				f2->SetLineColor(4);
				f2->Draw("same");
				//f3->SetLineColor(1);
				//f3->Draw("same");

				TF1 *fgaus = new TF1("fgaus","gaus",1.8,1.95);
				fgaus->SetParameter(0, f2->GetParameter(0));
				fgaus->SetParameter(1, f2->GetParameter(1));
				fgaus->SetParameter(2, f2->GetParameter(2));
				fgaus->Draw("same");

				TF1 *fbkg = new TF1("fbkg","pol2",1.7,2.2);
				fbkg->SetParameter(0, f2->GetParameter(3));
				fbkg->SetParameter(1, f2->GetParameter(4));
				fbkg->SetParameter(2, f2->GetParameter(5));
				fbkg->SetLineColor(1);
				fbkg->SetLineStyle(2);
				fbkg->Draw("same");

				float nsig = fgaus->Integral(1.8,1.95)/hmass[ii][0][ipt]->GetBinWidth(1);
				float nsig_err = nsig * (f2->IntegralError(1.8,1.95)/f2->Integral(1.8,1.95));

				hnsig_sub1[ii]->SetBinContent(ipt+1, nsig);
				hnsig_sub1[ii]->SetBinError(ipt+1, nsig_err);

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
					leg->AddEntry("",Form("N sig.=%d#pm%d",int(nsig),int(nsig_err)),"");
					//leg->AddEntry(f3,"Gaus sig. + pol3 bkg.","L");
					leg->Draw();
				}
			}//

			//continue;

			{
				c4[ii]->cd(ipt+1);
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
				float par2 = 0.007 + 0.001*ipt;

				TF1 *f2 = new TF1("f2","gaus(0) + pol2(3)",1.7,2.2);
				f2->SetParameter(0, par0);
				f2->SetParameter(1, par1);
				f2->SetParameter(2, par2);

				TF1 *f3 = new TF1("f3","gaus(0) + pol3(3)",1.7,2.2);
				f3->SetParameter(0, par0);
				f3->SetParameter(1, par1);
				f3->SetParameter(2, par2);
				if ( bFIX ){
					f3->FixParameter(2, par_width[ipt]);
				}

				hmass_sub2[ii][ipt]->Fit(f2,"R0Q");
				hmass_sub2[ii][ipt]->Fit(f3,"R0Q");
				f2->SetLineColor(4);
				//f2->Draw("same");
				f3->SetLineColor(1);
				f3->Draw("same");

				TF1 *fgaus = new TF1("fgaus","gaus",1.8,1.95);
				fgaus->SetParameter(0, f3->GetParameter(0));
				fgaus->SetParameter(1, f3->GetParameter(1));
				fgaus->SetParameter(2, f3->GetParameter(2));
				fgaus->Draw("same");

				TF1 *fbkg = new TF1("fbkg","pol3",1.7,2.2);
				fbkg->SetParameter(0, f3->GetParameter(3));
				fbkg->SetParameter(1, f3->GetParameter(4));
				fbkg->SetParameter(2, f3->GetParameter(5));
				fbkg->SetParameter(3, f3->GetParameter(6));
				fbkg->SetLineColor(4);
				fbkg->SetLineStyle(2);
				fbkg->Draw("same");

				float nsig = fgaus->Integral(1.8,1.95)/hmass[ii][0][ipt]->GetBinWidth(1);
				float nsig_err = nsig * (f3->IntegralError(1.8,1.95)/f3->Integral(1.8,1.95));

				hnsig_sub2[ii]->SetBinContent(ipt+1, nsig);
				hnsig_sub2[ii]->SetBinError(ipt+1, nsig_err);

				{
					TLegend *leg;
					if ( ipt<2 ){
						leg = new TLegend(0.5,0.25,0.94,0.25+0.05*4);
					}else{
						leg = new TLegend(0.5,0.90-0.05*4,0.94,0.90);
					}
					leg->SetFillStyle(0);
					leg->SetBorderSize(0);
					leg->SetTextSize(0.045);
					if ( ii==0 ){
						leg->AddEntry("",Form("w/o d_{0} cut"),"");
					}else{
						leg->AddEntry("",Form("w/ d_{0} cut"),"");
					}
					leg->AddEntry("",Form("%g<p_{T}<%g GeV/c",ptbin[ipt],ptbin[ipt+1]),"");
					//leg->AddEntry(f2,"Gaus sig. + pol2 bkg.","L");
					leg->AddEntry(f3,"Gaus sig. + pol3 bkg.","L");
					leg->AddEntry("",Form("N sig.=%d#pm%d",int(nsig),int(nsig_err)),"");
					leg->Draw();
				}
			}//

		}//ipt

		hnsig_sub1[ii]->Divide(heff[0]);
		hnsig_sub2[ii]->Divide(heff[0]);

		hnsig_sub1[ii]->Scale(1,"width");
		hnsig_sub2[ii]->Scale(1,"width");

		float factor = (1/2.0)*(1/1.0)*(1/BR);
		float norm = 1.1*(hntrig->GetBinContent(1));

		hnsig_sub1[ii]->Scale(factor/norm*visibleXection);
		hnsig_sub2[ii]->Scale(factor/norm*visibleXection);

	}//ii


	TCanvas *c10 = new TCanvas("c10","c10",1.2*500,500);
	SetPadStyle();
	gPad->SetLogy();

	htmp = (TH1D*)gPad->DrawFrame(0, 5e-2, 24, 5e2);
	SetHistoStyle("p_{T} (GeV/c)","d^{2}#sigma/dp_{T}dy (#mub GeV^{-1} c)","",24,20);

	hnsig_sub1[1]->SetMarkerStyle(24);
	hnsig_sub1[1]->SetMarkerColor(1);
	hnsig_sub1[1]->SetLineColor(1);
	hnsig_sub1[1]->Draw("p same");

	hnsig_sub2[1]->SetMarkerStyle(25);
	hnsig_sub2[1]->SetMarkerColor(4);
	hnsig_sub2[1]->SetLineColor(4);
	hnsig_sub2[1]->Draw("p same");

	{
		TLegend *leg = new TLegend(0.5,0.90-0.05*4,0.94,0.90);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.045);
		leg->AddEntry("","pp 13 TeV","");
		leg->AddEntry("","Inclusive D^{0}, |y|<0.5","");
		leg->AddEntry(hnsig_sub1[1],"Like-sign background","P");
		leg->AddEntry(hnsig_sub2[1],"Mixed-event background","P");
		leg->Draw();
	}


}
