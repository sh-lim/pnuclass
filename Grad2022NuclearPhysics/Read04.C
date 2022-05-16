#define NMAX_TRK 500

#include <TTree.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TMath.h>

#include <iostream>
#include <fstream>

using namespace std;

void Read04(const char *listname="file.lst"){

	gRandom = new TRandom3(0);

	const float const_pi = TMath::Pi();

	const int nz = 10;

	float const_mass_pion = 0.139; 
	float const_mass_kaon = 0.493;

	ifstream flist;
	flist.open(listname);

	char fname[500];

	bool _IsMB;
	bool _IsHMV0;
	float _VertexZ;
	float _Centrality;

	int _nTrack; 
	float _TrackPt[NMAX_TRK];
	float _TrackEta[NMAX_TRK];
	float _TrackPhi[NMAX_TRK];
	int _TrackBayesianPID[NMAX_TRK];
	int _TrackCharge[NMAX_TRK];
	float _TrackDCAXY[NMAX_TRK];

	const int npt = 7;
	const double ptbin[npt+1] = {0, 1, 2, 4, 6, 8, 12, 24};

	TH2D *hmass_pt_cut0[2];
	TH2D *hmass_pt_cut1[2];

	TH2D *hmass_pt_rot_cut0[2];
	TH2D *hmass_pt_rot_cut1[2];

	TH2D *hmass_pt_mix_cut0[2];
	TH2D *hmass_pt_mix_cut1[2];

	TH2D *hd01_pt[2];
	TH2D *hd02_pt[2];

	for (int ii=0; ii<2; ii++){
		hmass_pt_cut0[ii] = new TH2D(Form("hmass_pt_cut0_chg%d",ii), "", 200, 1.5, 2.5, npt, ptbin);
		hmass_pt_cut1[ii] = new TH2D(Form("hmass_pt_cut1_chg%d",ii), "", 200, 1.5, 2.5, npt, ptbin);

		hmass_pt_rot_cut0[ii] = new TH2D(Form("hmass_pt_rot_cut0_chg%d",ii), "", 200, 1.5, 2.5, npt, ptbin);
		hmass_pt_rot_cut1[ii] = new TH2D(Form("hmass_pt_rot_cut1_chg%d",ii), "", 200, 1.5, 2.5, npt, ptbin);

		hmass_pt_mix_cut0[ii] = new TH2D(Form("hmass_pt_mix_cut0_chg%d",ii), "", 200, 1.5, 2.5, npt, ptbin);
		hmass_pt_mix_cut1[ii] = new TH2D(Form("hmass_pt_mix_cut1_chg%d",ii), "", 200, 1.5, 2.5, npt, ptbin);

		hd01_pt[ii] = new TH2D(Form("hd01_pt_chg%d",ii), "", 200, -1e3, 1e3, npt, ptbin);
		hd02_pt[ii] = new TH2D(Form("hd02_pt_chg%d",ii), "", 200, -1e6, 1e6, npt, ptbin);
	}

	vector<float> vec_pt[nz];
	vector<float> vec_eta[nz];
	vector<float> vec_phi[nz];
	vector<float> vec_d0[nz];
	vector<int> vec_chg[nz];
	vector<int> vec_pid[nz];

	TH1D *htrig = new TH1D("htrig","",5,0,5);

	while ( flist >> fname ){

		TFile *infile = new TFile(fname,"read");
		if ( !infile->IsOpen() ){
			continue;
		}

		cout << "OPEN: " << fname << endl;

		TTree *T = (TTree*)infile->Get("tree");
		if ( !T ){
			infile->Close();
			delete infile;
		}

		T->SetBranchAddress("IsMB",&_IsMB);
		T->SetBranchAddress("IsHMV0",&_IsHMV0);
		T->SetBranchAddress("VertexZ",&_VertexZ);
		T->SetBranchAddress("Centrality",&_Centrality);
		T->SetBranchAddress("nTrack",&_nTrack);
		T->SetBranchAddress("TrackPt",_TrackPt);
		T->SetBranchAddress("TrackEta",_TrackEta);
		T->SetBranchAddress("TrackPhi",_TrackPhi);
		T->SetBranchAddress("TrackBayesianPID",_TrackBayesianPID);
		T->SetBranchAddress("TrackCharge",_TrackCharge);
		T->SetBranchAddress("TrackDCAXY",_TrackDCAXY);

		int nentries = T->GetEntries();

		//event loop
		for (int ien=0; ien<nentries; ien++){
		//for (int ien=0; ien<0; ien++){

			T->GetEntry(ien);

			if ( fabs(_VertexZ)>10.0 ) continue;

			int zbin = int((_VertexZ+10.0)/2);
			if ( zbin<0 || zbin>=nz ) continue;

			//if ( !_IsHMV0 ) continue;

			if ( _IsMB ){
				htrig->Fill(0);
			}else if ( _IsHMV0 ){
				htrig->Fill(1);
				continue;
			}else{
				htrig->Fill(2);
				continue;
			}

			for (int itrk=0; itrk<_nTrack; itrk++){

				if ( fabs(_TrackEta[itrk])>0.8 ) continue;
				if ( _TrackPt[itrk]<0.5 ) continue;

				int ipid = _TrackBayesianPID[itrk];
				if ( !(ipid==0 || ipid==1) ) continue;

				TLorentzVector ilvec;
				TLorentzVector ilvec_rot[10];
				if ( ipid==0 ){
					ilvec.SetPtEtaPhiM(_TrackPt[itrk], _TrackEta[itrk], _TrackPhi[itrk], const_mass_pion); 
				}else if ( ipid==1 ){
					ilvec.SetPtEtaPhiM(_TrackPt[itrk], _TrackEta[itrk], _TrackPhi[itrk], const_mass_kaon); 
				}

				for (int iter=0; iter<10; iter++){
					ilvec_rot[iter] = ilvec;
					float rot_ang = const_pi + (gRandom->Rndm()-0.5)*const_pi/3;
					ilvec_rot[iter].RotateZ(rot_ang);
				}

				for (int jtrk=itrk+1; jtrk<_nTrack; jtrk++){

					if ( fabs(_TrackEta[jtrk])>0.8 ) continue;
					if ( _TrackPt[jtrk]<0.5 ) continue;

					int jpid = _TrackBayesianPID[jtrk];
					if ( !(jpid==0 || jpid==1) ) continue;
					if ( ipid==jpid ) continue;

					TLorentzVector jlvec;
					TLorentzVector jlvec_rot[10];
					if ( jpid==0 ){
						jlvec.SetPtEtaPhiM(_TrackPt[jtrk], _TrackEta[jtrk], _TrackPhi[jtrk], const_mass_pion); 
					}else if ( jpid==1 ){
						jlvec.SetPtEtaPhiM(_TrackPt[jtrk], _TrackEta[jtrk], _TrackPhi[jtrk], const_mass_kaon); 
					}

					for (int jter=0; jter<10; jter++){
						jlvec_rot[jter] = jlvec;
						float rot_ang = const_pi + (gRandom->Rndm()-0.5)*const_pi/3;
						jlvec_rot[jter].RotateZ(rot_ang);
					}

					TLorentzVector plvec = ilvec + jlvec;

					float p_mass = plvec.M();
					float p_pt = plvec.Pt();

					if ( p_mass>1.5 && p_mass<2.5 ){

						bool bgood = true;

						if ( p_pt<2 ){
							//do nothing
						}else if ( p_pt<12 ){
							if ( _TrackPt[itrk]<0.7 || _TrackPt[jtrk]<0.7 ) bgood = false;
						}else{
							if ( _TrackPt[itrk]<0.6 || _TrackPt[jtrk]<0.6 ) bgood = false;
						}

						float i_d0 = _TrackDCAXY[itrk]*1e4;
						float j_d0 = _TrackDCAXY[jtrk]*1e4;

						if ( bgood ){
							if ( _TrackCharge[itrk]==_TrackCharge[jtrk] ){
								hmass_pt_cut0[1]->Fill(p_mass, p_pt);
								hd01_pt[1]->Fill(i_d0, p_pt);
								hd01_pt[1]->Fill(j_d0, p_pt);
								hd02_pt[1]->Fill(i_d0*j_d0, p_pt);
							}else{
								hmass_pt_cut0[0]->Fill(p_mass, p_pt);
								hd01_pt[0]->Fill(i_d0, p_pt);
								hd01_pt[0]->Fill(j_d0, p_pt);
								hd02_pt[0]->Fill(i_d0*j_d0, p_pt);
							}
						}//bgood

						if ( p_pt<2 ){
							if ( i_d0*j_d0>-2000 ) bgood = false;
						}else if ( p_pt<4 ){
							if ( i_d0*j_d0>-20000 ) bgood = false;
						}else if ( p_pt<12 ){
							if ( i_d0*j_d0>-5000 ) bgood = false;
						}else{
							//do nothing
						}

						if ( bgood ){
							if ( _TrackCharge[itrk]==_TrackCharge[jtrk] ){
								hmass_pt_cut1[1]->Fill(p_mass, p_pt);
							}else{
								hmass_pt_cut1[0]->Fill(p_mass, p_pt);
							}
						}//bgood
					}//p_mass

					for (int jter=0; jter<10; jter++){

						TLorentzVector plvec_rot;
						if ( ipid==1 ){
							plvec_rot = ilvec_rot[jter] + jlvec;
						}else{
							plvec_rot = ilvec + jlvec_rot[jter];
						}

						float p_rot_mass = plvec_rot.M();
						float p_rot_pt = plvec_rot.Pt();

						if ( p_rot_mass>1.5 && p_rot_mass<2.5 ){

							bool bgood = true;

							if ( p_rot_pt<2 ){
								//do nothing
							}else if ( p_rot_pt<12 ){
								if ( _TrackPt[itrk]<0.7 || _TrackPt[jtrk]<0.7 ) bgood = false;
							}else{
								if ( _TrackPt[itrk]<0.6 || _TrackPt[jtrk]<0.6 ) bgood = false;
							}

							float i_d0 = _TrackDCAXY[itrk]*1e4;
							float j_d0 = _TrackDCAXY[jtrk]*1e4;

							if ( bgood ){
								if ( _TrackCharge[itrk]==_TrackCharge[jtrk] ){
									hmass_pt_rot_cut0[1]->Fill(p_rot_mass, p_rot_pt);
								}else{
									hmass_pt_rot_cut0[0]->Fill(p_rot_mass, p_rot_pt);
								}
							}//bgood

							if ( p_rot_pt<2 ){
								if ( i_d0*j_d0>-2000 ) bgood = false;
							}else if ( p_rot_pt<4 ){
								if ( i_d0*j_d0>-20000 ) bgood = false;
							}else if ( p_rot_pt<12 ){
								if ( i_d0*j_d0>-5000 ) bgood = false;
							}else{
								//do nothing
							}

							if ( bgood ){
								if ( _TrackCharge[itrk]==_TrackCharge[jtrk] ){
									hmass_pt_rot_cut1[1]->Fill(p_rot_mass, p_rot_pt);
								}else{
									hmass_pt_rot_cut1[0]->Fill(p_rot_mass, p_rot_pt);
								}
							}//bgood
						}//p_mass

					}//jter

				}//jtrk

				if ( vec_pt[zbin].size()>0 ){

					for (unsigned int jj=0; jj<vec_pt[zbin].size(); jj++){

						int jpid = vec_pid[zbin][jj];
						if ( ipid==jpid ) continue;

						TLorentzVector jlvec;
						if ( jpid==0 ){
							jlvec.SetPtEtaPhiM(vec_pt[zbin][jj], vec_eta[zbin][jj], vec_phi[zbin][jj], const_mass_pion);
						}else if ( jpid==1 ){
							jlvec.SetPtEtaPhiM(vec_pt[zbin][jj], vec_eta[zbin][jj], vec_phi[zbin][jj], const_mass_kaon);
						}

						TLorentzVector plvec = ilvec + jlvec;

						float p_mass = plvec.M();
						float p_pt = plvec.Pt();

						if ( p_mass>1.5 && p_mass<2.5 ){

							bool bgood = true;

							if ( p_pt<2 ){
								//do nothing
							}else if ( p_pt<12 ){
								if ( _TrackPt[itrk]<0.7 || vec_pt[zbin][jj]<0.7 ) bgood = false;
							}else{
								if ( _TrackPt[itrk]<0.6 || vec_pt[zbin][jj]<0.6 ) bgood = false;
							}

							float i_d0 = _TrackDCAXY[itrk]*1e4;
							float j_d0 = vec_d0[zbin][jj]*1e4;

							if ( bgood ){
								if ( _TrackCharge[itrk]==vec_chg[zbin][jj] ){
									hmass_pt_mix_cut0[1]->Fill(p_mass, p_pt);
								}else{
									hmass_pt_mix_cut0[0]->Fill(p_mass, p_pt);
								}
							}//bgood

							if ( p_pt<2 ){
								if ( i_d0*j_d0>-2000 ) bgood = false;
							}else if ( p_pt<4 ){
								if ( i_d0*j_d0>-20000 ) bgood = false;
							}else if ( p_pt<12 ){
								if ( i_d0*j_d0>-5000 ) bgood = false;
							}else{
								//do nothing
							}

							if ( bgood ){
								if ( _TrackCharge[itrk]==vec_chg[zbin][jj] ){
									hmass_pt_mix_cut1[1]->Fill(p_mass, p_pt);
								}else{
									hmass_pt_mix_cut1[0]->Fill(p_mass, p_pt);
								}
							}//bgood
						}//mass

					}//jj

				}//mixed_event

			}//itrk

			//mixed event pool
			vec_pt[zbin].clear();
			vec_eta[zbin].clear();
			vec_phi[zbin].clear();
			vec_d0[zbin].clear();
			vec_chg[zbin].clear();
			vec_pid[zbin].clear();
			for (int itrk=0; itrk<_nTrack; itrk++){

				if ( fabs(_TrackEta[itrk])>0.8 ) continue;
				if ( _TrackPt[itrk]<0.5 ) continue;

				int ipid = _TrackBayesianPID[itrk];
				if ( !(ipid==0 || ipid==1) ) continue;

				vec_pt[zbin].push_back(_TrackPt[itrk]);
				vec_eta[zbin].push_back(_TrackEta[itrk]);
				vec_phi[zbin].push_back(_TrackPhi[itrk]);
				vec_d0[zbin].push_back(_TrackDCAXY[itrk]);
				vec_chg[zbin].push_back(_TrackCharge[itrk]);
				vec_pid[zbin].push_back(ipid);
			}//itrk


		}//ien

		infile->Close();
		delete infile;

	}//

	TFile *outfile = new TFile("outfile_Read04.root","recreate");

	htrig->Write();

	for (int ii=0; ii<2; ii++){
		hmass_pt_cut0[ii]->Write();
		hmass_pt_cut1[ii]->Write();

		hmass_pt_rot_cut0[ii]->Write();
		hmass_pt_rot_cut1[ii]->Write();

		hmass_pt_mix_cut0[ii]->Write();
		hmass_pt_mix_cut1[ii]->Write();

		hd01_pt[ii]->Write();
		hd02_pt[ii]->Write();
	}

	outfile->Close();


}
