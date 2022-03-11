#define NMAX_TRK 500

void Read00(const char *listname="file.lst"){

	ifstream flist;
	flist.open(listname);

	char fname[500];

	bool _IsMB;
	bool _IsHMV0;
	float _VertexZ;

	int _nTrack; 
	float _Centrality;
	float _TrackPt[NMAX_TRK];
	float _TrackEta[NMAX_TRK];

	TH1D *hcent = new TH1D("hcent","",100,0,100);
	TH2D *hcent_pt = new TH2D("hcent_pt","",100,0,100,50,0,5);


	while ( flist >> fname ){

		TFile *infile = new TFile(fname,"read");
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

		int nentries = T->GetEntries();

		//event loop
		for (int ien=0; ien<nentries; ien++){

			T->GetEntry(ien);

			if ( !_IsMB ) continue;

			hcent->Fill(_Centrality);

			//track loop
			for (int itrk=0; itrk<_nTrack; itrk++){

				if ( _TrackPt[itrk]>5.0 ) continue;
				if ( fabs(_TrackEta[itrk])>0.5 ) continue;

				hcent_pt->Fill(_Centrality, _TrackPt[itrk]);

			}//itrk

		}//ien

		infile->Close();
		delete infile;

	}//

	TFile *outfile = new TFile("outfile_Read00.root","recreate");
	hcent->Write();
	hcent_pt->Write();
	outfile->Close();


}
