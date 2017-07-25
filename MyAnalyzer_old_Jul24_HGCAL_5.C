// Analyser that can be used for any pdg id
// compile using, for example: g++ -o exampleAnalyser exampleAnalyser.cxx `root-config --cflags --libs`
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
using std::string;
using std::vector;
void MyAnalyzer()
{
  // get file and tree
  //  string fileName = "/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_SinglePhoton_35GeV_20170523/NTUP/partGun_PDGid22_x120_Pt35.0To35.0_NTUP_1.root"; // specified sample
  //  TFile *inputFile = new TFile( fileName.c_str() );
  //  TTree *theTree = (TTree*)inputFile->Get("ana/hgc");
  char fileName[500];
  char location[500];
  // define histograms

  char histRho[500];
  char histRhoPhi[500];
  char histEta[500];
  char histPhi[500];
  char histNumClusterNoThreshold[500];
  char histNumCluster[500];
  char histNumMultiClusPerPhoton[500];
  char histClusterEnergy[500];
  char histClusterPt[500];
  char EneProfForLayer[500];
  int EneClusterIndex=1;

  sprintf(histRho,"hDelta_rho_E%i",EneClusterIndex);
  sprintf(histRhoPhi,"hDelta_rhophi_E%i",EneClusterIndex);
  sprintf(histEta,"hDelta_eta_E%i",EneClusterIndex);
  sprintf(histPhi,"hDelta_phi_E%i",EneClusterIndex);
  sprintf(histNumCluster,"hNumMultiClusters");
  sprintf(histNumMultiClusPerPhoton,"hNumMultiClusPerPhoton");
  sprintf(histNumClusterNoThreshold,"hNumMultiClustersNoThreshold");
  sprintf(histClusterEnergy,"hLeadingMultiClusterEnergy");
  sprintf(histClusterPt,"hLeadingMultiClusterPt");
  std::vector<TH1*> histEneProfForLayer(50);

  TH1* hNumMultiClusPerNonIntPhoton=new TH1F(histNumMultiClusPerPhoton,"Number of multiclusters with energy above 1GeV for noninteracted photon",15,0,15);
  TH1* hNumClusForStrangePhotons=new TH1F("hNumClusForStrangePhotons","Number of 2D layer clusters corresponding to strange photon(no multicluster) with E>0",50,0,50);
  TH1* hTotEneOfClusForStrangePhotons=new TH1F("hTotEneOfClusForStrangePhotons","Total energy of 2D layer clusters with E>0 corresponding to strange photon(no multicluster)",400,-2,2);
  TH2F* hNumClus2DVsTotEneForStrangePhotons=new TH2F("hNumClus2DVsTotEneForStrangePhotons","Number of 2D layer clusters with E>0GeV versus total energy of such clusters for a strange photon",50,0,50,400,-2,2);
   TH1* hEtaForStrangePhotons=new TH1F("hEtaForStrangePhotons","Eta of a strange photon(no multicluster)",200,-5,5);
   TH1* hEtaForNormalPhotons=new TH1F("hEtaForNormalPhotons","Eta of a normal photon(nonzero multicluster)",200,-5,5);
   TH1* hEtaForNormalPhotonsOnStrangeRange=new TH1F("hEtaForNormalPhotonsOnStrangeRange","Eta of a normal photon(nonzero multicluster) for abs(eta) On Strange Range",200,-5,5);
   TH1* hPhiForStrangePhotons=new TH1F("hPhiForStrangePhotons","Phi of a strange photon(no multicluster)",200,-5,5);

  //  sprintf(location,/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_SinglePhoton_35GeV_20170523/NTUP/)
  int NumInputFiles=25;
  int TotNumOfEvents=0;
  int TotPhotonsOfInterest=0;
  int Num_DeltaRho_Plus=0;
  int Num_DeltaRho_Minus=0;
  int Num_DeltaEta_Plus=0;
  int Num_DeltaEta_Minus=0;
  int NonMatch_Photon_Multicluster=0;
  int NonMatch_MorePhoton_LessMulticluster=0;
  int NonMatch_LessPhoton_MoreMulticluster=0;
  int NonMatch_MorePhoton_ZeroMulticluster=0;
  int Match_Photon_Multicluster=0;
  int NoMatchingMultiClusToPhoton=0;
  for(int j=1;j<=NumInputFiles;j++){
    sprintf(fileName,"/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_SingleParticlePtGun_PID22_SinglePhoton_Pt35GeV_20170623/NTUP/partGun_PDGid22_x80_Pt35.0To35.0_NTUP_%i.root",j);
    TFile *inputFile = TFile::Open(fileName);
    std::cout<<" inputFile "<<fileName<<endl;
    TTree *theTree = (TTree*)inputFile->Get("ana/hgc");
    // define quantities in ntuple
    vector<bool> *genpart_reachedEE = 0;
    vector<float> *genpart_pt = 0;
    vector<float> *genpart_eta = 0;
    vector<float> *genpart_phi = 0;
    vector<int> *genpart_pid = 0;
    vector<float> *genpart_dvx = 0;
    vector<float> *genpart_dvy = 0;
    vector<float> *cluster2d_energy=0;
    vector<float> *cluster2d_x=0;
    vector<float> *cluster2d_y=0;
    vector<float> *cluster2d_z=0;
    vector<float> *cluster2d_phi=0;
    vector<float> *cluster2d_eta=0;
    vector<int> *cluster2d_layer=0;
    vector<int> *cluster2d_multicluster=0;
    vector<float> *multiclus_energy=0;
    vector<float> *multiclus_phi=0;
    vector<float> *multiclus_eta=0;
    vector<float> *multiclus_pt=0;
    vector<float> *multiclus_slopeX=0;
    vector<float> *multiclus_slopeY=0;
    vector<float> *multiclus_z=0;

    theTree->SetBranchAddress("multiclus_pt", &multiclus_pt);
    theTree->SetBranchAddress("multiclus_eta", &multiclus_eta);
    theTree->SetBranchAddress("multiclus_phi", &multiclus_phi);
    theTree->SetBranchAddress("multiclus_z", &multiclus_z);
    theTree->SetBranchAddress("multiclus_slopeX", &multiclus_slopeX);
    theTree->SetBranchAddress("multiclus_slopeY", &multiclus_slopeY);
    theTree->SetBranchAddress("multiclus_energy", &multiclus_energy);
    theTree->SetBranchAddress("cluster2d_energy", &cluster2d_energy);
    theTree->SetBranchAddress("cluster2d_eta", &cluster2d_eta);
    theTree->SetBranchAddress("cluster2d_layer", &cluster2d_layer);
    theTree->SetBranchAddress("cluster2d_multicluster", &cluster2d_multicluster);

    theTree->SetBranchAddress("genpart_pt", &genpart_pt);
    theTree->SetBranchAddress("genpart_eta", &genpart_eta);
    theTree->SetBranchAddress("genpart_phi", &genpart_phi);
    theTree->SetBranchAddress("genpart_dvx", &genpart_dvx);
    theTree->SetBranchAddress("genpart_dvy", &genpart_dvy);
    theTree->SetBranchAddress("genpart_reachedEE", &genpart_reachedEE);
    theTree->SetBranchAddress("genpart_pid", &genpart_pid);

    // loop over events
    uint nEntries = theTree->GetEntries();
    //    std::cout<<" nEntries "<<nEntries<<endl;
    for( int evtIndex = 0; evtIndex < nEntries; evtIndex++ )
      {
	TotNumOfEvents++;
	// get values
	theTree->GetEntry( evtIndex );
	//	std::cout<<" evtIndex "<<evtIndex<<" Number of multiclusters "<<multiclus_pt->size()<<endl;
	std::vector<int> leadingcluster_index;
	std::vector<float> leadingcluster_x;
	std::vector<float> leadingcluster_y;
	std::vector<float> leadingcluster_phi;
	std::vector<float> leadingcluster_eta;
	int PhotonsOfInterest=0;
	int NumMultiClusters=0;
	int NumMultiClustersNoThreshold=0;

	//loop over gen particles to check whether a gen particle is a photon having pt=35 GeV and is non interacted in a tracker.
	for( uint genIndex = 0; genIndex < genpart_pt->size(); genIndex++ ) {
	  float photon_Eta=0;
	  float photon_Phi=0;
	  float photon_X=0;
	  float photon_Y=0;
	  if(genpart_pid->at(genIndex)==22 && genpart_reachedEE->at(genIndex)==1 && genpart_pt->at(genIndex)==35.){
	    //	    std::cout<<" evt "<<evtIndex<<" gen particles "<<genpart_pt->size()<<endl;
	    std::vector<float> leadingcluster_energy;
	    float CumulativeEnergy=0;
	    int MultiClusPerPhoton=0;
	    PhotonsOfInterest++;
	    TotPhotonsOfInterest++;
	    photon_Eta=genpart_eta->at(genIndex);
	    photon_Phi=genpart_phi->at(genIndex);
	    photon_X=genpart_dvx->at(genIndex);
	    photon_Y=genpart_dvy->at(genIndex);
	    int jx=0;

	    //Determines number of multiclusters above 1GeV.
	      int startIndex=-1;
	      float MaxFound=9999999.;
	      for(int multiclusEneIndex = 0; multiclusEneIndex < multiclus_energy->size(); multiclusEneIndex++ ) {
		  if(multiclus_energy->at(multiclusEneIndex)>1.0 && ((photon_Eta>0 && multiclus_eta->at(multiclusEneIndex)>0) || (photon_Eta<0 && multiclus_eta->at(multiclusEneIndex)<0))){
		    NumMultiClusters++;
		    MultiClusPerPhoton++;
		  }
	      }
	      hNumMultiClusPerNonIntPhoton->Fill(MultiClusPerPhoton);

		//Finds leading energy 2D cluster per layer and fills into a vector leadingcluster_energy if no multicluster is matched to a noninteracted photon
		if(MultiClusPerPhoton==0){
		  int NumClusForStrangePhoton=0;
		  float TotEneClusForStrangePhoton=0;
		  for(int LIndex=0;LIndex<50;LIndex++){ 
		    int start2DIndex=-1;
		    float Max2DFound=9999999.;
		    float Max_cluster2D_energy=0;
		    for(int clusEneIndex = 0; clusEneIndex < cluster2d_energy->size(); clusEneIndex++ ) {
		      if((photon_Eta>0 && cluster2d_eta->at(clusEneIndex)>0 && cluster2d_layer->at(clusEneIndex)==(LIndex+1)) || (photon_Eta<0 && cluster2d_eta->at(clusEneIndex)<0 && cluster2d_layer->at(clusEneIndex)==(LIndex+1))){
			if(clusEneIndex!=start2DIndex && cluster2d_energy->at(clusEneIndex)<=Max2DFound){
			  if(cluster2d_energy->at(clusEneIndex)>Max_cluster2D_energy){
			    Max_cluster2D_energy=cluster2d_energy->at(clusEneIndex);
			    start2DIndex=clusEneIndex;
			  }
			}
		      }
		    } //end of loop over clusEneIndex
		    Max2DFound=Max_cluster2D_energy;
		    if(Max2DFound>0.){		   
		      leadingcluster_energy.push_back(Max2DFound);
		      CumulativeEnergy +=Max2DFound;
		    }
		  }
		  hNumClusForStrangePhotons->Fill(leadingcluster_energy.size());
		  hTotEneOfClusForStrangePhotons->Fill(CumulativeEnergy);
		  hNumClus2DVsTotEneForStrangePhotons->Fill(leadingcluster_energy.size(),CumulativeEnergy);
		  hEtaForStrangePhotons->Fill(photon_Eta);
		  hPhiForStrangePhotons->Fill(photon_Phi);
		}
		else{
		  hEtaForNormalPhotons->Fill(photon_Eta);
		  if((1.4<abs(photon_Eta) && abs(photon_Eta)<1.6) || (2.9<abs(photon_Eta) && abs(photon_Eta)<3.1))
		    hEtaForNormalPhotonsOnStrangeRange->Fill(photon_Eta);
		}
	  }
	}
      }
  }
 // draw histograms
  TFile *ouputFile = new TFile("BasicPlots.root","RECREATE");
  hNumMultiClusPerNonIntPhoton->Write();
  hNumClusForStrangePhotons->Write();
  hTotEneOfClusForStrangePhotons->Write();
  hNumClus2DVsTotEneForStrangePhotons->Write();
  hEtaForStrangePhotons->Write();
  hPhiForStrangePhotons->Write();
  hEtaForNormalPhotons->Write();
  hEtaForNormalPhotonsOnStrangeRange->Write();
  ouputFile->Close();
  
}

