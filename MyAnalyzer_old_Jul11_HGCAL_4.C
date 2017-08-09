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
void MyAnalyzer_old_Jul11_HGCAL_4()
{
  // get file and tree
  //  string fileName = "/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_SinglePhoton_35GeV_20170523/NTUP/partGun_PDGid22_x120_Pt35.0To35.0_NTUP_1.root"; // specified sample
  //  TFile *inputFile = new TFile( fileName.c_str() );
  //  TTree *theTree = (TTree*)inputFile->Get("ana/hgc");
  char fileName[500];
  char location[500];
  // define histograms
  TH1* hMultiPt = new TH1F( "hMultiPt", "Pt of a leading multicluster", 200, 0., 50. );
  TH1* hDeltaEta_Plus =new TH1F("hDeltaEta_Plus"," Number of events with (eta_cluster-eta_gen)>0 ",2000,0,2000); 
  TH1* hDeltaEta_Minus =new TH1F("hDeltaEta_Minus"," Number of events with (eta_cluster-eta_gen)<0 ",2000,0,2000); 
  TH1* hDeltaRho_Plus =new TH1F("hDeltaRho_Plus"," Number of events with (rho_cluster-rho_gen)>0 ",2000,0,2000); 
  TH1* hDeltaRho_Minus =new TH1F("hDeltaRho_Minus"," Number of events with (rho_cluster-rho_gen)<0 ",2000,0,2000); 

  TH1* hNonInteractedPhotons = new TH1F( "hNonInteractedPhotons","Number of NonInteractedPhotons",5,0,5);
  TH1* hDiff_etaC_etaM = new TH1F("hDiff_etaC_etaM","Difference between eta_measured and eta_calculated for a leading energy cluster",80,-0.02,0.02);
  TH1* hDelta_rhoM_etaPlus= new TH1F("hDelta_rhoM_etaPlus","Difference(rho_measured-rho_true) for a leading cluster with eta>0 and rho measured ",240,-3.,3.);
  TH1* hDelta_rhoM_etaMinus= new TH1F("hDelta_rhoM_etaMinus","Difference(rho_measured-rho_true) for a leading cluster with eta<0 and rho measured",240,-3.,3.);
  TH1* hDelta_rhoC_etaPlus= new TH1F("hDelta_rhoC_etaPlus","Difference(rho_measured-rho_true) for a leading cluster with eta>0 and rho calculated ",240,-3.,3.);
  TH1* hDelta_rhoC_etaMinus= new TH1F("hDelta_rhoC_etaMinus","Difference(rho_measured-rho_true) for a leading cluster with eta<0 and rho calculated",240,-3.,3.);
  TH2* hScatter_Delta_eta_Delta_rho_multicluster=new TH2F("hScatter_Delta_eta_Delta_rho_multicluster","Scatter plot of Delta_eta versus Delta_rho(cm) for a leading multicluster",160,-0.04,0.04,240,-3.,3.); 
  TH2* hScatter_Delta_eta_Delta_rho_lowerGenEta=new TH2F("hScatter_Delta_eta_Delta_rho_lowerGenEta","Scatter plot of Delta_eta versus Delta_rho for photons with abs(eta)<2.2",240,-0.06,0.06,400,-5.,5.); 
  TH2* hScatter_Delta_eta_Delta_rho_higherGenEta=new TH2F("hScatter_Delta_eta_Delta_rho_higherGenEta","Scatter plot of Delta_eta versus Delta_rho for photons with abs(eta)>2.2",240,-0.06,0.06,400,-5.,5.); 
  TH2* hScatter_Delta_theta_Delta_rho_lowerGenEta=new TH2F("hScatter_Delta_theta_Delta_rho_lowerGenEta","Scatter plot of Delta_theta versus Delta_rho for photons with abs(eta)<2.2",240,-0.06,0.06,400,-5.,5.); 
  TH2* hScatter_Delta_theta_Delta_rho_higherGenEta=new TH2F("hScatter_Delta_theta_Delta_rho_higherGenEta","Scatter plot of Delta_theta versus Delta_rho for photons with abs(eta)>2.2",240,-0.06,0.06,400,-5.,5.);
  TProfile* hEneProfileVsLayers = new TProfile("hEneProfileVsLayers", "Energy deposition as a function of layers", 50, 0., 50.,0.,1.);
  TProfile* hNumRechitsProfileVsLayers = new TProfile("hNumRechitsProfileVsLayers", "Profile of number of rechits for a leading 2D cluster as a function of layers", 50, 0., 50.,0.,1000.);

  char histRho[500];
  char histRhoPhi[500];
  char histEta[500];
  char histPhi[500];
  char histNumClusterNoThreshold[500];
  char histNumCluster[500];
  char histClusterEnergy[500];
  char histClusterPt[500];
  char EneProfForLayer[500];
  int EneClusterIndex=1;

  sprintf(histRho,"hDelta_rho_E%i",EneClusterIndex);
  sprintf(histRhoPhi,"hDelta_rhophi_E%i",EneClusterIndex);
  sprintf(histEta,"hDelta_eta_E%i",EneClusterIndex);
  sprintf(histPhi,"hDelta_phi_E%i",EneClusterIndex);
  sprintf(histNumCluster,"hNumMultiClusters");
  sprintf(histNumClusterNoThreshold,"hNumMultiClustersNoThreshold");
  sprintf(histClusterEnergy,"hLeadingMultiClusterEnergy");
  sprintf(histClusterPt,"hLeadingMultiClusterPt");
  TH1* hDelta_rho_multicluster= new TH1F(histRho,"Difference(rho_measured-rho_true) for a leading energy multicluster(cm)",240,-1.0,1.0);
  TH1* hDelta_rhophi_multicluster = new TH1F(histRhoPhi,"Difference[rho_true*(phi_measured-phi_true)] for a leading energy multicluster",240,-1.0,1.0);
  TH1* hDelta_eta_multicluster= new TH1F(histEta,"Difference(eta_measured-eta_true) for a leading energy multicluster",320,-0.02,0.02);
  TH1* hDelta_phi_multicluster= new TH1F(histPhi,"Difference(phi_measured-phi_true) for a leading energy multicluster",320,-0.02,0.02);
  TH1* hNumMultiClusters= new TH1F(histNumCluster,"Number of multiclusters with energy above 1GeV",15,0,15);
  TH1* hNumMultiClusters_NoThreshold= new TH1F(histNumClusterNoThreshold,"Number of multiclusters",15,0,15);
  TH1* hLeadingMultiClusterEnergy= new TH1F(histClusterEnergy,"Leading multicluster Energies",102,-5,505);
  TH1* hLeadingMultiClusterPt= new TH1F(histClusterPt,"Leading multicluster Pt",300,20,50);
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
    vector<int> *rechit_cluster2d=0;
    vector<float> *rechit_eta=0;
    vector<int> *rechit_layer=0;
    vector<float> *multiclus_energy=0;
    vector<float> *multiclus_phi=0;
    vector<float> *multiclus_eta=0;
    vector<float> *multiclus_pt=0;
    vector<float> *multiclus_slopeX=0;
    vector<float> *multiclus_slopeY=0;
    vector<float> *multiclus_z=0;
    vector<float> *rechit_energy=0;
    
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
    theTree->SetBranchAddress("rechit_cluster2d",&rechit_cluster2d);
    theTree->SetBranchAddress("rechit_energy",&rechit_energy);
    theTree->SetBranchAddress("rechit_eta",&rechit_eta);
    theTree->SetBranchAddress("rechit_layer",&rechit_layer);
    theTree->SetBranchAddress("genpart_pt", &genpart_pt);
    theTree->SetBranchAddress("genpart_eta", &genpart_eta);
    theTree->SetBranchAddress("genpart_phi", &genpart_phi);
    theTree->SetBranchAddress("genpart_dvx", &genpart_dvx);
    theTree->SetBranchAddress("genpart_dvy", &genpart_dvy);
    theTree->SetBranchAddress("genpart_reachedEE", &genpart_reachedEE);
    theTree->SetBranchAddress("genpart_pid", &genpart_pid);

    // loop over events
    int nEntries = theTree->GetEntries();
    //    std::cout<<" nEntries "<<nEntries<<endl;
    for( int evtIndex = 0; evtIndex < nEntries; evtIndex++ )
      {
	TotNumOfEvents++;
	// get values
	theTree->GetEntry( evtIndex );
	std::vector<int> leadingcluster_index;
	std::vector<float> leadingcluster_x;
	std::vector<float> leadingcluster_y;
	std::vector<float> leadingcluster_phi;
	std::vector<float> leadingcluster_eta;
	int PhotonsOfInterest=0;
	int NumMultiClusters=0;
	int NumMultiClustersNoThreshold=0;
	
	//loop over gen particles to check whether a gen particle is a photon having pt=35 GeV and is non interacted in a tracker.
	for( int genIndex = 0; genIndex < genpart_pt->size(); genIndex++ ) {
	  float photon_Eta=0;
	  float photon_Phi=0;
	  float photon_X=0;
	  float photon_Y=0;
	  if(genpart_pid->at(genIndex)==22 && genpart_reachedEE->at(genIndex)==1 && genpart_pt->at(genIndex)==35.){
	    std::vector<float> leadingcluster_energy;
	    std::vector<int> leadingcluster_NumRechits;
	    float CumulativeEnergy=0;
	    PhotonsOfInterest++;
	    TotPhotonsOfInterest++;
	    photon_Eta=genpart_eta->at(genIndex);
	    photon_Phi=genpart_phi->at(genIndex);
	    photon_X=genpart_dvx->at(genIndex);
	    photon_Y=genpart_dvy->at(genIndex);
	    int jx=0;
	    //Determines number of multiclusters above 0GeV and 1GeV.
	    int startIndex=-1;
	    float MaxFound=9999999.;
	    for(int multiclusEneIndex = 0; multiclusEneIndex < multiclus_energy->size(); multiclusEneIndex++ ) {
	      if(multiclus_energy->at(multiclusEneIndex)>0.0 && ((photon_Eta>0 && multiclus_eta->at(multiclusEneIndex)>0) || (photon_Eta<0 && multiclus_eta->at(multiclusEneIndex)<0)))
		NumMultiClustersNoThreshold++;
	      if(multiclus_energy->at(multiclusEneIndex)>1.0 && ((photon_Eta>0 && multiclus_eta->at(multiclusEneIndex)>0) || (photon_Eta<0 && multiclus_eta->at(multiclusEneIndex)<0))){
		NumMultiClusters++;
	      }
	    }
	    
	    //Finds the first leading energy multicluster with energy above 1 GeV
	    for(int clusOrderIndex=0;clusOrderIndex<1;clusOrderIndex++){   
	      double rhoT=0;
	      double rhoM=0;
	      double thetaT=0;
	      double thetaM=0;
	      double rhoM_rhoT=0;
	      double thetaM_thetaT=0;
	      double rhoC=0;
	      double rhoC_rhoT=0;
	      double rhophiM_rhophiT=0;
	      double etaM_etaT=0;
	      double phiM_phiT=0;
	      float Max_multiclus_energy=0;
	      float Max_multiclus_x=0;
	      float Max_multiclus_y=0;
	      float Max_multiclus_z=0;
	      float Max_multiclus_rhoC=0;
	      float Max_multiclus_thetaC=0;
	      float Max_multiclus_etaC=0;	
	      float Max_multiclus_phi=0;
	      float Max_multiclus_pt=0;
	      float Max_multiclus_eta=0;
	      for(int multiclusEneIndex = 0; multiclusEneIndex < multiclus_energy->size(); multiclusEneIndex++ ) {
		if(multiclus_energy->at(multiclusEneIndex)>1.0 && ((photon_Eta>0 && multiclus_eta->at(multiclusEneIndex)>0) || (photon_Eta<0 && multiclus_eta->at(multiclusEneIndex)<0))){
		  if(multiclusEneIndex!=startIndex && multiclus_energy->at(multiclusEneIndex)<=MaxFound){
		    if(multiclus_energy->at(multiclusEneIndex)>Max_multiclus_energy){
		      Max_multiclus_energy=multiclus_energy->at(multiclusEneIndex);
		      startIndex=multiclusEneIndex;
		      Max_multiclus_z=multiclus_z->at(multiclusEneIndex);
		      Max_multiclus_x=multiclus_slopeX->at(multiclusEneIndex);
		      Max_multiclus_y=multiclus_slopeY->at(multiclusEneIndex);
		      Max_multiclus_phi=multiclus_phi->at(multiclusEneIndex);
		      Max_multiclus_eta=multiclus_eta->at(multiclusEneIndex);
		      Max_multiclus_pt=multiclus_pt->at(multiclusEneIndex);
		    }
		  }
		}
	      } //end of loop over multiclusEneIndex
	      MaxFound=Max_multiclus_energy; 
	      etaM_etaT=abs(Max_multiclus_eta)-abs(photon_Eta);
	      phiM_phiT=Max_multiclus_phi-photon_Phi;
	      rhoT=Max_multiclus_z*tan(2*atan(exp(-photon_Eta))); //
	      rhoM=sqrt(Max_multiclus_x*Max_multiclus_x+Max_multiclus_y*Max_multiclus_y);
	      rhoM_rhoT=rhoM-rhoT;
	      rhophiM_rhophiT=rhoT*(Max_multiclus_phi-photon_Phi);
	      
	      //Finds leading energy 2D cluster per layer which belongs to leading energy multicluster having energy above 1GeV  and fills into a vector leadingcluster_energy.
	      if(Max_multiclus_energy>1.0){
		hMultiPt->Fill(Max_multiclus_pt);
		for(int LIndex=0;LIndex<50;LIndex++){ 
		  int start2DIndex=-1;
		  float Max2DFound=9999999.;
		  float Max_cluster2D_energy=0;
		  for(int clusEneIndex = 0; clusEneIndex < cluster2d_energy->size(); clusEneIndex++ ) {
		    if((photon_Eta>0 && cluster2d_eta->at(clusEneIndex)>0 && cluster2d_multicluster->at(clusEneIndex)==startIndex  && cluster2d_layer->at(clusEneIndex)==(LIndex+1)) || (photon_Eta<0 && cluster2d_eta->at(clusEneIndex)<0 && cluster2d_multicluster->at(clusEneIndex)==startIndex && cluster2d_layer->at(clusEneIndex)==(LIndex+1))){
		      if(clusEneIndex!=start2DIndex && cluster2d_energy->at(clusEneIndex)<=Max2DFound){
			if(cluster2d_energy->at(clusEneIndex)>Max_cluster2D_energy){
			  Max_cluster2D_energy=cluster2d_energy->at(clusEneIndex);
			  start2DIndex=clusEneIndex;
			}
		      }
		    }
		  } //end of loop over clusEneIndex
		  Max2DFound=Max_cluster2D_energy;
		  leadingcluster_energy.push_back(Max2DFound);
		  CumulativeEnergy +=Max2DFound;	  
		  int RecHitsFromCluster=0;
		  //If leading energy 2D cluster is found for a given layer, finds number of rechits belonging to that cluster and fills into a vector leadingcluster_NumRechits
		  if(Max_cluster2D_energy>0.){
		    for(int hitIndex = 0; hitIndex < rechit_energy->size(); hitIndex++){
		      if((photon_Eta>0 && rechit_eta->at(hitIndex)>0 && rechit_layer->at(hitIndex)==(LIndex+1) && rechit_cluster2d->at(hitIndex)==start2DIndex) || (photon_Eta<0 && rechit_eta->at(hitIndex)<0 && rechit_layer->at(hitIndex)==(LIndex+1) && rechit_cluster2d->at(hitIndex)==start2DIndex)){
			RecHitsFromCluster++;
		      }
		    }
		  }
		  leadingcluster_NumRechits.push_back(RecHitsFromCluster);
		}
		//Whenever 2D clusters(E>0) are found associated to a leading multicluster, plot profile of energy fractions and number of rechits belonging to leading cluster as a function of layers
		if(CumulativeEnergy!=0.){		 
		  for(int LIndex=0;LIndex<50;LIndex++){
		    hEneProfileVsLayers->Fill(LIndex+1,leadingcluster_energy[LIndex]/CumulativeEnergy,1);
		    hNumRechitsProfileVsLayers->Fill(LIndex+1,leadingcluster_NumRechits[LIndex],1);
		  }
		}
		hScatter_Delta_eta_Delta_rho_multicluster->Fill(etaM_etaT,rhoM_rhoT);
		hDelta_eta_multicluster->Fill(etaM_etaT);
		hDelta_phi_multicluster->Fill(phiM_phiT);
		hDelta_rho_multicluster->Fill(rhoM_rhoT);
		hDelta_rhophi_multicluster->Fill(rhophiM_rhophiT);
		hLeadingMultiClusterEnergy->Fill(MaxFound);
		hLeadingMultiClusterPt->Fill(Max_multiclus_pt);
		
	      }
	      else{
		NoMatchingMultiClusToPhoton++;
	      }
	      jx++;
	    }
	  }
	}
	
	hNonInteractedPhotons->Fill(PhotonsOfInterest);
	hNumMultiClusters->Fill(NumMultiClusters);
	hNumMultiClusters_NoThreshold->Fill(NumMultiClustersNoThreshold);
	if(PhotonsOfInterest!=NumMultiClusters){
	  NonMatch_Photon_Multicluster++;
	  if(PhotonsOfInterest>NumMultiClusters)
	    NonMatch_MorePhoton_LessMulticluster++;
	  if(PhotonsOfInterest<NumMultiClusters)
	    NonMatch_LessPhoton_MoreMulticluster++;
	  if(PhotonsOfInterest>NumMultiClusters && NumMultiClusters==0)
	    NonMatch_MorePhoton_ZeroMulticluster++;
	}
	else
	  Match_Photon_Multicluster++;
      }
  }
  
  std::cout<<" TotNumOfEvents "<<TotNumOfEvents<<" TotNoninteractedPhotons "<<TotPhotonsOfInterest<<" NoMatchingMultiClusToPhoton "<<NoMatchingMultiClusToPhoton<<" NonMatch_Photon_Multicluster "<<NonMatch_Photon_Multicluster<<" NonMatch_LessPhoton_MoreMulticluster "<<NonMatch_LessPhoton_MoreMulticluster<<" NonMatch_MorePhoton_LessMulticluster "<<NonMatch_MorePhoton_LessMulticluster<<" NonMatch_MorePhoton_ZeroMulticluster "<<NonMatch_MorePhoton_ZeroMulticluster<<" Match_Photon_Multicluster "<<Match_Photon_Multicluster<<endl;
  // draw histograms
  TFile *ouputFile = new TFile("BasicPlots_old_Jul11_HGCAL_4.root","RECREATE");
  hNonInteractedPhotons->Write();
  hScatter_Delta_eta_Delta_rho_multicluster->Write();
  hEneProfileVsLayers->Write();
  hMultiPt->Write();
  hNumRechitsProfileVsLayers->Write();
  hDelta_rho_multicluster->Write();
  hDelta_rhophi_multicluster->Write();
  hDelta_eta_multicluster->Write();
  hDelta_phi_multicluster->Write();
  hLeadingMultiClusterEnergy->Write();
  hLeadingMultiClusterPt->Write();
  hNumMultiClusters->Write();
  hNumMultiClusters_NoThreshold->Write();
  ouputFile->Close();
  
}

