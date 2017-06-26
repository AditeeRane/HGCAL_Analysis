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
  TH1* hMultiPt = new TH1F( "hMultiPt", "multicluster pt", 50, 0., 50. );
  TH1* hNonInteractedPhotons = new TH1F( "hNonInteractedPhotons","Number of NonInteractedPhotons",5,0,5);
  std::vector<TH1*> hDelta_rho(15);
  std::vector<TH1*> hDelta_rhophi(15);
  std::vector<TH1*> hDelta_eta(15);
  std::vector<TH1*> hDelta_phi(15);
  std::vector<TH1*> hNumClusters(5);
  std::vector<TH1*> hLeadingClusterEnergy(15);
  char histRho[500];
  char histRhoPhi[500];
  char histEta[500];
  char histPhi[500];
  char histNumCluster[500];
  char histClusterEnergy[500];
  for(int j=0;j<15;j++){
    int layerIndex=int(j/3)+6;
    int EneClusterIndex=(j%3)+1;
    sprintf(histRho,"hDelta_rho_L%i_E%i",layerIndex,EneClusterIndex);
    sprintf(histRhoPhi,"hDelta_rhophi_L%i_E%i",layerIndex,EneClusterIndex);
    sprintf(histEta,"hDelta_eta_L%i_E%i",layerIndex,EneClusterIndex);
    sprintf(histPhi,"hDelta_phi_L%i_E%i",layerIndex,EneClusterIndex);
    sprintf(histClusterEnergy,"hLeadingClusterEnergy_L%i_E%i",layerIndex,EneClusterIndex);
    hDelta_rho.at(j)= new TH1F(histRho,"Difference(rho_measured-rho_true)",60,-10,10);
    hDelta_rhophi.at(j) = new TH1F(histRhoPhi,"Difference[rho_true*(phi_measured-phi_true)]",60,-10,10);
    hDelta_eta.at(j)= new TH1F(histEta,"Difference(eta_measured-eta_true)",40,-0.1,0.1);
    hDelta_phi.at(j)= new TH1F(histPhi,"Difference(phi_measured-phi_true)",40,-0.1,0.1);
    hLeadingClusterEnergy.at(j)= new TH1F(histClusterEnergy,"2D cluster Energies",120,-5,35);
    if(j%3==0){
      sprintf(histNumCluster,"hNumClusters_L%i",layerIndex);
      hNumClusters.at(j/3)= new TH1F(histNumCluster,"Number of 2D clusters",15,0,15);
    }
  }
  //  sprintf(location,/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_SinglePhoton_35GeV_20170523/NTUP/)
  int NumInputFiles=9;
  for(int j=1;j<=NumInputFiles;j++){
    sprintf(fileName,"/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_SinglePhoton_35GeV_20170523/NTUP/partGun_PDGid22_x120_Pt35.0To35.0_NTUP_%i.root",j);
    TFile *inputFile = TFile::Open(fileName);
    std::cout<<" inputFile "<<fileName<<endl;
    TTree *theTree = (TTree*)inputFile->Get("ana/hgc");
    // define quantities in ntuple
    vector<float> *multiclus_pt = 0;
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
    vector<float> *multiclus_energy=0;
    vector<float> *multiclus_phi=0;
   
    theTree->SetBranchAddress("multiclus_pt", &multiclus_pt);
    theTree->SetBranchAddress("genpart_pt", &genpart_pt);
    theTree->SetBranchAddress("genpart_eta", &genpart_eta);
    theTree->SetBranchAddress("genpart_phi", &genpart_phi);
    theTree->SetBranchAddress("genpart_dvx", &genpart_dvx);
    theTree->SetBranchAddress("genpart_dvy", &genpart_dvy);
    theTree->SetBranchAddress("genpart_reachedEE", &genpart_reachedEE);
    theTree->SetBranchAddress("genpart_pid", &genpart_pid);
    theTree->SetBranchAddress("cluster2d_layer", &cluster2d_layer);
    theTree->SetBranchAddress("cluster2d_energy", &cluster2d_energy);
    theTree->SetBranchAddress("cluster2d_x", &cluster2d_x);
    theTree->SetBranchAddress("cluster2d_y", &cluster2d_y);
    theTree->SetBranchAddress("cluster2d_z", &cluster2d_z);
    theTree->SetBranchAddress("cluster2d_phi", &cluster2d_phi);
    theTree->SetBranchAddress("cluster2d_eta", &cluster2d_eta);
    //    theTree->SetBranchAddress("simcluster_pt",&simcluster_pt);
    // theTree->SetBranchAddress("simcluster_eta",&simcluster_eta);

    // loop over events
    uint nEntries = theTree->GetEntries();
    //    std::cout<<" nEntries "<<nEntries<<endl;
    for( int evtIndex = 0; evtIndex < nEntries; evtIndex++ )
      {
	// get values
	theTree->GetEntry( evtIndex );
	//	std::cout<<" evtIndex "<<evtIndex<<" Number of multiclusters "<<multiclus_pt->size()<<endl;
	std::vector<int> leadingcluster_index;
	std::vector<float> leadingcluster_energy;
	std::vector<float> leadingcluster_x;
	std::vector<float> leadingcluster_y;
	std::vector<float> leadingcluster_phi;
	std::vector<float> leadingcluster_eta;
	int PhotonsOfInterest=0;

	/*	for( uint genIndex = 0; genIndex < genpart_pt->size(); genIndex++ ) {
	  if(fabs(genpart_pid->at(genIndex))==11)
	    break;
	}
	break; */

	//	std::cout<<"*******New Event***********"<<" evt "<<evtIndex<<endl;
	//loop over gen particles to check whether a gen particle is a photon having pt=35 GeV and is non interacted in a tracker.
	for( uint genIndex = 0; genIndex < genpart_pt->size(); genIndex++ ) {
	  float photon_Eta=0;
	  float photon_Phi=0;
	  float photon_X=0;
	  float photon_Y=0;
	  if(genpart_pid->at(genIndex)==22 && genpart_reachedEE->at(genIndex)==1 && genpart_pt->at(genIndex)==35.){
	    //	    std::cout<<" evt "<<evtIndex<<" gen particles "<<genpart_pt->size()<<endl;
	    PhotonsOfInterest++;
	    photon_Eta=genpart_eta->at(genIndex);
	    photon_Phi=genpart_phi->at(genIndex);
	    photon_X=genpart_dvx->at(genIndex);
	    photon_Y=genpart_dvy->at(genIndex);
	    if(photon_Eta>0)
	      std::cout<<"photon + "<<endl;
	    else
	      std::cout<<"photon - "<<endl;
	    bool ConsiderCluster=false;	    
	    int jx=0;
	    //Determines 3 leading energy clusters for layers 5-10(which most probably corresponds to maximum of a shower) and finds deviations with respect to gen particle position in terms of delta_eta, delta_phi, delta_rho and delta_rhophi.
	    for(int LIndex=5;LIndex<10;LIndex++){
	      int startIndex=-1;
	      float MaxFound=9999999.;
	      int NumClusters=0;
	      for(int clusEneIndex = 0; clusEneIndex < cluster2d_energy->size(); clusEneIndex++ ) {
		if((photon_Eta>0 && cluster2d_eta->at(clusEneIndex)>0) || (photon_Eta<0 && cluster2d_eta->at(clusEneIndex)<0)){
		  if(cluster2d_layer->at(clusEneIndex)==(LIndex+1)){
		    NumClusters++;
		  }
		}
	      }
	      
	      for(int clusOrderIndex=0;clusOrderIndex<3;clusOrderIndex++){   
		//std::cout<<" LIndex "<<LIndex<<" clusOrderIndex "<<clusOrderIndex<<endl;
		double rhoT=0;
		double rhoM_rhoT=0;
		double rhophiM_rhophiT=0;
		double etaM_etaT=0;
		double phiM_phiT=0;
		float Max_cluster2D_energy=0;
		float Max_cluster2D_x=0;
		float Max_cluster2D_y=0;
		float Max_cluster2D_z=0;
		float Max_cluster2D_phi=0;
		float Max_cluster2D_pt=0;
		float Max_cluster2D_eta=0;
		for(int clusEneIndex = 0; clusEneIndex < cluster2d_energy->size(); clusEneIndex++ ) {
		  if((photon_Eta>0 && cluster2d_eta->at(clusEneIndex)>0) || (photon_Eta<0 && cluster2d_eta->at(clusEneIndex)<0)){
		    if(clusEneIndex!=startIndex && cluster2d_layer->at(clusEneIndex)==(LIndex+1) && cluster2d_energy->at(clusEneIndex)<=MaxFound){
		      if(cluster2d_energy->at(clusEneIndex)>Max_cluster2D_energy){
			Max_cluster2D_energy=cluster2d_energy->at(clusEneIndex);
			startIndex=clusEneIndex;
			Max_cluster2D_x=cluster2d_x->at(clusEneIndex);
			Max_cluster2D_y=cluster2d_y->at(clusEneIndex);
			Max_cluster2D_z=cluster2d_z->at(clusEneIndex);
			Max_cluster2D_phi=cluster2d_phi->at(clusEneIndex);
			Max_cluster2D_eta=cluster2d_eta->at(clusEneIndex);
		      }
		    }
		  }
		} //end of loop over clusEneIndex
		MaxFound=Max_cluster2D_energy; 
		etaM_etaT=Max_cluster2D_eta-photon_Eta;
		phiM_phiT=Max_cluster2D_phi-photon_Phi;
		rhoT=Max_cluster2D_z*tan(2*atan(exp(-photon_Eta)));
		//		rhoM_rhoT=sqrt(Max_cluster2D_x*Max_cluster2D_x+Max_cluster2D_y*Max_cluster2D_y)-sqrt(photon_X*photon_X+photon_Y*photon_Y);
		rhoM_rhoT=sqrt(Max_cluster2D_x*Max_cluster2D_x+Max_cluster2D_y*Max_cluster2D_y)-rhoT;
		//		rhophiM_rhophiT=sqrt(photon_X*photon_X+photon_Y*photon_Y)*(Max_cluster2D_phi-photon_Phi);
		rhophiM_rhophiT=rhoT*(Max_cluster2D_phi-photon_Phi);

		//	      	  std::cout<<" LIndex "<<LIndex+1<<" photon_eta "<<photon_Eta<<" NumClusters "<<NumClusters<<" clusOrderIndex "<<clusOrderIndex<<" MaxFound "<<MaxFound<<" rhoM_rhoT "<<rhoM_rhoT<<endl;
		if(Max_cluster2D_energy!=0){
		  hDelta_eta.at(jx)->Fill(etaM_etaT);
		  hDelta_phi.at(jx)->Fill(phiM_phiT);
		  hDelta_rho.at(jx)->Fill(rhoM_rhoT);
		  hDelta_rhophi.at(jx)->Fill(rhophiM_rhophiT);
		  hLeadingClusterEnergy.at(jx)->Fill(MaxFound);
		}
		if(jx%3==0)
		  hNumClusters.at(jx/3)->Fill(NumClusters);
		jx++;
	      }
	    }
	  }
	}
	hNonInteractedPhotons->Fill(PhotonsOfInterest);
      }
  }
  /*
	std::cout<<" evt "<<evtIndex<<" nPhotons "<<PhotonsOfInterest<<endl;
	    //	std::cout<<" photon_plusEta "<<photon_plusEta<<" photon_minusEta "<<photon_minusEta<<" photon_plusPhi "<<photon_plusPhi<<" photon_minusPhi "<<photon_minusPhi<<endl;
	    //Determine highest energy clusters corresponding to layer_1

	if(photon_minusEta!=0 && photon_minusX!=0 && photon_minusY!=0 && photon_plusEta!=0 && photon_plusX!=0 && photon_plusY!=0){	 
	  int jx=0;
	  for(int LIndex=5;LIndex<10;LIndex++){
	    float MaxFound=9999999.;
	    int NumClusters=0;
	   
	    for(int clusEneIndex = 0; clusEneIndex < cluster2d_energy->size(); clusEneIndex++ ) {
	      if(photon_Eta>0 && cluster2d_eta->at(clusEneIndex)>0 || photon_Eta<0 && cluster2d_eta->at(clusEneIndex)<0){
		if(cluster2d_layer->at(clusEneIndex)==(LIndex+1))
		  NumClusters++;
	      }
	    }
	    //	    std::cout<<" evt "<<evtIndex<<" LIndex "<<LIndex<<" NumClusters "<<NumClusters<<endl;
	    for(int clusOrderIndex=0;clusOrderIndex<3;clusOrderIndex++){   
	      //std::cout<<" LIndex "<<LIndex<<" clusOrderIndex "<<clusOrderIndex<<endl;
	      double rhoM_rhoT=0;
	      double rhophiM_rhophiT=0;
	      double etaM_etaT=0;
	      double phiM_phiT=0;
	      float Max_cluster2D_energy=0;
	      float Max_cluster2D_x=0;
	      float Max_cluster2D_y=0;
	      float Max_cluster2D_phi=0;
	      float Max_cluster2D_pt=0;
	      float Max_cluster2D_eta=0;
	      for(int clusEneIndex = 0; clusEneIndex < cluster2d_energy->size(); clusEneIndex++ ) {
		if(photon_Eta>0 && cluster2d_eta->at(clusEneIndex)>0 || photon_Eta<0 && cluster2d_eta->at(clusEneIndex)<0){
		  if(clusEneIndex!=startIndex && cluster2d_layer->at(clusEneIndex)==(LIndex+1) && cluster2d_energy->at(clusEneIndex)<=MaxFound){
		  if(cluster2d_energy->at(clusEneIndex)>Max_cluster2D_energy){
		    Max_cluster2D_energy=cluster2d_energy->at(clusEneIndex);
		    startIndex=clusEneIndex;
		    Max_cluster2D_x=cluster2d_x->at(clusEneIndex);
		    Max_cluster2D_y=cluster2d_y->at(clusEneIndex);
		    Max_cluster2D_phi=cluster2d_phi->at(clusEneIndex);
		    Max_cluster2D_eta=cluster2d_eta->at(clusEneIndex);
		  }
		  }
		}
	      } //end of loop over clusEneIndex
	      MaxFound=Max_cluster2D_energy; 
	      etaM_etaT=Max_cluster2D_eta-photon_Eta;
	      phiM_phiT=Max_cluster2D_phi-photon_Phi;
	      rhoM_rhoT=sqrt(Max_cluster2D_x*Max_cluster2D_x+Max_cluster2D_y*Max_cluster2D_y)-sqrt(photon_X*photon_X+photon_Y*photon_Y);
	      rhophiM_rhophiT=sqrt(photon_X*photon_X+photon_Y*photon_Y)*(Max_cluster2D_phi-photon_Phi);
	      if(rhoM_rhoT!=0 && rhophiM_rhophiT!=0 && etaM_etaT!=0 && phiM_phiT!=0){
		hDelta_eta.at(jx)->Fill(etaM_etaT);
		hDelta_phi.at(jx)->Fill(phiM_phiT);
		hDelta_rho.at(jx)->Fill(rhoM_rhoT);
		hDelta_rhophi.at(jx)->Fill(rhophiM_rhophiT);
	      }
	      if(jx%3==0)
		hNumClusters.at(jx/3)->Fill(NumClusters);
	      h3LeadingClusterEnergy.at(jx)->Fill(MaxFound);
	      jx++;
	    }
	  }
	}
	  
      }
  }
*/	
  // draw histograms
  TFile *ouputFile = new TFile("BasicPlots.root","RECREATE");
  hMultiPt->Write();
  hNonInteractedPhotons->Write();
  //  hNumClusterPerLayer->Write();
  for(int j=0;j<15;j++){
    int layerIndex=int(j/3)+6;
    int EneClusterIndex=(j%3)+1;
    hDelta_rho.at(j)->Write();
    hDelta_rhophi.at(j)->Write();
    hDelta_eta.at(j)->Write();
    hDelta_phi.at(j)->Write();
    hLeadingClusterEnergy.at(j)->Write();
    if(j%3==0){
      hNumClusters.at(j/3)->Write();
    }
    //    TCanvas *c0 = new TCanvas(histRho,histRho,500,500);   
    //hDelta_rho.at(j)->Draw();
    //TCanvas *c1 = new TCanvas(histRhoPhi,histRhoPhi,500,500);
    //hDelta_rhophi.at(j)->Draw();
  }
  ouputFile->Close();
  
  //  TCanvas *c0 = new TCanvas();
  //hMultiPt->Draw("hist");
  //c0->Print( "hMultiPt.pdf" );
  //c0->Print( "hMultiPt.png" );
  
}

