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

  //  TH1* hDelta_rho_etaPlus= new TH1F("hDelta_rho_etaPlus","Difference(rho_measured-rho_true) for a leading cluster with eta>0 ",240,-6.,6.);
  // TH1* hDelta_rho_etaMinus= new TH1F("hDelta_rho_etaMinus","Difference(rho_measured-rho_true) for a leading cluster with eta<0 ",240,-6.,6.);

  TH2* hScatter_Delta_eta_Delta_rho=new TH2F("hScatter_Delta_eta_Delta_rho","Scatter plot of Delta_eta versus Delta_rho",160,-0.04,0.04,240,-3.,3.); 
  TH2* hScatter_Delta_eta_Delta_rho_lowerGenEta=new TH2F("hScatter_Delta_eta_Delta_rho_lowerGenEta","Scatter plot of Delta_eta versus Delta_rho for photons with abs(eta)<2.2",240,-0.06,0.06,400,-5.,5.); 
  TH2* hScatter_Delta_eta_Delta_rho_higherGenEta=new TH2F("hScatter_Delta_eta_Delta_rho_higherGenEta","Scatter plot of Delta_eta versus Delta_rho for photons with abs(eta)>2.2",240,-0.06,0.06,400,-5.,5.); 
  TH2* hScatter_Delta_theta_Delta_rho_lowerGenEta=new TH2F("hScatter_Delta_theta_Delta_rho_lowerGenEta","Scatter plot of Delta_theta versus Delta_rho for photons with abs(eta)<2.2",240,-0.06,0.06,400,-5.,5.); 
  TH2* hScatter_Delta_theta_Delta_rho_higherGenEta=new TH2F("hScatter_Delta_theta_Delta_rho_higherGenEta","Scatter plot of Delta_theta versus Delta_rho for photons with abs(eta)>2.2",240,-0.06,0.06,400,-5.,5.); 

  char histRho[500];
  char histRhoPhi[500];
  char histEta[500];
  char histPhi[500];
  char histNumCluster[500];
  char histClusterEnergy[500];
  int layerIndex=6;
  int EneClusterIndex=1;
  sprintf(histRho,"hDelta_rho_L%i_E%i",layerIndex,EneClusterIndex);
  sprintf(histRhoPhi,"hDelta_rhophi_L%i_E%i",layerIndex,EneClusterIndex);
  sprintf(histEta,"hDelta_eta_L%i_E%i",layerIndex,EneClusterIndex);
  sprintf(histPhi,"hDelta_phi_L%i_E%i",layerIndex,EneClusterIndex);
  sprintf(histClusterEnergy,"hLeadingClusterEnergy_L%i_E%i",layerIndex,EneClusterIndex);
  TH1* hDelta_rho= new TH1F(histRho,"Difference(rho_measured-rho_true)",240,-3.,3.);
  TH1* hDelta_rhophi = new TH1F(histRhoPhi,"Difference[rho_true*(phi_measured-phi_true)]",240,-3.,3.);
  TH1* hDelta_eta= new TH1F(histEta,"Difference(eta_measured-eta_true)",160,-0.04,0.04);
  TH1* hDelta_phi= new TH1F(histPhi,"Difference(phi_measured-phi_true)",160,-0.04,0.04);
  TH1* hLeadingClusterEnergy= new TH1F(histClusterEnergy,"2D cluster Energies",120,-5,35);
  sprintf(histNumCluster,"hNumClusters_L%i",layerIndex);
  TH1* hNumClusters= new TH1F(histNumCluster,"Number of 2D clusters",15,0,15);
  

  //  sprintf(location,/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_SinglePhoton_35GeV_20170523/NTUP/)
  int NumInputFiles=25;
  int Num_DeltaRho_Plus=0;
  int Num_DeltaRho_Minus=0;
  int Num_DeltaEta_Plus=0;
  int Num_DeltaEta_Minus=0;
  for(int j=1;j<=NumInputFiles;j++){
    sprintf(fileName,"/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_SingleParticlePtGun_PID22_SinglePhoton_Pt35GeV_20170623/NTUP/partGun_PDGid22_x80_Pt35.0To35.0_NTUP_%i.root",j);
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
	    bool ConsiderCluster=false;	    
	    int jx=0;
	    //Determines 3 leading energy clusters for layers 5-10(which most probably corresponds to maximum of a shower) and finds deviations with respect to gen particle position in terms of delta_eta, delta_phi, delta_rho and delta_rhophi.
	    for(int LIndex=5;LIndex<6;LIndex++){
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
	      
	      for(int clusOrderIndex=0;clusOrderIndex<1;clusOrderIndex++){   
		//std::cout<<" LIndex "<<LIndex<<" clusOrderIndex "<<clusOrderIndex<<endl;
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
		float Max_cluster2D_energy=0;
		float Max_cluster2D_x=0;
		float Max_cluster2D_y=0;
		float Max_cluster2D_z=0;
		float Max_cluster2D_rhoC=0;
		float Max_cluster2D_thetaC=0;
		float Max_cluster2D_etaC=0;	
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
		Max_cluster2D_rhoC=sqrt(Max_cluster2D_x*Max_cluster2D_x+Max_cluster2D_y*Max_cluster2D_y); //Max_cluster2D_rhoC returned this way is always positive
		Max_cluster2D_thetaC=atan(Max_cluster2D_rhoC/Max_cluster2D_z); //if z is -ve, theta is negative,eta has to be negative 

		Max_cluster2D_etaC=-log(tan(fabs(Max_cluster2D_thetaC/2)));
		if(Max_cluster2D_thetaC<0)
		  Max_cluster2D_etaC *=-1;
		/*
		if(photon_Eta>0 && Max_cluster2D_z!=0)
		  std::cout<<" evtIndex "<<evtIndex<<" genIndex "<<genIndex<<" photon + "<<" z_cluster "<<Max_cluster2D_z<<" theta_cluster "<<Max_cluster2D_thetaC<<"  eta_clusterC "<<Max_cluster2D_etaC<<"  eta_clusterM "<<Max_cluster2D_eta<<" diff_etaC_etaM "<< fabs(Max_cluster2D_etaC)-fabs(Max_cluster2D_eta)<<endl;
		else if(photon_Eta<0 && Max_cluster2D_z!=0)
		  std::cout<<" evtIndex "<<evtIndex<<" genIndex "<<genIndex<<" photon - "<<" z_cluster "<<Max_cluster2D_z<<" theta_cluster "<<Max_cluster2D_thetaC<<"  eta_clusterC "<<Max_cluster2D_etaC<<"  eta_clusterM "<<Max_cluster2D_eta<<" diff_etaC_etaM "<< fabs(Max_cluster2D_etaC)-fabs(Max_cluster2D_eta)<<endl;
*/
		etaM_etaT=abs(Max_cluster2D_eta)-abs(photon_Eta);
		/*
		if(photon_Eta>0 && Max_cluster2D_z!=0 && abs(etaM_etaT)>0.015 && abs(etaM_etaT)<0.035)
		  std::cout<<" evtIndex "<<evtIndex<<" genIndex "<<genIndex<<" photon + "<<" eta_cluster "<<Max_cluster2D_eta<<" eta_cluster_abs "<<abs(Max_cluster2D_eta)<<" eta_photon "<<photon_Eta<<" eta_photon_abs "<<abs(photon_Eta)<<" etaM_etaT "<<etaM_etaT<<endl;
		if(photon_Eta<0 && Max_cluster2D_z!=0 && abs(etaM_etaT)>0.015 && abs(etaM_etaT)<0.035)
		  std::cout<<" evtIndex "<<evtIndex<<" genIndex "<<genIndex<<" photon - "<<" eta_cluster "<<Max_cluster2D_eta<<" eta_cluster_abs "<<abs(Max_cluster2D_eta)<<" eta_photon "<<photon_Eta<<" eta_photon_abs "<<abs(photon_Eta)<<" etaM_etaT "<<etaM_etaT<<endl;
*/			
		phiM_phiT=Max_cluster2D_phi-photon_Phi;
		rhoT=Max_cluster2D_z*tan(2*atan(exp(-photon_Eta))); //
		thetaT=tan(2*atan(exp(-photon_Eta)));
		thetaM=tan(2*atan(exp(-Max_cluster2D_eta)));
		thetaM_thetaT=abs(thetaM)-abs(thetaT);
	//here also it is found that for both +ve and -ve photon_Eta, rhoT is always true.
		/*
		if(photon_Eta>0 && Max_cluster2D_z!=0)
		  std::cout<<" evtIndex "<<evtIndex<<" genIndex "<<genIndex<<" photon + "<<" eta_photon "<<photon_Eta<<" tan "<<tan(2*atan(exp(-photon_Eta)))<<" z_cluster "<<Max_cluster2D_z<<" rho_photon "<<rhoT<<endl;
		if(photon_Eta<0 && Max_cluster2D_z!=0)
		  std::cout<<" evtIndex "<<evtIndex<<" genIndex "<<genIndex<<" photon - "<<" eta_photon "<<photon_Eta<<" tan "<<tan(2*atan(exp(-photon_Eta)))<<" z_cluster "<<Max_cluster2D_z<<" rho_photon "<<rhoT<<endl;
*/
		rhoC=Max_cluster2D_z*tan(2*atan(exp(-Max_cluster2D_eta)));
		rhoM=sqrt(Max_cluster2D_x*Max_cluster2D_x+Max_cluster2D_y*Max_cluster2D_y);
		rhoC_rhoT=rhoC-rhoT;

		//		rhoM_rhoT=sqrt(Max_cluster2D_x*Max_cluster2D_x+Max_cluster2D_y*Max_cluster2D_y)-sqrt(photon_X*photon_X+photon_Y*photon_Y);
		rhoM_rhoT=rhoM-rhoT;

		//		rhophiM_rhophiT=sqrt(photon_X*photon_X+photon_Y*photon_Y)*(Max_cluster2D_phi-photon_Phi);
		rhophiM_rhophiT=rhoT*(Max_cluster2D_phi-photon_Phi);

		//	      	  std::cout<<" LIndex "<<LIndex+1<<" photon_eta "<<photon_Eta<<" NumClusters "<<NumClusters<<" clusOrderIndex "<<clusOrderIndex<<" MaxFound "<<MaxFound<<" rhoM_rhoT "<<rhoM_rhoT<<endl;
		if(Max_cluster2D_energy!=0){
		  if(etaM_etaT>0)
		    Num_DeltaEta_Plus++;
		  if(etaM_etaT<0)
		    Num_DeltaEta_Minus++;
		  if(rhoM_rhoT>0)
		    Num_DeltaRho_Plus++;
		  if(rhoM_rhoT<0)
		    Num_DeltaRho_Minus++;
		  
		  if(Max_cluster2D_eta>0)
		    hDelta_rhoM_etaPlus->Fill(rhoM_rhoT);
		  if(Max_cluster2D_eta<0)		
		    hDelta_rhoM_etaMinus->Fill(rhoM_rhoT);
		  if(Max_cluster2D_eta>0)
		    hDelta_rhoC_etaPlus->Fill(rhoC_rhoT);
		  if(Max_cluster2D_eta<0)		
		    hDelta_rhoC_etaMinus->Fill(rhoC_rhoT);
		  hScatter_Delta_eta_Delta_rho->Fill(etaM_etaT,rhoM_rhoT);
		  if(abs(photon_Eta)<2.2){
		    hScatter_Delta_eta_Delta_rho_lowerGenEta->Fill(etaM_etaT,rhoM_rhoT);
		    hScatter_Delta_theta_Delta_rho_lowerGenEta->Fill(thetaM_thetaT,rhoM_rhoT);
		  }
		  if(abs(photon_Eta)>2.2){
		    hScatter_Delta_eta_Delta_rho_higherGenEta->Fill(etaM_etaT,rhoM_rhoT);
		    hScatter_Delta_theta_Delta_rho_higherGenEta->Fill(thetaM_thetaT,rhoM_rhoT);
		  }
		  hDelta_eta->Fill(etaM_etaT);
		  hDelta_phi->Fill(phiM_phiT);
		  hDelta_rho->Fill(rhoM_rhoT);
		  hDelta_rhophi->Fill(rhophiM_rhophiT);
		  hLeadingClusterEnergy->Fill(MaxFound);
		  hDiff_etaC_etaM->Fill(abs(Max_cluster2D_etaC)-abs(Max_cluster2D_eta));		
		}
		//if(jx%3==0)
		hNumClusters->Fill(NumClusters);
		jx++;
	      }
	    }
	  }
	}
	hNonInteractedPhotons->Fill(PhotonsOfInterest);
      }
  }
  std::cout<<" Num_DeltaEta_Plus "<<Num_DeltaEta_Plus<<" Num_DeltaEta_Minus "<<Num_DeltaEta_Minus<<" Num_DeltaRho_Plus "<<Num_DeltaRho_Plus<<" Num_DeltaRho_Minus "<<Num_DeltaRho_Minus<<endl;
  hDeltaEta_Plus->Fill(Num_DeltaEta_Plus);
  hDeltaEta_Minus->Fill(Num_DeltaEta_Minus);
  hDeltaRho_Plus->Fill(Num_DeltaRho_Plus);
  hDeltaRho_Minus->Fill(Num_DeltaRho_Minus);

  // draw histograms
  TFile *ouputFile = new TFile("BasicPlots.root","RECREATE");
  hMultiPt->Write();
  hDeltaEta_Plus->Write(); 
  hDeltaEta_Minus->Write();
  hDeltaRho_Plus->Write();
  hDeltaRho_Minus->Write();
  hNonInteractedPhotons->Write();
  hDiff_etaC_etaM->Write();
  hDelta_rhoC_etaPlus->Write();
  hDelta_rhoC_etaMinus->Write();
  hDelta_rhoM_etaPlus->Write();
  hDelta_rhoM_etaMinus->Write();
  hScatter_Delta_eta_Delta_rho->Write();
  hScatter_Delta_eta_Delta_rho_lowerGenEta->Write();
  hScatter_Delta_eta_Delta_rho_higherGenEta->Write();
  hScatter_Delta_theta_Delta_rho_lowerGenEta->Write();
  hScatter_Delta_theta_Delta_rho_higherGenEta->Write();

  //  hNumClusterPerLayer->Write();
 
  hDelta_rho->Write();
  hDelta_rhophi->Write();
  hDelta_eta->Write();
  hDelta_phi->Write();
  hLeadingClusterEnergy->Write();
  hNumClusters->Write();

  //    TCanvas *c0 = new TCanvas(histRho,histRho,500,500);   
  //hDelta_rho.at(j)->Draw();
  //TCanvas *c1 = new TCanvas(histRhoPhi,histRhoPhi,500,500);
  //hDelta_rhophi.at(j)->Draw();
  ouputFile->Close();

  //  TCanvas *c0 = new TCanvas();
  //hMultiPt->Draw("hist");
  //c0->Print( "hMultiPt.pdf" );
  //c0->Print( "hMultiPt.png" );
  
}

