// main02.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the pT_Z spectrum at the Tevatron.

#include "Pythia8/Pythia.h"
using namespace Pythia8; 
#define protonMass  0.93827204621
#define neutronMass 0.93956537921
int main() {

  Pythia* pythia;
  pythia = new Pythia("/home/guest/workspace/Pythia8/xmldoc");
  int idA, idB, pythiaInitializationFlag = 0;
  double pxA , pyA , pzA ,pEA;
  double pxB , pyB , pzB ,pEB;
  char ch;
  int motherLine1 ;
  int  idDaughter1, idDaughter2,lineDaughter1, lineDaughter12, lineDaughter2,lineDaughter22  , motherLine2;

  int numOfEvents =100;

  Particle* particleA;
  Particle* particleB;
  Particle* particleC;
  Vec4 vecA(0), vecCoordinate(0);
  double energyLab = 0;
  energyLab = 800.;
  vecA.pz(800.);
  vecA.px(00.);
  vecA.py(00.);
 // vecA.e(800.940);
  particleA = new Particle();
  particleB = new Particle();
  particleC = new Particle();
  particleA->p(vecA);
  particleA->id(-13);
  particleB->id(2112);
  int _verbose = 0;
  bool cutMass = true;
  pxA = particleA->px();//getVector()->px();
  pxB = particleB->px();//getVector()->px();
  pyA = particleA->py();//getVector()->py();
  pyB = particleB->py();//getVector()->py();
  pzA = particleA->pz();//getVector()->pz();
  pzB = particleB->pz();//getVector()->pz();
  idA = particleA->id();//getId();
  //idB = particleB->id();//getId();
  idB = particleB->id();
  pEA = particleA->e();
  pEB = particleB->e();
  if(_verbose)cout<<" id B particle = "<<particleB->id()<<" A particle momentum in initialization "<<particleA->p()<<endl;

//pythia->settings.resetAll();
		if(!_verbose){
  			pythia->readString("Init:showProcesses = off");
  			pythia->readString("Init:showMultipartonInteractions = off");
  			pythia->readString("Init:showChangedSettings = off");
  			pythia->readString("Init:showAllSettings = off");
  			pythia->readString("Init:showChangedParticleData = off");
  			pythia->readString("Init:showChangedResonanceData = off");

  			pythia->readString("Init:showAllParticleData = off");

  			pythia->readString("Main:showChangedSettings = off");
  			pythia->readString("Main:showAllSettings = off");
  			pythia->readString("Main:showChangedParticleData = off");

  			pythia->readString("Main:showChangedResonanceData = off");
  			pythia->readString("Main:showAllParticleData = off");

  			pythia->readString("Stat:showProcessLevel = off");
  		}

  	//	pythia->readString("PhaseSpace:pTHatMin = 3.0");

  	//	pythia->readString("PhaseSpace:mHatMin = 80.");
  	//	pythia->readString("PhaseSpace:mHatMax = 120.");




















		 // Generator.
		 // Pythia pythia;

		  // Allow no substructure in e+- beams: normal for corrected LEP data.
		//  pythia->readString("PDF:lepton = off");
		  // Process selection.
		  //pythia->readString("WeakSingleBoson:ffbar2gmZ = on");
		  // Switch off all Z0 decays and then switch back on those to quarks.
		 // pythia->readString("23:onMode = off");
		//  pythia->readString("23:onIfAny = 1 2 3 4 5");
			pythia->readString("SoftQCD:all = on");
			pythia->readString("HardQCD:all = on");
			//pythia->readString("ProcessLevel:all = on");
		//	pythia->readString("HadronLevel:all = on");
		//	pythia->readString("PartonLevel:all = on");

		  // LEP1 initialization at Z0 mass.
		  pythia->readString("Beams:idA =  -13");
		  pythia->readString("Beams:idB = 2212");
		//  double mZ = pythia->particleData.m0(2212);
		  //cout<<"mZ = "<<mZ<<endl;
		  pythia->settings.parm("Beams:eCM", 100);
		  pythia->init();

		  // Histograms.
		  Hist nCharge("charged multiplicity", 100, -0.5, 99.5);
		  Hist spheri("Sphericity", 100, 0., 1.);
		  Hist linea("Linearity", 100, 0., 1.);
		  Hist thrust("thrust", 100, 0.5, 1.);
		  Hist oblateness("oblateness", 100, 0., 1.);
		  Hist sAxis("cos(theta_Sphericity)", 100, -1., 1.);
		  Hist lAxis("cos(theta_Linearity)", 100, -1., 1.);
		  Hist tAxis("cos(theta_Thrust)", 100, -1., 1.);
		  Hist nLund("Lund jet multiplicity", 40, -0.5, 39.5);
		  Hist nJade("Jade jet multiplicity", 40, -0.5, 39.5);
		  Hist nDurham("Durham jet multiplicity", 40, -0.5, 39.5);
		  Hist eDifLund("Lund e_i - e_{i+1}", 100, -5.,45.);
		  Hist eDifJade("Jade e_i - e_{i+1}", 100, -5.,45.);
		  Hist eDifDurham("Durham e_i - e_{i+1}", 100, -5.,45.);

		  // Set up Sphericity, "Linearity", Thrust and cluster jet analyses.
		  Sphericity sph;
		  Sphericity lin(1.);
		  Thrust thr;
		  ClusterJet lund("Lund");
		  ClusterJet jade("Jade");
		  ClusterJet durham("Durham");

		  // Begin event loop. Generate event. Skip if error. List first few.
		  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
		    if (!pythia->next()) continue;
		    pythia->info.list();
		    // Find and histogram charged multiplicity.
		    int nCh = 0;
		    for (int i = 0; i < pythia->event.size(); ++i)
		      if (pythia->event[i].isFinal() && pythia->event[i].isCharged()) ++nCh;
		    nCharge.fill( nCh );

		    // Find and histogram sphericity.
		    if (sph.analyze( pythia->event )) {
		      if (iEvent < 3) sph.list();
		      spheri.fill( sph.sphericity() );
		      sAxis.fill( sph.eventAxis(1).pz() );
		      double e1 = sph.eigenValue(1);
		      double e2 = sph.eigenValue(2);
		      double e3 = sph.eigenValue(3);
		      if (e2 > e1 || e3 > e2) cout << "eigenvalues out of order: "
		      << e1 << "  " << e2 << "  " << e3 << endl;
		    }

		    // Find and histogram linearized sphericity.
		    if (lin.analyze( pythia->event )) {
		      if (iEvent < 3) lin.list();
		      linea.fill( lin.sphericity() );
		      lAxis.fill( lin.eventAxis(1).pz() );
		      double e1 = lin.eigenValue(1);
		      double e2 = lin.eigenValue(2);
		      double e3 = lin.eigenValue(3);
		      if (e2 > e1 || e3 > e2) cout << "eigenvalues out of order: "
		      << e1 << "  " << e2 << "  " << e3 << endl;
		    }

		    // Find and histogram thrust.
		    if (thr.analyze( pythia->event )) {
		      if (iEvent < 3) thr.list();
		      thrust.fill( thr.thrust() );
		      oblateness.fill( thr.oblateness() );
		      tAxis.fill( thr.eventAxis(1).pz() );
		      if ( abs(thr.eventAxis(1).pAbs() - 1.) > 1e-8
		        || abs(thr.eventAxis(2).pAbs() - 1.) > 1e-8
		        || abs(thr.eventAxis(3).pAbs() - 1.) > 1e-8
		        || abs(thr.eventAxis(1) * thr.eventAxis(2)) > 1e-8
		        || abs(thr.eventAxis(1) * thr.eventAxis(3)) > 1e-8
		        || abs(thr.eventAxis(2) * thr.eventAxis(3)) > 1e-8 ) {
		        cout << " suspicious Thrust eigenvectors " << endl;
		        thr.list();
		      }
		    }

		    // Find and histogram cluster jets: Lund, Jade and Durham distance.
		    if (lund.analyze( pythia->event, 0.01, 0.)) {
		      if (iEvent < 3) lund.list();
		      nLund.fill( lund.size() );
		      for (int j = 0; j < lund.size() - 1; ++j)
		        eDifLund.fill( lund.p(j).e() - lund.p(j+1).e() );
		    }
		    if (jade.analyze( pythia->event, 0.01, 0.)) {
		      if (iEvent < 3) jade.list();
		      nJade.fill( jade.size() );
		      for (int j = 0; j < jade.size() - 1; ++j)
		        eDifJade.fill( jade.p(j).e() - jade.p(j+1).e() );
		    }
		    if (durham.analyze( pythia->event, 0.01, 0.)) {
		      if (iEvent < 3) durham.list();
		      nDurham.fill( durham.size() );
		      for (int j = 0; j < durham.size() - 1; ++j)
		        eDifDurham.fill( durham.p(j).e() - durham.p(j+1).e() );
		    }

		  // End of event loop. Statistics. Output histograms.
		  }
		  pythia->stat();
		  cout << nCharge << spheri << linea << thrust << oblateness
		       << sAxis << lAxis << tAxis
		       << nLund << nJade << nDurham
		       << eDifLund << eDifJade << eDifDurham;

		  // Done.
		  return 0;

























		pythia->readString("SoftQCD:all = on");
		pythia->readString("HardQCD:all = on");
		pythia->readString("WeakSingleBoson:ffbar2gmZ = on");

  			pythia->readString("23: onMode = off");
  			pythia->readString("23:onIfAny = 1 2 3 4 5");
  	/*		pythia->readString("22: onMode = off");
  			pythia->readString("23: onIfAny = 13 -13");
  			pythia->readString("22: onIfAny = 13 -13");
  	  	   // pythia->readString("PhaseSpace:mHatMin = 4.0");
  	  	   // pythia->readString("PhaseSpace:mHatMax = 8.4");
  	  	    pythia->readString("23:mMin = 4.0");
  	  	    pythia->readString("23:mMax = 8.4");

  	*/
			//if(cutMass){
			//	pythia->readString("23:mMin = 2.0");
			//	pythia->readString("23:mMax = 4.1");
			//}


  	//		pythia->settings.mode("tune:pp",3);


  		if(_verbose){
  			cout<<"ida = "<<idA<<endl;
  			cout<<"idb = "<<idB<<endl;
  			cout<<"pxa = "<<pxA<<endl;
  			cout<<"pxb = "<<pxB<<endl;
  			cout<<"pya = "<<pyA<<endl;
  			cout<<"pyb = "<<pyB<<endl;
  			cout<<"pza = "<<pzA<<endl;
  			cout<<"pzb = "<<pzB<<endl;
  		}
  	  pythia->readString("PDF:lepton = off");
  		 pythia->init(13, -13, energyLab, pythia->particleData.m0(2212));
  		pythia->particleData.m0(23);



  		for(int iev = 0; iev < numOfEvents; iev++ ){
  			if (!pythia->next()) continue;
  			pythia->info.list();

  			for (int i = 3; i < pythia->event.size(); ++i){
  				if(pythia->event[i].isHadron()&&pythia->event[i].isFinal())cout<<pythia->event[i].id()<<" "<<pythia->event[i].p();
  			}

  		}






  	return 1;

  		// cin>>ch;
  		//pythiaInitializationFlag = pythia->init(idA,idB,pxA,pyA,pzA,pxB,pyB,pzB);
  		//return pythiaInitializationFlag;

  		ofstream fileDrellYan;
  		fileDrellYan.open("p_p_800GeV_DY_Mcut_10^6.txt");
  		bool flagQuit = false;
  	  // Begin event loop. Generate event. Skip if error. List first one.
  	  for (int iEvent = 0; iEvent < numOfEvents; ++iEvent) {
  		  if(iEvent % 1000 == 0){
  			//cout<<"number of event "<< iEvent <<endl;
  		  }
//cin>>ch;
  	    if (!pythia->next()) continue;
//cin>>ch;
  	  flagQuit = false;
  	    for (int i = 3; i < pythia->event.size(); ++i){
  	    	//if(pythia->event[i].idAbs() == 13 && pythia->event[i].isFinal()){



  	    	if(flagQuit)break;




  	    		if(abs(pythia->event.at(i+3).id()) == 23 ){

  	    								if(_verbose)cout<<"Z0 WAS BORNED !!!"<<endl;
  	    								/*int id = pythia->event.at(pythia->event.at(i_pyEv+3).mother1()).id();
  	    								if(id == 23){
  	    									cout<<"true leptone"<<endl;
  	    								}else{
  	    									id = pythia->event.at(id).mother1();
  	    								}*/
  	    								motherLine1 = i + 3;


  	    								//find line of Z0 which decays on mu+mu-
  	    								do{
  	    									if(_verbose)cout<<" Z0 cycle "<<endl;
  	    									lineDaughter1 = pythia->event.at(motherLine1).daughter1();
  	    									lineDaughter2 = pythia->event.at(motherLine1).daughter2();
  	    									idDaughter1 = pythia->event.at(lineDaughter1).id();
  	    									idDaughter2 = pythia->event.at(lineDaughter2).id();
  	    									if(_verbose){
  	    										cout<<" lineDaughter1 =  "<<lineDaughter1<<endl;
  	    										cout<<" idDaughter1 =  "<<idDaughter1<<endl;
  	    										cout<<" lineDaughter2 =  "<<lineDaughter2<<endl;
  	    										cout<<" idDaughter2 =  "<<idDaughter2<<endl;
  	    									}

  	    									if(idDaughter1 == idDaughter2 && idDaughter1 == 23){
  	    										motherLine1 = lineDaughter1;
  	    										if(_verbose)cout<< "motherLine = "<<motherLine1<<endl;
  	    									}
  	    									if(_verbose)cout<<"pythia->event.at(motherLine1).id() = "<<pythia->event.at(motherLine1).id()<<endl;
  	    									//cin>>lll;

  	    								}while(idDaughter1 == 23);

  	    								motherLine2 = motherLine1;
  	    								do{
  	    									if(_verbose)cout<<" lepton cycle "<<endl;
  	    									lineDaughter1 = pythia->event.at(motherLine1).daughter1();
  	    									lineDaughter2 = pythia->event.at(motherLine2).daughter2();
  	    									if(pythia->event.at(motherLine1).id() != 23){

  	    										lineDaughter12 = pythia->event.at(motherLine1).daughter2();
  	    										if(lineDaughter1 == lineDaughter12){
  	    											if(_verbose)cout<<" all right lineDaughter1 == lineDaughter12 "<<endl;
  	    										}else{

  	    											if (abs(pythia->event.at(lineDaughter1).id()) == 13){
  	    												if(_verbose)cout<<" pythia->event.at(lineDaughter1).id() = "<<pythia->event.at(lineDaughter1).id()<<endl;
  	    											}

  	    											if (abs(pythia->event.at(lineDaughter12).id()) == 13){
  	    												lineDaughter1 = lineDaughter12;
  	    												if(_verbose)cout<<" pythia->event.at(lineDaughter12).id() = "<<pythia->event.at(lineDaughter12).id()<<endl;
  	    											}


  	    										}


  	    										lineDaughter22 = pythia->event.at(motherLine2).daughter1();
  	    										if(lineDaughter2 == lineDaughter22){
  	    											if(_verbose)cout<<" all right lineDaughter2 == lineDaughter22 "<<endl;
  	    										}else{

  	    											if (abs(pythia->event.at(lineDaughter2).id()) == 13){
  	    												if(_verbose)cout<<" pythia->event.at(lineDaughter2).id() = "<<pythia->event.at(lineDaughter2).id()<<endl;
  	    											}

  	    											if (abs(pythia->event.at(lineDaughter22).id()) == 13){
  	    												lineDaughter2 = lineDaughter22;
  	    												if(_verbose)cout<<" pythia->event.at(lineDaughter22).id() = "<<pythia->event.at(lineDaughter22).id()<<endl;
  	    											}


  	    										}


  	    									}
  	    									idDaughter1 = pythia->event.at(lineDaughter1).id();
  	    									idDaughter2 = pythia->event.at(lineDaughter2).id();
  	    									if(_verbose){
  	    										cout<<" lineDaughter1 =  "<<lineDaughter1<<endl;
  	    										cout<<" idDaughter1 =  "<<idDaughter1<<endl;
  	    										cout<<" lineDaughter2 =  "<<lineDaughter2<<endl;
  	    										cout<<" idDaughter2 =  "<<idDaughter2<<endl;
  	    									}
  	    									//if()
  	    									motherLine1 = lineDaughter1;
  	    									motherLine2 = lineDaughter2;
  	    							//		cin>>lll;
  	    								}while(pythia->event.at(motherLine1).daughter1() != 0 && pythia->event.at(motherLine2).daughter1() != 0);

  	    								if(_verbose)cout<<" pythia->event.at(motherLine1) id "<<pythia->event.at(motherLine1).id()<<endl;
  	    								if(_verbose)cout<<" pythia->event.at(motherLine2) id "<<pythia->event.at(motherLine2).id()<<endl;
  	    								if(pythia->event.at(motherLine1).idAbs()==13
  	    										&& pythia->event.at(motherLine2).idAbs()==13){



  	    									if(_verbose)cout<<" partA "<<pythia->event.at(motherLine1).id()<<" "<<pythia->event.at(motherLine1).p();//<<" partA x1 = "<<pythia->info.x1()<<" partA x2 = "<<pythia->info.x2()<<endl;

  	    									if(_verbose)cout<<" partB "<<pythia->event.at(motherLine2).id()<<" "<<pythia->event.at(motherLine2).p()<<" partB x1 = "<<pythia->info.x1()<<" partB x2 = "<<pythia->info.x2()<<endl;






  	    									particleC->p(pythia->event.at(motherLine1).p() + pythia->event.at(motherLine2).p());

  	    									//fileDrellYan<<iEvent<<" "<<pythia->event.at(motherLine1).p();
  	    									//fileDrellYan<<iEvent<<" "<<pythia->event.at(motherLine2).p();
  	    									fileDrellYan<<iEvent<<" "<<particleC->p();
  	    									//cout<<iEvent<<" "<<particleC->p();
  	    									if(_verbose)cout<<" partC "<<" "<<particleC->p();
  	    									if(_verbose)cout<<" ecm = "<<pythia->info.eCM()<<endl;
  	    									if(_verbose)cout<<"XF = "<<particleC->pz()/particleA->e()<<"               "<<pythia->info.x1() - pythia->info.x2()<<endl;
  	    								//	cout<<"XF = "<<2*particleC->pz()/pythia->info.eCM()<<"               "<<pythia->info.x1() - pythia->info.x2()<<endl;
  	    					//				cin>>ch;
  	    									flagQuit = true;
  	    								}else{
  	    									cout<<" pythia->event.at(motherLine1) "<<pythia->event.at(motherLine1).id()<<" dauther = "<<pythia->event.at(motherLine1).daughter1()<<endl;
  	    									cout<<" pythia->event.at(motherLine2) "<<pythia->event.at(motherLine2).id()<<" dauther = "<<pythia->event.at(motherLine2).daughter1()<<endl;
  	    									//cin>>lll;
  	    								}





  	    							}// end of  Drell- Yan part





  	    				//break;





  	    		//fileDrellYan<<iEvent<<" "<<i<<" "<< pythia->event[i].id()<<" "<<pythia->event[i].p();
  	    		//cout<<iEvent<<" "<<i<<" "<<pythia->event[i].id()<<" "<<pythia->event[i].p();

  	    //	}
  	    }


  	  // End of event loop. Statistics. Histogram. Done.
  	  }
  	  pythia->stat();
  	fileDrellYan.close();


  return 0;
}
