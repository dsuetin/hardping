///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <stdlib.h>
//#include "timer.h"
//#include "PythiaConstants.h"


//#include "Pythia.h"

#include "../Pythia8/Pythia.h"
#include "../Pythia8/nucleus.h"
using namespace std;
using namespace Pythia8;
extern Pythia* pythia8;
extern ofstream coordinateSoftOutput;
extern ofstream coordinateHardOutput;
extern  double getRandomFromFile();

extern ofstream probabilityOutput;
extern ofstream impactParameterFileAfter;
extern ofstream impactParameterFileBefore;
//using namespace starlightConstants;
Timer time1;
nucleus::nucleus(const int Z,const int A)
:	_Z(Z),
 	_A(A),
	_precisionOfNuclearThicknessCalculation(0.000000001)
{
	switch (_A) {
		case 1:
			{
				_r0 = 1.4*pow(_A,1./3.);
			}
			break;
		case 2:
			{
				_r0   = 1.95;
			}
			break;
		default:
			cout << "density not defined for projectile with Z = " << _Z << ". using defaults." << endl;
		//The density of the nucleus is
			_r0 = 1.19*pow(_A,1./3.) - 1.61*pow(_A, -1./3.);
			cout << " r = "<< _r0<<endl;
		break;
		}

		//_r0 = 1.16 * (1. - 1.16 * pow(_A, -2. / 3.));  // for FRITIOF and FormFactor.
		cout<<"here we are "<<" A = "<<_A<<endl;
}

//______________________________________________________________________________
Hardping::Hardping(nucleus projectileNucleus,
				   nucleus targetNucleus     )
	: //_Z(targetNucleus.Z()),
	 // _A(targetNucleus.A()),
	  //_precisionOfNuclearThicknessCalculation(0.001),
	 // _maxZCoordinate(2*targetNucleus.getNuclearRadius()),
	   _BQ( 2./0.7),
	   _BP (2./1.3),
	   _BN (2./0.1),
	  _randomSeed(19780503),
	  _energyCut(5.),
	  _firstCall(true),
	  _isScattering(false),
	  _kEnergyLoss(3.0),
	  _projectileNucleus(projectileNucleus),
	  _targetNucleus(targetNucleus),
	  _indexParticle(0),
	  _scatteringOnProton(0),
	  _scatteringOnNeutron(0),
	  _softInteractionFlag(false),
	  _hardInteractionFlag(false),
	  _hardInteraction(false),
	  _softToHard(false),
	  _energyLab(800.0),
	  _hardInteractionCount(0),
	  _softInteractionCount(0),
	  _fortranHardping(false),
	  _verbose(0),
	  _cutMass(true)
	 // _randomFile("/home/dsuetin/workspace/Pythia8180/Debug/randomNumbersFile.txt"),
	//  _file(0),
//	  _time(),
	  //_random->init(_randomSeed),


{

	if(_targetNucleus.A() == 2 ) _BN = 2./0.01;
	if(_verbose)cout<<"before py"<<endl;
	//time1.start();
	//pythia = new Pythia("/home/dsuetin/workspace/Pythia8/xmldoc");
	pythia = pythia8;
	//double ts = time1.stop();

	  _BQ = 2./0.7;
	  _BP = 2./1.3;
	  _BN = 2./0.1;


	//cout<<"time stop = "<<ts<<endl;
	//time1.printTime(ts);
	char ch;


	//pythiaHard = new Pythia("/home/dsuetin/workspace/Pythia8/xmldoc");
	if(_verbose)cout<<"after py "<<endl;
	//_random->init(_randomSeed);
	/*randomGeneratorInitialize(_randomSeed);*/
//	pythia->readString("Random:setSeed = on");


	switch (targetNucleus.A()){
	case 9:
		{
			_maxZCoordinate = 8.1812924800334574;
		}
		break;
	case 184:
		{
			_maxZCoordinate = 12.93;
		}
		break;
	case 1:
		{
			_maxZCoordinate = 6.06;
		}
		break;
	default:
		cout << "using defaults _maxZCoordinate." << endl;
	//The density of the nucleus is
		_maxZCoordinate = 10.;

	break;
	}
///  impact parameter///////////////////
	_impactParameterMax = min(3*_targetNucleus.getNuclearRadius(),_maxZCoordinate + 6.06); // from Fortran version. 6.06 fm - is proton characteristic
	_impactParameterMin = 0;
////////////////////////////////////////

    _generation = NULL; // Initialization with NULL pointer.
	_generations = new vector <generation>;
	_finalState = NULL;
	_finalState =  new vector <hardpingParticle>;
	_notInit = NULL;
	_notInit = new vector <hardpingParticle>;
	_cutMassive = NULL;
	_cutMassive = new vector <hardpingParticle>;
	_index = NULL;
	_index =  new vector <index>;
	_indexBadInitializations = NULL;
	_indexBadInitializations =  new vector <unsigned int>;
	_outOfNucleus = NULL;
	_outOfNucleus = new vector <hardpingParticle>;
	_tempParticle = NULL;
	_vertexOfInteraction = NULL;
	_vertexOfInteraction =  new vector <unsigned int>;
	_indexSoftToHard = NULL;
	_indexSoftToHard = new vector <unsigned int>;
	_indexSoftToHardBadInit = NULL;
	_indexSoftToHardBadInit = new vector <unsigned int>;
	_vertexOfBadHardInteraction = NULL;
	_vertexOfBadHardInteraction = new vector <unsigned int>;
//cin>>ch;

//	int enter;
//	cout <<"proj = "<< _projectileNucleus.A()<<" "<<_projectileNucleus.Z()<< " target "<<_targetNucleus.A()<<" "<<_targetNucleus.Z()<<endl;
	//cin >>enter;
}


//______________________________________________________________________________
Hardping::~Hardping()
{
	//cout<<" iam in destructor"<<endl;
//	pythia->rndm.dumpState("randomState_Be");
//	delete pythia;
//	delete pythiaHard;
	delete _generations;
	delete _finalState;
	delete _notInit;
	delete _cutMassive;
	delete _index;
	delete _indexBadInitializations;
	delete _outOfNucleus;
	delete _vertexOfInteraction;
	delete _indexSoftToHard;
	delete _indexSoftToHardBadInit;
	delete _vertexOfBadHardInteraction;
	_randomFile.close();
	//pythia->rndm.dumpState("randomState");
}
void
Hardping::setVaribles(){
	_softInteractionCount = 0;
	_vertexOfInteraction->resize(0);
	_vertexOfBadHardInteraction->resize(0);
	_indexSoftToHardBadInit->resize(0);
	_outOfNucleus->resize(0);
	_indexBadInitializations->resize(0);
	_cutMassive->resize(0);
	_generations->resize(0);
	 _firstCall = true;
	 _isScattering = false;
	 _indexParticle = 0;
	_scatteringOnProton = 0;
	_scatteringOnNeutron = 0;
	_softInteractionFlag = false;
	_hardInteractionFlag = false;
	_hardInteraction = false;
	_hardInteractionCount = 0;
//	_softInteractionCount = 0;
	//_softCollisionIterator. = NULL;
	_softToHard = false;
	_indexSoftToHard->resize(0);
	return;




/*


	  _maxZCoordinate = 15.;
	  _firstCall =  true;
	  _isScattering = false;
	  _randomSeed = 19780503;
	  _indexParticle = 0;
	  _scatteringOnProton = 0;
	  _scatteringOnNeutron = 0;
	  _softInteractionFlag = false;
	  _hardInteractionFlag = false;
	  _hardInteraction = false;
	  _softToHard = false;
	  _hardInteractionCount = 0;
	  _softInteractionCount = 0;



	pythia = new Pythia("/home/dsuetin/workspace/Pythia8/xmldoc");
	pythiaHard = new Pythia("/home/dsuetin/workspace/Pythia8/xmldoc");
	cout<<"after py "<<endl;
	//_random->init(_randomSeed);

//	pythia->readString("Random:setSeed = on");
	pythia->rndm.init(_randomSeed);
	 _generation = NULL; // Initialization with NULL pointer.
	_generations = new vector <generation>;
	_finalState = NULL;
	_finalState =  new vector <hardpingParticle>;
	_notInit = NULL;
	_notInit = new vector <hardpingParticle>;
	_cutMassive = NULL;
	_cutMassive = new vector <hardpingParticle>;
	_index = NULL;
	_index =  new vector <index>;
	_indexBadInitializations = NULL;
	_indexBadInitializations =  new vector <unsigned int>;
	_outOfNucleus = NULL;
	_outOfNucleus = new vector <hardpingParticle>;
	_tempParticle = NULL;
	_vertexOfInteraction = NULL;
	_vertexOfInteraction =  new vector <unsigned int>;
	_indexSoftToHard = NULL;
	_indexSoftToHard = new vector <unsigned int>;
	_indexSoftToHardBadInit = NULL;
	_indexSoftToHardBadInit = new vector <unsigned int>;
	_vertexOfBadHardInteraction = NULL;
	_vertexOfBadHardInteraction = new vector <unsigned int>;




*/






}

//______________________________________________________________________________
double
nucleus::nuclearThicknessGauss5(const double impactParameter, double zMin, double zMax) const
{
	//    JS      This code calculates the nuclear thickness function as per Eq. 4 in
	//    Klein and Nystrand, PRC 60.
	//    former DOUBLE PRECISION FUNCTION T(b)
                                                  
	// data for Gauss integration
	double impactParameter2 = impactParameter*impactParameter;
	double nuclearThicness = 0.;
	const unsigned int nmbPoints         = 5;
	const double       xg[nmbPoints + 1] = {0., 0.1488743390, 0.4333953941, 0.6794095683,
	                                        0.8650633667, 0.9739065285};
	const double       ag[nmbPoints + 1] = {0., 0.2955242247, 0.2692667193, 0.2190863625,
	                                        0.1494513492, 0.0666713443};
  
	const double zRange = 0.5 * (zMax - zMin); 
	const double zMean  = 0.5 * (zMax + zMin); 
	double       sum    = 0;
	for(unsigned int i = 1; i <= nmbPoints; ++i) {
		double zsp    = zRange * xg[i] + zMean;
		double radius = sqrt(impactParameter2 + zsp * zsp);
		sum          += ag[i] * getNuclearDensity(radius);
		zsp           = zRange * (-xg[i]) + zMean;
		radius        = sqrt(impactParameter2 + zsp * zsp);
		sum          += ag[i] * getNuclearDensity(radius);
	}
	nuclearThicness = zRange * sum;
	//_lastResultOfThicknessCalculation = nuclearThicness;
	return nuclearThicness;
}
double
nucleus::nuclearThicknessGauss12(double impactParameter, double zMin, double zMax) const
{
   //  Return Integral of function between a and b.
   const double impactParameter2 = impactParameter*impactParameter;
  // cout<<"impactParameter2 = "<<impactParameter2;
   const double kHF = 0.5;
   const double kCST = 5./1000;

   double x[12] = { 0.96028985649753623,  0.79666647741362674,
                      0.52553240991632899,  0.18343464249564980,
                      0.98940093499164993,  0.94457502307323258,
                      0.86563120238783174,  0.75540440835500303,
                      0.61787624440264375,  0.45801677765722739,
                      0.28160355077925891,  0.09501250983763744};

   double w[12] = { 0.10122853629037626,  0.22238103445337447,
                      0.31370664587788729,  0.36268378337836198,
                      0.02715245941175409,  0.06225352393864789,
                      0.09515851168249278,  0.12462897125553387,
                      0.14959598881657673,  0.16915651939500254,
                      0.18260341504492359,  0.18945061045506850};

   double nuclearThicness, aconst, bb, aa, c1, c2, u, s8, s16, f1, f2;
   double xx[1];
   int i;
   int countInfityLoop1 = 0, countInfityLoop2 = 0;


   nuclearThicness = 0;
   if (zMax == zMin) return nuclearThicness;
   aconst = kCST/std::abs(zMax-zMin);
   bb = zMin;
CASE1:
   aa = bb;
   bb = zMax;
CASE2:
   c1 = kHF*(bb+aa);
   c2 = kHF*(bb-aa);
   s8 = 0;
   for (i=0;i<4;i++) {
      u     = c2*x[i];
      xx[0] = c1+u;
      xx[0] = sqrt(xx[0]*xx[0] + impactParameter2);
      f1    = getNuclearDensity(xx[0]);
  //    cout<<setprecision(10)<<" x = "<<xx[0]<<" rws = "<< f1<<endl;
   //   if (fgAbsValue) f1 = std::abs(f1);
      xx[0] = c1-u;
      xx[0] = sqrt(xx[0]*xx[0] + impactParameter2);
      f2    = getNuclearDensity(xx[0]);
 //     cout<<" x = "<<xx[0]<<" rws = "<< f2<<endl;
    //  if (fgAbsValue) f2 = std::abs(f2);
      s8   += w[i]*(f1 + f2);
   }
   s16 = 0;
   for (i=4;i<12;i++) {
      u     = c2*x[i];
      xx[0] = c1+u;
      xx[0] = sqrt(xx[0]*xx[0] + impactParameter2);
      f1    = getNuclearDensity(xx[0]);
    //  cout<<" x = "<<xx[0]<<" rws = "<< f1<<endl;
      //if (fgAbsValue) f1 = std::abs(f1);
      xx[0] = c1-u;
      xx[0] = sqrt(xx[0]*xx[0] + impactParameter2);
      f2    = getNuclearDensity(xx[0]);
  //    cout<<" x = "<<xx[0]<<" rws = "<< f2<<endl;
     // if (fgAbsValue) f2 = std::abs(f2);
      s16  += w[i]*(f1 + f2);
   }
   s16 = c2*s16;
   if (std::abs(s16-c2*s8) <= _precisionOfNuclearThicknessCalculation*(1. + std::abs(s16))) {
      nuclearThicness += s16;
      if(bb != zMax) {
    	 // cout<<"infinity loop 1"<<endl;
    	  countInfityLoop1++;
    	  if(countInfityLoop1 > 10000)return 0;
    	  goto CASE1;
      }
   } else {
      bb = c1;
      if(1. + aconst*std::abs(c2) != 1){
    	//  cout<<"infinity loop 2"<<endl;
    	  countInfityLoop2++;
    	  if(countInfityLoop2 > 10000)return 0;
    	  goto CASE2;
      }
      nuclearThicness = s8;  //this is a crude approximation (cernlib function returned 0 !)
   }

   //_lastResultOfThicknessCalculation = h;
   //_lastErrorOfThicknessCalculation = std::abs(s16-c2*s8);

   return nuclearThicness;
}
double
nucleus::renormalizedNuclearThicknessGauss5(double impactParameter, double zMin, double zMax)const{
	return nuclearThicknessGauss5(impactParameter, zMin, zMax)*(_A-1)/_A;
}
double
nucleus::renormalizedNuclearThicknessGauss12(double impactParameter, double zMin, double zMax)const{
	return nuclearThicknessGauss12(impactParameter, zMin, zMax)*(_A-1)/_A;
}
double
Hardping::probabilityOfNSoftScattering(double impactParameter, double zMin, double zMax, int numberOfCollisions)const{
 //   DMYSIN=((0.1D0*sigma*TA_(b,z,HIPR1(35)))**(n))
 //  &    *exp(-0.1D0*sigma*TA_(b,z,HIPR1(35)))
	double probabilityOfNSoftScattering = 0;
	if (numberOfCollisions > 12) return probabilityOfNSoftScattering;
	double quarkNucleonCrossSection = 10.; // in mb
	quarkNucleonCrossSection = 0.1*quarkNucleonCrossSection; // convert mb in fb
	double temp = quarkNucleonCrossSection*_targetNucleus.renormalizedNuclearThicknessGauss12(impactParameter,zMin,_maxZCoordinate);
	probabilityOfNSoftScattering = pow(temp,numberOfCollisions)*exp(-temp)/factorial(numberOfCollisions); // mul nuclear density?
	if(probabilityOfNSoftScattering>1)cerr<<"probability above 1"<<endl;
	//cout<<"probability = "<< probabilityOfNSoftScattering<<endl;
	return probabilityOfNSoftScattering;

}
double
Hardping::integralProbabilityOfNSoftScattering(int numberOfCollisions) const{
	double deltaZ, deltaImpactParameter, zMin = -15.,impactParameter = 0., z;
	double totalProbability = 0.;
	int nBins =100;
	deltaImpactParameter = (10. - 0.)/nBins;
	deltaZ = (_maxZCoordinate - zMin)/nBins;
	for(int ib = 1; ib <= nBins; ib++){
		impactParameter = ib*deltaImpactParameter;
		for(int iz = 0; iz < nBins; iz++){
			z = zMin + iz*deltaZ;
			totalProbability += 2*M_PIl*impactParameter*deltaZ*deltaImpactParameter*probabilityOfNSoftScattering(impactParameter, z, _maxZCoordinate, numberOfCollisions)*_targetNucleus.getNuclearDensity(sqrt(impactParameter*impactParameter+z*z));
		}
	}
	return totalProbability;
}


bool
Hardping::pathInNucleus2( hardpingParticle * particleA , double &zCoordinateOfCollision){
	_isScattering = false;

	int numOfCollisions = 0;
	int numCollisions = 0;
//	cout<<"parA "<<particleA->vProd();
	double crossSection1 = 10.0, crossSection2 = 10.0;
	double impactParameter;
	double path = 0.;
	double totalPath =0;
    int count = 0;
    int countInfinityLoop = 0;
    double xCoord = 0, yCoord = 0, zCoord = 0;
    double HIPR135 = 0;
    char ch;
    HIPR135 = _maxZCoordinate;
    //cout<<" HIPR135 = "<<HIPR135<<endl;
    bool isScattering = false;


	xCoord = particleA->vProd().px(); //DXC
	yCoord = particleA->vProd().py(); //DYC
	zCoord = particleA->vProd().pz();
	impactParameter = particleA->vProd().pT();
//	cout<<"coord  ggg= "<<particleA->vProd()<<"zCoordinateOfCollision = "<<zCoordinateOfCollision<<endl;
	//zCoord = particleA->vProd().pz(); // -HIPR1(35)
	//
	//if(_firstCall)zCoord = - HIPR135;


	//zCoord = particleA->vProd().pz();

/*
	particleA->setAngles();
	impactParameter = particleA->vProd().pT();
	particleA->rotateHardping();
	impactParameterFileBefore<<particleA->vProd();
	particleA->rotateBackHardping();

	*/
	int INOFCOLL = 0;
	int	INCOLL = 0;
	int	ICOUNT = 0;
	double KQ = 0.7;
	double KP = 1.3;
	double KN = 0.1;
	if(_targetNucleus.A() == 2 ) KN = 0.01;
    double BQ=2./KQ;
    double BP=2./KP;
    double BN=2./KN;
    double xMaxP0 = HIPR135;
    double yMaxP0 = 1.;
    double P0 = 0.;
    double DIFSC = 1.;
    double xMaxP, xMinP, yMaxP, z;
    double X,Y, YF;
   // double r = sqrt(_zPartonCoord*_zPartonCoord+impactParameter*impactParameter);
    double r = sqrt(zCoord*zCoord+impactParameter*impactParameter);
  //  cout<<"impactParameter = "<<impactParameter<<" r = "<<r<<endl;
   // cin>>ch;
    do{
    	countInfinityLoop++;

        do{

        	if(_verbose)cout<<" r "<<r<<" xMaxP0 = "<<xMaxP0<<" impactParameter = "<<impactParameter<<endl;
            if ( true || zCoord < 0){
            	//cout<<"r < xMaxP0 || zCoord < 0"<<endl;
            	if(numOfCollisions == 0){
            	//	cout<<"impactParameter "<<impactParameter<< " zCoord = "<<zCoord<<" HIPR135 = "<<HIPR135<<endl;
            	//	cout<<"huinya = "<<_targetNucleus.renormalizedNuclearThicknessGauss12(impactParameter,zCoord,HIPR135)<<endl;
            		 P0 = exp(-0.1*crossSection1*_targetNucleus.renormalizedNuclearThicknessGauss12(impactParameter,zCoord,HIPR135));//probability of zero scattering
            		 //probabilityOutput<<P0<<endl;
            		 cout.precision(16);
            		// cout<< "P01 = "<<P0<<endl;
            	//	 cin>>ch;
            	}else{
           // 		cout<<"2impactParameter "<<impactParameter<< " zCoord = "<<zCoord<<" HIPR135 = "<<HIPR135<<endl;
            //		cout<<"2huinya = "<<_targetNucleus.renormalizedNuclearThicknessGauss12(impactParameter,zCoord,HIPR135)<<endl;
            		 P0 = exp(-0.1*crossSection2*_targetNucleus.renormalizedNuclearThicknessGauss12(impactParameter,zCoord,HIPR135));
         //   		 cout<< "P02 = "<<P0<<endl;
            		// cin>>ch;
            	}
            }
            if(_verbose)cout<<" p0 = "<<P0<<endl;
      //      cout<<"nsoft = "<<_softInteractionCount<<endl;
            if(P0 == 0){

            	cout<<"p0= 0 ;impactParameter = "<<impactParameter<<" zCoord = "<<zCoord<<endl;
            	cin>>ch;
            	return false;
            }



     ///////////////////////////////////////////////////////////////////
    //// this block correspond to exist scattering or not///////////////
            do{
         //   	cout<< "in cycle "<<endl;
                X= xMaxP0*getRandom();//getRandomFromFile();//getRandom();		// x is NOT coordinate
                Y= yMaxP0*getRandom();//getRandomFromFile();//getRandom();
                if(_verbose)cout<<"xMaxP0 = "<<xMaxP0<<" X = "<<X<<" yMaxP0 = "<<yMaxP0<<" Y = "<<Y<<endl;
              //  cin>>ch;
                if ( X > xMaxP0/2.){
                	YF = 1. - P0;
                	_isScattering = true;
                	DIFSC = 1.;
                }else{
                	YF = P0;
                	_isScattering = false;
                	DIFSC = 0.;
                }
                ICOUNT = ICOUNT + 1;
                //count++;
                if(/*count*/ ICOUNT > 10){
                	yMaxP0 = yMaxP0/2.;
                	ICOUNT = 0;
                	//count = 0;
                }
            }while(Y > YF);
            isScattering = _isScattering;



            //change suetin 11.07.14 ///
            if(0&&_firstCall){
            	isScattering = true;
            	DIFSC = 1;
            	_isScattering = true;
            }

    ///////////////////end of block/////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
            //count = 0;
         //   if(isScattering)cin>>ch;
            ICOUNT = 0;
            //count = 0;

            xMaxP = HIPR135;

           // xMaxP = -HIPR135/2.;


            xMinP = zCoord;
       //     cout<<"xMaxP = "<<xMaxP<<" xMinP = "<<xMinP<<endl;
            //cin>>ch;
        //    cout<<" xMaxP = "<<xMaxP<<endl;
         //   cout<<" xMinP = "<<xMinP<<endl;
            yMaxP = 1. - P0;
            //if(yMaxP < pow(10.,-8))return 0;
    //       cout<<" yMaxP = "<<yMaxP<<endl;
    //    /    cin>>ch;

            if(_verbose)cout<<"xMinP =  "<<xMinP<<" yMaxP = "<<yMaxP<<endl;
            // if first collision do not happened we set new impact parameter for incident particle and try again
        	if(/*!_isScattering && _firstCall*/0){
        		getImpactParameter(xCoord,yCoord);
        		particleA->vProd().px(xCoord);
        		particleA->vProd().py(yCoord);
        		impactParameter = particleA->vProd().pT();
        		r = sqrt(zCoord*zCoord+impactParameter*impactParameter);
        	}
        	/////////////////////////////////////////////////////////////////////////////////////////////////////
        }while(/*!_isScattering && _firstCall*/0); //this condition correspond to case that collision of first incident particle happened always

        if(/*_isScattering DIFSC && INOFCOLL < _targetNucleus.A() */  _isScattering && _softInteractionCount < _targetNucleus.A()){
        	//numOfCollisions++;
        	INOFCOLL++;
        	do{

            	X = xMinP + (xMaxP -xMinP)*getRandom();//*getRandom();
      //      	cout<<" X = "<<X<<endl;
            	Y = yMaxP*getRandom();//getRandomFromFile();//*getRandom();
            	if(_verbose)cout<<"xMinP = "<<xMinP<<" xMaxP "<<xMaxP<<" yMaxP = "<<yMaxP<<" p0 "<<P0<<endl;
            	if(_verbose)cout<<" Xpath = "<<X<<"  Ypath = "<<Y<<" yMaxP = "<<yMaxP<<endl;
            	if(yMaxP < pow(10.,-17)){

            		cout<<"yMaxP < 10^ 17 "<<endl;
            	//	cin>>ch;
            		return 0;
            	}
            	z = X;
            	ICOUNT++;
            	//count++;
            	if(/*numCollisions*/INCOLL == 0){
            //		cout<<"im in numCollisions = 0"<<endl;

            		path = exp(-0.1*crossSection1*_targetNucleus.renormalizedNuclearThicknessGauss12(impactParameter,zCoord,z)) - P0;
            		 probabilityOutput<<path<<endl;
         //   		 cout<<"path1 = "<<path<<endl;
            		 if(path< -0.1)cin>>ch;
            	//	cout<<"path = "<<path<<endl;
            	}else{
          //  		cout<<"im in numCollisions != 0"<<endl;
            		path = exp(-0.1*crossSection2*_targetNucleus.renormalizedNuclearThicknessGauss12(impactParameter,zCoord,z)) - P0;
            		 probabilityOutput<<path<<endl;
           // 		 cout<<"path2 = "<<path<<endl;
            		 if(path< -0.1)cin>>ch;
            	}
            	//if(_verbose)cout<<"path = "<<path<<endl;
            	if(ICOUNT > 10){
            		yMaxP = yMaxP/2.0;
            		//count = 0;
            		ICOUNT = 0;
            	}
           // 	cout<<"infinity loop "<<endl;
            	if(_verbose)cout<<"Y = "<<Y<<" path = "<<path<<endl;
        	}while(Y >= path);
        	//numCollisions++;
        	INCOLL = INCOLL+1;
        	totalPath = totalPath + (z - zCoord);
        	zCoord = z;
       // 	cout<<" z1 = "<<z<<endl;
       // 	cin>>ch;
        	/*
        	particleA->rotateHardping();
        	impactParameterFileAfter<<particleA->vProd();
        	particleA->rotateBackHardping();
        	 */
        	zCoordinateOfCollision = z;
        	if(_verbose)cout<<" z =  "<<z<<endl;
        //	cin>>ch;
      //  //	zCoordinateOfCollision = totalPath;
        	//return isScattering;
        	return DIFSC;
        	//return totalPath;

 //       	return (z - zCoord);
        	//_zPartonCoord = z;
        	//cout<<" return z = "<< z<<endl;
        	//return z;
        }else{
        	if(_verbose)cout<< "return zero "<<endl;

        	if(!_isScattering){
        		if(_verbose)cout<<" z =  "<<z<<endl;
        	//	cin>>ch;
        	}

        	return false;//_isScattering;
        }
        //totalPath = totalPath +(z - _zPartonCoord);


    }while(countInfinityLoop < 10000000000);
  // cerr<<"infinity loop happened in func, 0.000 is returned"<<endl;

   return 0.;
}
double hardpingParticle::getMaxTransverseMomentum(int type){
	return 1;
}


/*double hardpingParticle::getNewPtInitialState(int type){
	double absolutePartonMomentum = this->pAbs();
	double maxTransverseMomentum  = this->getMaxTransverseMomentum(type);
	double newTransverseMomentum  = 0.;
	this->readRandomStage("random");
	_random->flat();
	cout<<"random number = "<<_random->flat()<<endl;
	double B = 0.;
	double _BP, _BN ;
	switch (type) {
			case 1:
				{
					B = _BP;
				}
				break;
			case 2:
				{
					B   = _BN;
				}
				break;
			default:
				cout << "in getNewPtInitialState() . Type = " << type << " not found." << endl;
				break;
		}

	do{
	//	randomNumber =getRandom();
		//newTransverseMomentum = -1./B*log(randomNumber*randomNumber);
	}while(newTransverseMomentum > absolutePartonMomentum
		|| newTransverseMomentum > maxTransverseMomentum);
//	anglePhi = 2.*M_PIl*getRandom();
	this->saveRandomStage("randomState");
	return 1;
}
*/

double
Hardping::getNewPtInitialState(hardpingParticle * particleA ,int type){

	char ch;
	cout.precision(12);
	double absoluteNucleonMomentum = particleA->pAbs();
	double newTransverseMomentum =0.;
	double B = 0.;
	switch (type) {
			case 1:
				{
					B = _BQ;
				}
				break;
			case 2:
				{
					B   = _BN;
				}
				break;
			default:
				cout << "in getNewPtInitialState() . Type = " << type << " not found." << endl;
				break;
		}

	double randomNumber = 0., anglePhi = 0., tempRandom1,tempRandom2;


	do{
		//randomNumber =getRandomFromFile();//getRandom();
		tempRandom1 = getRandom();
		tempRandom2 = getRandom();
		//cout<<"tempRandom1 = "<<tempRandom1<<" tempRandom2 = "<<tempRandom2<<" B "<<B<<endl;
		newTransverseMomentum = -1./B*log(tempRandom1*tempRandom2);//(randomNumber*randomNumber);
		//cout<<"newTransverseMomentum = "<<newTransverseMomentum<<endl;
	}while(newTransverseMomentum > absoluteNucleonMomentum);

	double theta, phi;//, newTheta, newPhi;

	theta = particleA->theta();
	phi = particleA->phi();

	//if(_verbose)cout<<"theta = "<<theta<<endl;
	//if(_verbose)cout<<"phi   = "<<phi<<endl;
	if(_verbose)cout<<" momentum  before  rotate   = "<<particleA->p()<<endl;
	//if(_verbose)cout<<" 3coordinate before rotation = "<<particleA->vProd()<<endl;
//	particleA->setAngles();




//	particleA->rotateBackHardping();



	//particleA->rot(0,-phi);
	//particleA->rot(-theta,0);
	//if(_verbose)cout<<" 3momentum  after rotate_    = "<<particleA->p()<<endl;
	//if(_verbose)cout<<" 3coordinate after rotation_ = "<<particleA->vProd()<<endl;
	//particleA->rotateHardping();
	//if(_verbose)cout<<" 000000000000000000000000000    = "<<particleA->p()<<endl;
	//particleA->rotateBackHardping();
	anglePhi = 2.*M_PIl*getRandom();//getRandom();
	if(_verbose)cout<<" phi = "<<anglePhi<<" new momentum = "<<newTransverseMomentum<<endl;
	//cout<<" anglePhi "<<anglePhi<<" newTransverseMomentum ="<<newTransverseMomentum<<endl;
	double px = 0, py = 0, pz = 0;
	double xc = 0, yc = 0, zc = 0;
	double newTheta = 0, newPhi = 0;
	px = newTransverseMomentum * cos(anglePhi);
	py = newTransverseMomentum * sin(anglePhi);

	pz = sqrt(absoluteNucleonMomentum * absoluteNucleonMomentum - newTransverseMomentum * newTransverseMomentum);
	newTheta = atan2(sqrt(px*px + py*py), pz);
	newPhi   = atan2(py,px);
	//cout<<"px = "<<px<<" py = "<<py<<" pz = "<<pz<<endl;
	//rotate correspond new pt

	double cosPhi = py/sqrt(px*px+py*py);
	double sinPhi = px/sqrt(px*px+py*py);
	double sinTheta = sqrt(px*px+py*py)/absoluteNucleonMomentum;
	double cosTheta = pz/absoluteNucleonMomentum;
	Vec4 tempVector(px,py,pz,0);
	Vec4 tempVectorCoord(0);

//	cout<<"vec "<<tempVector<<endl;
//	cout.precision(10);
	if(_verbose)cout<<"cosPhi = "<<cosPhi<<" sinPhi = "<<sinPhi<<endl;
	if(_verbose)cout<<" sinTheta = "<<sinTheta<<" cosTheta = "<<cosTheta<<endl;

//	px = particleA->p().px();
//	py = particleA->p().py();
//	pz = particleA->p().pz();
//	cout<<"px = "<<px<<" py "<<py<<" pz = "<<pz<<endl;
	//px = px;
	// momentum rotation;
	px = particleA->p().px();
    py = particleA->p().py()*cosTheta + particleA->p().pz()*sinTheta;
    pz =-particleA->p().py()*sinTheta + particleA->p().pz()*cosTheta;
    // coordinate rotation;
    xc = particleA->vProd().px();
    yc = particleA->vProd().py()*cosTheta + particleA->vProd().pz()*sinTheta;
    zc = -particleA->vProd().py()*sinTheta + particleA->vProd().pz()*cosTheta;
    //tempVector.reset();
    //cout<<"after px = "<<px<<" py "<<py<<" pz = "<<pz<<endl;
    tempVector.p(px,py,pz,0);
    tempVectorCoord.p(xc,yc,zc,0);
    particleA->p(tempVector);
    particleA->vProd(tempVectorCoord);
    if(_verbose)cout<<"  particleA 0 "<<particleA->p();
    cout.precision(10);
//	px = particleA->p().px();
//	py = particleA->p().py();
//	pz = particleA->p().pz();
    px = particleA->p().px()*cosPhi + particleA->p().py()*sinPhi;
    py = -particleA->p().px()*sinPhi + particleA->p().py()*cosPhi;


    //cord rotation
    xc = particleA->vProd().px()*cosPhi + particleA->vProd().py()*sinPhi;
    yc = -particleA->vProd().px()*sinPhi + particleA->vProd().py()*cosPhi;

    //cout<< " px = "<<px<<" py = "<<py<<" pz = "<<pz<<endl;
    //pz =-py*cosTheta + pz*cosTheta;
    tempVector.p(px,py,pz,sqrt(px*px+py*py+pz*pz+particleA->m2()));
    particleA->p(tempVector);
    tempVectorCoord.p(xc,yc,zc,0);
    particleA->vProd(tempVectorCoord);
	if(_verbose)cout<<" particleA 1 "<<particleA->p();

    //DPX1=DPX
    //DPY1=DPY*DCOSTHETA+DPZ*DSINTHETA
    //DPZ1=-DPY*DSINTHETA+DPZ*DCOSTHETA
    //WRITE(*,*)'DPX1 = ',DPX1, ' DPY1 = ',DPY1,' DPZ1 = ',DPZ1
    //DPX2=DPX1*DCOSFI+DPY1*DSINFI
    //DPY2=-DPX1*DSINFI+DPY1*DCOSFI
    //DPZ2=DPZ1
//	particleA->rot(newTheta+M_PIl/2.,0);
//	particleA->rot(0,newTheta+M_PIl*3/2.);
//	if(_verbose)cout<<" 1 particleA "<<particleA->p();
//	particleA->rot(0,newPhi);
//	particleA->rot(-newPhi,0);
//	if(_verbose)cout<<" 2 particleA"<<particleA->p();
	//rotate correspond to initial momentum


//	particleA->rotateHardping();





//	if(_verbose)cout<<" 3 particleA"<<particleA->p();
//	particleA->setAngles();

/*
	if(_verbose)cout<<" parA pt after     = "<<particleA->pT()<<endl;

	if(_verbose)cout<<"particleA momentum = "<<particleA->p();

	if(_verbose)cout<<"4 coordinate after rotation_ = "<<particleA->vProd()<<endl;
	if(_verbose)cout<<"4 momentum  after rotate_    = "<<particleA->p()<<endl;
	if(_verbose)cout<<" 77777777777777777777777777777777777777777"<<endl;
*/
	return 0;
}

double
Hardping::getMaxTransverseMomentum(hardpingParticle* particleA, int type){
	double kT = particleA->pAbs();
	double maxTransverseMomentum = 0;
	cout<<"kT = "<<kT;
			//getAbsPartonMomentum();
	double B = 0.;
	switch (type) {
		case 1:
			{
				B = _BP;
				cout<<" 1 case B = "<<B<<endl;
			}
			break;
		case 2:
			{
				B   = _BN;
				cout<<" 2 case B = "<<B<<endl;
			}
			break;
		default:
			cout << "in getMaxTransverseMomentum() . Type = " << type << " not found. return 0." << endl;
			break;
	}

	maxTransverseMomentum = B*B*kT;
	cout<<"maxTransverseMomentum = "<<maxTransverseMomentum<<endl;
	maxTransverseMomentum *= exp(-B*kT);
	double chislo = -B*kT;
	cout<<"chislo = "<<chislo<<endl;
	cout.precision(180);
	cout<<"e chislo = "<<exp(chislo)<<endl;
	return B*B*kT*exp(-B*kT);
}

void Hardping::hardping(){
	//cin.ignore();

char ch;
//cin>>ch;
	pythia->process.clear();
	//cout<<"pythia8->event.size() = "<<pythia->event.size()<<endl;
	//cin>>ch;
	//pythia->event.clear();
	int numOfPythiaProductedParticles;
	//producedParticlesInformation eventsInformations;
	//bool _firstCall = true;
	bool HardQCD = false;
	bool SoftQCD = true;
	double energy;
	useFortranMethod(true);
	int nProducedParticles = 0;
	int nProducedParticlesOld = 0;
	int numberOfGeneration = 0;
	int iMax = 0;
	int jMax = 0;
	int coutnFinalPy = 0;
	int coutnFinal = 0;
	int coutnFinalPyOld = 0;
	int countCut = 0;
	pythia->event.clear();
	std::vector <int> *vectorNProducedParticles;


	hardpingParticle *particleA, *particleB;
	Vec4 vecA(0), vecCoordinate(0),vecMomentum(0);
	vecA.pz(_energyLab);
	vecA.px(00.);
	vecA.py(00.);
	energy = sqrt( protonMass*protonMass + vecA.px()*vecA.px() + vecA.py()*vecA.py() + vecA.pz()*vecA.pz() );
	vecA.e(energy);


	particleA = new hardpingParticle();//(vecA,2212);
	particleA->p(vecA);
	particleA->vProd(vecCoordinate);
	particleA->id(2212);
	_initialParticle = *particleA;
	if(_verbose)cout<<"_maxZCoordinate = "<<_maxZCoordinate<<endl;
//todo add for _initialParticle idB
	//cout<<"vec4 "<<particleA->getVector()->pz()<<" id = "<<particleA->getId()<<endl;
	particleB = new hardpingParticle();
	particleB->id(2212);
	int pythiaNextFlag = 0;

	int countfinal=0;
	int newGenerationCount = 0;
	int newGenerationCountOld = 0;
	double phi = 0, theta = 0;
	double px = 0, py =0, pz = 0;
	double impactParameter = 0;
	int numberOfParticlesAtPreviousGeneration = 0;
	int isumMax = 0;
	double probability;
	int massTemp [10000];
	int indexCount = 0;
	int countFinalParticles = 0;
	int countInit = 0;
	bool flagInit = true;

	int sing = 1;
	int lll;
	unsigned int randomPossitionOfHardScattering;
	double deltaPath = 0;
	double zCoordinateOfCollisions = 0;
	double zCoordinateOfCollisionsTemp = 0;
	double deltaP = 0, deltaE = 0;
	double impactParameterMax = 0;
	double impactParameterMin = 0;
	double yImpact = 0., xImpact = 0., impactPar = 0.,phiImpact = 0.;
	int isNotAdsorbed = 0;
	bool isScattering = false;
	bool softToHardFlag = false;
	double x1 = 0;

	if(true)do{
//		/cout<<" begin " <<endl;
		//cin>>ch;
		if(numberOfGeneration == 1){
		//	cout<<"numberOfGeneration = 1"<<endl;
		}
		numberOfParticlesAtPreviousGeneration = 0;
		//	not execute in the first call////////////////////
		if(numberOfGeneration)if(_verbose)cout<<"size of previous generation =  "<<_generations->at(numberOfGeneration-1).getMatrix()->size()<<" number = "<<numberOfGeneration-1<<endl;
		///////////////////////////////////////////////////


			//if(_verbose)if(!_firstCall)cout<<" number of incident particle in previous generation "<<isumMax<<endl;

			numberOfParticlesAtPreviousGeneration = calculateNumberOfParticlesAtPreviousGeneration(numberOfGeneration);

			//////////////end of  calculating number of particles at previous generation/////////////////////////

			if(_verbose)if(!_firstCall)cout<<" number of incident particle in this generation "<<numberOfParticlesAtPreviousGeneration<<endl;

			 particleIndexation(numberOfGeneration, numberOfParticlesAtPreviousGeneration);
			//////////////////////////////////////////////////
			(_firstCall)? initGeteration(numberOfGeneration,_projectileNucleus.A()) : initGeteration(numberOfGeneration,numberOfParticlesAtPreviousGeneration) ;

			addGeneration();

			(_firstCall)? iMax = _projectileNucleus.A() : iMax = numberOfParticlesAtPreviousGeneration;
			for(int i_init = 0; i_init < iMax ; i_init++){
				// i_init - index number of particle in current generation

				if(_firstCall){
					setInitinalImpactAndIndex(particleA);
				}else{
					particleA->setNumberOfCurrentGeneration(numberOfGeneration);
					particleA->setInexNumber(i_init);
					getParticleFromPreviousGeneration(particleA);
					particleA->setHard(false);
					particleA->setSoft(false);

				//	cout<<"coord A = "<<particleA->vProd();
				//	cin>>ch;
					//	cout<<"p A = "<<particleA->p();
					// get vector of coordinate of particle after rotation
					vecCoordinate = particleA->vProd();
				}

			//	cout<< "part A"<<particleA->p();
				//cin>>ch;

				if(	!softToHardFlag){
				//	cout<<"after pA 1"<<endl;
					// return calculate z coordinate of collisions, if they occurred
					zCoordinateOfCollisions = 0;
					zCoordinateOfCollisionsTemp = 0;

				//	cout << "particleA->theta() 1 "<<particleA->theta()<<" particleA->phi() "<<particleA->phi()<<endl;

					particleA->setAngles();
				//	cout << "particleA->theta() 2 "<<particleA->theta()<<" particleA->phi() "<<particleA->phi()<<endl;
				//	theta = particleA->theta();
				//	phi =   particleA->phi();

					// rotate to z' system
					//cout<<"particleA 1 = "<<particleA->p();
					cout.precision(12);
					if(0){

						cout<<"pxa = "<<particleA->p().px()<<endl;

						cout<<"pya = "<<particleA->p().py()<<endl;

						cout<<"pza = "<<particleA->p().pz()<<endl;

					}
					//particleA->rot(0,-phi);
					//particleA->rot(-theta,0);
					particleA->rotateBackHardping();
					if(0){

						cout<<"pxa 2= "<<particleA->p().px()<<endl;

						cout<<"pya 2= "<<particleA->p().py()<<endl;

						cout<<"pza 2 = "<<particleA->p().pz()<<endl;

					}
					//cout<<"particleA 2 = "<<particleA->p();


				//	cout<<"impactMax = "<<_impactParameterMax<<" impactMin"<<_impactParameterMin<<endl;
					vecCoordinate = particleA->vProd();



					isScattering = pathInNucleus2(particleA,zCoordinateOfCollisions);
				//	cout<<"particleA "<<particleA->vProd();
					cout.precision(12);
			//		cout<<"particleA pt = "<<particleA->vProd().pT()<<endl;

				//	cout<<" z coord of col = "<<zCoordinateOfCollisions<<" isScattering "<<isScattering<<endl;
				//	cin>>ch;
					if(isScattering){
//						cout<<" zCoordinateOfCollisions = "<<zCoordinateOfCollisions<<endl;
						ch = 'a';
					}

					// for initial particle remember coordinate of first collision
					if(_firstCall){
						if(!isScattering)return;
						vecCoordinate.p(particleA->vProd());
						vecCoordinate.pz(zCoordinateOfCollisions);
						_initialParticle.vProd(vecCoordinate);
						_initialParticle.getHistory()->clear();
						_initialParticle.getHistory()->push_back(particleA->getHistory()->back());
					//	cout<<"partA size "<<particleA->getHistory()->size()<<endl;
					//	cout<<"history = "<<particleA->getHistory()->back()<<endl;
						//_initialParticle.vProd(particleA->vProd());
						//_initialParticle.vProd().pz(zCoordinateOfCollisions);
					}else{
					//	vecCoordinate.p(particleA->vProd()); //may be not nesesary
						//cout<<"zCoordinateOfCollisions "<<zCoordinateOfCollisions<<endl;
						if(isScattering)vecCoordinate.pz(zCoordinateOfCollisions);
					}

					if(isScattering){
						//if interaction occurred beyond boundary of nucleus hardron does not lose energy
						//if gotten coordinate above nuclear radius we assume that particle loosing energy until it reach nuclear radius
						zCoordinateOfCollisionsTemp = (zCoordinateOfCollisions >= _targetNucleus.getNuclearRadius()) ? _targetNucleus.getNuclearRadius() : zCoordinateOfCollisions;



						// energy loss will occur up to zCoordinateOfCollisionsTemp
						if(!energyLoss(particleA,zCoordinateOfCollisionsTemp))continue;
						//if(!energyLoss(particleA,zCoordinateOfCollisions))continue;
						//////////////////////////////////////////////////////////

						// set z coordinate of particle

						//particleA->vProd().pz(zCoordinateOfCollisions);

						particleA->vProd(vecCoordinate);
						//cout<<"coord a "<<particleA->vProd();


						 if(!checkEnergyCut(particleA))continue;


					}else{
						//cout<<" we have not collision "<<endl;


						isNotAdsorbed = energyLoss(particleA,_targetNucleus.getNuclearRadius());
						// isNotAdsorbed = 0 - particle is adsorbed
						// isNotAdsorbed = 1 - all right
						if(isNotAdsorbed){
							particleA->setOut();
							//dispose momentum of particle along initial beam direction




							particleA->rotateHardping();



							_outOfNucleus->push_back(*particleA);
							// massive with index of particles which not initialize next wave
							_indexBadInitializations->push_back(particleA->getIndexNumber());
						}
						continue;
					}



//					cout<<particleA->vProd();

					if(true/*numberOfGeneration*/){


						if(_verbose)cout<<"2 coordinate before rotation = "<<particleA->vProd()<<endl;
						if(_verbose)cout<<"2 momentum  before  rotate   = "<<particleA->p()<<endl;

						//dispose momentum of particle along initial beam direction
						particleA->rotateHardping();
						////////////////////////////////////////////////////////////

						if(_verbose)cout<<"2 coordinate after rotation_ = "<<particleA->vProd()<<endl;
						if(_verbose)cout<<"2 momentum  after rotate_    = "<<particleA->p()<<endl;


					}
					particleB->id(getIdTargetNucleon());
					particleA->scatteringOnParticle(particleB->id());



					if( false ){

						particleA->setHard();
						particleA->setSoft(false);
						coordinateHardOutput<<particleA->vProd();
					}else{
						particleA->setSoft();
						particleA->setHard(false);

						coordinateSoftOutput<<particleA->vProd();
					//	cout<<"impact parameter after writing = "<<particleA->vProd().pT()<<endl;
						//cin>>ch;
					}

				}else{//end of if (!_softToHard)
					softToHard(particleA, particleB);
				}// end of else (!_softToHard)

			//	cout<<"after pA 2"<<endl;
				if(particleA->isHard()){
				//	particleA->rotateBackHardping();
					//particleA->rotateHardping();
					if(_verbose)cout<<"before intitialization particle A momentum = "<<particleA->p();
					if(_verbose)cout<<" coordinate                                = "<<particleA->vProd();
				//1	cin>>ch;
				}

				if(_verbose)cout<<"momentum before rotate "<<particleA->p();

				if(!_softToHard&&_fortranHardping){


					//	if(1||_verbose)cout<<"momentum after 0 rotate "<<particleA->p();
						getNewPtInitialState(particleA,1);
					//	if(1||_verbose)cout<<"momentum after 1 rotate "<<particleA->p();
						getNewPtInitialState(particleA,2);
					//	if(1||_verbose)cout<<"momentum after 2 rotate "<<particleA->p();



				}
				particleA->setAngles();
			//	cout<<"momentum after rotate "<<particleA->p();

				if(_firstCall){
					_initialParticle.scatteringOnParticle(particleB->id());
					_initialParticle.p(particleA->p());
					_initialParticle.id(particleA->id());
				}
				if(!pythiaInitialization(particleA ,particleB)){

					if(_verbose)cout << "bad pythia initialization"<<particleA->p()<<" id  = "<<particleA->id()<<endl;

				}else{

					if(_verbose)cout << "good pythia initialization = "<<particleA->p()<<" id  = "<<particleA->id()<<endl;
				//	cin>>ch;

				}


				if(particleA->isSoft() ){
					if(_verbose)cout<<" iam in particleA->isSoft() && _fortranHardping"<<endl;
					if(_fortranHardping){
						pythiaNextFlag = 1;
					}else{
						pythiaNextFlag = pythia->next();
					//	pythia->info.list();
						if(_verbose)cout<<"name of processes "<< pythia->info.name()<<endl;
						if(_verbose)cout<<"code of processes "<< pythia->info.code()<<endl;
						//pythia->event.list();
					}



				}else{
					pythiaNextFlag = pythia->next();
					particleA->setXBjorkenProjectile( pythia->info.x1());
					//pythiaNextFlag = 1;
				}
				if(_verbose)cout<<"name of processes "<< pythia->info.name()<<endl;
				if(_verbose)cout<<"code of processes "<< pythia->info.code()<<endl;
				if(_verbose)pythia->event.list();// don't deleting processes f fbar -> gamma*/Z0 after re initialization in soft collision


				if(_verbose)cout<<"after next "<<particleA->p();
			//	cin>>ch;
				if(!pythiaNextFlag){
					notPythiaNext(particleA);
					continue;
				}// end of if(!pythia->next()


				//////////////////////////////////////////////////////
				// part - if pythia can initialize and run correctly//
				//////////////////////////////////////////////////////


				if(particleA->isHard()){
					_hardInteractionCount++;
					//cout<<"particleA history size = "<<particleA->getHistory()->size()<<endl;
					if(_verbose)cout<<" particleA = "<<particleA->getHistory()->back()<<" have a hard collision "<<endl;



/*
					vecMomentum = particleA->p();

			//		cout<<"mom before "<<particleA->p();
					vecMomentum.pz(vecMomentum.pz()*x1);
					vecMomentum.e(sqrt(vecMomentum.pT2()+vecMomentum.pz()*vecMomentum.pz()+protonMass*protonMass));
					particleA->p(vecMomentum);
			//		cout<<"mom after "<<particleA->p();
					particleA->setAngles();
*/

				}
				if(particleA->isSoft()){
					_softInteractionCount++;
					//cout<<"particleA history size = "<<particleA->getHistory()->size()<<endl;
					if(_verbose)cout<<" particleA = "<<particleA->getHistory()->back()<<" have a soft collision "<<endl;
				}

				///////////////////////////////////////////////////////////
				//  part for replacement soft collisions to hard,  ////////
				//and deleting particles correspond to the soft collisions/
				if(softToHardFlag)deleteParticleFromSoftCollision();
				///////////////////////////////////////////////////////////

				//save energetic and space particle parameters after energy loss, before new collision

//8888888888888888888888 this is for hard scattering 88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888u7

				saveParticle4VectorsAfterEnergyLoss(particleA);
				cout.precision(12);
				//particleA->rotateBackHardping();
	//			cout<<"after saveParticle4VectorsAfterEnergyLoss px "<<particleA->px()<<endl;
		//		cout<<"after saveParticle4VectorsAfterEnergyLoss py "<<particleA->py()<<endl;
			//	cout<<"after saveParticle4VectorsAfterEnergyLoss pz "<<particleA->pz()<<endl;
				if(_verbose)cout << " pythia after save"<<endl;
			//	cout<<particleA->p();
				//cin>>ch;
//888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
				if(_verbose)cout << " pythia cycle begin"<<endl;
			//	int particleIndex = particleA->getHistory()->back();
				for (int i_pyEv = 0; i_pyEv < pythia->event.size() - 3; i_pyEv++){
					if(_verbose)cout << "in pythia cycle "<<i_pyEv<<endl;

/////////////////////////   Drell- Yan part    ////////////////////////////////////////////////////////
					if(findDrellYanPairs( i_pyEv, particleA))return;

				//	cout<<"after findDrellYanPairs px "<<particleA->px()<<endl;
				//	cout<<"after findDrellYanPairs py "<<particleA->py()<<endl;
				//	cout<<"after findDrellYanPairs pz "<<particleA->pz()<<endl;
///////////////////////////////////////////////////////////////////////////////////////////////////////

					prepareNewGeneration(particleA, i_pyEv);

				}//end of for (unsigned int i_pyEv = 0; i_pyEv < pythia->event.size() - 3; i_pyEv++) pythia cycle

				if(!addVertexOfInteraction(particleA))continue;

				if(_verbose)cout <<"pythia->event.size() = "<<pythia->event.size() -3<<endl;
				if(_verbose)cout<< "size generation = "<<_generations->at(numberOfGeneration).getMatrix()->at(i_init).size()<<" at i_init = "<<i_init<<endl;

			} // ind i_init cicle
			if(_verbose)cout<< "_generations size before = "<<_generations->at(numberOfGeneration).getMatrix()->size()<<endl;
			if(_verbose)cout<<"_indexBadInitializations->size() = "<<_indexBadInitializations->size()<<endl;

			///// deleting from new generation all particles which can no initialize new wave///////
			deleteBadInitializationParticlesFromNewGeneration(numberOfGeneration);
			////////////////////////////////////////////////////////////////////////////////////////
		//	cin>>lll;
			if(_verbose)cout<< " _generations size after = "<<_generations->at(numberOfGeneration).getMatrix()->size()<<endl;


			// cout all branches
			// here delete zero rabnches
			if(_firstCall)_firstCall = false;




		if(_verbose)cout<< "NUMBER OF GENERATION "<<numberOfGeneration<<endl;
		numberOfGeneration++;

	//	cout<<"before ifNoHardCollisionHappened "<<endl;
		softToHardFlag = ifNoHardCollisionHappened(numberOfGeneration);
	//	cout<<"after ifNoHardCollisionHappened "<<endl;
	//	cin>>ch;
	}while(_generations->at(numberOfGeneration-1).getMatrix()->size()!=0);





	if(_verbose)cout<< "out of huge cycle"<<endl;

	finalOutput();



}

int
Hardping::pythiaInitialization( hardpingParticle * particleA ,hardpingParticle * particleB){
		int idA, idB, pythiaInitializationFlag = 0;
		double pxA , pyA , pzA ,pEA;
		double pxB , pyB , pzB ,pEB;
		char ch;
		int fff;
		Timer time;
		Particle partA;
	//	partA.pz(800.);
	//	pythia->event.reset();
		if(_verbose)cout<<"pythia->event.size() "<< pythia->event.size()<<endl;
		double _phiHardping = 0;
		double _thetaHardping = 0;
		_phiHardping = particleA->phi();
		_thetaHardping = particleA->theta();
		particleA->setAngles();

	//	particleA->setAngles();
	//	cout<<" A particle momentum in initialization1 "<<particleA->p();
		particleA->rotateBackHardping();
	//	cout<<" A particle momentum in initialization2 "<<particleA->p();




		if(_verbose)cout<<" particle momentum before initialization = "<<particleA->p();
		//cin>>ch;

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
		if(_verbose)cout<<" A particle momentum in initialization "<<particleA->p()<<" id B particle = "<<particleB->id()<<endl;
		if(_verbose)cout<<" hard processes "<<particleA->isHard()<<"  soft processes "<<particleA->isSoft()<<endl;

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


/*
			pythia->readString("Next:numberShowLHA = 0");
			pythia->readString("Next:numberShowInfo = 0");
			pythia->readString("Next:numberShowProcess = 0");
			pythia->readString("Next:numberShowEvent = 0");

*/



		}

	//	pythia->readString("PhaseSpace:pTHatMin = 3.0");

	//	pythia->readString("PhaseSpace:mHatMin = 80.");
	//	pythia->readString("PhaseSpace:mHatMax = 120.");


		if(particleA->isHard()){

			// pythia->readString("SoftQCD:all = off");

			//////////	soft            ////////
			  pythia->readString("SoftQCD:nonDiffractive  = off");
			  pythia->readString("SoftQCD:elastic = off");
			  pythia->readString("SoftQCD:singleDiffractive = off");
			  pythia->readString("SoftQCD:doubleDiffractive = off");
			  pythia->readString("SoftQCD:centralDiffractive = off");
			  pythia->readString("SoftQCD:inelastic = off");
			  if(_verbose)cout<<"HARDQCD is "<<particleA->isHard()<<" hard processes turn on"<<endl;
			  pythia->readString("HardQCD:qqbar2qqbarNew = off");
			  pythia->readString("HardQCD:gg2gg = off");
			  pythia->readString("HardQCD:gg2qqbar = off");
			  pythia->readString("HardQCD:qg2qg = off");
			  pythia->readString("HardQCD:qq2qq = off");
			  pythia->readString("HardQCD:qqbar2gg = off");

			  pythia->readString("WeakSingleBoson:ffbar2gmZ = on");
			  pythia->readString("23: onMode = off");
			  pythia->readString("22: onMode = off");
			  pythia->readString("23: onIfAny = 13 -13");
			  pythia->readString("22: onIfAny = 13 -13");
				if(_cutMass){
					pythia->readString("23:mMin = 4.0");
					pythia->readString("23:mMax = 8.4");
				}


/*

  			pythia->readString("WeakSingleBoson:ffbar2gmZ = on");
  			pythia->readString("23: onMode = off");
  			pythia->readString("22: onMode = off");
  			pythia->readString("23: onIfAny = 13 -13");
  			pythia->readString("22: onIfAny = 13 -13");
  	  	   // pythia->readString("PhaseSpace:mHatMin = 4.0");
  	  	   // pythia->readString("PhaseSpace:mHatMax = 8.4");
  	  	    pythia->readString("23:mMin = 4.0");
  	  	    pythia->readString("23:mMax = 8.4");
*/
		}
		if(particleA->isSoft()){
			  //pythia->readString("HardQCD:all = off");
			  if(_verbose)cout<<"SoftQCD is "<<particleA->isSoft()<<" soft processes turn on "<<endl;



		 	  pythia->readString("SoftQCD:nonDiffractive  = off");
		 	  //	if(particleA->isHard())pythia->readString("SoftQCD:nonDiffractive = on");
		 	  pythia->readString("SoftQCD:elastic = on");
			  pythia->readString("SoftQCD:singleDiffractive = off");
			  pythia->readString("SoftQCD:doubleDiffractive = off");
		 	  pythia->readString("SoftQCD:centralDiffractive = off");
		 	  pythia->readString("SoftQCD:inelastic = off");


			  pythia->readString("HardQCD:qqbar2qqbarNew = off");
			  pythia->readString("HardQCD:gg2gg = off");
			  pythia->readString("HardQCD:gg2qqbar = off");
			  pythia->readString("HardQCD:qg2qg = off");
			  pythia->readString("HardQCD:qq2qq = off");
			  pythia->readString("HardQCD:qqbar2gg = off");

			  pythia->readString("WeakSingleBoson:ffbar2gmZ = off");




			  if(_verbose)cout<<"pythia->event.size() "<< pythia->event.size()<<endl;

			//	getNewPtInitialState(particleA,1);
			//	getNewPtInitialState(particleA,2);

			  	//particleA->rotateHardping();
		//	  cout<<"_fortranHardping "<< _fortranHardping<<endl;
			  if(_fortranHardping){

					partA.p(particleA->p());
					//partA.vProd().p(particleA->vProd());
					partA.id(particleA->id());
					partA.vProd(particleA->vProd());
					//particleA->rotateBackHardping();
					//partA.theta();

					//partA.phi();
					partA.status(1);

					pythia->event.clear();
					pythia->event.append(0);
					pythia->event.append(0);
					pythia->event.append(0);
					fff = pythia->event.append(partA);
					if(_verbose)	cout<<" pz = "<< pythia->event.at(fff).pz()<<endl;
					if(_verbose)cout<<"pythia->event.size() ="<< pythia->event.size()<<endl;

			  }
		//		cin>>ch;
		}
	//	pythia->settings.mode("tune:pp",3);


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

		if(/*_fortranHardping && */particleA->isHard()){
			if(_verbose) cout<<"hard fortran collision occured "<<endl;
		//	 time.start();
			 //pythiaInitializationFlag = pythia->init(idA,idB,pxA,pyA,pzA,pxB,pyB,pzB);

			 pythiaInitializationFlag = pythia->init(idA,idB,pEA,pEB);
		//	 pythiaInitializationFlag = 1;

		//	 cout<<"time of pythia initializations"<<endl;
		//	 time.printTime(time.stop());
		//	 cin>>ch;

		}else{

			if(_verbose)cout<<"soft fortran collision occured "<<endl;
			if(_fortranHardping){
				pythiaInitializationFlag = 1;
			}else{
				pythiaInitializationFlag = pythia->init(idA,idB,pEA,pEB);
			}


		}




		return pythiaInitializationFlag;


	}
bool Hardping::prepareNewGeneration(hardpingParticle* particleA,int i_pyEv){
		char ch;
		//cout<<" in prepareNewGeneration beagin"<<endl;
		//cin>>ch;
		Vec4  vecCoordinate(particleA->vProd());

		hardpingParticle* tempHardpingParticle;

		if(pythia->event.at(i_pyEv+3).isFinal()){
			unsigned int numberOfGeneration = particleA->getNumberOfCurrentGeneration();
			unsigned int i_init =  particleA->getIndexNumber();

			if (/*false*/pythia->event.at(i_pyEv+3).e() < _energyCut) {

				//_tempParticle = new hardpingParticle(pythia->event.at(i_pyEv+3));
				//_cutMassive->push_back(*_tempParticle);

				//_tempParticle = new hardpingParticle(pythia->event.at(i_pyEv+3));
			/*
				_cutMassive->push_back(pythia->event.at(i_pyEv+3));
				_cutMassive->back().getHistory()->clear();
				_cutMassive->back().getHistory()->assign(particleA->getHistory()->begin(),particleA->getHistory()->end());
				_cutMassive->back().getHistory()->push_back(_indexParticle);
				cout<< "index of cut off particles =  ";
				for(unsigned int i = 0; i < _cutMassive->back().getHistory()->size(); i++){
						cout<<_cutMassive->back().getHistory()->at(i)<<" ";
				}
				cout<<endl;
				_indexParticle++;
				*/

											//if(pythia->event.at(i_pyEv+3).statusHepMC()==2 || pythia->event.at(i_pyEv+3).statusHepMC()==1)cout<<"final = "<<pythia->event.at(i_pyEv+3).e()<<" id = "<<pythia->event.at(i_pyEv+3).id()<<endl;
			}else{
							//_generations->at(numberOfGeneration).getMatrix()->at(i_init).push_back(pythia->event.at(i_pyEv+3));// = *particle;
				tempHardpingParticle = new hardpingParticle(pythia->event.at(i_pyEv+3));

				tempHardpingParticle->setPhiHardping(particleA->getPhiHardping());
				tempHardpingParticle->setThetaHardping(particleA->getThetaHardping());
				//.todo add coordinate
				//cout<<"before "<<tempHardpingParticle->p();

				tempHardpingParticle->rotateHardping();
				tempHardpingParticle->vProd(vecCoordinate);
				//cout<<"after "<<tempHardpingParticle->p();
				//cout<<"our part "<<endl;
				//particleA->rotateBackHardping();
				//cout<<"part a "<<particleA->p();

				//cout<<tempHardpingParticle->p();
				//cout<<tempHardpingParticle->vProd();
				//cout<<tempHardpingParticle->getPhiHardping()<<endl;
				//cout<<tempHardpingParticle->getThetaHardping()<<endl;

				tempHardpingParticle->getHistory()->assign(particleA->getHistory()->begin(),particleA->getHistory()->end());
				tempHardpingParticle->getHistory()->push_back(_indexParticle);
				if(/*pythia->event.size() != 4*/1){
					//cout<<"here i am "<<endl;



		//			cout<<"jhdsj "<<tempHardpingParticle->p();
			//		tempHardpingParticle->rotateHardping();
					tempHardpingParticle->setAngles();



					//cin>>ch;
				}
				/*
				tempHardpingParticle->rotateBackHardping();
				cout<<tempHardpingParticle->p();
				cout<<tempHardpingParticle->vProd();
				cout<<tempHardpingParticle->getPhiHardping()<<endl;
				cout<<"history of particle ";
				for(unsigned int i = 0 ; i < tempHardpingParticle->getHistory()->size();i++){
					cout<<	tempHardpingParticle->getHistory()->at(i)<<" ";
				}
				cout<<endl;
				tempHardpingParticle->setAngles();
				cout<<tempHardpingParticle->getThetaHardping()<<endl;
				//tempHardpingParticle->rotateBackHardping();
				cout<<tempHardpingParticle->p();

				cin>>ch;
				*/
				if(_verbose){
					cout<<"momentum of particle in prepare new gen "<<tempHardpingParticle->p();//<<endl;
					cout<<"coord of particle in prepare new gen "<<tempHardpingParticle->vProd();//<<endl;
					cout<<"history of particle ";
					for(unsigned int i = 0 ; i < tempHardpingParticle->getHistory()->size();i++){
						cout<<	tempHardpingParticle->getHistory()->at(i)<<" ";
					}
					cout<<endl;
		//			cin>>ch;
				}


			//	cout<<" 111"<< tempHardpingParticle->p();
				_generations->at(numberOfGeneration).getMatrix()->at(i_init).push_back(*tempHardpingParticle);// = *particle;

			//	_generations->at(numberOfGeneration).getMatrix()->at(i_init).back().id(particleA->id());
				if(_verbose)cout<<" producted particle "<<_generations->at(numberOfGeneration).getMatrix()->at(i_init).back().id()<<" at index "<<_generations->at(numberOfGeneration).getMatrix()->at(i_init).back().getHistory()->back()<<endl;

				_indexParticle++;
				//todo may be save all hardpingParticle characteristics
			/*	phi = _generations->at(numberOfGeneration).getMatrix()->at(i_init).back().phi();
				theta = _generations->at(numberOfGeneration).getMatrix()->at(i_init).back().theta();
				cout<<" coord 1=   "<<_generations->at(numberOfGeneration).getMatrix()->at(i_init).back().vProd()<<endl;
				cout<<"icoord 1=   "<<_generations->at(numberOfGeneration).getMatrix()->at(i_init).back().p();
				_generations->at(numberOfGeneration).getMatrix()->at(i_init).back().rot(0,-phi);
				_generations->at(numberOfGeneration).getMatrix()->at(i_init).back().rot(-theta,0);
				cout<<"icoord 2=   "<<_generations->at(numberOfGeneration).getMatrix()->at(i_init).back().p()<<endl;
				cout<<" coord 2=   "<<_generations->at(numberOfGeneration).getMatrix()->at(i_init).back().vProd()<<endl;*/
				delete tempHardpingParticle;
				return true;
			}// end else (pythia->event.at(i_pyEv+3).e() < _energyCut)

		}// end if(pythia->event.at(i_pyEv+3).isFinal())
		//cout<<" in prepareNewGeneration end"<<endl;
		//cin>>ch;
		return false;
	}
bool Hardping::findDrellYanPairs(int i_pyEv, hardpingParticle* particleA){
		char ch;
		//cout<<" in findDrellYanPairs beagin"<<endl;
		//cin>>ch;
		Vec4 vecMomentum(0);
		hardpingParticle* tempHardpingParticle;
		tempHardpingParticle = new hardpingParticle();
		tempHardpingParticle->setPhiHardping(particleA->getPhiHardping());
		tempHardpingParticle->setThetaHardping(particleA->getThetaHardping());
		double pt_old = 0, pt_new = 0, px_old =0, py_old = 0,pz_old =0,pe_old =0, px_new = 0, py_new =0, deltaPt = 0;

		Vec4 vecMomentumZ0(0);



		px_old = particleA->px();
		py_old = particleA->py();
		pz_old = particleA->pz();
		cout.precision(10);
	//	cout<<"findDrellYanPairs px  "<<px_old<<" findDrellYanPairs py  = "<<py_old<<" findDrellYanPairsp pz  = "<<pz_old<<endl;

		//cout << "particleA->theta() "<<particleA->theta()<<" particleA->phi() "<<particleA->phi()<<endl;

		double ptMother = particleA->pT();
	//	cout<<"phi = "<<tempHardpingParticle->getPhiHardping()<<" theta = "<<tempHardpingParticle->getThetaHardping()<<endl;
		if(abs(pythia->event.at(i_pyEv+3).id()) == 23 ){
			//cout<<"Z0 = "<<pythia->event.at(i_pyEv+3).p();

		//	cout<<"_softInteractionCount = "<<_softInteractionCount<<endl;

/*
			for(unsigned int i = 0; i <= _softInteractionCount; i++){
			//	cout<<"incident particle1 "<<particleA->p();
				getNewPtInitialState(particleA,1);
			//	cout<<"incident particle2 "<<particleA->p();
				getNewPtInitialState(particleA,2);
			//	cout<<"incident particle3 "<<particleA->p();
			}

*/

		//	cout<<"incident particle3 "<<particleA->p();
			vecMomentum = particleA->p();
		//	cout<<"x1 "<<particleA->getXBjorkenProjectile()<<endl;
	//		cout<<"mom before "<<particleA->p();
			vecMomentum.pz(vecMomentum.pz()*particleA->getXBjorkenProjectile());
			vecMomentum.e(sqrt(vecMomentum.pT2()+vecMomentum.pz()*vecMomentum.pz()+protonMass*protonMass));
			particleA->p(vecMomentum);
	//		cout<<"mom after "<<particleA->p();
			particleA->setAngles();
	//		cout<<"phi = "<<particleA->getPhiHardping()<<" theta = "<<particleA->getThetaHardping()<<endl;
/*
			if(_targetNucleus.A()==184){
				double tempZ =sqrt(particleA->pz()*particleA->pz()-5*5);
				Vec4 Wmomentum(5.,0.,tempZ,particleA->e());
				particleA->p(Wmomentum);
			}
*/
			px_old = particleA->px();
			py_old = particleA->py();
			pz_old = particleA->pz();
			//cout.precision(10);
	//		cout<<"px A = "<<px_old<<" py A = "<<py_old<<" pz A = "<<pz_old<<endl;


//			cout<<"getting pt "<<particleA->pT()-ptMother<<" pt mother "<<ptMother<<endl;
			//cin>>ch;


//			particleA->setAngles();

			cout.precision(10);
/*			cout<<"phi = "<<particleA->getPhiHardping()<<" theta = "<<particleA->getThetaHardping()<<endl;
			cout<<" cos phi = "<<cos(particleA->getPhiHardping())<<"sin phi = "<<sin(particleA->getPhiHardping())<<endl;
			cout<<" cos theta = "<<cos(particleA->getThetaHardping())<<"sin theta = "<<sin(particleA->getThetaHardping())<<endl;
*/
			//cin>>ch;
		/*	cout<<"p = "<<particleA->p();
			particleA->rotateBackHardping();
			cout<<"p2 = "<<particleA->p();
			particleA->rotateHardping();
			cout<<"p3 = "<<particleA->p();*/
	//		cin>>ch;

		//	cout<<"_softInteractionCount = "<<_softInteractionCount<<endl;
			if(_verbose)cout<<"Z0 WAS BORNED !!!"<<endl;
			/*int id = pythia->event.at(pythia->event.at(i_pyEv+3).mother1()).id();
			if(id == 23){
				cout<<"true leptone"<<endl;
			}else{
				id = pythia->event.at(id).mother1();
			}*/
			int motherLine1 = i_pyEv + 3;
			int  idDaughter1, idDaughter2,lineDaughter1, lineDaughter12, lineDaughter2,lineDaughter22  , motherLine2;


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

				px_old = pythia->event.at(motherLine1).px() + pythia->event.at(motherLine2).px();
				py_old = pythia->event.at(motherLine1).py() + pythia->event.at(motherLine2).py();
				pz_old = pythia->event.at(motherLine1).pz() + pythia->event.at(motherLine2).pz();
				pe_old = pythia->event.at(motherLine1).e()  + pythia->event.at(motherLine2).e();



			//	pz_old = sqrt(px_old*px_old +py_old*py_old);

				vecMomentumZ0.px(px_old);
				vecMomentumZ0.py(py_old);
				vecMomentumZ0.pz(pz_old);
				vecMomentumZ0.e(pe_old);

				tempHardpingParticle->p(vecMomentumZ0);
				tempHardpingParticle->setPhiHardping(particleA->getPhiHardping());
				tempHardpingParticle->setThetaHardping(particleA->getThetaHardping());

			//	cout<<"phi = "<<tempHardpingParticle->getPhiHardping()<<" theta = "<<tempHardpingParticle->getThetaHardping()<<endl;
















	//			cout<<" px_old = "<<px_old<<" py_old = "<<py_old<<" pz_old = "<<pz_old<<endl;
			//	cin>>ch;
				//*tempHardpingParticle = pythia->event.at(motherLine1);



/*
				cout.precision(10);
				px_old = pythia->event.at(motherLine1).px();
				py_old = pythia->event.at(motherLine1).py() ;
				pz_old = pythia->event.at(motherLine1).pz() ;
	//			cout<<"initial state first lepton px = "<<px_old<<" py = "<<py_old<<" pz = "<<pz_old<<endl;
				//cout<<"temp p1 ="<<tempHardpingParticle->p();


				tempHardpingParticle->setPhiHardping(particleA->getPhiHardping());
				tempHardpingParticle->setThetaHardping(particleA->getThetaHardping());
			//	tempHardpingParticle->setAngles();
				//cout<<tempHardpingParticle->p()<<tempHardpingParticle->pT()<<endl;
				cout.precision(10);
*/
/*
				cout<<"phi = "<<tempHardpingParticle->getPhiHardping()<<endl;
				cout<<"theta = "<<tempHardpingParticle->getThetaHardping()<<endl;
				cout<<"cos phi lep "<<cos(tempHardpingParticle->getPhiHardping())<<endl;
				cout<<"cos theta lep "<<cos(tempHardpingParticle->getThetaHardping())<<endl;
				cout<<"sin phi lep "<<sin(tempHardpingParticle->getPhiHardping())<<endl;
				cout<<"sin theta lep "<<sin(tempHardpingParticle->getThetaHardping())<<endl;
*/







				tempHardpingParticle->rotateHardping();


				cout.precision(10);
				px_old = tempHardpingParticle->px();
				py_old = tempHardpingParticle->py() ;
				pz_old = tempHardpingParticle->pz() ;


/*
				cout<<"px2 = "<<px_old<<" py2 = "<<py_old<<" pz2 = "<<pz_old<<endl;
				cout<<"first lepton pt = "<<tempHardpingParticle->pT()<<endl;
				cout<<" theta1 "<<tempHardpingParticle->getThetaHardping()<<" phi1 "<<tempHardpingParticle->getPhiHardping()<<endl;

*/

/*				px_new = tempHardpingParticle->px();
				py_new = tempHardpingParticle->py();
*/
				_finalState->push_back(*tempHardpingParticle);
		//		cout<<"111111111111111"<<tempHardpingParticle->p();

//				cout<<"final state size = "<<_finalState->back().p();
			//	cin>>ch;
	/*			*tempHardpingParticle = pythia->event.at(motherLine2);
				cout.precision(10);
				px_old = tempHardpingParticle->px();
				py_old = tempHardpingParticle->py() ;
				pz_old = tempHardpingParticle->pz() ;
	//			cout<<"initial state second lepton  px = "<<px_old<<" py = "<<py_old<<" pz = "<<pz_old<<endl;
			//	cout<<"temp p2 ="<<tempHardpingParticle->p();
				tempHardpingParticle->setPhiHardping(particleA->getPhiHardping());
				tempHardpingParticle->setThetaHardping(particleA->getThetaHardping());





				tempHardpingParticle->rotateHardping();


				cout.precision(10);
				px_old = tempHardpingParticle->px();
				py_old = tempHardpingParticle->py() ;
				pz_old = tempHardpingParticle->pz() ;

		//	cout<<"second lepton = "<<"px2 = "<<px_old<<" py2 = "<<py_old<<" pz2 = "<<pz_old<<endl;
		//	cout<<"second lepton pt = "<<tempHardpingParticle->pT()<<endl;


		//		cin>>ch;
				//cout<<" theta "<<tempHardpingParticle->getThetaHardping()<<" phi "<<tempHardpingParticle->getPhiHardping()<<endl;
		//		cin>>ch;



				//getNewPtInitialState(particleA,1);

			//	getNewPtInitialState(particleA,2);










				px_new += tempHardpingParticle->px();
				py_new += tempHardpingParticle->py();
				pt_new = sqrt(px_new*px_new +py_new*py_new);
				_finalState->push_back(*tempHardpingParticle);
				deltaPt = pt_new - pt_old;
				//cout<<"deltaPt = "<<deltaPt<<endl;
				deltaPtOutput<<deltaPt<<endl;
				//_finalState->push_back(pythia->event.at(motherLine1));
				//_finalState->push_back(pythia->event.at(motherLine2));
				*/
			}else{
				if(_verbose)cout<<" pythia->event.at(motherLine1) "<<pythia->event.at(motherLine1).id()<<" dauther = "<<pythia->event.at(motherLine1).daughter1()<<endl;
				if(_verbose)cout<<" pythia->event.at(motherLine2) "<<pythia->event.at(motherLine2).id()<<" dauther = "<<pythia->event.at(motherLine2).daughter1()<<endl;
				//cin>>lll;
			}


//		write in file total number of collisions
			//change 22.07.14
		//	softCollisionsNumberOutput<<_softInteractionCount +_hardInteractionCount<<endl;
			//change 22.07.14


		//	cout<<"num of col = "<<_softInteractionCount +_hardInteractionCount<<endl;


			//cin>>lll;
		//	cout<<" in findDrellYanPairs end"<<endl;
		//	cin>>ch;
			delete tempHardpingParticle;
			return true;

		}// end of  Drell- Yan part

		if(_hardInteractionCount){
		//	softCollisionsNumberOutput<<_softInteractionCount<<endl;

		//	return true;
		}
		return false;
	}
