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


#ifndef NUCLEUS_H
#define NUCLEUS_H

#define protonMass 0.93827204621
//#include <cmath>
#include "Pythia.h"
#include "timer.h"
#include <vector>
extern ofstream pathInNucleiOutput;
extern ofstream softCollisionsNumberOutput;
extern ofstream deltaPtOutput;
extern  double getRandomFromFile();
 //class Pythia8::Pythia;
//This class holds the information for a target nucleus
namespace Pythia8{
struct index {
  int i;
  int j;
} ;
//Pythia* pythia;
class hardpingParticle : public Particle{
public:
	hardpingParticle():Particle(0){
		_history = NULL;
		_hardInteraction = false;
		_softInteraction = false;
		_outOfNucleus = false;
		_scatteringOnparticle = 0;
		_history = new vector <unsigned int>;
		_numberOfCurrentGeneration = 0;
		_indexInGeneration = 0;
		_thetaHardping = 0;
		_phiHardping = 0;



	}
	hardpingParticle(Particle &particle):Particle(particle){
		_history = NULL;
		_hardInteraction = false;
		_softInteraction = false;
		_outOfNucleus = false;
		_scatteringOnparticle = 0;
		_history = new vector <unsigned int>;
		_numberOfCurrentGeneration = 0;
		_indexInGeneration = 0;
		_thetaHardping = 0;
		_phiHardping = 0;
	}
	vector <unsigned int> * getHistory(void){
		return _history;
	}
	bool isHard(void){	return _hardInteraction;}
	bool isSoft(void){	return _softInteraction;}
	bool isOut(void){	return _outOfNucleus;}
	void setSoft(bool b = true){ _softInteraction = b; }
	void setHard(bool b = true){ _hardInteraction = b; }
	void setOut (bool b = true){ _outOfNucleus = b; }
	void scatteringOnParticle(int id){ _scatteringOnparticle = id ;}
	int  getIdscatteringParticle(void){ return _scatteringOnparticle;}
	double getMaxTransverseMomentum(int type);
	double getNewPtInitialState(int type);
	void setThetaHardping(double theta){_thetaHardping = theta;}
	void setPhiHardping(double phi){_phiHardping = phi;}
	double getThetaHardping(){return _thetaHardping;}
	double getPhiHardping(){return _phiHardping;}

	void setAngles(void){
		_thetaHardping = this->theta();
		_phiHardping = this->phi();
		char ch;
/*
		_thetaHardping = acos(this->pz()/this->pAbs());
		if(this->pT() == 0){
			_phiHardping = 0;
		}else{
			_phiHardping   = acos(this->py()/this->pT());
		//	cout<<"_phiHardping = "<<_phiHardping<<endl;
		//	cout<<"asin phi"<<asin(this->px()/this->pT())+M_PIl<<endl;
		}
*/
	/*	cout.precision(18);
		cout<<"this->theta() = "<<this->theta()<<" this->phi() = "<<this->phi()<<endl;
		cout<<"thetaHardping = "<<_thetaHardping<<" phi = "<<_phiHardping<<endl;
		cout<<"sin1 "<<sin(this->phi())<<" sin2 "<<sin(_phiHardping+M_PIl/2.)<<endl;
		cout<<"cos1 "<<cos(this->phi())<<" cos2 "<<cos(_phiHardping+M_PIl/2.)<<endl;

		cout<<"tsin1 "<<sin(this->theta())<<" tsin2 "<<sin(_thetaHardping)<<endl;
		cout<<"tcos1 "<<cos(this->theta())<<" tcos2 "<<cos(_thetaHardping)<<endl;
	*/	//cin>>ch;

	}
	void rotateHardping(void){
		//if(_verbose)cout<<"thetaHardping = "<<_thetaHardping<<" phi = "<<_phiHardping<<endl;
		//dispose momentum of particle along initial beam direction
		this->rot(_thetaHardping,0);
		this->rot(0,_phiHardping);
	//	cout<<"_thetaHardping = "<<_thetaHardping<<endl;
	//	cout<<"_phiHardping = "<<_phiHardping<<endl;
		//this->p().rotHardpingTest(_thetaHardping,0);
		//this->p().rotHardpingTest(0,_phiHardping);
	//	cout<<"0 rotation "<<this->p();
		//rotateAroundX(_thetaHardping);
	//	cout<<"1 rotation "<<this->p();
		//rotateAroundZ(_phiHardping);
	//	cout<<"2 rotation "<<this->p();
	//	this->rot(0,_thetaHardping);
	//	this->rot(_phiHardping,0);
	}
	void rotateBackHardping(void){
		//if(_verbose)cout<<"theta = "<<_thetaHardping<<" phi = "<<_phiHardping<<endl;
		//set particle momentum along z' direction
		this->rot(0,-_phiHardping);
		this->rot(-_thetaHardping,0);
	//	rotateAroundZ(-_phiHardping);
	//	rotateAroundX(-_thetaHardping);

		//this->p().rotHardpingTest(0,-_phiHardping);
		//this->p().rotHardpingTest(-_thetaHardping,0);

//		this->rot(-_phiHardping,0);
//		this->rot(0,-_thetaHardping);

	}



	void rotateAroundX(double thetaIn, double phiIn = 0){
		// phi assumed equal zero//
	  double xx = this->px();
	  double yy = this->py();
	  double zz = this->pz();
	  double tt = this->e();

//	  cout<<"xx = "<<xx<<" yy = "<<yy<<" zz = "<<zz<<endl;
	  double cthe = cos(thetaIn);
	  double sthe = sin(thetaIn);
	  double cphi = cos(phiIn);
	  double sphi = sin(phiIn);
	  double tmpx =   cphi* xx  +  cthe* sphi * yy + sthe * sphi * zz;
	  double tmpy =  -sphi* xx  +  cphi* cthe * yy + sthe * cphi * zz;
	  double tmpz =  -sthe* yy  +  cthe* zz;
	  xx          = tmpx;
	  yy          = tmpy;
	  zz          = tmpz;
	  Vec4 tempVector(xx, yy, zz, tt);
	  this->p(tempVector);
	//  cout<<"xx2 = "<<xx<<" yy2 = "<<yy<<" zz2 = "<<zz<<endl;
	}
	void rotateAroundZ(double phiIn, double thetaIn = 0){
		// theta assumed equal zero//
	  double xx = this->px();
	  double yy = this->py();
	  double zz = this->pz();
	  double tt = this->e();

	  //cout<<"xx = "<<xx<<" yy = "<<yy<<" zz = "<<zz<<endl;
	  double cthe = cos(thetaIn);
	  double sthe = sin(thetaIn);
	  double cphi = cos(phiIn);
	  double sphi = sin(phiIn);
	  double tmpx =   cphi* xx  +  cthe* sphi * yy + sthe * sphi * zz;
	  double tmpy =  -sphi* xx  +  cphi* cthe * yy + sthe * cphi * zz;
	  double tmpz =  -sthe* yy  +  cthe* zz;
	  xx          = tmpx;
	  yy          = tmpy;
	  zz          = tmpz;
	  Vec4 tempVector(xx, yy, zz, tt);
	  this->p(tempVector);
	  //cout<<"xx2 = "<<xx<<" yy2 = "<<yy<<" zz2 = "<<zz<<endl;
	}
	void setInexNumber(unsigned int i){
		_indexInGeneration = i;
	}
	unsigned int getIndexNumber(void){
		return _indexInGeneration;
	}
	void setNumberOfCurrentGeneration(unsigned int i){
		_numberOfCurrentGeneration = i;
	}
	unsigned int getNumberOfCurrentGeneration(void){
		return _numberOfCurrentGeneration;
	}
	double getXBjorkenProjectile(void){
		return _x1;
	}
	void setXBjorkenProjectile(double x1){
		 _x1 = x1;
	}
private:
	vector <unsigned int> * _history;
	bool _hardInteraction;
	bool _softInteraction;
	bool _outOfNucleus;
	int _scatteringOnparticle;
	double _thetaHardping;
	double _phiHardping;
	unsigned int _indexInGeneration;// i_init - index number of particle in current generation
	unsigned int _numberOfCurrentGeneration;
	double _x1;
	//Rndm * _random;

};
class nucleus{
public:
	nucleus(const int Z, const int A);
	~nucleus(){};
	int    Z              () const { return _Z;                     }  ///< returns atomic number of nucleus
	int    A              () const { return _A;                     }  ///< returns nucleon number of nucleus

	double getNuclearDensity(const double r) const
	{ return (_A < 10) ? nuclearGaussDensity(r) : nuclearWoodSaxonDensity(r);}
	double renormalizedNuclearThicknessGauss12(double impactParameter, double zMin, double zMax)const;
	double renormalizedNuclearThicknessGauss5(double impactParameter, double zMin, double zMax)const;
	void setPrecisionOfNuclearThicknessCalculation(double precision){_precisionOfNuclearThicknessCalculation = precision;}
	double getNuclearRadius  () const { return _r0 ; }
private:
	double woodSaxonSkinDepth() const { return 0.53;  }  ///< returns surface (0.53 fm for Au)
	double woodSaxonRadius() const { return 1.2 * pow(_A, 1. / 3.); }  ///< returns Wood-Saxon nuclear radius [fm] (Fermi model)

	double nuclearThicknessGauss5 (const double impactParameter, double zMin, double zMax) const;  ///< calculates nuclear thickness function for given distance b in impact parameter space (Eq. 4 in KN, PRC 60)
	double nuclearThicknessGauss12(double impactParameter, double zMin, double zMax)const;

	//	Pythia8::Pythia* _pythia;

	int    _Z;                      ///< atomic number of nucleus
	int    _A;                      ///< nucleon number of nucleus
	double _r0;				// nuclear radius
	double _rho0;
	double _precisionOfNuclearThicknessCalculation;
	double nuclearWoodSaxonDensity(const double r) const
	{ return 3./4.*_A/M_PIl/getNuclearRadius()/getNuclearRadius()/getNuclearRadius()/(1+M_PIl*M_PIl*woodSaxonSkinDepth()*woodSaxonSkinDepth()/getNuclearRadius()/getNuclearRadius())/(1. + exp((r - getNuclearRadius()) / woodSaxonSkinDepth())); } ///< Wood-Saxon nuclear density
	double nuclearGaussDensity(const double r) const
	{ return exp(-r*r/getNuclearRadius()/getNuclearRadius())*_A/M_PIl/getNuclearRadius()/getNuclearRadius()/getNuclearRadius()/2*M_2_SQRTPIl;}
};
class producedParticlesInformation{
public:
	producedParticlesInformation():nProducedParticles(NULL),nParticleInFinalState(NULL),nParticleForNewGeneration(NULL){
		nProducedParticles = new vector <int>;
		nParticleInFinalState = new vector <int>;
		nParticleForNewGeneration  = new vector <int>;

	}
	~producedParticlesInformation(){}
	int getNProducedParticles(int nGeneration){
		return nProducedParticles->at(nGeneration);
	}
	void setNProducedParticles(int nGeneration, int nParticles){
		nProducedParticles->at(nGeneration) = nParticles;
	}
	int getNParticleInFinalState(int nGeneration){
		return nParticleInFinalState->at(nGeneration);
	}
	void setNParticleInFinalState(int nGeneration, int nParticles){
		nParticleInFinalState->at(nGeneration) = nParticles;
	}
	int getNParticleForNewGeneration(int nGeneration){
		return nParticleForNewGeneration->at(nGeneration);
	}
	void setNParticleForNewGeneration(int nGeneration, int nParticles){
		nParticleForNewGeneration->at(nGeneration) = nParticles;
	}
	void increaseNParticleInFinalState(int nGeneration){nParticleInFinalState->at(nGeneration)++;}
	std::vector <int> * getVectorNProducedParticlesPointer(){return nProducedParticles;}
	std::vector <int> * getVectorNParticleInFinalStatePointer(){return nParticleInFinalState;}
	std::vector <int> * getVectorNParticleForNewGenerationPointer(){return nParticleForNewGeneration;}
protected:

	std::vector <int> *nProducedParticles;
	std::vector <int> *nParticleInFinalState;
	std::vector <int> *nParticleForNewGeneration;
};
class vectorHardping : protected Vec4{
public:
	/*vectorHardping(){
		cout<<"in constructor vectorHardping 0 "<<endl;
		Vec4 * pVec = new Vec4(0);
		setVector(*pVec);
		//setVector(pVec);
		_idParticle = 0 ;
		delete pVec;
		cout<<"here!!!"<<endl;
	}*/
	vectorHardping(int id = 0):_idParticle(id),_vector4Hardping(0),_theta(0),_phi(0),_impactParameterX(0),_impactParameterY(0){}
	vectorHardping(Vec4 & vector, int id){
		//cout<<"in constructor vectorHardping 1"<<endl;
		setVector (vector);
		_idParticle = id;
	}
	vectorHardping( Particle * particle):
		_idParticle(particle->id()),
		_phi(particle->phi()),
		_theta(particle->theta()),
		_vector4HardpingProduction(NULL)
	{
		//particle->
	//	cout<<"in constructor vectorHardping 2"<<endl;
		Vec4 * pVec = new Vec4(particle->p());
		setVector(*pVec);
		Vec4 * pVecProd = new Vec4(particle->vProd());
		setVectorProduction(*pVecProd);
		delete pVec;

	}
	~vectorHardping(){
		delete _vector4Hardping;
	}
	void setVector(Vec4 &v){
		//cout<<"in set vector"<<endl;
	//	delete _vectorHardping;
		_vector4Hardping = new Vec4(v);
	}
	void setVectorProduction(Vec4 &v){
		_vector4HardpingProduction = new Vec4(v);
	}
	void setHardpingVector(vectorHardping &v){
	//	cout<<"in setHardpingVector"<<endl;
	//	delete _vectorHardping;
		//v.getVector()
		_vector4Hardping = new Vec4(*(v.getVector()));
		_idParticle = v.getId();
		_impactParameterX = v.getImpactParameterX();
		_impactParameterY = v.getImpactParameterY();
		_theta = v.getTheta();
		_phi = v.getPhi();


	}
	Vec4 * getVector(void){
		return _vector4Hardping;
	}
	void setId(int id){_idParticle = 1;}
	int getId(void){return _idParticle;}
	void setTheta(double theta){_theta = theta;}
	double getTheta(void){return _theta;}
	void setPhi(double phi){_phi = phi;}
	double getPhi(void){return _phi;}
	void setImpactParameter(double bx, double by){
		_impactParameterX = bx;
		_impactParameterY = by;
	}
	double getAbsImpactParameter(void){
		return sqrt(_impactParameterX*_impactParameterX+_impactParameterY*_impactParameterY);
	}
	double getImpactParameterX(void){
		return _impactParameterX;
	}
	double getImpactParameterY(void){
		return _impactParameterY;
	}

//	vector->

private:
	Vec4 * _vector4Hardping;
	Vec4 * _vector4HardpingProduction;
	int _idParticle;
	double _theta;
	double _phi;
	double _impactParameterX;
	double _impactParameterY;

};

class generation{
	typedef std::vector< std::vector<hardpingParticle> > vectorConstruction;
public:
	generation(const int numberOfGeneration, int iSize, int jSize):
	_iSize(iSize),
	_jSize(jSize),
	//_matrix(0),
	_numberOfGeneration(numberOfGeneration)
	{
		/*
		_matrix = vectorConstruction(0);
		if(_matrix.size() == 0){
			cout<< "_matrix is not initialized. Initialization matrix ..."<<endl;
			_matrix = vectorConstruction(_iSize,std::vector<Vec4>(_jSize));
		}else{
			cout<< "matrix already initialized"<<endl;
		}
		cout<<"matrix initialized "<<endl;
		*/
		initializeNewMatrix();
		//cout<<

	};
	generation():
		_iSize(0),
		_jSize(0),
		//_matrix(0),
		_numberOfGeneration(0)
		{
			/*
			_matrix = vectorConstruction(0);
			if(_matrix.size() == 0){
				cout<< "_matrix is not initialized. Initialization matrix ..."<<endl;
				_matrix = vectorConstruction(_iSize,std::vector<Vec4>(_jSize));
			}else{
				cout<< "matrix already initialized"<<endl;
			}
			cout<<"matrix initialized "<<endl;
			*/
			initializeNewMatrix();
			//cout<<

		};
	~generation(){};

	vectorConstruction* initializeNewMatrix(){ //return pointer to the matrix
		//cout<<"matrix size "<<_matrix.size()<<endl;
		if(_matrix.size() == 0){
		//	cout<< "_matrix is not initialized. Initialization matrix ..."<<endl;
			_matrix = vectorConstruction(_iSize,std::vector<hardpingParticle>(_jSize));

		}else{
			cout<< "matrix already initialized"<<endl;
		}

		/*for (unsigned int i =0; i< _matrix.size(); i++){
			for (unsigned int j =0; j< _matrix[i].size(); j++){
				_matrix[i][j] = 0;
			}

		}*/
		//cout<<"matrix initialized "<<endl;

		return &_matrix;
	}
	vectorConstruction* getMatrix(){ //return pointer to the matrix
				return &_matrix;
	}
	int* initialize_int(){
		cout<<"int initialized "<<endl;

		return &_iSize;
	}
	/*
	void printMatrix(void){
		for(int i = 0; i < _matrix.size(); i++){
			for(int j=0;j< _matrix[i].size();j++){
				cout<<(_matrix[i][j])<<" ";
			}
			cout<<endl;
		}
	}
*/
	//vectorConstruction matrix(10,std::vector <int>(10));
	//std::vector< std::vector< int > >  m(1, std::vector<int>(1));
private:
	vectorConstruction _matrix;
	int _iSize;
	int _jSize;
	int _numberOfGeneration;
	Vec4 tempVec4;
};// end class generation

class Pythia;
class Hardping
{
	//Pythia8::Pythia * pythia;

public:

	Hardping(nucleus projectileNucleus, nucleus targetNucleus);

	~Hardping();
 
//	int    Z              () const { return _Z;                     }  ///< returns atomic number of nucleus
//	int    A              () const { return _A;                     }  ///< returns nucleon number of nucleus
//	double woodSaxonRadius() const { return 1.2 * pow(_A, 1. / 3.); }  ///< returns Wood-Saxon nuclear radius [fm] (Fermi model)
/*	double nuclearThicknessGauss5 (const double impactParameter, double zMin, double zMax) const;  ///< calculates nuclear thickness function for given distance b in impact parameter space (Eq. 4 in KN, PRC 60)
	double nuclearThicknessGauss12(double impactParameter, double zMin, double zMax)const;
	double renormalizedNuclearThicknessGauss12(double impactParameter, double zMin, double zMax)const;
	double renormalizedNuclearThicknessGauss5(double impactParameter, double zMin, double zMax)const;*/
	double probabilityOfNSoftScattering(double impactParameter, double zMin, double zMax, int numberOfCollisions)const;
	double integralProbabilityOfNSoftScattering(int numberOfCollisions) const;
	double rho0() const { return _rho0; }
	double getr0() const { return _r0;}
	int factorial(int i)const{return _factorials[i];}
	//void setPrecisionOfNuclearThicknessCalculation(double precision){_precisionOfNuclearThicknessCalculation = precision;}
	void randomGeneratorInitialize(int seed){
		_randomSeed = seed;
		//_random = new Rndm(_randomSeed);
		_random->init(seed);
	}
	double getRandom(void){ return pythia->rndm.flat();}
/*	double nuclearWoodSaxonDensity(const double r) const
	{ return 3./4.*_targetNucleus.A()/M_PIl/getNuclearRadius()/getNuclearRadius()/getNuclearRadius()/(1+M_PIl*M_PIl*woodSaxonSkinDepth()*woodSaxonSkinDepth()/getNuclearRadius()/getNuclearRadius())/(1. + exp((r - getNuclearRadius()) / woodSaxonSkinDepth())); } ///< Wood-Saxon nuclear density
	double nuclearGaussDensity(const double r) const
	{ 	//std::cout<< "i ma in nuclearGaussDensity"<<std::endl;
		return exp(-r*r/getNuclearRadius()/getNuclearRadius())*_targetNucleus.A()/M_PIl/getNuclearRadius()/getNuclearRadius()/getNuclearRadius()/2*M_2_SQRTPIl;}
	double getNuclearDensity(const double r) const
	{ return (_targetNucleus.A() < 10) ? nuclearGaussDensity(r) : nuclearWoodSaxonDensity(r);}*/
	double pathInNucleus(hardpingParticle * particleA);
	bool pathInNucleus2(hardpingParticle * particleA, double &zCoordinateOfCollision);
	Pythia8::Pythia * pythia;
	double getAbsPartonMomentum(void){
		_absPartonMomentum = sqrt(_xPartonMomentum*_xPartonMomentum+_yPartonMomentum*_yPartonMomentum+_zPartonMomentum*_zPartonMomentum);
		return _absPartonMomentum;
	}
	double getNewPtInitialState(hardpingParticle* particleA, int type);
//	Pythia* pythia;
//	Pythia* pythiaHard;
	Pythia8::generation * _generation;
	void initGeteration(int nOfGeneration,int nInitialParcicle){
		_generation = new generation(nOfGeneration,nInitialParcicle,0);
	}
	int pythiaInitialization(hardpingParticle * particleA ,hardpingParticle * particleB);
	void hardping(void);
	double getMaxTransverseMomentum(hardpingParticle* particleA,int type);
	void setEnergy(double Elab){
		_energyLab = Elab;
	}
	unsigned int calculateNumberOfParticlesAtPreviousGeneration(unsigned int numberOfGeneration){
		unsigned int firstWave = 0;
		unsigned int isumMax = 0;
		isumMax = (numberOfGeneration) ? _generations->at(numberOfGeneration-1).getMatrix()->size() : 0;
		if(_verbose)if(!_firstCall)cout<<" number of incident particle in previous generation "<<isumMax<<endl;
		////////////// calculating number of particles at previous generation/////////////////////////
		for(int isum = 0; isum < isumMax ;isum++){
			//cout<<" in summ !!"<<endl;
			if(_verbose)cout<<" in Generation "<<numberOfGeneration-1<<" at i = "<<isum<< " size equal = "<<_generations->at(numberOfGeneration-1).getMatrix()->at(isum).size()<<endl;
			firstWave +=  _generations->at(numberOfGeneration-1).getMatrix()->at(isum).size();

		}
		return firstWave;
	}
	void particleIndexation(unsigned int numberOfGeneration,unsigned int numberOfParticlesAtPreviousGeneration){

		unsigned int isumMax = (numberOfGeneration) ? _generations->at(numberOfGeneration-1).getMatrix()->size() : 0;
		_index->resize(numberOfParticlesAtPreviousGeneration);
		int	indexCount = 0;


		// indexation //////////////////////////////////////
		for(int i_index = 0; i_index < isumMax ;i_index++){
			for(unsigned int j_index = 0; j_index < _generations->at(numberOfGeneration-1).getMatrix()->at(i_index).size() ;j_index++){
				//cout<<" in summ !!"<<endl;

				_index->at(indexCount).i = i_index;
				_index->at(indexCount).j = j_index;
				indexCount++;
			}

		}
	}
	double getEnergy(void){return _energyLab;}
	void deleteParticleFromSoftCollision(void){
		if(_hardInteractionCount && _softToHard){


			//change 22.07.14
			//_softInteractionCount--;
			//change 22.07.14

			_softToHard = false;
			_vertexOfBadHardInteraction->resize(0);
			//todo delete particle from corresponding soft collisions
			if(_indexSoftToHard->size())do{
				if(_verbose)cout<<" _indexSoftToHard back=  "<<_indexSoftToHard->back()<<" _indexSoftToHard->size() "<<_indexSoftToHard->size()<<endl;
				_outOfNucleus->erase(_outOfNucleus->begin()+_indexSoftToHard->back());
				_indexSoftToHard->pop_back();
			}while(_indexSoftToHard->size()!= 0);



			if(_indexSoftToHardBadInit->size())do{
				if(_verbose)cout<<" _indexSoftToHard back=  "<<_indexSoftToHardBadInit->back()<<" _indexSoftToHardBadInit->size() "<<_indexSoftToHardBadInit->size()<<endl;
				_notInit->erase(_notInit->begin()+_indexSoftToHardBadInit->back());
				_indexSoftToHardBadInit->pop_back();
			}while(_indexSoftToHardBadInit->size()!= 0);

		}
	}
	std::vector <generation> *_generations;
	std::vector <index> * _index;
	std::vector <hardpingParticle> * _finalState;
	std::vector <hardpingParticle> * _notInit;
	std::vector <hardpingParticle> * _cutMassive;
	std::vector <hardpingParticle> * _outOfNucleus;
	std::vector <unsigned int>* _indexBadInitializations;
	std::vector <unsigned int>* _indexSoftToHard;
	std::vector <unsigned int>* _indexSoftToHardBadInit;
	std::vector <unsigned int>* _vertexOfInteraction;
	std::vector <unsigned int>* _vertexOfBadHardInteraction;
	double _impactParameterMax;
	double _impactParameterMin;
	hardpingParticle _initialParticle;
	void addGeneration(/*	Pyt hia8::generation * gen*/){
		//generation gener(1,1,1);
		if(_verbose)cout<<"i am in addGeneration"<<endl;
		if(_generation==NULL){
			cerr<<"pointer '_generation' is not initialized. nothing happened "<<endl;
			return;
		}
		_generations->push_back(*_generation);
		return;
	}
	int energyLoss(hardpingParticle* particleA, double zCoordinateOfCollisions){
		// procedure return 0 if particle is adsorbed and 1 if not
//		zCoordinateOfCollisions - this is a coordinate of particle up to which particle loses energy
		bool isNucleon = false;
		if(particleA->idAbs() == 2212 || particleA->idAbs() == 2112)isNucleon = true;
		unsigned int i_init = particleA->getIndexNumber();
		double	deltaPath = zCoordinateOfCollisions - particleA->vProd().pz();
		if(particleA->vProd().pz() == -_maxZCoordinate) return 1;
//if interaction occurred beyond boundary of nucleus hardron does not lose energy
		pathInNucleiOutput<<deltaPath<<endl;
		if(deltaPath > 0 && (particleA->isHadron() || isNucleon)){
		//pathInNucleiOutput<<deltaPath<<endl;
		double	deltaE = deltaPath*_kEnergyLoss;
			if(_verbose)cout<<"zCoordinateOfCollisions = "<<zCoordinateOfCollisions<<" zCoordinateOfCollisionsTemp = "<<zCoordinateOfCollisions<<" vecCoordinate.pz() = "<<particleA->vProd().pz()<<" deltaPath  = "<<deltaPath<<endl;
		//
			if( particleA->e() < deltaE || (particleA->e() - deltaE)*(particleA->e() - deltaE) < particleA->m2() ){ // second condition for correct calculation particle momentum
				_indexBadInitializations->push_back(i_init);
				if(_verbose)cout << "particle "<<particleA->getHistory()->back()<<" is adsorbed "<<endl;
				//continue;
				return 0;
				//todo think may be do not escape from cycle
			}
	    double	deltaP = sqrt(particleA->e()*particleA->e() - particleA->m2()) - sqrt((particleA->e() - deltaE)*(particleA->e() - deltaE)- particleA->m2());
			if(_verbose)cout<<"deltaE = "<<deltaE<<" deltaP = "<<deltaP<<endl;
			particleA->e(particleA->e() - deltaE);
			particleA->pz(particleA->pz() - deltaP);

		}else{

			if(_verbose)	cout<<"energy lose do not occurred, deltapath = "<<deltaPath<<" id particle = "<<particleA->id()<<endl;
			deltaPath = 0;
		}// end of deltaPath > 0 && particleA->isHadron()
		return 1;
	}// end energyLoss
	int getIdTargetNucleon(void){
		double probability = double(_targetNucleus.Z() - _scatteringOnProton)/(_targetNucleus.A() - _scatteringOnProton - _scatteringOnNeutron);
		if(_verbose)cout<< "probability = "<<probability<<endl;
		if(getRandom() < double(_targetNucleus.Z() - _scatteringOnProton)/(_targetNucleus.A() - _scatteringOnProton - _scatteringOnNeutron)){
			if(_verbose)cout<< " in proton collision "<<endl;
			//particleB->id(2212);
			_scatteringOnProton++;
			return 2212;

		}else{
			if(_verbose)cout<< " in neutron collision "<<endl;
			//particleB->id(2112);
			_scatteringOnNeutron++;
			return 2112;
		}
	}
	double getEnergyCut(void){return _energyCut;}
	void setVaribles(void);
	bool isBad(unsigned int index){
		for(unsigned int ibad = 0; ibad < _vertexOfBadHardInteraction->size(); ibad++){
			if(index == _vertexOfBadHardInteraction->at(ibad))return true;
		}
		return false;
	}

	void saveParticle4VectorsAfterEnergyLoss(hardpingParticle* particleA){
		char ch;
		//cout<<" in saveParticle4VectorsAfterEnergyLoss beagin"<<endl;
		//cin>>ch;
	//	cout<<"saveParticle4VectorsAfterEnergyLoss1 = "<<particleA->p();

		particleA->rotateHardping();

	//	cout<<"saveParticle4VectorsAfterEnergyLoss2 = "<<particleA->p();
	//	cout.precision(12);
	//	cout<<"1particleA->theta() "<<particleA->theta()<<" 1particleA->phi() "<<particleA->phi()<<endl;
		particleA->setAngles();

		int numberOfGeneration = 0;
		numberOfGeneration = particleA->getNumberOfCurrentGeneration();
		if(numberOfGeneration)	{
		//	cout<<" particle momentum before energy loss "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).p();
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).p(particleA->p());
//		/	cout<<"saveParticle4VectorsAfterEnergyLoss3 "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).p()<<endl;
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).vProd(particleA->vProd());
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).id(particleA->id());
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).setPhiHardping(particleA->getPhiHardping());
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).setThetaHardping(particleA->getThetaHardping());
			if(_verbose)cout<<" if sca part "<<particleA->getIdscatteringParticle()<<endl;
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).scatteringOnParticle(particleA->getIdscatteringParticle());
		//	cout<<" particle momentum after  energy loss "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).p();
		}else{
			//todo extrapolate to the nucleus
		//	_generations->at(0).getMatrix()->at(0).at(0).scatteringOnParticle(particleA->getIdscatteringParticle());
		}
		if(numberOfGeneration&&0){
			cout.precision(12);
			//particleA->rotateBackHardping();
			cout<<"px "<<particleA->px()<<endl;
			cout<<"py "<<particleA->py()<<endl;
			cout<<"pz "<<particleA->pz()<<endl;
			cout<<"saveParticle4VectorsAfterEnergyLoss px "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).px()<<endl;
			cout<<"saveParticle4VectorsAfterEnergyLoss py "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).py()<<endl;
			cout<<"saveParticle4VectorsAfterEnergyLoss pz "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).pz()<<endl;
			cout<<"saveParticle4VectorsAfterEnergyLoss theta "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).theta()<<endl;
			cout<<"saveParticle4VectorsAfterEnergyLoss phi "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).phi()<<endl;
		}






	//	particleA->rotateBackHardping();
//		cout<<"saveParticle4VectorsAfterEnergyLoss3 = "<<particleA->p();
		//cout<<" in saveParticle4VectorsAfterEnergyLoss end"<<endl;
	//	cin>>ch;
	}
	bool findDrellYanPairs(int i_pyEv, hardpingParticle* particleA);


	bool addVertexOfInteraction(hardpingParticle* particleA){
		char ch;
		if(_generations->at(particleA->getNumberOfCurrentGeneration()).getMatrix()->at(particleA->getIndexNumber()).size() == 0){
			//this for situation when energy of all produced particles below energyCut
			if(_verbose)cout<<"getMatrix()->at(i_init).size() == 0"<<endl;
			_indexBadInitializations->push_back(particleA->getIndexNumber());
			return false;

		}else{
			// add vertex of interaction
			/*
			cout<<" history back particleA = "<<particleA->getHistory()->back()<<endl;


			cout<<"history particleA : ";
			for(unsigned int i = 0 ; i< particleA->getHistory()->size(); i++){
				cout<<particleA->getHistory()->at(i)<<" ";
			//	cin>>ch;
			}
			cout<<endl;
//			cin>>ch;

			 */
			_vertexOfInteraction->push_back(particleA->getHistory()->back());
			return true;
		}
	}


bool prepareNewGeneration(hardpingParticle* particleA,int i_pyEv);

void deleteBadInitializationParticlesFromNewGeneration(int numberOfGeneration){
	char ch;
	//cout<<" in deleteBadInitializationParticlesFromNewGeneration beagin"<<endl;
	//cin>>ch;
	if(_indexBadInitializations->size()){
	//	for(int i=0;i<_indexBadInitializations->size();i++){
	//
	//			cout<<"i = "<<i<<"  _indexBadInitializations = "<<_indexBadInitializations->at(i)<<endl;
	//	}

		do{
				//cout<<"  _indexBadInitializations = "<<_indexBadInitializations->back()<<endl;
	//		cout<<"bad size "<<_generations->at(numberOfGeneration).getMatrix()->at(_indexBadInitializations->back()).size()<<endl;//
	//		cout<<" val = "<<_generations->at(numberOfGeneration).getMatrix()->at(_indexBadInitializations->back()).at(0).p()<<endl;
			if(_verbose)cout<<" _indexBadInitializations =  "<<_indexBadInitializations->back()<<endl;
			_generations->at(numberOfGeneration).getMatrix()->erase(_generations->at(numberOfGeneration).getMatrix()->begin()+_indexBadInitializations->back());
			_indexBadInitializations->pop_back();
		}while(_indexBadInitializations->size()!= 0);
	}

	for(unsigned int i_out = 0; i_out <_generations->at(numberOfGeneration).getMatrix()->size(); i_out++){

		for(unsigned int j_out = 0; j_out < _generations->at(numberOfGeneration).getMatrix()->at(i_out).size();j_out++){

			if(_verbose)cout<<" i = "<<i_out<<" j = "<<j_out<<" val = "<<_generations->at(numberOfGeneration).getMatrix()->at(i_out).at(j_out).p()<<" id = "<<_generations->at(numberOfGeneration).getMatrix()->at(i_out).at(j_out).id()<<" index = "<< _generations->at(numberOfGeneration).getMatrix()->at(i_out).at(j_out).getHistory()->back()<<endl;
		}
	}
	//cout<<" in deleteBadInitializationParticlesFromNewGeneration end"<<endl;
	//cin>>ch;
}
unsigned int getIndexHardOfHardCollision(void){

	unsigned int randomPossitionOfHardScattering = 0;
	char ch;
	//cout<<"_vertexOfInteraction size = "<<_vertexOfInteraction->size()<<endl;
	//cout<<"back = "<<_vertexOfInteraction->back()<<endl;

	//_softCollisionIterator = _vertexOfInteraction->begin();

	_indexOfHardCollosion = _vertexOfInteraction->back();

//	cout<<"_indexOfHardCollosion = "<<_indexOfHardCollosion<<endl;
//	cin>>ch;
	// we choose last of soft collisions. if hardping cant't initialize we get random position from _vertexOfInteraction

	if(isBad(_indexOfHardCollosion)){
		if(_vertexOfInteraction->size() == 1 || _vertexOfInteraction->size() == _vertexOfBadHardInteraction->size() +1){
				if(_verbose)cout<< " _vertexOfInteraction->size() == 1 "<<endl;
				_indexOfHardCollosion = 0;
				/// todo add condition for initial particle
			//	particleA->p(_initialParticle.p());
			//	particleA->vProd(_initialParticle.vProd());
				//todo may be add some new propertes

			}else{
				do{

					_softCollisionIterator = _vertexOfInteraction->begin();
					randomPossitionOfHardScattering = (unsigned int)(getRandom()*(_vertexOfInteraction->size() - 1)) +1;
					_softCollisionIterator = _softCollisionIterator + randomPossitionOfHardScattering;
					if(_verbose)cout<< "randomPossitionOfHardScattering = "<< randomPossitionOfHardScattering<<" index of hard scattering = "<< *_softCollisionIterator<<endl;
					_indexOfHardCollosion = *_softCollisionIterator;
			//		cout<<"_indexOfHardCollosion = "<<_indexOfHardCollosion<<endl;
			//		cin>>ch;
					if(_verbose)cout<< " is bad = "<< isBad(_indexOfHardCollosion)<<endl;
				}while(isBad(_indexOfHardCollosion));
			}//end of else _vertexOfInteraction->size() == 1
	}

	return _indexOfHardCollosion;
}
class pythia6Event{
public:

	pythia6Event(){
		pythia6Px = new vector<double>;
		pythia6Py = new vector<double>;
		pythia6Pz = new vector<double>;
		pythia6E = new vector<double>;
		pythia6Id = new vector<int>;
		leptonPosition = new vector<int>;
		pythia6Particle = new vector<hardpingParticle>;
	}
	~pythia6Event(){
		delete pythia6Px;
		delete pythia6Py;
		delete pythia6Pz;
		delete pythia6E;
		delete pythia6Id;
		delete pythia6Particle;
	}
	std::vector<double>* pythia6Px;
	std::vector<double>* pythia6Py;
	std::vector<double>* pythia6Pz;
	std::vector<double>* pythia6E;
	std::vector<  int >* pythia6Id;
	std::vector<  int >* leptonPosition;
	std::vector<hardpingParticle>* pythia6Particle;

};
double getQ2LeptonHadron(pythia6Event* py6Ev){

		double Q2 = 0;


		// in position 0 and 1 locate initial particles

 		for(int i = 3; i <= py6Ev->pythia6Particle->size(); i++ ){
			if(py6Ev->pythia6Particle->at(0).id() == py6Ev->pythia6Particle->at(i).id())py6Ev->leptonPosition->push_back(i);
		}
		Vec4 deltaLeptonP(0);

		//py6Ev->leptonPosition->at(py6Ev->leptonPosition->)
		deltaLeptonP = py6Ev->pythia6Particle->at(0).p() - py6Ev->pythia6Particle->at(py6Ev->leptonPosition->at(0)).p();

		double nu = 0; //energy of virtual photon
		nu = (deltaLeptonP * py6Ev->pythia6Particle->at(1).p())/py6Ev->pythia6Particle->at(1).m();

		Q2 = deltaLeptonP*deltaLeptonP;

		double bl = 0;

		bl = nu/_kEnergyLoss;

		return Q2;
/*        IF(((ABS(IHNT2(5)).LT.11).OR.(ABS(IHNT2(5)).GT.18))) THEN
            z=2.0D0*sqrt(PATT(i,1)**2.0D0+PATT(i,2)**2.0D0)
&              /sqrt(PARI(33)*PARI(34)*PARI(12))
            DZA(I)=2.0D0*DSQRT(PATT(i,1)**2.0D0+PATT(i,2)**2.0D0)
&                   /DSQRT(DXFREL1*DXFR2*HINT1(1)*HINT1(1))
            z=DZA(I)
            IF(z.GT.1.0D0) THEN
                z=0.99D0
            ENDIF
            DZA(I)=z
            snu=dphadron/z
        ENDIF
        */
}
double getQ2HadronHadron(pythia6Event* py6Ev){

}
bool findHardParticle(hardpingParticle* particleA){
	bool itIsOk = false;
	bool exit = false;
	char ch;
	if(_verbose)cout<<"_generations->size() = "<<_generations->size()<<endl;
	for(unsigned int indexOfGeneration = 0; indexOfGeneration < _generations->size(); indexOfGeneration++ ){
		if (indexOfGeneration == 3){
			//cout<<"_generations->at(indexOfGeneration).getMatrix()->size() = "<<_generations->at(indexOfGeneration).getMatrix()->size()<<endl;
	//		cin>>ch;
		}
		//cout<<"indexOfGeneration = "<<indexOfGeneration<<" _generations->size() = "<<_generations->size()<<endl;
		if(exit)break;
		for(unsigned int initialParticleIndex = 0; initialParticleIndex < _generations->at(indexOfGeneration).getMatrix()->size(); initialParticleIndex++){
			//cout<<"at indexOfGeneration = "<<indexOfGeneration<< "_generations->at(indexOfGeneration).getMatrix()->size() = "<<_generations->at(indexOfGeneration).getMatrix()->size()<<endl;
			if(exit)break;
			for(unsigned int producedParticleIndex = 0; producedParticleIndex < _generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).size(); producedParticleIndex++){
				if(exit)break;
				//cout<<"indexOfGeneration = "<<(int)indexOfGeneration<<" initialParticleIndex = "<<initialParticleIndex<<" producedParticleIndex = "<<producedParticleIndex<<endl;
				//cout<<"_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->size() = "<<_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->size()<<endl;
				//todo some particle havent history
				for(int i =0; i< _generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->size(); i++){
				//	cout<<_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->at(i)<<endl;
				}
				if(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->back() ==
						_indexOfHardCollosion){
				//	cin>>ch;

					if(_verbose)	cout<<" diff = "<<_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).p();
					particleA->p(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).p());
					particleA->vProd(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).vProd());
					particleA->scatteringOnParticle(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getIdscatteringParticle());
					particleA->setPhiHardping(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).phi());
					particleA->setThetaHardping(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).theta());
					particleA->getHistory()->clear();
					particleA->getHistory()->assign(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->begin(),_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->end());
					particleA->id(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).id());
/*
					for(unsigned int i = 0 ; i < particleA->getHistory()->size(); i++){
						cout<<particleA->getHistory()->at(i)<<" ";
					}
					cout<<endl;
*/
					//todo may be add some new parameter
					if(_verbose)cout<<" particleA = "<<particleA->p()<<" idb "<<particleA->getIdscatteringParticle()<<endl;
					return itIsOk;
					//exit = true;
				//	cin>>ch;
					break;
				}// end of if getHistory()->back() ==_indexOfHardCollosion
			//	cout<<"1"<<endl;
			}// end producedParticleIndex cycle
		//	cout<<"2"<<endl;
		}// end initialParticleIndex cycle
	//	cout<<"3"<<endl;
	}// end indexOfGeneration cycle
//	cout<<"4"<<endl;
	return itIsOk;
}

bool ifNoHardCollisionHappened(unsigned int numberOfGeneration/*, hardpingParticle * particleA*/){
	char ch;
	//cout<<" in ifNoHardCollisionHappened beagin"<<endl;
	//cin>>ch;
	hardpingParticle* particleA;
	particleA = new hardpingParticle();

	bool softToHard = false;
	//unsigned int randomPossitionOfHardScattering = 0 ;

//	cout<<"mat size _generations->at(numberOfGeneration-1).getMatrix()->size() = "<<_generations->at(numberOfGeneration-1).getMatrix()->size()<<" _hardInteractionCount = "<<_hardInteractionCount<<endl;
//	cin>>ch;
	if(_generations->at(numberOfGeneration-1).getMatrix()->size() == 0 && _hardInteractionCount == 0){


		if(_verbose)cout<<"_vertexOfInteraction->size() = "<<_vertexOfInteraction->size()<<" _vertexOfBadHardInteraction->size() = "<<_vertexOfBadHardInteraction->size()<<endl;
			if( _vertexOfInteraction->size() == _vertexOfBadHardInteraction->size())return softToHard;//condition for exit

				_softToHard = true;

				if(_verbose)cout<< "!!!!!!!!!!!!!!!!!!hard to soft part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
				if(_verbose)for(unsigned int iv = 0; iv < _vertexOfInteraction->size(); iv ++){
					cout<<_vertexOfInteraction->at(iv)<<" ";
				}
				if(_verbose)cout<<endl;
				if(_verbose)cout<<" vertexOfInteraction->size() "<<_vertexOfInteraction->size()<<" vertexOfBadHardInteraction->size() "<<_vertexOfBadHardInteraction->size()<<endl;


				if(getIndexHardOfHardCollision()){
					//when we get particle she is oriented along initial z direction
					findHardParticle(particleA);


				}else{
						particleA->p(_initialParticle.p());
						particleA->vProd(_initialParticle.vProd());
						particleA->id(_initialParticle.id());
						particleA->scatteringOnParticle(_initialParticle.getIdscatteringParticle());
						particleA->getHistory()->clear();
						for(unsigned int i_history = 0; i_history < _initialParticle.getHistory()->size(); i_history++){
							particleA->getHistory()->push_back(_initialParticle.getHistory()->at(i_history));
						}
				}
				findParticlesWithIndexOfHardCollisions();

				_generations->at(numberOfGeneration-1).getMatrix()->resize(1);

				_generations->at(numberOfGeneration-1).getMatrix()->at(0).push_back(*particleA);//at(0).push_back(*particleA);
				if(_verbose)cout<<" new generation = "<< _generations->at(numberOfGeneration-1).getMatrix()->at(0).at(0).p();
			}// end of if(_generations->at(numberOfGeneration-1).getMatrix()->size() == 0 && _hardInteractionCount == 0)
	delete particleA;
	return _softToHard;
	//cout<<" in ifNoHardCollisionHappened end"<<endl;
	//cin>>ch;
}
void findParticlesWithIndexOfHardCollisions(void){
	// find all particle which produced in soft collision with index _indexOfHardCollosion and put them into massive _indexSoftToHard
				for(unsigned int isoft = 0; isoft < _outOfNucleus->size(); isoft++){
					for(unsigned int jsoft = 0; jsoft <_outOfNucleus->at(isoft).getHistory()->size() ; jsoft++){
						if(_outOfNucleus->at(isoft).getHistory()->at(jsoft) == _indexOfHardCollosion) _indexSoftToHard->push_back(isoft);
					}
				}

				for(unsigned int isoftBadInit = 0; isoftBadInit < _notInit->size(); isoftBadInit++){
					for(unsigned int jsoftBadInit = 0; jsoftBadInit <_notInit->at(isoftBadInit).getHistory()->size() ; jsoftBadInit++){
						if(_notInit->at(isoftBadInit).getHistory()->at(jsoftBadInit) == _indexOfHardCollosion) _indexSoftToHardBadInit->push_back(isoftBadInit);
					}
				}
}
void getParticleFromPreviousGeneration(hardpingParticle* particleA){
	// get particle from previous generation
	unsigned int numberOfGeneration = 0;
	unsigned int i_init = 0;
	numberOfGeneration = particleA->getNumberOfCurrentGeneration();
	i_init = particleA->getIndexNumber();


	//particleA->getHistory()->clear();
	cout.precision(12);
	/*cout<<"hui "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).px()<<endl;
	cout<<"hui "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).py()<<endl;
	cout<<"hui "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).pz()<<endl;





cout<<"px "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).px()<<endl;
cout<<"py "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).py()<<endl;
cout<<"pz "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).pz()<<endl;
cout<<"theta "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).theta()<<endl;
cout<<"phi "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).phi()<<endl;
*/
	*particleA = _generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j);
	// set index number of particle in current generation
	particleA->setInexNumber(i_init);
	//and number of current generation
	particleA->setNumberOfCurrentGeneration(numberOfGeneration);


//	cout<<"theta "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).theta()<<" phi "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).phi()<<endl;
//	cout << "particleA->theta() 0 "<<particleA->theta()<<" particleA->phi() "<<particleA->phi()<<endl;
//	cout<<"p  getParticleFromPreviousGeneration= "<<particleA->p();
//	cout.precision(12);
//	cout << "particleA->theta() 0 "<<particleA->theta()<<" particleA->phi() "<<particleA->phi()<<endl;
	//particleA->rotateBackHardping();
	//cout<<"p  getParticleFromPreviousGeneration 2= "<<particleA->p();
	////////////////////////////////////////
//	cout<<"particleA id = "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).id()<<endl;
//	cout<<"particleA id = "<<particleA->id()<<endl;
	if(	!_softToHard){
		// todo may be we need remember angles after particle production immediately


		/// save angles relatively initial beam direction///
//		particleA->setAngles();
		//////////////////////////////////////////////////
	//	cout.precision(10);

		if(_verbose)cout<<" momentum  before  rotate   = "<<particleA->p()<<endl;
		if(_verbose)cout<<" coordinate before rotation = "<<particleA->vProd()<<endl;

		//dispose momentum of particle along z direction
	//	particleA->rotateBackHardping();
		//////////////////////////////////////////////////

		if(_verbose)cout<<" momentum  after rotate_    = "<<particleA->p()<<endl;
		if(_verbose)cout<<" coordinate after rotation_ = "<<particleA->vProd()<<endl;



		//////////////////////////////////////////////////////////
	}

}
void setInitinalImpactAndIndex(hardpingParticle* particleA){

	char ch;
	////////// calculating impact parameter of incident particles//////////////////////////////////////////
	double xImpact = 0, yImpact = 0;
	getImpactParameter(xImpact,yImpact);
	Vec4 vecCoordinate(0);
	//////////end of calculating impact parameter of incident particles//////////////////////////////////////////

	/////// set coordinate of incident particle////////////////////////////////////

	vecCoordinate.px(xImpact);
	vecCoordinate.py(yImpact);
	vecCoordinate.pz(-_maxZCoordinate);
	particleA->vProd(vecCoordinate);

	//////// end of set coordinate of incident particle ////////////////////

	// set sequentially number for each incident particle. It done for usable trace of particle in program

	particleA->getHistory()->clear();
	particleA->getHistory()->push_back(_indexParticle);
//	cout<<" initial index particle "<<particleA->getHistory()->back()<<endl;
	//cin>>ch;
	_indexParticle++;

	//////////////////////////////////////////////////////////////////////////////////////////////////////
}
void softToHard(hardpingParticle* particleA,hardpingParticle* particleB){
	char ch;
	//cout<<" in softToHard beagin"<<endl;
	//cin>>ch;
	particleA->setHard();
	particleA->setSoft(false);
	//coordinateHardOutput<<particleA->vProd();
	if(particleA->getIdscatteringParticle() == 0){
		_scatteringOnNeutron = 0;
		_scatteringOnProton = 0;
		particleB->id(getIdTargetNucleon());
	}//end if(particleA->getIdscatteringParticle() == 0)
	else{
		if(_verbose)cout<< " else of (particleA->getIdscatteringParticle() == 0) "<<endl;
		particleB->id(particleA->getIdscatteringParticle());
		if(_verbose)cout<<"id particle B = "<<particleB->id()<<endl;
	} // else if(particleA->getIdscatteringParticle() == 0)
	//particleB->id(particleA->getIdscatteringParticle());
	//cout<<"here i am "<<particleB->id()<<endl;
//	cout<<" in softToHard end"<<endl;
	//cin>>ch;
}
bool checkEnergyCut(hardpingParticle* particleA){
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// if after energy loss particle not have enough energy for initializing new collision.
	// if particle set out of nucleus put her into energy cut massive,
	// if not decrease her energy correspond to residual path in nucleus.
	// if energy became bellow zero or equal assume that particle is adsorbed if not add her into energy cut massive.

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//part for particles which can't initialize new event. for this particles - if z coordinate position more that nuclear radius assumed that particle set out of nucleus
	//	else we propose particle propagates in nuclear matter along the z direction and keeps losing energy until it reaches the boundary of nucleus (i.e. nuclear radius).
	bool isNotAdsorbed = true;
	if(particleA->e() < _energyCut){

		isNotAdsorbed = energyLoss(particleA,_targetNucleus.getNuclearRadius());
		// isNotAdsorbed = 0 - particle is adsorbed
		// isNotAdsorbed = 1 - all right
		if(isNotAdsorbed){
			particleA->setOut();
			//dispose momentum of particle along initial beam direction
			particleA->rotateHardping();
			// put particle in massive
			_outOfNucleus->push_back(*particleA);
			// massive with index of particles which not initialize next wave
			_indexBadInitializations->push_back(particleA->getIndexNumber());
		}

	}//end of if (particleA->e() < _energyCut)
	return isNotAdsorbed;
}
void finalOutput(void){
	if(_verbose)for(unsigned int i=0;i < _finalState->size() ; i++){
		//_finalState->at(i).setHardpingVector(*particleA);
		cout<<"i = "<<i<<"  id = "<<_finalState->at(i).id()<<" energy = "<<_finalState->at(i).e()<<endl;
		cout<< "index =  ";
		for(unsigned int i = 0; i < _finalState->at(i).getHistory()->size(); i++){
			cout<<_finalState->at(i).getHistory()->at(i)<<" ";
		}
	}


	if(_verbose)cout<<"final state size = "<<_finalState->size()<<endl;

	for(unsigned int i = 0 ; i < _cutMassive->size() ; i++ ){
			//_finalState->at(i).setHardpingVector(*particleA);
		if(_verbose)cout<<"i = "<<i<<"  id = "<<_cutMassive->at(i).id()<<" energy = "<<_cutMassive->at(i).e()<<endl;
		if(_verbose)cout<< "index =  ";
			for(unsigned int j = 0; j < _cutMassive->at(i).getHistory()->size(); j++){
					cout<<_cutMassive->at(i).getHistory()->at(j)<<" ";
			}
			if(_verbose)cout<<endl;
		}
	if(_verbose)cout<<"_cutMassive state size = "<<_cutMassive->size()<<endl;
	for(unsigned int ii = 0; ii<	_indexSoftToHard->size(); ii++){
		if(_verbose)cout<< " index of soft to hard  = "<<_indexSoftToHard->at(ii)<<" _indexSoftToHard->size() = "<<_indexSoftToHard->size()<<endl;
	}

	if(_verbose)for(unsigned int i = 0 ; i < _notInit->size() ; i++ ){
		//_finalState->at(i).setHardpingVector(*particleA);
		cout<<"i = "<<i<<"  id = "<<_notInit->at(i).id()<<" energy = "<<_notInit->at(i).e()<<endl;
		cout<< "index =  ";
		for(unsigned int j = 0; j < _notInit->at(i).getHistory()->size(); j++){
				cout<<_notInit->at(i).getHistory()->at(j)<<" ";
		}
		cout<<endl;

	}
	if(_verbose)cout<<"_notInit state size = "<<_notInit->size()<<endl;

	if(_verbose)for(unsigned int i = 0 ; i < _outOfNucleus->size() ; i++ ){
		//_finalState->at(i).setHardpingVector(*particleA);
		cout<<"i = "<<i<<"  id = "<<_outOfNucleus->at(i).id()<<" energy = "<<_outOfNucleus->at(i).e()<<endl;
		cout<< "index =  ";
		if(_verbose)for(unsigned int j = 0; j < _outOfNucleus->at(i).getHistory()->size(); j++){
				cout<<_outOfNucleus->at(i).getHistory()->at(j)<<" ";
		}
		if(_verbose)cout<<endl;
	}
	if(_verbose)cout<<"_outOfNucleus state size = "<<_outOfNucleus->size()<<endl;
	if(_verbose)cout<<"_vertexOfInteraction = "<<endl;
	if(_verbose)for(unsigned int iv = 0; iv < _vertexOfInteraction->size(); iv ++){
			cout<<_vertexOfInteraction->at(iv)<<" ";
	}
	if(_verbose)cout<< "number of hard scattering = "<<_hardInteractionCount<<endl;
	if(_verbose)cout<< "number of soft scattering = "<<_softInteractionCount<<endl;
}
	void notPythiaNext(hardpingParticle* particleA){
		if(_verbose)cout<<"in !pythiaNextFlag = "<<particleA->id()<<"  "<<particleA->p()<<endl;

		int isNotAdsorbed = 0;
		if(_softToHard){
			_indexBadInitializations->push_back(particleA->getIndexNumber());
			_vertexOfBadHardInteraction->push_back(_indexOfHardCollosion);
			_indexSoftToHard->resize(0);
			return;
		}
		if(particleA->getNumberOfCurrentGeneration() && particleA->isHadron()){
			//////////////////////////////////////////////////////////////////////////////
			//if pythia can not initialize new event we assumed that projectile hardron //
			//lose energy during propagation in nucleus along the z direction			//
			//until it reaches the boundary of nucleus (i.e. nuclear radius).			//
			//////////////////////////////////////////////////////////////////////////////

			if(_verbose)cout<<"1before particleA->isHadron() = "<<particleA->p();
			//set particle momentum along z direction
			//particleA->setAngles();
			//particleA->rotateBackHardping();
			if(_verbose)cout<<"1after particleA->isHadron() = "<<particleA->p();

		//	cin>>lll;

			isNotAdsorbed = energyLoss(particleA,_targetNucleus.getNuclearRadius());
			// isNotAdsorbed = 0 - particle is adsorbed
			// isNotAdsorbed = 1 - all right
			if(isNotAdsorbed){
				particleA->setOut();
				//dispose momentum of particle along initial beam direction
				particleA->rotateHardping();
			}else{return;}
		}// end of if(numberOfGeneration && particleA->isHadron())
		// put particle in massive
		_notInit->push_back(*particleA);
		// massive with index of particles which not initialize next wave
		_indexBadInitializations->push_back(particleA->getIndexNumber());
		return;

	}
	void useFortranMethod(bool value){
		_fortranHardping = value;
	}


	//double ::getRandomFromFile();
	void getImpactParameter(double &impactX,double &impactY){
	//	double impactParameterMax = min(3*_targetNucleus.getNuclearRadius(),_maxZCoordinate + 6.06); // from Fortran version. 6.06 fm - is proton characteristic

		//BB=DSQRT(BMIN**2+RAN(NSEED)*(BMAX**2-BMIN**2))



		double impactPar = sqrt(_impactParameterMin*_impactParameterMin + getRandom()*(_impactParameterMax*_impactParameterMax - _impactParameterMin*_impactParameterMin));
		//impactPar = getRandom()*impactParameterMax;


		double phiImpact = 2*M_PIl*getRandom();
		//cout<<" phi1 = "<<phiImpact<<endl;
		impactX = sin(phiImpact)*impactPar;
		impactY = cos(phiImpact)*impactPar;

	}


	/*
	double getRandomFromFile(){
		double a;
		_randomFile>>a;
		return a;
	}

	*/

	int getNumberOfSoftInteraction(void){
		return _softInteractionCount;
	}


	bool   _firstCall;
	unsigned int _softInteractionCount;
private:

	double woodSaxonSkinDepth() const { return 0.53;  }  ///< returns surface (0.53 fm for Au)

//	double getNuclearRadius  () const { return _r0 ; }

//	Pythia8::Pythia* _pythia;


	double _BQ;
	double _BP;
	double _BN;
	int    _Z;                      ///< atomic number of nucleus
	int    _A;                      ///< nucleon number of nucleus
	double _energyLab;
	int _randomSeed;
	Rndm * _random;
	double _xPartonCoord;
	double _yPartonCoord;
	double _zPartonCoord;
	double _xPartonMomentum;
	double _yPartonMomentum;
	double _zPartonMomentum;
	double _absPartonMomentum;
	Vec4 vek;
//	std::vector *a;

//	std::vector< std::vector< Vec4 > > * m(1, std::vector<Vec4>(1));
	double _energyCut;
	double _kEnergyLoss;
	double _r0;				// nuclear radius
	double _rho0;
	double _Q0;
	double _lastResultOfThicknessCalculation;
	double _lastErrorOfThicknessCalculation;
	unsigned int _scatteringOnProton;
	unsigned int _scatteringOnNeutron;
	//double _precisionOfNuclearThicknessCalculation;
	//bool   _firstCall;
	double _maxZCoordinate; //Maximum radial coordinator for target nucleons
	unsigned int _factorials[13];//={1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800};
	nucleus _targetNucleus;
	nucleus _projectileNucleus;
	bool 	_isScattering;
	hardpingParticle* _tempParticle;
	unsigned int _indexParticle;
	bool _softInteractionFlag;
	bool _hardInteractionFlag;
	bool _hardInteraction; // correspond to situation when hard collision do not occurred in regular regime of running
	unsigned int _hardInteractionCount;
	//unsigned int _softInteractionCount;
	unsigned int _indexOfHardCollosion;
	int _verbose;
	bool _softToHard;
	bool _fortranHardping;
	Timer _time;
	bool _cutMass;
	ifstream _randomFile;

	std::vector<unsigned int>::iterator _softCollisionIterator;
}; // end nucleus class

}// end of namespace pythia8
#endif  // NUCLEUS_H
