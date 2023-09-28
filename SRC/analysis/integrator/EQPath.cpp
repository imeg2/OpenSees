/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// Written: M. Salehi opensees.net@gmail.com
// Created: 02/19
// Revision: A

//refs
//Structural Engineeringand Mechanics   Volume 48, Number 6, December25 2013, pages 849 - 878
//DOI: https://doi.org/10.12989/sem.2013.48.6.849	
//
//Comprehensive evaluation of structural geometrical nonlinear solution techniques Part I : Formulation and characteristics of the methods
//M.Rezaiee - Pajand, M.Ghalishooyan and M.Salehi - Ahmadabad
//FULLTEXT : https://www.researchgate.net/publication/264146397_Comprehensive_evaluation_of_structural_geometrical_nonlinear_solution_techniques_Part_I_Formulation_and_characteristics_of_the_methods


//Structural Engineeringand Mechanics   Volume 48, Number 6, December25 2013, pages 879 - 914
//DOI: https://doi.org/10.12989/sem.2013.48.6.879	
//
//Comprehensive evaluation of structural geometrical nonlinear solution techniques Part II : Comparing efficiencies of the methods
//M.Rezaiee - Pajand, M.Ghalishooyan and M.Salehi - Ahmadabad
//FULLTEXT : https://www.researchgate.net/publication/263361974_Comprehensive_evaluation_of_structural_geometrical_nonlinear_solution_techniques_Part_II_Comparing_efficiencies_of_the_methods

#include <EQPath.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <math.h>
#include "ElementIter.h"
#include "NodeIter.h"
#include "Domain.h"
#include "Node.h"
#include "Element.h"
#include "DOF_Group.h"

EQPath::EQPath(double arcLen, int method, ID* ids, ID* dofIDs, int gmrdType)
	:StaticIntegrator(INTEGRATOR_TAGS_EQPath),
	arclen(arcLen),
	uq(0), ur(0), du(0), du0(0), uq_n1(0), uq_n0(0), q(0), type(method), dl(0),
	nitr(0), sign_n0(1), theIDs(ids), theDofIDs(dofIDs), step(0), iter(0)
{

	GRMDType = gmrdType;
	theTaggedObjectDic = {};
}

EQPath::~EQPath()
{
	// delete any vector object created
	if (uq != 0)
		delete uq;
	if (uq_n0 != 0)
		delete uq_n0;
	if (uq_n1 != 0)
		delete uq_n1;
	if (ur != 0)
		delete ur;
	if (du != 0)
		delete du;
	if (du0 != 0)
		delete du0;
	/*if (theDofs != 0)
		delete theDofs;
	if (theNodeIDs != 0)
		delete theNodeIDs;
	if (theElementIDs != 0)
		delete theElementIDs;*/
	if (du0 != 0)
		delete du0;
	if (q != 0)
		delete q;
}

int
EQPath::newStep(void)
{
    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING EQPath::newStep() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    // get the current load factor
    double currentLambda = theModel->getCurrentDomainTime();


    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*q);
    if (theLinSOE->solve() < 0) {
      opserr << "EQPath::newStep(void) - failed in solver\n";
      return -1;
    }

	int ndofs = q->Size();


	if (step == 0) {
		// initialize uq vectors
		uq_n0 = new Vector(ndofs);
		uq_n1 = new Vector(ndofs);

		(*uq_n0) = theLinSOE->getX();
		(*uq_n1) = theLinSOE->getX();

		sign_n0 = 1;
	}
	else
	{
		(*uq_n1) = (*uq_n0);
		(*uq_n0) = theLinSOE->getX();

		//sign_n0 = (*du) ^ (*uq_n0);
		sign_n0 = ((*uq_n0) ^ (*uq_n1)) > 0 ? sign_n0 : -sign_n0;
	}

	du->Zero();

	double dLambda = sign_n0 * arclen / uq_n0->Norm();

	(*du) = dLambda * (*uq_n0);

	du0 = new Vector(ndofs);
	(*du0) = (*du);

	currentLambda += dLambda;
	dl += dLambda;

	// update model with delta lambda and delta U
	theModel->incrDisp(*du);
	theModel->applyLoadDomain(currentLambda);
	if (theModel->updateDomain() < 0) {
		opserr << "EQPath::newStep - model failed to update for new dU\n";
		return -1;
	}

	//iter = 0;
	int flag = 1;

	step++;
	nitr = 0;
	return 0;
}

int
EQPath::update(const Vector &dU)
{
	AnalysisModel* theModel = this->getAnalysisModel();
	LinearSOE* theLinSOE = this->getLinearSOE();
	Domain* domain = theModel->getDomainPtr();
	if (theModel == 0 || theLinSOE == 0) {
		opserr << "WARNING EQPath::update() ";
		opserr << "No AnalysisModel or LinearSOE has been set\n";
		return -1;
	}

	(*ur) = dU;

	theLinSOE->setB(*q);
	theLinSOE->solve();
	(*uq) = theLinSOE->getX();

	double dLambda = 0;
	iter++;

	if (type == EQPath_Method_MRD) // minimum residual disp
	{
		double a = (*ur) ^ (*uq);
		double b = (*uq) ^ (*uq);
		if (b == 0) {
			opserr << "EQPath::update() - zero denominator\n";
			return -1;
		}

		//opserr << "MRD  A:" << a << " B:" << b << endln;
		dLambda = -a / b;
	}
	else if (type == EQPath_Method_NP) // normal plain
	{
		double a = (*du0) ^ (*ur);
		double b = (*du0) ^ (*uq);
		if (b == 0) {
			opserr << "EQPath::update() - zero denominator\n";
			return -1;
		}

		dLambda = -a / b;
	}
	else if (type == EQPath_Method_UNP) // update normal plain
	{
		double a = (*du) ^ (*ur);
		double b = (*du) ^ (*uq);
		if (b == 0) {
			opserr << "EQPath::update() - zero denominator\n";
			return -1;
		}

		dLambda = -a / b;
	}
	else if (type == EQPath_Method_CYL) // cylindrical arc-length
	{
		double A = (*uq) ^ (*uq);
		double B = 2 * (((*du) + (*ur)) ^ (*uq));
		double c1 = ((*du) + (*ur)) ^ (*du);
		double c2 = ((*du) + (*ur)) ^ (*ur);
		double C = c1 + c2 - arclen * arclen;
		double delta = B * B - 4 * A * C;

		dLambda = 0;
		if (delta < 0)
		{
			opserr << "EQPath::update() - negetive denominator\n";
			return -1;
		}
		else if (delta == 0)
		{
			dLambda = -B / 2 / A;
		}
		else
		{
			double sl1 = (-B + pow(delta, 0.5)) / 2 / A;
			double sl2 = (-B - pow(delta, 0.5)) / 2 / A;
			double aa1 = (*du) ^ (*ur);
			double aa2 = (*du) ^ (*du);
			double aa3 = (*du) ^ (*uq);
			double costl1 = aa1 + aa2 + sl1 * aa3;
			double costl2 = aa1 + aa2 + sl2 * aa3;
			dLambda = sl1;
			if (costl2 > costl1)
				dLambda = sl2;
		}

	}
	else if (type == EQPath_Method_MNP) // modified normal flow
	{
		/*double a = (*ur) ^ (*uq);
		double b = (*uq) ^ (*uq);
		if (b == 0) {
			opserr << "EQPath::update() - zero denominator\n";
			return -1;
		}

		dLambda = -a / b;*/
		opserr << "WARNING EQPath::update() ";
		opserr << "Removed update method has been set\n";
		return -1;
	}
	else if (type == EQPath_Method_GDC) // GDC
	{
		double a, b;
		a = (*ur) ^ (*uq_n1);
		b = (*uq) ^ (*uq_n1);


		if (b == 0) {
			opserr << "EQPath::update() - zero denominator\n";
			return -1;
		}
		dLambda = -a / b;

	}
	else if (type == EQPath_Method_MGDC) // MGDC
	{
		double a, b;
		a = (*ur) ^ (*uq_n0);
		b = (*uq) ^ (*uq_n0);

		if (b == 0) {
			opserr << "EQPath::update() - zero denominator\n";
			return -1;
		}
		dLambda = -a / b;
	}
	else if (type == EQPath_Method_MUNP) // Modified Update Normal Plane
	{
		double p1 = (*uq) ^ (*uq);
		double p2 = (*du) ^ (*uq);
		double p3 = (*ur) ^ (*uq);
		double p4 = (*ur) ^ (*du);
		double p5 = (*ur) ^ (*ur);
		double dl0 = -p3 / p1; //Chan constraint


		double A = p1;//(*uq)^(*uq);
		double B = p2 + 2 * p3;//(((*du) + (*ur)+ (*ur)) ^ (*uq));
		double C = p4 + p5;//((*du) + (*ur)) ^ (*ur);
		double delta = B * B - 4 * A * C;
		//opserr << "negative denominator - alpha = " << 0 <<"\n";
		dLambda = 0;
		if (delta < 0)
		{
			Vector* v1 = new Vector(ur->Size());
			Vector* v2 = new Vector(ur->Size());
			(*v2) = (*ur);
			v2->addVector(1, *uq, dl0);
			(*v1) = (*du);
			v1->addVector(1, *v2, 1);
			double l1 = v1->Norm();
			double l2 = v2->Norm();
			double alpha = (C - B * B / 4 / A) / l1 / l2;
			alpha += 0.1 * (1 - alpha);
			//opserr << "negative denominator - alpha = " << alpha <<"\n";
			delta = B * B - 4 * A * (C - alpha * l1 * l2);
			//dLambda=dl0;
		}

		if (delta == 0)
		{
			dLambda = -B / 2 / A;
		}
		else
		{
			double sl1 = (-B + pow(delta, 0.5)) / 2 / A;
			double sl2 = (-B - pow(delta, 0.5)) / 2 / A;

			double* roots = new double[2] { sl1, sl2 };
			int index = chooseBestDLambda(roots, 2, du, uq, ur);
			if (index == -1)
			{
				opserr << "EQPath::ChooseBestRoot() - failed to choose best root\n";
				free(roots);
				return -1;
			}
			dLambda = roots[index];

		}
	}
	
	else if (type == EQPath_Method_GMRD) // generalized minimum residual displacement
	{
		AnalysisModel* am = this->getAnalysisModel();
		Domain* theDomain = theModel->getDomainPtr();
		Element* elePtr;
		double a = 0, b = 0;
		if (theTaggedObjectDic.size() == 0
			&& fillTaggedObjectDic(theDomain, GRMDType, theIDs) == -1)
		{
			opserr << "WARNING EQPath::update()" << endln;
			opserr << "error in fill tagged object dictionary" << endln;
			return -1;
		}

		Vector* uqb = 0, * urb = 0;

		if (GRMDType == GMRD_TYPE_NODE) {
			// iterate over tagged objects
			for (std::map<int, TaggedObject*>::iterator it = theTaggedObjectDic.begin(); it != theTaggedObjectDic.end(); ++it)
			{
				Node* nodePtr = (Node*)it->second;
				if (updateGMRDParameters(nodePtr, uq, ur, uqb, urb, &a, &b) == -1) return -1;
			}
		}
		else if (GRMDType == GMRD_TYPE_ELEMENT) {
			for (std::map<int, TaggedObject*>::iterator it = theTaggedObjectDic.begin(); it != theTaggedObjectDic.end(); ++it)
			{
				Element* elementPtr = (Element*)it->second;
				int numOfNodes = elementPtr->getNumExternalNodes();
				Node** nodes = elementPtr->getNodePtrs();
				int size = 0;
				for (int i = 0; i < numOfNodes; i++)
				{
					Node* nodePtr = nodes[i];
					if (updateGMRDParameters(nodePtr, uq, ur, uqb, urb, &a, &b) == -1) return -1;
				}
			}
		}

		if (b == 0) {
			opserr << "EQPath::update() - divide by zero [b == 0]\n";
			return -1;
		}

		//opserr << "GMRD A:" << a << " B:" << b << endln;
		dLambda = -a / b;
	}
	else if (type == EQPath_Method_PEP) // pep
	{
		double a = (*uq) ^ (*uq);
		double b = (*du + *ur + *ur) ^ (*uq);
		double c = (*du + *ur) ^ (*ur);
		double delta = b * b - 4 * a * c;
		if (delta < 0) {
			opserr << "EQPath::update() - unable to find root [delta < 0]\n";
			return -1;
		}
		else {
			if (delta == 0) {
				dLambda = -b / 2 / a;
			}
			else {

			}
		}

		dLambda = -a / b;
	}
	else {
		opserr << "WARNING EQPath::update() ";
		opserr << "Unknown update method has been set\n";
		return -1;
	}



	Vector* su = new Vector(ur->Size());
	// determine delta U(i)
	(*su) = (*ur);
	su->addVector(1.0, *uq, dLambda);

	(*du) += (*su);
	dl += dLambda;

	double currentLambda = theModel->getCurrentDomainTime();
	currentLambda += dLambda;
	// update the model
	theModel->incrDisp(*su);
	theModel->applyLoadDomain(currentLambda);
	int flag = 1;

	if (theModel->updateDomain() < 0) {
		opserr << "EQPath::update - model failed to update for new dU\n";
		return -1;
	}

	// set the X soln in linearSOE to be deltaU for convergence Test
	theLinSOE->setX(*su);
	nitr++;

	return 0;
}

int
EQPath::updateGMRDParameters(Node* nodePtr, Vector* uq, Vector* ur, Vector* uqb, Vector* urb, double* a, double* b) {

	DOF_Group* dofGroup = nodePtr->getDOF_GroupPtr();
	if (dofGroup == 0) {
		opserr << "EQPath::update() - fe graph of node " << nodePtr->getTag() << " not found\n";
		return -1;
	}

	const ID& id = dofGroup->getID();
	int size = theDofIDs != 0 ? theDofIDs->Size() : id.Size();

	// resize vectors if needed
	if (uqb == 0 || uqb->Size() != size || urb == 0 || urb->Size() != size)
	{
		if (uqb != 0)
			delete uqb;
		if (urb != 0)
			delete urb;

		uqb = new Vector(size);
		urb = new Vector(size);
	}

	uqb->Zero();
	urb->Zero();

	for (int i = 0; i < size; i++)
	{
		int index = theDofIDs != 0 ? id(theDofIDs->operator()(i) - 1) : id(i);
		if (index == -1)
			continue;
		uqb->operator[](i) = uq->operator[](index);
		urb->operator[](i) = ur->operator[](index);
	}

	(*a) += (*urb) ^ (*uqb);
	(*b) += (*uqb) ^ (*uqb);
	return 0;
}

int
EQPath::fillTaggedObjectDic(Domain* theDomain, GMRD_TYPE type, ID* ids) {
	theTaggedObjectDic.clear();
	if (type == GMRD_TYPE_NODE) {
		if (ids != 0)
			for (int i = 0; i < ids->Size(); i++)
			{
				int tag = ids->operator[](i);
				Node* nodePtr = theDomain->getNode(tag);
				if (nodePtr != 0)
					theTaggedObjectDic[tag] = nodePtr;
			}
		else
		{
			NodeIter& nodeIter = theDomain->getNodes();
			Node* nodePtr = 0;
			while ((nodePtr = nodeIter()) != 0)
			{
				theTaggedObjectDic[nodePtr->getTag()] = nodePtr;
			}
		}
		return 0;
	}
	else if (type == GMRD_TYPE_ELEMENT) {
		if (ids != 0)
			for (int i = 0; i < ids->Size(); i++)
			{
				int tag = ids->operator[](i);
				Element* elementPtr = theDomain->getElement(tag);
				if (elementPtr != 0)
					theTaggedObjectDic[tag] = elementPtr;
			}
		else
		{
			ElementIter& elementIter = theDomain->getElements();
			Element* elementPtr = 0;
			while ((elementPtr = elementIter()) != 0)
			{
				theTaggedObjectDic[elementPtr->getTag()] = elementPtr;
			}
		}
		return 0;
	}
	else
	{
		opserr << "GMRD type not implemented - " << type;
		return -1;
	}
}

int
EQPath::chooseBestDLambda(double* dlambdas, int count, Vector* du, Vector* uq, Vector* ur) {
	if (dlambdas == 0 || count == 0) {
		opserr << "EQPath::ChooseBestRoot() - failed to send the data\n";
		return -1;
	}
	double aa1 = (*du) ^ (*ur);
	double aa2 = (*du) ^ (*du);
	double aa3 = (*du) ^ (*uq);

	int index = 0;
	double bestRoot = dlambdas[0];
	double maxCos = aa1 + aa2 + dlambdas[0] * aa3;
	for (int i = 1; i < count; i++)
	{
		double cos = aa1 + aa2 + dlambdas[i] * aa3;
		if (cos > maxCos)
		{
			index = i;
			maxCos = cos;
		}
	}
	return index;
}

int 
EQPath::domainChanged(void)
{
	// we first create the Vectors needed
	AnalysisModel* theModel = this->getAnalysisModel();
	LinearSOE* theLinSOE = this->getLinearSOE();
	if (theModel == 0 || theLinSOE == 0) {
		opserr << "WARNING EQPath::update() ";
		opserr << "No AnalysisModel or LinearSOE has been set\n";
		return -1;
	}

	int size = theModel->getNumEqn(); // ask model in case N+1 space

	if (uq == 0 || uq->Size() != size) { // create new Vector
		if (uq != 0)
			delete uq;   // delete the old
		uq = new Vector(size);
		if (uq == 0 || uq->Size() != size) { // check got it
			opserr << "FATAL EQPath::domainChanged() - ran out of memory for";
			opserr << " uq Vector of size " << size << endln;
			exit(-1);
		}
	}

	if (du == 0 || du->Size() != size) { // create new Vector
		if (du != 0)
			delete du;   // delete the old
		du = new Vector(size);
		if (du == 0 || du->Size() != size) { // check got it
			opserr << "FATAL EQPath::domainChanged() - ran out of memory for";
			opserr << " du Vector of size " << size << endln;
			exit(-1);
		}
	}


	if (ur == 0 || ur->Size() != size) { // create new Vector
		if (ur != 0)
			delete ur;   // delete the old
		ur = new Vector(size);
		if (ur == 0 || ur->Size() != size) { // check got it
			opserr << "FATAL EQPath::domainChanged() - ran out of memory for";
			opserr << " deltaU Vector of size " << size << endln;
			exit(-1);
		}
	}

	if (q == 0 || q->Size() != size) {
		if (q != 0)
			delete q;
		q = new Vector(size);
		if (q == 0 || q->Size() != size) {
			opserr << "FATAL EQPath::domainChanged() - ran out of memory for";
			opserr << " q Vector of size " << size << endln;
			exit(-1);
		}
	}

	// now we have to determine phat
	// do this by incrementing lambda by 1, applying load
	// and getting phat from unbalance.
	double currentLambda = theModel->getCurrentDomainTime();
	currentLambda += 1.0;
	theModel->applyLoadDomain(currentLambda);
	this->formUnbalance(); // NOTE: this assumes unbalance at last was 0
	(*q) = theLinSOE->getB();
	currentLambda -= 1.0;
	theModel->setCurrentDomainTime(currentLambda);


	// check there is a reference load
	int haveLoad = 0;
	for (int i = 0; i < size; i++)
		if ((*q)(i) != 0.0) {
			haveLoad = 1;
			i = size;
		}

	if (haveLoad == 0) {
		opserr << "WARNING ArcLength::domainChanged() - zero reference load";
		return -1;
	}


	return 0;
}

int
EQPath::sendSelf(int cTag,
		    Channel &theChannel)
{
	Vector data(3);
	data(0) = arclen;
	data(1) = dl;
	data(2) = sign_n0;

	if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
		opserr << "EQPath::sendSelf() - failed to send the data\n";
		return -1;
	}
	return 0;
}


int
EQPath::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	Vector data(3);
	if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
		opserr << "EQPath::sendSelf() - failed to send the data\n";
		return -1;
	}

	// set the data

	arclen = data(0);
	dl = data(1);
	sign_n0 = data(2);

	return 0;
}

void
EQPath::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double cLambda = theModel->getCurrentDomainTime();
	s << "\t EQPath - currentLambda: " << cLambda <<"\n";
	s << "\t EQPath - arcLength: " << arclen <<"\n";
	s << "\t EQPath - sign: " << sign_n0 <<"\n";
    } else 
	s << "\t EQPath - no associated AnalysisModel\n";
}

