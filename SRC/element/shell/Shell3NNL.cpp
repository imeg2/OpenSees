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
                                                                        
// $Revision: 1.8 $
// $Date: 2007-02-02 01:30:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Shell3NNL.cpp,v $
                                                                        
// Written: msa // msa_imeg@yahoo.com 
// Created: 11/11
//
// Description: This file contains the implementation for the Shell3NNL class.
//
// What: "@(#) Shell3NNL.cpp, revA"

#include "Shell3NNL.h"
#include <Information.h>
#include <ElementResponse.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include "Vector.h"
#include "FileStream.h"
#include "DomainComponent.h"

static int NUM_NODE =3;
static int NUM_DOF  =18;

Matrix Shell3NNL::theMatrix(18,18);
Vector Shell3NNL::theVector(18);

#ifdef _WIN32

extern "C" int STF_SHELLNL(double *kt0, double *rt0,double *e,
				    double *nu,double *psig0,double *th);

extern "C" int RST_SHELLNL(double *f0,double *rt0,double *e,
		  double *nu,double *psig0,double *th);

extern "C" int UPDATE_FORCES(double *rl0,double *p0,double *e,
		  double *nu,double *psig0,double *th);

#else

extern "C" void stf_shellnl_(double *kt0, double *rt0,double *e,
					   double *nu,double *psig0,double *th);

extern "C" void rst_shellnl_(double *f0,double *rt0,double *e,
					   double *nu,double *psig0,double *th);

extern "C" void update_forces_(double *rl0,double *p0,double *e,
						double *nu,double *psig0,double *th);

#endif

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_NewShell3NNLElement()
{
    Element *theElement = 0;

    int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs != 7) {
	   opserr << "Invalid Args want: element Shell3NNL $tag $iNode $jNode $kNode $E $nu $t \n";
	   exit(-1);	
    }

    int    iData[4];
    double E = 0.0;
    double nu = 0.0;
    double t = 0.0;
    
    int numData = 4;
    if (OPS_GetInt(&numData, iData) != 0) {
	   opserr << "WARNING invalid integer (tag, iNode, jNode, kNode) in element Shell3NNL " << endln;
	   return 0;
    }

    numData = 1;
    if (OPS_GetDouble(&numData, &E) != 0) {
	   opserr << "WARNING: Invalid A: element Shell3NNL " << iData[0] << 
		  " Invalid Modulus of elasticity\n";
	   return 0;	
    }

    numData = 1;
    if (OPS_GetDouble(&numData, &nu) != 0) {
	   opserr << "WARNING: Invalid matTag: element Shell3NNL " << iData[0] << 
		  " Invalid poisson's ratio\n";
	   return 0;
    }

    numData = 1;
    if (OPS_GetDouble(&numData, &t) != 0) {
	   opserr << "WARNING: Invalid matTag: element Shell3NNL " << iData[0] << 
		  " Invalid thickness\n";
	   return 0;
    }

    theElement = new Shell3NNL(iData[0],iData[1],iData[2],iData[3],E,nu,t);

    if (theElement == 0) {
	   opserr << "WARNING: out of memory: element Shell3NNL " << iData[0] << 
		  " $tag $iNode $jNode $kNode $E $nu $t\n";
    }

    return theElement;
}



// constructors:
Shell3NNL::Shell3NNL(int tag,int n1,int n2,int n3,double E,double nu1,double t)
 :Element(tag,ELE_TAG_Shell3NNL),     
  connectedExternalNodes(3),e(E),nu(nu1),th(t),lu(0),cu(0),rt(0),rl(0),rlf(0),psit(0),psih(0)
{
    //opserr.setPrecision(20);

    if (connectedExternalNodes.Size()!=3)
    {
	   opserr << "FATAL Shell3NNL::Shell3NNL - " <<  tag << "failed to create an ID of size 3\n";
	   exit(-1);
    }
    
    updated=false;

    connectedExternalNodes(0)=n1;
    connectedExternalNodes(1)=n2;
    connectedExternalNodes(2)=n3;

    // set node pointers to NULL
    for (int i=0; i<3; i++)
	   theNodes[i] = 0;
    
    this->lu = new double[18];
    this->cu = new double[18];
    this->psih = new double[12];
    this->psit = new double[12];
    this->rl = new double[9];
    this->rlf = new double[9];
    this->rt = new double[9];
    for(int i=0;i<18;i++)
    {
	   lu[i]=0.0;
	   cu[i]=0.0;
    }

    for(int i=0;i<12;i++)
    {
	   psih[i]=0.0;
	   psit[i]=0.0;
    }

    for(int i=0;i<9;i++)
    {
	   rt[i]=0.0;
	   rl[i]=0.0;
	   rlf[i]=0.0;
    }
}

Shell3NNL::Shell3NNL()
 :Element(0,ELE_TAG_Shell3NNL),     
  connectedExternalNodes(2),lu(0),cu(0),rt(0),rl(0),rlf(0),psit(0),psih(0)
{

}

//  destructor:
Shell3NNL::~Shell3NNL()
{
    
    if(lu!=0)
	   delete [] lu;
    if(cu!=0)
	   delete [] cu;
    if(rt!=0)
	   delete [] rt;
    if(rl!=0)
	   delete [] rl;
    if(rlf!=0)
	   delete [] rlf;
    if(psit!=0)
	   delete [] psit;
    if(psih!=0)
	   delete [] psih;
}

int
Shell3NNL::getNumExternalNodes(void) const
{
    return NUM_NODE;
}

const ID &
Shell3NNL::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
Shell3NNL::getNodePtrs(void)
{
  return theNodes;
}

int
Shell3NNL::getNumDOF(void) 
{
    return NUM_DOF;
}


void
Shell3NNL::setDomain(Domain *theDomain)
{
    if (theDomain == 0) {
	   theNodes[0] = 0;
	   theNodes[1] = 0;
	   theNodes[1] = 0;
	   return;
    }
    
    int *n=new int[3];
    n[0]=this->connectedExternalNodes(0);
    n[1]=this->connectedExternalNodes(1);
    n[2]=this->connectedExternalNodes(2);

    this->theNodes[0]=theDomain->getNode(n[0]);
    this->theNodes[1]=theDomain->getNode(n[1]);
    this->theNodes[2]=theDomain->getNode(n[2]);

    for (int i=0;i<3;i++)
    {
	   if (theNodes[i] == 0) 
	   {
		  opserr <<"Shell3NNL::setDomain() - Shell3NNL " << this->getTag() << " node " << n[i] <<
			 "does not exist in the model\n"<< endln;
		  exit(-1);
	   }
    }
    
    // now determine the number of dof and the dimesnion    
    
    for (int i=0;i<3;i++)
    {
	   if (theNodes[i]->getNumberDOF() != 6) 
	   {
		  opserr <<"WARNING Shell3NNL::setDomain(): nodes " << theNodes[i]->getTag() << 
			 "have invalid dof ("<<theNodes[i]->getNumberDOF() <<" != 6) " << endln;
		  exit(-1);
	   }
    }

    this->DomainComponent::setDomain(theDomain);

}   	 


int
Shell3NNL::commitState()
{
    
    
  int retVal = 0;

  // call the base class method
  retVal = this->Element::commitState();
  if (retVal < 0) {
    opserr << "Shell3NNL::commitState() - failed in base class\n";
    return retVal;
  }

  for(int i=0;i<12;i++)
  {
	 psih[i]=psit[i];
  }

  update();
  return retVal;
}

int
Shell3NNL::revertToLastCommit()
{
    for(int i=0;i<12;i++)
    {
	   psit[i]=psih[i];
    }
  return 0;
}

int
Shell3NNL::revertToStart()
{
    
    for(int i=0;i<18;i++)
    {
	   lu[i]=0.0;
    }

    for(int i=0;i<12;i++)
    {
	   psih[i]=0.0;
	   psit[i]=0.0;
    }

    for(int i=0;i<9;i++)
    {
	   rt[i]=0.0;
	   rl[i]=0.0;
	   rlf[i]=0.0;
    }
  return 0;
}


int
Shell3NNL::update(void)
{
    updated=false;
    for(int i=0;i<18;i++)
	   lu[i]=cu[i];

    for (int i = 0; i < 3; i++)
    {
	   for (int j = 0; j < 6; j++)
	   {
		  int index = i * 6 + j;
		  cu[index]=(theNodes[i]->getTrialDisp())(j);
	   }
    }
    
    
    
    
   
    return 0;

}


const Matrix &
Shell3NNL::getTangentStiff(void)
{
    for (int i = 0; i < 3; i++)
    {
	   for (int j = 0; j < 3; j++)
	   {
		  int index = i * 3 + j;
		  int index1 = i * 6 + j;
		  rl[index]=(theNodes[i]->getCrds())(j) + lu[index1];
		  rt[index]=(theNodes[i]->getCrds())(j) + cu[index1];
	   }
    }
    double *p0 = new double[18];
    for (int i = 0; i < 3; i++)
    {
	   for (int j = 0; j < 6; j++)
	   {
		  int index = i * 6 + j;
		  p0[index] = cu[index] - lu[index];
	   }
    }

    if(updated==false)
    {
	   UPDATE_FORCES(rl,p0,&e,&nu,psit,&th);   
	   updated=true;
    }

    
   

    //double *Kt=new double[NUM_DOF*NUM_DOF];
    //STF_SHELLNL(Kt,rt, &e, &nu,psit, &th);
    //theMatrix.setData(Kt,NUM_DOF,NUM_DOF); 
    double *Kt=new double[NUM_DOF*NUM_DOF];
    STF_SHELLNL(Kt,rt, &e, &nu,psit, &th);

    for(int i=0;i<NUM_DOF*NUM_DOF;i++)
    {
	   int row=i/NUM_DOF;
	   int col=i-row*NUM_DOF;
	   theMatrix(col,row)=Kt[i];
    }
    delete [] Kt;
    delete [] p0;


    return theMatrix;
}

const Vector &
Shell3NNL::getResistingForce()
{

    for (int i = 0; i < 3; i++)
    {
	   for (int j = 0; j < 3; j++)
	   {
		  int index = i * 3 + j;
		  int index1 = i * 6 + j;
		  rl[index]=(theNodes[i]->getCrds())(j) + lu[index1];
		  rt[index]=(theNodes[i]->getCrds())(j) + cu[index1];
	   }
    }
    double *p0 = new double[18];
    for (int i = 0; i < 3; i++)
    {
	   for (int j = 0; j < 6; j++)
	   {
		  int index = i * 6 + j;
		  p0[index] = cu[index] - lu[index]; // 
	   }
    }
    if(updated==false)
    {
	   UPDATE_FORCES(rl,p0,&e,&nu,psit,&th);   
	   updated=true;
    }
    
    
    double *f=new double[NUM_DOF];
    RST_SHELLNL(f,rt, &e, &nu,psit, &th);
    
    for(int i=0;i<NUM_DOF;i++)
    {

	   theVector(i)=f[i];
    }
    delete [] f;
    delete [] p0;

    //double *f=new double[NUM_DOF];
    //RST_SHELLNL(f,rt, &e, &nu,psit, &th);
    //theVector.setData(f,NUM_DOF);
    return theVector;
}



const Matrix &
Shell3NNL::getInitialStiff(void)
{
    double *rl1,*rt1,*lu1,*psit1;

    for (int i = 0; i < 3; i++)
	   for (int j = 0; j < 3; j++)
	   {
		  int index = i * 3 + j;
		  int index1 = i * 6 + j;
		  rl1[index]=(theNodes[i]->getCrds())(j) + lu1[index1]*0;
		  rt1[index]=(theNodes[i]->getCrds())(j) + (theNodes[i]->getDisp()+theNodes[i]->getTrialDisp())(j)*0;
	   }

    double *p0 = new double[18];
    for (int i = 0; i < 3; i++)
		for (int j = 0; j < 6; j++)
		{
			    int index = i * 6 + j;
			    p0[index] = (theNodes[i]->getDisp()+theNodes[i]->getTrialDisp())(j)*0 - lu1[index]*0;
			    lu1[index]=(theNodes[i]->getDisp()+theNodes[i]->getTrialDisp())(j)*0;
		}
    
	
	UPDATE_FORCES(rl1,p0,&e,&nu,psit1,&th);

	double *Kt=new double[NUM_DOF*NUM_DOF];
	STF_SHELLNL(Kt,rt1, &e, &nu,psit1, &th);

	for(int i=0;i<NUM_DOF*NUM_DOF;i++)
	{
	    int row=i/NUM_DOF;
	    int col=i-row*NUM_DOF;
	    theMatrix(row,col)=Kt[i];
	}
	delete [] Kt;
	delete [] p0,rl1,rt1,lu1,psit1;;

	return theMatrix;

}
    
void 
Shell3NNL::zeroLoad(void)
{
  return;
}


int 
Shell3NNL::addLoad(const Vector &addP)
{
  return 0;
}


int 
Shell3NNL::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}



const Vector &
Shell3NNL::getResistingForceIncInertia()
{	
    return Shell3NNL::getResistingForce();
}


int
Shell3NNL::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
Shell3NNL::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}


int
Shell3NNL::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
  return 0;
}


void
Shell3NNL::Print(OPS_Stream &s, int flag)
{
  return;
}


Response*
Shell3NNL::setResponse(const char **argv, int argc, OPS_Stream &S)
{
	return this->Element::setResponse(argv, argc, S);
}


int 
Shell3NNL::getResponse(int responseID, Information &eleInfo)
{
	return this->Element::getResponse(responseID, eleInfo);
}


int
Shell3NNL::setParameter(const char **argv, int argc, Parameter &param)
{
  return 0;
}
    

int
Shell3NNL::updateParameter(int parameterID, Information &info)
{
  return -1;
}
