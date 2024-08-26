
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


#ifndef EQPath_h
#define EQPath_h

#include <StaticIntegrator.h>
#include <map>
#include <iostream>
#include <fstream>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class ID;
class Node;
class TaggedObject;

#define SIGN_LAST_STEP      1
#define CHANGE_DETERMINANT  2

#define GMRD_TYPE int
#define GMRD_TYPE_NODE 1
#define GMRD_TYPE_ELEMENT 2


#define EQPath_Method int
#define EQPath_Method_MRD 1 
#define EQPath_Method_NP 2
#define EQPath_Method_UNP 3
#define EQPath_Method_CYL 4
#define EQPath_Method_MNP 5
#define EQPath_Method_GDC 6
#define EQPath_Method_MUNP 7
#define EQPath_Method_MGDC 8
#define EQPath_Method_GMRD 9
#define EQPath_Method_PEP 10
#define EQPath_Method_MRW 11 
#define EQPath_Method_ZRW 12 

class EQPath : public StaticIntegrator
{
  public:
    EQPath(double arcLeng,int type, ID* ids, ID* dofIDs, int gmrdType);

    ~EQPath();

    int newStep(void);    
    int update(const Vector &deltaU);
    int domainChanged(void);
    
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:
      int updateGMRDParameters(Node* nodePtr, Vector* uq, Vector* ur, Vector* uqb, Vector* urb, double* a, double* b);
      int chooseBestDLambda(double* dlambdas, int count, Vector* du, Vector* uq, Vector* ur);
      int fillTaggedObjectDic(Domain* theDomain, GMRD_TYPE type, ID* ids);
      ID* theDofIDs;
      ID* theIDs;
      GMRD_TYPE GRMDType;
    double arclen,dl;
    double sign_n0; // current sign of increment
    int type, nitr;
    Vector* du, * du0, * uq, * uq_n0, * uq_n1, * ur;
    Vector* q;
    std::map<int, TaggedObject*> theTaggedObjectDic;
    int step, iter;
    
};

#endif


