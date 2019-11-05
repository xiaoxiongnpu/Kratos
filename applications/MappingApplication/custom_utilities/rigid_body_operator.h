//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

//#if !defined(KRATOS_BEAM_MAPPER_H_INCLUDED )
//#define  KRATOS_BEAM_MAPPER_H_INCLUDED
//
//// System includes
//#include <tuple>
//
//// External includes
//
//// Project includes
//#include "includes/define.h"
//
//namespace Kratos
//{
///@name Kratos Classes
///@{

// This is an auxiliary class for the Beam Mapper
//class RigidBodyOperator 
//{
//public:
//
//    typedef array_1d<double, 3>  PointType;
//
//    KinematicMotion();
//
//    virtual ~KinematicMotion();
//    
//    void resetTranslation();
//
//    void resetRotation();
//
//    const double* getTranslationVector() const;
//
//    const double* getRotationMatrix() const;
//
//    void getRotationVector(double* rot) const;
//
//    void addKinematicMotion(const RigidBodyOperator* rigidBodyOperator);
//
//    RigidBodyOperator* newInverse();
//
//    void addTranslation(const double *_translationVector);
//
//    void addRotation(const double* axis, bool normalized, double angle);
//
//    void addRotation(const double* vec1, const double* vec2, bool normalized);
//
//    void addRotation(const double* xAxisNew, const double *yAxisNew, const double *zAxisNew, bool normalized);
//
//    void addRotation(const double* _rotationMatrix);
//
//    void apply(PointType* Coordinates) const;
//
//    void checkRotationCorrectness() const;
//
//    void print() const;
//
//protected:
//    double *rotationMatrix;
//    double *translationVector;
//
//    static void normalizeVector(double* vector);
//
//    static void matrixProduct(const double* A, double* B);
//
//    static void copyMatrix(const double* M, double* MCopy);
//
//    static void copyVector(const double *v, double* vCopy);
//
//    static double vectorProduct(const double* v1, const double* v2);
//
//    static void matrixVectorProduct(const double* M, double* v):
//
//    static void vectorAddition(double* v1, const double* v2);
//
//    static void matrixTranspose(double* M);
//
//    static void vectorInverse(double* v);
//
//    static void printMatrix(const double* M);
//
//    static void printVector(const double* v);
//};
//
//
//}; // Class BeamMapper
//
/////@} addtogroup block
//}  // namespace Kratos.
//
//#endif // KRATOS_BEAM_MAPPER_H_INCLUDED  defined