#pragma once
#include "Particle.h"
class RigidBody
{
public:
	// Constructors
	RigidBody();
	RigidBody(Particle Collection[50], double MomentOfInertia[3][3]);
	// Retrievers

	double GetVelocity(int i);
	double GetQuaternion(int i);
	double GetCenterOfMass(int i);
	double GetBodyRate(int i);
	double GetBodyTorque(int i);
	double GetQuaternionVelocity(int i);
	double GetMass();
	double GetParticleMass(int i);
	double GetMomentOfInertia(int i, int j);




	// Editors


	void SetVelocity(double NewVelocity[3], double OldVelocity[3]);
	void SetQuaternion(double NewQuaternion[4], double OldQuaternion[4]);
	void SetCenterOfMass(double NewCenterOfMass[3], double OldCenterOfMass[3]);
	void SetBodyRate(double NewBodyRate[3], double OldBodyRate[3]);
	void SetBodyTorque(double NewBodyTorque[3], double OldBodyTorque[3]);
	void SetQuaternionVelocity(double NewQVelocity[4], double OldQVelocity[4]);





	//
	// Mathematical Operations
	//

	// Swapping Back and Forth Between Quaternions, Vectors, and Operators
	void VectorToQuaternion(double Vector[3], double Quaternion[4]);
	void QuaternionToVector(double Quaternion[4], double OutputVector[3]);

	// Scaling Vectors, Operators, and Quaternions
	void QScalar(double a, double Q[4]);
	void VScalar(double a, double V[3]);
	void VScalar(double a, double V[3], double OutV[3]);
	void QScalar(double a, double Q[4], double OutQ[4]);
	void MScalar(double a, double M[3][3]);
	void MScalar(double a, double M[3][3], double OutM[3][3]);
	void VNormalize(double Vector[3]);
	void QNormalize(double Vector[4]);



	// Sums Of Vectors Quaternions And Operators
	void VSum(double Vec1[3], double Vec2[3]);
	void QSum(double Q1[4], double Q2[4]);
	void VSum(double V1[3], double V2[3], double ReturnV[3]);
	void QSum(double Q1[4], double Q2[4], double ReturnQ[4]);
	void MixSum(double Q[4], double V[3]);
	void MixSum(double Q[4], double V[3], double OutQ[4]);



	//Vector, Quaternion, and Operator Products
	double InnerProduct(double Vec1[3], double Vec2[3]);
	void CrossProduct(double Vec1[3], double Vec2[3], double ReturnVec[3]);
	void QProduct(double Q1[4], double Q2[4], double OutputQ[4]);
	void QProduct(double Q1[4], double Q2[4]);
	void MixProduct(double Q[4], double V[3], double OutQ[4]);
	void MixProduct(double Q[4], double V[3]);
	void MatrixOperation(double M[3][3], double V[3]);
	void MatrixOperation(double M[3][3], double V[3], double OutV[3]);
	void OperatorComposition(double M1[3][3], double M2[3][3], double OutM[3][3]);
	void OperatorComposition(double M1[3][3], double M2[3][3]);
	double Det(double M1[3][3]);
	void AdjugateMatrix(double M1[3][3], double M2[3][3]);
	double CoFactor(double M[3][3], int i, int j);

	//
	// Other Helpful and Necessary Functions
	//

	void InvertMatrix(double M1[3][3], double M2[3][3]);
	void CalcRFromQ(double Quat[4], double R[3][3]);
	void CalcRStar(double R[3][3], double RS[3][3]);



	// The Big Kahuna


	void Timestep(double dt, double COM[3], double V[3], double VDot[3], double Q[4], double QDot[4], double BR[3], double BRDot[3], double BT[3], double MOI[3][3], double Force[3], double m, double MoiInverse[3][3]);
	void Timestep(double dt);
private:

	// Particles that "Make Up" the Rigid body
	Particle Particles[50];
	Particle CenteredParticles[50];

	//Constants For the RigidBody
	double BodyMass;
	double MomentOfInertia[3][3];
	double MomentOfInertiaInverse[3][3];

	//Force 
	double TotalForce[3]; // F


	// StateVariables (What we will print out to Document)
	double CenterOfMass[3]; // r  
	double Velocity[3]; // r dot
	double Quaternion[4]; // Quaternion q
	double BodyRate[3]; // Omega Star


	// Quantities involved in State Variable Derivatives
	double CenterAcceleration[3]; // r double dot
	double QuaternionVelocity[4]; // q dot
	double BodyRateDerivative[3]; // Omega star dot
	double BodyTorque[3]; // tau star 
	double RotationMatrix[3][3]; // 
	double RStar[3][3];
};

