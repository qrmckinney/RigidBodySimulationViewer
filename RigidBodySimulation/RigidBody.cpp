
#include <iostream>
#include "RigidBody.h"

using namespace std;

// Constructors and such

RigidBody::RigidBody()
{

}
RigidBody::RigidBody(Particle Collection[50], double Moment[3][3])
{
	BodyMass = 0;
	for (int i = 0; i < 50; i++)
	{
		BodyMass += Collection[i].GetMass();
	}
	cout << "The Body's Total Mass is " << BodyMass << endl; 
	
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			MomentOfInertia[i][j] = Moment[i][j];
		}

		cout << "Please Enter The Desired ";
		if (i == 0)
		{
			cout << "X ";
		}
		else if (i == 1)
		{
			cout << "Y ";
		}
		else
		{
			cout << "Z ";
		}

		

		cout << "Coordinate of the Body's Initial " << endl;
		cout << "Center Of Mass, Velocity, and Body Rate, in that order, separated by a space." << endl;
		cin >> CenterOfMass[i] >> Velocity[i] >> BodyRate[i];
		BodyTorque[i] = 0;
		TotalForce[i] = 0;
		Quaternion[i] = 0;
		

	}
	TotalForce[2] = -9.8*BodyMass;
	
	InvertMatrix(MomentOfInertia, MomentOfInertiaInverse);
	Quaternion[3] = 0;
	Quaternion[0] = 1;
	VScalar(1 / BodyMass, TotalForce, CenterAcceleration);
	QProduct(Quaternion, BodyRate, QuaternionVelocity);
	QScalar(.5, QuaternionVelocity);
	double NewDude[3];
	double OtherNewDude[3];
	MatrixOperation(MomentOfInertia, BodyRate, NewDude);
	CrossProduct(BodyRate, NewDude, OtherNewDude);
	VScalar(-1, OtherNewDude);
	VSum(BodyTorque, OtherNewDude);
	MatrixOperation(MomentOfInertiaInverse, OtherNewDude, BodyRateDerivative);
	
}

// Retrievers and Such

double RigidBody::GetVelocity(int i)
{
	return Velocity[i];
}
double RigidBody::GetQuaternion(int i)
{
	return Quaternion[i];
}
double RigidBody::GetCenterOfMass(int i)
{
	return CenterOfMass[i];
}
double RigidBody::GetBodyRate(int i)
{
	return BodyRate[i];
}
double RigidBody::GetBodyTorque(int i)
{
	return BodyTorque[i];
}
double RigidBody::GetQuaternionVelocity(int i)
{
	return QuaternionVelocity[i];
}
double RigidBody::GetMass()
{
	return BodyMass;
}
double RigidBody::GetParticleMass(int i)
{
	return Particles[i].GetMass();
}
double RigidBody::GetMomentOfInertia(int i, int j)
{
	return MomentOfInertia[i][j];
}

// Editors

void RigidBody::SetVelocity(double NewVelocity[3], double OldVelocity[3])
{
	for (int i = 0; i < 3; i++)
	{
		OldVelocity[i] = NewVelocity[i];
	}
}
void RigidBody::SetQuaternion(double NewQuaternion[4], double OldQuaternion[4])
{
	for (int i = 0; i < 4; i++)
	{
		OldQuaternion[i] = NewQuaternion[i];
	}
}
void RigidBody::SetCenterOfMass(double NewCenterOfMass[3], double OldCenterOfMass[3])
{
	for (int i = 0; i < 3; i++)
	{
		OldCenterOfMass[i] = NewCenterOfMass[i];
	}
}
void RigidBody::SetBodyRate(double NewBodyRate[3], double OldBodyRate[3])
{
	for (int i = 0; i < 3; i++)
	{
		OldBodyRate[i] = NewBodyRate[i];
	}
}
void RigidBody::SetBodyTorque(double NewBodyTorque[3], double OldBodyTorque[3])
{
	for (int i = 0; i < 3; i++)
	{
		OldBodyTorque[i] = NewBodyTorque[3];
	}
}
void RigidBody::SetQuaternionVelocity(double NewQVelocity[4], double OldQVelocity[4])
{
	for (int i = 0; i < 4; i++)
	{
		OldQVelocity[i] = NewQVelocity[i];
	}
}

//
//
// Basic Mathematical Operations
//
//

// Sum Of Quaterniions and Vectors
void RigidBody::QSum(double Q1[4], double Q2[4], double ReturnQ[4])
{
	for (int i = 0; i < 4; i++)
	{
		ReturnQ[i] = Q1[i] + Q2[i];

	}
}
void RigidBody::QSum(double Q1[4], double Q2[4])
{
	for (int i = 0; i < 4; i++)
	{
		Q1[i] += Q2[i];
	}
}
void RigidBody::VSum(double V1[3], double V2[3], double ReturnV[3])
{
	for (int i = 0; i < 3; i++)
	{
		ReturnV[i] = V1[i] + V2[i];
	}
}
void RigidBody::VSum(double V1[3], double V2[3])
{
	for (int i = 0; i < 3; i++)
	{
		V1[i] += V2[i];
	}
}
void RigidBody::MixSum(double Q[4], double V[3])
{
	for (int i = 0; i < 3; i++)
	{
		Q[i + 1] += V[i];
	}
}
void RigidBody::MixSum(double Q[4], double V[3], double OutQ[4])
{
	for (int i = 0; i < 3; i++)
	{
		OutQ[i + 1] = Q[i + 1] + V[i];
	}
	OutQ[0] = Q[0];
}


// Scalaing Vectors, Operators, and Quaternions
void RigidBody::QScalar(double a, double Q[4], double OutQ[4])
{
	for (int i = 0; i < 4; i++)
	{
		OutQ[i] = a * Q[i];
	}
}
void RigidBody::VScalar(double a, double V[3], double OutV[3])
{
	for (int i = 0; i < 3; i++)
	{
		OutV[i] = a * V[i];
	}
}
void RigidBody::VScalar(double a, double V[3])
{
	for (int i = 0; i < 3; i++)
	{
		V[i] *= a;
	}
}
void RigidBody::QScalar(double a, double Q[4])
{
	for (int i = 0; i < 4; i++)
	{
		Q[i] *= a;
	}
}
void RigidBody::MScalar(double a, double M[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			M[i][j] *= a;
		}
	}
}
void RigidBody::MScalar(double a, double M[3][3], double OutM[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			OutM[i][j] = M[i][j] * a;
		}
	}
}
void RigidBody::VNormalize(double Vector[3])
{
	double a = InnerProduct(Vector, Vector);
	double b = sqrt(a);
	VScalar(1/b, Vector);
}
void RigidBody::QNormalize(double Q[4])
{
	double F[3];
	QuaternionToVector(Q, F);
	double a = Q[0] * Q[0] + InnerProduct(F, F);
	double b = sqrt(a);
	QScalar(1/b, Q);
}


// Swapping between Vectors and Quaternions
void RigidBody::VectorToQuaternion(double Vector[3], double Quaternion[4])
{
	Quaternion[0] = 0;
	for (int i = 0; i < 3; i++)
	{
		Quaternion[i + 1]=Vector[i];

	}
}
void RigidBody::QuaternionToVector(double Quaternion[4], double Vector[3])
{
	for (int i = 0; i < 3; i++)
	{
		Vector[i] = Quaternion[i + 1];
	}
}



// Products Of Various Types
void RigidBody::QProduct(double Q1[4], double Q2[4], double OutputQ[4])
{
	double V1[3];
	double V2[3];
	double Cross[3];
	QuaternionToVector(Q1, V1);
	QuaternionToVector(Q2, V2);
	CrossProduct(V1, V2, Cross);
	VScalar(Q1[0], V2);
	VScalar(Q2[0], V1);


	for (int i = 0; i < 4; i++)
	{
		OutputQ[i] = 0;
	}
	OutputQ[0] = Q1[0] * Q2[0] - InnerProduct(V1, V2);
	MixSum(OutputQ, V1);
	MixSum(OutputQ, V2);
	MixSum(OutputQ, Cross);
	


}
void RigidBody::QProduct(double Q1[4], double Q2[4])
{
	double LocalDude[4];
	QProduct(Q1, Q2, LocalDude);
	for (int i = 0; i < 4; i++)
	{
		Q1[i] = LocalDude[i];

	}
}
void RigidBody::MixProduct(double Q[4], double V[3])
{
	double LocalVariable[4];
	VectorToQuaternion(V, LocalVariable);
	QProduct(Q, LocalVariable);
	for (int i = 0; i < 4; i++)
	{
		Q[i] = LocalVariable[i];
	}

}
void RigidBody::MixProduct(double Q[4], double V[3], double OutQ[4])
{
	double LocalVariable[4];
	VectorToQuaternion(V, LocalVariable);
	QProduct(Q, LocalVariable, OutQ);
}
double RigidBody::InnerProduct(double Vec1[3], double Vec2[3])
{
	double a = 0;
	for (int i = 0; i < 3; i++)
	{
		a += Vec1[i] * Vec2[i];

	}
	return a;
}
void RigidBody::CrossProduct(double Vec1[3], double Vec2[3], double ReturnVec[3])
{
	ReturnVec[0] = Vec1[1] * Vec2[2] - Vec1[2] * Vec2[1];
	ReturnVec[1] = Vec1[2] * Vec2[0] - Vec1[0] * Vec2[2];
	ReturnVec[2] = Vec1[0] * Vec2[1] - Vec1[1] * Vec2[0];
}
void RigidBody::MatrixOperation(double M[3][3], double V[3], double OutV[3])
{
	for (int i = 0; i < 3; i++)
	{
		OutV[i] = 0;
		for (int j = 0; j < 3; j++)
		{
			OutV[i] += M[i][j] * V[j];
		}
	}
}
void RigidBody::MatrixOperation(double M[3][3], double V[3])
{
	double LocalVariable[3];
	MatrixOperation(M, V, LocalVariable);
	for (int i = 0; i < 3; i++)
	{
		V[i] = LocalVariable[i];
	}
}
void RigidBody::OperatorComposition(double M1[3][3], double M2[3][3], double OutM[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			OutM[i][j] = 0;
			for (int k = 0; k < 3; k++)
			{
				OutM[i][j] += M1[i][k] * M2[k][j];
			}

		}
	}
}
void RigidBody::OperatorComposition(double M1[3][3], double M2[3][3])
{
	double LocalVariable[3][3];
	OperatorComposition(M1, M2, LocalVariable);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			M1[i][j] = LocalVariable[i][j];
		}
	}
}
void RigidBody::InvertMatrix(double M1[3][3], double  M2[3][3])
{
	double DummyMatrix[3][3];
	AdjugateMatrix(M1, DummyMatrix);
	MScalar(1 / Det(M1), DummyMatrix, M2);

}
void RigidBody::AdjugateMatrix(double M1[3][3], double M2[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			M2[j][i] = CoFactor(M1, i, j);
		}
	}
}
double RigidBody::CoFactor(double M[3][3], int i, int j)
{
	int a = 1;
	for (int k = 0; k < i + j; k++)
	{
		a *= -1;
	}

	double NewMatrix[2][2];
	for (int m = 0; m < 2; m++)
	{
		for (int n = 0; n < 2; n++)
		{
			if (m < i&&n < j)
			{
				NewMatrix[m][n] = M[m][n];
			}
			else if (m < i)
			{
				NewMatrix[m][n] = M[m][n + 1];
			}
			else if (n < j)
			{
				NewMatrix[m][n] = M[m + 1][n];
			}
			else
			{
				NewMatrix[m][n] = M[m + 1][n + 1];
			}
		}
	}
	double b = NewMatrix[0][0] * NewMatrix[1][1] - NewMatrix[0][1] * NewMatrix[1][0];
	double c = a * b;
	return c;

}
double RigidBody::Det(double M[3][3])
{
	double a = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0]*M[2][2]-M[2][0]*M[1][2]) + M[0][2] * (M[1][0]*M[2][1]-M[2][0]*M[1][1]);

	return a;
}

//
// Physicsey Methods
//

void RigidBody::CalcRFromQ(double Q[4], double R[3][3])
{
	R[0][0] = Q[0] * Q[0] + Q[1] * Q[1] - Q[2] * Q[2] - Q[3] * Q[3];
	R[0][1] = 2 * (Q[1] * Q[2] - Q[0] * Q[3]);
	R[0][2] = 2 * (Q[1] * Q[3] + Q[0] * Q[2]);
	R[1][0] = 2 * (Q[1] * Q[2] + Q[0] * Q[3]);
	R[1][1] = Q[0] * Q[0] - Q[1] * Q[1] + Q[2] * Q[2] - Q[3] * Q[3];
	R[1][2] = 2 * (Q[2] * Q[3] - Q[0] * Q[1]);
	R[2][0] = 2 * (Q[1] * Q[3] - Q[0] * Q[2]);
	R[2][1] = 2 * (Q[2] * Q[3] + Q[0] * Q[1]);
	R[2][2] = Q[0] * Q[0] - Q[1] * Q[1] - Q[2] * Q[2] + Q[3] * Q[3];
}
void RigidBody::CalcRStar(double R[3][3], double RS[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			RS[i][j] = R[j][i];
		}
	}
}


// The Big Kahuna Dude
void RigidBody::Timestep(double dt)
{
	Timestep(dt, CenterOfMass, Velocity, CenterAcceleration, Quaternion, QuaternionVelocity, BodyRate, BodyRateDerivative, BodyTorque, MomentOfInertia, TotalForce, BodyMass, MomentOfInertiaInverse);
}

void RigidBody::Timestep(double dt, double COM[3], double V[3], double VDot[3], double Q[4], double QDot[4], double BR[3], double BRDot[3], double BT[3], double MOI[3][3], double Force[3], double m, double MoiInverse[3][3])
{
	// Calculating New State Vector
	
	// Position
	double dr[3];
	VScalar(dt, V, dr);
	VSum(COM, dr);
	//Velocity Return to this one
	double dv[3];
	VScalar(dt, VDot, dv);
	VSum(V, dv);
	// Quaternion
	double dq[4];
	QScalar(dt, QDot, dq);
	QSum(Q, dq);
	QNormalize(Q);
	// Body Rate
	double dw[3];
	VScalar(dt, BRDot, dw);
	VSum(BR, dw);


	// Calculating Derivative Of New State Vector (Except velocity, which is part of State Vector)

	//VDot
	VScalar(1 / m, Force, VDot);

	//QDot
	QProduct(Q, BR, QDot);
	QScalar(.5, QDot);


	// BRDot
	double NewDude[3];
	double OtherNewDude[3];
	MatrixOperation(MOI, BR, NewDude);
	CrossProduct(BR, NewDude, OtherNewDude);
	VScalar(-1, OtherNewDude);
	VSum(OtherNewDude, BT);
	MatrixOperation(MoiInverse, OtherNewDude, BRDot);

}