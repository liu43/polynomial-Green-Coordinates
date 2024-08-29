#pragma once
#include <vector>
#include <cmath>
#include <map>
#include <iostream>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace std;
class BezierCurve
{
protected:
	int m_Degree;
	vector<OpenMesh::Vec3d> m_CtrlPoints;

public:
	BezierCurve():m_Degree(0),m_CtrlPoints(vector<OpenMesh::Vec3d>(1)){}
	BezierCurve(const vector<OpenMesh::Vec3d>& p_CtrlPoints);
	

public:
	vector<OpenMesh::Vec3d> GetCtrlPoints() const { return m_CtrlPoints; }

public:
	OpenMesh::Vec3d Evaluate(const double& para) const;
};



