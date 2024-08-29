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


namespace MVC {
	struct Point2D {
		Point2D(double x_ = 0.00, double y_ = 0.00) {
			x = x_;
			y = y_;
		}
		double x;
		double y;
	};
	
	Point2D operator + (const Point2D& a, const Point2D& b);
	Point2D operator - (const Point2D& a, const Point2D& b);
	Point2D operator * (const Point2D& a, const Point2D& b);
	Point2D operator * (const Point2D& a, double t);
	Point2D operator / (const Point2D& a, double t);
	Point2D operator / (const Point2D& a, const Point2D& b);
	

	void cubicMVCs(const std::vector<Point2D>& poly, const Point2D& p, std::vector<double>& vCoords, std::vector<double>& gnCoords, std::vector<double>& gtCoords);
	void cubicMVCs(const std::vector<OpenMesh::Vec3d>& poly, const OpenMesh::Vec3d& p, std::vector<double>& vCoords, std::vector<double>& gnCoords, std::vector<double>& gtCoords);
}