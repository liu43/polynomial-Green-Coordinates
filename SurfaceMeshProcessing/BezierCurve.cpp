#include "BezierCurve.h"

BezierCurve::BezierCurve(const vector<OpenMesh::Vec3d>& p_CtrlPoints) {
	m_CtrlPoints = p_CtrlPoints;
	m_Degree = p_CtrlPoints.size() - 1;
}

OpenMesh::Vec3d BezierCurve::Evaluate(const double& para) const {
	double t = para;
	assert(t >= 0.0 && t <= 1.0);
	//de Casteljau Ëã·¨
	vector<OpenMesh::Vec3d> temp_points = m_CtrlPoints;
	
	for (uint i = 1; i <= m_Degree; ++i) {
		for (uint j = 0; j <= m_Degree - i; ++j) {
			
			temp_points[j] = temp_points[j] * (1.0 - t)  + temp_points[j + 1] * t ;
			
		}
	}
	return temp_points[0];
}