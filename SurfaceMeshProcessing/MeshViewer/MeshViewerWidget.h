#pragma once
#include <QString>
#include <QEvent>
#include <QMouseEvent>
#include "QGLViewerWidget.h"
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include "MeshDefinition.h"

#include "BezierCurve.h"
#include "MVC.h"

class MeshViewerWidget : public QGLViewerWidget
{
	Q_OBJECT
public:
	enum SelectMode {
		NoSelect,
		SelectAdjust,
		SelectCustom,
		Move
	};
	struct Triangle {
		std::array<OpenMesh::Vec3f, 3> vertices;
		std::array<OpenMesh::Vec2f, 3> texcoords;
		float zCenter;
	};
	MeshViewerWidget(QWidget* parent = 0);
	virtual ~MeshViewerWidget(void);
	bool LoadMesh(const std::string & filename);
	void Clear(void);
	void UpdateMesh(void);
	bool SaveMesh(const std::string & filename);
	bool ScreenShot(void);
	void SetDrawBoundingBox(bool b);
	void SetDrawBoundary(bool b);
	void EnableLighting(bool b);
	void EnableDoubleSide(bool b);
	void ResetView(void);
	void ViewCenter(void);
	void CopyRotation(void);
	void LoadRotation(void);
signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
	void SelectSMAdjust(void);
	void SetSMCustom(void);
	void SetSMMove(void);
	void SetSMUpdegree(void);
	void SetSMChangedegree(void);
	void SetSMNodrawpoint(void);
	void SetSMAddpoints(void);
	void SetSMNoSelect(void);
	void ClearSelected(void);
protected:
	virtual bool event(QEvent* _event) override;
	virtual void mouseDoubleClickEvent(QMouseEvent* _event) override;
	virtual void DrawScene(void) override;
	void DrawSceneMesh(void);
	bool NearestVertex(OpenMesh::Vec3d objCor, OpenMesh::VertexHandle& minVh);
	bool NearestCageVertex(OpenMesh::Vec3d objCor, OpenMesh::VertexHandle& minVh);
	void Bezier1Poly1(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage2poly);
	void Bezier2Poly2(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage2poly);
	void Bezier2Poly3(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage2poly);
	void Bezier2Bezier7(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage7);
	void Bezier2Bezier227(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage3);
	void Bezier2Bezier223(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage3);
	void Bezier2Bezier221(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage7);
	void Bezier2Bezier321(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage7);
	void Bezier2Bezier322(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage7);
	void Bezier2Poly7(std::vector<std::vector<Mesh::Point>> curvecage7, std::vector<std::vector<Mesh::Point>>& curvecage72poly);
	void calculate_specialgreen_weight123(std::vector<int>sp_id);
	void calculate_green_weight123(void);
	void calculate_green_weight121(void);
	void calculate_green_weight222(void);
	void calculate_green_weight221(void);
	void calculate_green_weight223(void);
	void calculate_green_weight227(void);
	void calculate_green_weight321(void);
	void calculate_green_weight322(void);
	void calculate_green_weight323(void);
	void calculate_green_weight327(void);
	void ccPoints2BezierSegs(void);
	void deform_mesh_from_cc(Mesh& deformedmesh);
	void Set_Texture_coord();

private:
	void DrawPoints(void);
	void DrawWireframe(void);
	void DrawHiddenLines(void);
	void DrawFlatLines(void);
	void DrawFlat(void);
	void DrawSmooth(void);
	void DrawBezierCurve(const BezierCurve& bezierCurve, float r = 1.0f, float g = 0.0f, float b = 0.0f, float linewidth = 1.0);
	void DrawTexture();
	void DrawCageWireframe(void);
	void DrawCagePoints(void);
	void DrawCurveCage(void);
	void DrawBoundingBox(void) const;
	void DrawBoundary(void) const;
protected:
	Mesh mesh;
	Mesh CC_mesh;//控制点组成的网格，由CC_points初始化，在move过程中顶点位置发生改变
	int degree = 3;
	int todegree = 7;
	bool highdegree = false;
	bool usecvm = false;
	bool drawpoints = true;
	std::vector<std::vector<Mesh::Point>> curvecage2;//若干段bezier曲线的控制点，曲线按逆时针顺序排序
	std::vector<Mesh::Point > CC_points;//初始化：将curvecage2中的控制点逆时针排序进一个数组内，改变：costume在move时
	std::vector<Mesh::Point> mvcGn;
	std::vector<Mesh::Point> mvcGt;
	std::vector<double> mvcL;
	std::vector<std::vector<double>> weights;
	QString strMeshFileName;
	QString strMeshBaseName;
	QString strMeshPath;
	Mesh::Point ptMin;
	Mesh::Point ptMax;
	double avgEdgeLength;
	SelectMode selectMode;
	bool isMovable;
	double moveDepth;
	OpenMesh::Vec3d lastObjCor;
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool isDrawBoundingBox;
	bool isDrawBoundary;

public slots:
	void CubicMVC_Test(void);
	void CurveCage_Test(void);
	void PolyGC_Test(void);
	void CalculateWeight_Test(void);

protected:
	bool mode_drawCurveCage = false;
};
