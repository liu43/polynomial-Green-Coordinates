#include <QtCore>
#include <complex.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshViewerWidget.h"
#include <math.h>
#include <cmath>
#define _USE_MATH_DEFINES

void DrawPoints3d(const vector<OpenMesh::Vec3d>& Points, float r = 0.0f, float g = 0.0f, float b = 1.0f, float pointsize = 5.0f);
Mesh createMeshFromCurveCage(const std::vector<Mesh::Point>& curvecage2);
std::vector<std::vector<Mesh::Point>> CCpoints_fromCCmesh(const Mesh CC_mesh,int degree);
std::vector<double> mvc(const Mesh::Point p, const std::vector<Mesh::Point> vts);
//计算积分，分子t的m次方,分母是t的1次方
double F1_n(const Mesh::Point eta, const Mesh::Point c0, const Mesh::Point c1,  int m);
//计算积分，分子t的m次方,分母是t的2次方
double F2_n(const Mesh::Point eta, const Mesh::Point c0, const Mesh::Point c1, const Mesh::Point c2, int m);
void Cardano(complex<double> a, complex<double> b, complex<double> c, complex<double>d, complex<double>& x1, complex<double>& x2, complex<double>& x3);
//计算积分，分子t的m次方，分母是t的3次方
double F3_n(const Mesh::Point eta, const Mesh::Point c0, const Mesh::Point c1, const Mesh::Point c2, const Mesh::Point c3, int m);
Mesh::Point evaluate_coor(const std::vector<double> weight, const std::vector<Mesh::Point> ctps);
//lsb关灯
MeshViewerWidget::MeshViewerWidget(QWidget* parent)
	: QGLViewerWidget(parent),
	ptMin(0.0),
	ptMax(0.0),
	isEnableLighting(false),
	isTwoSideLighting(false),
	isDrawBoundingBox(false),
	isDrawBoundary(false)
{
}

MeshViewerWidget::~MeshViewerWidget(void)
{
}

bool MeshViewerWidget::LoadMesh(const std::string & filename)
{
	Clear();
	bool read_OK = MeshTools::ReadMesh(mesh, filename);
	std::cout << "Load mesh from file " << filename << std::endl;
	if (read_OK)
	{
		strMeshFileName = QString::fromStdString(filename);
		QFileInfo fi(strMeshFileName);
		strMeshPath = fi.path();
		strMeshBaseName = fi.baseName();
		auto vertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(mesh, "vertexState");
		for (auto vh : mesh.vertices())
		{
			vertexState[vh] = NotSelected;
		}
		selectMode = NoSelect;
		UpdateMesh();
		update();
		return true;
	}
	return false;
}

void MeshViewerWidget::Clear(void)
{
	mesh.clear();
}

void MeshViewerWidget::UpdateMesh(void)
{
	mesh.update_normals();
	if (mesh.vertices_empty())
	{
		std::cerr << "ERROR: UpdateMesh() No vertices!" << std::endl;
		return;
	}
	ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
	ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;
	for (const auto& vh : mesh.vertices())
	{
		ptMin.minimize(mesh.point(vh));
		ptMax.maximize(mesh.point(vh));
	}

	double avelen = 0.0;
	double maxlen = 0.0;
	double minlen = DBL_MAX;
	for (const auto& eh : mesh.edges())
	{
		double len = mesh.calc_edge_length(eh);
		maxlen = len > maxlen ? len : maxlen;
		minlen = len < minlen ? len : minlen;
		avelen += len;
	}
	avgEdgeLength = avelen / mesh.n_vertices();

	SetScenePosition((ptMin + ptMax)*0.5, (ptMin - ptMax).norm()*0.5);
	std::cout << "Information of the input mesh:" << std::endl;
	std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
	std::cout << "  Edge Length: [" << minlen << ", " << maxlen << "]; AVG: " << avelen / mesh.n_edges() << std::endl;
}

bool MeshViewerWidget::SaveMesh(const std::string & filename)
{
	return MeshTools::WriteMesh(mesh, filename, DBL_DECIMAL_DIG);
}

bool MeshViewerWidget::ScreenShot()
{
	update();
	QString filename = strMeshPath + "/" + QDateTime::currentDateTime().toString("yyyyMMddHHmmsszzz") + QString(".png");
	QImage image = grabFramebuffer();
	image.save(filename);
	std::cout << "Save screen shot to " << filename.toStdString() << std::endl;
	return true;
}

void MeshViewerWidget::SetDrawBoundingBox(bool b)
{
	isDrawBoundingBox = b;
	update();
}
void MeshViewerWidget::SetDrawBoundary(bool b)
{
	isDrawBoundary = b;
	update();
}
void MeshViewerWidget::EnableLighting(bool b)
{
	isEnableLighting = b;
	update();
}
void MeshViewerWidget::EnableDoubleSide(bool b)
{
	isTwoSideLighting = b;
	update();
}

void MeshViewerWidget::ResetView(void)
{
	ResetModelviewMatrix();
	ViewCenter();
	update();
}

void MeshViewerWidget::ViewCenter(void)
{
	if (!mesh.vertices_empty())
	{
		UpdateMesh();
	}
	update();
}

void MeshViewerWidget::CopyRotation(void)
{
	CopyModelViewMatrix();
}

void MeshViewerWidget::LoadRotation(void)
{
	LoadCopyModelViewMatrix();
	update();
}

void MeshViewerWidget::PrintMeshInfo(void)
{
	std::cout << "Mesh Info:\n";
	std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;

	if (degree == 2)
		calculate_green_weight222();
	else if (degree == 3)
		calculate_green_weight323();
}

void MeshViewerWidget::SelectSMAdjust(void)
{
	selectMode = SelectAdjust;
}

void MeshViewerWidget::SetSMCustom(void)
{
	selectMode = SelectCustom;
}

void MeshViewerWidget::SetSMMove(void)
{
	selectMode = Move;
}
//升阶
void MeshViewerWidget::SetSMUpdegree(void)
{
	std::cout << "升阶!" << std::endl;
	calculate_green_weight327();//先计算3次cage到7次控制点的权重
	auto curvecage7 = curvecage2;
	Bezier2Bezier7(curvecage2, curvecage7);
	curvecage2 = curvecage7;
	//通过curvecage2计算cc_points;
	CC_points.clear();
	for (int i = 0; i < curvecage2.size(); i++) {
		for (int j = 0; j < todegree; j++) {
			CC_points.push_back(curvecage2[i][j]);
		}
	}
	CC_mesh = createMeshFromCurveCage(CC_points);
	auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
	for (auto vh : CC_mesh.vertices())
	{
		CagevertexState[vh] = NotSelected;
	}
	highdegree = true;
	update();
}

//变阶
void MeshViewerWidget::SetSMChangedegree(void)
{
	todegree = 2;
	highdegree = true;
	std::cout << "变阶到" << todegree << std::endl;
	auto curvecage1 = curvecage2;
	if (degree == 2) {
		if (todegree == 1) {
			calculate_green_weight221();//先计算2次cage到1次控制点的权重
			Bezier2Bezier221(curvecage2, curvecage1);
			curvecage2 = curvecage1;
		}
		if (todegree == 3) {
			calculate_green_weight223();//先计算2次cage到3次控制点的权重
			Bezier2Bezier223(curvecage2, curvecage1);
			curvecage2 = curvecage1;
		}
		if (todegree == 7) {
			calculate_green_weight227();//先计算2次cage到7次控制点的权重
			Bezier2Bezier227(curvecage2, curvecage1);
			curvecage2 = curvecage1;
		}
	}
	if (degree == 3) {
		if (todegree == 1) {
			calculate_green_weight321();//先计算2次cage到1次控制点的权重
			{
				CC_points = { OpenMesh::Vec3d(0.060138, -2.34968, 0),
	OpenMesh::Vec3d(0.167533, -1.16113, 0),
	OpenMesh::Vec3d(0.304982, 0.160104, 0),
	OpenMesh::Vec3d(0.231798, 1.08177, 0),
	OpenMesh::Vec3d(0.61252, 2.32314, 0),
	OpenMesh::Vec3d(0.771864, 3.13087, 0),
	OpenMesh::Vec3d(0.855257, 4.5897, 0),
	OpenMesh::Vec3d(3.38975, 5.46352, 0),
	OpenMesh::Vec3d(4.41819, 7.10405, 0),
	OpenMesh::Vec3d(8.67744, 7.8645, 0),
	OpenMesh::Vec3d(4.59444, 11.1053, 0),
	OpenMesh::Vec3d(-2.40215, 10.8965, 0),
	OpenMesh::Vec3d(-7.26403, 7.93585, 0),
	OpenMesh::Vec3d(-4.30306, 6.87478, 0),
	OpenMesh::Vec3d(-2.76766, 5.30325, 0),
	OpenMesh::Vec3d(-0.552003, 4.00657, 0),
	OpenMesh::Vec3d(-0.0860061, 2.21267, 0),
	OpenMesh::Vec3d(0.0878497, 0.616274, 0),
	OpenMesh::Vec3d(-0.0636332, -0.969635, 0),
	OpenMesh::Vec3d(-1.42388, 1.41961, 0),
	OpenMesh::Vec3d(-3.71789, 1.96421, 0),
	OpenMesh::Vec3d(-7.39052, 1.52015, 0),
	OpenMesh::Vec3d(-5.55652, -0.704119, 0),
	OpenMesh::Vec3d(-3.16591, -1.64514, 0),
	OpenMesh::Vec3d(-0.29709, -1.64989, 0),
	OpenMesh::Vec3d(-0.284548, -4.19343, 0),
	OpenMesh::Vec3d(-4.56897, -5.40077, 0),
	OpenMesh::Vec3d(-3.58096, -9.71999, 0),
	OpenMesh::Vec3d(-3.51223, -4.61304, 0),
	OpenMesh::Vec3d(-0.378727, -4.82229, 0),
	OpenMesh::Vec3d(0.125836, -2.84884, 0),
	OpenMesh::Vec3d(1.10361, -2.36519, 0),
	OpenMesh::Vec3d(3.49, -2.80699, 0),
	OpenMesh::Vec3d(7.08069, -3.99408, 0),
	OpenMesh::Vec3d(4.17284, -0.459165, 0),
	OpenMesh::Vec3d(0.512756, -0.0850323, 0),
				};
			}
			CC_mesh = createMeshFromCurveCage(CC_points);
			curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
			Bezier2Bezier321(curvecage2, curvecage1);
			curvecage2 = curvecage1;
		}
		if (todegree == 2) {
			calculate_green_weight322();//先计算2次cage到1次控制点的权重
			{
				CC_points = { OpenMesh::Vec3d(0.060138, -2.34968, 0),
	OpenMesh::Vec3d(0.167533, -1.16113, 0),
	OpenMesh::Vec3d(0.304982, 0.160104, 0),
	OpenMesh::Vec3d(0.231798, 1.08177, 0),
	OpenMesh::Vec3d(0.61252, 2.32314, 0),
	OpenMesh::Vec3d(0.771864, 3.13087, 0),
	OpenMesh::Vec3d(0.855257, 4.5897, 0),
	OpenMesh::Vec3d(3.38975, 5.46352, 0),
	OpenMesh::Vec3d(4.41819, 7.10405, 0),
	OpenMesh::Vec3d(8.67744, 7.8645, 0),
	OpenMesh::Vec3d(4.59444, 11.1053, 0),
	OpenMesh::Vec3d(-2.40215, 10.8965, 0),
	OpenMesh::Vec3d(-7.26403, 7.93585, 0),
	OpenMesh::Vec3d(-4.30306, 6.87478, 0),
	OpenMesh::Vec3d(-2.76766, 5.30325, 0),
	OpenMesh::Vec3d(-0.552003, 4.00657, 0),
	OpenMesh::Vec3d(-0.0860061, 2.21267, 0),
	OpenMesh::Vec3d(0.0878497, 0.616274, 0),
	OpenMesh::Vec3d(-0.0636332, -0.969635, 0),
	OpenMesh::Vec3d(-1.42388, 1.41961, 0),
	OpenMesh::Vec3d(-3.71789, 1.96421, 0),
	OpenMesh::Vec3d(-7.39052, 1.52015, 0),
	OpenMesh::Vec3d(-5.55652, -0.704119, 0),
	OpenMesh::Vec3d(-3.16591, -1.64514, 0),
	OpenMesh::Vec3d(-0.29709, -1.64989, 0),
	OpenMesh::Vec3d(-0.284548, -4.19343, 0),
	OpenMesh::Vec3d(-4.56897, -5.40077, 0),
	OpenMesh::Vec3d(-3.58096, -9.71999, 0),
	OpenMesh::Vec3d(-3.51223, -4.61304, 0),
	OpenMesh::Vec3d(-0.378727, -4.82229, 0),
	OpenMesh::Vec3d(0.125836, -2.84884, 0),
	OpenMesh::Vec3d(1.10361, -2.36519, 0),
	OpenMesh::Vec3d(3.49, -2.80699, 0),
	OpenMesh::Vec3d(7.08069, -3.99408, 0),
	OpenMesh::Vec3d(4.17284, -0.459165, 0),
	OpenMesh::Vec3d(0.512756, -0.0850323, 0),
				};
			}
			CC_mesh = createMeshFromCurveCage(CC_points);
			curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
			Bezier2Bezier322(curvecage2, curvecage1);
			curvecage2 = curvecage1;
		}
		if (todegree == 7) {
			calculate_green_weight327();//先计算2次cage到1次控制点的权重
			{
				CC_points = { OpenMesh::Vec3d(0.060138, -2.34968, 0),
	OpenMesh::Vec3d(0.167533, -1.16113, 0),
	OpenMesh::Vec3d(0.304982, 0.160104, 0),
	OpenMesh::Vec3d(0.231798, 1.08177, 0),
	OpenMesh::Vec3d(0.61252, 2.32314, 0),
	OpenMesh::Vec3d(0.771864, 3.13087, 0),
	OpenMesh::Vec3d(0.855257, 4.5897, 0),
	OpenMesh::Vec3d(3.38975, 5.46352, 0),
	OpenMesh::Vec3d(4.41819, 7.10405, 0),
	OpenMesh::Vec3d(8.67744, 7.8645, 0),
	OpenMesh::Vec3d(4.59444, 11.1053, 0),
	OpenMesh::Vec3d(-2.40215, 10.8965, 0),
	OpenMesh::Vec3d(-7.26403, 7.93585, 0),
	OpenMesh::Vec3d(-4.30306, 6.87478, 0),
	OpenMesh::Vec3d(-2.76766, 5.30325, 0),
	OpenMesh::Vec3d(-0.552003, 4.00657, 0),
	OpenMesh::Vec3d(-0.0860061, 2.21267, 0),
	OpenMesh::Vec3d(0.0878497, 0.616274, 0),
	OpenMesh::Vec3d(-0.0636332, -0.969635, 0),
	OpenMesh::Vec3d(-1.42388, 1.41961, 0),
	OpenMesh::Vec3d(-3.71789, 1.96421, 0),
	OpenMesh::Vec3d(-7.39052, 1.52015, 0),
	OpenMesh::Vec3d(-5.55652, -0.704119, 0),
	OpenMesh::Vec3d(-3.16591, -1.64514, 0),
	OpenMesh::Vec3d(-0.29709, -1.64989, 0),
	OpenMesh::Vec3d(-0.284548, -4.19343, 0),
	OpenMesh::Vec3d(-4.56897, -5.40077, 0),
	OpenMesh::Vec3d(-3.58096, -9.71999, 0),
	OpenMesh::Vec3d(-3.51223, -4.61304, 0),
	OpenMesh::Vec3d(-0.378727, -4.82229, 0),
	OpenMesh::Vec3d(0.125836, -2.84884, 0),
	OpenMesh::Vec3d(1.10361, -2.36519, 0),
	OpenMesh::Vec3d(3.49, -2.80699, 0),
	OpenMesh::Vec3d(7.08069, -3.99408, 0),
	OpenMesh::Vec3d(4.17284, -0.459165, 0),
	OpenMesh::Vec3d(0.512756, -0.0850323, 0),
				};
			}
			CC_mesh = createMeshFromCurveCage(CC_points);
			curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
			Bezier2Bezier7(curvecage2, curvecage1);
			curvecage2 = curvecage1;
		}
	}
	//通过curvecage2计算cc_points;
	CC_points.clear();
	for (int i = 0; i < curvecage2.size(); i++) {
		for (int j = 0; j < todegree; j++) {
			CC_points.push_back(curvecage2[i][j]);
		}
	}
	CC_mesh = createMeshFromCurveCage(CC_points);
	auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
	for (auto vh : CC_mesh.vertices())
	{
		CagevertexState[vh] = NotSelected;
	}
	//highdegree = true;
	update();
}

void MeshViewerWidget::SetSMNodrawpoint(void)
{
	std::cout << "no select!" << std::endl;
	if (drawpoints)
		drawpoints = false;
	else
		drawpoints = true;
	update();
}

void MeshViewerWidget::SetSMAddpoints(void)
{
	std::cout << "add points!" << std::endl;
	auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
	int vertex_index;
	for (auto vh : CC_mesh.vertices())
	{
		if (CagevertexState[vh] == Custom) {
			vertex_index = vh.idx();
		}
		
	}
	// Get the points at vertex_index and vertex_index + 1
	Mesh::Point p1 = CC_points[vertex_index];
	Mesh::Point p2 = CC_points[vertex_index + 1];

	// Calculate the three division points
	Mesh::Point p1_third = p1 + (p2 - p1) / 4.0f;
	Mesh::Point p2_third = p1 + 2.0f * (p2 - p1) / 4.0f;
	Mesh::Point p3_third = p1 + 3.0f * (p2 - p1) / 4.0f;

	// Insert the three points into the vector
	CC_points.insert(CC_points.begin() + vertex_index + 1, p3_third); // Insert third division point
	CC_points.insert(CC_points.begin() + vertex_index + 1, p2_third); // Insert second division point
	CC_points.insert(CC_points.begin() + vertex_index + 1, p1_third); // Insert first division point
	
	CC_mesh = createMeshFromCurveCage(CC_points);
	CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
	for (auto vh : CC_mesh.vertices())
	{
		CagevertexState[vh] = NotSelected;
	}

	update();
}

void MeshViewerWidget::SetSMNoSelect(void)
{
	//selectMode = NoSelect;
	std::cout << "no select!" << std::endl;
	drawpoints = false;
	update();
}

void MeshViewerWidget::ClearSelected(void)
{
	/*auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(mesh, "vertexState");
	for (auto vh : mesh.vertices())
	{
		vertexState[vh] = NotSelected;
	}*/
	Mesh deformedMesh;
	deformedMesh.assign(mesh);

	CC_points = { //deer
OpenMesh::Vec3d(8.01426, -9.42736, 0),
OpenMesh::Vec3d(8.33338, -9.22516, 0),
OpenMesh::Vec3d(8.57722, -9.0101, 0),
OpenMesh::Vec3d(8.8489, -8.78959, 0),
OpenMesh::Vec3d(7.10063, -6.55714, 0),
OpenMesh::Vec3d(6.71007, -5.29668, 0),
OpenMesh::Vec3d(6.70587, -3.83298, 0),
OpenMesh::Vec3d(5.83477, -2.3301, 0),
OpenMesh::Vec3d(7.05116, -0.928495, 0),
OpenMesh::Vec3d(6.61916, 0.630525, 0),
OpenMesh::Vec3d(8.02456, 2.09566, 0),
OpenMesh::Vec3d(8.65167, 3.70863, 0),
OpenMesh::Vec3d(6.40582, 5.04416, 0),
OpenMesh::Vec3d(6.37169, 3.87665, 0),
OpenMesh::Vec3d(5.80833, 2.98216, 0),
OpenMesh::Vec3d(5.62158, 2.09523, 0),
OpenMesh::Vec3d(3.23511, 2.54016, 0),
OpenMesh::Vec3d(0.750374, 1.38221, 0),
OpenMesh::Vec3d(-1.69483, 2.49143, 0),
OpenMesh::Vec3d(-3.65702, 4.17105, 0),
OpenMesh::Vec3d(-2.17201, 6.0254, 0),
OpenMesh::Vec3d(-3.05934, 7.31317, 0),
OpenMesh::Vec3d(-1.89982, 8.15725, 0),
OpenMesh::Vec3d(-1.63298, 9.67336, 0),
OpenMesh::Vec3d(-2.82579, 10.398, 0),
OpenMesh::Vec3d(-3.92575, 10.5965, 0),
OpenMesh::Vec3d(-2.76592, 8.26476, 0),
OpenMesh::Vec3d(-5.78046, 7.07922, 0),
OpenMesh::Vec3d(-6.48863, 6.34431, 0),
OpenMesh::Vec3d(-5.30909, 5.86309, 0),
OpenMesh::Vec3d(-7.14945, 4.31455, 0),
OpenMesh::Vec3d(-6.72589, 3.45316, 0),
OpenMesh::Vec3d(-5.35703, 4.33777, 0),
OpenMesh::Vec3d(-4.45691, 4.19383, 0),
OpenMesh::Vec3d(-3.65641, 2.89317, 0),
OpenMesh::Vec3d(-4.13472, -0.366618, 0),
OpenMesh::Vec3d(-2.87246, -1.14708, 0),
OpenMesh::Vec3d(-3.61957, -1.90331, 0),
OpenMesh::Vec3d(-2.36707, -2.13026, 0),
OpenMesh::Vec3d(-2.28264, -2.59782, 0),
OpenMesh::Vec3d(-4.87634, -4.58312, 0),
OpenMesh::Vec3d(-3.97765, -6.59246, 0),
OpenMesh::Vec3d(-1.14421, -9.59584, 0),
OpenMesh::Vec3d(-0.906403, -9.58534, 0),
OpenMesh::Vec3d(-0.639413, -9.58154, 0),
OpenMesh::Vec3d(-0.305103, -9.56834, 0),
OpenMesh::Vec3d(-2.6836, -6.57468, 0),
OpenMesh::Vec3d(-3.98979, -4.66978, 0),
OpenMesh::Vec3d(-0.925819, -3.20408, 0),
OpenMesh::Vec3d(0.770266, -1.45608, 0),
OpenMesh::Vec3d(2.14696, -1.64475, 0),
OpenMesh::Vec3d(3.40852, -2.90979, 0),
OpenMesh::Vec3d(4.44002, -4.23804, 0),
OpenMesh::Vec3d(5.05317, -4.38792, 0),
OpenMesh::Vec3d(6.7502, -6.72501, 0),
OpenMesh::Vec3d(7.59131, -7.41388, 0),
OpenMesh::Vec3d(7.20539, -8.66758, 0),
	};


//	CC_points = { //ball 变形
//		OpenMesh::Vec3d(4.93469, -4.95787, 0),
//OpenMesh::Vec3d(4.00952, -2.92598, 0),
//OpenMesh::Vec3d(3.73102, 1.06033, 0),
//OpenMesh::Vec3d(4.89974, 5.07107, 0),
//OpenMesh::Vec3d(1.80626, 2.01981, 0),
//OpenMesh::Vec3d(-2.8923, 2.12251, 0),
//OpenMesh::Vec3d(-5.78844, 4.26243, 0),
//OpenMesh::Vec3d(-4.26235, 1.69334, 0),
//OpenMesh::Vec3d(-4.06703, -2.70747, 0),
//OpenMesh::Vec3d(-4.89242, -5.17084, 0),
//OpenMesh::Vec3d(-2.12912, -3.0845, 0),
//OpenMesh::Vec3d(1.4489, -2.94237, 0),
//	};
//	CC_points = {//changjinlu
//OpenMesh::Vec3d(4.26544, -2.08319, 0),
//OpenMesh::Vec3d(2.14808, -1.76166, 0),
//OpenMesh::Vec3d(1.05683, 0.58523, 0),
//OpenMesh::Vec3d(-0.959219, 1.60384, 0),
//OpenMesh::Vec3d(-2.12218, 6.64445, 0),
//OpenMesh::Vec3d(-7.00934, 10.7648, 0),
//OpenMesh::Vec3d(-3.68572, 13.2197, 0),
//OpenMesh::Vec3d(-5.01245, 15.6318, 0),
//OpenMesh::Vec3d(-7.78227, 15.2923, 0),
//OpenMesh::Vec3d(-9.72932, 12.9165, 0),
//OpenMesh::Vec3d(-8.66079, 7.13327, 0),
//OpenMesh::Vec3d(-4.31018, 5.64754, 0),
//OpenMesh::Vec3d(-4.30619, 1.46692, 0),
//OpenMesh::Vec3d(-5.6163, -8.3802, 0),
//OpenMesh::Vec3d(-7.78844, -18.0231, 0),
//OpenMesh::Vec3d(2.65463, -13.4286, 0),
//OpenMesh::Vec3d(0.741755, -4.48237, 0),
//OpenMesh::Vec3d(2.73214, -4.76554, 0),
//OpenMesh::Vec3d(3.64475, -3.15402, 0),
//OpenMesh::Vec3d(5.22527, -2.66165, 0),
//OpenMesh::Vec3d(6.60758, -3.9218, 0),
//OpenMesh::Vec3d(6.79063, -7.53675, 0),
//OpenMesh::Vec3d(8.52689, -0.989955, 0),
//OpenMesh::Vec3d(5.22416, -1.9387, 0),
//	};

//	CC_points = {//paris tower
//		OpenMesh::Vec3d(2.21312, -3.64502, 0),
//OpenMesh::Vec3d(0.0729413, -1.26159, 0),
//OpenMesh::Vec3d(0.71968, 2.04631, 0),
//OpenMesh::Vec3d(-3.67721, 2.8085, 0),
//OpenMesh::Vec3d(-0.533977, 2.06356, 0),
//OpenMesh::Vec3d(0.020478, -1.26554, 0),
//OpenMesh::Vec3d(-2.13172, -3.62273, 0),
//OpenMesh::Vec3d(-1.5717, -3.65431, 0),
//OpenMesh::Vec3d(-0.789005, -3.66778, 0),
//OpenMesh::Vec3d(-0.017268, -3.66608, 0),
//OpenMesh::Vec3d(0.997374, -3.65183, 0),
//OpenMesh::Vec3d(1.73898, -3.63311, 0),
//	};

//	CC_points = { //xiyi
//		OpenMesh::Vec3d(4.51664, -0.509111, 0),
//OpenMesh::Vec3d(3.51089, 0.0120201, 0),
//OpenMesh::Vec3d(2.58026, 0.369127, 0),
//OpenMesh::Vec3d(0.899758, 1.92868, 0),
//OpenMesh::Vec3d(-6.36326, 3.14212, 0),
//OpenMesh::Vec3d(-6.99906, 0.160185, 0),
//OpenMesh::Vec3d(-3.73234, -0.87026, 0),
//OpenMesh::Vec3d(-4.44834, -2.77133, 0),
//OpenMesh::Vec3d(-1.61845, -2.46282, 0),
//OpenMesh::Vec3d(1.3063, -2.28422, 0),
//OpenMesh::Vec3d(4.22574, -1.92609, 0),
//OpenMesh::Vec3d(4.29148, -1.36913, 0),
//OpenMesh::Vec3d(5.32878, -1.19727, 0),
//OpenMesh::Vec3d(5.83854, -1.11686, 0),
//OpenMesh::Vec3d(6.02733, -1.16366, 0),
//OpenMesh::Vec3d(6.34854, -1.15323, 0),
//OpenMesh::Vec3d(6.61149, -1.13321, 0),
//OpenMesh::Vec3d(6.64945, -1.1902, 0),
//OpenMesh::Vec3d(6.79034, -1.02896, 0),
//OpenMesh::Vec3d(6.37798, -1.20285, 0),
//OpenMesh::Vec3d(6.18146, -0.962081, 0),
//OpenMesh::Vec3d(6.05019, -1.01034, 0),
//OpenMesh::Vec3d(5.83587, -0.977423, 0),
//OpenMesh::Vec3d(5.29081, -0.736826, 0),
//	};
//CC_points = { //yu
//	OpenMesh::Vec3d(1.92641, -1.5622, 0),
//OpenMesh::Vec3d(1.92093, -1.42981, 0),
//OpenMesh::Vec3d(1.91418, -1.27198, 0),
//OpenMesh::Vec3d(1.90741, -1.11983, 0),
//OpenMesh::Vec3d(1.90298, -0.988665, 0),
//OpenMesh::Vec3d(1.90126, -0.842423, 0),
//OpenMesh::Vec3d(1.89499, -0.723717, 0),
//OpenMesh::Vec3d(1.88534, -0.589824, 0),
//OpenMesh::Vec3d(0.918267, 2.37635, 0),
//OpenMesh::Vec3d(0.189832, -5.62026, 0),
//OpenMesh::Vec3d(-1.08374, 2.40682, 0),
//OpenMesh::Vec3d(-1.3776, 0.930637, 0),
//OpenMesh::Vec3d(-2.33054, -3.02618, 0),
//OpenMesh::Vec3d(-2.78674, 1.10829, 0),
//OpenMesh::Vec3d(-3.24116, -0.577623, 0),
//OpenMesh::Vec3d(-3.23842, -0.686008, 0),
//OpenMesh::Vec3d(-3.23358, -0.780003, 0),
//OpenMesh::Vec3d(-3.22776, -0.859145, 0),
//OpenMesh::Vec3d(-3.22375, -0.936334, 0),
//OpenMesh::Vec3d(-3.22635, -1.01901, 0),
//OpenMesh::Vec3d(-3.22532, -1.11623, 0),
//OpenMesh::Vec3d(-3.2268, -1.23156, 0),
//OpenMesh::Vec3d(-2.75638, 0.473092, 0),
//OpenMesh::Vec3d(-2.32616, -3.6897, 0),
//OpenMesh::Vec3d(-1.35225, 0.368065, 0),
//OpenMesh::Vec3d(-1.1074, 2.01805, 0),
//OpenMesh::Vec3d(0.203981, -5.85821, 0),
//OpenMesh::Vec3d(0.949495, 1.71095, 0),
//};


//CC_points = {//kuzi
//	OpenMesh::Vec3d(11.4493, 14.8504, 0),
//OpenMesh::Vec3d(6.50769, 22.2384, 0),
//OpenMesh::Vec3d(-3.7327, 14.8504, 0),
//OpenMesh::Vec3d(-11.3237, 14.8504, 0),
//OpenMesh::Vec3d(-10.855, -2.94221, 0),
//OpenMesh::Vec3d(-13.0942, -0.373684, 0),
//OpenMesh::Vec3d(-21.0616, -0.669252, 0),
//OpenMesh::Vec3d(-20.9948, -3.2209, 0),
//OpenMesh::Vec3d(-21.0322, -5.72048, 0),
//OpenMesh::Vec3d(-21.382, -7.49102, 0),
//OpenMesh::Vec3d(-14.1957, -7.37668, 0),
//OpenMesh::Vec3d(-0.656309, -9.60571, 0),
//OpenMesh::Vec3d(-0.812533, 1.34013, 0),
//OpenMesh::Vec3d(-0.232666, 1.34382, 0),
//OpenMesh::Vec3d(0.3472, 1.34751, 0),
//OpenMesh::Vec3d(0.927067, 1.3512, 0),
//OpenMesh::Vec3d(11.3977, 1.806, 0),
//OpenMesh::Vec3d(13.901, 1.74005, 0),
//OpenMesh::Vec3d(14.4254, -15.146, 0),
//OpenMesh::Vec3d(17.9292, -15.146, 0),
//OpenMesh::Vec3d(21.4329, -15.146, 0),
//OpenMesh::Vec3d(24.9366, -15.146, 0),
//OpenMesh::Vec3d(24.7283, 0.459509, 0),
//OpenMesh::Vec3d(11.4493, 4.86897, 0),
//};
 

//CC_points = {//car
//OpenMesh::Vec3d(18.207, -4.48427, 0),
//OpenMesh::Vec3d(17.0602, -0.865747, 0),
//OpenMesh::Vec3d(16.9603, 3.13988, 0),
//OpenMesh::Vec3d(13.1047, 2.52286, 0),
//OpenMesh::Vec3d(14.1094, 21.6961, 0),
//OpenMesh::Vec3d(4.01465, 16.2523, 0),
//OpenMesh::Vec3d(-2.6146, 4.44862, 0),
//OpenMesh::Vec3d(-6.64821, 6.31885, 0),
//OpenMesh::Vec3d(-12.1076, 5.34944, 0),
//OpenMesh::Vec3d(-12.7923, -2.21603, 0),
//OpenMesh::Vec3d(-5.40276, -4.51037, 0),
//OpenMesh::Vec3d(8.38458, 0.846978, 0),
//};

//CC_points = {//tuzi
//	OpenMesh::Vec3d(1.2338, -0.872146, 0),
//OpenMesh::Vec3d(0.848741, -0.488863, 0),
//OpenMesh::Vec3d(0.994741, -0.428605, 0),
//OpenMesh::Vec3d(1.03727, -0.118332, 0),
//OpenMesh::Vec3d(1.48523, -0.403536, 0),
//OpenMesh::Vec3d(0.959169, -0.893975, 0),
//OpenMesh::Vec3d(1.47613, -1.00747, 0),
//OpenMesh::Vec3d(1.97429, -0.710009, 0),
//OpenMesh::Vec3d(1.54183, 0.144788, 0),
//OpenMesh::Vec3d(0.507068, 0.378323, 0),
//OpenMesh::Vec3d(-0.335255, 0.525864, 0),
//OpenMesh::Vec3d(-0.896345, 0.022467, 0),
//OpenMesh::Vec3d(-0.979956, -0.381595, 0),
//OpenMesh::Vec3d(-1.03288, -0.89268, 0),
//OpenMesh::Vec3d(-1.25456, -0.83442, 0),
//OpenMesh::Vec3d(-1.21265, -1.61744, 0),
//OpenMesh::Vec3d(-0.690139, -1.62269, 0),
//OpenMesh::Vec3d(-0.28394, -1.64832, 0),
//OpenMesh::Vec3d(0.223731, -1.65626, 0),
//OpenMesh::Vec3d(0.416271, -1.6679, 0),
//OpenMesh::Vec3d(0.533081, -1.66271, 0),
//OpenMesh::Vec3d(0.753259, -1.66561, 0),
//OpenMesh::Vec3d(0.811291, -1.37615, 0),
//OpenMesh::Vec3d(1.19484, -1.34047, 0),
//};
 

//CC_points = { //daizi
//	OpenMesh::Vec3d(-45.4715, -4.44834, 0),
//OpenMesh::Vec3d(-45.1077, -3.19904, 0),
//OpenMesh::Vec3d(-43.3737, -2.21489, 0),
//OpenMesh::Vec3d(-44.2633, -0.503097, 0),
//OpenMesh::Vec3d(-44.2644, 0.337876, 0),
//OpenMesh::Vec3d(-41.3244, 1.5609, 0),
//OpenMesh::Vec3d(-41.3281, 2.27688, 0),
//OpenMesh::Vec3d(-41.3336, 3.49199, 0),
//OpenMesh::Vec3d(-72.1444, -2.66667, 0),
//OpenMesh::Vec3d(-77.7012, -26.7162, 0),
//OpenMesh::Vec3d(-60.2857, -59.3187, 0),
//OpenMesh::Vec3d(-22.2499, -59.0844, 0),
//OpenMesh::Vec3d(-7.31521, -26.3709, 0),
//OpenMesh::Vec3d(-16.4601, -3.0048, 0),
//OpenMesh::Vec3d(-45.0842, 3.50565, 0),
//OpenMesh::Vec3d(-45.0786, 2.25173, 0),
//OpenMesh::Vec3d(-45.0749, 1.22532, 0),
//OpenMesh::Vec3d(-42.7613, 0.304172, 0),
//OpenMesh::Vec3d(-42.7602, -0.536801, 0),
//OpenMesh::Vec3d(-43.7393, -2.27171, 0),
//OpenMesh::Vec3d(-41.8222, -3.21588, 0),
//OpenMesh::Vec3d(-41.8222, -4.45255, 0),
//OpenMesh::Vec3d(-20.0514, -7.14788, 0),
//OpenMesh::Vec3d(-14.1207, -26.0952, 0),
//OpenMesh::Vec3d(-25.508, -54.5302, 0),
//OpenMesh::Vec3d(-56.9787, -54.7188, 0),
//OpenMesh::Vec3d(-71.292, -26.2217, 0),
//OpenMesh::Vec3d(-66.6177, -7.2701, 0),
//};
// 

//CC_points = { //fish
//	OpenMesh::Vec3d(8.43404, -9.00585, 0),
//OpenMesh::Vec3d(7.47304, -1.84145, 0),
//OpenMesh::Vec3d(7.60398, -0.916129, 0),
//OpenMesh::Vec3d(10.8396, 7.86985, 0),
//OpenMesh::Vec3d(8.99175, 1.64069, 0),
//OpenMesh::Vec3d(4.69127, -0.627544, 0),
//OpenMesh::Vec3d(2.66809, -0.0643164, 0),
//OpenMesh::Vec3d(2.72136, 0.590422, 0),
//OpenMesh::Vec3d(1.37385, 1.77557, 0),
//OpenMesh::Vec3d(0.240029, 1.49932, 0),
//OpenMesh::Vec3d(-2.49334, 2.50637, 0),
//OpenMesh::Vec3d(-6.33127, 3.18449, 0),
//OpenMesh::Vec3d(-8.90656, 2.49647, 0),
//OpenMesh::Vec3d(-10.3053, 3.09553, 0),
//OpenMesh::Vec3d(-10.0761, 2.09196, 0),
//OpenMesh::Vec3d(-9.00522, 1.38062, 0),
//OpenMesh::Vec3d(-8.98836, 1.07329, 0),
//OpenMesh::Vec3d(-8.99063, 0.702774, 0),
//OpenMesh::Vec3d(-8.75359, 0.482605, 0),
//OpenMesh::Vec3d(-8.63689, 0.763358, 0),
//OpenMesh::Vec3d(-8.78799, 1.02137, 0),
//OpenMesh::Vec3d(-8.72299, 1.21815, 0),
//OpenMesh::Vec3d(-6.77407, -2.07211, 0),
//OpenMesh::Vec3d(-3.72537, -3.12714, 0),
//OpenMesh::Vec3d(-1.61215, -3.01141, 0),
//OpenMesh::Vec3d(-1.16986, -3.56759, 0),
//OpenMesh::Vec3d(0.232353, -4.09102, 0),
//OpenMesh::Vec3d(0.0920223, -2.80436, 0),
//OpenMesh::Vec3d(0.365635, -2.75175, 0),
//OpenMesh::Vec3d(0.651136, -2.71796, 0),
//OpenMesh::Vec3d(0.895707, -2.68846, 0),
//OpenMesh::Vec3d(1.40535, -1.90493, 0),
//OpenMesh::Vec3d(3.27536, -4.52581, 0),
//OpenMesh::Vec3d(2.61347, -2.47371, 0),
//OpenMesh::Vec3d(5.24685, -1.85539, 0),
//OpenMesh::Vec3d(6.49604, -3.83112, 0),
//};
 
//CC_points = { //choufish
//	OpenMesh::Vec3d(2.27477, -5.47137, 0),
//    OpenMesh::Vec3d(1.8946, -4.29797, 0),
//    OpenMesh::Vec3d(1.76649, -3.4866, 0),
//    OpenMesh::Vec3d(2.33876, -3.72963, 0),
//    OpenMesh::Vec3d(4.26167, -3.63543, 0),
//    OpenMesh::Vec3d(3.85949, -1.66655, 0),
//    OpenMesh::Vec3d(5.77752, -1.50034, 0),
//    OpenMesh::Vec3d(4.98695, 0.732251, 0),
//    OpenMesh::Vec3d(3.27073, 0.866982, 0),
//    OpenMesh::Vec3d(2.85876, 2.93215, 0),
//    OpenMesh::Vec3d(0.10979, 3.08279, 0),
//    OpenMesh::Vec3d(-1.88859, 4.21388, 0),
//    OpenMesh::Vec3d(0.0442605, 5.69605, 0),
//    OpenMesh::Vec3d(-0.172236, 5.83505, 0),
//    OpenMesh::Vec3d(-0.488855, 6.00424, 0),
//    OpenMesh::Vec3d(-0.872024, 6.13353, 0),
//    OpenMesh::Vec3d(-3.33598, 4.58416, 0),
//    OpenMesh::Vec3d(-1.63781, 2.39221, 0),
//    OpenMesh::Vec3d(0.310345, 2.27181, 0),
//    OpenMesh::Vec3d(-0.884167, 1.77491, 0),
//    OpenMesh::Vec3d(-2.55503, 1.80748, 0),
//    OpenMesh::Vec3d(-3.35906, 1.4997, 0),
//    OpenMesh::Vec3d(-5.28685, 1.27267, 0),
//    OpenMesh::Vec3d(-5.5253, -1.89937, 0),
//    OpenMesh::Vec3d(-3.81445, -2.608, 0),
//    OpenMesh::Vec3d(-2.6613, -1.81394, 0),
//    OpenMesh::Vec3d(-2.12056, -2.29342, 0),
//    OpenMesh::Vec3d(-2.20238, -3.11821, 0),
//    OpenMesh::Vec3d(-1.63304, -3.92693, 0),
//    OpenMesh::Vec3d(-0.393606, -4.68184, 0),
//    OpenMesh::Vec3d(0.503098, -2.97635, 0),
//    OpenMesh::Vec3d(0.916171, -3.35227, 0),
//    OpenMesh::Vec3d(1.32529, -3.13754, 0),
//    OpenMesh::Vec3d(1.78622, -2.87407, 0),
//    OpenMesh::Vec3d(1.21669, -3.3626, 0),
//    OpenMesh::Vec3d(1.27813, -4.24631, 0),
//};


//CC_points= { //zomm fish
//	OpenMesh::Vec3d(10.5207, 1.26963, 0),
//OpenMesh::Vec3d(9.80995, 20.178, 0),
//OpenMesh::Vec3d(-10.8162, 20.0895, 0),
//OpenMesh::Vec3d(-11.1068, 1.26963, 0),
//OpenMesh::Vec3d(-11.1068, 0.22077, 0),
//OpenMesh::Vec3d(-11.1068, -0.82809, 0),
//OpenMesh::Vec3d(-11.1068, -1.87695, 0),
//OpenMesh::Vec3d(-7.5541, 10.7143, 0),
//OpenMesh::Vec3d(7.14498, 10.7733, 0),
//OpenMesh::Vec3d(10.5207, -1.87695, 0),
//OpenMesh::Vec3d(10.5207, -0.82809, 0),
//OpenMesh::Vec3d(10.5207, 0.22077, 0),
//};


//CC_points = { //deformed circle
//	OpenMesh::Vec3d(3.45289, 1.63941, 0),
//OpenMesh::Vec3d(4.30953, -0.696405, 0),
//OpenMesh::Vec3d(8.83846, -0.334435, 0),
//OpenMesh::Vec3d(7.30388, 2.62878, 0),
//OpenMesh::Vec3d(6.28855, -2.62518, 0),
//OpenMesh::Vec3d(3.90656, 6.20451, 0),
//OpenMesh::Vec3d(-1.69168, 3.85054, 0),
//OpenMesh::Vec3d(1.09155, 5.00703, 0),
//OpenMesh::Vec3d(1.30912, 9.30382, 0),
//OpenMesh::Vec3d(-2.54111, 7.34319, 0),
//OpenMesh::Vec3d(3.17407, 6.91966, 0),
//OpenMesh::Vec3d(-5.90614, 4.59064, 0),
//OpenMesh::Vec3d(-3.65812, -1.22062, 0),
//OpenMesh::Vec3d(-4.82526, 1.79114, 0),
//OpenMesh::Vec3d(-8.67249, 0.318421, 0),
//OpenMesh::Vec3d(-7.36877, -2.80688, 0),
//OpenMesh::Vec3d(-6.62193, 2.72038, 0),
//OpenMesh::Vec3d(-4.64774, -5.87751, 0),
//OpenMesh::Vec3d(1.34646, -3.2845, 0),
//OpenMesh::Vec3d(-1.61716, -4.57378, 0),
//OpenMesh::Vec3d(-1.13283, -8.95429, 0),
//OpenMesh::Vec3d(2.75132, -7.28116, 0),
//OpenMesh::Vec3d(-3.34859, -6.70248, 0),
//OpenMesh::Vec3d(5.46868, -4.0388, 0),
//};

//CC_points = { //deformed circle new
//	OpenMesh::Vec3d(3.45289, 1.63941, 0),
//OpenMesh::Vec3d(6.76375, 1.57117, 0),
//OpenMesh::Vec3d(6.50663, 6.68737, 0),
//OpenMesh::Vec3d(2.59476, 7.52767, 0),
//OpenMesh::Vec3d(6.52016, 3.41986, 0),
//OpenMesh::Vec3d(3.51282, 5.94974, 0),
//OpenMesh::Vec3d(-1.69168, 3.85054, 0),
//OpenMesh::Vec3d(-0.739426, 7.46023, 0),
//OpenMesh::Vec3d(-5.80258, 7.03216, 0),
//OpenMesh::Vec3d(-6.85908, 3.35323, 0),
//OpenMesh::Vec3d(-3.01995, 6.31121, 0),
//OpenMesh::Vec3d(-5.5124, 3.82632, 0),
//OpenMesh::Vec3d(-3.65812, -1.22062, 0),
//OpenMesh::Vec3d(-7.57777, -1.29563, 0),
//OpenMesh::Vec3d(-6.16779, -6.64343, 0),
//OpenMesh::Vec3d(-2.63291, -6.80906, 0),
//OpenMesh::Vec3d(-6.03714, -3.68413, 0),
//OpenMesh::Vec3d(-3.69563, -5.43646, 0),
//OpenMesh::Vec3d(1.34646, -3.2845, 0),
//OpenMesh::Vec3d(1.52332, -6.63453, 0),
//OpenMesh::Vec3d(5.49032, -7.07563, 0),
//OpenMesh::Vec3d(6.97355, -2.60171, 0),
//OpenMesh::Vec3d(2.96145, -6.57085, 0),
//OpenMesh::Vec3d(5.27339, -3.66767, 0),
//};


//hunman
//CC_points = { //human 2次
//	OpenMesh::Vec3d(2.2434, -4.09831, 0),
//OpenMesh::Vec3d(1.33897, -1.2763, 0),
//OpenMesh::Vec3d(0.168209, 0.44022, 0),
//OpenMesh::Vec3d(1.4422, 0.351278, 0),
//OpenMesh::Vec3d(2.71398, 0.333548, 0),
//OpenMesh::Vec3d(2.71397, 1.68434, 0),
//OpenMesh::Vec3d(3.46302, 4.87632, 0),
//OpenMesh::Vec3d(2.02594, 5.163, 0),
//OpenMesh::Vec3d(0.502351, 5.27545, 0),
//OpenMesh::Vec3d(0.555983, 4.54789, 0),
//OpenMesh::Vec3d(0.618441, 3.76385, 0),
//OpenMesh::Vec3d(1.34768, 3.70515, 0),
//OpenMesh::Vec3d(2.12808, 3.63751, 0),
//OpenMesh::Vec3d(1.15359, 2.4928, 0),
//OpenMesh::Vec3d(0.929425, 1.2952, 0),
//OpenMesh::Vec3d(0.312689, 3.28967, 0),
//OpenMesh::Vec3d(-0.901575, 4.85807, 0),
//OpenMesh::Vec3d(-1.23977, 3.5368, 0),
//OpenMesh::Vec3d(-1.9027, 4.60764, 0),
//OpenMesh::Vec3d(-2.60159, 4.08557, 0),
//OpenMesh::Vec3d(-3.67377, 2.94282, 0),
//OpenMesh::Vec3d(-4.19734, 1.3431, 0),
//OpenMesh::Vec3d(-2.67894, 0.569752, 0),
//OpenMesh::Vec3d(-4.41661, 0.875907, 0),
//OpenMesh::Vec3d(-5.35647, 1.84359, 0),
//OpenMesh::Vec3d(-4.28499, -1.33886, 0),
//OpenMesh::Vec3d(-2.25433, -1.02231, 0),
//OpenMesh::Vec3d(-2.9142, -1.57256, 0),
//OpenMesh::Vec3d(-3.11938, -1.99522, 0),
//OpenMesh::Vec3d(-2.32945, -2.77795, 0),
//OpenMesh::Vec3d(-2.39632, -3.96978, 0),
//OpenMesh::Vec3d(-0.151571, -4.57271, 0),
//};

//CC_points = { //human 1次
//	OpenMesh::Vec3d(2.2434, -4.09831, 0),
//OpenMesh::Vec3d(0.168209, 0.44022, 0),
//OpenMesh::Vec3d(2.59141, 0.159537, 0),
//OpenMesh::Vec3d(2.75152, 5.52769, 0),
//OpenMesh::Vec3d(0.311402, 5.24667, 0),
//OpenMesh::Vec3d(0.0708214, 3.50104, 0),
//OpenMesh::Vec3d(2.1565, 3.58798, 0),
//OpenMesh::Vec3d(0.929425, 1.2952, 0),
//OpenMesh::Vec3d(-0.901575, 4.85807, 0),
//OpenMesh::Vec3d(-1.9027, 4.60764, 0),
//OpenMesh::Vec3d(-3.67377, 2.94282, 0),
//OpenMesh::Vec3d(-2.67894, 0.569752, 0),
//OpenMesh::Vec3d(-4.13832, -0.216033, 0),
//OpenMesh::Vec3d(-2.25433, -1.02231, 0),
//OpenMesh::Vec3d(-3.11938, -1.99522, 0),
//OpenMesh::Vec3d(-2.39632, -3.96978, 0),
//};
//CC_points = { //三次
//	OpenMesh::Vec3d(2.2434, -4.09831, 0),
//OpenMesh::Vec3d(-0.769148, -3.04015, 0),
//OpenMesh::Vec3d(2.3222, -0.889475, 0),
//OpenMesh::Vec3d(0.168209, 0.44022, 0),
//OpenMesh::Vec3d(1.009, -0.045892, 0),
//OpenMesh::Vec3d(1.84906, -0.38449, 0),
//OpenMesh::Vec3d(2.68837, -0.575573, 0),
//OpenMesh::Vec3d(1.55687, 1.20745, 0),
//OpenMesh::Vec3d(4.88948, 3.42979, 0),
//OpenMesh::Vec3d(3.94164, 4.94133, 0),
//OpenMesh::Vec3d(3.07853, 4.98889, 0),
//OpenMesh::Vec3d(2.06426, 5.05415, 0),
//OpenMesh::Vec3d(0.591506, 5.43161, 0),
//OpenMesh::Vec3d(0.584249, 4.90679, 0),
//OpenMesh::Vec3d(0.512881, 4.38422, 0),
//OpenMesh::Vec3d(0.210941, 3.88952, 0),
//OpenMesh::Vec3d(0.623509, 3.67733, 0),
//OpenMesh::Vec3d(1.0465, 3.47584, 0),
//OpenMesh::Vec3d(1.45429, 3.28504, 0),
//OpenMesh::Vec3d(2.04615, 2.5743, 0),
//OpenMesh::Vec3d(0.742298, 2.41874, 0),
//OpenMesh::Vec3d(0.929425, 1.2952, 0),
//OpenMesh::Vec3d(0.518268, 2.62485, 0),
//OpenMesh::Vec3d(-0.0920657, 3.81247, 0),
//OpenMesh::Vec3d(-0.901575, 4.85807, 0),
//OpenMesh::Vec3d(-1.12704, 3.97722, 0),
//OpenMesh::Vec3d(-1.46075, 3.89375, 0),
//OpenMesh::Vec3d(-1.9027, 4.60764, 0),
//OpenMesh::Vec3d(-2.36863, 4.25959, 0),
//OpenMesh::Vec3d(-2.95898, 3.70465, 0),
//OpenMesh::Vec3d(-3.67377, 2.94282, 0),
//OpenMesh::Vec3d(-4.02282, 1.87634, 0),
//OpenMesh::Vec3d(-3.69121, 1.08532, 0),
//OpenMesh::Vec3d(-2.67894, 0.569752, 0),
//OpenMesh::Vec3d(-4.13322, 1.15321, 0),
//OpenMesh::Vec3d(-4.05724, -1.83736, 0),
//OpenMesh::Vec3d(-5.46007, 0.228285, 0),
//OpenMesh::Vec3d(-3.21245, -2.86571, 0),
//OpenMesh::Vec3d(-3.14689, -0.405828, 0),
//OpenMesh::Vec3d(-2.25433, -1.02231, 0),
//OpenMesh::Vec3d(-2.69424, -1.38914, 0),
//OpenMesh::Vec3d(-2.98259, -1.71345, 0),
//OpenMesh::Vec3d(-3.11938, -1.99522, 0),
//OpenMesh::Vec3d(-2.59276, -2.51704, 0),
//OpenMesh::Vec3d(-2.35174, -3.17523, 0),
//OpenMesh::Vec3d(-2.39632, -3.96978, 0),
//OpenMesh::Vec3d(-0.899821, -4.37173, 0),
//OpenMesh::Vec3d(0.646753, -4.41458, 0),
//};

//CC_points= { //human 7次
//	OpenMesh::Vec3d(2.2434, -4.09831, 0),
//OpenMesh::Vec3d(1.98499, -3.29202, 0),
//OpenMesh::Vec3d(1.7139, -2.53838, 0),
//OpenMesh::Vec3d(1.43013, -1.83737, 0),
//OpenMesh::Vec3d(1.13367, -1.18901, 0),
//OpenMesh::Vec3d(0.824533, -0.593291, 0),
//OpenMesh::Vec3d(0.502712, -0.0502143, 0),
//OpenMesh::Vec3d(0.168209, 0.44022, 0),
//OpenMesh::Vec3d(0.528549, 0.231886, 0),
//OpenMesh::Vec3d(0.888784, 0.044626, 0),
//OpenMesh::Vec3d(1.24891, -0.121561, 0),
//OpenMesh::Vec3d(1.60894, -0.266674, 0),
//OpenMesh::Vec3d(1.96885, -0.390714, 0),
//OpenMesh::Vec3d(2.32866, -0.49368, 0),
//OpenMesh::Vec3d(2.68837, -0.575573, 0),
//OpenMesh::Vec3d(2.36958, 0.574717, 0),
//OpenMesh::Vec3d(4.04595, -0.611701, 0),
//OpenMesh::Vec3d(3.37615, 1.16304, 0),
//OpenMesh::Vec3d(4.94639, 0.0309733, 0),
//OpenMesh::Vec3d(4.59013, 1.75041, 0),
//OpenMesh::Vec3d(5.29918, 1.4732, 0),
//OpenMesh::Vec3d(5.28407, 2.03122, 0),
//OpenMesh::Vec3d(5.2983, 2.51805, 0),
//OpenMesh::Vec3d(5.30008, 3.00193, 0),
//OpenMesh::Vec3d(5.28941, 3.48284, 0),
//OpenMesh::Vec3d(5.2663, 3.9608, 0),
//OpenMesh::Vec3d(5.23073, 4.4358, 0),
//OpenMesh::Vec3d(5.18272, 4.90783, 0),
//OpenMesh::Vec3d(5.12226, 5.37691, 0),
//OpenMesh::Vec3d(4.68562, 5.47027, 0),
//OpenMesh::Vec3d(4.2691, 5.53468, 0),
//OpenMesh::Vec3d(3.87268, 5.57015, 0),
//OpenMesh::Vec3d(3.49637, 5.57668, 0),
//OpenMesh::Vec3d(3.14017, 5.55426, 0),
//OpenMesh::Vec3d(2.80408, 5.50289, 0),
//OpenMesh::Vec3d(2.4881, 5.42258, 0),
//OpenMesh::Vec3d(2.41248, 5.1286, 0),
//OpenMesh::Vec3d(2.35664, 4.84529, 0),
//OpenMesh::Vec3d(2.32059, 4.57266, 0),
//OpenMesh::Vec3d(2.30431, 4.3107, 0),
//OpenMesh::Vec3d(2.30782, 4.05942, 0),
//OpenMesh::Vec3d(2.3311, 3.81881, 0),
//OpenMesh::Vec3d(2.37417, 3.58887, 0),
//OpenMesh::Vec3d(2.70708, 3.03427, 0),
//OpenMesh::Vec3d(2.26469, 3.32868, 0),
//OpenMesh::Vec3d(2.8365, 1.6402, 0),
//OpenMesh::Vec3d(1.43069, 2.817, 0),
//OpenMesh::Vec3d(2.21378, 1.10076, 0),
//OpenMesh::Vec3d(0.599588, 2.35945, 0),
//OpenMesh::Vec3d(0.929425, 1.2952, 0),
//OpenMesh::Vec3d(0.753215, 1.86505, 0),
//OpenMesh::Vec3d(0.548551, 2.41461, 0),
//OpenMesh::Vec3d(0.315433, 2.94388, 0),
//OpenMesh::Vec3d(0.0538616, 3.45286, 0),
//OpenMesh::Vec3d(-0.236164, 3.94155, 0),
//OpenMesh::Vec3d(-0.554642, 4.40996, 0),
//OpenMesh::Vec3d(-0.901575, 4.85807, 0),
//OpenMesh::Vec3d(-0.998202, 4.48056, 0),
//OpenMesh::Vec3d(-1.11029, 4.21697, 0),
//OpenMesh::Vec3d(-1.23785, 4.06728, 0),
//OpenMesh::Vec3d(-1.38087, 4.03151, 0),
//OpenMesh::Vec3d(-1.53935, 4.10964, 0),
//OpenMesh::Vec3d(-1.71329, 4.30169, 0),
//OpenMesh::Vec3d(-1.9027, 4.60764, 0),
//OpenMesh::Vec3d(-2.10238, 4.45848, 0),
//OpenMesh::Vec3d(-2.31984, 4.27976, 0),
//OpenMesh::Vec3d(-2.55508, 4.07148, 0),
//OpenMesh::Vec3d(-2.80809, 3.83365, 0),
//OpenMesh::Vec3d(-3.07887, 3.56626, 0),
//OpenMesh::Vec3d(-3.36743, 3.26932, 0),
//OpenMesh::Vec3d(-3.67377, 2.94282, 0),
//OpenMesh::Vec3d(-3.82336, 2.48576, 0),
//OpenMesh::Vec3d(-3.87572, 2.06805, 0),
//OpenMesh::Vec3d(-3.83083, 1.68968, 0),
//OpenMesh::Vec3d(-3.68872, 1.35067, 0),
//OpenMesh::Vec3d(-3.44936, 1.05102, 0),
//OpenMesh::Vec3d(-3.11277, 0.790709, 0),
//OpenMesh::Vec3d(-2.67894, 0.569752, 0),
//OpenMesh::Vec3d(-3.72316, 1.23398, 0),
//OpenMesh::Vec3d(-3.70137, 0.279413, 0),
//OpenMesh::Vec3d(-3.43965, -0.725036, 0),
//OpenMesh::Vec3d(-4.43595, 1.305, 0),
//OpenMesh::Vec3d(-4.93327, -0.948555, 0),
//OpenMesh::Vec3d(-5.54271, 1.48429, 0),
//OpenMesh::Vec3d(-5.86143, 0.0619877, 0),
//OpenMesh::Vec3d(-5.18066, 1.04357, 0),
//OpenMesh::Vec3d(-5.04011, -1.7729, 0),
//OpenMesh::Vec3d(-4.45911, 0.61759, 0),
//OpenMesh::Vec3d(-3.2615, -1.3835, 0),
//OpenMesh::Vec3d(-3.24322, -0.137838, 0),
//OpenMesh::Vec3d(-2.92926, -1.4303, 0),
//OpenMesh::Vec3d(-2.25433, -1.02231, 0),
//OpenMesh::Vec3d(-2.44286, -1.17952, 0),
//OpenMesh::Vec3d(-2.60975, -1.33066, 0),
//OpenMesh::Vec3d(-2.75498, -1.47573, 0),
//OpenMesh::Vec3d(-2.87856, -1.61471, 0),
//OpenMesh::Vec3d(-2.98048, -1.74762, 0),
//OpenMesh::Vec3d(-3.06076, -1.87446, 0),
//OpenMesh::Vec3d(-3.11938, -1.99522, 0),
//OpenMesh::Vec3d(-2.89369, -2.21886, 0),
//OpenMesh::Vec3d(-2.70879, -2.46198, 0),
//OpenMesh::Vec3d(-2.5647, -2.72457, 0),
//OpenMesh::Vec3d(-2.4614, -3.00665, 0),
//OpenMesh::Vec3d(-2.39891, -3.30822, 0),
//OpenMesh::Vec3d(-2.37721, -3.62926, 0),
//OpenMesh::Vec3d(-2.39632, -3.96978, 0),
//OpenMesh::Vec3d(-1.75496, -4.14205, 0),
//OpenMesh::Vec3d(-1.10645, -4.26301, 0),
//OpenMesh::Vec3d(-0.450789, -4.33267, 0),
//OpenMesh::Vec3d(0.212028, -4.35103, 0),
//OpenMesh::Vec3d(0.881999, -4.31809, 0),
//OpenMesh::Vec3d(1.55912, -4.23385, 0),
//};

//CC_points = {//flower3?
//  OpenMesh::Vec3d(0.060138, -2.34968, 0),
//OpenMesh::Vec3d(0.167533, -1.16113, 0),
//OpenMesh::Vec3d(0.304982, 0.160104, 0),
//OpenMesh::Vec3d(0.231798, 1.08177, 0),
//OpenMesh::Vec3d(0.61252, 2.32314, 0),
//OpenMesh::Vec3d(0.771864, 3.13087, 0),
//OpenMesh::Vec3d(0.855257, 4.5897, 0),
//OpenMesh::Vec3d(3.38975, 5.46352, 0),
//OpenMesh::Vec3d(4.41819, 7.10405, 0),
//OpenMesh::Vec3d(8.67744, 7.8645, 0),
//OpenMesh::Vec3d(4.59444, 11.1053, 0),
//OpenMesh::Vec3d(-2.40215, 10.8965, 0),
//OpenMesh::Vec3d(-7.26403, 7.93585, 0),
//OpenMesh::Vec3d(-4.30306, 6.87478, 0),
//OpenMesh::Vec3d(-2.76766, 5.30325, 0),
//OpenMesh::Vec3d(-0.552003, 4.00657, 0),
//OpenMesh::Vec3d(-0.0860061, 2.21267, 0),
//OpenMesh::Vec3d(0.0878497, 0.616274, 0),
//OpenMesh::Vec3d(-0.0636332, -0.969635, 0),
//OpenMesh::Vec3d(-1.42388, 1.41961, 0),
//OpenMesh::Vec3d(-3.71789, 1.96421, 0),
//OpenMesh::Vec3d(-7.39052, 1.52015, 0),
//OpenMesh::Vec3d(-5.55652, -0.704119, 0),
//OpenMesh::Vec3d(-3.16591, -1.64514, 0),
//OpenMesh::Vec3d(-0.29709, -1.64989, 0),
//OpenMesh::Vec3d(-0.284548, -4.19343, 0),
//OpenMesh::Vec3d(-4.56897, -5.40077, 0),
//OpenMesh::Vec3d(-3.58096, -9.71999, 0),
//OpenMesh::Vec3d(-3.51223, -4.61304, 0),
//OpenMesh::Vec3d(-0.378727, -4.82229, 0),
//OpenMesh::Vec3d(0.125836, -2.84884, 0),
//OpenMesh::Vec3d(1.10361, -2.36519, 0),
//OpenMesh::Vec3d(3.49, -2.80699, 0),
//OpenMesh::Vec3d(7.08069, -3.99408, 0),
//OpenMesh::Vec3d(4.17284, -0.459165, 0),
//OpenMesh::Vec3d(0.512756, -0.0850323, 0),
//};


//CC_points = { //flower 2
//	OpenMesh::Vec3d(0.060138, -2.34968, 0),
//OpenMesh::Vec3d(0.236257, -0.500513, 0),
//OpenMesh::Vec3d(0.231798, 1.08177, 0),
//OpenMesh::Vec3d(0.692192, 2.72701, 0),
//OpenMesh::Vec3d(0.855257, 4.5897, 0),
//OpenMesh::Vec3d(5.47551, 4.08825, 0),
//OpenMesh::Vec3d(8.67744, 7.8645, 0),
//OpenMesh::Vec3d(1.09614, 11.0009, 0),
//OpenMesh::Vec3d(-7.26403, 7.93585, 0),
//OpenMesh::Vec3d(-4.36735, 4.74858, 0),
//OpenMesh::Vec3d(-0.552003, 4.00657, 0),
//OpenMesh::Vec3d(0.0009218, 1.41447, 0),
//OpenMesh::Vec3d(-0.0636332, -0.969635, 0),
//OpenMesh::Vec3d(-2.57089, 1.69191, 0),
//OpenMesh::Vec3d(-7.39052, 1.52015, 0),
//OpenMesh::Vec3d(-4.36121, -1.17463, 0),
//OpenMesh::Vec3d(-0.29709, -1.64989, 0),
//OpenMesh::Vec3d(-4.32186, -5.28243, 0),
//OpenMesh::Vec3d(-3.58096, -9.71999, 0),
//OpenMesh::Vec3d(-3.37836, -5.11055, 0),
//OpenMesh::Vec3d(0.125836, -2.84884, 0),
//OpenMesh::Vec3d(2.94391, -4.43496, 0),
//OpenMesh::Vec3d(7.08069, -3.99408, 0),
//OpenMesh::Vec3d(2.3428, -0.272099, 0),
//};

//CC_points = {//flower 7
//	OpenMesh::Vec3d(0.060138, -2.34968, 0),
//OpenMesh::Vec3d(0.106164, -1.8403, 0),
//OpenMesh::Vec3d(0.156484, -1.31197, 0),
//OpenMesh::Vec3d(0.204221, -0.779887, 0),
//OpenMesh::Vec3d(0.242497, -0.259265, 0),
//OpenMesh::Vec3d(0.264437, 0.234689, 0),
//OpenMesh::Vec3d(0.725381, 0.565134, 0),
//OpenMesh::Vec3d(1.27787, 1.1061, 0),
//OpenMesh::Vec3d(0.394965, 1.61379, 0),
//OpenMesh::Vec3d(0.526506, 2.08385, 0),
//OpenMesh::Vec3d(0.630576, 2.52296, 0),
//OpenMesh::Vec3d(0.711332, 2.96211, 0),
//OpenMesh::Vec3d(0.772927, 3.43229, 0),
//OpenMesh::Vec3d(0.819517, 3.96449, 0),
//OpenMesh::Vec3d(0.855257, 4.5897, 0),
//OpenMesh::Vec3d(1.94147, 4.96419, 0),
//OpenMesh::Vec3d(3.73668, 4.93421, 0),
//OpenMesh::Vec3d(4.52793, 5.48071, 0),
//OpenMesh::Vec3d(4.45056, 6.55665, 0),
//OpenMesh::Vec3d(5.02614, 7.61761, 0),
//OpenMesh::Vec3d(6.38999, 8.06925, 0),
//OpenMesh::Vec3d(8.67744, 7.8645, 0),
//OpenMesh::Vec3d(6.24266, 10.6718, 0),
//OpenMesh::Vec3d(3.86555, 7.90151, 0),
//OpenMesh::Vec3d(2.47202, 11.7662, 0),
//OpenMesh::Vec3d(-0.517681, 9.43055, 0),
//OpenMesh::Vec3d(-3.73159, 11.8117, 0),
//OpenMesh::Vec3d(-3.7706, 7.81967, 0),
//OpenMesh::Vec3d(-7.26403, 7.93585, 0),
//OpenMesh::Vec3d(-5.99504, 7.48111, 0),
//OpenMesh::Vec3d(-4.92971, 6.95344, 0),
//OpenMesh::Vec3d(-3.63227, 7.001, 0),
//OpenMesh::Vec3d(-2.79374, 6.3948, 0),
//OpenMesh::Vec3d(-2.35396, 5.15727, 0),
//OpenMesh::Vec3d(-1.50157, 4.56229, 0),
//OpenMesh::Vec3d(-0.552003, 4.00657, 0),
//OpenMesh::Vec3d(-0.35229, 3.23776, 0),
//OpenMesh::Vec3d(-1.09442, 2.13225, 0),
//OpenMesh::Vec3d(-0.079016, 1.77943, 0),
//OpenMesh::Vec3d(-0.00735184, 1.07923, 0),
//OpenMesh::Vec3d(0.0197323, 0.391214, 0),
//OpenMesh::Vec3d(0.00128804, -0.28996, 0),
//OpenMesh::Vec3d(-0.0636332, -0.969635, 0),
//OpenMesh::Vec3d(-0.646596, 0.0543271, 0),
//OpenMesh::Vec3d(-1.36295, 0.814769, 0),
//OpenMesh::Vec3d(-2.22542, 1.33615, 0),
//OpenMesh::Vec3d(-3.24669, 1.64292, 0),
//OpenMesh::Vec3d(-4.4395, 1.75954, 0),
//OpenMesh::Vec3d(-5.81654, 1.71046, 0),
//OpenMesh::Vec3d(-7.39052, 1.52015, 0),
//OpenMesh::Vec3d(-6.60452, 0.566892, 0),
//OpenMesh::Vec3d(-5.739, -0.203045, 0),
//OpenMesh::Vec3d(-4.79621, -0.799575, 0),
//OpenMesh::Vec3d(-3.77839, -1.23261, 0),
//OpenMesh::Vec3d(-2.68776, -1.51207, 0),
//OpenMesh::Vec3d(-1.52658, -1.64785, 0),
//OpenMesh::Vec3d(-0.29709, -1.64989, 0),
//OpenMesh::Vec3d(-1.72703, -2.95892, 0),
//OpenMesh::Vec3d(0.632428, -4.75824, 0),
//OpenMesh::Vec3d(-0.365144, -5.35037, 0),
//OpenMesh::Vec3d(-3.78949, -4.69212, 0),
//OpenMesh::Vec3d(-3.67462, -6.46236, 0),
//OpenMesh::Vec3d(-1.37705, -8.74468, 0),
//OpenMesh::Vec3d(-3.53231, -9.71999, 0),
//OpenMesh::Vec3d(-0.924157, -8.40708, 0),
//OpenMesh::Vec3d(-3.08422, -6.10206, 0),
//OpenMesh::Vec3d(-3.26623, -4.53686, 0),
//OpenMesh::Vec3d(-0.00293407, -5.54073, 0),
//OpenMesh::Vec3d(0.85041, -5.34761, 0),
//OpenMesh::Vec3d(-0.0904053, -3.6946, 0),
//OpenMesh::Vec3d(0.125836, -2.84884, 0),
//OpenMesh::Vec3d(0.544882, -2.64156, 0),
//OpenMesh::Vec3d(1.16516, -2.56649, 0),
//OpenMesh::Vec3d(1.98083, -2.61848, 0),
//OpenMesh::Vec3d(2.98605, -2.79238, 0),
//OpenMesh::Vec3d(4.175, -3.08304, 0),
//OpenMesh::Vec3d(5.54182, -3.48533, 0),
//OpenMesh::Vec3d(7.08069, -3.99408, 0),
//OpenMesh::Vec3d(5.83447, -2.47912, 0),
//OpenMesh::Vec3d(4.48079, -1.41569, 0),
//OpenMesh::Vec3d(3.13277, -0.788896, 0),
//OpenMesh::Vec3d(1.90357, -0.58381, 0),
//OpenMesh::Vec3d(0.906306, -0.785522, 0),
//OpenMesh::Vec3d(0.254117, -1.37912, 0),
//};



//CC_points = { //flower 7,new
// OpenMesh::Vec3d(0.060138, -2.34968, 0),
//OpenMesh::Vec3d(0.106164, -1.8403, 0),
//OpenMesh::Vec3d(0.156484, -1.31197, 0),
//OpenMesh::Vec3d(0.204221, -0.779887, 0),
//OpenMesh::Vec3d(0.242497, -0.259265, 0),
//OpenMesh::Vec3d(0.264437, 0.234689, 0),
//OpenMesh::Vec3d(0.263163, 0.68677, 0),
//OpenMesh::Vec3d(0.316944, 1.0988, 0),
//OpenMesh::Vec3d(1.13816, 1.9118, 0),
//OpenMesh::Vec3d(0.526506, 2.08385, 0),
//OpenMesh::Vec3d(0.630576, 2.52296, 0),
//OpenMesh::Vec3d(2.33153, 2.20675, 0),
//OpenMesh::Vec3d(2.39312, 2.67693, 0),
//OpenMesh::Vec3d(0.819517, 3.96449, 0),
//OpenMesh::Vec3d(0.855257, 4.5897, 0),
//OpenMesh::Vec3d(2.70697, 4.46184, 0),
//OpenMesh::Vec3d(3.57803, 4.94586, 0),
//OpenMesh::Vec3d(3.60378, 5.99472, 0),
//OpenMesh::Vec3d(4.45056, 6.55665, 0),
//OpenMesh::Vec3d(5.4882, 7.08696, 0),
//OpenMesh::Vec3d(6.85205, 7.53859, 0),
//OpenMesh::Vec3d(8.67744, 7.8645, 0),
//OpenMesh::Vec3d(4.12873, 7.55496, 0),
//OpenMesh::Vec3d(1.96264, 8.45108, 0),
//OpenMesh::Vec3d(2.32342, 10.5728, 0),
//OpenMesh::Vec3d(-0.242403, 10.5431, 0),
//OpenMesh::Vec3d(-3.77254, 12.1138, 0),
//OpenMesh::Vec3d(-3.86467, 8.29567, 0),
//OpenMesh::Vec3d(-7.26403, 7.93585, 0),
//OpenMesh::Vec3d(-5.99504, 7.48111, 0),
//OpenMesh::Vec3d(-4.30774, 7.59933, 0),
//OpenMesh::Vec3d(-3.38589, 7.02118, 0),
//OpenMesh::Vec3d(-3.16933, 5.76909, 0),
//OpenMesh::Vec3d(-3.14338, 4.7506, 0),
//OpenMesh::Vec3d(-2.29099, 4.15562, 0),
//OpenMesh::Vec3d(-0.552003, 4.00657, 0),
//OpenMesh::Vec3d(-0.35229, 3.23776, 0),
//OpenMesh::Vec3d(-1.11388, 2.51905, 0),
//OpenMesh::Vec3d(0.369823, 2.075, 0),
//OpenMesh::Vec3d(0.14591, 1.53901, 0),
//OpenMesh::Vec3d(-0.661432, 0.618268, 0),
//OpenMesh::Vec3d(-0.679876, -0.062905, 0),
//OpenMesh::Vec3d(-0.0636332, -0.969635, 0),
//OpenMesh::Vec3d(-0.646596, 0.0543271, 0),
//OpenMesh::Vec3d(-1.36295, 0.814769, 0),
//OpenMesh::Vec3d(-2.22542, 1.33615, 0),
//OpenMesh::Vec3d(-3.24669, 1.64292, 0),
//OpenMesh::Vec3d(-4.4395, 1.75954, 0),
//OpenMesh::Vec3d(-5.81654, 1.71046, 0),
//OpenMesh::Vec3d(-7.39052, 1.52015, 0),
//OpenMesh::Vec3d(-6.60452, 0.566892, 0),
//OpenMesh::Vec3d(-5.739, -0.203045, 0),
//OpenMesh::Vec3d(-4.79621, -0.799575, 0),
//OpenMesh::Vec3d(-3.77839, -1.23261, 0),
//OpenMesh::Vec3d(-2.68776, -1.51207, 0),
//OpenMesh::Vec3d(-1.52658, -1.64785, 0),
//OpenMesh::Vec3d(-0.29709, -1.64989, 0),
//OpenMesh::Vec3d(-1.28548, -3.01893, 0),
//OpenMesh::Vec3d(1.314, -5.34777, 0),
//OpenMesh::Vec3d(-0.157957, -5.71244, 0),
//OpenMesh::Vec3d(-4.15521, -5.37328, 0),
//OpenMesh::Vec3d(-3.67462, -6.46236, 0),
//OpenMesh::Vec3d(-1.05795, -7.9212, 0),
//OpenMesh::Vec3d(-3.58096, -9.71999, 0),
//OpenMesh::Vec3d(-0.605064, -7.5836, 0),
//OpenMesh::Vec3d(-3.08422, -6.10206, 0),
//OpenMesh::Vec3d(-3.63195, -5.21803, 0),
//OpenMesh::Vec3d(0.204253, -5.9028, 0),
//OpenMesh::Vec3d(1.53198, -5.93714, 0),
//OpenMesh::Vec3d(-0.0904053, -3.6946, 0),
//OpenMesh::Vec3d(0.125836, -2.84884, 0),
//OpenMesh::Vec3d(0.544882, -2.64156, 0),
//OpenMesh::Vec3d(1.16516, -2.56649, 0),
//OpenMesh::Vec3d(1.98083, -2.61848, 0),
//OpenMesh::Vec3d(2.98605, -2.79238, 0),
//OpenMesh::Vec3d(4.175, -3.08304, 0),
//OpenMesh::Vec3d(5.54182, -3.48533, 0),
//OpenMesh::Vec3d(7.08069, -3.99408, 0),
//OpenMesh::Vec3d(5.83447, -2.47912, 0),
//OpenMesh::Vec3d(4.48079, -1.41569, 0),
//OpenMesh::Vec3d(3.13277, -0.788896, 0),
//OpenMesh::Vec3d(1.90357, -0.58381, 0),
//OpenMesh::Vec3d(0.906306, -0.785522, 0),
//OpenMesh::Vec3d(0.254117, -1.37912, 0),
//};

//CC_points = { //octcpos
//OpenMesh::Vec3d(14.1163, -2.13572, 0),
//OpenMesh::Vec3d(11.6986, -2.58835, 0),
//OpenMesh::Vec3d(9.93025, -2.04599, 0),
//OpenMesh::Vec3d(6.50689, -3.68536, 0),
//OpenMesh::Vec3d(5.42654, -1.87831, 0),
//OpenMesh::Vec3d(4.92811, 1.19895, 0),
//OpenMesh::Vec3d(6.11625, 1.44333, 0),
//OpenMesh::Vec3d(6.71954, -0.124826, 0),
//OpenMesh::Vec3d(7.67902, 0.404205, 0),
//OpenMesh::Vec3d(8.47221, 1.76796, 0),
//OpenMesh::Vec3d(9.05846, 3.72884, 0),
//OpenMesh::Vec3d(8.84838, 5.3717, 0),
//OpenMesh::Vec3d(6.49071, 5.0889, 0),
//OpenMesh::Vec3d(7.77231, 9.87716, 0),
//OpenMesh::Vec3d(5.30888, 14.7788, 0),
//OpenMesh::Vec3d(-0.841746, 14.3125, 0),
//OpenMesh::Vec3d(-5.90767, 11.3197, 0),
//OpenMesh::Vec3d(-8.86451, 3.55846, 0),
//OpenMesh::Vec3d(-0.246632, 4.78125, 0),
//OpenMesh::Vec3d(1.06736, 1.00327, 0),
//OpenMesh::Vec3d(-4.23325, 0.219585, 0),
//OpenMesh::Vec3d(-11.8515, 4.19812, 0),
//OpenMesh::Vec3d(-12.0648, 2.31314, 0),
//OpenMesh::Vec3d(-10.1253, 0.317715, 0),
//OpenMesh::Vec3d(-6.06773, 0.106026, 0),
//OpenMesh::Vec3d(-7.90343, -1.09784, 0),
//OpenMesh::Vec3d(-12.4609, -0.711102, 0),
//OpenMesh::Vec3d(-10.8011, -3.7199, 0),
//OpenMesh::Vec3d(-11.7545, -4.62186, 0),
//OpenMesh::Vec3d(-11.8664, -5.95728, 0),
//OpenMesh::Vec3d(-11.2416, -8.89297, 0),
//OpenMesh::Vec3d(-9.92311, -4.87563, 0),
//OpenMesh::Vec3d(-4.59648, -5.33218, 0),
//OpenMesh::Vec3d(-2.05936, -4.13197, 0),
//OpenMesh::Vec3d(-3.09482, -7.04246, 0),
//OpenMesh::Vec3d(-13.654, -11.0081, 0),
//OpenMesh::Vec3d(-6.70591, -12.2449, 0),
//OpenMesh::Vec3d(-7.3868, -10.8175, 0),
//OpenMesh::Vec3d(-8.22787, -9.63252, 0),
//OpenMesh::Vec3d(-6.34981, -9.67327, 0),
//OpenMesh::Vec3d(-6.18194, -10.5259, 0),
//OpenMesh::Vec3d(-5.16694, -10.9844, 0),
//OpenMesh::Vec3d(-4.43513, -10.705, 0),
//OpenMesh::Vec3d(-2.69772, -8.66487, 0),
//OpenMesh::Vec3d(-1.78497, -7.37903, 0),
//OpenMesh::Vec3d(0.488619, -4.78304, 0),
//OpenMesh::Vec3d(-0.0659694, -6.98257, 0),
//OpenMesh::Vec3d(-0.438017, -8.1058, 0),
//OpenMesh::Vec3d(2.47188, -10.8048, 0),
//OpenMesh::Vec3d(1.71534, -8.54231, 0),
//OpenMesh::Vec3d(1.34671, -7.00428, 0),
//OpenMesh::Vec3d(1.89121, -5.32397, 0),
//OpenMesh::Vec3d(3.61222, -8.87285, 0),
//OpenMesh::Vec3d(5.47686, -9.79907, 0),
//OpenMesh::Vec3d(4.93164, -12.9952, 0),
//OpenMesh::Vec3d(8.7507, -9.71153, 0),
//OpenMesh::Vec3d(2.57694, -5.41496, 0),
//OpenMesh::Vec3d(4.12527, -0.841089, 0),
//OpenMesh::Vec3d(4.43461, -3.18152, 0),
//OpenMesh::Vec3d(5.09308, -5.72131, 0),
//OpenMesh::Vec3d(7.13664, -5.11346, 0),
//OpenMesh::Vec3d(9.72134, -3.74938, 0),
//OpenMesh::Vec3d(11.762, -4.0687, 0),
//OpenMesh::Vec3d(14.823, -3.17095, 0),
//OpenMesh::Vec3d(15.2188, -3.22438, 0),
//OpenMesh::Vec3d(14.9258, -2.34444, 0),
//};

//CC_points = {//sheji
//	OpenMesh::Vec3d(-15.1181, -7.27403, 0),
//OpenMesh::Vec3d(-7.47435, -7.27403, 0),
//OpenMesh::Vec3d(-1.31281, -7.27403, 0),
//OpenMesh::Vec3d(1.437, -7.27403, 0),
//OpenMesh::Vec3d(1.39394, -6.40042, 0),
//OpenMesh::Vec3d(1.437, 1.82655, 0),
//OpenMesh::Vec3d(1.437, 5.21889, 0),
//OpenMesh::Vec3d(-4.08329, 5.21889, 0),
//OpenMesh::Vec3d(-10.0778, 5.21889, 0),
//OpenMesh::Vec3d(-15.1181, 5.21889, 0),
//OpenMesh::Vec3d(-15.0535, 1.26594, 0),
//OpenMesh::Vec3d(-15.075, -6.50998, 0),
//};
//CC_points = {//sheji polyGC
//	OpenMesh::Vec3d(-15.1181, -7.31403, 0),
//	OpenMesh::Vec3d(-13.1181, -7.36403, 0),
//	OpenMesh::Vec3d(-11.1181, -7.41403, 0),
//OpenMesh::Vec3d(-7.47435, -7.28, 0),
//OpenMesh::Vec3d(-3.31281, -7.35403, 0),
//OpenMesh::Vec3d(-1.31281, -7.59403, 0),
//OpenMesh::Vec3d(1.437, -7.71403, 0),
//OpenMesh::Vec3d(1.39394, -6.40042, 0),
//OpenMesh::Vec3d(1.437, 1.82655, 0),
//OpenMesh::Vec3d(1.437, 5.21889, 0),
//OpenMesh::Vec3d(-4.08329, 5.21889, 0),
//OpenMesh::Vec3d(-10.0778, 5.21889, 0),
//OpenMesh::Vec3d(-15.1181, 5.21889, 0),
//OpenMesh::Vec3d(-15.0535, 1.26594, 0),
//OpenMesh::Vec3d(-15.075, -6.50998, 0),
//};
//CC_points = { //shejiCVM
// OpenMesh::Vec3d(-12.729, -6.31988, 0),
//OpenMesh::Vec3d(-10.997, -6.3684, 0),
//OpenMesh::Vec3d(-8.98264, -6.41056, 0),
//OpenMesh::Vec3d(-6.73546, -6.27695, 0),
//OpenMesh::Vec3d(-4.10206, -6.35682, 0),
//OpenMesh::Vec3d(-1.83219, -6.59565, 0),
//OpenMesh::Vec3d(1.47363, -6.71038, 0),
//OpenMesh::Vec3d(1.91409, -1.48812, 0),
//OpenMesh::Vec3d(1.91601, 1.91747, 0),
//OpenMesh::Vec3d(2.24009, 5.19837, 0),
//OpenMesh::Vec3d(-4.52556, 4.91904, 0),
//OpenMesh::Vec3d(-9.04531, 4.60278, 0),
//OpenMesh::Vec3d(-12.3981, 4.54636, 0),
//OpenMesh::Vec3d(-12.5084, 0.92428, 0),
//OpenMesh::Vec3d(-12.6187, -2.6978, 0),
//OpenMesh::Vec3d(-12.729, -6.31988, 0),
//};

//CC_points = { //pangxie
//	OpenMesh::Vec3d(7.71977, -2.05287, 0),
//OpenMesh::Vec3d(14.0104, 0.411218, 0),
//OpenMesh::Vec3d(19.374, -4.0279, 0),
//OpenMesh::Vec3d(22.4267, -11.7043, 0),
//OpenMesh::Vec3d(22.6234, -4.69254, 0),
//OpenMesh::Vec3d(16.5132, 3.32898, 0),
//OpenMesh::Vec3d(7.68061, 1.27629, 0),
//OpenMesh::Vec3d(11.8108, 5.72265, 0),
//OpenMesh::Vec3d(24.3071, 0.0673024, 0),
//OpenMesh::Vec3d(25.5395, -8.6918, 0),
//OpenMesh::Vec3d(24.2683, 1.50686, 0),
//OpenMesh::Vec3d(15.8142, 7.55389, 0),
//OpenMesh::Vec3d(7.82502, 4.58309, 0),
//OpenMesh::Vec3d(13.0287, 7.31829, 0),
//OpenMesh::Vec3d(20.2624, 6.36308, 0),
//OpenMesh::Vec3d(25.4787, 9.17267, 0),
//OpenMesh::Vec3d(18.4321, 10.7327, 0),
//OpenMesh::Vec3d(13.6405, 8.00646, 0),
//OpenMesh::Vec3d(4.85222, 5.92983, 0),
//OpenMesh::Vec3d(10.4038, 10.2729, 0),
//OpenMesh::Vec3d(14.7017, 15.1746, 0),
//OpenMesh::Vec3d(12.0969, 19.6625, 0),
//OpenMesh::Vec3d(11.6911, 15.6914, 0),
//OpenMesh::Vec3d(8.19465, 10.9663, 0),
//OpenMesh::Vec3d(2.16398, 7.18829, 0),
//OpenMesh::Vec3d(-0.0248164, 9.29763, 0),
//OpenMesh::Vec3d(-1.17408, 9.21065, 0),
//OpenMesh::Vec3d(-3.33054, 6.62634, 0),
//OpenMesh::Vec3d(-6.47158, 9.13021, 0),
//OpenMesh::Vec3d(-13.0069, 12.1429, 0),
//OpenMesh::Vec3d(-13.8588, 18.6794, 0),
//OpenMesh::Vec3d(-18.8441, 13.6221, 0),
//OpenMesh::Vec3d(-12.6782, 8.694, 0),
//OpenMesh::Vec3d(-5.13643, 5.07839, 0),
//OpenMesh::Vec3d(-5.27052, 4.97835, 0),
//OpenMesh::Vec3d(-5.46014, 4.84129, 0),
//OpenMesh::Vec3d(-5.64977, 4.70424, 0),
//OpenMesh::Vec3d(-16.6586, 6.71951, 0),
//OpenMesh::Vec3d(-22.5467, 6.17214, 0),
//OpenMesh::Vec3d(-27.2511, 11.2113, 0),
//OpenMesh::Vec3d(-24.0825, 3.01043, 0),
//OpenMesh::Vec3d(-15.7064, 4.55985, 0),
//OpenMesh::Vec3d(-7.0058, 3.14714, 0),
//OpenMesh::Vec3d(-17.113, 3.54583, 0),
//OpenMesh::Vec3d(-25.8122, -0.336797, 0),
//OpenMesh::Vec3d(-28.889, -9.77593, 0),
//OpenMesh::Vec3d(-21.7434, -1.11048, 0),
//OpenMesh::Vec3d(-15.7677, 1.49321, 0),
//OpenMesh::Vec3d(-7.4344, 0.51362, 0),
//OpenMesh::Vec3d(-18.7285, 0.0406011, 0),
//OpenMesh::Vec3d(-23.595, -7.64864, 0),
//OpenMesh::Vec3d(-19.4429, -13.2397, 0),
//OpenMesh::Vec3d(-19.1456, -5.33267, 0),
//OpenMesh::Vec3d(-13.9446, -1.15383, 0),
//OpenMesh::Vec3d(-7.53781, -3.00233, 0),
//OpenMesh::Vec3d(-16.6865, -2.97738, 0),
//OpenMesh::Vec3d(-18.0268, -14.6281, 0),
//OpenMesh::Vec3d(-2.51255, -20.4804, 0),
//OpenMesh::Vec3d(-4.32991, -14.8028, 0),
//OpenMesh::Vec3d(-7.41556, -7.88393, 0),
//OpenMesh::Vec3d(-11.3053, -6.94069, 0),
//OpenMesh::Vec3d(-7.82506, -6.74206, 0),
//OpenMesh::Vec3d(-4.0724, -6.50565, 0),
//OpenMesh::Vec3d(0.270762, -12.5408, 0),
//OpenMesh::Vec3d(4.64637, -10.1421, 0),
//OpenMesh::Vec3d(7.19738, -5.76408, 0),
//OpenMesh::Vec3d(11.4721, -5.97088, 0),
//OpenMesh::Vec3d(9.8597, -10.8824, 0),
//OpenMesh::Vec3d(3.22116, -19.9134, 0),
//OpenMesh::Vec3d(-0.202762, -25.0171, 0),
//OpenMesh::Vec3d(1.82509, -36.1817, 0),
//OpenMesh::Vec3d(10.8714, -18.1249, 0),
//OpenMesh::Vec3d(15.5892, -11.8939, 0),
//OpenMesh::Vec3d(18.353, -8.94576, 0),
//OpenMesh::Vec3d(15.939, -2.21534, 0),
//};
	CC_mesh = createMeshFromCurveCage(CC_points);
	auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
	for (auto vh : CC_mesh.vertices())
	{
		CagevertexState[vh] = NotSelected;
	}
	if (!highdegree)
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
	else
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, todegree);
	deform_mesh_from_cc(deformedMesh);//通过重心坐标改变mesh的顶点位置
	MeshTools::AssignPoints(mesh, deformedMesh);
	update();
}

bool MeshViewerWidget::event(QEvent* _event)
{
	if (_event->type() == QEvent::MouseButtonPress)
	{
		if (selectMode == Move)
		{		    
			auto e = static_cast<QMouseEvent*>(_event);
			QPoint winCor = e->pos();
			double depth;
			OpenMesh::Vec3d objCor;
			WinCor2ObjCor(winCor.x(), winCor.y(), objCor, depth);
			objCor = OpenMesh::Vec3d{ objCor[0],objCor[1],0 };

			OpenMesh::VertexHandle minVh;
			if (NearestCageVertex(objCor, minVh))
			{
				auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
				if (vertexState[minVh] == Custom)
				{
					isMovable = true;
					moveDepth = depth;
					lastObjCor = objCor;
					return true;
				}
			}
		}
		if (selectMode == SelectAdjust)
		{
			auto e = static_cast<QMouseEvent*>(_event);
			QPoint winCor = e->pos();
			double depth;
			OpenMesh::Vec3d objCor;
			WinCor2ObjCor(winCor.x(), winCor.y(), objCor, depth);
			objCor = OpenMesh::Vec3d{ objCor[0],objCor[1],0 };

			OpenMesh::VertexHandle minVh;
			if (NearestCageVertex(objCor, minVh))
			{
				auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
				if (vertexState[minVh] == Custom)
				{
					isMovable = true;
					moveDepth = depth;
					lastObjCor = objCor;
					return true;
				}
			}
		}
	}
	else if (_event->type() == QEvent::MouseMove)
	{
		
		if (selectMode == Move && isMovable)
		{
			
			auto e = static_cast<QMouseEvent*>(_event);
			QPoint winCor = e->pos();
			OpenMesh::Vec3d objCor;
			WinCor2ObjCor(winCor.x(), winCor.y(), moveDepth, objCor);
			objCor = OpenMesh::Vec3d{ objCor[0],objCor[1],0 };

			auto moveVec = objCor - lastObjCor;
			lastObjCor = objCor;
			Mesh deformedMesh;
			deformedMesh.assign(mesh);
			auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(mesh, "vertexState");
			for (auto vh_ : deformedMesh.vertices())
			{
				if (vertexState[mesh.vertex_handle(vh_.idx())] != Custom)
				{
					deformedMesh.set_point(vh_, deformedMesh.point(vh_) + moveVec);
				}
			}
			
			MeshTools::AssignPoints(mesh, deformedMesh);
			//for cage mesh
			//Mesh deformedCCMesh;
			//deformedCCMesh.assign(CC_mesh);
			auto CCvertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
			for (auto vh_ : CC_mesh.vertices())
			{
				if (CCvertexState[CC_mesh.vertex_handle(vh_.idx())] == Custom)
				{
					CC_mesh.set_point(vh_, CC_mesh.point(vh_) + moveVec);
					int vertex_id = vh_.idx();
					CC_points[vertex_id] = CC_points[vertex_id] + moveVec;
				}
			}
			
			
			//MeshTools::AssignPoints(CC_mesh, deformedCCMesh);
			if (!highdegree)
			{
				auto deformedcurvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
				curvecage2 = deformedcurvecage2;
			}
			else
			{
				auto deformedcurvecage2 = CCpoints_fromCCmesh(CC_mesh, todegree);
				curvecage2 = deformedcurvecage2;
			}
			deform_mesh_from_cc(deformedMesh);//通过重心坐标改变mesh的顶点位置
			MeshTools::AssignPoints(mesh, deformedMesh);
			update();
			return true;
		}

		if (selectMode == SelectAdjust && isMovable)
		{

			auto e = static_cast<QMouseEvent*>(_event);
			QPoint winCor = e->pos();
			OpenMesh::Vec3d objCor;
			WinCor2ObjCor(winCor.x(), winCor.y(), moveDepth, objCor);
			objCor = OpenMesh::Vec3d{ objCor[0],objCor[1],0 };

			auto moveVec = objCor - lastObjCor;
			lastObjCor = objCor;
			moveVec = moveVec / 2;
			//for cage mesh
			//Mesh deformedCCMesh;
			//deformedCCMesh.assign(CC_mesh);
			auto CCvertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
			for (auto vh_ : CC_mesh.vertices())
			{
				if (CCvertexState[CC_mesh.vertex_handle(vh_.idx())] == Custom )
				{
					CC_mesh.set_point(vh_, CC_mesh.point(vh_) + moveVec);
					int vertex_id = vh_.idx();
					CC_points[vertex_id] = CC_points[vertex_id] + moveVec;
				}
			}


			//MeshTools::AssignPoints(CC_mesh, deformedCCMesh);
			auto deformedcurvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
			curvecage2 = deformedcurvecage2;
			
			update();
			return true;
		}
	}
	else if (_event->type() == QEvent::MouseButtonRelease)
	{
		
		if (selectMode == Move && isMovable)
		{
			std::cout << "{";
			for (const auto& point : CC_points) {
				std::cout << "OpenMesh::Vec3d(" << point[0] << ", " << point[1] << ", " << point[2] << ")," << std::endl;
			}
			std::cout << "}" << std::endl;
			isMovable = false;
			return true;
		}
		if (selectMode == SelectAdjust)
		{
			std::cout << "{";
			for (const auto& point : CC_points) {
				std::cout << "OpenMesh::Vec3d(" << point[0] << ", " << point[1] << ", " << point[2] << ")," << std::endl;
			}
			std::cout << "}" << std::endl;
			isMovable = false;
			return true;
		}
	}
	return QGLViewerWidget::event(_event);
}

void MeshViewerWidget::mouseDoubleClickEvent(QMouseEvent* _event)
{
	switch (selectMode)
	{
	case NoSelect:
		break;
	case SelectAdjust:
	case SelectCustom: 
	{
		QPoint winCor = _event->pos();
		//winCor.setX(609);
		//winCor.setY(183);
		double depth;
		OpenMesh::Vec3d objCor;
		WinCor2ObjCor(winCor.x(), winCor.y(), objCor, depth);
		OpenMesh::VertexHandle minVh;

		if (NearestCageVertex(objCor,minVh))
		{
			auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
			if (selectMode == SelectAdjust)
			{
				if (vertexState[minVh] == NotSelected)
				{
					vertexState[minVh] = Custom;
				}
				else if (vertexState[minVh] == Custom)
				{
					vertexState[minVh] = NotSelected;
				}
			}
			else if (selectMode == SelectCustom)
			{
				if (vertexState[minVh] == NotSelected)
				{
					vertexState[minVh] = Custom;
				}
				else if (vertexState[minVh] == Custom)
				{
					vertexState[minVh] = NotSelected;
				}
			}
			update();
		}
		break;
	}
	case Move:
		break;
	default:
		break;
	}
}

void MeshViewerWidget::DrawScene(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(&projectionmatrix[0]);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(&modelviewmatrix[0]);
	//DrawAxis();
	if (isDrawBoundingBox) DrawBoundingBox();
	if (isDrawBoundary) DrawBoundary();
	if (isEnableLighting) glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, isTwoSideLighting);
	DrawSceneMesh();
	if (isEnableLighting) glDisable(GL_LIGHTING);
}

void MeshViewerWidget::DrawSceneMesh(void)
{
	if (mesh.n_vertices() == 0) { return; }
	SetMaterial();
	switch (drawmode)
	{
	case POINTS:
		DrawPoints();
		break;
	case WIREFRAME:
		DrawWireframe();
		break;
	case HIDDENLINES:
		DrawHiddenLines();
		break;
	case FLATLINES:
		DrawFlatLines();
		break;
	case FLAT:
		glColor3d(0.8, 0.8, 0.8);
		DrawFlat();
		break;
	case SMOOTH:
		DrawSmooth();
		break;
	case CURVECAGE:
		DrawCurveCage();
		break;
	default:
		break;
	}
}

void MeshViewerWidget::CurveCage_Test(void)
{
	drawmode = CURVECAGE;
	degree = 3;

	if (degree==2)//if mode
	{
		
		int N = 16;
		curvecage2.resize(N);
		std::vector<OpenMesh::Vec3d> square_vertices = {
			OpenMesh::Vec3d(2, -2, 0),
			OpenMesh::Vec3d(2, 2, 0),
			OpenMesh::Vec3d(-2, 2, 0),
			OpenMesh::Vec3d(-2, -2, 0)
		};
		int segments_per_side = N / 4;
		for (int i = 0; i < 4; i++) {
			OpenMesh::Vec3d start = square_vertices[i];
			OpenMesh::Vec3d end = square_vertices[(i + 1) % 4];

			for (int j = 0; j < segments_per_side; j++) {
				double t1 = (double)j / segments_per_side;
				double t2 = (double)(j + 1) / segments_per_side;

				// Compute the coordinates for the control points
				OpenMesh::Vec3d p1 = start + (end - start) * t1;
				OpenMesh::Vec3d p4 = start + (end - start) * t2;

				// Intermediate control points (simple linear interpolation for demonstration)
				OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
				OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

				curvecage2[i * segments_per_side + j] = { p1, p2,  p4 };
			}
		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}
		CC_points = { //human
			OpenMesh::Vec3d(2.2434, -4.09831, 0),
OpenMesh::Vec3d(1.33897, -1.2763, 0),
OpenMesh::Vec3d(0.168209, 0.44022, 0),
OpenMesh::Vec3d(1.4294, -0.288948, 0),
OpenMesh::Vec3d(2.68837, -0.575573, 0),
OpenMesh::Vec3d(3.60324, 1.3499, 0),
OpenMesh::Vec3d(5.97756, 1.52252, 0),
OpenMesh::Vec3d(6.02737, 3.22644, 0),
OpenMesh::Vec3d(5.81575, 4.86821, 0),
OpenMesh::Vec3d(4.28753, 5.19497, 0),
OpenMesh::Vec3d(3.18159, 4.91388, 0),
OpenMesh::Vec3d(2.91693, 3.88495, 0),
OpenMesh::Vec3d(3.06766, 3.08017, 0),
OpenMesh::Vec3d(1.80562, 2.91953, 0),
OpenMesh::Vec3d(0.929425, 1.2952, 0),
OpenMesh::Vec3d(0.312689, 3.28967, 0),
OpenMesh::Vec3d(-0.901575, 4.85807, 0),
OpenMesh::Vec3d(-1.23977, 3.5368, 0),
OpenMesh::Vec3d(-1.9027, 4.60764, 0),
OpenMesh::Vec3d(-2.60159, 4.08557, 0),
OpenMesh::Vec3d(-3.67377, 2.94282, 0),
OpenMesh::Vec3d(-4.19734, 1.3431, 0),
OpenMesh::Vec3d(-2.67894, 0.569752, 0),
OpenMesh::Vec3d(-4.59146, 0.273326, 0),
OpenMesh::Vec3d(-5.02472, -2.5631, 0),
OpenMesh::Vec3d(-3.32428, -1.2692, 0),
OpenMesh::Vec3d(-2.25433, -1.02231, 0),
OpenMesh::Vec3d(-2.9142, -1.57256, 0),
OpenMesh::Vec3d(-3.11938, -1.99522, 0),
OpenMesh::Vec3d(-2.32945, -2.77795, 0),
OpenMesh::Vec3d(-2.39632, -3.96978, 0),
OpenMesh::Vec3d(-0.151571, -4.57271, 0),
		};

		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_green_weight222();//计算2次cage到2次控制点的权重
		
	}
	else if (degree == 3) {
		/*{
			complex<double>x1, x2, x3;
			Cardano(complex<double>(1, 1), complex<double>(2, 1), complex<double>(3, 1), complex<double>(4, 1), x1, x2, x3);
			std::cout << "testing" << std::endl;
			exit(-1);
		}*/
		int N = 4;//deer16球4，giraffe8,tower4,xiyi8,yu4,kuzi 8,car 4,tuzi 8.daizi 4,fish12,choufish12,circle 8,human 16,flower12,octopus16,sheji 4
		curvecage2.resize(N);
		/*curvecage2[0] = { OpenMesh::Vec3d(1.5,0.5,0),OpenMesh::Vec3d(1,1.5,0), OpenMesh::Vec3d(0,1.5,0),OpenMesh::Vec3d(-0.5,0.5,0) };
		curvecage2[1] = { OpenMesh::Vec3d(-0.5,0.5,0),OpenMesh::Vec3d(0,-0.5,0),OpenMesh::Vec3d(1,-0.5,0), OpenMesh::Vec3d(1.5,0.5,0) };*/
		//月亮obj的简单曲边cage
		/*curvecage2[0] = { OpenMesh::Vec3d(0.688794, - 5.48505, 0),OpenMesh::Vec3d(11.2451, - 0.744474,0), OpenMesh::Vec3d(-0.186055, 10.7086,0),OpenMesh::Vec3d(-5.20355, 1.80194,0) };
		curvecage2[1] = { OpenMesh::Vec3d(-5.20355, 1.80194,0),OpenMesh::Vec3d(-0.151066, 2.44716,0),OpenMesh::Vec3d(2.57412, - 0.470955,0), OpenMesh::Vec3d(0.688794, -5.48505,0) };*/
		// 定义正方形的四个顶点
		std::vector<OpenMesh::Vec3d> square_vertices = {
			OpenMesh::Vec3d(2, -2, 0),
			OpenMesh::Vec3d(2, 2, 0),
			OpenMesh::Vec3d(-2, 2, 0),
			OpenMesh::Vec3d(-2, -2, 0)
		};
		int segments_per_side = N / 4;
		for (int i = 0; i < 4; i++) {
			OpenMesh::Vec3d start = square_vertices[i];
			OpenMesh::Vec3d end = square_vertices[(i + 1) % 4];

			for (int j = 0; j < segments_per_side; j++) {
				double t1 = (double)j / segments_per_side;
				double t2 = (double)(j + 1) / segments_per_side;

				// Compute the coordinates for the control points
				OpenMesh::Vec3d p1 = start + (end - start) * t1;
				OpenMesh::Vec3d p4 = start + (end - start) * t2;

				// Intermediate control points (simple linear interpolation for demonstration)
				OpenMesh::Vec3d p2 = p1 + (p4 - p1) /3;
				OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2/3;

				curvecage2[i * segments_per_side + j] = { p1, p2, p3, p4 };
			}
		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}
		CC_points = { //deer
			OpenMesh::Vec3d(4.5276, -10.0723, 0),
OpenMesh::Vec3d(4.89459, -10.0739, 0),
OpenMesh::Vec3d(5.18646, -10.0632, 0),
OpenMesh::Vec3d(5.55932, -10.0582, 0),
OpenMesh::Vec3d(6.25252, -6.75231, 0),
OpenMesh::Vec3d(6.71007, -5.29668, 0),
OpenMesh::Vec3d(6.70587, -3.83298, 0),
OpenMesh::Vec3d(5.83477, -2.3301, 0),
OpenMesh::Vec3d(7.05116, -0.928495, 0),
OpenMesh::Vec3d(6.61916, 0.630525, 0),
OpenMesh::Vec3d(8.02456, 2.09566, 0),
OpenMesh::Vec3d(8.65167, 3.70863, 0),
OpenMesh::Vec3d(6.40582, 5.04416, 0),
OpenMesh::Vec3d(6.37169, 3.87665, 0),
OpenMesh::Vec3d(5.80833, 2.98216, 0),
OpenMesh::Vec3d(5.62158, 2.09523, 0),
OpenMesh::Vec3d(3.23511, 2.54016, 0),
OpenMesh::Vec3d(0.750374, 1.38221, 0),
OpenMesh::Vec3d(-1.69483, 2.49143, 0),
OpenMesh::Vec3d(-3.65702, 4.17105, 0),
OpenMesh::Vec3d(-2.17201, 6.0254, 0),
OpenMesh::Vec3d(-3.05934, 7.31317, 0),
OpenMesh::Vec3d(-1.89982, 8.15725, 0),
OpenMesh::Vec3d(-1.63298, 9.67336, 0),
OpenMesh::Vec3d(-2.82579, 10.398, 0),
OpenMesh::Vec3d(-3.92575, 10.5965, 0),
OpenMesh::Vec3d(-2.76592, 8.26476, 0),
OpenMesh::Vec3d(-5.78046, 7.07922, 0),
OpenMesh::Vec3d(-6.48863, 6.34431, 0),
OpenMesh::Vec3d(-5.30909, 5.86309, 0),
OpenMesh::Vec3d(-7.14945, 4.31455, 0),
OpenMesh::Vec3d(-6.72589, 3.45316, 0),
OpenMesh::Vec3d(-5.35703, 4.33777, 0),
OpenMesh::Vec3d(-4.45691, 4.19383, 0),
OpenMesh::Vec3d(-3.65641, 2.89317, 0),
OpenMesh::Vec3d(-4.13472, -0.366618, 0),
OpenMesh::Vec3d(-3.00663, -1.36788, 0),
OpenMesh::Vec3d(-3.10867, -2.51138, 0),
OpenMesh::Vec3d(-2.47899, -3.00126, 0),
OpenMesh::Vec3d(-2.19824, -3.3226, 0),
OpenMesh::Vec3d(-1.97104, -5.31759, 0),
OpenMesh::Vec3d(-1.18862, -7.33088, 0),
OpenMesh::Vec3d(-2.64391, -10.14, 0),
OpenMesh::Vec3d(-2.4061, -10.1295, 0),
OpenMesh::Vec3d(-2.13911, -10.1257, 0),
OpenMesh::Vec3d(-1.8048, -10.1125, 0),
OpenMesh::Vec3d(-1.04107, -7.68947, 0),
OpenMesh::Vec3d(-1.53657, -5.56524, 0),
OpenMesh::Vec3d(-0.549384, -3.05351, 0),
OpenMesh::Vec3d(0.598589, -3.14832, 0),
OpenMesh::Vec3d(1.97528, -3.33699, 0),
OpenMesh::Vec3d(3.40852, -2.90979, 0),
OpenMesh::Vec3d(4.44002, -4.23804, 0),
OpenMesh::Vec3d(6.07989, -4.12357, 0),
OpenMesh::Vec3d(5.53946, -7.39734, 0),
OpenMesh::Vec3d(5.59344, -8.20186, 0),
OpenMesh::Vec3d(4.89491, -9.35694, 0),
		};
	


//		CC_points = {//ball
//			OpenMesh::Vec3d(4.77039, -5.28648, 0),
//OpenMesh::Vec3d(7.50588, -2.64831, 0),
//OpenMesh::Vec3d(8.28221, 1.94099, 0),
//OpenMesh::Vec3d(4.89974, 5.07107, 0),
//OpenMesh::Vec3d(1.63539, 8.05301, 0),
//OpenMesh::Vec3d(-3.23734, 7.6858, 0),
//OpenMesh::Vec3d(-5.78844, 4.26243, 0),
//OpenMesh::Vec3d(-7.71271, 1.69334, 0),
//OpenMesh::Vec3d(-7.84599, -2.21456, 0),
//OpenMesh::Vec3d(-4.89242, -5.17084, 0),
//OpenMesh::Vec3d(-1.80051, -8.01358, 0),
//OpenMesh::Vec3d(2.59902, -7.54285, 0),
//		};

//		CC_points = {//giraffe
//		OpenMesh::Vec3d(4.26544, -2.08319, 0),
//OpenMesh::Vec3d(2.14808, -1.76166, 0),
//OpenMesh::Vec3d(1.05683, 0.58523, 0),
//OpenMesh::Vec3d(-0.959219, 1.60384, 0),
//OpenMesh::Vec3d(-2.08942, 6.26007, 0),
//OpenMesh::Vec3d(1.37963, 11.5202, 0),
//OpenMesh::Vec3d(5.01997, 9.59617, 0),
//OpenMesh::Vec3d(5.45668, 12.4244, 0),
//OpenMesh::Vec3d(2.73499, 14.4063, 0),
//OpenMesh::Vec3d(0.408011, 13.9141, 0),
//OpenMesh::Vec3d(-2.38072, 10.0441, 0),
//OpenMesh::Vec3d(-4.10297, 7.24421, 0),
//OpenMesh::Vec3d(-4.30619, 1.46692, 0),
//OpenMesh::Vec3d(-5.6163, -8.3802, 0),
//OpenMesh::Vec3d(-7.78844, -18.0231, 0),
//OpenMesh::Vec3d(2.65463, -13.4286, 0),
//OpenMesh::Vec3d(0.741755, -4.48237, 0),
//OpenMesh::Vec3d(2.73214, -4.76554, 0),
//OpenMesh::Vec3d(3.64475, -3.15402, 0),
//OpenMesh::Vec3d(8.02143, -2.49543, 0),
//OpenMesh::Vec3d(4.51096, -0.319015, 0),
//OpenMesh::Vec3d(4.58592, 2.2061, 0),
//OpenMesh::Vec3d(1.09773, 0.39865, 0),
//OpenMesh::Vec3d(5.22416, -1.9387, 0),
//		};

//CC_points = { //paris tower
//	OpenMesh::Vec3d(2.21312, -3.64502, 0),
//OpenMesh::Vec3d(0.0729413, -1.26159, 0),
//OpenMesh::Vec3d(0.71968, 2.04631, 0),
//OpenMesh::Vec3d(-0.00294265, 4.30579, 0),
//OpenMesh::Vec3d(-0.533977, 2.06356, 0),
//OpenMesh::Vec3d(0.020478, -1.26554, 0),
//OpenMesh::Vec3d(-2.13172, -3.62273, 0),
//OpenMesh::Vec3d(-1.5717, -3.65431, 0),
//OpenMesh::Vec3d(-0.789005, -3.66778, 0),
//OpenMesh::Vec3d(-0.017268, -3.66608, 0),
//OpenMesh::Vec3d(0.997374, -3.65183, 0),
//OpenMesh::Vec3d(1.73898, -3.63311, 0),
//};

//CC_points = { //xiyi
//	
//OpenMesh::Vec3d(4.87927, 0.893929, 0),
//OpenMesh::Vec3d(5.09543, -0.452398, 0),
//OpenMesh::Vec3d(2.58026, 0.369127, 0),
//OpenMesh::Vec3d(0.899758, 1.92868, 0),
//OpenMesh::Vec3d(-6.36326, 3.14212, 0),
//OpenMesh::Vec3d(-6.99906, 0.160185, 0),
//OpenMesh::Vec3d(-3.73234, -0.87026, 0),
//OpenMesh::Vec3d(-4.44834, -2.77133, 0),
//OpenMesh::Vec3d(-1.61845, -2.46282, 0),
//OpenMesh::Vec3d(1.3063, -2.28422, 0),
//OpenMesh::Vec3d(7.03627, -1.80929, 0),
//OpenMesh::Vec3d(5.87444, 2.51706, 0),
//OpenMesh::Vec3d(4.08634, 1.48513, 0),
//OpenMesh::Vec3d(3.23814, 0.190855, 0),
//OpenMesh::Vec3d(4.99338, 0.133408, 0),
//OpenMesh::Vec3d(4.78922, 0.868629, 0),
//OpenMesh::Vec3d(4.60864, 1.23347, 0),
//OpenMesh::Vec3d(4.20604, 1.08959, 0),
//OpenMesh::Vec3d(4.34064, 0.844649, 0),
//OpenMesh::Vec3d(4.51878, 1.00467, 0),
//OpenMesh::Vec3d(4.6284, 0.795768, 0),
//OpenMesh::Vec3d(4.43999, 0.739965, 0),
//OpenMesh::Vec3d(3.96868, 0.757685, 0),
//OpenMesh::Vec3d(4.47345, 1.60985, 0),
//};

//CC_points = { //yu
//	OpenMesh::Vec3d(1.55761, 0.298392, 0),
//OpenMesh::Vec3d(1.36297, 0.283486, 0),
//OpenMesh::Vec3d(1.1833, 0.280326, 0),
//OpenMesh::Vec3d(0.741254, 0.259984, 0),
//OpenMesh::Vec3d(0.29268, -2.05241, 0),
//OpenMesh::Vec3d(-1.99653, 0.0359225, 0),
//OpenMesh::Vec3d(-0.382603, 1.17903, 0),
//OpenMesh::Vec3d(-0.410516, 1.27791, 0),
//OpenMesh::Vec3d(-0.449854, 1.40897, 0),
//OpenMesh::Vec3d(-0.477108, 1.49262, 0),
//OpenMesh::Vec3d(-3.41672, -0.230562, 0),
//OpenMesh::Vec3d(1.02345, -3.44974, 0),
//
//};

//CC_points = {//kuzi
//	OpenMesh::Vec3d(0.677801, 1.2003, 0),
//OpenMesh::Vec3d(0.780243, -3.25194, 0),
//OpenMesh::Vec3d(0.81439, -9.3091, 0),
//OpenMesh::Vec3d(0.81439, -15.1614, 0),
//OpenMesh::Vec3d(3.80227, -15.1916, 0),
//OpenMesh::Vec3d(6.17551, -15.2559, 0),
//OpenMesh::Vec3d(11.5537, -15.2348, 0),
//OpenMesh::Vec3d(11.5672, -5.14217, 0),
//OpenMesh::Vec3d(11.6845, 9.43102, 0),
//OpenMesh::Vec3d(11.6459, 14.9842, 0),
//OpenMesh::Vec3d(-2.88771, 15.2544, 0),
//OpenMesh::Vec3d(-7.15668, 14.971, 0),
//OpenMesh::Vec3d(-11.383, 15.0709, 0),
//OpenMesh::Vec3d(-11.5878, -8.03373, 0),
//OpenMesh::Vec3d(-11.5878, -11.3671, 0),
//OpenMesh::Vec3d(-11.4171, -15.2809, 0),
//OpenMesh::Vec3d(-9.16339, -15.3532, 0),
//OpenMesh::Vec3d(-6.09014, -15.2718, 0),
//OpenMesh::Vec3d(-0.62658, -15.2758, 0),
//OpenMesh::Vec3d(-0.776266, -11.2293, 0),
//OpenMesh::Vec3d(-0.738142, -5.45842, 0),
//OpenMesh::Vec3d(-0.751239, 1.18322, 0),
//OpenMesh::Vec3d(-0.422864, 1.2003, 0),
//OpenMesh::Vec3d(-0.0432679, 1.2003, 0),
//};

//CC_points = {//car
//	OpenMesh::Vec3d(17.0636, -3.96522, 0),
//OpenMesh::Vec3d(17.0602, -0.865747, 0),
//OpenMesh::Vec3d(16.9603, 3.13988, 0),
//OpenMesh::Vec3d(13.1047, 2.52286, 0),
//OpenMesh::Vec3d(11.3294, 9.4102, 0),
//OpenMesh::Vec3d(0.0590485, 6.34927, 0),
//OpenMesh::Vec3d(-2.5501, 2.84747, 0),
//OpenMesh::Vec3d(-12.414, 2.59491, 0),
//OpenMesh::Vec3d(-16.5235, 1.10555, 0),
//OpenMesh::Vec3d(-17.0582, -3.72901, 0),
//OpenMesh::Vec3d(-4.71722, -8.88474, 0),
//OpenMesh::Vec3d(6.7265, -7.13234, 0),
//};

//CC_points = {//tuzi
//	OpenMesh::Vec3d(1.2338, -0.872146, 0),
//OpenMesh::Vec3d(0.848741, -0.488863, 0),
//OpenMesh::Vec3d(0.900976, -0.193908, 0),
//OpenMesh::Vec3d(0.789363, 0.186805, 0),
//OpenMesh::Vec3d(1.20989, 0.53253, 0),
//OpenMesh::Vec3d(1.40575, 1.25384, 0),
//OpenMesh::Vec3d(0.776231, 1.62091, 0),
//OpenMesh::Vec3d(-0.294722, 1.55362, 0),
//OpenMesh::Vec3d(0.952108, 0.857486, 0),
//OpenMesh::Vec3d(0.145927, 0.50355, 0),
//OpenMesh::Vec3d(-0.335255, 0.525864, 0),
//OpenMesh::Vec3d(-0.896345, 0.022467, 0),
//OpenMesh::Vec3d(-0.979956, -0.381595, 0),
//OpenMesh::Vec3d(-1.03288, -0.89268, 0),
//OpenMesh::Vec3d(-1.25456, -0.83442, 0),
//OpenMesh::Vec3d(-1.21265, -1.61744, 0),
//OpenMesh::Vec3d(-0.690139, -1.62269, 0),
//OpenMesh::Vec3d(-0.28394, -1.64832, 0),
//OpenMesh::Vec3d(0.223731, -1.65626, 0),
//OpenMesh::Vec3d(0.416271, -1.6679, 0),
//OpenMesh::Vec3d(0.533081, -1.66271, 0),
//OpenMesh::Vec3d(0.753259, -1.66561, 0),
//OpenMesh::Vec3d(0.811291, -1.37615, 0),
//OpenMesh::Vec3d(1.19484, -1.34047, 0),
//};

//CC_points = { //daizi
//	OpenMesh::Vec3d(31.4678, -3.60352, 0),
//OpenMesh::Vec3d(31.4678, -0.717967, 0),
//OpenMesh::Vec3d(31.4678, 0.615366, 0),
//OpenMesh::Vec3d(31.4549, 3.54117, 0),
//OpenMesh::Vec3d(0.653732, 3.59247, 0),
//OpenMesh::Vec3d(-0.679602, 3.59247, 0),
//OpenMesh::Vec3d(-31.6356, 3.56382, 0),
//OpenMesh::Vec3d(-31.6226, 0.638025, 0),
//OpenMesh::Vec3d(-31.6226, -0.695309, 0),
//OpenMesh::Vec3d(-31.6226, -3.58087, 0),
//OpenMesh::Vec3d(-0.666667, -3.55222, 0),
//OpenMesh::Vec3d(0.666667, -3.55222, 0),
//};
//CC_points = { //fish
//	OpenMesh::Vec3d(7.88917, -4.19279, 0),
//OpenMesh::Vec3d(7.47304, -1.84145, 0),
//OpenMesh::Vec3d(7.60398, -0.916129, 0),
//OpenMesh::Vec3d(10.1018, -0.0975394, 0),
//OpenMesh::Vec3d(8.99175, 1.64069, 0),
//OpenMesh::Vec3d(4.941, 0.121658, 0),
//OpenMesh::Vec3d(3.34918, 1.5249, 0),
//OpenMesh::Vec3d(3.69759, 2.6337, 0),
//OpenMesh::Vec3d(2.35008, 3.81885, 0),
//OpenMesh::Vec3d(0.921122, 3.08854, 0),
//OpenMesh::Vec3d(-1.81225, 4.09559, 0),
//OpenMesh::Vec3d(-5.65018, 4.77371, 0),
//OpenMesh::Vec3d(-8.65683, 3.24567, 0),
//OpenMesh::Vec3d(-10.3053, 3.09553, 0),
//OpenMesh::Vec3d(-10.0761, 2.09196, 0),
//OpenMesh::Vec3d(-9.00522, 1.38062, 0),
//OpenMesh::Vec3d(-8.98836, 1.07329, 0),
//OpenMesh::Vec3d(-8.99063, 0.702774, 0),
//OpenMesh::Vec3d(-8.75359, 0.482605, 0),
//OpenMesh::Vec3d(-8.63689, 0.763358, 0),
//OpenMesh::Vec3d(-8.78799, 1.02137, 0),
//OpenMesh::Vec3d(-8.72299, 1.21815, 0),
//OpenMesh::Vec3d(-7.24225, -0.0052856, 0),
//OpenMesh::Vec3d(-4.19355, -1.06031, 0),
//OpenMesh::Vec3d(-2.08033, -0.944584, 0),
//OpenMesh::Vec3d(-1.63804, -1.50076, 0),
//OpenMesh::Vec3d(-0.235829, -2.02419, 0),
//OpenMesh::Vec3d(-0.37616, -0.737527, 0),
//OpenMesh::Vec3d(-0.207549, -0.699108, 0),
//OpenMesh::Vec3d(0.0779521, -0.665316, 0),
//OpenMesh::Vec3d(0.427525, -0.621636, 0),
//OpenMesh::Vec3d(1.40535, -1.90493, 0),
//OpenMesh::Vec3d(2.01614, -2.07862, 0),
//OpenMesh::Vec3d(3.86555, -0.802197, 0),
//OpenMesh::Vec3d(5.24685, -1.85539, 0),
//OpenMesh::Vec3d(6.49604, -3.83112, 0),
//};

//CC_points = { //choufish
//	OpenMesh::Vec3d(0.621947, -3.00053, 0),
//	OpenMesh::Vec3d(1.19422, -3.24356, 0),
//	OpenMesh::Vec3d(1.76649, -3.4866, 0),
//	OpenMesh::Vec3d(2.33876, -3.72963, 0),
//	OpenMesh::Vec3d(3.48501, -2.98653, 0),
//	OpenMesh::Vec3d(4.63127, -2.24344, 0),
//	OpenMesh::Vec3d(5.77752, -1.50034, 0),
//	OpenMesh::Vec3d(4.8046, -0.0228433, 0),
//	OpenMesh::Vec3d(3.83168, 1.45465, 0),
//	OpenMesh::Vec3d(2.85876, 2.93215, 0),
//	OpenMesh::Vec3d(0.10979, 3.08279, 0),
//	OpenMesh::Vec3d(-2.63918, 3.23342, 0),
//	OpenMesh::Vec3d(-5.38815, 3.38406, 0),
//	OpenMesh::Vec3d(-5.43681, 3.13371, 0),
//	OpenMesh::Vec3d(-5.48547, 2.88337, 0),
//	OpenMesh::Vec3d(-5.53413, 2.63302, 0),
//	OpenMesh::Vec3d(-3.58597, 2.51262, 0),
//	OpenMesh::Vec3d(-1.63781, 2.39221, 0),
//	OpenMesh::Vec3d(0.310345, 2.27181, 0),
//	OpenMesh::Vec3d(-0.884167, 1.77491, 0),
//	OpenMesh::Vec3d(-2.07868, 1.27801, 0),
//	OpenMesh::Vec3d(-3.27319, 0.781111, 0),
//	OpenMesh::Vec3d(-3.36141, -0.046006, 0),
//	OpenMesh::Vec3d(-3.44963, -0.873123, 0),
//	OpenMesh::Vec3d(-3.53785, -1.70024, 0),
//	OpenMesh::Vec3d(-3.00611, -2.34053, 0),
//	OpenMesh::Vec3d(-2.47437, -2.98081, 0),
//	OpenMesh::Vec3d(-1.94263, -3.6211, 0),
//	OpenMesh::Vec3d(-1.12739, -3.40618, 0),
//	OpenMesh::Vec3d(-0.312145, -3.19127, 0),
//	OpenMesh::Vec3d(0.503098, -2.97635, 0),
//	OpenMesh::Vec3d(0.930805, -2.94226, 0),
//	OpenMesh::Vec3d(1.35851, -2.90816, 0),
//	OpenMesh::Vec3d(1.78622, -2.87407, 0),
//	OpenMesh::Vec3d(1.39813, -2.91622, 0),
//	OpenMesh::Vec3d(1.01004, -2.95838, 0),
//};

//CC_points = { //circle
//	OpenMesh::Vec3d(3.45289, 1.63941, 0),
//OpenMesh::Vec3d(4.48355, 0.137393, 0),
//OpenMesh::Vec3d(4.01564, -1.4928, 0),
//OpenMesh::Vec3d(3.44068, -3.85162, 0),
//OpenMesh::Vec3d(6.55334, -0.345064, 0),
//OpenMesh::Vec3d(3.90656, 6.20451, 0),
//OpenMesh::Vec3d(-1.69168, 3.85054, 0),
//OpenMesh::Vec3d(-0.0258138, 5.08684, 0),
//OpenMesh::Vec3d(1.91392, 4.70734, 0),
//OpenMesh::Vec3d(3.83488, 3.89295, 0),
//OpenMesh::Vec3d(1.07448, 7.12438, 0),
//OpenMesh::Vec3d(-5.90614, 4.59064, 0),
//OpenMesh::Vec3d(-3.65812, -1.22062, 0),
//OpenMesh::Vec3d(-4.84896, 1.04082, 0),
//OpenMesh::Vec3d(-4.17104, 2.38932, 0),
//OpenMesh::Vec3d(-3.86462, 4.22689, 0),
//OpenMesh::Vec3d(-6.53473, 1.87253, 0),
//OpenMesh::Vec3d(-4.64774, -5.87751, 0),
//OpenMesh::Vec3d(1.34646, -3.2845, 0),
//OpenMesh::Vec3d(-0.105171, -4.56524, 0),
//OpenMesh::Vec3d(-2.44788, -4.06556, 0),
//OpenMesh::Vec3d(-4.21109, -3.22378, 0),
//OpenMesh::Vec3d(-1.36415, -6.68839, 0),
//OpenMesh::Vec3d(5.46868, -4.0388, 0),
//};

//CC_points = { //flower
//	OpenMesh::Vec3d(0.060138, -2.34968, 0),
//OpenMesh::Vec3d(0.167533, -1.16113, 0),
//OpenMesh::Vec3d(0.304982, 0.160104, 0),
//OpenMesh::Vec3d(0.231798, 1.08177, 0),
//OpenMesh::Vec3d(0.61252, 2.32314, 0),
//OpenMesh::Vec3d(0.771864, 3.13087, 0),
//OpenMesh::Vec3d(0.855257, 4.5897, 0),
//OpenMesh::Vec3d(2.09554, 4.91859, 0),
//OpenMesh::Vec3d(2.47687, 5.80984, 0),
//OpenMesh::Vec3d(2.67184, 7.06251, 0),
//OpenMesh::Vec3d(1.05076, 8.34497, 0),
//OpenMesh::Vec3d(-1.00333, 8.56513, 0),
//OpenMesh::Vec3d(-3.01162, 7.65609, 0),
//OpenMesh::Vec3d(-2.46392, 6.1255, 0),
//OpenMesh::Vec3d(-1.33722, 5.03078, 0),
//OpenMesh::Vec3d(-0.552003, 4.00657, 0),
//OpenMesh::Vec3d(-0.0860061, 2.21267, 0),
//OpenMesh::Vec3d(0.0878497, 0.616274, 0),
//OpenMesh::Vec3d(-0.0636332, -0.969635, 0),
//OpenMesh::Vec3d(-0.469408, 0.4811, 0),
//OpenMesh::Vec3d(-1.36141, 0.898247, 0),
//OpenMesh::Vec3d(-2.47641, 1.32553, 0),
//OpenMesh::Vec3d(-2.85447, -0.680847, 0),
//OpenMesh::Vec3d(-1.57427, -0.75945, 0),
//OpenMesh::Vec3d(-0.29709, -1.64989, 0),
//OpenMesh::Vec3d(-0.284548, -4.19343, 0),
//OpenMesh::Vec3d(-0.449725, -7.19035, 0),
//OpenMesh::Vec3d(-0.0458234, -9.10946, 0),
//OpenMesh::Vec3d(0.607011, -6.40262, 0),
//OpenMesh::Vec3d(-0.378727, -4.82229, 0),
//OpenMesh::Vec3d(0.125836, -2.84884, 0),
//OpenMesh::Vec3d(0.941884, -1.5343, 0),
//OpenMesh::Vec3d(2.07567, -2.26778, 0),
//OpenMesh::Vec3d(2.4943, -2.94274, 0),
//OpenMesh::Vec3d(4.17284, -0.459165, 0),
//OpenMesh::Vec3d(0.512756, -0.0850323, 0),
//};

//CC_points = { //octopus
// OpenMesh::Vec3d(8.57042, 3.297, 0),
//OpenMesh::Vec3d(5.09424, 1.93482, 0),
//OpenMesh::Vec3d(8.81644, -0.913078, 0),
//OpenMesh::Vec3d(6.50689, -3.68536, 0),
//OpenMesh::Vec3d(5.42654, -1.87831, 0),
//OpenMesh::Vec3d(4.92811, 1.19895, 0),
//OpenMesh::Vec3d(5.93302, 3.06514, 0),
//OpenMesh::Vec3d(6.39359, 2.98059, 0),
//OpenMesh::Vec3d(6.76971, 2.94039, 0),
//OpenMesh::Vec3d(7.1323, 3.39599, 0),
//OpenMesh::Vec3d(7.34691, 3.96736, 0),
//OpenMesh::Vec3d(7.36178, 4.91438, 0),
//OpenMesh::Vec3d(6.49071, 5.0889, 0),
//OpenMesh::Vec3d(7.77231, 9.87716, 0),
//OpenMesh::Vec3d(5.30888, 14.7788, 0),
//OpenMesh::Vec3d(-0.841746, 14.3125, 0),
//OpenMesh::Vec3d(-5.90767, 11.3197, 0),
//OpenMesh::Vec3d(-8.86451, 3.55846, 0),
//OpenMesh::Vec3d(-0.246632, 4.78125, 0),
//OpenMesh::Vec3d(1.06736, 1.00327, 0),
//OpenMesh::Vec3d(-4.23325, 0.219585, 0),
//OpenMesh::Vec3d(-8.03561, 3.02382, 0),
//OpenMesh::Vec3d(-8.08751, 2.08127, 0),
//OpenMesh::Vec3d(-7.52196, 0.833628, 0),
//OpenMesh::Vec3d(-6.06773, 0.106026, 0),
//OpenMesh::Vec3d(-7.90343, -1.09784, 0),
//OpenMesh::Vec3d(-9.28057, -2.77283, 0),
//OpenMesh::Vec3d(-7.6208, -5.78163, 0),
//OpenMesh::Vec3d(-8.57419, -6.68359, 0),
//OpenMesh::Vec3d(-8.68614, -8.01901, 0),
//OpenMesh::Vec3d(-8.06132, -10.9547, 0),
//OpenMesh::Vec3d(-6.74283, -6.93736, 0),
//OpenMesh::Vec3d(-4.59648, -5.33218, 0),
//OpenMesh::Vec3d(-2.05936, -4.13197, 0),
//OpenMesh::Vec3d(-3.09482, -7.04246, 0),
//OpenMesh::Vec3d(-9.21992, -12.7162, 0),
//OpenMesh::Vec3d(-2.27186, -13.953, 0),
//OpenMesh::Vec3d(-2.95275, -12.5256, 0),
//OpenMesh::Vec3d(-3.79382, -11.3406, 0),
//OpenMesh::Vec3d(-3.42902, -9.89676, 0),
//OpenMesh::Vec3d(-2.6323, -11.6362, 0),
//OpenMesh::Vec3d(-1.81066, -12.1583, 0),
//OpenMesh::Vec3d(-0.00107587, -12.4131, 0),
//OpenMesh::Vec3d(-2.69772, -8.66487, 0),
//OpenMesh::Vec3d(-1.78497, -7.37903, 0),
//OpenMesh::Vec3d(0.488619, -4.78304, 0),
//OpenMesh::Vec3d(-0.0659694, -6.98257, 0),
//OpenMesh::Vec3d(-0.438017, -8.1058, 0),
//OpenMesh::Vec3d(0.945469, -9.78723, 0),
//OpenMesh::Vec3d(1.71534, -8.54231, 0),
//OpenMesh::Vec3d(1.34671, -7.00428, 0),
//OpenMesh::Vec3d(1.89121, -5.32397, 0),
//OpenMesh::Vec3d(3.61222, -8.87285, 0),
//OpenMesh::Vec3d(2.05656, -8.95106, 0),
//OpenMesh::Vec3d(1.51134, -12.1472, 0),
//OpenMesh::Vec3d(7.13949, -10.1638, 0),
//OpenMesh::Vec3d(2.57694, -5.41496, 0),
//OpenMesh::Vec3d(4.12527, -0.841089, 0),
//OpenMesh::Vec3d(4.43461, -3.18152, 0),
//OpenMesh::Vec3d(5.09308, -5.72131, 0),
//OpenMesh::Vec3d(7.13664, -5.11346, 0),
//OpenMesh::Vec3d(10.4172, -2.70598, 0),
//OpenMesh::Vec3d(6.50459, 2.41615, 0),
//OpenMesh::Vec3d(9.27712, 2.26177, 0),
//OpenMesh::Vec3d(9.67284, 2.20834, 0),
//OpenMesh::Vec3d(9.37983, 3.08828, 0),
//};

//CC_points = { //sheji
//	OpenMesh::Vec3d(-13.903, -6.69777, 0),
//OpenMesh::Vec3d(1.05228, -8.39994, 0),
//OpenMesh::Vec3d(7.48085, -1.03554, 0),
//OpenMesh::Vec3d(5.9094, 9.53619, 0),
//OpenMesh::Vec3d(2.46234, 8.10173, 0),
//OpenMesh::Vec3d(-1.18305, 8.99174, 0),
//OpenMesh::Vec3d(-4.21192, 10.1169, 0),
//OpenMesh::Vec3d(-4.81756, 4.98018, 0),
//OpenMesh::Vec3d(-7.6985, 4.82465, 0),
//OpenMesh::Vec3d(-12.3981, 4.54636, 0),
//OpenMesh::Vec3d(-12.0119, 0.141044, 0),
//OpenMesh::Vec3d(-12.2488, -2.39566, 0),
//};
//CC_points = {//pangxie
//	OpenMesh::Vec3d(7.71977, -2.05287, 0),
//OpenMesh::Vec3d(16.1661, 2.5924, 0),
//OpenMesh::Vec3d(22.5491, -0.870733, 0),
//OpenMesh::Vec3d(15.6577, -8.58216, 0),
//OpenMesh::Vec3d(26.8967, -7.59486, 0),
//OpenMesh::Vec3d(21.084, 11.3199, 0),
//OpenMesh::Vec3d(7.68061, 1.27629, 0),
//OpenMesh::Vec3d(11.5959, 6.63227, 0),
//OpenMesh::Vec3d(24.6401, 7.52288, 0),
//OpenMesh::Vec3d(21.7966, -5.4645, 0),
//OpenMesh::Vec3d(31.1961, 2.17972, 0),
//OpenMesh::Vec3d(15.6549, 13.7466, 0),
//OpenMesh::Vec3d(7.82502, 4.58309, 0),
//OpenMesh::Vec3d(11.8374, 11.0923, 0),
//OpenMesh::Vec3d(20.9685, 10.4649, 0),
//OpenMesh::Vec3d(24.9294, 1.87976, 0),
//OpenMesh::Vec3d(27.8247, 9.55625, 0),
//OpenMesh::Vec3d(14.2738, 16.6223, 0),
//OpenMesh::Vec3d(5.4765, 7.35184, 0),
//OpenMesh::Vec3d(8.64119, 11.865, 0),
//OpenMesh::Vec3d(14.7017, 15.1746, 0),
//OpenMesh::Vec3d(23.1277, 9.7121, 0),
//OpenMesh::Vec3d(19.424, 19.1598, 0),
//OpenMesh::Vec3d(6.83002, 15.4013, 0),
//OpenMesh::Vec3d(2.16398, 7.18829, 0),
//OpenMesh::Vec3d(-0.0248164, 9.29763, 0),
//OpenMesh::Vec3d(-1.17408, 9.21065, 0),
//OpenMesh::Vec3d(-3.33054, 6.62634, 0),
//OpenMesh::Vec3d(-8.85175, 14.7474, 0),
//OpenMesh::Vec3d(-24.3841, 17.0936, 0),
//OpenMesh::Vec3d(-23.6175, 7.39743, 0),
//OpenMesh::Vec3d(-18.8441, 13.6221, 0),
//OpenMesh::Vec3d(-10.5837, 11.5502, 0),
//OpenMesh::Vec3d(-5.13643, 5.07839, 0),
//OpenMesh::Vec3d(-5.27052, 4.97835, 0),
//OpenMesh::Vec3d(-5.46014, 4.84129, 0),
//OpenMesh::Vec3d(-5.64977, 4.70424, 0),
//OpenMesh::Vec3d(-12.9059, 15.3053, 0),
//OpenMesh::Vec3d(-28.4601, 9.29942, 0),
//OpenMesh::Vec3d(-24.863, 0.0668226, 0),
//OpenMesh::Vec3d(-22.7747, 8.29837, 0),
//OpenMesh::Vec3d(-14.3986, 9.84779, 0),
//OpenMesh::Vec3d(-7.0058, 3.14714, 0),
//OpenMesh::Vec3d(-18.5345, 12.5865, 0),
//OpenMesh::Vec3d(-30.6453, -0.223078, 0),
//OpenMesh::Vec3d(-21.5541, -7.27411, 0),
//OpenMesh::Vec3d(-24.3589, 0.197288, 0),
//OpenMesh::Vec3d(-18.8381, 7.34974, 0),
//OpenMesh::Vec3d(-7.4344, 0.51362, 0),
//OpenMesh::Vec3d(-24.585, 8.96755, 0),
//OpenMesh::Vec3d(-25.4714, -12.1974, 0),
//OpenMesh::Vec3d(-13.6432, -10.1693, 0),
//OpenMesh::Vec3d(-22.614, -4.87779, 0),
//OpenMesh::Vec3d(-17.811, 3.67923, 0),
//OpenMesh::Vec3d(-7.53781, -3.00233, 0),
//OpenMesh::Vec3d(-23.8508, 3.73205, 0),
//OpenMesh::Vec3d(-18.0268, -14.6281, 0),
//OpenMesh::Vec3d(0.581672, -14.6728, 0),
//OpenMesh::Vec3d(-1.23569, -8.99517, 0),
//OpenMesh::Vec3d(-10.3154, -6.40558, 0),
//OpenMesh::Vec3d(-14.2051, -5.46234, 0),
//OpenMesh::Vec3d(-10.7249, -5.26371, 0),
//OpenMesh::Vec3d(-4.0724, -6.50565, 0),
//OpenMesh::Vec3d(0.270762, -12.5408, 0),
//OpenMesh::Vec3d(6.57959, -7.9814, 0),
//OpenMesh::Vec3d(9.1306, -3.60342, 0),
//OpenMesh::Vec3d(13.4053, -3.81022, 0),
//OpenMesh::Vec3d(11.9635, -6.90228, 0),
//OpenMesh::Vec3d(4.07802, -8.25055, 0),
//OpenMesh::Vec3d(0.654098, -13.3543, 0),
//OpenMesh::Vec3d(3.50707, -16.7318, 0),
//OpenMesh::Vec3d(8.3762, -13.8882, 0),
//OpenMesh::Vec3d(12.9466, -9.91377, 0),
//OpenMesh::Vec3d(19.5164, -2.8494, 0),
//OpenMesh::Vec3d(19.745, 1.90087, 0),
//};
		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_green_weight323();//计算3次cage到3次控制点的权重
	}

	Set_Texture_coord();//设置纹理
}


void MeshViewerWidget::CubicMVC_Test(void)
{
	std::cout << "CubicMVC!" << std::endl;
	drawmode = CURVECAGE;
	usecvm = true;
	auto arthono = [](const OpenMesh::Vec3d& p) -> OpenMesh::Vec3d {//left( -y, x)
		return OpenMesh::Vec3d(p[1],- p[0], p[2]);
	};
	//正儿八经的MVC
	if (1) {
		int N;
		std::vector<OpenMesh::Vec3d> polygon_vertices;
		//四边形
		if (0) {
			N = 4;
			curvecage2.resize(N);
			
			for (int v_it = 0; v_it <  mesh.n_vertices(); v_it++) {
				auto vh = mesh.vertex_handle(v_it);
				auto point = mesh.point(vh);
				Mesh::Point pointuni;
				pointuni[0] = (point[0] - ptMin[0]) / (ptMax[0] - ptMin[0]); // Normalize x
				pointuni[1] =  (point[1] - ptMin[1]) / (ptMax[1] - ptMin[1]); // Normalize y
				pointuni[0] *= (1 + pointuni[1]);
				pointuni[2] = 0;
				//std::cout << pointuni[0]<<" " << pointuni[1] << std::endl;
				mesh.set_point(vh, pointuni);
			}
			// 逆时针定义8个顶点
			polygon_vertices = {
				OpenMesh::Vec3d(2,1,0),
				OpenMesh::Vec3d(0,1,0),
				OpenMesh::Vec3d(0,0,0),
				OpenMesh::Vec3d(1,0,0)
			};
		}
		//choufish
		if (0) {
			N = 12;
			curvecage2.resize(N);
			polygon_vertices = {
				OpenMesh::Vec3d(0.621947, -3.00053, 0),
                OpenMesh::Vec3d(2.33876, -3.72963, 0),
                OpenMesh::Vec3d(5.77752, -1.50034, 0),
                OpenMesh::Vec3d(2.85876, 2.93215, 0),
                OpenMesh::Vec3d(-5.38815, 3.38406, 0),
                OpenMesh::Vec3d(-5.53413, 2.63302, 0),
                OpenMesh::Vec3d(0.310345, 2.27181, 0),
                OpenMesh::Vec3d(-3.27319, 0.781111, 0),
                OpenMesh::Vec3d(-3.53785, -1.70024, 0),
                OpenMesh::Vec3d(-1.94263, -3.6211, 0),
                OpenMesh::Vec3d(0.503098, -2.97635, 0),
                OpenMesh::Vec3d(1.78622, -2.87407, 0),
			};
			for (int v_it = 0; v_it < mesh.n_vertices(); v_it++) {
				auto vh = mesh.vertex_handle(v_it);
				auto point = mesh.point(vh);
				Mesh::Point pointuni;
				pointuni[0] = (point[0]); // Normalize x
				pointuni[1] = (point[1]); // Normalize y

				pointuni[2] = point[0] / 10000;
				//std::cout << pointuni[0]<<" " << pointuni[1] << std::endl;
				mesh.set_point(vh, pointuni);
			}
		}
		//kuzi
		if (0) {
			N = 8;//kuzi 8
			curvecage2.resize(N);
			// 逆时针定义8个顶点
			polygon_vertices = {
				OpenMesh::Vec3d(11.4393 + 0.01, 14.8404 + 0.01, 0),
				OpenMesh::Vec3d(-11.3137 - 0.01, 14.8404 + 0.01, 0),
				OpenMesh::Vec3d(-11.3137 - 0.01, -15.0839 - 0.01, 0),
				OpenMesh::Vec3d(-0.822533 + 0.01, -15.0839 - 0.01, 0),
				OpenMesh::Vec3d(-0.822533 + 0.01, 1.35013 - 0.01, 0),
				OpenMesh::Vec3d(0.937067 - 0.01, 1.36120 - 0.01, 0),
				OpenMesh::Vec3d(0.948133 - 0.01, -15.0839 - 0.01, 0),
				OpenMesh::Vec3d(11.4393 + 0.01, -15.0839 - 0.01, 0)
			};
			for (int v_it = 0; v_it < mesh.n_vertices(); v_it++) {
				auto vh = mesh.vertex_handle(v_it);
				auto point = mesh.point(vh);
				Mesh::Point pointuni;
				pointuni[0] = (point[0]); // Normalize x
				pointuni[1] = (point[1]); // Normalize y
				
				pointuni[2] = point[0]/10000;
				//std::cout << pointuni[0]<<" " << pointuni[1] << std::endl;
				mesh.set_point(vh, pointuni);
			}

		}
		//zoomfish
		if (0) {
			N = 4;//zoom fish 4
			curvecage2.resize(N);
			
			// 逆时针定义4个顶点
			polygon_vertices = {
				OpenMesh::Vec3d(15.9365 + 0.01,  1.25963 + 0.01, 0),
				OpenMesh::Vec3d(-16.3751 - 0.01,  1.25963 + 0.01, 0),
				OpenMesh::Vec3d(-16.3751 - 0.01, -1.86695 - 0.01, 0),
				OpenMesh::Vec3d(15.9365 + 0.01, -1.86695 - 0.01, 0),

			};
			

		}
		//sheji
		if (1) {
			N = 5;//circle 4
			curvecage2.resize(N);

			// 逆时针定义3个顶点
			polygon_vertices = {
				OpenMesh::Vec3d(-12.729, -6.31988, 0),
					OpenMesh::Vec3d(5.51534, -6.20705, 0),
					OpenMesh::Vec3d(6.0044, 9.48202, 0),
					OpenMesh::Vec3d(-4.21192, 10.1169, 0),
					OpenMesh::Vec3d(-12.3981, 4.54636, 0),
			};
		}


		for (int i = 0; i < N; i++) {
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}

		std::vector<OpenMesh::Vec3d> t_e(N);
		std::vector<OpenMesh::Vec3d> n_e(N);
		mvcGn.resize(2 * N);
		mvcGt.resize(2 * N);
		mvcL.resize(N);
		
		std::vector<double> vc(N);
		std::vector<double> gnc(2 * N);
		std::vector<double> gtc(2 * N);
		OpenMesh::Vec3d ptest = OpenMesh::Vec3d(1.1, 1, 0);
		MVC::cubicMVCs(polygon_vertices, ptest, vc, gnc, gtc);
		OpenMesh::Vec3d r_test = OpenMesh::Vec3d(0, 0, 0);

		std::vector<double> L(N);
		for (int i = 0; i < N; i++)
		{
			mvcGt[2 * i] = 3*(curvecage2[i][1] - curvecage2[i][0])/norm(curvecage2[i][3] - curvecage2[i][0]);
			mvcGt[2 * i + 1] = 3*(curvecage2[i][2] - curvecage2[i][3]) / norm(curvecage2[i][3] - curvecage2[i][0]);
			mvcL[i] = (curvecage2[i][3] - curvecage2[i][0]).norm();//每条边的长度
			n_e[i] = arthono((curvecage2[i][3] - curvecage2[i][0])).normalized();//每条边的单位法向
			t_e[i] = (curvecage2[i][3] - curvecage2[i][0]).normalized();//每条边的单位切向
		}
		auto f = polygon_vertices;
		auto df = mvcGt;
		auto h = mvcGt;
		//求解h
		for (int i = 0; i < N; i++)
		{
			int index = (i - 1) % N;
			if (index < 0) {
				index += N;
			}
			Eigen::Matrix2d A;
			Eigen::Matrix2d B;
			OpenMesh::Vec3d Bx = -mvcGt[2 * index + 1][0] * t_e[index] - mvcGt[2 * i][0] * t_e[i];
			OpenMesh::Vec3d By = -mvcGt[2 * index + 1][1] * t_e[index] - mvcGt[2 * i][1] * t_e[i];
			// 为矩阵A和B赋值
			A << n_e[i][0], -n_e[index][0],
				n_e[i][1], -n_e[index][1];
			B << Bx[0], By[0],
				Bx[1], By[1];
			//solve gn
			Eigen::Matrix2d X = A.inverse() * B;
			//std::cout << X << std::endl;
			mvcGn[2 * index + 1] = OpenMesh::Vec3d(X(1, 0), X(1, 1), 0);
			mvcGn[2 * i] = OpenMesh::Vec3d(X(0, 0), X(0, 1), 0);
			/*std::cout <<"x的梯度" << gt[2 * i][0] * t_e[i] + gn[2 * i][0] * n_e[i] << std::endl;
			std::cout << "y的梯度" << gt[2 * i][1] * t_e[i] + gn[2 * i][1] * n_e[i] << std::endl;
			std::cout << gt[2 * i][0] * t_e[i] + gn[2 * i][0] * n_e[i] - (-gt[2 * index + 1][0] * t_e[index] + gn[2 * index + 1][0] * n_e[index]) << std::endl;
			std::cout << gt[2 * i][1] * t_e[i] + gn[2 * i][1] * n_e[i] - (-gt[2 * index + 1][1] * t_e[index] + gn[2 * index + 1][1] * n_e[index]) << std::endl;*/

		}
		
		for (int i = 0; i < N; i++)
		{
			r_test += vc[i] * polygon_vertices[i];
			r_test += gtc[2 * i] * mvcGt[2 * i];
			r_test += gtc[2 * i + 1] * mvcGt[2 * i + 1];
			r_test += gnc[2 * i] * mvcGn[2 * i];
			r_test += gnc[2 * i + 1] * mvcGn[2 * i + 1];
		}
		std::cout << r_test - ptest << std::endl;
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}

		//calculate_green_weight123();//计算3次cage到3次控制点的权重
		weights.resize(mesh.n_vertices(), std::vector<double>(5 * curvecage2.size(), 0.0));
		for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {
			N = curvecage2.size();
			std::vector<double> vc(N);
			std::vector<double> gnc(2 * N);
			std::vector<double> gtc(2 * N);
			auto vh = mesh.vertex_handle(v_id);
			Mesh::Point eta = mesh.point(vh);
			MVC::cubicMVCs(polygon_vertices, eta, vc, gnc, gtc);
			for (int i = 0; i < curvecage2.size(); i++)
			{
				weights[v_id][5 * i] = vc[i];
				weights[v_id][5 * i + 1] = gtc[2 * i];
				weights[v_id][5 * i + 2] = gtc[2 * i + 1];
				weights[v_id][5 * i + 3] = gnc[2 * i];
				weights[v_id][5 * i + 4] = gnc[2 * i + 1];
			}

		}
		Set_Texture_coord();//设置纹理

	}
}


//仅针对kuzi测试
void MeshViewerWidget::PolyGC_Test(void)
{
	drawmode = CURVECAGE;
	double test = F2_n(Mesh::Point(0, 0, 0), Mesh::Point(0, 2, 0), Mesh::Point(-2, -2, 0), Mesh::Point(1, 0, 0), 1);

	//正儿八经的GC
	if (0) {
		int N = 8;//kuzi 8
		curvecage2.resize(N);

		// 逆时针定义8个顶点
		std::vector<OpenMesh::Vec3d> square_vertices = {
			OpenMesh::Vec3d(11.4393 + 0.01, 14.8404 + 0.01, 0),
			OpenMesh::Vec3d(-11.3137 - 0.01, 14.8404 + 0.01, 0),
			OpenMesh::Vec3d(-11.3137 - 0.01, -15.0839 - 0.01, 0),
			OpenMesh::Vec3d(-0.822533 + 0.01, -15.0839 - 0.01, 0),
			OpenMesh::Vec3d(-0.822533 + 0.01, 1.35013 - 0.01, 0),
			OpenMesh::Vec3d(0.937067 - 0.01, 1.36120 - 0.01, 0),
			OpenMesh::Vec3d(0.948133 - 0.01, -15.0839 - 0.01, 0),
			OpenMesh::Vec3d(11.4393 + 0.01, -15.0839 - 0.01, 0)
		};
		for (int i = 0; i < 8; i++) {
			OpenMesh::Vec3d start = square_vertices[i];
			OpenMesh::Vec3d end = square_vertices[(i + 1) % 8];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_green_weight123();//计算3次cage到3次控制点的权重
	}
	//choufish的GC
	if (0) {
		int N = 12;
		curvecage2.resize(N);
		std::vector<OpenMesh::Vec3d> polygon_vertices = {
			OpenMesh::Vec3d(0.621947, -3.00053, 0),
			OpenMesh::Vec3d(2.33876, -3.72963, 0),
			OpenMesh::Vec3d(5.77752, -1.50034, 0),
			OpenMesh::Vec3d(2.85876, 2.93215, 0),
			OpenMesh::Vec3d(-5.38815, 3.38406, 0),
			OpenMesh::Vec3d(-5.53413, 2.63302, 0),
			OpenMesh::Vec3d(0.310345, 2.27181, 0),
			OpenMesh::Vec3d(-3.27319, 0.781111, 0),
			OpenMesh::Vec3d(-3.53785, -1.70024, 0),
			OpenMesh::Vec3d(-1.94263, -3.6211, 0),
			OpenMesh::Vec3d(0.503098, -2.97635, 0),
			OpenMesh::Vec3d(1.78622, -2.87407, 0),
		};
		for (int i = 0; i < N; i++) {
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_green_weight123();//计算3次cage到3次控制点的权重
	}
	//第一条边特殊的GCkuzi
	if (0)
	{
		int N = 8;//kuzi 8
		curvecage2.resize(N);

		// 逆时针定义8个顶点
		std::vector<OpenMesh::Vec3d> square_vertices = {
			OpenMesh::Vec3d(11.4393 + 0.01, 14.8404 + 0.01, 0),
			OpenMesh::Vec3d(-11.3137 - 0.01, 14.8404 + 0.01, 0),
			OpenMesh::Vec3d(-11.3137 - 0.01, -15.0839 - 0.01, 0),
			OpenMesh::Vec3d(-0.822533 + 0.01, -15.0839 - 0.01, 0),
			OpenMesh::Vec3d(-0.822533 + 0.01, 1.35013 - 0.01, 0),
			OpenMesh::Vec3d(0.937067 - 0.01, 1.36120 - 0.01, 0),
			OpenMesh::Vec3d(0.948133 - 0.01, -15.0839 - 0.01, 0),
			OpenMesh::Vec3d(11.4393 + 0.01, -15.0839 - 0.01, 0)
		};
		{//0特殊设置控制点
			int i = 0;
			OpenMesh::Vec3d start = square_vertices[i];
			OpenMesh::Vec3d end = square_vertices[(i + 1) % 8];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 2;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) /3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		for (int i = 1; i < 8; i++) {
			OpenMesh::Vec3d start = square_vertices[i];
			OpenMesh::Vec3d end = square_vertices[(i + 1) % 8];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_specialgreen_weight123({0});//计算3次cage到3次控制点的权重
	}

	//第3,5条边特殊的GCfish,chousifh
	if (0)
	{
		int N = 12;
		curvecage2.resize(N);
		std::vector<OpenMesh::Vec3d> polygon_vertices = {
			OpenMesh::Vec3d(0.621947, -3.00053, 0),
			OpenMesh::Vec3d(2.33876, -3.72963, 0),
			OpenMesh::Vec3d(5.77752, -1.50034, 0),
			OpenMesh::Vec3d(2.85876, 2.93215, 0),
			OpenMesh::Vec3d(-5.38815, 3.38406, 0),
			OpenMesh::Vec3d(-5.53413, 2.63302, 0),
			OpenMesh::Vec3d(0.310345, 2.27181, 0),
			OpenMesh::Vec3d(-3.27319, 0.781111, 0),
			OpenMesh::Vec3d(-3.53785, -1.70024, 0),
			OpenMesh::Vec3d(-1.94263, -3.6211, 0),
			OpenMesh::Vec3d(0.503098, -2.97635, 0),
			OpenMesh::Vec3d(1.78622, -2.87407, 0),
		};
		if (1)
		{//3特殊设置控制点
			int i = 3;
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 1.2;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) / 8;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		if(1)
		{//3特殊设置控制点
			int i = 7;
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 1.2;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) / 2;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		{//5特殊设置控制点
			int i =5;
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 1.2;//2
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) / 8;//3

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		for (int i = 0; i < N; i++) {
			if (i==3||i == 7 || i == 5)
				continue;
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_specialgreen_weight123({3,5,7});//计算3次cage到3次控制点的权重
	}

	//正儿八经的GC,zoom fish
	if (0) {
		int N = 4;//zoom_fish 4
		curvecage2.resize(N);

		// 逆时针定义4个顶点
		std::vector<OpenMesh::Vec3d> square_vertices = {
			OpenMesh::Vec3d(15.9365 + 0.01,  1.25963 + 0.01, 0),
			OpenMesh::Vec3d(-16.3751 - 0.01,  1.25963 + 0.01, 0),
			OpenMesh::Vec3d(-16.3751 - 0.01, -1.86695 - 0.01, 0),
			OpenMesh::Vec3d(15.9365 + 0.01, -1.86695 - 0.01, 0),
			
		};
		for (int i = 0; i < N; i++) {
			OpenMesh::Vec3d start = square_vertices[i];
			OpenMesh::Vec3d end = square_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_green_weight123();//计算3次cage到3次控制点的权重
	}
	//第0条边特殊的GCfish,zoomfifh
	if (0)
	{
		int N = 4;//zoom_fish 4
		curvecage2.resize(N);

		// 逆时针定义4个顶点
		std::vector<OpenMesh::Vec3d> polygon_vertices = {
			OpenMesh::Vec3d(15.9365 + 0.01,  1.25963 + 0.01, 0),
			OpenMesh::Vec3d(-16.3751 - 0.01,  1.25963 + 0.01, 0),
			OpenMesh::Vec3d(-16.3751 - 0.01, -1.86695 - 0.01, 0),
			OpenMesh::Vec3d(15.9365 + 0.01, -1.86695 - 0.01, 0),

		};
		if (0)
		{//0特殊设置控制点
			int i = 0;
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) *2/ 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		if (0)
		{//0特殊设置控制点
			int i = 0;
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) * 0.6;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) *0.9;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		if (1)
		{//0特殊设置控制点
			int i = 0;
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) * 0.2;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 0.8;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		if(0)
		{//5特殊设置控制点
			int i = 5;
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 1.2;//2
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) / 8;//3

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		for (int i = 0; i < N; i++) {
			if (i == 0)
				continue;
			OpenMesh::Vec3d start = polygon_vertices[i];
			OpenMesh::Vec3d end = polygon_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_specialgreen_weight123({ 0 });//计算3次cage到3次控制点的权重
	}

	//正儿八经的GC,circle
	if (0) {
		int N = 8;//circle 4
		curvecage2.resize(N);

		// 逆时针定义8个顶点
		std::vector<OpenMesh::Vec3d> square_vertices = {
			OpenMesh::Vec3d(4.96702, 1.64565, 0),
			OpenMesh::Vec3d(2.85338, 4.94866, 0),
			OpenMesh::Vec3d(-0.779829, 5.71698, 0),
OpenMesh::Vec3d(-4.53453, 3.86186, 0),
OpenMesh::Vec3d(-5.43369, -0.358081, 0),
OpenMesh::Vec3d(-3.93307, -4.14067, 0),
OpenMesh::Vec3d(0.472388, -5.10482, 0),
OpenMesh::Vec3d(4.43436, -3.03273, 0),


		};
		for (int i = 0; i < N; i++) {
			OpenMesh::Vec3d start = square_vertices[i];
			OpenMesh::Vec3d end = square_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_green_weight123();//计算3次cage到3次控制点的权重
	}
	//线性cage
	if (0) {
		int N = 8*3;//circle 
		degree = 1;
		curvecage2.resize(N);
		std::vector<OpenMesh::Vec3d>square_vertices = {
			OpenMesh::Vec3d(3.87407, 0.191556, 0),
OpenMesh::Vec3d(4.09041, -1.15695, 0),
OpenMesh::Vec3d(3.90557, -1.76057, 0),
OpenMesh::Vec3d(3.44068, -3.85162, 0),
OpenMesh::Vec3d(5.42266, -1.35024, 0),
OpenMesh::Vec3d(3.33319, 4.06603, 0),
OpenMesh::Vec3d(-0.485222, 4.40195, 0),
OpenMesh::Vec3d(0.747137, 4.73402, 0),
OpenMesh::Vec3d(2.05847, 4.4998, 0),
OpenMesh::Vec3d(3.83488, 3.89295, 0),
OpenMesh::Vec3d(1.39062, 6.04187, 0),
OpenMesh::Vec3d(-3.94472, 3.75008, 0),
OpenMesh::Vec3d(-4.09774, 0.129278, 0),
OpenMesh::Vec3d(-4.37443, 1.31151, 0),
OpenMesh::Vec3d(-4.17104, 2.38932, 0),
OpenMesh::Vec3d(-3.86462, 4.22689, 0),
OpenMesh::Vec3d(-5.74373, 1.56731, 0),
OpenMesh::Vec3d(-3.57387, -3.48392, 0),
OpenMesh::Vec3d(0.194383, -3.81429, 0),
OpenMesh::Vec3d(-0.758251, -4.06885, 0),
OpenMesh::Vec3d(-2.37286, -3.93278, 0),
OpenMesh::Vec3d(-4.21109, -3.22378, 0),
OpenMesh::Vec3d(-1.77591, -5.4711, 0),
OpenMesh::Vec3d(3.56099, -3.33112, 0),
		};
		/*{OpenMesh::Vec3d(3.87407, 0.191556, 0),
			OpenMesh::Vec3d(5.79932, -0.8039, 0),
			OpenMesh::Vec3d(6.87269, -0.0378885, 0),
			OpenMesh::Vec3d(7.54488, 1.72362, 0),
			OpenMesh::Vec3d(5.75912, 0.547401, 0),
			OpenMesh::Vec3d(3.33319, 4.06603, 0),
			OpenMesh::Vec3d(-0.485222, 4.40195, 0),
			OpenMesh::Vec3d(0.42163, 6.0035, 0),
			OpenMesh::Vec3d(0.190249, 7.34517, 0),
			OpenMesh::Vec3d(-1.22957, 7.85349, 0),
			OpenMesh::Vec3d(-0.500449, 6.15799, 0),
			OpenMesh::Vec3d(-3.94472, 3.75008, 0),
			OpenMesh::Vec3d(-4.09774, 0.129278, 0),
			OpenMesh::Vec3d(-5.63014, 1.21918, 0),
			OpenMesh::Vec3d(-7.03331, 0.985885, 0),
			OpenMesh::Vec3d(-7.98259, -1.01753, 0),
			OpenMesh::Vec3d(-5.89146, 0.256206, 0),
			OpenMesh::Vec3d(-3.57387, -3.48392, 0),
			OpenMesh::Vec3d(0.194383, -3.81429, 0),
			OpenMesh::Vec3d(-0.924447, -5.26916, 0),
			OpenMesh::Vec3d(-0.6001, -6.86891, 0),
			OpenMesh::Vec3d(1.23645, -7.72955, 0),
			OpenMesh::Vec3d(0.14458, -5.74809, 0),
			OpenMesh::Vec3d(3.56099, -3.33112, 0),
		}*/
		for (int i = 0; i < N; i++) {
			OpenMesh::Vec3d start = square_vertices[i];
			OpenMesh::Vec3d end = square_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			

			curvecage2[i] = { p1, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_green_weight121();//计算3次cage到3次控制点的权重
	}
	//paristower
	if (0) {
		int N = 3;//circle 4
		curvecage2.resize(N);

		// 逆时针定义3个顶点
		std::vector<OpenMesh::Vec3d> square_vertices = {
			
			
			OpenMesh::Vec3d(2.21312, -3.64502, 0),
			OpenMesh::Vec3d(-0.00294265, 4.30579, 0),
			OpenMesh::Vec3d(-2.13172, -3.62273, 0),
		};
		for (int i = 0; i < N; i++) {
			OpenMesh::Vec3d start = square_vertices[i];
			OpenMesh::Vec3d end = square_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_green_weight123();//计算3次cage到3次控制点的权重
	}
	//paristower
	if (1) {
		int N = 5;//circle 4
		curvecage2.resize(N);

		// 逆时针定义3个顶点
		std::vector<OpenMesh::Vec3d> square_vertices = {
			OpenMesh::Vec3d(-12.729, -6.31988, 0),
				OpenMesh::Vec3d(5.51534, -6.20705, 0),
                OpenMesh::Vec3d(6.0044, 9.48202, 0),
				OpenMesh::Vec3d(-4.21192, 10.1169, 0),
				OpenMesh::Vec3d(-12.3981, 4.54636, 0),
		};
		for (int i = 0; i < N; i++) {
			OpenMesh::Vec3d start = square_vertices[i];
			OpenMesh::Vec3d end = square_vertices[(i + 1) % N];

			OpenMesh::Vec3d p1 = start;
			OpenMesh::Vec3d p4 = end;

			// Intermediate control points (simple linear interpolation for demonstration)
			OpenMesh::Vec3d p2 = p1 + (p4 - p1) / 3;
			OpenMesh::Vec3d p3 = p1 + (p4 - p1) * 2 / 3;

			curvecage2[i] = { p1, p2, p3, p4 };

		}
		//通过curvecage2计算cc_points;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}


		CC_mesh = createMeshFromCurveCage(CC_points);
		curvecage2 = CCpoints_fromCCmesh(CC_mesh, degree);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
		calculate_green_weight123();//计算3次cage到3次控制点的权重
	}
	Set_Texture_coord();//设置纹理
}

void MeshViewerWidget::CalculateWeight_Test(void)
{
	drawmode = CURVECAGE;
	if (true)//if mode
	{
		weights.resize(mesh.n_vertices());
		// 遍历每个顶点，将每个内向量的大小调整为 CC_points.size()
		for (auto& weight : weights) {
			weight.resize(CC_points.size());
		}
		int N = 4;
		curvecage2.resize(N);
		curvecage2[0] = { OpenMesh::Vec3d(1,0,0),OpenMesh::Vec3d(1,1,0), OpenMesh::Vec3d(0,1,0) };
		curvecage2[1] = { OpenMesh::Vec3d(0,1,0),OpenMesh::Vec3d(-1,1,0), OpenMesh::Vec3d(-1,0,0) };
		curvecage2[2] = { OpenMesh::Vec3d(-1,0,0),OpenMesh::Vec3d(-1,-1,0), OpenMesh::Vec3d(0,-1,0) };
		curvecage2[3] = { OpenMesh::Vec3d(0,-1,0),OpenMesh::Vec3d(1,-1,0), OpenMesh::Vec3d(1,0,0) };
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < degree; j++) {
				CC_points.push_back(curvecage2[i][j]);
			}
		}

		CC_mesh = createMeshFromCurveCage(CC_points);
		auto CagevertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
		for (auto vh : CC_mesh.vertices())
		{
			CagevertexState[vh] = NotSelected;
		}
	}
}

void DrawPoints3d(const vector<OpenMesh::Vec3d>& Points, float r, float g, float b, float pointsize) {
	if (!Points.empty())
	{
		glEnable(GL_POINT_SMOOTH);
		glPointSize(pointsize);
		glColor3f(r, g, b);
		
		glBegin(GL_POINTS);
		for (const auto& p : Points)
		{
			glVertex3d(p[0], p[1], p[2]);
		}
		glEnd();
	}
}

Mesh createMeshFromCurveCage(const std::vector<Mesh::Point>& curvecage2) {
	// 创建一个网格对象
	Mesh mesh;

	// 存储顶点句柄的向量
	std::vector<Mesh::VertexHandle> vhandles;

	// 将 curvecage2 中的点添加到网格的顶点
	for (const auto& point : curvecage2) {
		vhandles.push_back(mesh.add_vertex(point));
	}
	
	//// 按照顺序用边连接这些顶点
	//for (size_t i = 0; i < vhandles.size(); ++i) {
	//	Mesh::VertexHandle vh0 = vhandles[i];
	//	Mesh::VertexHandle vh1 = vhandles[i];
	//	Mesh::VertexHandle vh2 = vhandles[(i + 1) % vhandles.size()];
	//	std::vector<Mesh::VertexHandle> face_vhandles = { vh0, vh1, vh2 };
	//	mesh.add_face(face_vhandles);
	//}

	return mesh;
}

std::vector<std::vector<Mesh::Point>> CCpoints_fromCCmesh(const Mesh CC_mesh,int degree) {
	std::vector<std::vector<Mesh::Point>> CCpoints;
	int segs = CC_mesh.n_vertices() / degree;
	CCpoints.resize(segs);
	for (int i = 0; i < segs; i++)
	{
		for (int j = 0; j < degree+1; j++)
		{
			if (i * degree + j == CC_mesh.n_vertices()) {
				auto vh = CC_mesh.vertex_handle(0);
				CCpoints[i].push_back(CC_mesh.point(vh));
			}
			else {
				auto vh = CC_mesh.vertex_handle(i * degree + j);
				CCpoints[i].push_back(CC_mesh.point(vh));
			}
			
		}
	}
	return CCpoints;
}

std::vector<double> mvc(const Mesh::Point p, const std::vector<Mesh::Point> vts) {
	int n_v = vts.size();
	std::vector<double> w(n_v);
	double sum = 0;
	for (int i = 0; i < n_v; i++) {
		Mesh::Point v_left = vts[(i - 1) % (n_v )];
		Mesh::Point v = vts[i];
		Mesh::Point v_right = vts[(i+1)%n_v];
		double alpha_i_1 = acos((v - p).dot(v_left - p) / ((v - p).norm() * (v_left - p).norm()));
		double alpha_i = acos((v - p).dot(v_right - p) / ((v - p).norm() * (v_right - p).norm()));
		double wi = (tan(alpha_i_1 / 2) + tan(alpha_i / 2)) / (v - p).norm();
		w[i] = wi;
		sum += wi;

	}
	for (int i = 0; i < n_v; i++) {
		w[i] /= sum;
	}
	return w;
}

//计算积分，分子t的m次方,分母是t的1次方
double  F1_n(const Mesh::Point eta, const Mesh::Point c0, const Mesh::Point c1, int m)
{

	// 定义lambda函数
	auto Un = [](complex<double> x, int n) -> complex<double> {
		complex<double> sum = 0;
		for (int k = 1; k < n; k++) {
			sum += std::pow(x, k) / static_cast<double>(n-k);
		}
		return sum;
	};

	complex<double> com_c0(c0[0] - eta[0], c0[1] - eta[1]);
	complex<double> com_c1(c1[0], c1[1]);
	
	complex<double> w = - com_c0/com_c1;
	
	double result = imag(std::pow(w, m) * (std::log(1.0 - w) - std::log(-w)) + Un(w, m));
	result /= imag(w);
	result /= (2 * M_PI * norm(com_c1));//norm已经是平方
	return result;
}

//计算积分，分子t的m次方,分母是t的2次方
double F2_n(const Mesh::Point eta, const Mesh::Point c0, const Mesh::Point c1, const Mesh::Point c2, int m)
{

	// 定义lambda函数
	auto Hyper = [](complex<double> x, int n) -> complex<double> {
		complex<double> sum = 0;
		for (int m = 1; m <= n; ++m) {
			sum += std::pow(x, m) / static_cast<double>(m);
		}
		return complex<double>(n + 1) / std::pow(x, n + 1) * (-std::log(1.0 - x) - sum);
	};

	complex<double> com_c0(c0[0]-eta[0], c0[1]-eta[1]);
	complex<double> com_c1(c1[0], c1[1]);
	complex<double> com_c2(c2[0], c2[1]);
	complex<double> w1 = (-com_c1 - sqrt(com_c1 * com_c1 - com_c0 * com_c2 * 4.0)) / (2.0 * com_c2);
	complex<double> w2 = (-com_c1 + sqrt(com_c1 * com_c1 - com_c0 * com_c2 * 4.0)) / (2.0 * com_c2);

	double a = w1.real(); double b = abs(w1.imag());
	double c = w2.real(); double d = abs(w2.imag());
	w1 = complex<double>(a, b);
	w2 = complex<double>(c, d);
	complex<double> w1_c = conj(w1);
	complex<double> w2_c = conj(w2);
	double result = imag(w1 * d * norm(w2) * (w1 * w1 - 2. * c * w1 + c * c + d * d) * Hyper(1. / w1_c, m) +
		                          w2 * b * norm(w1) * (w2 * w2 - 2. * a * w2 + a * a + b * b) * Hyper(1. / w2_c, m));
	result /= ( b * norm(w1) * ((b - d) * (b - d) + (a - c) * (a - c)) * norm(w2) * d * ((b + d) * (b + d) + (a - c) * (a - c)) * (m + 1));
	result /= (2 * M_PI*norm(com_c2));
	return result;
}
//输入顺序为3，2，1，0次项系数
void Cardano(complex<double> a, complex<double> b, complex<double> c, complex<double>d, complex<double>& x1, complex<double>& x2, complex<double>& x3)
{
	complex<double> u = (9.0 * a * b * c - 27.0 * a * a * d - 2.0 * b * b * b) / (54.0 * a * a * a);
	complex<double> v = sqrt(3.0 * (4.0 * a * c * c * c - b * b * c * c - 18.0 * a * b * c * d + 27.0 * a * a * d * d + 4.0 * b * b * b * d)) / (18.0 * a * a);

	complex<double> m;
	if (norm(u + v) >= norm(u - v))
	{
		m = pow(u + v, 1.0 / 3.0);
	}
	else
	{
		m = pow(u - v, 1.0 / 3.0);
	}

	complex<double> n;
	if (norm(m) != 0)
	{
		n = (b * b - 3.0 * a * c) / (9.0 * a * a * m);
	}
	else
	{
		n = 0;
	}

	complex<double> omega1 = complex<double>(-0.5, sqrt(3.0) / 2.0);
	complex<double> omega2 = complex<double>(-0.5, -sqrt(3.0) / 2.0);

	x1 = m + n - b / (3.0 * a);
	x2 = omega1 * m + omega2 * n - b / (3.0 * a);
	x3 = omega2 * m + omega1 * n - b / (3.0 * a);
}

//计算积分，分子t的m次方,分母是t的3次方
double F3_n(const Mesh::Point eta, const Mesh::Point c0, const Mesh::Point c1, const Mesh::Point c2, const Mesh::Point c3, int m)
{

	auto accumulateSum =[](complex<double> w,int m) {
		complex<double> sum = 0.0;
		for (int k = 1; k <= m; ++k) {
			sum += 1.0 / static_cast<double>(k) * pow(w, m - k);
		}
		return sum;
	};


	complex<double> com_c0(c0[0] - eta[0], c0[1] - eta[1]);
	complex<double> com_c1(c1[0], c1[1]);
	complex<double> com_c2(c2[0], c2[1]);
	complex<double> com_c3(c3[0], c3[1]);
	complex<double> w1;
	complex<double> w2;
	complex<double> w3;
	Cardano(com_c3, com_c2, com_c1, com_c0, w1, w2, w3);
	

	double a = w1.real(); double b = abs(w1.imag());
	double c = w2.real(); double d = abs(w2.imag());
	double e = w3.real(); double f = abs(w3.imag());
	w1 = complex<double>(a, b);
	w2 = complex<double>(c, d);
	w3 = complex<double>(e, f);
	complex<double> w1_c = conj(w1);
	complex<double> w2_c = conj(w2);
	complex<double> w3_c = conj(w3);

	complex<double> term1 =- pow(w1_c, m) * log(1.0 - 1.0 / w1_c) - accumulateSum(w1_c, m);
	complex<double> denominator1 = b  * (w1_c - w2_c) * (w1_c - w2) * (w1_c - w3_c) * (w1_c - w3);

	complex<double> term2 = -pow(w2_c, m) * log(1.0 - 1.0 / w2_c) - accumulateSum(w2_c, m);
	complex<double> denominator2 = d * (w2_c - w1_c) * (w2_c - w1) * (w2_c - w3_c) * (w2_c - w3);

	complex<double> term3 = -pow(w3_c, m) * log(1.0 - 1.0 / w3_c) - accumulateSum(w3_c, m);
	complex<double> denominator3 = f * (w3_c - w1_c) * (w3_c - w1) * (w3_c - w2_c) * (w3_c - w2);

	double result = imag(term1 / denominator1) + imag(term2 / denominator2) + imag(term3 / denominator3);
	
	
	result /= (2 * M_PI * norm(com_c3));
	return result;
}


Mesh::Point evaluate_coor(const std::vector<double> weight, const std::vector<Mesh::Point> ctps)
{
	assert(weight.size() == ctps.size());
	Mesh::Point result(0,0,0);
	for (int i = 0; i < weight.size(); i++)
	{
		result += weight[i] * ctps[i];
	}
	//std::cout << result << std::endl;
	return result;
}

bool MeshViewerWidget::NearestVertex(OpenMesh::Vec3d objCor, OpenMesh::VertexHandle& minVh)
{
	double maxAllowedDis = avgEdgeLength * 0.5;
	double minDis = INFINITY;
	OpenMesh::VertexHandle mv;

	for (auto vh : mesh.vertices())
	{
		auto vp = mesh.point(vh);
		double dis = (objCor - vp).norm();
		if (dis < minDis)
		{
			minDis = dis;
			mv = vh;
		}
	}
	//std::cout << minDis << "," << maxAllowedDis << endl;
	if (minDis <= maxAllowedDis)
	{
		minVh = mv;
		return true;
	}
	else
	{
		return false;
	}
}

bool MeshViewerWidget::NearestCageVertex(OpenMesh::Vec3d objCor, OpenMesh::VertexHandle& minVh)
{
	//cal平均边长的一半
	double totalLength = 0.0;
	size_t numEdges = 0;
	// 获取所有顶点句柄
	std::vector<Mesh::VertexHandle> vhandles;
	for (const auto& vh : CC_mesh.vertices()) {
		vhandles.push_back(vh);
	}

	// 按顺序计算相邻顶点之间的距离
	for (size_t i = 0; i < vhandles.size(); ++i) {
		auto vh0 = vhandles[i];
		auto vh1 = vhandles[(i + 1) % vhandles.size()];  // 循环连接最后一个顶点和第一个顶点

		// 获取两个顶点的位置
		Mesh::Point p0 = CC_mesh.point(vh0);
		Mesh::Point p1 = CC_mesh.point(vh1);

		// 计算欧几里得距离
		double length = std::sqrt(
			std::pow(p1[0] - p0[0], 2) +
			std::pow(p1[1] - p0[1], 2) +
			std::pow(p1[2] - p0[2], 2)
		);

		totalLength += length;
		numEdges++;
	}

	double maxAllowedDis = totalLength/ numEdges * 0.5;
	double minDis = INFINITY;
	OpenMesh::VertexHandle mv;

	for (auto vh : CC_mesh.vertices())
	{
		auto vp = CC_mesh.point(vh);
		double dis = (objCor - vp).norm();
		if (dis < minDis)
		{
			minDis = dis;
			mv = vh;
		}
	}
	//std::cout << minDis << "," << maxAllowedDis << endl;
	if (minDis <= maxAllowedDis)
	{
		minVh = mv;
		return true;
	}
	else
	{
		return false;
	}
}
//将1次bezier曲线的控制点转为多项式基的控制点
void MeshViewerWidget::Bezier1Poly1(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage2poly)
{
	curvecage2poly = curvecage2;
	for (int i = 0; i < curvecage2.size(); i++)
	{
		const Mesh::Point& P0 = curvecage2[i][0];
		const Mesh::Point& P1 = curvecage2[i][1];
		

		Mesh::Point a0 = P0;
		Mesh::Point a1 =(P1 - P0);
		

		curvecage2poly[i][0] = a0;
		curvecage2poly[i][1] = a1;
	}
}

//将2次bezier曲线的控制点转为多项式基的控制点
void MeshViewerWidget::Bezier2Poly2(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage2poly)
{
	curvecage2poly = curvecage2;
	for (int i = 0; i < curvecage2.size(); i++)
	{
		const Mesh::Point& P0 = curvecage2[i][0];
		const Mesh::Point& P1 = curvecage2[i][1];
		const Mesh::Point& P2 = curvecage2[i][2];

		Mesh::Point a0 = P0;
		Mesh::Point a1 = 2 * (P1 - P0);
		Mesh::Point a2 = P2 - 2 * P1 + P0;

		curvecage2poly[i][0] = a0;
		curvecage2poly[i][1] = a1;
		curvecage2poly[i][2] = a2;
	}
}

//将3次bezier曲线的控制点转为多项式基的控制点
void MeshViewerWidget::Bezier2Poly3(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage2poly)
{
	curvecage2poly = curvecage2;
	for (int i = 0; i < curvecage2.size(); i++)
	{
		const Mesh::Point& P0 = curvecage2[i][0];
		const Mesh::Point& P1 = curvecage2[i][1];
		const Mesh::Point& P2 = curvecage2[i][2];
		const Mesh::Point& P3 = curvecage2[i][3];

		Mesh::Point a0 = P0;
		Mesh::Point a1 = 3 * (P1 - P0);
		Mesh::Point a2 = 3 * (P2 - 2 * P1 + P0);
		Mesh::Point a3 = P3 - 3 * P2 + 3 * P1 - P0;

		curvecage2poly[i][0] = a0;
		curvecage2poly[i][1] = a1;
		curvecage2poly[i][2] = a2;
		curvecage2poly[i][3] = a3;
	}
}


//将3次bezier曲线的控制点转为7次bezier的控制点,用于可视化
void MeshViewerWidget::Bezier2Bezier7(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage7)
{
	int n = 3, m = 7;
	curvecage7.resize(curvecage2.size());
	auto curvecage = curvecage2;
	auto shengjie = [](std::vector<Mesh::Point> points, std::vector<Mesh::Point> &points1) -> void {
		int n = points.size() - 1;
		points1.resize(points.size() + 1);
		points1[0] = points[0];
		points1[points.size()] = points[points.size() - 1];
		for (int i = 1; i < points.size(); i++)
		{
			points1[i] = (double)i / (n + 1.) * points[i - 1] + (1. - (double)i / (n + 1)) * points[i];
		}
	};
	for (int seg_i = 0; seg_i < curvecage2.size(); seg_i++) {
		for (int i = 0; i < m - n; i++)
		{
			shengjie(curvecage[seg_i], curvecage7[seg_i]);
			curvecage[seg_i] = curvecage7[seg_i];
		}
	}
	
}

//将2次bezier曲线的控制点转为1次bezier的控制点,用于可视化
void MeshViewerWidget::Bezier2Bezier221(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage1)
{
	curvecage1.clear();
	curvecage1.resize(curvecage2.size());
	for (int seg_i = 0; seg_i < curvecage2.size(); seg_i++) {
		curvecage1[seg_i].push_back(curvecage2[seg_i][0]);
		curvecage1[seg_i].push_back(curvecage2[seg_i][2]);
	}

}
//将3次bezier曲线的控制点转为1次bezier的控制点,用于可视化
void MeshViewerWidget::Bezier2Bezier321(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage1)
{
	curvecage1.clear();
	curvecage1.resize(curvecage2.size());
	for (int seg_i = 0; seg_i < curvecage2.size(); seg_i++) {
		curvecage1[seg_i].push_back(curvecage2[seg_i][0]);
		curvecage1[seg_i].push_back(curvecage2[seg_i][3]);
	}

}

//将3次bezier曲线的控制点转为2次bezier的控制点,用于可视化
void MeshViewerWidget::Bezier2Bezier322(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage1)
{
	curvecage1.clear();
	curvecage1.resize(curvecage2.size());
	for (int seg_i = 0; seg_i < curvecage2.size(); seg_i++) {
		curvecage1[seg_i].push_back(curvecage2[seg_i][0]);
		curvecage1[seg_i].push_back((curvecage2[seg_i][1] + curvecage2[seg_i][2]) / 2);
		curvecage1[seg_i].push_back(curvecage2[seg_i][3]);
	}

}

//将2次bezier曲线的控制点转为3次bezier的控制点,用于可视化
void MeshViewerWidget::Bezier2Bezier223(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage3)
{
	int n = 2, m = 3;
	curvecage3.resize(curvecage2.size());
	auto curvecage = curvecage2;
	auto shengjie = [](std::vector<Mesh::Point> points, std::vector<Mesh::Point>& points1) -> void {
		int n = points.size() - 1;
		points1.resize(points.size() + 1);
		points1[0] = points[0];
		points1[points.size()] = points[points.size() - 1];
		for (int i = 1; i < points.size(); i++)
		{
			points1[i] = (double)i / (n + 1.) * points[i - 1] + (1. - (double)i / (n + 1)) * points[i];
		}
	};
	for (int seg_i = 0; seg_i < curvecage2.size(); seg_i++) {
		for (int i = 0; i < m - n; i++)
		{
			shengjie(curvecage[seg_i], curvecage3[seg_i]);
			curvecage[seg_i] = curvecage3[seg_i];
		}
	}

}

//将2次bezier曲线的控制点转为7次bezier的控制点,用于可视化
void MeshViewerWidget::Bezier2Bezier227(std::vector<std::vector<Mesh::Point>> curvecage2, std::vector<std::vector<Mesh::Point>>& curvecage3)
{
	int n = 2, m = 7;
	curvecage3.resize(curvecage2.size());
	auto curvecage = curvecage2;
	auto shengjie = [](std::vector<Mesh::Point> points, std::vector<Mesh::Point>& points1) -> void {
		int n = points.size() - 1;
		points1.resize(points.size() + 1);
		points1[0] = points[0];
		points1[points.size()] = points[points.size() - 1];
		for (int i = 1; i < points.size(); i++)
		{
			points1[i] = (double)i / (n + 1.) * points[i - 1] + (1. - (double)i / (n + 1)) * points[i];
		}
	};
	for (int seg_i = 0; seg_i < curvecage2.size(); seg_i++) {
		for (int i = 0; i < m - n; i++)
		{
			shengjie(curvecage[seg_i], curvecage3[seg_i]);
			curvecage[seg_i] = curvecage3[seg_i];
		}
	}

}


//将7次bezier曲线的控制点转为多项式基的控制点
void MeshViewerWidget::Bezier2Poly7(std::vector<std::vector<Mesh::Point>> curvecage7, std::vector<std::vector<Mesh::Point>>& curvecage7poly)
{
	int degree = 7;
	curvecage7poly.resize(curvecage7.size());
	// Precompute binomial coefficients
	std::vector<std::vector<double>> binomials(degree + 1, std::vector<double>(degree + 1, 0));
	for (int n = 0; n <= degree; ++n) {
		binomials[n][0] = 1;
		for (int k = 1; k <= n; ++k) {
			binomials[n][k] = binomials[n - 1][k - 1] + binomials[n - 1][k];
		}
	}

	for (int c = 0; c < curvecage7.size(); ++c) {
		std::vector<Mesh::Point> poly(degree + 1, { 0, 0, 0 });
		// Calculate polynomial coefficients
		for (int k = 0; k <= degree; ++k) {
			Mesh::Point sum = { 0, 0, 0 };
			for (int i = 0; i <= degree; ++i) {
				for (int j = 0; j <= degree - i; ++j) {
					if (i + j == k) {
						double coeff = binomials[degree][i] * binomials[degree - i][j] * (j % 2 == 0 ? 1 : -1);
						sum = sum + curvecage7[c][i] * coeff;
					}
				}
			}
			poly[k] = sum;
		}
		curvecage7poly[c] = poly;
	}
}

//2到2的权重计算,次数m=2
void MeshViewerWidget::calculate_green_weight222(void)
{
	assert(degree == 2);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};
	
	weights.resize(mesh.n_vertices());
	for (auto& weight : weights) {
		weight.resize(curvecage2.size()*(2*degree+1));//使用论文中的每段曲线使用幂基来计算
	}
	auto poly_Ctps = curvecage2;
	Bezier2Poly2(curvecage2, poly_Ctps);
	double max_err=0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1]; auto c2 = poly_Ctp[2];
			double F0,F1,F2, F3, F4, F5;
			F0 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 0);
			F1 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 1);
			F2 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 2);
			F3 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 3);
			F4 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 4);
			F5 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 5);
			weights[v_id][5 * i + 0] = (c0 - eta).dot(arthono(c1)) * F0 + 2 * (c0 - eta).dot(arthono(c2)) * F1 + (c1).dot(arthono(c2)) * F2;
			weights[v_id][5 * i + 1] = (c0 - eta).dot(arthono(c1)) * F1 + 2 * (c0 - eta).dot(arthono(c2)) * F2 + (c1).dot(arthono(c2)) * F3;
			weights[v_id][5 * i + 2] = (c0 - eta).dot(arthono(c1)) * F2 + 2 * (c0 - eta).dot(arthono(c2)) * F3 + (c1).dot(arthono(c2)) * F4;
			weights[v_id][5 * i + 3] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F2 + 3 * (c1).dot(c2) * F3 + 2 * c2.dot(c2) * F4 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			weights[v_id][5 * i + 4] = (c0 - eta).dot(c1) * F2 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F3 + 3 * (c1).dot(c2) * F4 + 2 * c2.dot(c2) * F5 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			test_eta += weights[v_id][5 * i + 0] * poly_Ctp[0];
			test_eta += weights[v_id][5 * i + 1] * poly_Ctp[1];
			test_eta += weights[v_id][5 * i + 2] * poly_Ctp[2];
			test_eta += weights[v_id][5 * i + 3] * arthono(poly_Ctp[1]);
			test_eta += weights[v_id][5 * i + 4] * arthono(poly_Ctp[2]);
		}
		//check
		/*std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 0) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 1) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 2) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 5) << std::endl;*/
		//std::cout << eta << endl;
		//std::cout << "check:" << test_eta-eta << std::endl;
		if ((test_eta - eta).norm() > max_err)
			max_err = (test_eta - eta).norm();
	}
	std::cout << "max norm err:"<<max_err << std::endl;
	

}

//2到1的权重计算,次数m=2
void MeshViewerWidget::calculate_green_weight221(void)
{
	assert(degree == 2);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};

	weights.resize(mesh.n_vertices());
	int N = 2 * todegree + 1;
	for (auto& weight : weights) {
		weight.resize(curvecage2.size() * N);//使用论文中的每段曲线使用幂基来计算
	}
	auto poly_Ctps = curvecage2;
	Bezier2Poly2(curvecage2, poly_Ctps);
	double max_err = 0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1]; auto c2 = poly_Ctp[2];
			double F0, F1, F2, F3, F4, F5;
			F0 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 0);
			F1 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 1);
			F2 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 2);
			F3 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 3);
			F4 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 4);
			F5 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 5);
			weights[v_id][N * i + 0] = (c0 - eta).dot(arthono(c1)) * F0 + 2 * (c0 - eta).dot(arthono(c2)) * F1 + (c1).dot(arthono(c2)) * F2;
			weights[v_id][N * i + 1] = (c0 - eta).dot(arthono(c1)) * F1 + 2 * (c0 - eta).dot(arthono(c2)) * F2 + (c1).dot(arthono(c2)) * F3;
			
			weights[v_id][N * i + 2] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F2 + 3 * (c1).dot(c2) * F3 + 2 * c2.dot(c2) * F4 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			
		}
		//check
	
		//std::cout << eta << endl;
		//std::cout << "check:" << test_eta-eta << std::endl;
		/*if ((test_eta - eta).norm() > max_err)
			max_err = (test_eta - eta).norm();*/
	}
	std::cout << "max norm err:" << max_err << std::endl;


}
//2到3的权重计算,次数m=2
void MeshViewerWidget::calculate_green_weight223(void)
{
	assert(degree == 2);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};

	weights.resize(mesh.n_vertices());
	int N = 2 * todegree + 1;
	for (auto& weight : weights) {
		weight.resize(curvecage2.size() * N);//使用论文中的每段曲线使用幂基来计算
	}
	auto poly_Ctps = curvecage2;
	Bezier2Poly2(curvecage2, poly_Ctps);
	double max_err = 0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1]; auto c2 = poly_Ctp[2];
			double F0, F1, F2, F3, F4, F5, F6;
			F0 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 0);
			F1 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 1);
			F2 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 2);
			F3 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 3);
			F4 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 4);
			F5 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 5);
			F6 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 6);
			weights[v_id][N * i + 0] = (c0 - eta).dot(arthono(c1)) * F0 + 2 * (c0 - eta).dot(arthono(c2)) * F1 + (c1).dot(arthono(c2)) * F2;
			weights[v_id][N * i + 1] = (c0 - eta).dot(arthono(c1)) * F1 + 2 * (c0 - eta).dot(arthono(c2)) * F2 + (c1).dot(arthono(c2)) * F3;
			weights[v_id][N * i + 2] = (c0 - eta).dot(arthono(c1)) * F2 + 2 * (c0 - eta).dot(arthono(c2)) * F3 + (c1).dot(arthono(c2)) * F4;
			weights[v_id][N * i + 3] = (c0 - eta).dot(arthono(c1)) * F3 + 2 * (c0 - eta).dot(arthono(c2)) * F4 + (c1).dot(arthono(c2)) * F5;
			weights[v_id][N * i + 4] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F2 + 3 * (c1).dot(c2) * F3 + 2 * c2.dot(c2) * F4 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			weights[v_id][N * i + 5] = (c0 - eta).dot(c1) * F2 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F3 + 3 * (c1).dot(c2) * F4 + 2 * c2.dot(c2) * F5 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			weights[v_id][N * i + 6] = (c0 - eta).dot(c1) * F3 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F4 + 3 * (c1).dot(c2) * F5 + 2 * c2.dot(c2) * F6 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			
		}
		
	}
	std::cout << "max norm err:" << max_err << std::endl;


}

//2到7的权重计算,次数m=2
void MeshViewerWidget::calculate_green_weight227(void)
{
	assert(degree == 2);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};

	weights.resize(mesh.n_vertices());
	int N = 2 * todegree + 1;
	for (auto& weight : weights) {
		weight.resize(curvecage2.size() * N);//使用论文中的每段曲线使用幂基来计算
	}
	auto poly_Ctps = curvecage2;
	Bezier2Poly2(curvecage2, poly_Ctps);
	double max_err = 0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1]; auto c2 = poly_Ctp[2];
			double F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, F10;
			F0 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 0);
			F1 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 1);
			F2 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 2);
			F3 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 3);
			F4 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 4);
			F5 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 5);
			F6 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 6);
			F7 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 7);
			F8 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 8);
			F9 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 9);
			F10 = F2_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], 10);
			weights[v_id][N * i + 0] = (c0 - eta).dot(arthono(c1)) * F0 + 2 * (c0 - eta).dot(arthono(c2)) * F1 + (c1).dot(arthono(c2)) * F2;
			weights[v_id][N * i + 1] = (c0 - eta).dot(arthono(c1)) * F1 + 2 * (c0 - eta).dot(arthono(c2)) * F2 + (c1).dot(arthono(c2)) * F3;
			weights[v_id][N * i + 2] = (c0 - eta).dot(arthono(c1)) * F2 + 2 * (c0 - eta).dot(arthono(c2)) * F3 + (c1).dot(arthono(c2)) * F4;
			weights[v_id][N * i + 3] = (c0 - eta).dot(arthono(c1)) * F3 + 2 * (c0 - eta).dot(arthono(c2)) * F4 + (c1).dot(arthono(c2)) * F5;
			weights[v_id][N * i + 4] = (c0 - eta).dot(arthono(c1)) * F4 + 2 * (c0 - eta).dot(arthono(c2)) * F5 + (c1).dot(arthono(c2)) * F6;
			weights[v_id][N * i + 5] = (c0 - eta).dot(arthono(c1)) * F5 + 2 * (c0 - eta).dot(arthono(c2)) * F6 + (c1).dot(arthono(c2)) * F7;
			weights[v_id][N * i + 6] = (c0 - eta).dot(arthono(c1)) * F6 + 2 * (c0 - eta).dot(arthono(c2)) * F7 + (c1).dot(arthono(c2)) * F8;
			weights[v_id][N * i + 7] = (c0 - eta).dot(arthono(c1)) * F7 + 2 * (c0 - eta).dot(arthono(c2)) * F8 + (c1).dot(arthono(c2)) * F9;
			weights[v_id][N * i + 8] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F2 + 3 * (c1).dot(c2) * F3 + 2 * c2.dot(c2) * F4 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			weights[v_id][N * i + 9] = (c0 - eta).dot(c1) * F2 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F3 + 3 * (c1).dot(c2) * F4 + 2 * c2.dot(c2) * F5 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			weights[v_id][N * i + 10] = (c0 - eta).dot(c1) * F3 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F4 + 3 * (c1).dot(c2) * F5 + 2 * c2.dot(c2) * F6 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			weights[v_id][N * i + 11] = (c0 - eta).dot(c1) * F4 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F5 + 3 * (c1).dot(c2) * F6 + 2 * c2.dot(c2) * F7 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			weights[v_id][N * i + 12] = (c0 - eta).dot(c1) * F5 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F6 + 3 * (c1).dot(c2) * F7 + 2 * c2.dot(c2) * F8 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			weights[v_id][N * i + 13] = (c0 - eta).dot(c1) * F6 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F7 + 3 * (c1).dot(c2) * F8 + 2 * c2.dot(c2) * F9 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);
			weights[v_id][N * i + 14] = (c0 - eta).dot(c1) * F7 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F8 + 3 * (c1).dot(c2) * F9 + 2 * c2.dot(c2) * F10 - log(norm(c0 + c1 + c2 - eta)) / (2 * M_PI);


		}

	}
	std::cout << "max norm err:" << max_err << std::endl;


}


//特殊的1到3的权重计算，第一条边是三次的
void MeshViewerWidget::calculate_specialgreen_weight123(std::vector<int>sp_id)
{
	assert(degree == 3);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};

	weights.resize(mesh.n_vertices());
	for (auto& weight : weights) {
		weight.resize(curvecage2.size() * (2 * degree + 1));//使用论文中的每段曲线使用幂基来计算
	}
	std::vector<std::vector<OpenMesh::Vec3d>> poly_Ctps;
	poly_Ctps.resize(curvecage2.size());
	std::vector<std::vector<OpenMesh::Vec3d>> poly_Ctps_pro = curvecage2;
	Bezier2Poly3(curvecage2, poly_Ctps_pro);
	for(auto spi: sp_id)
		poly_Ctps[spi] = poly_Ctps_pro[spi];
	//poly_Ctps[0] = poly_Ctps_pro[0];
	for (int i = 0; i < poly_Ctps.size(); i++)//从分为两拨
	{
		bool isspec = false;
		for (int spn = 0; spn < sp_id.size(); spn++) {
			if (i == sp_id[spn]) {
				isspec = true;
				poly_Ctps[i] = poly_Ctps_pro[i];
			}
		}
		if (isspec == false) {
			poly_Ctps[i].push_back(curvecage2[i][0]);
			poly_Ctps[i].push_back(curvecage2[i][3] - curvecage2[i][0]);
		}
	}

	double max_err = 0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for(auto i:sp_id)
		{//i表示边界曲线的段数,现在对于第一条边求三次的权重
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1]; auto c2 = poly_Ctp[2]; auto c3 = poly_Ctp[3];
			double F0, F1, F2, F3, F4, F5, F6, F7, F8;
			F0 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 0);
			F1 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 1);
			F2 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 2);
			F3 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 3);
			F4 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 4);
			F5 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 5);
			F6 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 6);
			F7 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 7);
			F8 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 8);
			weights[v_id][(2 * degree + 1) * i + 0] = (c0 - eta).dot(arthono(c1)) * F0 + 2 * (c0 - eta).dot(arthono(c2)) * F1
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F2 + 2 * c1.dot(arthono(c3)) * F3 + c2.dot(arthono(c3)) * F4;

			weights[v_id][(2 * degree + 1) * i + 1] = (c0 - eta).dot(arthono(c1)) * F1 + 2 * (c0 - eta).dot(arthono(c2)) * F2
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F3 + 2 * c1.dot(arthono(c3)) * F4 + c2.dot(arthono(c3)) * F5;

			weights[v_id][(2 * degree + 1) * i + 2] = (c0 - eta).dot(arthono(c1)) * F2 + 2 * (c0 - eta).dot(arthono(c2)) * F3
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F4 + 2 * c1.dot(arthono(c3)) * F5 + c2.dot(arthono(c3)) * F6;

			weights[v_id][(2 * degree + 1) * i + 3] = (c0 - eta).dot(arthono(c1)) * F3 + 2 * (c0 - eta).dot(arthono(c2)) * F4
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F5 + 2 * c1.dot(arthono(c3)) * F6 + c2.dot(arthono(c3)) * F7;

			weights[v_id][(2 * degree + 1) * i + 4] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F2 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F3
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F4 + 5 * c2.dot(c3) * F5 + 3 * c3.dot(c3) * F6 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][(2 * degree + 1) * i + 5] = (c0 - eta).dot(c1) * F2 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F3 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F4
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F5 + 5 * c2.dot(c3) * F6 + 3 * c3.dot(c3) * F7 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][(2 * degree + 1) * i + 6] = (c0 - eta).dot(c1) * F3 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F4 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F5
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F6 + 5 * c2.dot(c3) * F7 + 3 * c3.dot(c3) * F8 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			test_eta += weights[v_id][(2 * degree + 1) * i + 0] * poly_Ctp[0];
			test_eta += weights[v_id][(2 * degree + 1) * i + 1] * poly_Ctp[1];
			test_eta += weights[v_id][(2 * degree + 1) * i + 2] * poly_Ctp[2];
			test_eta += weights[v_id][(2 * degree + 1) * i + 3] * poly_Ctp[3];
			test_eta += weights[v_id][(2 * degree + 1) * i + 4] * arthono(poly_Ctp[1]);
			test_eta += weights[v_id][(2 * degree + 1) * i + 5] * arthono(poly_Ctp[2]);
			test_eta += weights[v_id][(2 * degree + 1) * i + 6] * arthono(poly_Ctp[3]);
		}
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			if (poly_Ctps[i].size() == 2) {
				auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
				auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1];
				double F0, F1, F2, F3, F4, F5, F6, F7, F8;
				F0 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 0);
				F1 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 1);
				F2 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 2);
				F3 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 3);
				F4 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 4);
				F5 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 5);
				F6 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 6);
				F7 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 7);
				F8 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 8);
				weights[v_id][(2 * degree + 1) * i + 0] = (c0 - eta).dot(arthono(c1)) * F0;

				weights[v_id][(2 * degree + 1) * i + 1] = (c0 - eta).dot(arthono(c1)) * F1;

				weights[v_id][(2 * degree + 1) * i + 2] = (c0 - eta).dot(arthono(c1)) * F2;

				weights[v_id][(2 * degree + 1) * i + 3] = (c0 - eta).dot(arthono(c1)) * F3;

				weights[v_id][(2 * degree + 1) * i + 4] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1)) * F2 - log(norm(c0 + c1 - eta)) / (2 * M_PI);

				weights[v_id][(2 * degree + 1) * i + 5] = (c0 - eta).dot(c1) * F2 + (c1.dot(c1)) * F3 - log(norm(c0 + c1 - eta)) / (2 * M_PI);

				weights[v_id][(2 * degree + 1) * i + 6] = (c0 - eta).dot(c1) * F3 + (c1.dot(c1)) * F4 - log(norm(c0 + c1 - eta)) / (2 * M_PI);

				test_eta += weights[v_id][(2 * degree + 1) * i + 0] * poly_Ctp[0];
				test_eta += weights[v_id][(2 * degree + 1) * i + 1] * poly_Ctp[1];
				test_eta += weights[v_id][(2 * degree + 1) * i + 2] * (OpenMesh::Vec3d(0, 0, 0));
				test_eta += weights[v_id][(2 * degree + 1) * i + 3] * (OpenMesh::Vec3d(0, 0, 0));
				test_eta += weights[v_id][(2 * degree + 1) * i + 4] * arthono(poly_Ctp[1]);
				test_eta += weights[v_id][(2 * degree + 1) * i + 5] * (OpenMesh::Vec3d(0, 0, 0));
				test_eta += weights[v_id][(2 * degree + 1) * i + 6] * (OpenMesh::Vec3d(0, 0, 0));
			}
		}
		//check
		/*std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 0) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 1) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 2) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 5) << std::endl;*/
		//std::cout << eta << endl;
		//std::cout << "check:" << test_eta - eta << std::endl;
		if ((test_eta - eta).norm() > max_err)
			max_err = (test_eta - eta).norm();
	}
	std::cout << "max norm err:" << max_err << std::endl;

}


//1到3的权重计算
void MeshViewerWidget::calculate_green_weight123(void)
{
	assert(degree == 3);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};

	weights.resize(mesh.n_vertices());
	for (auto& weight : weights) {
		weight.resize(curvecage2.size() * (2 * degree + 1));//使用论文中的每段曲线使用幂基来计算
	}
	std::vector<std::vector<OpenMesh::Vec3d>> poly_Ctps;
	poly_Ctps.resize(curvecage2.size());
	for (int i = 0; i < poly_Ctps.size(); i++)
	{
		poly_Ctps[i].push_back(curvecage2[i][0]);
		poly_Ctps[i].push_back(curvecage2[i][3]-curvecage2[i][0]);
	}
	
	double max_err = 0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1]; 
			double F0, F1, F2, F3, F4, F5, F6, F7, F8;
			F0 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 0);
			F1 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 1);
			F2 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 2);
			F3 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 3);
			F4 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 4);
			F5 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 5);
			F6 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 6);
			F7 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 7);
			F8 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 8);
			weights[v_id][(2 * degree + 1) * i + 0] = (c0 - eta).dot(arthono(c1)) * F0 ;

			weights[v_id][(2 * degree + 1) * i + 1] = (c0 - eta).dot(arthono(c1)) * F1 ;

			weights[v_id][(2 * degree + 1) * i + 2] = (c0 - eta).dot(arthono(c1)) * F2 ;

			weights[v_id][(2 * degree + 1) * i + 3] = (c0 - eta).dot(arthono(c1)) * F3 ;

			weights[v_id][(2 * degree + 1) * i + 4] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1)) * F2  - log(norm(c0 + c1 - eta)) / (2 * M_PI);

			weights[v_id][(2 * degree + 1) * i + 5] = (c0 - eta).dot(c1) * F2 + (c1.dot(c1)) * F3  - log(norm(c0 + c1 - eta)) / (2 * M_PI);

			weights[v_id][(2 * degree + 1) * i + 6] = (c0 - eta).dot(c1) * F3 + (c1.dot(c1)) * F4  - log(norm(c0 + c1 - eta)) / (2 * M_PI);

			test_eta += weights[v_id][(2 * degree + 1) * i + 0] * poly_Ctp[0];
			test_eta += weights[v_id][(2 * degree + 1) * i + 1] * poly_Ctp[1];
			test_eta += weights[v_id][(2 * degree + 1) * i + 2] * (OpenMesh::Vec3d(0,0,0));
			test_eta += weights[v_id][(2 * degree + 1) * i + 3] * (OpenMesh::Vec3d(0, 0, 0));
			test_eta += weights[v_id][(2 * degree + 1) * i + 4] * arthono(poly_Ctp[1]);
			test_eta += weights[v_id][(2 * degree + 1) * i + 5] * (OpenMesh::Vec3d(0, 0, 0));
			test_eta += weights[v_id][(2 * degree + 1) * i + 6] * (OpenMesh::Vec3d(0, 0, 0));
		}
		//check
		/*std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 0) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 1) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 2) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 5) << std::endl;*/
		//std::cout << eta << endl;
		//std::cout << "check:" << test_eta - eta << std::endl;
		if ((test_eta - eta).norm() > max_err)
			max_err = (test_eta - eta).norm();
	}
	std::cout << "max norm err:" << max_err << std::endl;

}

//1到1的权重计算
void MeshViewerWidget::calculate_green_weight121(void)
{
	assert(degree == 1);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};

	weights.resize(mesh.n_vertices());
	for (auto& weight : weights) {
		weight.resize(curvecage2.size() * (2 * degree + 1));//使用论文中的每段曲线使用幂基来计算
	}
	std::vector<std::vector<OpenMesh::Vec3d>> poly_Ctps;
	poly_Ctps.resize(curvecage2.size());
	for (int i = 0; i < poly_Ctps.size(); i++)
	{
		poly_Ctps[i].push_back(curvecage2[i][0]);
		poly_Ctps[i].push_back(curvecage2[i][degree] - curvecage2[i][0]);
	}

	double max_err = 0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1];
			double F0, F1, F2, F3, F4, F5, F6, F7, F8;
			F0 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 0);
			F1 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 1);
			F2 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 2);
			F3 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 3);
			F4 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 4);
			F5 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 5);
			F6 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 6);
			F7 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 7);
			F8 = F1_n(eta, poly_Ctp[0], poly_Ctp[1], 8);
			weights[v_id][(2 * degree + 1) * i + 0] = (c0 - eta).dot(arthono(c1)) * F0;

			weights[v_id][(2 * degree + 1) * i + 1] = (c0 - eta).dot(arthono(c1)) * F1;

			weights[v_id][(2 * degree + 1) * i + 2] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1)) * F2 - log(norm(c0 + c1 - eta)) / (2 * M_PI);


			test_eta += weights[v_id][(2 * degree + 1) * i + 0] * poly_Ctp[0];
			test_eta += weights[v_id][(2 * degree + 1) * i + 1] * poly_Ctp[1];
			test_eta += weights[v_id][(2 * degree + 1) * i + 2] * arthono(poly_Ctp[1]);
			
		}
		//check
		/*std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 0) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 1) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 2) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 5) << std::endl;*/
		//std::cout << eta << endl;
		//std::cout << "check:" << test_eta - eta << std::endl;
		if ((test_eta - eta).norm() > max_err)
			max_err = (test_eta - eta).norm();
	}
	std::cout << "max norm err:" << max_err << std::endl;

}
//3到1的权重计算
void MeshViewerWidget::calculate_green_weight321(void)
{
	assert(degree == 3);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};

	weights.resize(mesh.n_vertices());
	int N = 2 * todegree + 1;
	for (auto& weight : weights) {
		weight.resize(curvecage2.size() * N);//使用论文中的每段曲线使用幂基来计算
	}
	auto poly_Ctps = curvecage2;
	Bezier2Poly3(curvecage2, poly_Ctps);
	double max_err = 0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1]; auto c2 = poly_Ctp[2]; auto c3 = poly_Ctp[3];
			double F0, F1, F2, F3, F4, F5, F6, F7, F8;
			F0 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 0);
			F1 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 1);
			F2 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 2);
			F3 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 3);
			F4 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 4);
			F5 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 5);
			F6 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 6);
			F7 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 7);
			F8 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 8);
			weights[v_id][N * i + 0] = (c0 - eta).dot(arthono(c1)) * F0 + 2 * (c0 - eta).dot(arthono(c2)) * F1
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F2 + 2 * c1.dot(arthono(c3)) * F3 + c2.dot(arthono(c3)) * F4;

			weights[v_id][N * i + 1] = (c0 - eta).dot(arthono(c1)) * F1 + 2 * (c0 - eta).dot(arthono(c2)) * F2
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F3 + 2 * c1.dot(arthono(c3)) * F4 + c2.dot(arthono(c3)) * F5;

			weights[v_id][N * i + 2] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F2 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F3
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F4 + 5 * c2.dot(c3) * F5 + 3 * c3.dot(c3) * F6 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			test_eta += weights[v_id][N * i + 0] * poly_Ctp[0];
			test_eta += weights[v_id][N * i + 1] * poly_Ctp[1];
			test_eta += weights[v_id][N * i + 2] * arthono(poly_Ctp[1]);
			
		}
		
	}
	std::cout << "max norm err:" << max_err << std::endl;

}

//3到2的权重计算
void MeshViewerWidget::calculate_green_weight322(void)
{
	assert(degree == 3);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};

	weights.resize(mesh.n_vertices());
	int N = 2 * todegree + 1;
	for (auto& weight : weights) {
		weight.resize(curvecage2.size() * N);//使用论文中的每段曲线使用幂基来计算
	}
	auto poly_Ctps = curvecage2;
	Bezier2Poly3(curvecage2, poly_Ctps);
	double max_err = 0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1]; auto c2 = poly_Ctp[2]; auto c3 = poly_Ctp[3];
			double F0, F1, F2, F3, F4, F5, F6, F7, F8;
			F0 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 0);
			F1 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 1);
			F2 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 2);
			F3 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 3);
			F4 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 4);
			F5 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 5);
			F6 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 6);
			F7 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 7);
			F8 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 8);
			weights[v_id][N * i + 0] = (c0 - eta).dot(arthono(c1)) * F0 + 2 * (c0 - eta).dot(arthono(c2)) * F1
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F2 + 2 * c1.dot(arthono(c3)) * F3 + c2.dot(arthono(c3)) * F4;

			weights[v_id][N * i + 1] = (c0 - eta).dot(arthono(c1)) * F1 + 2 * (c0 - eta).dot(arthono(c2)) * F2
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F3 + 2 * c1.dot(arthono(c3)) * F4 + c2.dot(arthono(c3)) * F5;

			weights[v_id][N * i + 2] = (c0 - eta).dot(arthono(c1)) * F2 + 2 * (c0 - eta).dot(arthono(c2)) * F3
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F4 + 2 * c1.dot(arthono(c3)) * F5 + c2.dot(arthono(c3)) * F6;

			weights[v_id][N * i + 3] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F2 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F3
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F4 + 5 * c2.dot(c3) * F5 + 3 * c3.dot(c3) * F6 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][N * i + 4] = (c0 - eta).dot(c1) * F2 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F3 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F4
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F5 + 5 * c2.dot(c3) * F6 + 3 * c3.dot(c3) * F7 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			test_eta += weights[v_id][N * i + 0] * poly_Ctp[0];
			test_eta += weights[v_id][N * i + 1] * poly_Ctp[1];
			test_eta += weights[v_id][N * i + 2] * poly_Ctp[2];
			test_eta += weights[v_id][N * i + 3] * arthono(poly_Ctp[1]);
			test_eta += weights[v_id][N * i + 4] * arthono(poly_Ctp[2]);

		}

	}
	std::cout << "max norm err:" << max_err << std::endl;

}
//3到3的权重计算
void MeshViewerWidget::calculate_green_weight323(void)
{
	assert(degree == 3);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};

	weights.resize(mesh.n_vertices());
	for (auto& weight : weights) {
		weight.resize(curvecage2.size() * (2 * degree + 1));//使用论文中的每段曲线使用幂基来计算
	}
	auto poly_Ctps = curvecage2;
	Bezier2Poly3(curvecage2, poly_Ctps);
	double max_err = 0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1]; auto c2 = poly_Ctp[2]; auto c3 = poly_Ctp[3];
			double F0, F1, F2, F3, F4, F5, F6, F7, F8;
			F0 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 0);
			F1 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 1);
			F2 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 2);
			F3 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 3);
			F4 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 4);
			F5 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 5);
			F6 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 6);
			F7 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 7);
			F8 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 8);
			weights[v_id][(2 * degree + 1) * i + 0] = (c0 - eta).dot(arthono(c1)) * F0 + 2 * (c0 - eta).dot(arthono(c2)) * F1
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F2 + 2 * c1.dot(arthono(c3)) * F3 + c2.dot(arthono(c3)) * F4;

			weights[v_id][(2 * degree + 1) * i + 1] = (c0 - eta).dot(arthono(c1)) * F1 + 2 * (c0 - eta).dot(arthono(c2)) * F2
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F3 + 2 * c1.dot(arthono(c3)) * F4 + c2.dot(arthono(c3)) * F5;

			weights[v_id][(2 * degree + 1) * i + 2] = (c0 - eta).dot(arthono(c1)) * F2 + 2 * (c0 - eta).dot(arthono(c2)) * F3
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F4 + 2 * c1.dot(arthono(c3)) * F5 + c2.dot(arthono(c3)) * F6;

			weights[v_id][(2 * degree + 1) * i + 3] = (c0 - eta).dot(arthono(c1)) * F3 + 2 * (c0 - eta).dot(arthono(c2)) * F4
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F5 + 2 * c1.dot(arthono(c3)) * F6 + c2.dot(arthono(c3)) * F7;

			weights[v_id][(2 * degree + 1) * i + 4] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F2 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F3
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F4 + 5 * c2.dot(c3) * F5 + 3 * c3.dot(c3) * F6 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][(2 * degree + 1) * i + 5] = (c0 - eta).dot(c1) * F2 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F3 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F4
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F5 + 5 * c2.dot(c3) * F6 + 3 * c3.dot(c3) * F7 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][(2 * degree + 1) * i + 6] = (c0 - eta).dot(c1) * F3 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F4 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F5
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F6 + 5 * c2.dot(c3) * F7 + 3 * c3.dot(c3) * F8 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);
			
			test_eta += weights[v_id][(2 * degree + 1) * i + 0] * poly_Ctp[0];
			test_eta += weights[v_id][(2 * degree + 1) * i + 1] * poly_Ctp[1];
			test_eta += weights[v_id][(2 * degree + 1) * i + 2] * poly_Ctp[2];
			test_eta += weights[v_id][(2 * degree + 1) * i + 3] * poly_Ctp[3];
			test_eta += weights[v_id][(2 * degree + 1) * i + 4] * arthono(poly_Ctp[1]);
			test_eta += weights[v_id][(2 * degree + 1) * i + 5] * arthono(poly_Ctp[2]);
			test_eta += weights[v_id][(2 * degree + 1) * i + 6] * arthono(poly_Ctp[3]);
		}
		//check
		/*std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 0) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 1) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 2) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 5) << std::endl;*/
		//std::cout << eta << endl;
		//std::cout << "check:" << test_eta - eta << std::endl;
		if ((test_eta - eta).norm() > max_err)
			max_err = (test_eta - eta).norm();
	}
	std::cout << "max norm err:" << max_err << std::endl;

}
//3到7的权重计算
void MeshViewerWidget::calculate_green_weight327(void)
{
	assert(degree == 3);
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};

	weights.resize(mesh.n_vertices());
	for (auto& weight : weights) {
		weight.resize(curvecage2.size() * (2 * todegree + 1));//使用论文中的每段曲线使用幂基来计算
	}
	auto poly_Ctps = curvecage2;
	Bezier2Poly3(curvecage2, poly_Ctps);
	double max_err = 0;
	for (int v_id = 0; v_id < mesh.n_vertices(); v_id++) {//v_id表示mesh中的点数
		auto vh = mesh.vertex_handle(v_id);
		Mesh::Point eta = mesh.point(vh);
		//Mesh::Point eta(0, 0, 0);
		Mesh::Point test_eta(0, 0, 0);//check equality
		for (int i = 0; i < curvecage2.size(); i++)
		{//i表示边界曲线的段数
			auto poly_Ctp = poly_Ctps[i];//每一段曲线，原来的poly_ctp
			auto c0 = poly_Ctp[0]; auto c1 = poly_Ctp[1]; auto c2 = poly_Ctp[2]; auto c3 = poly_Ctp[3];
			double F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12;
			F0 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 0);
			F1 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 1);
			F2 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 2);
			F3 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 3);
			F4 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 4);
			F5 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 5);
			F6 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 6);
			F7 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 7);
			F8 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 8);
			F9 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 9);
			F10 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 10);
			F11 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 11);
			F12 = F3_n(eta, poly_Ctp[0], poly_Ctp[1], poly_Ctp[2], poly_Ctp[3], 12);
			weights[v_id][(2 * todegree + 1) * i + 0] = (c0 - eta).dot(arthono(c1)) * F0 + 2 * (c0 - eta).dot(arthono(c2)) * F1
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F2 + 2 * c1.dot(arthono(c3)) * F3 + c2.dot(arthono(c3)) * F4;

			weights[v_id][(2 * todegree + 1) * i + 1] = (c0 - eta).dot(arthono(c1)) * F1 + 2 * (c0 - eta).dot(arthono(c2)) * F2
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F3 + 2 * c1.dot(arthono(c3)) * F4 + c2.dot(arthono(c3)) * F5;

			weights[v_id][(2 * todegree + 1) * i + 2] = (c0 - eta).dot(arthono(c1)) * F2 + 2 * (c0 - eta).dot(arthono(c2)) * F3
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F4 + 2 * c1.dot(arthono(c3)) * F5 + c2.dot(arthono(c3)) * F6;

			weights[v_id][(2 * todegree + 1) * i + 3] = (c0 - eta).dot(arthono(c1)) * F3 + 2 * (c0 - eta).dot(arthono(c2)) * F4
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F5 + 2 * c1.dot(arthono(c3)) * F6 + c2.dot(arthono(c3)) * F7;

			weights[v_id][(2 * todegree + 1) * i + 4] = (c0 - eta).dot(arthono(c1)) * F4 + 2 * (c0 - eta).dot(arthono(c2)) * F5
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F6 + 2 * c1.dot(arthono(c3)) * F7 + c2.dot(arthono(c3)) * F8;

			weights[v_id][(2 * todegree + 1) * i + 5] = (c0 - eta).dot(arthono(c1)) * F5 + 2 * (c0 - eta).dot(arthono(c2)) * F6
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F7 + 2 * c1.dot(arthono(c3)) * F8 + c2.dot(arthono(c3)) * F9;

			weights[v_id][(2 * todegree + 1) * i + 6] = (c0 - eta).dot(arthono(c1)) * F6 + 2 * (c0 - eta).dot(arthono(c2)) * F7
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F8 + 2 * c1.dot(arthono(c3)) * F9 + c2.dot(arthono(c3)) * F10;

			weights[v_id][(2 * todegree + 1) * i + 7] = (c0 - eta).dot(arthono(c1)) * F7 + 2 * (c0 - eta).dot(arthono(c2)) * F8
				+ (3 * (c0 - eta).dot(arthono(c3)) + 4 / 3 * (c1).dot(arthono(c2)) + 1 / 3 * c2.dot(arthono(c3))) * F9 + 2 * c1.dot(arthono(c3)) * F10 + c2.dot(arthono(c3)) * F11;

			weights[v_id][(2 * todegree + 1) * i + 8] = (c0 - eta).dot(c1) * F1 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F2 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F3
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F4 + 5 * c2.dot(c3) * F5 + 3 * c3.dot(c3) * F6 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][(2 * todegree + 1) * i + 9] = (c0 - eta).dot(c1) * F2 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F3 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F4
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F5 + 5 * c2.dot(c3) * F6 + 3 * c3.dot(c3) * F7 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][(2 * todegree + 1) * i + 10] = (c0 - eta).dot(c1) * F3 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F4 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F5
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F6 + 5 * c2.dot(c3) * F7 + 3 * c3.dot(c3) * F8 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][(2 * todegree + 1) * i + 11] = (c0 - eta).dot(c1) * F4 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F5 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F6
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F7 + 5 * c2.dot(c3) * F8 + 3 * c3.dot(c3) * F9 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][(2 * todegree + 1) * i + 12] = (c0 - eta).dot(c1) * F5 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F6 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F7
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F8 + 5 * c2.dot(c3) * F9 + 3 * c3.dot(c3) * F10 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][(2 * todegree + 1) * i + 13] = (c0 - eta).dot(c1) * F6 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F7 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F8
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F9 + 5 * c2.dot(c3) * F10 + 3 * c3.dot(c3) * F11 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);

			weights[v_id][(2 * todegree + 1) * i + 14] = (c0 - eta).dot(c1) * F7 + (c1.dot(c1) + 2 * (c0 - eta).dot(c2)) * F8 + (3 * (c1).dot(c2) + 3 * (c0 - eta).dot(c3)) * F9
				+ (2 * c2.dot(c2) + 4 * c1.dot(c3)) * F10 + 5 * c2.dot(c3) * F11 + 3 * c3.dot(c3) * F12 - log(norm(c0 + c1 + c2 + c3 - eta)) / (2 * M_PI);
			std::vector<std::vector<Mesh::Point>> bezier_ctp7;
			std::vector<std::vector<Mesh::Point>> poly_Ctp7;
			Bezier2Bezier7({ curvecage2[i] }, bezier_ctp7);
			Bezier2Poly7(bezier_ctp7, poly_Ctp7);
			test_eta += weights[v_id][(2 * todegree + 1) * i + 0] * poly_Ctp7[0][0];
			test_eta += weights[v_id][(2 * todegree + 1) * i + 1] * poly_Ctp7[0][1];
			test_eta += weights[v_id][(2 * todegree + 1) * i + 2] * poly_Ctp7[0][2];
			test_eta += weights[v_id][(2 * todegree + 1) * i + 3] * poly_Ctp7[0][3];
			test_eta += weights[v_id][(2 * todegree + 1) * i + 4] * poly_Ctp7[0][4];
			test_eta += weights[v_id][(2 * todegree + 1) * i + 5] * poly_Ctp7[0][5];
			test_eta += weights[v_id][(2 * todegree + 1) * i + 6] * poly_Ctp7[0][6];
			test_eta += weights[v_id][(2 * todegree + 1) * i + 7] * poly_Ctp7[0][7];
			test_eta += weights[v_id][(2 * todegree + 1) * i + 8] * arthono(poly_Ctp7[0][1]);
			test_eta += weights[v_id][(2 * todegree + 1) * i + 9] * arthono(poly_Ctp7[0][2]);
			test_eta += weights[v_id][(2 * todegree + 1) * i + 10] * arthono(poly_Ctp7[0][3]);
			test_eta += weights[v_id][(2 * todegree + 1) * i + 11] * arthono(poly_Ctp7[0][4]);
			test_eta += weights[v_id][(2 * todegree + 1) * i + 12] * arthono(poly_Ctp7[0][5]);
			test_eta += weights[v_id][(2 * todegree + 1) * i + 13] * arthono(poly_Ctp7[0][6]);
			test_eta += weights[v_id][(2 * todegree + 1) * i + 14] * arthono(poly_Ctp7[0][7]);
		}
		//check
		/*std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 0) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 1) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 2) << std::endl;
		std::cout << F2_n(Mesh::Point{ 0,0,0 }, Mesh::Point{ 1,0,0 }, Mesh::Point{ 0,1,0 }, Mesh::Point{ -0.5,1,0 }, 5) << std::endl;*/
		//std::cout << eta << endl;
		//std::cout << "check:" << test_eta - eta << std::endl;
		if ((test_eta - eta).norm() > max_err)
			max_err = (test_eta - eta).norm();
	}
	std::cout << "max norm err:" << max_err << std::endl;

}


//根据degree将CCpoints放进curvecage2里
void MeshViewerWidget::ccPoints2BezierSegs(void)
{
	for (int i = 0; i < curvecage2.size(); i++)
	{
		for (int j = 0; j < degree + 1; j++)
		{
			curvecage2[i][j] = CC_points[i * degree + j];
		}
	}
}

//根据curvecage2和坐标确定deformedmesh的顶点位置
void MeshViewerWidget::deform_mesh_from_cc(Mesh &deformedmesh)
{
	auto arthono = [](const Mesh::Point& p) -> Mesh::Point {
		return Mesh::Point(p[1], -p[0], p[2]);
	};
	auto curvecage2poly = curvecage2;
	if (usecvm) {
		//cal mvcGn
		int N = curvecage2.size();
		if (1) {
			std::vector<OpenMesh::Vec3d> t_e(N);
			std::vector<OpenMesh::Vec3d> n_e(N);
			for (int i = 0; i < N; i++)
			{
				mvcGt[2 * i] = 3 * (curvecage2[i][1] - curvecage2[i][0])/mvcL[i];
				mvcGt[2 * i + 1] = 3 * (curvecage2[i][2] - curvecage2[i][3]) / mvcL[i];
				n_e[i] = arthono((curvecage2[i][3] - curvecage2[i][0])).normalized();//每条边的单位法向
				t_e[i] = (curvecage2[i][3] - curvecage2[i][0]).normalized();//每条边的单位切向
			}
			for (int i = 0; i < N; i++)
			{
				int index = (i - 1) % N;
				if (index < 0) {
					index += N;
				}
				Eigen::Matrix2d A;
				Eigen::Matrix2d B;
				OpenMesh::Vec3d Bx = -mvcGt[2 * index + 1][0] * t_e[index] - mvcGt[2 * i][0] * t_e[i];
				OpenMesh::Vec3d By = -mvcGt[2 * index + 1][1] * t_e[index] - mvcGt[2 * i][1] * t_e[i];
				// 为矩阵A和B赋值
				A << n_e[i][0], -n_e[index][0],
					n_e[i][1], -n_e[index][1];
				B << Bx[0], By[0],
					Bx[1], By[1];
				//solve gn
				Eigen::Matrix2d X = A.inverse() * B;
				//std::cout << X << std::endl;
				mvcGn[2 * index + 1] = OpenMesh::Vec3d(X(1, 0), X(1, 1), 0);
				mvcGn[2 * i] = OpenMesh::Vec3d(X(0, 0), X(0, 1), 0);
				std::cout << "x的梯度" << mvcGt[2 * i][0] * t_e[i] + mvcGn[2 * i][0] * n_e[i] << std::endl;
				std::cout << "y的梯度" << mvcGt[2 * i][1] * t_e[i] + mvcGn[2 * i][1] * n_e[i] << std::endl;
				std::cout << mvcGt[2 * i][0] * t_e[i] + mvcGn[2 * i][0] * n_e[i] - (-mvcGt[2 * index + 1][0] * t_e[index] + mvcGn[2 * index + 1][0] * n_e[index]) << std::endl;
				std::cout << mvcGt[2 * i][1] * t_e[i] + mvcGn[2 * i][1] * n_e[i] - (-mvcGt[2 * index + 1][1] * t_e[index] + mvcGn[2 * index + 1][1] * n_e[index]) << std::endl;

			}
		}
		for (int i = 0; i < deformedmesh.n_vertices(); i++)
		{
			auto vh = deformedmesh.vertex_handle(i);
			auto point = mesh.point(vh);
			auto weight = weights[i];
			OpenMesh::Vec3d new_point(0, 0, point[2]);
			for (int j = 0; j < curvecage2.size(); j++)
			{
				new_point += curvecage2[j][0] * weight[5*j];
				new_point += mvcGt[2 * j] * weight[5 * j + 1];
				new_point += mvcGt[2 * j+1] * weight[5 * j +2];
				new_point += -mvcGn[2 * j] * weight[5 * j + 3];
				new_point += -mvcGn[2 * j + 1] * weight[5 * j + 4];
			}
			
			deformedmesh.set_point(vh, new_point);
		}
		return;
	}
	if (!highdegree) {
		if (degree == 2)
			Bezier2Poly2(curvecage2, curvecage2poly);
		else if (degree == 3)
			Bezier2Poly3(curvecage2, curvecage2poly);
		else if (degree == 1)
			Bezier1Poly1(curvecage2, curvecage2poly);
		std::vector < Mesh::Point> cpts(curvecage2poly.size() * (2 * degree + 1));
		for (int i = 0; i < curvecage2poly.size(); i++)
		{
			for (int j = 0; j < degree + 1; j++)
			{
				cpts[i * (2 * degree + 1) + j] = curvecage2poly[i][j];
			}
			for (int j = 1; j < degree + 1; j++)
			{
				cpts[i * (2 * degree + 1) + degree + j] = arthono(curvecage2poly[i][j]);
			}
		}
		for (int i = 0; i < deformedmesh.n_vertices(); i++)
		{
			auto vh = deformedmesh.vertex_handle(i);
			auto weight = weights[i];
			auto new_point = evaluate_coor(weight, cpts);
			deformedmesh.set_point(vh, new_point);
		}
	}
	else {
		if (todegree == 1) {
			Bezier1Poly1(curvecage2, curvecage2poly);
		}
		if (todegree == 2) {
			Bezier2Poly2(curvecage2, curvecage2poly);
		}
		if (todegree == 3) {
			Bezier2Poly3(curvecage2, curvecage2poly);
		}
		if (todegree == 7) {
			Bezier2Poly7(curvecage2, curvecage2poly);
		}
		std::vector < Mesh::Point> cpts(curvecage2poly.size() * (2 * todegree + 1));
		for (int i = 0; i < curvecage2poly.size(); i++)
		{
			for (int j = 0; j < todegree + 1; j++)
			{
				cpts[i * (2 * todegree + 1) + j] = curvecage2poly[i][j];
			}
			for (int j = 1; j < todegree + 1; j++)
			{
				cpts[i * (2 * todegree + 1) + todegree + j] = arthono(curvecage2poly[i][j]);
			}
		}
		for (int i = 0; i < deformedmesh.n_vertices(); i++)
		{
			auto vh = deformedmesh.vertex_handle(i);
			auto weight = weights[i];
			auto new_point = evaluate_coor(weight, cpts);
			deformedmesh.set_point(vh, new_point);
		}
	}
}

void MeshViewerWidget::Set_Texture_coord()
{
	LoadTexture();
	ptMin = Mesh::Point(-10, -11.26, 0);
	ptMax = Mesh::Point(10, 11.26, 0);//deer
	//ptMin = Mesh::Point(-6.4, -6.4, 0);
	//ptMax = Mesh::Point(6.4, 6.4, 0);//ball
	//ptMin = Mesh::Point(-14, -14, 0);
	//ptMax = Mesh::Point(14, 14, 0);//giraffe
	//ptMin = Mesh::Point(-1.98, -3.6, 0);
	//ptMax = Mesh::Point(1.98, 3.6, 0);//paristower
	//ptMin = Mesh::Point(-5.405, -2.14, 0);
	//ptMax = Mesh::Point(5.405, 2.14, 0);//xiyi
	//ptMin = Mesh::Point(-1.80, -2.285, 0);
	//ptMax = Mesh::Point(1.80, 2.285, 0);//yu
	//ptMin = Mesh::Point(-15.01, 16.6, 0);
	//ptMax = Mesh::Point(15.01, 16.6, 0);//kuzi
	//ptMin = Mesh::Point(-16.91, -6.11, 0);
	//ptMax = Mesh::Point(16.91, 6.11, 0);//car
	//ptMin = Mesh::Point(-1.21, -1.635, 0);
	//ptMax = Mesh::Point(1.21, 1.635, 0);//tuzi
	//ptMin = Mesh::Point(-32.13, 5.25, 0);
	//ptMax = Mesh::Point(32.13, 5.25, 0);//daizi
	//ptMin = Mesh::Point(-9.655,-3.92, 0);
	//ptMax = Mesh::Point(9.655, 3.92, 0);//fish
	//ptMin = Mesh::Point(0,0, 0);
	//ptMax = Mesh::Point(1, 1, 0);//zhengfangxin
	//ptMin = Mesh::Point(-5.625,-5.625, 0);
	//ptMax = Mesh::Point(5.625, 5.625, 0);//choufish
	//ptMin = Mesh::Point(-16.87, -16.87, 0);
	//ptMax = Mesh::Point(16.87, 16.87, 0);//zoomfish
	//ptMin = Mesh::Point(-6.15, -5.855, 0);
	//ptMax = Mesh::Point(6.15, 5.855, 0);//circle
	//ptMin = Mesh::Point(-6.23, -6.23, 0);
	//ptMax = Mesh::Point(6.23, 6.23, 0);//human
	//ptMin = Mesh::Point(-5.475, -9.6, 0);
	//ptMax = Mesh::Point(5.475, 9.6, 0);//flower
	//ptMin = Mesh::Point(-10, -15, 0);
	//ptMax = Mesh::Point(10, 15, 0);//octopus
	//ptMin = Mesh::Point(-16.86, -16.86, 0);
	//ptMax = Mesh::Point(16.86, 16.86, 0);//sheji
	//ptMin = Mesh::Point(-25.885, -25.885, 0);
	//ptMax = Mesh::Point(25.885, 25.885, 0);//pangxie
	float L = max({ ptMax[0] - ptMin[0],ptMax[1] - ptMin[1] });
	float Lx = 2 * ptMax[0];
	float Ly = 2 * ptMax[1];
	std::cout << "bbMIn:   " << ptMin[0] << " " << ptMin[1] << std::endl;
	std::cout << "bbMax:   " << ptMax[0] << " " << ptMax[1] << std::endl;
	/*ptMin = Mesh::Point(-10, -11.26, 0);
	ptMax = Mesh::Point(10, 11.26, 0);*/
	mesh.request_vertex_texcoords2D();//打开纹理坐标
	Mesh::Point P;
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		P = mesh.point(*v_it) - ptMin;
		Mesh::TexCoord2D t = { P[0] / Lx,P[1] / Ly };
		mesh.set_texcoord2D(*v_it, t);
	}

	std::cout << "load coordinate ok!" << std::endl;
}

void MeshViewerWidget::DrawPoints(void)
{
	auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(mesh, "vertexState");
	glColor3d(0.2, 0.2, 0.2);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (const auto& vh : mesh.vertices())
	{
		if (vertexState[vh] == NotSelected)
		{
			glNormal3dv(mesh.normal(vh).data());
			glVertex3dv(mesh.point(vh).data());
		}
	}
	glEnd();
	
	glColor3d(0.7, 0.7, 0.0);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (const auto& vh : mesh.vertices())
	{
		if (vertexState[vh] == Custom)
		{
			glNormal3dv(mesh.normal(vh).data());
			glVertex3dv(mesh.point(vh).data());
		}
	}
	glEnd();
}

void MeshViewerWidget::DrawWireframe(void)
{
	glColor3d(0.2, 0.2, 0.2);
	glBegin(GL_LINES);
	for (const auto& eh : mesh.edges())
	{
		auto heh = mesh.halfedge_handle(eh, 0);
		auto vh0 = mesh.from_vertex_handle(heh);
		auto vh1 = mesh.to_vertex_handle(heh);
		glNormal3dv(mesh.normal(vh0).data());
		glVertex3dv(mesh.point(vh0).data());
		glNormal3dv(mesh.normal(vh1).data());
		glVertex3dv(mesh.point(vh1).data());
	}
	glEnd();
	DrawPoints();
}

void MeshViewerWidget::DrawHiddenLines()
{
	glLineWidth(1.0);
	float backcolor[4];
	glGetFloatv(GL_COLOR_CLEAR_VALUE, backcolor);
	glColor4fv(backcolor);
	glDepthRange(0.01, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawFlat();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawFlat();
	}
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor3d(.3, .3, .3);
	DrawFlat();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::DrawFlatLines(void)
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.5f, 2.0f);
	glShadeModel(GL_FLAT);
	//glColor3d(0.8, 0.8, 0.8);
	glColor3d(1.0, 1.0, 1.0);
	DrawFlat();
	glDisable(GL_POLYGON_OFFSET_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawWireframe();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawWireframe();
	}
}

void MeshViewerWidget::DrawFlat(void)
{
	glBegin(GL_TRIANGLES);
	for (const auto& fh : mesh.faces())
	{
		glNormal3dv(mesh.normal(fh).data());
		for (const auto& fvh : mesh.fv_range(fh))
		{
			glVertex3dv(mesh.point(fvh).data());
		}
	}
	glEnd();
	DrawPoints();
}

void MeshViewerWidget::DrawSmooth(void)
{
	glColor3d(0.8, 0.8, 0.8);
	glShadeModel(GL_SMOOTH);
	glLoadName(static_cast<GLuint>(mesh.n_vertices()));
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());
	for (const auto& fh : mesh.faces())
	{
		glBegin(GL_POLYGON);
		for (const auto& fvh : mesh.fv_range(fh))
		{
			glArrayElement(fvh.idx());
		}
		glEnd();
	}
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	DrawPoints();
}

void MeshViewerWidget::DrawBezierCurve(const BezierCurve& bezierCurve, float r, float g, float b, float linewidth) {
	//std::cout << r << g << b << std::endl;
	if (&bezierCurve == nullptr)
	{
		return;
	}
	if ((!bezierCurve.GetCtrlPoints().empty())) {
		//DrawPoints3d(bezierCurve.GetCtrlPoints());
		
		linewidth = 3.0;
		glLineWidth(linewidth);
		glColor3f(r, g, b);
		glBegin(GL_LINE_STRIP);
		OpenMesh::Vec3d p;
		double t = 0;
		double h = 1.0 / 1000;
		int i = 0;
		while (i <= 1000)
		{
			p = bezierCurve.Evaluate(t + i++ * h);
			glVertex3f(p[0],p[1],p[2]);
		}
		glEnd();
	}
	else {
		//cout << "Empty Bezier Curve\n";
		return;
	}
}

void MeshViewerWidget::DrawTexture()
{

	glShadeModel(GL_SMOOTH);

	glBindTexture(GL_TEXTURE_2D, glTextureID);

	glEnable(GL_TEXTURE_2D);//开启2D纹理功能

	std::vector<Triangle> triangles;

	// Step 1: Extract triangles and calculate their center z-coordinate
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		Triangle tri;
		int i = 0;
		for (auto fv_it = mesh.fv_begin(*f_it); fv_it != mesh.fv_end(*f_it); ++fv_it, ++i)
		{
			tri.vertices[i] = mesh.point(*fv_it);
			tri.texcoords[i] = mesh.texcoord2D(*fv_it);
		}
		tri.zCenter = (tri.vertices[0][2] + tri.vertices[1][2] + tri.vertices[2][2]) / 3.0f;
		triangles.push_back(tri);
	}

	// Step 2: Sort triangles based on their zCenter (back to front)
	std::sort(triangles.begin(), triangles.end(), [](const Triangle& a, const Triangle& b) {
		return a.zCenter > b.zCenter;
		});

	// Step 3: Draw sorted triangles
	glBegin(GL_TRIANGLES);
	for (const auto& tri : triangles)
	{
		for (int i = 0; i < 3; ++i)
		{
			float uvScale = 1.0f; // or 3.0f depending on the use case
			glTexCoord2f(tri.texcoords[i][0] * uvScale, tri.texcoords[i][1] * uvScale);
			glVertex3fv(tri.vertices[i].data());
		}
	}
	glEnd();

	glDisable(GL_TEXTURE_2D);

	//glBegin(GL_TRIANGLES);
	//{
	//	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	//	{
	//		for (auto fv_it = mesh.fv_begin(*f_it); fv_it != mesh.fv_end(*f_it); fv_it++)
	//		{
	//			auto uv = mesh.texcoord2D(*fv_it);
	//			float uvScale = 1.f;//怪兽的例子改成1.，其他改3.
	//			glTexCoord2f(uv[0] * uvScale, uv[1] * uvScale);
	//			glVertex3dv(mesh.point(*fv_it).data());

	//		}

	//	}


	//}
	//


	//glEnd();
	//glDisable(GL_TEXTURE_2D);//关闭纹理贴图功能
}

void MeshViewerWidget::DrawCageWireframe(void)
{
	glColor3d(0.2, 0.2, 0.2);
	glBegin(GL_LINES);
	// 获取所有顶点句柄
	std::vector<Mesh::VertexHandle> vhandles;
	for (const auto& vh : CC_mesh.vertices()) {
		vhandles.push_back(vh);
	}

	// 按顺序连接相邻的顶点
	for (size_t i = 0; i < vhandles.size(); ++i) {
		auto vh0 = vhandles[i];
		auto vh1 = vhandles[(i + 1) % vhandles.size()];  // 循环连接最后一个顶点和第一个顶点
		glVertex3dv(CC_mesh.point(vh0).data());
		glVertex3dv(CC_mesh.point(vh1).data());
	}
	glEnd();
	
}

void MeshViewerWidget::DrawCagePoints(void)
{	
	auto CC_vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(CC_mesh, "vertexState");
	glColor3d(0.2, 0.2, 0.2);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (const auto& vh : CC_mesh.vertices())
	{
		if (CC_vertexState[vh] == NotSelected)
		{
			
			glVertex3dv(CC_mesh.point(vh).data());
		}
	}
	glEnd();
	
	glColor3d(1.0, 0.0, 0.0);
	glPointSize(50);
	glBegin(GL_POINTS);
	for (const auto& vh : CC_mesh.vertices())
	{
		if (CC_vertexState[vh] == Custom)
		{
			
			glVertex3dv(CC_mesh.point(vh).data());
		}
	}
	glEnd();
}

void MeshViewerWidget::DrawCurveCage(void)
{
	int N = curvecage2.size();
	std::cout << N << std::endl;
	for (int i = 0; i < N; i++) {
		DrawBezierCurve(curvecage2[i],0,0,1);
		//DrawPoints3d(curvecage2[i], 0, 1, 0);
	}
	
	//DrawFlat();
	//DrawPoints();
	//cage的边和点
	if(drawpoints)
		DrawCagePoints();
	//DrawCageWireframe();
	//DrawCagePoints();
	glColor3d(1.0, 1.0, 1.0);
	DrawTexture();
	
	//DrawFlat();
}

void MeshViewerWidget::DrawBoundingBox(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(.3, .7, .3);
	glBegin(GL_LINES);
	for (const auto& i : { 0, 1 })
	{
		for (const auto& j : { 0, 1 })
		{
			for (const auto& k : { 0, 1 })
			{
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(~i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], ~j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], ~k ? ptMin[2] : ptMax[2]);
			}
		}
	}
	glEnd();
	glLineWidth(linewidth);
}

void MeshViewerWidget::DrawBoundary(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(0.1, 0.1, 0.1);
	glBegin(GL_LINES);
	for (const auto& eh : mesh.edges())
	{
		if (mesh.is_boundary(eh))
		{
			auto heh = mesh.halfedge_handle(eh, 0);
			auto vh0 = mesh.from_vertex_handle(heh);
			auto vh1 = mesh.to_vertex_handle(heh);
			glNormal3dv(mesh.normal(vh0).data());
			glVertex3dv(mesh.point(vh0).data());
			glNormal3dv(mesh.normal(vh1).data());
			glVertex3dv(mesh.point(vh1).data());
		}
	}
	glEnd();
	glLineWidth(linewidth);
}
