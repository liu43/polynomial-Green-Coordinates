#pragma once
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <QOpenGLWidget>
#include <GL/GLU.h>
class QOpenGLTexture;

class QGLViewerWidget : public QOpenGLWidget
{
	Q_OBJECT
public:
	QGLViewerWidget(QWidget* _parent = 0);
	virtual ~QGLViewerWidget(void);
private:
	void Init(void);
public:
	QSize minimumSizeHint(void) const override;
	QSize sizeHint(void) const override;
	const double & Radius(void) const;
	const OpenMesh::Vec3d & Center(void) const;
	const double* GetModelviewMatrix(void) const;
	void ResetModelviewMatrix(void);
	void CopyModelViewMatrix(void);
	void LoadCopyModelViewMatrix(void);
	const double* GetProjectionMatrix(void) const;

	enum ProjectionMode{ PERSPECTIVE, ORTHOGRAPHIC };
	void SetProjectionMode(const ProjectionMode &pm);
	const ProjectionMode & GetProjectionMode(void) const;

	enum DrawMode{ POINTS, WIREFRAME, HIDDENLINES, FLATLINES, FLAT, SMOOTH, CURVECAGE };
	void SetDrawMode(const DrawMode &dm);
	const DrawMode& GetDrawMode(void) const;

protected:
	enum MaterialType { MaterialDefault, MaterialGold, MaterialSilver, MaterialEmerald, MaterialTin };
	void SetMaterial(const MaterialType & mattype = MaterialGold) const;
	void SetDefaultLight(void) const;
	void initializeGL(void) override;
	void resizeGL(int w, int h) override;
	void paintGL(void) override;
	virtual void DrawScene(void);
	virtual void mousePressEvent(QMouseEvent*) override;
	virtual void mouseMoveEvent(QMouseEvent*) override;
	virtual void mouseReleaseEvent(QMouseEvent*) override;
	virtual void wheelEvent(QWheelEvent*) override;
	virtual void keyPressEvent(QKeyEvent*) override;
	virtual void keyReleaseEvent(QKeyEvent*) override;
private:
	void Translation(const QPoint & p);
	void Translate(const OpenMesh::Vec3d & trans);
	void Rotation(const QPoint & p);
	void Rotate(const OpenMesh::Vec3d & axis, const double & angle);
	bool MapToSphere(const QPoint & point, OpenMesh::Vec3d & result);
	void UpdateProjectionMatrix(void);
	
public:
	void SetScenePosition(const OpenMesh::Vec3d & c, const double & r);
	void ViewAll(void);
	void WinCor2ObjCor(int x, int y, OpenMesh::Vec3d& objCor, double& depth);
	void WinCor2ObjCor(int x, int y, double givenDepth, OpenMesh::Vec3d& objCor);
	// Load uv texture
	void LoadTexture();
protected:
	DrawMode drawmode;
	ProjectionMode projectionmode;
	double windowleft;
	double windowright;
	double windowtop;
	double windowbottom;
	Qt::MouseButton mousemode;
	OpenMesh::Vec3d center;
	double radius;
	std::vector<double> projectionmatrix;
	std::vector<double> modelviewmatrix;
	std::vector<double> copymodelviewmatrix;
	QPoint lastpoint2;
	OpenMesh::Vec3d lastpoint3;
	bool lastpointok;
	//texture by hzb
	unsigned int glTextureID;
private:
	static const double trackballradius;
};
