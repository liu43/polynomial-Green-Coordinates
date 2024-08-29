#pragma once

#include <QWidget>
#include <QtGui>
#include <QtWidgets>

class MeshParamWidget : public QWidget
{
	Q_OBJECT

public:
	MeshParamWidget(QWidget *parent = 0);
	~MeshParamWidget(void);
private:
	void CreateTabWidget(void);
	void CreateLayout(void);
signals:
	void PrintInfoSignal();
	void NoSelectSignal();
	void SelectAdjustSignal();
	void SelectCustomSignal();
	void MoveSignal();
	void UpdegreeSignal();
	void ChangedegreeSignal();
	void NodrawpointSignal();
	void AddpointsSignal();
	void ClearSignal();
private slots:
	void ResetCheck();
private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;
	QPushButton* noSelectBtn;
	QPushButton* SelectAdjustBtn;
	QPushButton* selectCustomBtn;
	QPushButton* moveBtn;
	QPushButton* clearBtn;
	QPushButton* updegreeBtn;
	QPushButton* changedegreeBtn;
	QPushButton* nodrawpointBtn;
	QPushButton* addpointsBtn;
	QButtonGroup* deformBtnGroup;
	QWidget* deformWidget;
};
