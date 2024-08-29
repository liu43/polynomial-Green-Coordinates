#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateLayout();
}

MeshParamWidget::~MeshParamWidget()
{
}

void MeshParamWidget::CreateTabWidget(void)
{
	pbPrintInfo = new QPushButton(tr("Cal Weight"));
	connect(pbPrintInfo, SIGNAL(clicked()), SIGNAL(PrintInfoSignal()));

	noSelectBtn = new QPushButton(tr("No Select"));
	SelectAdjustBtn = new QPushButton(tr("Adjust Cage"));
	selectCustomBtn = new QPushButton(tr("Select Custom"));
	moveBtn = new QPushButton(tr("Move"));
	clearBtn = new QPushButton(tr("Deform"));
	updegreeBtn = new QPushButton(tr("Updegree"));
	changedegreeBtn = new QPushButton(tr("Changedegree"));
	nodrawpointBtn = new QPushButton(tr("No drawpoint"));
	addpointsBtn = new QPushButton(tr("Add points"));

	
	connect(SelectAdjustBtn, SIGNAL(clicked()), SIGNAL(SelectAdjustSignal()));
	connect(selectCustomBtn, SIGNAL(clicked()), SIGNAL(SelectCustomSignal()));
	connect(moveBtn, SIGNAL(clicked()), SIGNAL(MoveSignal()));
	connect(clearBtn, SIGNAL(clicked()), SIGNAL(ClearSignal()));
	connect(clearBtn, SIGNAL(clicked()), SLOT(ResetCheck()));
	connect(updegreeBtn, SIGNAL(clicked()), SIGNAL(UpdegreeSignal()));
	connect(changedegreeBtn, SIGNAL(clicked()), SIGNAL(ChangedegreeSignal()));
	connect(nodrawpointBtn, SIGNAL(clicked()), SIGNAL(NodrawpointSignal()));
	connect(addpointsBtn, SIGNAL(clicked()), SIGNAL(AddpointsSignal()));

	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(pbPrintInfo);
	layout->addWidget(noSelectBtn);
	layout->addWidget(SelectAdjustBtn);
	layout->addWidget(selectCustomBtn);
	layout->addWidget(moveBtn);
	layout->addWidget(clearBtn);
	layout->addWidget(updegreeBtn);
	layout->addWidget(changedegreeBtn);
	layout->addWidget(nodrawpointBtn);
	layout->addWidget(addpointsBtn);
	layout->addStretch();
	wParam = new QWidget();
	wParam->setLayout(layout);
	saParam = new QScrollArea();
	saParam->setFocusPolicy(Qt::NoFocus);
	saParam->setFrameStyle(QFrame::NoFrame);
	saParam->setWidget(wParam);
	saParam->setWidgetResizable(true);
}


void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	
	this->setLayout(layout);
}


void MeshParamWidget::ResetCheck()
{
	noSelectBtn->setChecked(true);
	SelectAdjustBtn->setChecked(false);
	selectCustomBtn->setChecked(false);
	moveBtn->setChecked(false);
	updegreeBtn->setChecked(false);
	changedegreeBtn->setChecked(false);
	nodrawpointBtn->setChecked(false);
	addpointsBtn->setChecked(false);
}
