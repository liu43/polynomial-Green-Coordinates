/****************************************************************************
** Meta object code from reading C++ file 'MeshParamWidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.12.9)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../MeshParamWidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshParamWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.12.9. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_MeshParamWidget_t {
    QByteArrayData data[13];
    char stringdata0[188];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MeshParamWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MeshParamWidget_t qt_meta_stringdata_MeshParamWidget = {
    {
QT_MOC_LITERAL(0, 0, 15), // "MeshParamWidget"
QT_MOC_LITERAL(1, 16, 15), // "PrintInfoSignal"
QT_MOC_LITERAL(2, 32, 0), // ""
QT_MOC_LITERAL(3, 33, 14), // "NoSelectSignal"
QT_MOC_LITERAL(4, 48, 18), // "SelectAdjustSignal"
QT_MOC_LITERAL(5, 67, 18), // "SelectCustomSignal"
QT_MOC_LITERAL(6, 86, 10), // "MoveSignal"
QT_MOC_LITERAL(7, 97, 14), // "UpdegreeSignal"
QT_MOC_LITERAL(8, 112, 18), // "ChangedegreeSignal"
QT_MOC_LITERAL(9, 131, 17), // "NodrawpointSignal"
QT_MOC_LITERAL(10, 149, 15), // "AddpointsSignal"
QT_MOC_LITERAL(11, 165, 11), // "ClearSignal"
QT_MOC_LITERAL(12, 177, 10) // "ResetCheck"

    },
    "MeshParamWidget\0PrintInfoSignal\0\0"
    "NoSelectSignal\0SelectAdjustSignal\0"
    "SelectCustomSignal\0MoveSignal\0"
    "UpdegreeSignal\0ChangedegreeSignal\0"
    "NodrawpointSignal\0AddpointsSignal\0"
    "ClearSignal\0ResetCheck"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MeshParamWidget[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      11,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
      10,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   69,    2, 0x06 /* Public */,
       3,    0,   70,    2, 0x06 /* Public */,
       4,    0,   71,    2, 0x06 /* Public */,
       5,    0,   72,    2, 0x06 /* Public */,
       6,    0,   73,    2, 0x06 /* Public */,
       7,    0,   74,    2, 0x06 /* Public */,
       8,    0,   75,    2, 0x06 /* Public */,
       9,    0,   76,    2, 0x06 /* Public */,
      10,    0,   77,    2, 0x06 /* Public */,
      11,    0,   78,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      12,    0,   79,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,

       0        // eod
};

void MeshParamWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<MeshParamWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->PrintInfoSignal(); break;
        case 1: _t->NoSelectSignal(); break;
        case 2: _t->SelectAdjustSignal(); break;
        case 3: _t->SelectCustomSignal(); break;
        case 4: _t->MoveSignal(); break;
        case 5: _t->UpdegreeSignal(); break;
        case 6: _t->ChangedegreeSignal(); break;
        case 7: _t->NodrawpointSignal(); break;
        case 8: _t->AddpointsSignal(); break;
        case 9: _t->ClearSignal(); break;
        case 10: _t->ResetCheck(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (MeshParamWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MeshParamWidget::PrintInfoSignal)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (MeshParamWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MeshParamWidget::NoSelectSignal)) {
                *result = 1;
                return;
            }
        }
        {
            using _t = void (MeshParamWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MeshParamWidget::SelectAdjustSignal)) {
                *result = 2;
                return;
            }
        }
        {
            using _t = void (MeshParamWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MeshParamWidget::SelectCustomSignal)) {
                *result = 3;
                return;
            }
        }
        {
            using _t = void (MeshParamWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MeshParamWidget::MoveSignal)) {
                *result = 4;
                return;
            }
        }
        {
            using _t = void (MeshParamWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MeshParamWidget::UpdegreeSignal)) {
                *result = 5;
                return;
            }
        }
        {
            using _t = void (MeshParamWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MeshParamWidget::ChangedegreeSignal)) {
                *result = 6;
                return;
            }
        }
        {
            using _t = void (MeshParamWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MeshParamWidget::NodrawpointSignal)) {
                *result = 7;
                return;
            }
        }
        {
            using _t = void (MeshParamWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MeshParamWidget::AddpointsSignal)) {
                *result = 8;
                return;
            }
        }
        {
            using _t = void (MeshParamWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MeshParamWidget::ClearSignal)) {
                *result = 9;
                return;
            }
        }
    }
    Q_UNUSED(_a);
}

QT_INIT_METAOBJECT const QMetaObject MeshParamWidget::staticMetaObject = { {
    &QWidget::staticMetaObject,
    qt_meta_stringdata_MeshParamWidget.data,
    qt_meta_data_MeshParamWidget,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *MeshParamWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MeshParamWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MeshParamWidget.stringdata0))
        return static_cast<void*>(this);
    return QWidget::qt_metacast(_clname);
}

int MeshParamWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 11)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 11;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 11)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 11;
    }
    return _id;
}

// SIGNAL 0
void MeshParamWidget::PrintInfoSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 0, nullptr);
}

// SIGNAL 1
void MeshParamWidget::NoSelectSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}

// SIGNAL 2
void MeshParamWidget::SelectAdjustSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 2, nullptr);
}

// SIGNAL 3
void MeshParamWidget::SelectCustomSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 3, nullptr);
}

// SIGNAL 4
void MeshParamWidget::MoveSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 4, nullptr);
}

// SIGNAL 5
void MeshParamWidget::UpdegreeSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 5, nullptr);
}

// SIGNAL 6
void MeshParamWidget::ChangedegreeSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 6, nullptr);
}

// SIGNAL 7
void MeshParamWidget::NodrawpointSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 7, nullptr);
}

// SIGNAL 8
void MeshParamWidget::AddpointsSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 8, nullptr);
}

// SIGNAL 9
void MeshParamWidget::ClearSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 9, nullptr);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
