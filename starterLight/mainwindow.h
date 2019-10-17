#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QMessageBox>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();


    float compute_area(MyMesh *_mesh);
    float get_min_area(MyMesh *_mesh);
    float get_max_area(MyMesh *_mesh);
    float compute_face_area(MyMesh *_mesh, int n_face);
    void displayMesh(MyMesh *_mesh);
    void resetAllColorsAndThickness(MyMesh *_mesh);
    MyMesh::Point getNormalPoint (MyMesh* _mesh,VertexHandle vertexFromFace);
    bool checkGlobalNeighbours(MyMesh* _mesh);
    bool checkAllTriangularFace(MyMesh* _mesh);
    void normals_points(MyMesh * _mesh);
    bool checkOnlyPoint();
    float angle_vector(MyMesh::Point vecteur1, MyMesh::Point vecteur2);
    float moy_angle_vertice_faces(MyMesh * _mesh, VertexHandle v);
    void angles_normal_points(MyMesh * _mesh);
    void putVertexColor(QVector <float> angles, float degresMin, float degresMax);
    MyMesh::Point getBarycenterFromFace(VertexHandle vh ,FaceHandle fh,MyMesh* _mesh);
    MyMesh::Point getNormalFace (MyMesh* _mesh,VertexHandle v0,VertexHandle v1, VertexHandle v2);
    //std::vector<MyMesh::Point> getNormalFace (MyMesh* _mesh,VertexHandle v1, VertexHandle v2);// getNormalFace (MyMesh* _mesh,VertexHandle vertexFromFace, float barycentre);

    /**
     * @brief dihedralAngles
     * @param _mesh
     * @details Pour chaque arête calculer l'angle entre les 2 normales des 2 faces concourantes.
     */
    void dihedralAngles(MyMesh* _mesh);

private slots:

    void on_pushButton_chargement_clicked();
    void on_pushButton_barycentre_clicked();
    void on_getInformation_clicked();
    void on_boundingBox_clicked();
    void on_pushButton_area_clicked();
    void on_triangleSurface_proportion_clicked();
    void on_meshIsValid_clicked();
    void on_getValanceRing_clicked();
    void on_show_pts_norm_clicked();
    void on_pushButton_fv_angle_clicked();

    void on_getDiedreAngle_clicked();

    void on_pushButton_save_clicked();

private:

    bool modevoisinage;

    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
