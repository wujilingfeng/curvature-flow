
#include<igl\viewer\Viewer.h>
#include <igl/readOFF.h>
#include<igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include<OpenMesh\Core\IO\IOManager.hh>
#include<OpenMesh\Core\IO\MeshIO.hh>
#include<OpenMesh\Core\Mesh\PolyMesh_ArrayKernelT.hh>
#include<OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>
#include<nanogui\screen.h>
#include<nanogui\formhelper.h>
#include<nanogui.h>
#include <hedra/polygonal_edge_topology.h>
#include<hedra\polygonal_read_OFF.h>
#include<hedra\triangulate_mesh.h>
#include<iostream>
using namespace std;
using namespace Eigen;
using namespace hedra;
using namespace OpenMesh;
#ifndef Sgn
#define Sgn(y)(y>0?1:-1)
#endif

struct MyTraits :OpenMesh::DefaultTraits
{
	VertexTraits
	{
	private:
		double dual_length;
		double angle;
		double energy;
		double area;
		double curvature;
		double gradient;

		RowVector3d  gradientenergy;
	public:
		VertexT()
		{
			angle = 0;
			dual_length = 0;
			energy = 0;
			area = 0;
			curvature = 0;
			gradient = 0;
		}
		double &dl_get()
		{
			return dual_length;
		}
		double &angle_get()
		{
			return angle;
		}
		double &energy_get()
		{
			return energy;
		}
		double &area_get()
		{
			return area;
		}
		RowVector3d &gradiente_get()
		{
			return gradientenergy;

		}
		double &curvature_get()
		{
			return curvature;
		}
		double &gradient_get()
		{
			return gradient;
		}

	};

};
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> PolyMesh;
igl::viewer::Viewer viewer;
Eigen::MatrixXd V;
VectorXi D, innerEdges, TF;
Eigen::MatrixXi F, EF, EV, FE, EFi, T;
MatrixXd FEs;
double t = 0.001, energy = 0;
bool boolVariable = false, boolg = false;
class Willmoreflow
{
public:
	Willmoreflow(PolyMesh &mesh_) :mesh(mesh_)
	{
	}
	~Willmoreflow() {}
	void defaul()
	{
		int n = mesh.n_vertices();
		A.resize(n, n); k.resize(n, 1); g.resize(n, 3);
		A.setZero(); k.setZero(); g.setZero();
		C[0].resize(n, 1); C[1].resize(n, 1); C[2].resize(n, 1);
		C[0].setOnes(); C[1].setOnes(); C[2].setOnes();

	}
	void init()
	{
		HalfedgeHandle he;
		VertexHandle v;
		double  len = 0; k.setZero();
		for (PolyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			he = mesh.halfedge_handle(*vter);
			v = mesh.to_vertex_handle(he);
			mesh.point(v)[0] = V.coeff(v.idx(), 0);
			mesh.point(v)[1] = V.coeff(v.idx(), 1);
			mesh.point(v)[2] = V.coeff(v.idx(), 2);
			len = mesh.calc_edge_length(mesh.edge_handle(he));
			mesh.data(v).angle_get() = Sgn(mesh.calc_sector_angle(he))*M_PI - mesh.calc_sector_angle(he);
			mesh.data(v).area_get() = mesh.calc_sector_area(he)*Sgn(sin(mesh.data(v).angle_get()));
			//cout << "mesh.data(v).angle_get(): "<< mesh.data(v).angle_get() << endl;
			//cout << "mesh.calc_sector_area(he): " << mesh.data(v).area_get() << endl;

			len += mesh.calc_edge_length(mesh.edge_handle(mesh.next_halfedge_handle(he)));
			mesh.data(v).dl_get() = len / 2.0;
			mesh.data(v).curvature_get() = mesh.data(v).angle_get() / mesh.data(v).dl_get();
			mesh.data(v).gradient_get() = -2 * mesh.data(v).curvature_get();
			mesh.data(v).gradiente_get().setZero();
			A.coeffRef(v.idx(), v.idx()) = mesh.data(v).dl_get();
			k.coeffRef(v.idx(), 0) = mesh.data(v).gradient_get();
			g.setZero();
			C[1].coeffRef(v.idx(), 0) = V.coeff(v.idx(), 0);
			C[2].coeffRef(v.idx(), 0) = V.coeff(v.idx(), 1);
		}
	}
	void Schmid()
	{
		double temp1 = (C[1].transpose()*A*C[0]).coeff(0, 0) / (C[0].transpose()*A*C[0]).coeff(0, 0);
		C[1] = C[1] - temp1 * C[0];
		double temp2 = (C[2].transpose()*A*C[0]).coeff(0, 0) / (C[0].transpose()*A*C[0]).coeff(0, 0), temp3 = (C[2].transpose() * A*C[1]).coeff(0, 0) / (C[1].transpose() * A*C[1]).coeff(0, 0);
		C[2] = C[2] - temp2 * C[0] - temp3 * C[1];
		C[0] /= sqrt((C[0].transpose()*A*C[0]).coeff(0, 0));
		C[1] /= sqrt((C[1].transpose()*A*C[1]).coeff(0, 0));
		C[2] /= sqrt((C[2].transpose()*A*C[2]).coeff(0, 0));
	}
	void constraint()
	{
		for (int i = 0; i<3; i++)
		{
			k -= (k.transpose()*A*C[i]).coeff(0, 0)*C[i];
		}
	}
	void integrate()
	{
		int n = mesh.n_vertices();
		normalk();
		VectorXd pha = t * A * k;
		MatrixXd T(n, 3);
		T.setZero();
		PolyMesh::VertexIter vter = mesh.vertices_begin();
		VertexHandle v = *vter;
		HalfedgeHandle he = mesh.halfedge_handle(v);
		RowVector3d p, initp = V.row(v.idx()), tempp;
		tempp.setZero();
		double temp = 0;
		do
		{
			temp += pha.coeff(mesh.from_vertex_handle(he).idx());
			p = V.row(mesh.to_vertex_handle(he).idx()) - V.row(mesh.from_vertex_handle(he).idx());
			T.coeffRef(mesh.to_vertex_handle(he).idx(), 0) = (p[0] * cos(temp) - p[1] * sin(temp));
			T.coeffRef(mesh.to_vertex_handle(he).idx(), 1) = (p[1] * cos(temp) + p[0] * sin(temp));
			he = mesh.next_halfedge_handle(he);
		} while (mesh.from_vertex_handle(he) != v);

		do
		{
			tempp += T.row(mesh.to_vertex_handle(he).idx());
			V.row(mesh.to_vertex_handle(he).idx()) = (initp + tempp);
			he = mesh.next_halfedge_handle(he);
		} while (mesh.from_vertex_handle(he) != v);

		init();
		init();
		energy = (k.transpose()*A*k).coeff(0, 0);
		cout << energy << endl;
	}
	void computegradient()
	{
		PolyMesh::VertexIter vter = mesh.vertices_begin();

		VertexHandle v = *vter;

		HalfedgeHandle he = mesh.halfedge_handle(v);
		if (mesh.from_vertex_handle(he) != v)
		{
			cout << "cuowu" << endl;
			return;
		}
		RowVector3d U, VV;
		do
		{
			VertexHandle vv = mesh.from_vertex_handle(he), v1 = mesh.to_vertex_handle(he), v_1 = mesh.from_vertex_handle(mesh.prev_halfedge_handle(he));
			U = V.row(vv.idx()) - V.row(v_1.idx());
			VV = V.row(v1.idx()) - V.row(vv.idx());
			RowVector3d Uhat = U / U.norm(), Vhat = VV / VV.norm();
			RowVector3d Vper = (VV - VV.dot(Uhat)*Uhat), Uper = (U - Vhat.dot(U)*Vhat);
			double phi = mesh.data(vv).angle_get(), area = mesh.data(vv).area_get();
			mesh.data(v_1).gradiente_get() += phi / (mesh.data(v_1).dl_get() * mesh.data(vv).dl_get())*(Vper / area + phi * 0.5 / mesh.data(vv).dl_get()*Uhat);
			mesh.data(vv).gradiente_get() -= phi / (mesh.data(v1).dl_get() *mesh.data(vv).dl_get())*(Uper / area + phi * 0.5 / mesh.data(vv).dl_get()*Vhat);
			mesh.data(v1).gradiente_get() += phi / pow(mesh.data(vv).dl_get(), 2)*((Uper - Vper) / area + phi * 0.5 / mesh.data(vv).dl_get()*(Vhat - Uhat));
			he = mesh.next_halfedge_handle(he);
		} while (mesh.from_vertex_handle(he) != v);

	}
	void mycomputegradient()
	{
		PolyMesh::VertexIter vter = mesh.vertices_begin();
		VertexHandle v = *vter;
		HalfedgeHandle he = mesh.halfedge_handle(v);
		if (mesh.from_vertex_handle(he) != v)
		{
			cout << "cuowu" << endl;
			return;
		}
		RowVector3d U, VV;
		do
		{
			VertexHandle vv = mesh.from_vertex_handle(he), v1 = mesh.to_vertex_handle(he), v_1 = mesh.from_vertex_handle(mesh.prev_halfedge_handle(he));
			U = V.row(vv.idx()) - V.row(v_1.idx());
			VV = V.row(v1.idx()) - V.row(vv.idx());
			RowVector3d Uhat = U / U.norm(), Vhat = VV / VV.norm();
			RowVector3d Vper = (VV - VV.dot(Uhat)*Uhat), Uper = (U - Vhat.dot(U)*Vhat);
			double phi = fabs(mesh.data(vv).angle_get()), area = fabs(mesh.data(vv).area_get());
			mesh.data(v_1).gradiente_get() -= phi / mesh.data(vv).dl_get()*(Vper / area + phi * 0.5 / mesh.data(vv).dl_get()*Uhat);
			mesh.data(vv).gradiente_get() += phi / mesh.data(vv).dl_get()*((Vper - Uper) / area + phi * 0.5 / mesh.data(vv).dl_get()*(Uhat - Vhat));
			mesh.data(v1).gradiente_get() += phi / mesh.data(vv).dl_get()*((Uper) / area + phi * 0.5 / mesh.data(vv).dl_get()*(Vhat));
			he = mesh.next_halfedge_handle(he);
		} while (mesh.from_vertex_handle(he) != v);

	}
	void inte()
	{
		for (PolyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			g.row((*vter).idx()) = mesh.data(*vter).gradiente_get();
		}
		normalg();
		for (PolyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			V.row((*vter).idx()) += (t * g.row((*vter).idx()));

			//cout << mesh.data((*vter)).gradiente_get() << endl;
		}
		init();
		init();
	}
	void normalk()
	{
		double temp1 = k.maxCoeff(), temp2 = k.minCoeff();
		if (fabs(temp1)>fabs(temp2))
		{
			k /= fabs(temp1);
		}
		else
		{
			k /= fabs(temp2);
		}

	}
	void normalg()
	{

		double temp1 = g.maxCoeff(), temp2 = g.minCoeff();
		if (fabs(temp1)>fabs(temp2))
		{
			g /= fabs(temp1);
		}
		else
		{
			g /= fabs(temp2);
		}
	}
	void test()
	{
		//cout << k << endl;
		cout << "test: " << k.transpose()*A * C[1] << endl;
	}
private:
	MatrixXd A;
	VectorXd C[3], k;
	MatrixXd g;
	PolyMesh &mesh;
};
PolyMesh mesh;
Willmoreflow W(mesh);
