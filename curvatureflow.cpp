
#include "willmoreflow.h"
bool pre_draw(igl::viewer::Viewer &viewer)
{
	if (!viewer.core.is_animating)
	{
		return false;
	}
	else {
		if (boolVariable)
		{
			W.Schmid();
			W.constraint();
			W.integrate();

			viewer.data.clear();
			viewer.data.set_mesh(V, T);
			viewer.data.set_edges(V, EV, Eigen::RowVector3d(1, 1, 1));
			return false;
		}
		else if (boolg)
		{
			W.mycomputegradient();
			W.inte();
			viewer.data.clear();
			viewer.data.set_mesh(V, T);
			viewer.data.set_edges(V, EV, Eigen::RowVector3d(1, 1, 1));
			return false;
		}
	}
	return false;
}
int main(int argc, char *argv[])
{
	// Load a mesh in OFF format
	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, "e:/hedra/libhedra/examples/data/test.off", opt))
	{
		cout << "无法读取文件:" << "e: / hedra / libhedrxa / examples / data / test.off" << endl;
		return 1;
	}

	polygonal_read_OFF("e:/hedra/libhedra/examples/data/test.off", V, D, F);
	triangulate_mesh(D, F, T, TF);
	polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
	W.defaul();
	W.init();
	cout << "sign: " << Sgn(1230.3) << endl;
	viewer.callback_init = [&](igl::viewer::Viewer &viewer)
	{
		viewer.ngui->addButton("iter_gradient", [&]() {
			W.computegradient();
			W.inte();
			viewer.data.clear();
			viewer.data.set_mesh(V, T);
			viewer.data.set_edges(V, EV, Eigen::RowVector3d(1, 1, 1));
		});
		viewer.ngui->addButton("iter_k", [&]() {W.Schmid();
		W.constraint();
		W.integrate();
		viewer.data.clear();
		viewer.data.set_mesh(V, T);
		viewer.data.set_edges(V, EV, Eigen::RowVector3d(1, 1, 1));
		});

		viewer.ngui->addVariable("t", t);
		viewer.ngui->addVariable("energy", energy);
		//viewer.ngui->
		viewer.ngui->addVariable<bool>("boolk", [&](bool val) {boolVariable = val; viewer.core.is_animating = boolVariable;
		}, [&]() {return boolVariable; });
		viewer.ngui->addVariable<bool>("boolg", [&](bool val) {boolg = val; viewer.core.is_animating = boolg;
		}, [&]() {return boolg; });

		viewer.screen->performLayout();
		return false;
	};
	viewer.callback_pre_draw = &pre_draw;
	viewer.data.clear();
	viewer.data.set_mesh(V, T);
	viewer.data.set_edges(V, EV, Eigen::RowVector3d(1, 0, 0));

	viewer.launch();
}
