#include "model.h"
#include "project.h"
#include <iostream>

int main() {
    project prj;
    prj.name = "RectangularWG";
    prj.full_path_name = prj.data_path + prj.name;
    prj.task = project::LOAD_FES;
    prj.execute_task();
    // prj.model.msh.save_vtk_mesh(prj.name);
    std::cout << "Nodes  = " << prj.model.msh.n_nodes << "\n"
              << "Edges  = " << prj.model.msh.n_edges << "\n"
              << "Triangles  = " << prj.model.msh.n_faces << "\n"
              << "Tetrahedra = " << prj.model.msh.n_tetras << "\n";
    prj.task = project::ANALYZE;
    prj.execute_task();
    return 0;
}
