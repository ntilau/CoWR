#include "model.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

mdl_sld::mdl_sld() {}

mdl_sld::~mdl_sld() {}

double mdl_sld::get_geom_dim() {
    double dimx, dimy, dimz;
    std::vector<std::vector<double>> transpose(3, std::vector<double>(nodes.size()));
    for (unsigned int i = 0; i < nodes.size(); i++)
        for (unsigned int j = 0; j < 3; j++)
            transpose[j][i] = nodes[i][j];
    dimx = *std::max_element(transpose[0].begin(), transpose[0].end()) -
           *std::min_element(transpose[0].begin(), transpose[0].end());
    dimy = *std::max_element(transpose[1].begin(), transpose[1].end()) -
           *std::min_element(transpose[1].begin(), transpose[1].end());
    dimz = *std::max_element(transpose[2].begin(), transpose[2].end()) -
           *std::min_element(transpose[2].begin(), transpose[2].end());
    return dimx * dimy * dimz;
}

bool mdl_sld::read_stl_file(std::string &name) {
    clear();
    std::map<std::tuple<double, double, double>, size_t> nod_map;
    bool is_binary = false;
    std::string tmp_str, tmp_str_sld;
    std::string line;
    std::vector<double> tmp_pnt(3);
    std::vector<size_t> tmp_tri(3);
    size_t nod_cnt = 0, fac_cnt = 0;
    std::ifstream stl_file(std::string(name + ".stl").data());
    if (stl_file.is_open()) {
        getline(stl_file, line);
        std::istringstream iss;
        iss.str(line);
        iss >> tmp_str_sld;
        if (tmp_str_sld == "solid") {
            iss >> sld_name;
            std::cout << "Solid name: " << sld_name << "\n";
            dim = 3;
            do {
                getline(stl_file, line);
                iss.clear();
                iss.str(line);
                iss >> tmp_str; // facet
                if (tmp_str == "endsolid") {
                    getline(stl_file, line);
                    iss.clear();
                    iss.str(line);
                    iss >> tmp_str; // facet
                    if (tmp_str == "endsolid")
                        break;
                    else if (stl_file.eof())
                        break;
                }
                iss >> tmp_str;    // normal
                iss >> tmp_pnt[0]; // ni
                iss >> tmp_pnt[1]; // nj
                iss >> tmp_pnt[2]; // nk
                faces_normals.push_back(tmp_pnt);
                getline(stl_file, line);
                iss.clear();
                iss.str(line);
                iss >> tmp_str; // outer
                iss >> tmp_str; // loop
                getline(stl_file, line);
                iss.clear();
                iss.str(line);
                iss >> tmp_str;    // vertex
                iss >> tmp_pnt[0]; // v1x
                iss >> tmp_pnt[1]; // v1y
                iss >> tmp_pnt[2]; // v1z
                if (nod_map.find(std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])) ==
                    nod_map.end()) {
                    nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])] =
                        nod_cnt++;
                }
                tmp_tri[0] =
                    nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])];
                getline(stl_file, line);
                iss.clear();
                iss.str(line);
                iss >> tmp_str;    // vertex
                iss >> tmp_pnt[0]; // v2x
                iss >> tmp_pnt[1]; // v2y
                iss >> tmp_pnt[2]; // v2z
                if (nod_map.find(std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])) ==
                    nod_map.end()) {
                    nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])] =
                        nod_cnt++;
                }
                tmp_tri[1] =
                    nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])];
                getline(stl_file, line);
                iss.clear();
                iss.str(line);
                iss >> tmp_str;    // vertex
                iss >> tmp_pnt[0]; // v3x
                iss >> tmp_pnt[1]; // v3y
                iss >> tmp_pnt[2]; // v3z
                if (nod_map.find(std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])) ==
                    nod_map.end()) {
                    nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])] =
                        nod_cnt++;
                }
                tmp_tri[2] =
                    nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])];
                face_type fac_ply;
                fac_ply.polygons.push_back(tmp_tri);
                faces.push_back(fac_ply);
                faces_marker.push_back(std::vector<int>(1, 1));
                fac_cnt++;
                getline(stl_file, line); // endloop
                getline(stl_file, line); // endfacet
            } while (!stl_file.eof());
            stl_file.close();
            std::cout << "Nodes = " << nod_cnt << "\n";
            std::cout << "Faces = " << fac_cnt << "\n";
            nodes.resize(nod_cnt);
            for (std::map<std::tuple<double, double, double>, size_t>::iterator it =
                     nod_map.begin();
                 it != nod_map.end(); it++) {
                std::vector<double> node(3);
                node[0] = std::get<0>(it->first);
                node[1] = std::get<1>(it->first);
                node[2] = std::get<2>(it->first);
                nodes[it->second] = node;
            }
            max_faces_marker = 1;
            bc_markers.push_back(1);
        } else {
            is_binary = true;
        }
    } else {
        std::cout << name + ".stl not found\n";
        return false;
    }
    return true;
}

bool mdl_sld::read_poly_file(std::string &name) {
    clear();
    size_t tmp_int, n_pts;
    double tmp_dbl;
    std::string tmp_str;
    std::string line;
    std::ifstream poly_file(std::string(name + ".poly").data());
    if (poly_file.is_open()) {
        {
            do {
                getline(poly_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            std::istringstream iss(line);
            iss >> tmp_int;
            nodes.resize(tmp_int);
            iss >> n_pts;
            for (size_t i = 0; i < nodes.size(); i++) {
                do {
                    getline(poly_file, line);
                } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                std::istringstream iss(line);
                iss >> tmp_int;
                iss >> tmp_dbl;
                nodes[i].push_back(tmp_dbl);
                iss >> tmp_dbl;
                nodes[i].push_back(tmp_dbl);
                iss >> tmp_dbl;
                if (n_pts < 3)
                    nodes[i].push_back(0.0);
                else
                    nodes[i].push_back(tmp_dbl);
            }
        }
        dim = check_dim();
        std::cout << "Model is " << dim << "-dimensional\n";
        if (dim == 3) {
            do {
                getline(poly_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            std::istringstream iss(line);
            iss >> tmp_int;
            faces.resize(tmp_int);
            faces_marker.resize(tmp_int);
            unsigned int n_faces_marker;
            iss >> n_faces_marker;
            for (size_t i = 0; i < faces.size(); i++) {
                int polygons = 0, hole = 0, bmark = 0;
                do {
                    getline(poly_file, line);
                } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                std::istringstream iss(line);
                iss >> polygons;
                iss >> hole;
                iss >> bmark;
                max_faces_marker = std::max(max_faces_marker, bmark);
                if (std::find(bc_markers.begin(), bc_markers.end(), bmark) ==
                    bc_markers.end())
                    bc_markers.push_back(bmark);
                faces_marker[i].push_back(bmark);
                faces[i].polygons.resize(polygons);
                if (hole > 0)
                    faces[i].holes.resize(hole);
                for (size_t j = 0; j < polygons; j++) {
                    unsigned int nodes_nbr = 0;
                    do {
                        getline(poly_file, line);
                    } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                    std::istringstream iss(line);
                    iss >> nodes_nbr;
                    for (unsigned int k = 0; k < nodes_nbr; k++) {
                        iss >> tmp_int;
                        if (faces[i].polygons[j].size() > 0) {
                            if (faces[i].polygons[j].back() == (tmp_int - 1)) {
                                do {
                                    getline(poly_file, line);
                                } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                                iss.clear();
                                iss.str("");
                                iss.str(line);
                                iss >> tmp_int;
                            }
                        }
                        faces[i].polygons[j].push_back(tmp_int - 1);
                    }
                }
                for (size_t j = 0; j < hole; j++) {
                    unsigned int nodes_nbr = 0;
                    do {
                        getline(poly_file, line);
                    } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                    std::istringstream iss(line);
                    iss >> tmp_int;
                    iss >> tmp_dbl;
                    faces[i].holes[j].push_back(tmp_dbl);
                    iss >> tmp_dbl;
                    faces[i].holes[j].push_back(tmp_dbl);
                    iss >> tmp_dbl;
                    faces[i].holes[j].push_back(tmp_dbl);
                }
            }
        } else if (dim == 2) {
            do {
                getline(poly_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            std::istringstream iss(line);
            iss >> tmp_int;
            edges.resize(tmp_int);
            edges_marker.resize(tmp_int);
            iss >> tmp_int;
            size_t n_edges_marker = tmp_int;
            for (size_t i = 0; i < edges.size(); i++) {
                do {
                    getline(poly_file, line);
                } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                std::istringstream iss(line);
                iss >> tmp_int;
                iss >> tmp_int;
                edges[i].push_back(tmp_int - 1);
                iss >> tmp_int;
                edges[i].push_back(tmp_int - 1);
                for (unsigned int j = 0; j < n_edges_marker; j++) {
                    iss >> tmp_int;
                    int bmark = tmp_int;
                    max_edges_marker = std::max(max_edges_marker, bmark);
                    edges_marker[i].push_back(bmark);
                    if (std::find(bc_markers.begin(), bc_markers.end(), bmark) ==
                        bc_markers.end())
                        bc_markers.push_back(bmark);
                }
            }
        }
        // importing holes
        {
            do {
                getline(poly_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            std::istringstream iss(line);
            iss >> tmp_int;
            if (tmp_int > 0) {
                holes.resize(tmp_int);
                for (size_t i = 0; i < holes.size(); i++) {
                    do {
                        getline(poly_file, line);
                    } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                    std::istringstream iss(line);
                    iss >> tmp_int;
                    iss >> tmp_dbl;
                    holes[i].push_back(tmp_dbl);
                    iss >> tmp_dbl;
                    holes[i].push_back(tmp_dbl);
                    if (dim == 2)
                        holes[i].push_back(0.0);
                }
            }
        }
        // importing regions
        {
            do {
                getline(poly_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            std::istringstream iss(line);
            iss >> tmp_int;
            if (tmp_int > 0) {
                regions.resize(tmp_int);
                regions_marker.resize(tmp_int);
                for (size_t i = 0; i < regions.size(); i++) {
                    do {
                        getline(poly_file, line);
                    } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                    std::istringstream iss(line);
                    iss >> tmp_int;
                    iss >> tmp_dbl;
                    regions[i].push_back(tmp_dbl);
                    iss >> tmp_dbl;
                    regions[i].push_back(tmp_dbl);
                    if (dim == 3) {
                        iss >> tmp_dbl;
                        regions[i].push_back(tmp_dbl);
                    } else if (dim == 2)
                        regions[i].push_back(0.0);
                    iss >> regions_marker[i];
                    max_regions_marker = std::max(max_regions_marker, regions_marker[i]);
                    if (std::find(mtrl_markers.begin(), mtrl_markers.end(),
                                  regions_marker[i]) == mtrl_markers.end())
                        mtrl_markers.push_back(regions_marker[i]);
                }
            }
        }
        std::sort(bc_markers.begin(), bc_markers.end());
        std::sort(mtrl_markers.begin(), mtrl_markers.end());
        return true;
    } else {
        std::cout << name + ".poly not found\n";
        return false;
    }
}

void mdl_sld::clear() {
    std::vector<std::vector<double>>().swap(nodes);
    std::vector<std::vector<size_t>>().swap(edges);
    std::vector<std::vector<int>>().swap(edges_marker);
    std::vector<face_type>().swap(faces);
    std::vector<std::vector<int>>().swap(faces_marker);
    std::vector<std::vector<double>>().swap(holes);
    std::vector<std::vector<double>>().swap(regions);
    std::vector<int>().swap(regions_marker);
    std::vector<std::vector<double>>().swap(bounding_box);
    std::vector<int>().swap(bc_markers);
    std::vector<int>().swap(mtrl_markers);
    max_dimension = 0.0;
    dim = 0;
    max_edges_marker = -INT_MAX;
    max_faces_marker = -INT_MAX;
    max_regions_marker = -INT_MAX;
}

void mdl_sld::get_bounding_info() {
    bounding_box.resize(2);
    bounding_box[0].assign(3, DBL_MAX);
    bounding_box[1].assign(3, -DBL_MAX);
    for (size_t i = 0; i < nodes.size(); i++) {
        bounding_box[0][0] = std::min(bounding_box[0][0], nodes[i][0]);
        bounding_box[0][1] = std::min(bounding_box[0][1], nodes[i][1]);
        bounding_box[0][2] = std::min(bounding_box[0][2], nodes[i][2]);
        bounding_box[1][0] = std::max(bounding_box[1][0], nodes[i][0]);
        bounding_box[1][1] = std::max(bounding_box[1][1], nodes[i][1]);
        bounding_box[1][2] = std::max(bounding_box[1][2], nodes[i][2]);
    }
    double x_dim, y_dim, z_dim;
    x_dim = std::abs(bounding_box[1][0] - bounding_box[0][0]);
    y_dim = std::abs(bounding_box[1][1] - bounding_box[0][1]);
    z_dim = std::abs(bounding_box[1][2] - bounding_box[0][2]);
    max_dimension = std::max(max_dimension, x_dim);
    max_dimension = std::max(max_dimension, y_dim);
    max_dimension = std::max(max_dimension, z_dim);
}

unsigned int mdl_sld::check_dim() {
    double x, y, z;
    if (nodes.size() >= 1) {
        x = nodes[0][0];
        y = nodes[0][1];
        z = nodes[0][2];
    }
    double x_norm = 0, y_norm = 0, z_norm = 0;
    for (size_t i = 0; i < nodes.size(); i++) {
        x_norm += std::abs(nodes[i][0] * nodes[i][0]);
        y_norm += std::abs(nodes[i][1] * nodes[i][1]);
        z_norm += std::abs(nodes[i][2] * nodes[i][2]);
    }
    x_norm = std::sqrt(x_norm);
    y_norm = std::sqrt(y_norm);
    z_norm = std::sqrt(z_norm);
    if (x == x_norm || y == y_norm || z == z_norm)
        return 2;
    else
        return 3;
}

void mdl_sld::write_prj_file(std::string &name) {
    std::ofstream sld_out_file(std::string(name + ".fes").c_str(),
                               std::ios::out | std::ios::ate | std::ios::app);
    sld_out_file << "#Sld_Nodes " << nodes.size() << "\n";
    for (size_t i = 0; i < nodes.size(); i++) {
        sld_out_file << std::scientific << std::setprecision(16) << nodes[i][0]
                     << " " << std::setprecision(16) << nodes[i][1] << " "
                     << std::setprecision(16) << nodes[i][2] << "\n";
    }
    if (dim == 2) {
        sld_out_file << "#Sld_Edges " << edges.size() << " "
                     << edges_marker[0].size() << "\n";
        for (size_t i = 0; i < edges.size(); i++) {
            for (size_t j = 0; j < edges[i].size(); j++)
                sld_out_file << edges[i][j] << " ";
            for (size_t j = 0; j < edges_marker[i].size(); j++)
                sld_out_file << edges_marker[i][j] << " ";
            sld_out_file << "\n";
        }
    }
    if (dim == 3) {
        sld_out_file << "#Sld_Faces " << faces.size() << "\n";
        for (size_t i = 0; i < faces.size(); i++) {
            sld_out_file << faces[i].polygons.size() << " " << faces[i].holes.size()
                         << " " << faces_marker[i].back() << "\n";
            for (size_t j = 0; j < faces[i].polygons.size(); j++) {
                sld_out_file << faces[i].polygons[j].size() << " ";
                for (size_t k = 0; k < faces[i].polygons[j].size(); k++) {
                    sld_out_file << faces[i].polygons[j][k] << " ";
                }
                sld_out_file << "\n";
            }
            for (size_t j = 0; j < faces[i].holes.size(); j++) {
                sld_out_file << std::scientific << std::setprecision(16)
                             << faces[i].holes[j][0] << " " << std::setprecision(16)
                             << faces[i].holes[j][1] << " " << std::setprecision(16)
                             << faces[i].holes[j][2] << "\n";
            }
        }
    }
    sld_out_file << "#Sld_Holes " << holes.size() << "\n";
    for (size_t i = 0; i < holes.size(); i++) {
        sld_out_file << std::scientific << std::setprecision(16) << holes[i][0]
                     << " " << std::setprecision(16) << holes[i][1] << " "
                     << std::setprecision(16) << holes[i][2] << "\n";
    }
    sld_out_file << "#Sld_Regions " << regions.size() << "\n";
    for (size_t i = 0; i < regions.size(); i++) {
        sld_out_file << std::scientific << std::setprecision(16) << regions[i][0]
                     << " " << std::setprecision(16) << regions[i][1] << " "
                     << std::setprecision(16) << regions[i][2] << " "
                     << regions_marker[i] << "\n";
    }
    sld_out_file.close();
}

void mdl_sld::read_prj_file(std::string &name) {
    clear();
    std::ifstream sld_in_file(std::string(name + ".fes").c_str(), std::ios::in);
    std::string line;
    std::istringstream iss;
    unsigned int tmp_uint;
    double tmp_dbl;
    std::string tmp_str;
    if (sld_in_file.is_open()) {
        while (getline(sld_in_file, line)) {
            iss.clear();
            iss.str(line);
            iss >> tmp_str;
            if (strcmp(tmp_str.data(), "#Sld_Nodes") == 0) {
                iss >> tmp_uint;
                nodes.resize(tmp_uint);
                for (size_t i = 0; i < nodes.size(); i++) {
                    getline(sld_in_file, line);
                    iss.clear();
                    iss.str(line);
                    nodes[i].resize(3);
                    iss >> nodes[i][0];
                    iss >> nodes[i][1];
                    iss >> nodes[i][2];
                }
            }
            if (strcmp(tmp_str.data(), "#Sld_Edges") == 0) {
                dim = 2;
                iss >> tmp_uint;
                edges.resize(tmp_uint);
                edges_marker.resize(tmp_uint);
                iss >> tmp_uint;
                for (size_t i = 0; i < edges.size(); i++) {
                    edges[i].resize(2);
                    edges_marker[i].resize(tmp_uint);
                }
                for (size_t i = 0; i < edges.size(); i++) {
                    getline(sld_in_file, line);
                    iss.clear();
                    iss.str(line);
                    iss >> edges[i][0];
                    iss >> edges[i][1];
                    for (unsigned int j = 0; j < edges_marker[i].size(); j++) {
                        iss >> edges_marker[i][j];
                        max_edges_marker = std::max(max_edges_marker, edges_marker[i][j]);
                        if (std::find(bc_markers.begin(), bc_markers.end(),
                                      edges_marker[i][j]) == bc_markers.end())
                            bc_markers.push_back(edges_marker[i][j]);
                    }
                }
            }
            if (strcmp(tmp_str.data(), "#Sld_Faces") == 0) {
                dim = 3;
                iss >> tmp_uint;
                faces.resize(tmp_uint);
                faces_marker.resize(tmp_uint);
                for (size_t i = 0; i < faces.size(); i++) {
                    getline(sld_in_file, line);
                    iss.clear();
                    iss.str(line);
                    unsigned int polygons, hole;
                    int bmark;
                    iss >> polygons;
                    iss >> hole;
                    iss >> bmark;
                    faces[i].polygons.resize(polygons);
                    faces[i].holes.resize(hole);
                    faces_marker[i].push_back(bmark);
                    max_faces_marker = std::max(max_faces_marker, bmark);
                    if (std::find(bc_markers.begin(), bc_markers.end(), bmark) ==
                        bc_markers.end())
                        bc_markers.push_back(bmark);
                    for (unsigned int j = 0; j < faces[i].polygons.size(); j++) {
                        getline(sld_in_file, line);
                        iss.clear();
                        iss.str(line);
                        iss >> tmp_uint;
                        faces[i].polygons[j].resize(tmp_uint);
                        for (unsigned int k = 0; k < faces[i].polygons[j].size(); k++)
                            iss >> faces[i].polygons[j][k];
                    }
                    for (unsigned int j = 0; j < faces[i].holes.size(); j++) {
                        getline(sld_in_file, line);
                        iss.clear();
                        iss.str(line);
                        faces[i].holes[j].resize(3);
                        for (unsigned int k = 0; k < 3; k++)
                            iss >> faces[i].holes[j][k];
                    }
                }
            }
            if (strcmp(tmp_str.data(), "#Sld_Holes") == 0) {
                iss >> tmp_uint;
                holes.resize(tmp_uint);
                for (size_t i = 0; i < holes.size(); i++) {
                    getline(sld_in_file, line);
                    iss.clear();
                    iss.str(line);
                    nodes[i].resize(3);
                    iss >> nodes[i][0];
                    iss >> nodes[i][1];
                    iss >> nodes[i][2];
                }
            }
            if (strcmp(tmp_str.data(), "#Sld_Regions") == 0) {
                iss >> tmp_uint;
                regions.resize(tmp_uint);
                regions_marker.resize(tmp_uint);
                for (size_t i = 0; i < regions.size(); i++) {
                    getline(sld_in_file, line);
                    iss.clear();
                    iss.str(line);
                    regions[i].resize(3);
                    iss >> regions[i][0];
                    iss >> regions[i][1];
                    iss >> regions[i][2];
                    iss >> regions_marker[i];
                    max_regions_marker = std::max(max_regions_marker, regions_marker[i]);
                    if (std::find(mtrl_markers.begin(), mtrl_markers.end(),
                                  regions_marker[i]) == mtrl_markers.end())
                        mtrl_markers.push_back(regions_marker[i]);
                }
            }
        }
    }
    sld_in_file.close();
}

mdl_bc::mdl_bc() {
    label = 0;
}

mdl_bc::~mdl_bc() {
    std::vector<std::complex<double>>().swap(mode_beta);
    std::vector<std::vector<std::complex<double>>>().swap(mode_eig_vec);
    std::vector<std::vector<std::complex<double>>>().swap(mode_eig_vec_f);
    std::vector<std::vector<size_t>>().swap(mode_dof_map);
}

mdl_mtrl::mdl_mtrl() {}

mdl_mtrl::~mdl_mtrl() {}

void mdl_mtrl::upd_mtrl() {
    epsr2 = -tand * epsr;
}

void mdl_mtrl::upd_mtrl(double &freq) {
    epsr2 =
        -sigma / (2.0 * phys_const::pi * freq * phys_const::eps0) - tand * epsr;
}

double mdl_mtrl::calc_epsr2(double &freq) {
    return (-sigma / (2.0 * phys_const::pi * freq * phys_const::eps0) -
            tand * epsr);
}

mdl_msh::mdl_msh() {
    n_domains = 0;
    n_edges = 0;
    n_faces = 0;
    n_nodes = 0;
    n_tetras = 0;
}

void mdl_msh::write_prj_file(std::string &name) {
    std::ofstream msh_out_file(std::string(name + ".fes").c_str(),
                               std::ios::out | std::ios::ate | std::ios::app);
    msh_out_file << "#Mesh " << type << "\n";
    msh_out_file << "#Nodes " << n_nodes << "\n";
    for (size_t i = 0; i < n_nodes; i++) {
        msh_out_file << std::scientific << std::setprecision(16) << nod_pos[i][0]
                     << " " << std::setprecision(16) << nod_pos[i][1] << " "
                     << std::setprecision(16) << nod_pos[i][2] << "\n";
    }
    msh_out_file << "#Edges " << n_edges << "\n";
    for (size_t i = 0; i < n_edges; i++) {
        msh_out_file << edg_nodes[i][0] << " " << edg_nodes[i][1] << " "
                     << edg_lab[i] << "\n";
    }
    msh_out_file << "#Faces " << n_faces << "\n";
    for (size_t i = 0; i < n_faces; i++) {
        msh_out_file << fac_nodes[i][0] << " " << fac_nodes[i][1] << " "
                     << fac_nodes[i][2] << " " << fac_lab[i] << "\n";
    }
    msh_out_file << "#Tetras " << n_tetras << "\n";
    for (size_t i = 0; i < n_tetras; i++) {
        msh_out_file << tet_nodes[i][0] << " " << tet_nodes[i][1] << " "
                     << tet_nodes[i][2] << " " << tet_nodes[i][3] << " "
                     << tet_lab[i] << "\n";
    }
    msh_out_file.close();
}

void mdl_msh::read_prj_file(std::string &name) {
    clear();
    std::ifstream msh_in_file(std::string(name + ".fes").c_str(), std::ios::in);
    std::string line;
    std::istringstream iss;
    unsigned int tmp_uint;
    double tmp_dbl;
    std::string tmp_str;
    if (msh_in_file.is_open()) {
        while (getline(msh_in_file, line)) {
            iss.clear();
            iss.str(line);
            iss >> tmp_str;
            if (strcmp(tmp_str.data(), "#Mesh") == 0) {
                iss >> type;
            }
            if (strcmp(tmp_str.data(), "#Nodes") == 0) {
                iss >> n_nodes;
                nod_pos.resize(n_nodes);
                for (size_t i = 0; i < n_nodes; i++) {
                    getline(msh_in_file, line);
                    iss.clear();
                    iss.str(line);
                    nod_pos[i].resize(3);
                    iss >> nod_pos[i][0];
                    iss >> nod_pos[i][1];
                    iss >> nod_pos[i][2];
                }
            }
            if (strcmp(tmp_str.data(), "#Edges") == 0) {
                iss >> n_edges;
                edg_nodes.resize(n_edges);
                edg_lab.resize(n_edges);
                for (size_t i = 0; i < n_edges; i++) {
                    getline(msh_in_file, line);
                    iss.clear();
                    iss.str(line);
                    edg_nodes[i].resize(2);
                    iss >> edg_nodes[i][0];
                    iss >> edg_nodes[i][1];
                    iss >> edg_lab[i];
                }
            }
            if (strcmp(tmp_str.data(), "#Faces") == 0) {
                iss >> n_faces;
                fac_nodes.resize(n_faces);
                fac_lab.resize(n_faces);
                for (size_t i = 0; i < n_faces; i++) {
                    getline(msh_in_file, line);
                    iss.clear();
                    iss.str(line);
                    fac_nodes[i].resize(3);
                    iss >> fac_nodes[i][0];
                    iss >> fac_nodes[i][1];
                    iss >> fac_nodes[i][2];
                    iss >> fac_lab[i];
                }
            }
            if (strcmp(tmp_str.data(), "#Tetras") == 0) {
                iss >> n_tetras;
                tet_nodes.resize(n_tetras);
                tet_lab.resize(n_tetras);
                for (size_t i = 0; i < n_tetras; i++) {
                    getline(msh_in_file, line);
                    iss.clear();
                    iss.str(line);
                    tet_nodes[i].resize(4);
                    iss >> tet_nodes[i][0];
                    iss >> tet_nodes[i][1];
                    iss >> tet_nodes[i][2];
                    iss >> tet_nodes[i][3];
                    iss >> tet_lab[i];
                }
            }
        }
    }
    msh_in_file.close();
    regularize_mesh();
}

void mdl_msh::read_tetgen_files(std::string &name) {
    clear();
    type = "TETRA";
    unsigned int lvl = 1;
    std::ostringstream str_lvl;
    str_lvl << lvl;
    std::string line;
    std::istringstream iss;
    double tmp_dbl;
    int tmp_int;
    size_t tmp_uint;
    std::string tmp_str;
    std::ifstream tet_node_file(
        std::string(name + "." + str_lvl.str() + ".node").c_str());
    if (tet_node_file.is_open()) {
        std::cout << "Loading " << std::string(name + "." + str_lvl.str() + ".node")
                  << "\n";
        iss.clear();
        do {
            getline(tet_node_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_nodes;
        iss >> tmp_int;
        nod_pos.resize(n_nodes);
        for (size_t i = 0; i < n_nodes; i++) {
            iss.clear();
            do {
                getline(tet_node_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 3; j++) {
                iss >> tmp_dbl;
                nod_pos[i].push_back(tmp_dbl);
            }
        }
    }
    tet_node_file.close();
    std::ifstream tet_edge_file(
        std::string(name + "." + str_lvl.str() + ".edge").c_str());
    if (tet_edge_file.is_open()) {
        std::cout << "Loading " << std::string(name + "." + str_lvl.str() + ".edge")
                  << "\n";
        unsigned int dim;
        iss.clear();
        do {
            getline(tet_edge_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_edges;
        iss >> tmp_int;
        edg_nodes.resize(n_edges);
        edg_lab.assign(n_edges, 0);
        for (size_t i = 0; i < n_edges; i++) {
            iss.clear();
            do {
                getline(tet_edge_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 2; j++) {
                iss >> tmp_uint;
                edg_nodes[i].push_back(tmp_uint - 1);
            }
            iss >> tmp_int;
            edg_lab[i] = tmp_int;
        }
    }
    tet_edge_file.close();
    std::ifstream tet_face_file(
        std::string(name + "." + str_lvl.str() + ".face").c_str());
    if (tet_face_file.is_open()) {
        std::cout << "Loading " << std::string(name + "." + str_lvl.str() + ".face")
                  << "\n";
        unsigned int dim;
        iss.clear();
        do {
            getline(tet_face_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_faces;
        iss >> tmp_int;
        fac_nodes.resize(n_faces);
        fac_lab.assign(n_faces, 0);
        for (size_t i = 0; i < n_faces; i++) {
            iss.clear();
            do {
                getline(tet_face_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 3; j++) {
                iss >> tmp_uint;
                fac_nodes[i].push_back(tmp_uint - 1);
            }
            iss >> fac_lab[i];
        }
    }
    tet_face_file.close();
    std::ifstream tet_ele_file(
        std::string(name + "." + str_lvl.str() + ".ele").c_str());
    if (tet_ele_file.is_open()) {
        std::cout << "Loading " << std::string(name + "." + str_lvl.str() + ".ele")
                  << "\n";
        unsigned int dim;
        iss.clear();
        do {
            getline(tet_ele_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_tetras;
        iss >> tmp_int;
        tet_nodes.resize(n_tetras);
        tet_lab.assign(n_tetras, 0);
        for (size_t i = 0; i < n_tetras; i++) {
            iss.clear();
            do {
                getline(tet_ele_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 4; j++) {
                iss >> tmp_uint;
                tet_nodes[i].push_back(tmp_uint - 1);
            }
            iss >> tmp_int;
            tet_lab[i] = tmp_int;
        }
    }
    tet_ele_file.close();
    regularize_mesh();
}

void mdl_msh::read_triangle_files(std::string &name) {
    clear();
    type = "TRIA";
    unsigned int lvl = 1;
    std::ostringstream str_lvl;
    str_lvl << lvl;
    std::string line;
    std::istringstream iss;
    double tmp_dbl;
    int tmp_int;
    size_t tmp_uint;
    std::string tmp_str;
    std::ifstream tria_node_file(
        std::string(name + "." + str_lvl.str() + ".node").c_str());
    if (tria_node_file.is_open()) {
        std::cout << "Loading " << std::string(name + "." + str_lvl.str() + ".node")
                  << "\n";
        iss.clear();
        do {
            getline(tria_node_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_nodes;
        iss >> tmp_int;
        nod_pos.resize(n_nodes);
        for (size_t i = 0; i < n_nodes; i++) {
            iss.clear();
            do {
                getline(tria_node_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 3; j++) {
                iss >> tmp_dbl;
                nod_pos[i].push_back(tmp_dbl);
            }
        }
    }
    tria_node_file.close();
    std::ifstream tria_edge_file(
        std::string(name + "." + str_lvl.str() + ".edge").c_str());
    if (tria_edge_file.is_open()) {
        std::cout << "Loading " << std::string(name + "." + str_lvl.str() + ".edge")
                  << "\n";
        unsigned int dim;
        iss.clear();
        do {
            getline(tria_edge_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_edges;
        iss >> tmp_int;
        edg_nodes.resize(n_edges);
        edg_lab.assign(n_edges, 0);
        for (size_t i = 0; i < n_edges; i++) {
            iss.clear();
            do {
                getline(tria_edge_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 2; j++) {
                iss >> tmp_uint;
                edg_nodes[i].push_back(tmp_uint - 1);
            }
            iss >> tmp_int;
            edg_lab[i] = tmp_int;
        }
    }
    tria_edge_file.close();
    std::ifstream tria_ele_file(
        std::string(name + "." + str_lvl.str() + ".ele").c_str());
    if (tria_ele_file.is_open()) {
        std::cout << "Loading " << std::string(name + "." + str_lvl.str() + ".ele")
                  << "\n";
        unsigned int dim;
        iss.clear();
        do {
            getline(tria_ele_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_faces;
        iss >> tmp_int;
        fac_nodes.resize(n_faces);
        fac_lab.assign(n_faces, 0);
        for (size_t i = 0; i < n_faces; i++) {
            iss.clear();
            do {
                getline(tria_ele_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 3; j++) {
                iss >> tmp_uint;
                fac_nodes[i].push_back(tmp_uint - 1);
            }
            iss >> tmp_int;
            fac_lab[i] = tmp_int;
        }
    }
    tria_ele_file.close();
    regularize_mesh();
}

/// Refinement mapping
static unsigned int tet_tet_nodes_loc_map_id[] = {
    0, 4, 5, 6, 1, 4, 7, 8, 2, 5, 7, 9, 3, 6, 8, 9,
    4, 5, 6, 7, 4, 6, 7, 8, 5, 6, 7, 9, 6, 7, 8, 9
};

static unsigned int tet_fac_nodes_loc_map_id[] = {
    1, 7, 8, 2, 7, 9, 3, 8, 9, 7, 8, 9, 0, 5, 6, 2, 5, 9, 3, 6, 9, 5, 6, 9,
    0, 4, 6, 1, 4, 8, 3, 6, 8, 4, 6, 8, 0, 4, 5, 1, 4, 7, 2, 5, 7, 4, 5, 7,
    4, 5, 6, 4, 6, 7, 4, 7, 8, 5, 6, 7, 5, 7, 9, 6, 7, 8, 6, 7, 9, 6, 8, 9
};

static unsigned int fac_fac_nodes_loc_map_id[] = {0, 4, 5, 1, 3, 5,
                                                  2, 3, 4, 3, 4, 5
                                                 };

static unsigned int tet_edg_nodes_loc_map_id[] = {
    0, 4, 0, 5, 0, 6, 1, 4, 1, 7, 1, 8, 2, 5, 2, 7, 2, 9, 3, 6, 3, 8, 3, 9, 4,
    5, 4, 6, 4, 7, 4, 8, 5, 6, 5, 7, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9
};

static unsigned int fac_edg_nodes_loc_map_id[] = {1, 3, 2, 3, 0, 4, 2, 4, 0,
                                                  5, 1, 5, 3, 4, 3, 5, 4, 5
                                                 };

static int tet_fac_lab_loc_map_id[] = {0,  0,  0,  0,  1,  1,  1,  1,
                                       2,  2,  2,  2,  3,  3,  3,  3,
                                       -1, -1, -1, -1, -1, -1, -1, -1
                                       };

static int fac_edg_lab_loc_map_id[] = {0, 0, 1, 1, 2, 2, -1, -1, -1};

mdl_msh::mdl_msh(mdl_msh *msh) {
    msh->tet_nodes = tet_nodes;
    msh->tet_edges = tet_edges;
    msh->tet_faces = tet_faces;
    msh->fac_nodes = fac_nodes;
    msh->fac_edges = fac_edges;
    msh->edg_nodes = edg_nodes;
    msh->nod_pos = nod_pos;
    msh->edg_lab = edg_lab;
    msh->fac_lab = fac_lab;
    msh->tet_lab = tet_lab;
    msh->fac_adj_tet = fac_adj_tet;
    msh->edg_adj_fac = edg_adj_fac;
    msh->dom_tetras = dom_tetras;
    msh->dom_faces = dom_faces;
    msh->n_nodes = n_nodes;
    msh->n_edges = n_edges;
    msh->n_faces = n_faces;
    msh->n_tetras = n_tetras;
    msh->n_domains = n_domains;
}

mdl_msh::~mdl_msh() {
    clear();
}

void mdl_msh::get_mesh_statistics() {
    std::cout << "tet_nodes = " << tet_nodes.size() << "\n";
    std::cout << "tet_edges = " << tet_edges.size() << "\n";
    std::cout << "tet_faces = " << tet_faces.size() << "\n";
    std::cout << "fac_nodes = " << fac_nodes.size() << "\n";
    std::cout << "fac_edges = " << fac_edges.size() << "\n";
    std::cout << "edg_nodes = " << edg_nodes.size() << "\n";
    std::cout << "nod_pos = " << nod_pos.size() << "\n";
    std::cout << "fac_lab = " << fac_lab.size() << "\n";
    std::cout << "tet_lab = " << tet_lab.size() << "\n";
    std::cout << "fac_adj_tet = " << fac_adj_tet.size() << "\n";
    std::cout << "edg_adj_fac = " << edg_adj_fac.size() << "\n";
    std::cout << "dom_tetras = " << dom_tetras.size() << "\n";
    std::cout << "dom_faces = " << dom_faces.size() << "\n";
    std::cout << "n_nodes = " << n_nodes << "\n";
    std::cout << "n_edges = " << n_edges << "\n";
    std::cout << "n_faces = " << n_faces << "\n";
    std::cout << "n_tetras = " << n_tetras << "\n";
    std::cout << "n_domains = " << n_domains << "\n";
}

void mdl_msh::clear() {
    tet_nodes.clear();
    tet_edges.clear();
    tet_faces.clear();
    fac_nodes.clear();
    fac_edges.clear();
    edg_nodes.clear();
    nod_pos.clear();
    fac_lab.clear();
    tet_lab.clear();
    fac_adj_tet.clear();
    edg_adj_fac.clear();
    dom_tetras.clear();
    dom_faces.clear();
    n_nodes = 0;
    n_edges = 0;
    n_faces = 0;
    n_tetras = 0;
    n_domains = 0;
    max_edg_marker = -INT_MAX;
    max_fac_marker = -INT_MAX;
    max_tet_marker = -INT_MAX;
}

void mdl_msh::regularize_mesh() {
    max_edg_marker = -INT_MAX;
    max_fac_marker = -INT_MAX;
    max_tet_marker = -INT_MAX;
    std::map<std::pair<size_t, size_t>, size_t> edgesMap;
    std::map<std::tuple<size_t, size_t, size_t>, size_t> facesMap;
    if (strcmp(type.data(), "EDGE")) { // if not edges
        for (size_t i = 0; i < n_edges; i++) {
            std::sort(edg_nodes[i].begin(), edg_nodes[i].end());
            edgesMap[std::make_pair(edg_nodes[i][0], edg_nodes[i][1])] = i;
        }
        fac_edges.clear();
        edg_adj_fac.clear();
        fac_edges.resize(n_faces);
        edg_adj_fac.resize(n_edges);
        for (size_t i = 0; i < n_faces; i++) {
            std::sort(fac_nodes[i].begin(), fac_nodes[i].end());
            fac_edges[i].push_back(
                edgesMap[std::make_pair(fac_nodes[i][1], fac_nodes[i][2])]);
            fac_edges[i].push_back(
                edgesMap[std::make_pair(fac_nodes[i][0], fac_nodes[i][2])]);
            fac_edges[i].push_back(
                edgesMap[std::make_pair(fac_nodes[i][0], fac_nodes[i][1])]);
            facesMap[std::make_tuple(fac_nodes[i][0], fac_nodes[i][1],
                                                                      fac_nodes[i][2])] = i;
            for (size_t j = 0; j < 3; j++) {
                edg_adj_fac[fac_edges[i][j]].push_back(i);
            }
        }
    }
    if (strcmp(type.data(), "TETRA") == 0) {
        tet_edges.clear();
        tet_faces.clear();
        fac_adj_tet.clear();
        tet_edges.resize(n_tetras);
        tet_faces.resize(n_tetras);
        fac_adj_tet.resize(n_faces);
        for (size_t i = 0; i < n_tetras; i++) {
            std::sort(tet_nodes[i].begin(), tet_nodes[i].end());
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][0], tet_nodes[i][1])]);
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][0], tet_nodes[i][2])]);
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][0], tet_nodes[i][3])]);
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][1], tet_nodes[i][2])]);
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][1], tet_nodes[i][3])]);
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][2], tet_nodes[i][3])]);
            tet_faces[i].push_back(facesMap[std::make_tuple(
                                                tet_nodes[i][1], tet_nodes[i][2], tet_nodes[i][3])]);
            tet_faces[i].push_back(facesMap[std::make_tuple(
                                                tet_nodes[i][0], tet_nodes[i][2], tet_nodes[i][3])]);
            tet_faces[i].push_back(facesMap[std::make_tuple(
                                                tet_nodes[i][0], tet_nodes[i][1], tet_nodes[i][3])]);
            tet_faces[i].push_back(facesMap[std::make_tuple(
                                                tet_nodes[i][0], tet_nodes[i][1], tet_nodes[i][2])]);
            for (size_t j = 0; j < 4; j++) {
                fac_adj_tet[tet_faces[i][j]].push_back(i);
            }
        }
    }
    for (size_t i = 0; i < edg_lab.size(); i++)
        max_edg_marker = std::max(max_edg_marker, edg_lab[i]);
    for (size_t i = 0; i < fac_lab.size(); i++)
        max_fac_marker = std::max(max_fac_marker, fac_lab[i]);
    for (size_t i = 0; i < tet_lab.size(); i++)
        max_tet_marker = std::max(max_tet_marker, tet_lab[i]);
}

void mdl_msh::refine_homogeneous() {
    std::vector<std::vector<double>> new_nod_pos(n_nodes + n_edges,
                                  std::vector<double>(3));
    std::vector<std::vector<size_t>> new_tet_nodes(n_tetras * 8,
                                  std::vector<size_t>(4));
    std::vector<std::vector<size_t>> new_tet_edges(n_tetras * 8,
                                  std::vector<size_t>(6));
    std::vector<std::vector<size_t>> new_tet_faces(n_tetras * 8,
                                  std::vector<size_t>(4));
    std::vector<std::vector<size_t>> new_fac_nodes(n_faces * 4 + n_tetras * 8,
                                  std::vector<size_t>(3));
    std::vector<std::vector<size_t>> new_fac_edges(n_faces * 4 + n_tetras * 8,
                                  std::vector<size_t>(3));
    std::vector<std::vector<size_t>> new_edg_nodes(
                                      n_edges * 2 + n_faces * 3 + n_tetras, std::vector<size_t>(2));
    std::vector<int> new_edg_lab(n_edges * 2 + n_faces * 3 + n_tetras, 0);
    std::vector<int> new_fac_lab(n_faces * 4 + n_tetras * 8, 0);
    std::vector<int> new_tet_lab(n_tetras * 8, 0);
    std::vector<std::vector<size_t>> new_edg_adj_fac(n_edges * 2 + n_faces * 3 +
                                  n_tetras);
    std::vector<std::vector<size_t>> new_fac_adj_tet(n_faces * 4 + n_tetras * 8);
    //
    std::map<std::pair<size_t, size_t>, size_t> edgesMap;
    std::map<std::tuple<size_t, size_t, size_t>, size_t> facesMap;
    //
    for (size_t nid = 0; nid < n_nodes; nid++)
        new_nod_pos[nid] = nod_pos[nid];
    for (size_t eid = 0; eid < n_edges; eid++) {
        new_nod_pos[n_nodes + eid][0] =
            (nod_pos[edg_nodes[eid][0]][0] + nod_pos[edg_nodes[eid][1]][0]) / 2;
        new_nod_pos[n_nodes + eid][1] =
            (nod_pos[edg_nodes[eid][0]][1] + nod_pos[edg_nodes[eid][1]][1]) / 2;
        new_nod_pos[n_nodes + eid][2] =
            (nod_pos[edg_nodes[eid][0]][2] + nod_pos[edg_nodes[eid][1]][2]) / 2;
    }
    size_t tet_lvl = 0;
    size_t fac_lvl = 0;
    size_t edg_lvl = 0;
    if (strcmp(type.data(), "TRIA") == 0) {
        std::vector<size_t> nod_glob(6), edg_tmp(2), fac_tmp(3);
        for (size_t fid = 0; fid < n_faces; fid++) {
            for (unsigned int i = 0; i < 3; i++) {
                nod_glob[i] = fac_nodes[fid][i];
                nod_glob[3 + i] = n_nodes + fac_edges[fid][i];
            }
            // populating edges
            for (unsigned int i = 0; i < 9; i++) {
                edg_tmp[0] = nod_glob[fac_edg_nodes_loc_map_id[2 * i]];
                edg_tmp[1] = nod_glob[fac_edg_nodes_loc_map_id[2 * i + 1]];
                std::sort(edg_tmp.begin(), edg_tmp.end());
                if (edgesMap.find(std::make_pair(edg_tmp[0], edg_tmp[1])) ==
                    edgesMap.end()) {
                    edgesMap[std::make_pair(edg_tmp[0], edg_tmp[1])] = edg_lvl;
                    new_edg_nodes[edg_lvl][0] = edg_tmp[0];
                    new_edg_nodes[edg_lvl][1] = edg_tmp[1];
                    if (fac_edg_lab_loc_map_id[i] > -1) {
                        new_edg_lab[edg_lvl] =
                            edg_lab[fac_edges[fid][fac_edg_lab_loc_map_id[i]]];
                    }
                    ++edg_lvl;
                }
            }
            // populating faces
            for (unsigned int i = 0; i < 4; i++) {
                fac_tmp[0] = nod_glob[fac_fac_nodes_loc_map_id[3 * i]];
                fac_tmp[1] = nod_glob[fac_fac_nodes_loc_map_id[3 * i + 1]];
                fac_tmp[2] = nod_glob[fac_fac_nodes_loc_map_id[3 * i + 2]];
                std::sort(fac_tmp.begin(), fac_tmp.end());
                new_fac_nodes[fac_lvl][0] = fac_tmp[0];
                new_fac_nodes[fac_lvl][1] = fac_tmp[1];
                new_fac_nodes[fac_lvl][2] = fac_tmp[2];
                new_fac_edges[fac_lvl][0] =
                    edgesMap[std::make_pair(fac_tmp[1], fac_tmp[2])];
                new_fac_edges[fac_lvl][1] =
                    edgesMap[std::make_pair(fac_tmp[0], fac_tmp[2])];
                new_fac_edges[fac_lvl][2] =
                    edgesMap[std::make_pair(fac_tmp[0], fac_tmp[1])];
                new_fac_lab[fac_lvl] = fac_lab[fid];
                for (size_t j = 0; j < 3; j++) {
                    size_t eid = new_fac_edges[fac_lvl][j];
                    new_edg_adj_fac[eid].push_back(fac_lvl);
                }
                ++fac_lvl;
            }
        }
    }
    if (strcmp(type.data(), "TETRA") == 0) {
        std::vector<size_t> nod_glob(10), edg_tmp(2), fac_tmp(3), tet_tmp(4);
        for (size_t tid = 0; tid < n_tetras; tid++) {
            for (unsigned int i = 0; i < 4; i++)
                nod_glob[i] = tet_nodes[tid][i];
            for (unsigned int i = 0; i < 6; i++)
                nod_glob[4 + i] = n_nodes + tet_edges[tid][i];
            // populating edges
            for (unsigned int i = 0; i < 25; i++) {
                edg_tmp[0] = nod_glob[tet_edg_nodes_loc_map_id[2 * i]];
                edg_tmp[1] = nod_glob[tet_edg_nodes_loc_map_id[2 * i + 1]];
                std::sort(edg_tmp.begin(), edg_tmp.end());
                if (edgesMap.find(std::make_pair(edg_tmp[0], edg_tmp[1])) ==
                    edgesMap.end()) {
                    edgesMap[std::make_pair(edg_tmp[0], edg_tmp[1])] = edg_lvl;
                    new_edg_nodes[edg_lvl][0] = edg_tmp[0];
                    new_edg_nodes[edg_lvl][1] = edg_tmp[1];
                    ++edg_lvl;
                }
            }
            // populating faces
            for (unsigned int i = 0; i < 24; i++) {
                fac_tmp[0] = nod_glob[tet_fac_nodes_loc_map_id[3 * i]];
                fac_tmp[1] = nod_glob[tet_fac_nodes_loc_map_id[3 * i + 1]];
                fac_tmp[2] = nod_glob[tet_fac_nodes_loc_map_id[3 * i + 2]];
                std::sort(fac_tmp.begin(), fac_tmp.end());
                if (facesMap.find(std::make_tuple(fac_tmp[0], fac_tmp[1],
                                                  fac_tmp[2])) == facesMap.end()) {
                    facesMap[std::make_tuple(fac_tmp[0], fac_tmp[1], fac_tmp[2])] =
                        fac_lvl;
                    new_fac_nodes[fac_lvl][0] = fac_tmp[0];
                    new_fac_nodes[fac_lvl][1] = fac_tmp[1];
                    new_fac_nodes[fac_lvl][2] = fac_tmp[2];
                    new_fac_edges[fac_lvl][0] =
                        edgesMap[std::make_pair(fac_tmp[1], fac_tmp[2])];
                    new_fac_edges[fac_lvl][1] =
                        edgesMap[std::make_pair(fac_tmp[0], fac_tmp[2])];
                    new_fac_edges[fac_lvl][2] =
                        edgesMap[std::make_pair(fac_tmp[0], fac_tmp[1])];
                    if (tet_fac_lab_loc_map_id[i] > -1) {
                        new_fac_lab[fac_lvl] =
                            fac_lab[tet_faces[tid][tet_fac_lab_loc_map_id[i]]];
                    }
                    ++fac_lvl;
                }
            }
            for (unsigned int i = 0; i < 8; i++) {
                for (unsigned int j = 0; j < 4; j++) {
                    tet_tmp[j] = nod_glob[tet_tet_nodes_loc_map_id[i * 4 + j]];
                }
                std::sort(tet_tmp.begin(), tet_tmp.end());
                for (unsigned int j = 0; j < 4; j++) {
                    new_tet_nodes[tet_lvl][j] = tet_tmp[j];
                }
                new_tet_lab[tet_lvl] = tet_lab[tid];
                new_tet_edges[tet_lvl][0] =
                    edgesMap[std::make_pair(tet_tmp[0], tet_tmp[1])];
                new_tet_edges[tet_lvl][1] =
                    edgesMap[std::make_pair(tet_tmp[0], tet_tmp[2])];
                new_tet_edges[tet_lvl][2] =
                    edgesMap[std::make_pair(tet_tmp[0], tet_tmp[3])];
                new_tet_edges[tet_lvl][3] =
                    edgesMap[std::make_pair(tet_tmp[1], tet_tmp[2])];
                new_tet_edges[tet_lvl][4] =
                    edgesMap[std::make_pair(tet_tmp[1], tet_tmp[3])];
                new_tet_edges[tet_lvl][5] =
                    edgesMap[std::make_pair(tet_tmp[2], tet_tmp[3])];
                new_tet_faces[tet_lvl][0] =
                    facesMap[std::make_tuple(tet_tmp[1], tet_tmp[2], tet_tmp[3])];
                new_tet_faces[tet_lvl][1] =
                    facesMap[std::make_tuple(tet_tmp[0], tet_tmp[2], tet_tmp[3])];
                new_tet_faces[tet_lvl][2] =
                    facesMap[std::make_tuple(tet_tmp[0], tet_tmp[1], tet_tmp[3])];
                new_tet_faces[tet_lvl][3] =
                    facesMap[std::make_tuple(tet_tmp[0], tet_tmp[1], tet_tmp[2])];
                for (size_t j = 0; j < 4; j++) {
                    size_t fid = new_tet_faces[tet_lvl][j];
                    new_fac_adj_tet[fid].push_back(tet_lvl);
                }
                ++tet_lvl;
            }
        }
    }
    nod_pos = new_nod_pos;
    edg_nodes = new_edg_nodes;
    edg_lab = new_edg_lab;
    edg_adj_fac = new_edg_adj_fac;
    fac_nodes = new_fac_nodes;
    fac_edges = new_fac_edges;
    fac_lab = new_fac_lab;
    fac_adj_tet = new_fac_adj_tet;
    tet_nodes = new_tet_nodes;
    tet_edges = new_tet_edges;
    tet_faces = new_tet_faces;
    tet_lab = new_tet_lab;
    n_nodes = n_nodes + n_edges;
    n_edges = n_edges * 2 + n_faces * 3 + n_tetras;
    n_faces = n_faces * 4 + n_tetras * 8;
    n_tetras = n_tetras * 8;
}

void mdl_msh::save_vtk_mesh(std::string vtkMshName) {
    unsigned int n_bc_mtrl_faces = 0;
    for (size_t i = 0; i < n_faces; i++) {
        if (fac_lab[i] > -1)
            ++n_bc_mtrl_faces;
        else if (fac_adj_tet[i].size() > 1) {
            if (tet_lab[fac_adj_tet[i][0]] != tet_lab[fac_adj_tet[i][1]]) {
                --fac_lab[i];
                ++n_bc_mtrl_faces;
            }
        }
    }
    std::ofstream out_vol_msh(std::string(vtkMshName + "_volmsh.vtk").data());
    out_vol_msh << "# vtk DataFile Version 2.0\n";
    out_vol_msh << "Mesh data\n";
    out_vol_msh << "ASCII\n";
    out_vol_msh << "DATASET UNSTRUCTURED_GRID\n";
    out_vol_msh << "POINTS " << n_nodes << " double \n";
    for (size_t i = 0; i < n_nodes; i++) {
        out_vol_msh << std::setprecision(16) << nod_pos[i][0] << " ";
        out_vol_msh << std::setprecision(16) << nod_pos[i][1] << " ";
        out_vol_msh << std::setprecision(16) << nod_pos[i][2] << "\n";
    }
    out_vol_msh << "CELLS " << n_tetras << " " << 5 * n_tetras << "\n";
    for (size_t i = 0; i < n_tetras; i++) {
        out_vol_msh << 4 << " ";
        out_vol_msh << tet_nodes[i][0] << " ";
        out_vol_msh << tet_nodes[i][1] << " ";
        out_vol_msh << tet_nodes[i][2] << " ";
        out_vol_msh << tet_nodes[i][3] << "\n";
    }
    out_vol_msh << "CELL_TYPES " << n_tetras << "\n";
    for (size_t i = 0; i < n_tetras; i++) {
        out_vol_msh << 10 << "\n";
    }
    out_vol_msh << "CELL_DATA " << n_tetras << "\n";
    out_vol_msh << "SCALARS "
                << "Materials int 1\n";
    out_vol_msh << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_tetras; i++) {
        out_vol_msh << tet_lab[i] << "\n";
    }
    out_vol_msh.close();

    std::ofstream out_srf_msh(std::string(vtkMshName + "_srfmsh.vtk").data());
    out_srf_msh << "# vtk DataFile Version 2.0\n";
    out_srf_msh << "Mesh data\n";
    out_srf_msh << "ASCII\n";
    out_srf_msh << "DATASET UNSTRUCTURED_GRID\n";
    out_srf_msh << "POINTS " << n_nodes << " double \n";
    for (size_t i = 0; i < n_nodes; i++) {
        out_srf_msh << std::setprecision(16) << nod_pos[i][0] << " ";
        out_srf_msh << std::setprecision(16) << nod_pos[i][1] << " ";
        out_srf_msh << std::setprecision(16) << nod_pos[i][2] << "\n";
    }
    out_srf_msh << "CELLS " << n_bc_mtrl_faces << " " << 4 * n_bc_mtrl_faces
                << "\n";
    for (size_t i = 0; i < n_faces; i++) {
        if (fac_lab[i] != -1) {
            out_srf_msh << 3 << " ";
            out_srf_msh << fac_nodes[i][0] << " ";
            out_srf_msh << fac_nodes[i][1] << " ";
            out_srf_msh << fac_nodes[i][2] << "\n";
        }
    }
    out_srf_msh << "CELL_TYPES " << n_bc_mtrl_faces << "\n";
    for (size_t i = 0; i < n_faces; i++) {
        if (fac_lab[i] != -1)
            out_srf_msh << 5 << "\n";
    }
    out_srf_msh << "CELL_DATA " << n_bc_mtrl_faces << "\n";
    out_srf_msh << "SCALARS "
                << "Boundaries int 1\n";
    out_srf_msh << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_faces; i++) {
        if (fac_lab[i] != -1)
            out_srf_msh << fac_lab[i] << "\n";
    }
    out_srf_msh.close();
    for (size_t i = 0; i < n_faces; i++)
        if (fac_lab[i] < -1)
            fac_lab[i] = -1;
}

std::vector<std::vector<double>> mdl_msh::tet_geo(size_t id) {
    std::vector<std::vector<double>> cGeo(4, std::vector<double>(3));
    cGeo[0] = nod_pos[tet_nodes[id][0]];
    cGeo[1] = nod_pos[tet_nodes[id][1]];
    cGeo[2] = nod_pos[tet_nodes[id][2]];
    cGeo[3] = nod_pos[tet_nodes[id][3]];
    return cGeo;
}

std::vector<std::vector<double>> mdl_msh::fac_geo(size_t id) {
    std::vector<std::vector<double>> cGeo(3, std::vector<double>(3));
    cGeo[0] = nod_pos[fac_nodes[id][0]];
    cGeo[1] = nod_pos[fac_nodes[id][1]];
    cGeo[2] = nod_pos[fac_nodes[id][2]];
    return cGeo;
}

std::vector<std::vector<double>> mdl_msh::fac_geo2(size_t id) {
    std::vector<double> v0 = nod_pos[fac_nodes[id][0]];
    std::vector<double> v1 = nod_pos[fac_nodes[id][1]];
    std::vector<double> v2 = nod_pos[fac_nodes[id][2]];
    for (unsigned int i = 0; i < 3; i++) {
        v1[i] -= v0[i];
        v2[i] -= v0[i];
    }
    std::vector<double> u = v1;
    std::vector<double> n(3), v(3);
    n[0] = v1[1] * v2[2];
    n[1] = v1[2] * v2[0];
    n[2] = v1[0] * v2[1];
    double norm2_u = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    double norm2_n = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (unsigned int i = 0; i < 3; i++) {
        u[i] /= norm2_u;
        n[i] /= norm2_n;
    }
    v[0] = n[1] * u[2];
    v[1] = n[2] * u[0];
    v[2] = n[0] * u[1];
    std::vector<std::vector<double>> cGeo2(3, std::vector<double>(2, 0.0));
    cGeo2[1][0] = v1[0] * u[0] + v1[1] * u[1] + v1[2] * u[2];
    cGeo2[2][0] = v2[0] * u[0] + v2[1] * u[1] + v2[2] * u[2];
    cGeo2[2][1] = v2[0] * v[0] + v2[1] * v[1] + v2[2] * v[2];
    return cGeo2;
}

std::vector<std::vector<double>> mdl_msh::edg_geo(size_t id) {
    std::vector<std::vector<double>> cGeo(2, std::vector<double>(3));
    cGeo[0] = nod_pos[edg_nodes[id][0]];
    cGeo[1] = nod_pos[edg_nodes[id][1]];
    return cGeo;
}

std::vector<double> mdl_msh::int_node(size_t id) {
    std::vector<size_t> nfac = fac_nodes[id];
    std::vector<size_t> ntet = tet_nodes[fac_adj_tet[id][0]];
    size_t intid;
    for (unsigned int i = 0; i < 4; i++) {
        bool found = true;
        intid = ntet[i];
        for (unsigned int j = 0; j < 3; j++)
            if (ntet[i] == nfac[j]) {
                found = false;
            }
        if (found) {
            break;
        }
    }
    return nod_pos[intid];
}

std::vector<double> mdl_msh::int_node(size_t id, size_t &ref_face) {
    std::vector<size_t> nfac = fac_nodes[id];
    std::vector<size_t> ntet = tet_nodes[fac_adj_tet[id][ref_face]];
    size_t intid = 0;
    for (unsigned int i = 0; i < 4; i++) {
        bool found = true;
        intid = ntet[i];
        for (unsigned int j = 0; j < 3; j++)
            if (ntet[i] == nfac[j]) {
                found = false;
            }
        if (found) {
            ref_face = i;
            break;
        }
    }
    return nod_pos[intid];
}

// double mdl_msh::tet_vol(size_t id) {
////         1      [ax bx cx dx]T
////    V = --- det [ay by cy dy]
////         6      [az bz cz dz]
////                [ 1  1  1  1]
//    std::vector<std::vector<double> > matrix = tet_geo(id);
//    matrix.resize(4,4);
//    matrix.col(3).fill(1.0);
//    return std::abs(arma::det(matrix]]/6.0;
//}
//
// double mdl_msh::tet_mean_edg_length(size_t id) {
//    std::vector<std::vector<double> > matrix = tet_geo(id);
//    return (arma::norm(matrix[1)-matrix[0),2) +
//            arma::norm(matrix[2)-matrix[0),2) +
//            arma::norm(matrix[3)-matrix[0),2) +
//            arma::norm(matrix[2)-matrix[1),2) +
//            arma::norm(matrix[3)-matrix[1),2) +
//            arma::norm(matrix[3)-matrix[2),2]] / 6.0;
//}
//
// double mdl_msh::tet_max_edg_length(size_t id) {
//    std::vector<std::vector<double> > matrix = tet_geo(id);
//    double max_edge = 0.0;
//    max_edge = std::max(max_edge, arma::norm(matrix[1)-matrix[0),2]];
//    max_edge = std::max(max_edge, arma::norm(matrix[2)-matrix[0),2]];
//    max_edge = std::max(max_edge, arma::norm(matrix[3)-matrix[0),2]];
//    max_edge = std::max(max_edge, arma::norm(matrix[2)-matrix[1),2]];
//    max_edge = std::max(max_edge, arma::norm(matrix[3)-matrix[1),2]];
//    max_edge = std::max(max_edge, arma::norm(matrix[3)-matrix[2),2]];
//    return max_edge;
//}

void mdl_msh::get_bounding_info() {
    bounding_box.resize(2);
    bounding_box[0].assign(3, DBL_MAX);
    bounding_box[1].assign(3, -DBL_MAX);
    for (size_t i = 0; i < nod_pos.size(); i++) {
        bounding_box[0][0] = std::min(bounding_box[0][0], nod_pos[i][0]);
        bounding_box[0][1] = std::min(bounding_box[0][1], nod_pos[i][1]);
        bounding_box[0][2] = std::min(bounding_box[0][2], nod_pos[i][2]);
        bounding_box[1][0] = std::max(bounding_box[1][0], nod_pos[i][0]);
        bounding_box[1][1] = std::max(bounding_box[1][1], nod_pos[i][1]);
        bounding_box[1][2] = std::max(bounding_box[1][2], nod_pos[i][2]);
    }
    double x_dim, y_dim, z_dim;
    x_dim = std::abs(bounding_box[1][0] - bounding_box[0][0]);
    y_dim = std::abs(bounding_box[1][1] - bounding_box[0][1]);
    z_dim = std::abs(bounding_box[1][2] - bounding_box[0][2]);
    max_dimension = std::max(max_dimension, x_dim);
    max_dimension = std::max(max_dimension, y_dim);
    max_dimension = std::max(max_dimension, z_dim);
}

double mdl_msh::get_geom_dim() {
    double dim = 0;
    for (unsigned int i = 0; i < nod_pos.size(); i++) {
        dim = std::max(dim, std::abs(nod_pos[i][0]));
        dim = std::max(dim, std::abs(nod_pos[i][1]));
        dim = std::max(dim, std::abs(nod_pos[i][2]));
    }
    dim *= 1e-1;
    for (int i = -15; i < 15; i++) {
        if (dim < pow(10, double(i)))
            return pow(10, double(i));
    }
    return 0;
}

mdl_src::mdl_src() {}

mdl_src::~mdl_src() {}

mdl_frm::mdl_frm() {}

mdl_frm::~mdl_frm() {
    std::vector<mdl_bc>().swap(bcs);
    std::vector<mdl_mtrl>().swap(mtrls);
}

void mdl_frm::clear() {
    bcs.clear();
    mtrls.clear();
}

void mdl_frm::update_msh_info(mdl_msh &msh) {
    if (strcmp(type.data(), "TETRA") == 0) {
        for (size_t i = 0; i < bcs.size(); i++) {
            bcs[i].faces.clear();
            size_t cLab = bcs[i].label;
            for (size_t fid = 0; fid < msh.n_faces; fid++) {
                if (cLab == msh.fac_lab[fid]) {
                    bcs[i].faces.push_back(fid);
                }
            }
        }
        std::vector<size_t> mtrl_tets(mtrls.size(), 0), LabMap(mtrls.size(), 0);
        for (size_t tid = 0; tid < msh.n_tetras; tid++) {
            mtrl_tets[msh.tet_lab[tid]]++;
        }
        for (size_t i = 0; i < mtrls.size(); i++) {
            mtrls[i].tetras.resize(mtrl_tets[i]);
            mtrl_tets[i] = 0;
        }
        for (size_t tid = 0; tid < msh.n_tetras; tid++) {
            size_t cLab = msh.tet_lab[tid];
            mtrls[cLab].tetras[mtrl_tets[cLab]++] = tid;
        }
    }
}

void mdl_frm::write_prj_file(std::string &name) {
    std::ofstream prj_out_file(std::string(name + ".fes").c_str(),
                               std::ios::out | std::ios::ate);
    prj_out_file << "#Formulation " << type << "\n";
    prj_out_file << "#Materials " << mtrls.size() << "\n";
    for (size_t i = 0; i < mtrls.size(); i++) {
        prj_out_file << mtrls[i].label << " " << mtrls[i].type << " "
                     << mtrls[i].epsr << " " << mtrls[i].mur << " "
                     << mtrls[i].sigma << " " << mtrls[i].tand << " "
                     << mtrls[i].name << "\n";
    }
    prj_out_file << "#Boundaries " << bcs.size() << "\n";
    for (size_t i = 0; i < bcs.size(); i++) {
        prj_out_file << bcs[i].label << " " << bcs[i].name << " " << bcs[i].type;
        if (strcmp(bcs[i].type.data(), "WavePort") == 0) {
            prj_out_file << " " << bcs[i].num_modes;
        } else if (strcmp(bcs[i].type.data(), "Impedance") == 0) {
            prj_out_file << " " << bcs[i].surf_impedance.real() << " "
                         << bcs[i].surf_impedance.imag();
        } else if (strcmp(bcs[i].type.data(), "LumpedPort") == 0) {
            prj_out_file << " " << bcs[i].lumped_impedance.real() << " "
                         << bcs[i].lumped_impedance.imag();
        } else if (strcmp(bcs[i].type.data(), "LumpedRLC") == 0) {
            prj_out_file << " " << bcs[i].R << " " << bcs[i].L << " " << bcs[i].C;
        } else if (strcmp(bcs[i].type.data(), "Voltage") == 0) {
            prj_out_file << " " << bcs[i].voltage;
        } else if (strcmp(bcs[i].type.data(), "Current") == 0) {
            prj_out_file << " " << bcs[i].current;
        }
        prj_out_file << "\n";
    }
    prj_out_file.close();
}

void mdl_frm::read_prj_file(std::string &name) {
    clear();
    std::ifstream frm_in_file(std::string(name + ".fes").c_str(), std::ios::in);
    std::string line;
    std::istringstream iss;
    unsigned int tmp_uint;
    double tmp_dbl;
    std::string tmp_str;
    if (frm_in_file.is_open()) {
        while (getline(frm_in_file, line)) {
            iss.clear();
            iss.str(line);
            iss >> tmp_str;
            if (strcmp(tmp_str.data(), "#Formulation") == 0) {
                iss >> type;
            }
            if (strcmp(tmp_str.data(), "#Materials") == 0) {
                iss >> tmp_uint;
                mtrls.resize(tmp_uint);
                for (size_t i = 0; i < mtrls.size(); i++) {
                    getline(frm_in_file, line);
                    iss.clear();
                    iss.str(line);
                    iss >> mtrls[i].label;
                    iss >> mtrls[i].type;
                    iss >> mtrls[i].epsr;
                    iss >> mtrls[i].mur;
                    iss >> mtrls[i].sigma;
                    iss >> mtrls[i].tand;
                    iss >> mtrls[i].name;
                }
            }
            if (strcmp(tmp_str.data(), "#Boundaries") == 0) {
                iss >> tmp_uint;
                bcs.resize(tmp_uint);
                for (size_t i = 0; i < bcs.size(); i++) {
                    getline(frm_in_file, line);
                    iss.clear();
                    iss.str(line);
                    iss >> bcs[i].label;
                    iss >> bcs[i].name;
                    iss >> bcs[i].type;
                    if (strcmp(bcs[i].type.data(), "WavePort") == 0) {
                        iss >> bcs[i].num_modes;
                    } else if (strcmp(bcs[i].type.data(), "Impedance") == 0) {
                        double real;
                        double imag;
                        iss >> real;
                        iss >> imag;
                        bcs[i].surf_impedance = std::complex<double>(real, imag);
                    } else if (strcmp(bcs[i].type.data(), "LumpedPort") == 0) {
                        double real;
                        double imag;
                        iss >> real;
                        iss >> imag;
                        bcs[i].lumped_impedance = std::complex<double>(real, imag);
                    } else if (strcmp(bcs[i].type.data(), "LumpedRLC") == 0) {
                        iss >> bcs[i].R;
                        iss >> bcs[i].L;
                        iss >> bcs[i].C;
                    } else if (strcmp(bcs[i].type.data(), "Voltage") == 0) {
                        iss >> bcs[i].voltage;
                    } else if (strcmp(bcs[i].type.data(), "Current") == 0) {
                        iss >> bcs[i].current;
                    }
                }
            }
        }
    }
    frm_in_file.close();
}

void mdl_core::create_tri_mesh() {
    msh.clear();
    msh.type = "TRIA";
    msh.n_nodes = sld.nodes.size();
    msh.nod_pos = sld.nodes;
    msh.n_faces = sld.faces.size();
    msh.fac_nodes.resize(msh.n_faces);
    msh.fac_edges.resize(msh.n_faces);
    msh.fac_lab.assign(msh.n_faces, 1);
    std::map<std::pair<size_t, size_t>, size_t> edg_map;
    size_t edg_cnt = 0;
    for (size_t i = 0; i < msh.n_faces; i++) {
        msh.fac_nodes[i] = sld.faces[i].polygons[0];
        std::sort(msh.fac_nodes[i].begin(), msh.fac_nodes[i].end());
        //		msh.fac_lab[i] = sld.faces_marker[i][0];
        if (edg_map.find(std::make_pair(msh.fac_nodes[i][0],
                                        msh.fac_nodes[i][1])) == edg_map.end())
            edg_map[std::make_pair(msh.fac_nodes[i][0], msh.fac_nodes[i][1])] =
                edg_cnt++;
        if (edg_map.find(std::make_pair(msh.fac_nodes[i][0],
                                        msh.fac_nodes[i][1])) == edg_map.end())
            edg_map[std::make_pair(msh.fac_nodes[i][0], msh.fac_nodes[i][2])] =
                edg_cnt++;
        if (edg_map.find(std::make_pair(msh.fac_nodes[i][0],
                                        msh.fac_nodes[i][1])) == edg_map.end())
            edg_map[std::make_pair(msh.fac_nodes[i][1], msh.fac_nodes[i][2])] =
                edg_cnt++;
    }
    std::cout << edg_cnt << "\n";
    msh.n_edges = edg_cnt;
    msh.edg_nodes.resize(edg_cnt);
    for (std::map<std::pair<size_t, size_t>, size_t>::iterator it =
             edg_map.begin();
         it != edg_map.end(); it++) {
        std::vector<size_t> edge(2);
        edge[0] = std::get<0>(it->first);
        edge[1] = std::get<1>(it->first);
        msh.edg_nodes[it->second] = edge;
    }
    msh.edg_lab.assign(edg_cnt, 1);
    msh.max_edg_marker = 1;
    msh.get_mesh_statistics();
}
