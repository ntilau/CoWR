#include "mdl_core.h"


void mdl_core::create_tri_mesh()
{
    msh.clear();
    msh.type = "TRIA";
    msh.n_nodes = sld.nodes.size();
    msh.nod_pos = sld.nodes;
    msh.n_faces = sld.faces.size();
    msh.fac_nodes.resize(msh.n_faces);
    msh.fac_edges.resize(msh.n_faces);
    msh.fac_lab.assign(msh.n_faces, 1);
    map<pair<size_t, size_t>, size_t> edg_map;
    size_t edg_cnt = 0;
    for (size_t i = 0; i < msh.n_faces; i++)
    {
        msh.fac_nodes[i] = sld.faces[i].polygons[0];
        sort(msh.fac_nodes[i].begin(), msh.fac_nodes[i].end());
        //		msh.fac_lab[i] = sld.faces_marker[i][0];
        if (edg_map.find(make_pair(msh.fac_nodes[i][0],
                                        msh.fac_nodes[i][1])) == edg_map.end())
            edg_map[make_pair(msh.fac_nodes[i][0], msh.fac_nodes[i][1])] =
                edg_cnt++;
        if (edg_map.find(make_pair(msh.fac_nodes[i][0],
                                        msh.fac_nodes[i][1])) == edg_map.end())
            edg_map[make_pair(msh.fac_nodes[i][0], msh.fac_nodes[i][2])] =
                edg_cnt++;
        if (edg_map.find(make_pair(msh.fac_nodes[i][0],
                                        msh.fac_nodes[i][1])) == edg_map.end())
            edg_map[make_pair(msh.fac_nodes[i][1], msh.fac_nodes[i][2])] =
                edg_cnt++;
    }
    cout << edg_cnt << "\n";
    msh.n_edges = edg_cnt;
    msh.edg_nodes.resize(edg_cnt);
    for (map<pair<size_t, size_t>, size_t>::iterator it =
             edg_map.begin();
         it != edg_map.end(); it++)
    {
        vector<size_t> edge(2);
        edge[0] = get<0>(it->first);
        edge[1] = get<1>(it->first);
        msh.edg_nodes[it->second] = edge;
    }
    msh.edg_lab.assign(edg_cnt, 1);
    msh.max_edg_marker = 1;
    msh.get_mesh_statistics();
}
