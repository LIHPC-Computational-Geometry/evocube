#pragma once
#include <Eigen/Core>
#include <iostream>

#include "disjointset.h"
#include "logging.h"

int updateCharts(const Eigen::MatrixXi& TT, const Eigen::VectorXi& labeling, Eigen::VectorXi& charts){
    int n_fs = TT.rows();
    DisjointSet ds(n_fs);

    for (int f = 0; f < n_fs; f ++) {
		int dim = labeling(f) / 2;
		for (int side = 0; side < 3; side ++) {
            int neighbor = TT(f, side);
            if (labeling(neighbor) == labeling(f))
                ds.merge(f, neighbor);
        }
	}

	std::vector<int> idmap;
	int nb_charts = ds.get_sets_id(idmap);
    if (idmap.size() != TT.rows()) coloredPrint("Unexpected sizes in updateChart", "red");

    charts.resize(idmap.size());
    for (int i = 0; i < idmap.size(); i++){ // TODO optimize ?
        charts(i) = idmap[i];
    }
    return nb_charts;
}

Eigen::VectorXi perChartLabels(const Eigen::VectorXi& charts, const Eigen::VectorXi& labeling){
    Eigen::VectorXi res(charts.maxCoeff() + 1);
    for (int i=0; i<charts.rows(); i++){
        res(charts(i)) = labeling(i); 
    }
    return res;
}

std::vector<std::vector<int>> chartAdjacency(const Eigen::MatrixXi& TT, const Eigen::VectorXi& charts){
    if (TT.rows() != charts.rows()) coloredPrint("Error, wrong sizes in chartAdjacency", "red");

    std::vector<std::vector<int>> adj;
    for (int c=0; c<=charts.maxCoeff(); c++){
        adj.push_back({});
    }
    for (int i=0; i<charts.rows(); i++){
        for (int neigh=0; neigh<3; neigh++){
            if (charts(i) != charts(TT(i, neigh))){
                adj[charts(i)].push_back(charts(TT(i, neigh)));
                adj[charts(TT(i, neigh))].push_back(charts(i));
            }
        }
    }

    for (int c=0; c<adj.size(); c++){
        std::sort(adj[c].begin(), adj[c].end());
        adj[c].erase(std::unique(adj[c].begin(), adj[c].end() ), adj[c].end());
    }

    /*
    for (int c=0; c<adj.size(); c++){        
        std::cout << c << ": ";
        for (int i=0; i<adj[c].size(); i++){
            std::cout << adj[c][i] << " ";
        }
        std::cout << std::endl;
    }*/

    return adj;
}

void computeBoundaries(const Eigen::MatrixXi& TT, const Eigen::MatrixXi& F, 
                       const Eigen::VectorXi& charts,
                       std::vector<std::vector<int>>& ordered_borders,
                       std::vector<std::pair<int,int>>& patches_per_border){
    std::vector<std::pair<int,int>> edges;
    std::vector<std::pair<int,int>> patches_per_edge;

    for (int i=0; i<TT.rows(); i++){
        for (int neigh=0; neigh<3; neigh++){
            if (TT(i, neigh) < i || TT(i, neigh) == -1) continue; // process each edge only once
            if (charts(i) != charts(TT(i, neigh))){ // this edge is on a border
                int ei1 = F(i, neigh);
                int ei2 = F(i, (neigh + 1) % 3);
                edges.push_back(std::make_pair(ei1, ei2));
                patches_per_edge.push_back(std::make_pair(charts(i), charts(TT(i, neigh))));
            }
        }
    }

    auto equalPair = [](std::pair<int, int> p1, std::pair<int, int> p2){ // Can be replaced with something from std? Maybe std::set?
        return (p1 == p2 || (p1.first == p2.second && p1.second == p2.first));
    };

    DisjointSet ds(edges.size());

    for (int i=0; i<edges.size(); i++){ // possible optim here, O(E2)
        for (int j=i+1; j<edges.size(); j++){
            if (equalPair(patches_per_edge[i], patches_per_edge[j])){
                if (edges[i].first == edges[j].first || edges[i].first == edges[j].second ||
                    edges[i].second == edges[j].first || edges[i].second == edges[j].second){
                    ds.merge(i, j);
                }
            }
        }
    }

    std::vector<int> idmap;
	int n_borders = ds.get_sets_id(idmap);

    std::vector<std::vector<std::pair<int,int>>> boundaries;
    for (int border=0; border<n_borders; border++){
        boundaries.push_back({});
    }

    patches_per_border.clear();
    patches_per_border.resize(n_borders);
    for (int i=0; i<idmap.size(); i++){
        boundaries[idmap[i]].push_back(edges[i]);
        patches_per_border[idmap[i]] = patches_per_edge[i];
    }

    // ordering
    /*
    for (int border=0; border<n_borders; border++){ // possible optim here, O(#B * max boundary size 2)
        // count vertex occurrences
        //Eigen::VectorXi occ = Eigen::VectorXi::Zero(boundaries[border].size() + 1);
        std::map<int,int> occ;
        for (std::pair<int,int> p: boundaries[border]){
            if (occ.count(p.first) == 0) occ[p.first] = 1;
            else occ[p.first] = occ[p.first] + 1;
            if (occ.count(p.second) == 0) occ[p.second] = 1;
            else occ[p.second] = occ[p.second] + 1;
            
        }

        // find starting vertex
        int start = -1;
        for (auto it = occ.begin(); it != occ.end(); it++){
            if (it->second == 1){
                start = it->first;
                break;
            }
        }
        if (start == -1){ // it's a loop, start anywhere
            start = boundaries[border][0].first;
        }
    }*/

    // --- From flagging.cpp in hexercise: ---
    ordered_borders.clear();
    for (int border=0; border<n_borders; border++){
        std::vector<std::pair<int,int>> vec = boundaries[border];
        // find endings: they're only on one edge of this boundary
        std::set<int> singles; 
        for (int i=0; i<vec.size(); i++) {
            for (int k: {vec[i].first, vec[i].second}){
                if (singles.count(k)){
                    singles.erase(singles.find(k));
                }
                else {
                    singles.insert(k);
                }
            }
        }

        std::vector<int> v_ids = {};
        int last_elem;
        if (singles.size() == 2){
            for (int i: singles){
                if (v_ids.empty()) v_ids.push_back(i);
                else last_elem = i;
            }
        }
        else if (singles.size() == 0){
            // it's a loop: start and end with same id
            v_ids.push_back(vec[0].first);
            v_ids.push_back(vec[0].second);
            last_elem = vec[0].first;
            vec.erase(vec.begin());
        } else {
            std::cout << "singles.size() = " << singles.size() << std::endl;
            for (int elem: singles) std::cout << elem << " ";
            std::cout << std::endl;
            std::cout << "Error identifying boundaries: Bad mesh/flagging ?" << std::endl;

            // TODO find out why following snippet is needed...
            for (int i: singles){
                if (v_ids.empty()) v_ids.push_back(i);
                else last_elem = i;
            }
        }

        // Go through consecutive edges until end is reached
        int last_added = v_ids[v_ids.size()-1];
        while (last_added != last_elem){
            for (int i=0; i<vec.size(); i++){
                std::pair<int,int> pair = vec[i];
                if (pair.first == last_added){
                    v_ids.push_back(pair.second);
                    last_added = pair.second;
                    vec.erase(vec.begin()+i);
                    break;
                }
                if (pair.second == last_added){
                    v_ids.push_back(pair.first);
                    last_added = pair.first;
                    vec.erase(vec.begin()+i);
                    break;
                }
            }
        }

        ordered_borders.push_back(v_ids);
    }
}