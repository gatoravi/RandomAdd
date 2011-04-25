/*
 * Copyright (C) 2009 Andre Wehe
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#ifndef TREE_EDGE_WEIGHTS_H
#define TREE_EDGE_WEIGHTS_H

#include "common.h"
#include "tree_traversal.h"
#include "tree_IO.h"
#include <vector>
#include <limits.h>
#include <boost/foreach.hpp>

// uncomment the next line to use boost::unordered_map instead of std::map
#define MAP23

#ifdef MAP23
#include <map>
#else
#include <boost/unordered_map.hpp>
#endif

namespace aw {

using namespace std;

class EdgeWeights {
    public: EdgeWeights() { default_value = 1.0;}
    public: ~EdgeWeights() { this->free(); }
    // convert weights from nodes to their parent edges (rooted tree required)
    public: template<class TREE> inline bool create(TREE &tree, idx2weight &node_weights) {
        unsigned int parent = UINT_MAX;
        TREE_DFS2(v,tree) {
            if (v.direction == PREORDER) {
                if (parent != UINT_MAX) set(v.idx, parent, node_weights[v.idx][0]);
            }
            parent = v.idx;
        }
//         BOOST_FOREACH(edge2double::value_type i, weights) {
//             cout << '(' << i.first.u << ',' << i.first.v << ") " << get(i.first.u,i.first.v) << endl;
//         }
        return true;
    }
    // convert weights from nodes to their parent edges (explicit root required)
    public: template<class TREE> inline bool create(TREE &tree, idx2weight &node_weights, const unsigned int root) {
        const unsigned int savedroot = tree.root; tree.root = root;
        bool ret = create(tree, node_weights);
        tree.root = savedroot;
        return ret;
    }
    public: inline double get(const unsigned int u, const unsigned int v) {
        unsigned int l = (u < v) ? u : v;
        unsigned int h = (u < v) ? v : u;
        edge2double::iterator itr = weights.find(NodePair(l,h));
        if (itr == weights.end()) return default_value;
        return itr->second;
    }
    public: inline void set(const unsigned int u, const unsigned int v, double val) {
        unsigned int l = (u < v) ? u : v;
        unsigned int h = (u < v) ? v : u;
        weights[NodePair(l,h)] = val;
    }
    public: inline bool exist(const unsigned int u, const unsigned int v) {
        unsigned int l = (u < v) ? u : v;
        unsigned int h = (u < v) ? v : u;
        edge2double::iterator itr = weights.find(NodePair(l,h));
        return itr != weights.end();
    }
    public: template<class TREE> inline bool convert2_idx2weight(TREE &tree, idx2weight &node_weights) {
        if (tree.root == TREE::NOROOT) ERROR_return("rooted tree expected");
        unsigned int parent = UINT_MAX;
        TREE_DFS2(v,tree) {
            if (v.direction == PREORDER) {
                if (parent != UINT_MAX) {
                    if (exist(v.idx, parent)) node_weights[v.idx] = get(v.idx, parent);
                }
            }
            parent = v.idx;
        }
        return true;
    }
    protected: inline void free() {
        weights.clear();
    }
    public: class NodePair {
        public: unsigned int u,v;
        public: NodePair() {}
        public: NodePair(const unsigned int _u, const unsigned int _v) : u(_u), v(_v) {}
        public: inline bool operator<(const NodePair& o) const { return (u < o.u) || ((u == o.u) && (v < o.v)); }
        public: inline bool operator==(const NodePair& o) const { return (u == o.u) && (v == o.v); }
    };
    #ifdef MAP23
    protected: typedef std::map<NodePair,double> edge2double;
    #else
    friend std::size_t hash_value(NodePair const& p);
    protected: typedef boost::unordered_map<NodePair,double> edge2double;
    #endif
    protected: edge2double weights;
    public: double default_value;
};
#ifndef MAP23
std::size_t hash_value(EdgeWeights::NodePair const& p) {
    std::size_t seed = 0;
    boost::hash_combine(seed, p.u);
    boost::hash_combine(seed, p.v);
    return seed;
}
#endif

} // namespace end

#endif
