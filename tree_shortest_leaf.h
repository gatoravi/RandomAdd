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

#ifndef TREE_SHORTEST_LEAF_H
#define TREE_SHORTEST_LEAF_H

#include "common.h"
#include "tree_traversal.h"
#include "tree_edge_weights.h"
#include <limits.h>
#include <map>
#include <stack>
#include <boost/foreach.hpp>

namespace aw {

using namespace std;

class ShortestLeaf {
    public: const static unsigned int UNROOTED = UINT_MAX;
    public: ShortestLeaf() { }
    public: ~ShortestLeaf() { this->free(); }
    protected: inline void free() {
        closestleaves.clear();
    }
    // precompute shortest leaves
    public: template<class TREE> inline bool create(TREE &tree, EdgeWeights &weights) {
        this->free();
        TREE_FOREACHLEAF(v,tree) if (tree.degree(v) == 1) {
            const unsigned int p = *(tree.adjacent(v).begin());
            set(v,p,v,0.0);
        }
        TREE_FOREACHLEAF(v,tree) if (tree.degree(v) == 1) {
            std::stack<Subtree> todo;
            todo.push(Subtree(v,UNROOTED));
            while (!todo.empty()) {
                Subtree s = todo.top();
                unsigned int leaf = UINT_MAX;
                double len = UNROOTED;
                BOOST_FOREACH(const unsigned int &i, tree.adjacent(s.u)) if (i != s.v) {
                    const unsigned int &u = i, &v = s.u;
                    unsigned int leaf2;
                    double len2;
                    if (get(i,s.u,leaf2,len2)) {
                        len2 += weights.get(u,v);
                        if (len2 < len) {
                            leaf = leaf2;
                            len = len2;
                        }
                    } else {
                        todo.push(Subtree(u,v));
                    }
                }
                if (todo.top() == s) {
                    set(s.u,s.v,leaf,len);
                    todo.pop();
                }
            }
        }
//         BOOST_FOREACH(subtree2shortestleaf::value_type i, closestleaves) {
//             unsigned int leaf;
//             double len;
//             get(i.first.u,i.first.v,leaf,len);
//             cout << '(' << i.first.u << "->" << i.first.v << ") " << leaf << '[' << len << ']' << endl;
//         }
        return true;
    }
    // get shortest leaf and paths length
    // subtree root u, parent v (UNROOTED for no parent)
    public: inline bool get(const unsigned int u, const unsigned int v, unsigned int &leaf, double &len) {
        subtree2shortestleaf::iterator itr = closestleaves.find(Subtree(u,v));
        if (itr == closestleaves.end()) return false;
        leaf = itr->second.leaf;
        len = itr->second.len;
        return true;
    }
    // set shortest leaf and path length
    // subtree root u, parent v (UNROOTED for no parent)
    public: inline void set(const unsigned int u, const unsigned int v, const unsigned int leaf, const double len) {
        closestleaves[Subtree(u,v)] = ClosestLeaf(leaf,len);
    }
    protected: class ClosestLeaf {
        public: unsigned int leaf;
        public: double len;
        public: ClosestLeaf() {}
        public: ClosestLeaf(const unsigned int _leaf, const double _len) : leaf(_leaf), len(_len) {}
    };
    protected: class Subtree {
        public: unsigned int u,v;
        public: Subtree() {}
        public: Subtree(const unsigned int _u, const unsigned int _v) : u(_u), v(_v) {}
        public: inline bool operator<(const Subtree& o) const { return (u < o.u) || ((u == o.u) && (v < o.v)); }
        public: inline bool operator==(const Subtree& o) const { return (u == o.u) && (v == o.v); }
    };
    protected: typedef std::map<Subtree,ClosestLeaf> subtree2shortestleaf;
    protected: subtree2shortestleaf closestleaves;
};

} // namespace end

#endif
