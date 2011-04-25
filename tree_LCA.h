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

#ifndef TREE_LCA_H
#define TREE_LCA_H

#include "common.h"
namespace RMQ {
#include "rmq.h"
}
#include "tree_traversal.h"
#include <vector>
#include <iostream>

namespace aw {

using namespace std;

// preprocess the LCA computation (build RMQ)
class LCA {
    protected: RMQ::INT *R; // first occurences in sequence
    protected: RMQ::VAL *E, *L; // sequence - E:nodes, L:levels
    protected: struct RMQ::rmqinfo *ri; // lookup table
    protected: unsigned int R_size, E_size, L_size;

    public: LCA() {
        init();
    }
    public: LCA(const LCA &r) { // copy constructor
        init();
        copy(r);
    }
    public: LCA& operator=(const LCA& r) { // assign operator
        if (this != &r) {
            free();
            copy(r);
        }
        return *this;
    }
    protected: inline void init() {
        R = NULL; E = NULL; L = NULL; ri = NULL;
        R_size = E_size = L_size = 0;
    }
    protected: inline void copy(const LCA &r) {
        if (r.E_size != 0) {E_size = r.E_size; E = new RMQ::VAL[E_size]; memcpy(E, r.E, E_size * sizeof(RMQ::VAL));}
        if (r.L_size != 0) {L_size = r.L_size; L = new RMQ::VAL[L_size]; memcpy(L, r.L, L_size * sizeof(RMQ::VAL));}
        if (r.R_size != 0) {R_size = r.R_size; R = new RMQ::INT[R_size]; memcpy(R, r.R, R_size * sizeof(RMQ::INT));}
        if (r.L_size != 0) {ri = RMQ::rm_query_preprocess(L, L_size);}
    }
    public: ~LCA() { this->free(); }
    protected: inline void free() {
        if (R != NULL) delete [] R;
        if (E != NULL) delete [] E;
        if (L != NULL) delete [] L;
        if (ri != NULL) rm_free(ri);
        init();
    }

    public: template<class TREE> inline bool create(TREE &tree) {
        free();
        vector<unsigned int> seq, lvl; seq.reserve(tree.node_size()); lvl.reserve(tree.node_size());
        TREE_INORDER2(v, tree) {
            seq.push_back(v.idx);
            lvl.push_back(v.lvl);
        }
        const unsigned int n = seq.size();
        E_size = n; E = new RMQ::VAL[E_size]; for (unsigned int i=0;i<n;i++) E[i] = seq[i];
        L_size = n; L = new RMQ::VAL[L_size]; for (unsigned int i=0;i<n;i++) L[i] = lvl[i];
        R_size = tree.node_size(); R = new RMQ::INT[R_size]; for (unsigned int i=n; i>0; i--) R[E[i-1]] = i-1;
        ri = RMQ::rm_query_preprocess(L, n);
        return true;
    }
    public: inline unsigned int lca(const unsigned int u, const unsigned int v) {
        if (u == v) return u;
        if (u == NONODE) return v;
        if (v == NONODE) return u;
        const RMQ::INT uidx = R[u];
        const RMQ::INT vidx = R[v];
        return E[rm_query(ri, uidx, vidx)];
    }
    public: template<class T> unsigned int lca(T leaves) {
        unsigned int r = 0;
        bool first = true;
        BOOST_FOREACH(const unsigned int &v, leaves) if (first) { r = v; first=false; } else r = lca(r,v);
        return r;
    }

    public: void clear() {
        free();
    }
};

} // namespace end

#endif
