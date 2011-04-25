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

#ifndef TREE_DUPLICATIONS_H
#define TREE_DUPLICATIONS_H

#include "common.h"
#include "tree_traversal.h"
#include "tree_LCA.h"
#include "tree_name_map.h"
#include "tree_LCA_mapping.h"
#include "tree_subtree_info.h"
#include <limits.h>
#include <boost/random.hpp>

namespace aw {

using namespace std;

// compute the gene duplications induced by a gene tree
template<class TREE>
unsigned int compute_duplications(TREE &g_tree, LCAmapping &g_map) {
    unsigned int dups = 0;
    TREE_POSTORDER2(v,g_tree) {
        if (!g_tree.is_leaf(v.idx)) {
            unsigned int ch[2];
            g_tree.children(v.idx,v.parent,ch);
            const unsigned int &c_map_0 = g_map.mapping(ch[0]);
            const unsigned int &c_map_1 = g_map.mapping(ch[1]);
            if ((c_map_0 == NONODE) || (c_map_1 == NONODE)) continue;
            const unsigned int &s_map = g_map.mapping(v.idx);
            if ((c_map_0 == s_map) || (c_map_1 == s_map)) {
                ++dups;
            }
        }
    }
    return dups;
}

// compute the gene duplications induced by multiple gene trees
template<class TREE>
inline unsigned int compute_duplications(std::vector<TREE> &g_trees, std::vector<LCAmapping> &g_maps) {
    unsigned int dups = 0;
    for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
        dups += compute_duplications(g_trees[i],g_maps[i]);
    }
    return dups;
}

// compute the gene duplications induced by a gene tree
template<class STREE,class GTREE>
inline unsigned int compute_duplications(STREE &s_tree, idx2name &s_names, GTREE &g_tree, idx2name &g_names) {
    TaxaMap taxamap;
    taxamap.insert(s_names);
    taxamap.insert(g_names);
    TreetaxaMap s_nmap; s_nmap.create(s_names,taxamap);
    TreetaxaMap g_nmap; g_nmap.create(g_names,taxamap);
    LCA lca; lca.create(s_tree);
    LCAmapping lca_map; lca_map.create(lca,s_nmap,g_tree,g_nmap);
    const unsigned int dups = compute_duplications(g_tree,lca_map);
    return dups;
}

// compute the gene duplications induced by multiple gene trees
template<class STREE,class GTREE>
inline unsigned int compute_duplications(STREE &s_tree, TreetaxaMap &s_nmap, std::vector<GTREE> &g_trees, std::vector<TreetaxaMap> &g_nmaps) {
    LCA s_lca; s_lca.create(s_tree);
    unsigned int dups = 0;
    for (unsigned int i=0,iEE=g_nmaps.size(); i<iEE; ++i) {
        GTREE &g_tree = g_trees[i];
        TreetaxaMap &g_nmap = g_nmaps[i];
        LCAmapping lca_map; lca_map.create(s_lca,s_nmap,g_nmap,g_tree);
        dups += compute_duplications(g_tree,lca_map);
    }
    return dups;
}

// compute the gene duplications induced by multiple gene trees
template<class STREE,class GTREE>
inline unsigned int compute_duplications(STREE &s_tree, idx2name &s_names, std::vector<GTREE> &g_trees, std::vector<idx2name> &g_names) {
    TaxaMap taxamap;
    taxamap.insert(s_names);
    TreetaxaMap s_nmap; s_nmap.create(s_names,taxamap);
    std::vector<TreetaxaMap> g_nmaps(g_names.size());
    for (unsigned int i=0,iEE=g_names.size(); i<iEE; ++i) {
        idx2name &n = g_names[i];
        taxamap.insert(n);
        g_nmaps[i].create(n,taxamap);
    }
    const unsigned int dups = compute_duplications(s_tree,s_nmap,g_trees,g_nmaps);
    return dups;
}

#ifndef AW_RANDOMGEN
boost::mt19937 rng;
#endif

// find the best SPR move on the subtree - there can be multiple equal ones, then only one of them is returned
// return
//   location = edge(u,v) for the location with lowest duplications
//   duplications = lowest duplications
// not thread-safe
template<class STREE,class GTREE>
bool bestSPRlocation(
    const unsigned int subtree, const unsigned int subtree_parent, STREE &s_tree,
    std::vector<GTREE> &g_trees, std::vector<aw::LCAmapping> &g_lmaps,
    std::pair<unsigned int, unsigned int> &location, unsigned int &duplications
) {
    const unsigned int &subtree_left = subtree;
    const unsigned int current_root = s_tree.root;
    if (subtree_left == current_root) return false; // someone wants me to place the tree inside the tree -> can't do that
    const bool move = !(subtree_parent == current_root); // unless the subtree is already placed at the root ...
    std::vector<unsigned int> subtree_parent_adj;
    if (move) { // move the subtree to the root
        s_tree.children(subtree_parent,subtree_left,subtree_parent_adj);
        // if (subtree_parent_adj.size() != 2) ERROR_exit("something is wrong here");
        s_tree.spr_to_root(subtree_left,subtree_parent);
    }
    { // Bansal, Eulenstein, Wehe algorithm to determine best SPR for subtree_left
        // update current LCA mapping
        aw::LCA s_lca; s_lca.create(s_tree); // compute LCAs for the species tree
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) g_lmaps[i].update_LCA_internals(s_lca,g_trees[i]); // update the LCA mapping for internal nodes of the gene trees
        static aw::SubtreeInfoRooted<aw::Tree> s_info; s_info.create(s_tree);
        // determine locations and gene duplication changes
        const unsigned int &s_root = subtree_parent;
        const unsigned int subtree_right = *s_tree.children(subtree_parent,subtree_left).begin();
        std::vector<unsigned int> dups_inc(s_tree.node_size(),0);
        std::vector<unsigned int> dups_dec(s_tree.node_size(),0);
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            aw::Tree &g_tree = g_trees[i];
            static std::vector<unsigned int> g_lmap2; // store the secondary LCA mappings of a gene tree (reuseable)
            if (g_tree.node_size() > g_lmap2.size()) g_lmap2.resize(g_tree.node_size());
            aw::LCAmapping &g_lmap = g_lmaps[i];
            TREE_DFS2(v,g_tree) {
                switch (v.direction) {
                    case aw::PREORDER: {
                        const unsigned int &v_map = g_lmap.mapping(v.idx);
                        if (v_map != s_root) { // no need to traverse further into the subtree
                            v.skip();
                        }
                    } break;
                    case aw::POSTORDER: {
                        const unsigned int &v_map = g_lmap.mapping(v.idx);
                        if (v_map == s_root) {
                            unsigned int ch[2]; g_tree.children(v.idx,v.parent,ch);
                            const unsigned int ch_map[2] = { g_lmap.mapping(ch[0]),g_lmap.mapping(ch[1]) };
                            bool case0;
                            // missing lca mapping, usually caused by missing taxa in the species tree
                            // no gene duplications will ever be lost or gained
                            if ((case0 = ch_map[0] == NONODE) || (ch_map[1] == NONODE)) {
                                g_lmap2[v.idx] = g_lmap2[ch[case0 ? 1 : 0]];
                                break;
                            }
                            const bool to_root_0 = (ch_map[0] == s_root);
                            const bool to_root_1 = (ch_map[1] == s_root);
                            // parent and both children map to the root
                            // no gene duplications will every be lost or gained
                            if (to_root_0 && to_root_1) {
//                                 const unsigned int ch_map2[2] = { g_lmap2[ch[0]],g_lmap2[ch[1]] };
//                                 const unsigned int v_map2 = s_lca.lca(ch_map2[0],ch_map2[1]);
                                const unsigned int &ch_map2_0 = g_lmap2[ch[0]];
                                const unsigned int &ch_map2_1 = g_lmap2[ch[1]];
                                const unsigned int v_map2 = s_lca.lca(ch_map2_0,ch_map2_1);
                                g_lmap2[v.idx] = v_map2;
                                break;
                            }
                            const bool in_r_0 = s_info.is_contained(ch_map[0],subtree_right);
                            const bool in_r_1 = s_info.is_contained(ch_map[1],subtree_right);
                            // parent and one child map to the root, one child maps into the right subtree
                            if ((case0 = to_root_0 && in_r_1) || (to_root_1 && in_r_0)) {
//                                 const unsigned int ch_map2[2] = {
//                                     case0 ? g_lmap2[ch[0]] : ch_map[0],
//                                     case0 ? ch_map[1] : g_lmap2[ch[1]]
//                                 };
//                                 const unsigned int &s_l = ch_map2[case0 ? 0 : 1]; // virtually maps to the right subtree
//                                 const unsigned int &s_r = ch_map2[case0 ? 1 : 0]; // maps to the right subtree
//                                 const unsigned int &s_p = g_lmap2[v.idx] = s_lca.lca(ch_map2[0],ch_map2[1]); // virtually maps to the lca of both child mappings
                                const unsigned int &ch_map2_0 = case0 ? g_lmap2[ch[0]] : ch_map[0];
                                const unsigned int &ch_map2_1 = case0 ? ch_map[1] : g_lmap2[ch[1]];
                                const unsigned int &s_l = case0 ? ch_map2_0 : ch_map2_1; // virtually maps to the right subtree
                                const unsigned int &s_r = case0 ? ch_map2_1 : ch_map2_0; // maps to the right subtree
                                const unsigned int &s_p = g_lmap2[v.idx] = s_lca.lca(ch_map2_0,ch_map2_1); // virtually maps to the lca of both child mappings
                                const unsigned int &s_pp = s_info.parent(s_p); // parent of s_p
                                unsigned int s_pch[2]; s_tree.children(s_p, s_pp, s_pch); // children of s_p
                                if (s_r == s_p) { // duplication remains because of the right child mapping
                                } else
                                if (s_info.is_contained(s_l,s_pch[0])) {
                                    const unsigned int &s_ll = s_pch[0]; // left child of s_p
                                    ++dups_dec[s_ll];
                                } else
                                if (s_info.is_contained(s_l,s_pch[1])) {
                                    const unsigned int &s_ll = s_pch[1]; // left child of s_p
                                    ++dups_dec[s_ll];
                                }
                                break;
                            }
                            const bool in_l_0 = s_info.is_contained(ch_map[0],subtree_left);
                            const bool in_l_1 = s_info.is_contained(ch_map[1],subtree_left);
                            // parent maps to the root, one child maps into the left subtree, one child maps into the right subtree
                            // then gene duplication occurs when moving below the mapped node in the right subtree
                            if ((case0 = in_l_0 && in_r_1) || (in_l_1 && in_r_0)) {
                                const unsigned int &v_map2 = ch_map[case0 ? 1 : 0];
                                g_lmap2[v.idx] = v_map2;
                                ++dups_inc[v_map2];
                                break;
                            }
                            // parent and one child map to the root, one child maps into the left subtree
                            // then no gene duplications will every be lost
                            if ((case0 = to_root_0 && in_l_1) || (to_root_1 && in_l_0)) {
                                const unsigned int &v_map2 = g_lmap2[ch[case0 ? 0 : 1]];
                                g_lmap2[v.idx] = v_map2;
                                break;
                            }
                            // something is wrong
                            ERROR_exit("something is wrong");
                        }
                    } break;
                    default: break;
                }
            }
        }
        { // accumulate gene duplication changes
            unsigned int d = duplications = compute_duplications(g_trees,g_lmaps); // compute initial gene duplications
            // static std::vector<std::pair<unsigned int,unsigned int> > candidates;
            static util::vector<std::pair<unsigned int,unsigned int> > candidates; candidates.resize(s_tree.node_size());
            location = std::pair<unsigned int,unsigned int>(subtree_right,s_root);
            for (aw::Tree::iterator_dfs v=s_tree.begin_dfs(subtree_right,s_root),vEE=s_tree.end_dfs(); v!=vEE; ++v) {
                switch (v.direction) {
                    case aw::PREORDER: {
                        d -= dups_dec[v.idx];
                        if (d < duplications) {
                            duplications = d;
                            candidates.clear();
                        }
                        if (d == duplications) {
                            candidates.push_back(std::pair<unsigned int,unsigned int>(v.idx,v.parent));
                        }
                        // { // brute force verification of the algorithm
                        //     aw::Tree s_tree2 = s_tree;
                        //     s_tree2.spr_from_root(subtree_left,s_info.sibling_binary(subtree_left),v.idx,v.parent);
                        //     aw::LCA s_lca2; s_lca2.create(s_tree2);
                        //     std::vector<aw::LCAmapping> g_lmaps2 = g_lmaps;
                        //     for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) g_lmaps2[i].update_LCA_internals(s_lca2,g_trees[i]); // update the LCA mapping for internal nodes of the gene trees
                        //     unsigned int dd=compute_duplications(g_trees,g_lmaps2);
                        //     if (d != dd) std::cout << v.idx << ':' << d << '(' << dd << ") ";
                        // }
                        d += dups_inc[v.idx];
                    } break;
                    case aw::POSTORDER: {
                        d -= dups_inc[v.idx];
                        d += dups_dec[v.idx];
                    } break;
                    default: break;
                }
            }
            if (candidates.empty()) ERROR_exit("something is wrong here");
            boost::uniform_int<> range(0,candidates.size()-1);
            boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(aw::rng, range);
            location = candidates[die()];
            candidates.clear();
        }
    }
    if (move) { // move the subtree back to its original location
        s_tree.spr_from_root(subtree_left,current_root,subtree_parent_adj[0],subtree_parent_adj[1]);
    }
    return true;
}

} // namespace end

#endif
