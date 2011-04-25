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

#ifndef TREE_graphics_H
#define TREE_graphics_H

#include "common.h"
#include "input.h"
#include <iostream>
#include <stack>
#include <set>
#include <vector>
#include <boost/bimap.hpp>
#include <boost/foreach.hpp>
#include "Board.h"
#include "tree_IO.h"

namespace NS_rootedtree {

enum {SLANTED, RECTANGULAR};

typedef LibBoard::Color color;

using namespace LibBoard;
using namespace std;

template<class TREE>
class treegraphic_internal {
    protected: TREE &t;
    protected: bimap_str2val &name2idx;
    protected: color &tree_color, &font_color;
    public: Board board;
    protected: vector<double> node_x, node_y;
    protected: unsigned int dim_x, dim_y;
    public: treegraphic_internal(
        TREE &_t, bimap_str2val &_name2idx,
        color _tree_color=color(0,0,0), color _font_color=color(0,0,0)
    ) : t(_t), name2idx(_name2idx),
        tree_color(_tree_color), font_color(_font_color)
    { }
    // compute node positions in a rectangular cladogram
    public: inline bool rectangular_cladogram() {
        node_x.clear(); node_x.resize(t.node_size());
        node_y.clear(); node_y.resize(t.node_size());
        dim_y = 0;
        dim_x = 0;
        {
            vector<unsigned int> lvl2(t.node_size()); BOOST_FOREACH(unsigned int &v, lvl2) v = 0;
            unsigned int lvl = 0;
            TREE_POSTORDER(v,t) {
                if (t.nodes[v].is_leaf()) {
                    lvl = 0;
                } else {
                    BOOST_FOREACH(const unsigned int &i, t.nodes[v].adjacent) {
                        if (lvl2[i]+1 > lvl) lvl = lvl2[i]+1;
                    }
                }
                lvl2[v] = lvl;
                lvl++;
            }
            TREE_FOREACH(v,t) node_x[v] = lvl2[v];
        }
        {
            unsigned int leafcounter = 0;
            vector<unsigned int> l(t.node_size());
            vector<unsigned int> r(t.node_size());
            TREE_DFS2(v,t) {
                switch (v.direction) {
                    case PREORDER: l[v.idx] = leafcounter; break;
                    case INORDER: if (!t.nodes[v.idx].is_leaf()) leafcounter++; break;
                    case POSTORDER: r[v.idx] = leafcounter; break;
                }
            }
            TREE_FOREACH(v,t) {
                node_y[v] = double(l[v] + r[v]) / 2;
            }
            dim_y = leafcounter;
            dim_x = 0; BOOST_FOREACH(const double &v, node_x) if (dim_x < v) dim_x = v;
        }
        return true;
    }
    // compute node positions in a slanted cladogram
    public: inline bool slanted_cladogram() {
        node_x.clear(); node_x.resize(t.node_size());
        node_y.clear(); node_y.resize(t.node_size());
        dim_y = 0;
        dim_x = 0;
        {
            unsigned int leafcounter = 0;
            vector<unsigned int> l(t.node_size());
            vector<unsigned int> r(t.node_size());
            TREE_DFS2(v,t) {
                switch (v.direction) {
                    case PREORDER: l[v.idx] = leafcounter; break;
                    case INORDER: if (!t.nodes[v.idx].is_leaf()) leafcounter++; break;
                    case POSTORDER: r[v.idx] = leafcounter; break;
                }
            }
            TREE_FOREACH(v,t) {
                node_y[v] = double(l[v] + r[v]) / 2;
                node_x[v] = r[v] - l[v];
            }
            dim_y = leafcounter;
            dim_x = 0; BOOST_FOREACH(const double &v, node_x) if (dim_x < v) dim_x = v;
        }
        return true;
    }
    // draw tree with slanted edges
    public: inline bool draw_slanted() {
        const double font_height = 10;
        const double scale = 10;
        TREE_FOREACH(v,t) { // flip tree
            node_x[v] = dim_x - node_x[v];
        }
        TREE_FOREACH(v,t) { // scale tree
            node_x[v] *= scale;
            node_y[v] *= scale;
        } dim_x *= scale; dim_y *= scale;
        TREE_FOREACH(v,t) { // translate tree
            // node_x[v] += font_height/2;
            node_y[v] += font_height/2;
        }
        board.setLineCap(Shape::RoundCap);
        TREE_POSTORDER(v,t) if (!t.nodes[v].is_root()) {
            const unsigned int &p = t.nodes[v].parent;
            board.setLineWidth(1).setPenColor(tree_color).drawLine(node_x[p],node_y[p], node_x[v],node_y[v], 4);
        }
        unsigned int max_str = 0;
        TREE_FOREACH(v,t) {
            bimap_str2val::right_const_iterator itr = name2idx.right.find(v);
            if (itr != name2idx.right.end()) {
                const unsigned int h = font_height;
                const unsigned int w = itr->second.length();
                const unsigned int x = node_x[v] + font_height * 0.5;
                const unsigned int y = node_y[v] - font_height/3;
                if (w > max_str) max_str = w;
                // board.fillRectangle(x, y + *h, w * 4*h, 5*h, 3);
                board << Text(x, y, itr->second, Fonts::Helvetica, "'Bookman Old Style',Verdana", h, font_color, 2);
            }
        }
        {
            const double x = 0;
            const double y = dim_y + font_height;
            const double w = dim_x + (font_height*0.5) + (max_str*font_height*0.75);
            const double h = dim_y + font_height;
            board.setPenColorRGBi(0xff,0xff,0xff,0).fillRectangle(x-10, y+10, w+20, h+20, 5);
        }
        return true;
    }
    // draw tree with rectangular edges
    public: inline bool draw_rectangluar() {
        const double font_height = 10;
        const double scale = 10;
        TREE_FOREACH(v,t) { // flip tree
            node_x[v] = dim_x - node_x[v];
        }
        TREE_FOREACH(v,t) { // scale tree
            node_x[v] *= scale;
            node_y[v] *= scale;
        } dim_x *= scale; dim_y *= scale;
        TREE_FOREACH(v,t) { // translate tree
            node_x[v] += scale;
            node_y[v] += font_height/2;
        }
        board.setLineCap(Shape::RoundCap);
        TREE_POSTORDER(v,t) if (!t.nodes[v].is_root()) {
            const unsigned int &p = t.nodes[v].parent;
            board.setLineWidth(1).setPenColor(tree_color).drawLine(node_x[p],node_y[p], node_x[p],node_y[v], 4);
            board.setLineWidth(1).setPenColor(tree_color).drawLine(node_x[p],node_y[v], node_x[v],node_y[v], 4);
        } else {
            board.setLineWidth(1).setPenColor(tree_color).drawLine(node_x[v] - scale,node_y[v], node_x[v],node_y[v], 4);
        }
        unsigned int max_str = 0;
        TREE_FOREACH(v,t) {
            bimap_str2val::right_const_iterator itr = name2idx.right.find(v);
            if (itr != name2idx.right.end()) {
                const unsigned int h = font_height;
                const unsigned int w = itr->second.length();
                const unsigned int x = node_x[v] + font_height * 0.5;
                const unsigned int y = node_y[v] - font_height/3;
                if (w > max_str) max_str = w;
                // board.fillRectangle(x, y + *h, w * 4*h, 5*h, 3);
                board << Text(x, y, itr->second, Fonts::Helvetica, "'Bookman Old Style',Verdana", h, font_color, 2);
            }
        }
        {
            const double x = 0;
            const double y = dim_y + font_height;
            const double w = dim_x + (font_height*0.5) + (max_str*font_height*0.75) + scale;
            const double h = dim_y + font_height;
            board.setPenColorRGBi(0xff,0xff,0xff,0).fillRectangle(x-10, y+10, w+20, h+20, 5);
        }
        return true;
    }
};

template<class TREE>
bool inline tree2svg(std::string filename, TREE &t, bimap_str2val &name2idx_t, unsigned int style=SLANTED, color tree_color=color(0,0,200), color font_color=color(0,0,0)) {
    treegraphic_internal<TREE> tg(t, name2idx_t, tree_color, font_color);
    switch (style) {
        case SLANTED: if (!tg.slanted_cladogram() || !tg.draw_slanted()) return false; break;
        case RECTANGULAR: if (!tg.rectangular_cladogram() || !tg.draw_rectangluar()) return false; break;
    }
    tg.board.saveSVG(filename.c_str());
    return true;
}

template<class TREE>
bool inline tree2eps(std::string filename, TREE &t, bimap_str2val &name2idx_t, unsigned int style=SLANTED, color tree_color=color(0,0,200), color font_color=color(0,0,0)) {
    treegraphic_internal<TREE> tg(t, name2idx_t, tree_color, font_color);
    switch (style) {
        case SLANTED: if (!tg.slanted_cladogram() || !tg.draw_slanted()) return false; break;
        case RECTANGULAR: if (!tg.rectangular_cladogram() || !tg.draw_rectangluar()) return false; break;
    }
    tg.board.saveEPS(filename.c_str());
    return true;
}

template<class TREE>
bool inline tree2fig(std::string filename, TREE &t, bimap_str2val &name2idx_t, unsigned int style=SLANTED, color tree_color=color(0,0,200), color font_color=color(0,0,0)) {
    treegraphic_internal<TREE> tg(t, name2idx_t, tree_color, font_color);
    switch (style) {
        case SLANTED: if (!tg.slanted_cladogram() || !tg.draw_slanted()) return false; break;
        case RECTANGULAR: if (!tg.rectangular_cladogram() || !tg.draw_rectangluar()) return false; break;
    }
    tg.board.saveFIG(filename.c_str());
    return true;
}

} // namespace end

#endif
