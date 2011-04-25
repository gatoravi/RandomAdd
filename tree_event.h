/*
 * Copyright (C) 2010 Andre Wehe
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

#ifndef TREE_EVENT_H
#define TREE_EVENT_H

#include "common.h"
#include "tree_subtree_info.h"
#include "tree_node_distance.h"

namespace aw {

// tree if g is a gene duplication
inline bool IsDuplication(const unsigned int g_map, const unsigned int c0_map, const unsigned int c1_map) {
    return ((g_map == c0_map) || (g_map == c1_map));
}

// return the losses implied by a gene mapping to s
inline unsigned int GetLosses(
    const unsigned int s, const unsigned int s1, const unsigned int s2,
    aw::NodeDistance &nd2,
    aw::SubtreeInfoRooted<aw::Tree> &stree_info, aw::SubtreeInfoRooted<aw::Tree> &stree_info2,
    boost::unordered_map<unsigned int, unsigned int> &nodes
) {
    unsigned int loss = 0;
    if ((s == s1) && (s == s2)) {
        // no losses
    } else {
        const unsigned int &d1 = nd2.distance(s,s1);
        const unsigned int &d2 = nd2.distance(s,s2);
        unsigned int d = d1 + d2;
        if (d1 == 0) ++d; else --d;
        if (d2 == 0) ++d; else --d;
        loss += d;
        if (d1 != 0) {
            if (d2 != 0) {
                unsigned int v = s1;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == s) break;
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        ++nodes[v_sib];
                    }
                    v = pv;
                }
            } else {
                unsigned int v = s1;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        ++nodes[v_sib];
                    }
                    if (pv == s) break;
                    v = pv;
                }
            }
        }
        if (d2 != 0) {
            if (d1 != 0) {
                unsigned int v = s2;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == s) break;
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        ++nodes[v_sib];
                    }
                    v = pv;
                }
            } else {
                unsigned int v = s2;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        ++nodes[v_sib];
                    }
                    if (pv == s) break;
                    v = pv;
                }
            }
        }
    }
    return loss;
}

// return the losses implied by a gene mapping to s
inline unsigned int GetLosses(
    const unsigned int s, const unsigned int s1, const unsigned int s2,
    aw::NodeDistance &nd2,
    aw::SubtreeInfoRooted<aw::Tree> &stree_info, aw::SubtreeInfoRooted<aw::Tree> &stree_info2,
    std::vector<unsigned int> &nodes
) {
    unsigned int loss = 0;
    if ((s == s1) && (s == s2)) {
        // no losses
    } else {
        const unsigned int &d1 = nd2.distance(s,s1);
        const unsigned int &d2 = nd2.distance(s,s2);
        unsigned int d = d1 + d2;
        if (d1 == 0) ++d; else --d;
        if (d2 == 0) ++d; else --d;
        loss += d;
        if (d1 != 0) {
            if (d2 != 0) {
                unsigned int v = s1;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == s) break;
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        ++nodes[v_sib];
                    }
                    v = pv;
                }
            } else {
                unsigned int v = s1;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        ++nodes[v_sib];
                    }
                    if (pv == s) break;
                    v = pv;
                }
            }
        }
        if (d2 != 0) {
            if (d1 != 0) {
                unsigned int v = s2;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == s) break;
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        ++nodes[v_sib];
                    }
                    v = pv;
                }
            } else {
                unsigned int v = s2;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        ++nodes[v_sib];
                    }
                    if (pv == s) break;
                    v = pv;
                }
            }
        }
    }
    return loss;
}

// return the losses implied by a gene mapping to s
inline unsigned int GetLosses2(
    const unsigned int s, const unsigned int s1, const unsigned int s2,
    aw::NodeDistance &nd2,
    aw::SubtreeInfoRooted<aw::Tree> &stree_info, aw::SubtreeInfoRooted<aw::Tree> &stree_info2,
    std::vector<unsigned int> &nodes
) {
    unsigned int loss = 0;
    if ((s == s1) && (s == s2)) {
        // no losses
    } else {
        const unsigned int &d1 = nd2.distance(s,s1);
        const unsigned int &d2 = nd2.distance(s,s2);
        unsigned int d = d1 + d2;
        if (d1 == 0) ++d; else --d;
        if (d2 == 0) ++d; else --d;
        loss += d;
        if (d1 != 0) {
            if (d2 != 0) {
                unsigned int v = s1;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == s) break;
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        nodes.push_back(v_sib);
                    }
                    v = pv;
                }
            } else {
                unsigned int v = s1;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        nodes.push_back(v_sib);
                    }
                    if (pv == s) break;
                    v = pv;
                }
            }
        }
        if (d2 != 0) {
            if (d1 != 0) {
                unsigned int v = s2;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == s) break;
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        nodes.push_back(v_sib);
                    }
                    v = pv;
                }
            } else {
                unsigned int v = s2;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        nodes.push_back(v_sib);
                    }
                    if (pv == s) break;
                    v = pv;
                }
            }
        }
    }
    return loss;
}

// return the losses implied by a gene mapping to s
inline unsigned int GetLosses2(
    const unsigned int s, const unsigned int s1, const unsigned int s2,
    aw::NodeDistance &nd2,
    aw::SubtreeInfoRooted<aw::Tree> &stree_info, aw::SubtreeInfoRooted<aw::Tree> &stree_info2,
    std::vector<unsigned int> &nodes,
    std::vector<unsigned int> &edges
) {
    unsigned int loss = 0;
    if ((s == s1) && (s == s2)) {
        // no losses
    } else {
        const unsigned int &d1 = nd2.distance(s,s1);
        const unsigned int &d2 = nd2.distance(s,s2);
        unsigned int d = d1 + d2;
        if (d1 == 0) ++d; else --d;
        if (d2 == 0) ++d; else --d;
        loss += d;
        if (d1 != 0) {
            if (d2 != 0) {
                unsigned int v = s1;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    edges.push_back(v);
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == s) break;
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        nodes.push_back(v_sib);
                    }
                    v = pv;
                }
            } else {
                unsigned int v = s1;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    edges.push_back(v);
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        nodes.push_back(v_sib);
                    }
                    if (pv == s) break;
                    v = pv;
                }
            }
        }
        if (d2 != 0) {
            if (d1 != 0) {
                unsigned int v = s2;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    edges.push_back(v);
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == s) break;
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        nodes.push_back(v_sib);
                    }
                    v = pv;
                }
            } else {
                unsigned int v = s2;
                unsigned int x = stree_info2.parent(v);
                for (;;) {
                    edges.push_back(v);
                    const unsigned int &pv = stree_info.parent(v);
                    if (pv == x) {
                        x = stree_info2.parent(x);
                        const unsigned int &v_sib = stree_info.sibling_binary(v);
                        nodes.push_back(v_sib);
                    }
                    if (pv == s) break;
                    v = pv;
                }
            }
        }
    }
    return loss;
}

} // namespace end

#endif
