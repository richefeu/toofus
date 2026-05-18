#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace lv {

// ============================================================================
// Vec2
// ============================================================================

template <typename T> struct Vec2 {
  T x{}, y{};
};

template <typename T> inline bool eq(const Vec2<T> &a, const Vec2<T> &b) {
  constexpr T eps = std::is_floating_point<T>::value ? T(1e-9) : T(0);
  return std::abs(a.x - b.x) <= eps && std::abs(a.y - b.y) <= eps;
}

// ============================================================================
// LaguerreBuilder
// ============================================================================

template <typename T> struct LaguerreBuilder {

  struct Site {
    T x{}, y{}, r{};
    Site(T x, T y, T r) : x(x), y(y), r(r) {}
  };

  using Poly = std::vector<Vec2<T>>;

  std::vector<Site> sites;
  std::vector<Poly> cells;
  int n_real{0};

  void add_site(T x, T y, T r) {
    sites.emplace_back(x, y, r);
  }

  T bbox_xmin_{}, bbox_ymin_{}, bbox_xmax_{}, bbox_ymax_{};
  bool use_custom_bbox_{false};

  void set_bbox(T xmin, T ymin, T xmax, T ymax) {
    bbox_xmin_       = xmin;
    bbox_ymin_       = ymin;
    bbox_xmax_       = xmax;
    bbox_ymax_       = ymax;
    use_custom_bbox_ = true;
  }

  // Fit bbox to all sites (including their radius).
  void fit_bbox() {
    if (sites.empty()) {
      use_custom_bbox_ = false;
      return;
    }
    T xmin = sites[0].x - sites[0].r, ymin = sites[0].y - sites[0].r;
    T xmax = sites[0].x + sites[0].r, ymax = sites[0].y + sites[0].r;
    for (const auto &s : sites) {
      xmin = std::min(xmin, s.x - s.r);
      xmax = std::max(xmax, s.x + s.r);
      ymin = std::min(ymin, s.y - s.r);
      ymax = std::max(ymax, s.y + s.r);
    }
    set_bbox(xmin, ymin, xmax, ymax);
  }

  // Sutherland-Hodgman clip against half-plane ax*x + ay*y <= c.
  Poly clip(const Poly &poly, T ax, T ay, T c) {
    Poly out;
    int n = (int)poly.size();
    if (n == 0) return out;
    for (int i = 0; i < n; ++i) {
      const auto &a = poly[i];
      const auto &b = poly[(i + 1) % n];
      T da          = ax * a.x + ay * a.y - c;
      T db          = ax * b.x + ay * b.y - c;
      if (da <= 0) out.push_back(a);
      if ((da < 0) != (db < 0)) {
        const T denom = da - db;
        if (std::abs(denom) > std::numeric_limits<T>::epsilon() * 1000) {
          T t = da / denom;
          out.push_back({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)});
        }
      }
    }
    return out;
  }

  // Laguerre bisector coefficients between site i and site j.
  void bisector(int i, int j, T &ax, T &ay, T &c) const {
    ax = 2 * (sites[j].x - sites[i].x);
    ay = 2 * (sites[j].y - sites[i].y);
    c  = sites[j].x * sites[j].x + sites[j].y * sites[j].y - sites[j].r * sites[j].r - sites[i].x * sites[i].x -
         sites[i].y * sites[i].y + sites[i].r * sites[i].r;
  }

  // O(n²) reference implementation.
  void compute() {
    int n = (int)sites.size();
    if (n_real == 0) n_real = n;
    cells.resize(n);
    if (!use_custom_bbox_) fit_bbox();
    Poly bbox = {
        {bbox_xmin_, bbox_ymin_}, {bbox_xmax_, bbox_ymin_}, {bbox_xmax_, bbox_ymax_}, {bbox_xmin_, bbox_ymax_}};
    for (int i = 0; i < n; ++i) {
      Poly poly = bbox;
      for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        T ax, ay, c;
        bisector(i, j, ax, ay, c);
        poly = clip(poly, ax, ay, c);
        if (poly.empty()) break;
      }
      cells[i] = poly;
    }
    
  }

  // O(n log n) practical: sort neighbours by distance, prune early.
  void compute_fast() {
    int n = (int)sites.size();
    if (n_real == 0) n_real = n;
    cells.resize(n);
    if (!use_custom_bbox_) fit_bbox();

    T xmin = bbox_xmin_, ymin = bbox_ymin_;
    T xmax = bbox_xmax_, ymax = bbox_ymax_;
    T diag2 = (xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin);

    Poly bbox = {{xmin, ymin}, {xmax, ymin}, {xmax, ymax}, {xmin, ymax}};

    for (int i = 0; i < n; ++i) {
      Poly poly  = bbox;
      const T xi = sites[i].x, yi = sites[i].y;

      // Sort all other sites by squared distance to site i.
      std::vector<std::pair<T, int>> neighbors;
      neighbors.reserve(n - 1);
      for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        T dx = sites[j].x - xi, dy = sites[j].y - yi;
        neighbors.push_back({dx * dx + dy * dy, j});
      }
      std::sort(neighbors.begin(), neighbors.end());

      // cell_radius2: squared max distance from site i to any cell vertex.
      // Start with bbox diagonal so the first clip is never pruned.
      T cell_radius2 = diag2;

      for (auto &[dist2, j] : neighbors) {
        // A site j whose bisector with i lies beyond the current cell
        // cannot clip further. Since the bisector is at dist/2 from i,
        // prune when (dist/2)² > cell_radius².
        if (dist2 > 4 * cell_radius2) break;

        T ax, ay, c;
        bisector(i, j, ax, ay, c);
        poly = clip(poly, ax, ay, c);
        if (poly.empty()) break;

        cell_radius2 = T(0);
        for (const auto &v : poly) {
          T vdx = v.x - xi, vdy = v.y - yi;
          cell_radius2 = std::max(cell_radius2, vdx * vdx + vdy * vdy);
        }
      }
      cells[i] = poly;
    }
  }

  // O(n log n) practical: spatial grid + ring expansion + early pruning.
  void compute_fast2() {
    int n = (int)sites.size();
    if (n_real == 0) n_real = n;
    cells.resize(n);
    if (!use_custom_bbox_) fit_bbox();

    T xmin = bbox_xmin_, ymin = bbox_ymin_;
    T xmax = bbox_xmax_, ymax = bbox_ymax_;
    T diag2 = (xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin);

    Poly bbox_poly = {{xmin, ymin}, {xmax, ymin}, {xmax, ymax}, {xmin, ymax}};

    // Build a G×G spatial grid. G≈sqrt(n) gives ~1 site per cell on average.
    int G    = std::max(4, (int)std::sqrt((double)n * 0.7));
    T cell_w = (xmax - xmin) / G;
    T cell_h = (ymax - ymin) / G;

    std::vector<std::vector<int>> grid(G * G);

    auto grid_cell = [&](T x, T y) -> std::pair<int, int> {
      int gx = std::min(G - 1, (int)((x - xmin) / cell_w));
      int gy = std::min(G - 1, (int)((y - ymin) / cell_h));
      return {gx, gy};
    };

    for (int i = 0; i < n; ++i) {
      auto [gx, gy] = grid_cell(sites[i].x, sites[i].y);
      grid[gy * G + gx].push_back(i);
    }

    for (int i = 0; i < n; ++i) {
      Poly poly  = bbox_poly;
      const T xi = sites[i].x, yi = sites[i].y;
      const T ri2 = sites[i].r * sites[i].r;

      auto [igx, igy] = grid_cell(xi, yi);
      T cell_radius2  = diag2;

      std::vector<std::pair<T, int>> candidates;

      // Expand outward ring by ring. Stop when the nearest point of the
      // next ring is already beyond the current cell radius.
      for (int ring = 0; ring <= G; ++ring) {
        if (ring > 1) {
          T d_min = (ring - 1) * std::min(cell_w, cell_h);
          if (d_min * d_min > 4 * cell_radius2) break;
        }

        candidates.clear();
        int gx0 = std::max(0, igx - ring), gx1 = std::min(G - 1, igx + ring);
        int gy0 = std::max(0, igy - ring), gy1 = std::min(G - 1, igy + ring);

        for (int gy = gy0; gy <= gy1; ++gy) {
          for (int gx = gx0; gx <= gx1; ++gx) {
            // Only process the border of the ring (interior already handled).
            if (gx != gx0 && gx != gx1 && gy != gy0 && gy != gy1) continue;
            for (int j : grid[gy * G + gx]) {
              if (j == i) continue;
              T dx = sites[j].x - xi, dy = sites[j].y - yi;
              candidates.push_back({dx * dx + dy * dy, j});
            }
          }
        }

        std::sort(candidates.begin(), candidates.end());

        for (auto &[dist2, j] : candidates) {
          if (dist2 > 4 * cell_radius2) continue;
          T ax, ay, c;
          bisector(i, j, ax, ay, c);
          poly = clip(poly, ax, ay, c);
          if (poly.empty()) break;
          cell_radius2 = T(0);
          for (const auto &v : poly) {
            T vdx = v.x - xi, vdy = v.y - yi;
            cell_radius2 = std::max(cell_radius2, vdx * vdx + vdy * vdy);
          }
        }
        if (poly.empty()) break;
      }
      cells[i] = poly;
    }
  }
  
  
  
  // Pre-processing: add a ring of ghost sites of radius R around the cloud,
  // using a gift-wrapping approach.  Call this BEFORE compute*().
  //
  // Algorithm:
  //   - Start from the rightmost site S0; place G0 tangent to its right.
  //   - At each step, given the current ghost Gk and incoming direction θk,
  //     find among all original sites the one Si that, combined with Gk,
  //     yields a new ghost position Gk+1 (tangent to both Gk and Si, exterior
  //     side) that requires the smallest CCW turn from θk.
  //   - Reject candidates that would penetrate any existing site or ghost.
  //   - Stop when the next ghost would coincide with G0.
  //
  // All ghost sites are appended to this->sites.  To distinguish them from
  // real sites, record n_real = sites.size() before calling this function.
  void add_surrounding_sites(T R) {
    n_real = (int)sites.size();
    if (n_real == 0) return;
    
    constexpr T pi  = T(3.14159265358979323846);
    constexpr T eps = T(1e-9);

    // Penetration check: does a candidate ghost at (cx,cy) radius R
    // overlap any existing site or already-placed ghost?
    auto penetrates = [&](T cx, T cy, int skip_site = -1) -> bool {
      for (int k = 0; k < (int)sites.size(); ++k) {
        if (k == skip_site) continue;
        T dx = sites[k].x - cx, dy = sites[k].y - cy;
        T min_d = sites[k].r + R;
        if (dx * dx + dy * dy < min_d * min_d - eps) return true;
      }
      return false;
    };

    // Normalise angle to [0, 2π).
    auto norm_angle = [&](T a) -> T {
      while (a < 0)    a += 2 * pi;
      while (a >= 2 * pi) a -= 2 * pi;
      return a;
    };

    // --- Initialisation: rightmost site ---
    int i0 = 0;
    for (int k = 1; k < n_real; ++k)
      if (sites[k].x + sites[k].r > sites[i0].x + sites[i0].r) i0 = k;

    const Site &s0 = sites[i0];
    T g0x = s0.x + s0.r + R, g0y = s0.y;
    sites.emplace_back(g0x, g0y, R);

    // Incoming direction for G0: we arrived from the right (angle = 0 → π,
    // i.e. the previous ghost would have been further right).  We set the
    // reference direction to "downward" (−π/2 in trigo) so the first turn
    // searches CCW from there.
    T theta = norm_angle(-pi / 2);  // pointing down = start searching CCW

    T gx = g0x, gy = g0y;          // current ghost centre

    for (;;) {
      // For each real site Si, compute the ghost position tangent to
      // current ghost Gk (radius R) and to Si (radius ri), on the exterior.
      //
      // |Gk Gk+1| = 2R  and  |Si Gk+1| = ri + R
      // => Gk+1 lies at the intersection of two circles.
      // We keep the intersection that is most CCW from theta.

      T   best_angle = -1;          // smallest CCW turn (initialised invalid)
      T   best_gx = 0, best_gy = 0;
      int best_si = -1;

      for (int si = 0; si < n_real; ++si) {
        const Site &s = sites[si];
        T d_gg = 2 * R;             // ghost-ghost tangency distance
        T d_sg = s.r + R;           // site-ghost tangency distance

        T dx = s.x - gx, dy = s.y - gy;
        T d2 = dx * dx + dy * dy;
        T d  = std::sqrt(d2);

        // Intersection of two circles: radii d_gg (centred on Gk) and
        // d_sg (centred on Si).
        if (d < eps) continue;
        if (d > d_gg + d_sg + eps) continue;   // too far
        if (d < std::abs(d_gg - d_sg) - eps) continue; // one inside the other

        // a = distance from Gk to the radical line
        T a   = (d_gg * d_gg - d_sg * d_sg + d2) / (2 * d);
        T h2  = d_gg * d_gg - a * a;
        if (h2 < 0) h2 = 0;
        T h   = std::sqrt(h2);

        // Midpoint on the line Gk→Si
        T mx  = gx + a * dx / d;
        T my  = gy + a * dy / d;

        // Two candidate ghost positions (perpendicular offsets)
        T px  = -dy / d * h, py = dx / d * h;

        T cands[2][2] = {{mx + px, my + py}, {mx - px, my - py}};

        for (auto &cand : cands) {
          T cx = cand[0], cy = cand[1];

          // Compute CCW angle from theta to direction Gk→Gk+1
          T angle_to = norm_angle(std::atan2(cy - gy, cx - gx));
          T turn     = norm_angle(angle_to - theta);

          // Prefer smallest CCW turn (gift-wrapping).
          if (best_si != -1 && turn >= best_angle) continue;

          // Reject if it penetrates any site or existing ghost
          // (allow contact with the pivot site si itself).
          if (penetrates(cx, cy, si)) continue;

          best_angle = turn;
          best_gx    = cx;
          best_gy    = cy;
          best_si    = si;
        }
      }

      if (best_si == -1) break;  // no valid candidate — open contour, stop

      // Check if we are back to G0 (closure).
      T dxclose = best_gx - g0x, dyclose = best_gy - g0y;
      if ((int)sites.size() > n_real + 1 &&
          dxclose * dxclose + dyclose * dyclose < R * R * T(1e-4))
        break;

      // Safety: never place more ghosts than 6 * n_real (sanity cap).
      if ((int)sites.size() - n_real > 6 * n_real) break;

      // Place the new ghost.
      sites.emplace_back(best_gx, best_gy, R);

      // Update state.
      theta = norm_angle(std::atan2(best_gy - gy, best_gx - gx));
      gx    = best_gx;
      gy    = best_gy;
    }
  }
  
  
  
  
};

// ============================================================================
// DCEL helpers (hash types)
// ============================================================================

struct PairHash {
  size_t operator()(const std::pair<int, int> &p) const {
    return std::hash<long long>()((long long)p.first << 32 | (unsigned)p.second);
  }
};

// Quantized 2D key used to deduplicate DCEL vertices.
struct DCELKey {
  long long x, y;
  bool operator==(const DCELKey &o) const {
    return x == o.x && y == o.y;
  }
};

struct DCELKeyHash {
  size_t operator()(const DCELKey &k) const {
    return std::hash<long long>()(k.x * 2654435761LL ^ k.y * 2246822519LL);
  }
};

// ============================================================================
// DCEL (Doubly Connected Edge List)
// ============================================================================

template <typename T> struct DCEL {

  struct Vertex {
    Vec2<T> p;
  };

  // Halfedges are stored in pairs: index 2k and 2k+1 are always twins.
  struct HalfEdge {
    int origin = -1; // start vertex index
    int twin   = -1; // opposite halfedge index
    int face   = -1; // face to the left (-1 for boundary)
    int next   = -1; // next halfedge around the face
  };

  // One face per Laguerre cell.
  struct Face {
    int site = -1; // generating site index
    int edge = -1; // any halfedge on this face's boundary
  };

  std::vector<Vertex> vertices;
  std::vector<HalfEdge> halfedges;
  std::vector<Face> faces;

  std::unordered_map<DCELKey, int, DCELKeyHash> vmap;
  std::unordered_map<std::pair<int, int>, int, PairHash> emap;

  // Return (or create) the vertex index for point p.
  int get_vertex(const Vec2<T> &p) {
    constexpr long long scale = std::is_same<T, double>::value ? 1000000000000LL : 1000000000LL;
    DCELKey key{(long long)std::round(p.x * scale), (long long)std::round(p.y * scale)};
    auto it = vmap.find(key);
    if (it != vmap.end()) return it->second;
    int id = (int)vertices.size();
    vertices.push_back({p});
    vmap[key] = id;
    return id;
  }

  // Return the halfedge index for directed edge a→b, creating the pair if needed.
  int get_halfedge(int a, int b) {
    if (a == b) return -1;
    int lo = std::min(a, b), hi = std::max(a, b);
    auto it = emap.find({lo, hi});
    if (it != emap.end()) {
      int id = it->second; // id = lo→hi
      return (a == lo) ? id : id + 1;
    }
    int id = (int)halfedges.size();
    halfedges.push_back({}); // id   : lo→hi
    halfedges.push_back({}); // id+1 : hi→lo
    halfedges[id].origin     = lo;
    halfedges[id].twin       = id + 1;
    halfedges[id + 1].origin = hi;
    halfedges[id + 1].twin   = id;
    emap[{lo, hi}]           = id;
    return (a == lo) ? id : id + 1;
  }
};

// ============================================================================
// build_dcel — populate a DCEL from LaguerreBuilder cells
// ============================================================================


template <typename T>
inline void build_dcel(const LaguerreBuilder<T> & b, DCEL<T> &dcel) {
  for (size_t i = 0; i < b.cells.size(); ++i) {
    const auto &poly = b.cells[i];
    if (poly.size() < 3) continue;

    std::vector<int> vids;
    for (const auto &p : poly) vids.push_back(dcel.get_vertex(p));

    int m     = (int)vids.size();
    int first = -1, prev = -1;

    for (int k = 0; k < m; ++k) {
      int a = vids[k], b = vids[(k + 1) % m];
      if (a == b) continue;
      int e = dcel.get_halfedge(a, b);
      if (e == -1) continue;
      if (first == -1) first = e;
      if (prev != -1) dcel.halfedges[prev].next = e;
      dcel.halfedges[e].face = (int)i;
      prev                   = e;
    }

    if (first != -1) {
      dcel.halfedges[prev].next = first;
      typename DCEL<T>::Face f;
      f.site = (i < b.n_real) ? (int)i : -1;
      f.edge = first;
      dcel.faces.push_back(f);
    }
  }
}


// ============================================================================
// EdgeInfo + for_each_edge
// ============================================================================

template <typename T> struct EdgeInfo {
  Vec2<T> A, B; // endpoints
  Vec2<T> mid;  // midpoint
  T length;     // edge length
  int he;       // canonical halfedge index (always even)
  int site0;    // site on the left of he    (-1 if boundary)
  int site1;    // site on the left of twin  (-1 if boundary)
  bool border;  // true if one side has no face

  std::pair<int, int> sorted_sites() const {
    return {std::min(site0, site1), std::max(site0, site1)};
  }
};

// Iterate over each unique geometric edge once, calling fn(EdgeInfo).
template <typename T, typename Fn> void for_each_edge(const DCEL<T> &dcel, Fn &&fn) {
  for (size_t i = 0; i < dcel.halfedges.size(); i += 2) {
    const auto &e = dcel.halfedges[i];
    const auto &t = dcel.halfedges[i + 1];
    if (e.origin < 0 || t.origin < 0) continue;

    const auto &A = dcel.vertices[e.origin].p;
    const auto &B = dcel.vertices[t.origin].p;
    T dx = B.x - A.x, dy = B.y - A.y;

    EdgeInfo<T> info;
    info.he     = (int)i;
    info.A      = A;
    info.B      = B;
    info.mid    = {(A.x + B.x) * T(0.5), (A.y + B.y) * T(0.5)};
    info.length = std::sqrt(dx * dx + dy * dy);
    info.site0  = (e.face >= 0) ? dcel.faces[e.face].site : -1;
    info.site1  = (t.face >= 0) ? dcel.faces[t.face].site : -1;
    info.border = (info.site0 < 0 || info.site1 < 0);

    fn(info);
  }
}

// ============================================================================
// SVGOptions
// ============================================================================

struct SVGOptions {
  int width{800}, height{800}, margin{30};
  bool draw_cells{true};
  bool draw_sites{true};
  bool draw_site_labels{true};
  bool draw_radii{true};
  bool draw_limits{true};
  bool draw_dcel_labels{true};
};

// ============================================================================
// SVGContext — internal helper shared by both SVG writers
// ============================================================================

template <typename T> struct SVGContext {

  T xmin, ymin, xmax, ymax, scale;
  const SVGOptions &opt;
  std::ofstream f;

  static constexpr const char *palette[] = {"#AEC6E8", "#FFD9A0", "#B5EAD7", "#FFDDE1", "#D5AAFF",
                                            "#FFF4A3", "#C7E9B0", "#FFB7B2", "#B5D8EB", "#E2C9F5"};

  SVGContext(const SVGOptions &opt) : opt(opt) {}

  T tx(T x) const {
    return opt.margin + (x - xmin) * scale;
  }
  T ty(T y) const {
    return opt.margin + (ymax - y) * scale;
  }

  bool open(const std::string &filename, const LaguerreBuilder<T> &b) {
    f.open(filename);
    if (!f) return false;

    xmin = T(1e30);
    ymin = T(1e30);
    xmax = T(-1e30);
    ymax = T(-1e30);

    auto expand = [&](T x, T y) {
      xmin = std::min(xmin, x);
      xmax = std::max(xmax, x);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    };
    for (const auto &s : b.sites) expand(s.x, s.y);
    for (const auto &c : b.cells)
      for (const auto &p : c) expand(p.x, p.y);

    T w = xmax - xmin, h = ymax - ymin;
    xmin -= T(0.05) * w;
    xmax += T(0.05) * w;
    ymin -= T(0.05) * h;
    ymax += T(0.05) * h;
    scale = std::min((opt.width - 2 * opt.margin) / (xmax - xmin), (opt.height - 2 * opt.margin) / (ymax - ymin));

    f << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
      << "width=\"" << opt.width << "\" height=\"" << opt.height << "\">\n";
    return true;
  }

  void write_cells(const LaguerreBuilder<T> &b) {
    if (!opt.draw_cells) return;
    for (size_t i = 0; i < /*b.cells.size()*/ b.n_real; ++i) {
      const auto &poly = b.cells[i];
      if (poly.empty()) continue;
      f << "<polygon points=\"";
      for (const auto &v : poly) f << tx(v.x) << "," << ty(v.y) << " ";
      f << "\" fill=\"" << palette[i % 10] << "\" opacity=\"0.75\" stroke=\"#333\" stroke-width=\"1\" />\n";
    }
  }

  void write_sites(const LaguerreBuilder<T> &b) {
    if (!opt.draw_sites) return;
    for (size_t i = 0; i < /*b.sites.size()*/ b.n_real; ++i) {
      const auto &s = b.sites[i];
      f << "<circle cx=\"" << tx(s.x) << "\" cy=\"" << ty(s.y) << "\" r=\"4\" fill=\"#000\" />\n";
      if (opt.draw_site_labels) {
        f << "<text x=\"" << tx(s.x) + 5 << "\" y=\"" << ty(s.y) - 5 << "\" font-size=\"10\" fill=\"#000\">" << i
          << "</text>\n";  
      }
    }
  }

  void write_radii(const LaguerreBuilder<T> &b) {
    if (!opt.draw_radii) return;
    //for (const auto &s : b.sites) {
      for (size_t i = 0; i < /*b.sites.size()*/ b.n_real; ++i) {
        const auto &s = b.sites[i];
      f << "<circle cx=\"" << tx(s.x) << "\" cy=\"" << ty(s.y) << "\" r=\"" << s.r * scale
        << "\" fill=\"none\" stroke=\"#666\" stroke-dasharray=\"4 3\" />\n";
    }
  }

  void write_limits(const LaguerreBuilder<T> &b) {
    if (!opt.draw_limits) return;
    f << "<rect x=\"" << tx(b.bbox_xmin_) << "\" y=\"" << ty(b.bbox_ymax_) << "\" width=\""
      << (b.bbox_xmax_ - b.bbox_xmin_) * scale << "\" height=\"" << (b.bbox_ymax_ - b.bbox_ymin_) * scale
      << "\" fill=\"none\" stroke=\"#000\""
      << " stroke-width=\"1.5\" stroke-dasharray=\"6 4\" />\n";
  }

  void close() {
    f << "</svg>\n";
  }
};

template <typename T> constexpr const char *SVGContext<T>::palette[];

// ============================================================================
// write_svg — quick render of Laguerre cells
// ============================================================================

template <typename T>
inline bool write_svg(const LaguerreBuilder<T> &b, const std::string &filename, const SVGOptions &opt = {}) {
  SVGContext<T> ctx(opt);
  if (!ctx.open(filename, b)) return false;
  ctx.write_cells(b);
  ctx.write_sites(b);
  ctx.write_radii(b);
  ctx.write_limits(b);
  ctx.close();
  return true;
}

// ============================================================================
// write_svg_dcel — render with DCEL edge topology and labels
// ============================================================================

template <typename T>
inline bool write_svg_dcel(const DCEL<T> &dcel, const LaguerreBuilder<T> &b, const std::string &filename,
                           const SVGOptions &opt = {}) {
  SVGContext<T> ctx(opt);
  if (!ctx.open(filename, b)) return false;

  ctx.write_cells(b);

  for_each_edge<T>(dcel, [&](const EdgeInfo<T> &info) {
    auto [smin, smax] = info.sorted_sites();
    ctx.f << "<line x1=\"" << ctx.tx(info.A.x) << "\""
          << " y1=\"" << ctx.ty(info.A.y) << "\""
          << " x2=\"" << ctx.tx(info.B.x) << "\""
          << " y2=\"" << ctx.ty(info.B.y) << "\""
          << " stroke=\"" << (info.border ? "#f78de2" : "#8db6f7") << "\""
          << " stroke-width=\"2\" />\n";
    
    if (opt.draw_dcel_labels) {
      ctx.f << "<text x=\"" << ctx.tx(info.mid.x) << "\" y=\"" << ctx.ty(info.mid.y)
            << "\" font-size=\"9\" fill=\"#111\" text-anchor=\"middle\">";
      if (info.border) ctx.f << "B (" << info.length << ")";
      else ctx.f << smin << "-" << smax << " (" << info.length << ")";
      ctx.f << "</text>\n";
    }
    
  });

  ctx.write_sites(b);
  ctx.write_radii(b);
  ctx.write_limits(b);
  ctx.close();
  return true;
}

} // namespace lv