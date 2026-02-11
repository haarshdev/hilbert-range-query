#include <chrono>
// tests/serf_driver.cc  5D+pruning
// Lawder Hilbert DB range query driver implementing professor Step 1–3
// with an efficient "sphere -> sub-quadrants" mapping.
//
// Core idea for performance in 5D:
//   Use geometric pruning (minDist/maxDist) AND occupancy pruning
//   (skip boxes that contain no dataset points) to avoid exploring
//   billions of empty sub-quadrants.
//
// Flags:
//   --qnode <name> --rtt <ms> --horder <order>
//   [--rebuild] [--step 1|2|3] [--debug]
//   [--vec_already_ms] (do not multiply Vec by 1000)
//   [--json <path>]
//   [--fp_counts_json <path>]
//   [--out_summary_csv <path>]   (append Table A summary; default tableA_summary.csv)
//   [--out_nodes_csv <path>]     (append Table B nodes; default tableB_nodes.csv)
//   [--out_boxes_csv <path>]     (optional; append Table C box-level debug)
//
// Notes:
//   Hilbert order k implies 2^k cells per axis.
//   Therefore sub-quadrant (cell) size is latency_max / (2^ORDER).
//
// IMPORTANT ALIGNMENT WITH PROF EMAIL:
//   Step 2 is a SHIFTED coordinate system:
//     - (0,0,...) is the leftmost/bottommost (per-dimension minimum) node
//     - units are latency units (ms) and identical across dimensions
//   We do NOT apply any extra global-span scaling/normalization.
//
// UPDATES (requested):
//   (1) Professor experiment metrics: TP/FP/FN computed vs RTT truth set directly
//       (no intersection with coordinate-defined region).
//   (3) Repeat-offender FP tracking across multiple runs via --fp_counts_json.
//   (4) Append two CSVs per run: Table A (summary) + Table B (node lists).
//
// NEW UPDATE (this request):
//   (5) Box-level RTT lower-bound pruning BEFORE calling Lawder:
//       Keep existing prunes (Vec-sphere + occupancy) and add a third prune:
//         "Even in the best possible case, could any node in this box have RTT <= T?"
//       Best possible case means:
//         - Use the minimum possible Vec distance from query to this box (minDist).
//         - Use the minimum Height among nodes inside the box.
//         - Use the minimum (most negative) Adjustment among nodes inside the box.
//       If that optimistic lower bound already exceeds T_ms, the entire box cannot
//       contain any true answers -> prune it (reduces FPs without post filtering).
//
// Compatibility fix:
//   Accept JSON root as either an array (old) or an object containing "nodes" array (new).
#ifdef DEV
  #include "../db/db.h"
  #include "../hilbert/hilbert.h"
#else
  #include "db.h"
#endif
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <limits>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <sys/stat.h>
#include <cstdint>
#include <errno.h>
#include "../tests/third-party/json.hpp"
using json = nlohmann::json;
using namespace std;

// Box hit record used by optional Table C output.
struct BoxHit {
    std::string box_key;
    std::string node_name; // full name
};
// ---------------- helpers ----------------
static bool read_all(const string& path, string& out) {
    ifstream f(path.c_str(), ios::in | ios::binary);
    if (!f) return false;
    ostringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    return true;
}
static bool file_exists(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}
static bool ensure_dir_exists(const std::string& dir) {
    int rc = mkdir(dir.c_str(), 0755);
    if (rc == 0) return true;          // created successfully
    return errno == EEXIST;            // already exists
}
static bool read_json_file(const std::string& path, json& out) {
    std::string text;
    if (!read_all(path, text)) return false;
    try { out = json::parse(text); }
    catch (...) { return false; }
    return true;
}
static bool extract_nodes_array(const json& root, json& out_nodes) {
    if (root.is_array()) { out_nodes = root; return true; }
    if (root.is_object() && root.contains("nodes") && root["nodes"].is_array()) {
        out_nodes = root["nodes"]; return true;
    }
    return false;
}
static void usage(const char* prog) {
    std::cerr
      << "Usage:\n"
      << "  " << prog << " --qnode <name> --rtt <ms> --horder <order>\n"
      << "         [--rebuild] [--step 1|2|3] [--debug]\n"
      << "         [--vec_already_ms]\n"
      << "         [--json <path>]\n"
      << "         [--fp_counts_json <path>]\n"
      << "         [--out_summary_csv <path>]   (default tableA_summary.csv)\n"
      << "         [--out_nodes_csv <path>]     (default tableB_nodes.csv)\n"
      << "         [--out_boxes_csv <path>]     (optional; Table C box-level debug)\n"
      << "Example:\n"
      << "  " << prog << " --json cluster-status-15012026.json --qnode clab-nebula-serf1 --rtt 15 --horder 10 \\\n"
      << "      --out_summary_csv tableA_summary.csv --out_nodes_csv tableB_nodes.csv\n";
}
static string point_key(const array<PU_int,5>& p) {
    return to_string(p[0]) + "," + to_string(p[1]) + "," + to_string(p[2]) + "," +
           to_string(p[3]) + "," + to_string(p[4]);
}
static PU_int clamp_pu(long long v, PU_int lo, PU_int hi) {
    if (v < (long long)lo) return lo;
    if (v > (long long)hi) return hi;
    return (PU_int)v;
}
// Short name helper for CSV output (e.g., "clab-nebula-serf1" -> "serf1")
static std::string short_name(const std::string& full) {
    size_t pos = full.rfind('-');
    if (pos == std::string::npos) return full;
    if (pos + 1 >= full.size()) return full;
    return full.substr(pos + 1);
}
// CSV helpers
static bool file_empty_or_missing(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return true;
    return (st.st_size == 0);
}
static std::string csv_escape(const std::string& s) {
    // Escape for CSV: wrap in quotes if contains comma/quote/newline, double quotes inside.
    bool need = false;
    for (char c : s) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') { need = true; break; }
    }
    if (!need) return s;
    std::string out;
    out.reserve(s.size() + 8);
    out.push_back('"');
    for (char c : s) {
        if (c == '"') out.push_back('"');
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}
static std::string join_names_sorted(std::vector<std::string> v) {
    std::sort(v.begin(), v.end());
    std::string out;
    for (size_t i = 0; i < v.size(); i++) {
        if (i) out += " ";
        out += v[i];
    }
    return out;
}
static std::vector<std::string> to_short_names(const std::vector<std::string>& in) {
    std::vector<std::string> out;
    out.reserve(in.size());
    for (const auto& s : in) out.push_back(short_name(s));
    return out;
}
// ---------------- exact integer geometry ----------------
static inline bool isSingleCell(const array<long long,5>& lo,
                                const array<long long,5>& hi)
{
    for (int d=0; d<5; d++) if (lo[d] != hi[d]) return false;
    return true;
}
static inline long long dist2_PointToPoint(
    const array<long long,5>& a,
    const array<long long,5>& b)
{
    long long acc = 0;
    for (int d=0; d<5; d++) {
        long long t = a[d] - b[d];
        acc += t*t;
    }
    return acc;
}
// ---------------- ms-space geometry helpers ----------------
static inline double minDist2_PointToBox_ms(
    const array<double,5>& c,
    const array<double,5>& lo,
    const array<double,5>& hi)
{
    double acc = 0.0;
    for (int d=0; d<5; d++) {
        double x = c[d];
        if (x < lo[d]) { double t = lo[d] - x; acc += t*t; }
        else if (x > hi[d]) { double t = x - hi[d]; acc += t*t; }
    }
    return acc;
}
static inline double maxDist2_PointToBox_ms(
    const array<double,5>& c,
    const array<double,5>& lo,
    const array<double,5>& hi)
{
    double acc = 0.0;
    for (int d=0; d<5; d++) {
        double x = c[d];
        double a = fabs(x - lo[d]);
        double b = fabs(x - hi[d]);
        double m = std::max(a,b);
        acc += m*m;
    }
    return acc;
}
// ---------------- Lawder range query helper ----------------
struct QueryStats {
    uint64_t open_ok = 0;
    uint64_t open_fail = 0;
    uint64_t fetched_rows = 0;
    uint64_t matched_names = 0;

    // Lawder 'cell' instrumentation (accumulated across all submitted boxes)
    uint64_t lawder_box_queries = 0;
};
// open a range set and fetch at most 1 record.
// Used ONLY for occupancy pruning.
static bool box_has_any_point(DBASE* DB, const array<long long,5>& lo, const array<long long,5>& hi) {
    PU_int LB[5], UB[5], result[5];
    for (int d=0; d<5; d++) { LB[d]=(PU_int)lo[d]; UB[d]=(PU_int)hi[d]; }
    int set_id = -1;
    if (true != DB->db_range_open_set(LB, UB, &set_id)) return false;
    bool any = (true == DB->db_range_fetch_another(set_id, result));
    DB->db_close_set(set_id);
    return any;
}
static void run_one_box_query(
    DBASE* DB,
    const array<long long,5>& lo,
    const array<long long,5>& hi,
    const unordered_map<string, vector<string>>& point_to_names,
    const string& qname,
    const std::string& box_key,
    std::vector<BoxHit>* hits,
    unordered_map<string,int>& hilbert_set,
    vector<string>& hilbert_names,
    QueryStats& qs)
{
    PU_int LB[5], UB[5], result[5];
    for (int d=0; d<5; d++) { LB[d]=(PU_int)lo[d]; UB[d]=(PU_int)hi[d]; }
    int set_id = -1;
    qs.lawder_box_queries++;
    if (true == DB->db_range_open_set(LB, UB, &set_id)) {
        qs.open_ok++;
        while (true == DB->db_range_fetch_another(set_id, result)) {
            qs.fetched_rows++;
            array<PU_int,5> rp = {result[0],result[1],result[2],result[3],result[4]};
            string k = point_key(rp);
            auto itp = point_to_names.find(k);
            if (itp != point_to_names.end()) {
                for (const string& nm : itp->second) {
                    if (nm == qname) continue;
                    if (hilbert_set.emplace(nm,1).second) {
                        hilbert_names.push_back(nm);
                    }
                    if (hits) { hits->push_back(BoxHit{box_key, nm}); }
                    qs.matched_names++;
                }
            }
        }
        
        DB->db_close_set(set_id);
    } else {
        qs.open_fail++;
    }
}
// ---------------- Sphere cover recursion (exact + occupancy pruning) ----------------
struct CoverStats {
    uint64_t visited = 0;
    uint64_t pruned_outside = 0;
    uint64_t pruned_empty = 0;
    // NEW: pruned by RTT optimistic lower bound (before querying Lawder)
    uint64_t pruned_rtt_lb = 0;
    uint64_t accepted_inside = 0;
    uint64_t accepted_leaf = 0;

    // NEW (Table A only): how many dataset points fall inside boxes of each category
    uint64_t pruned_rtt_lb_points = 0;
    uint64_t accepted_leaf_points = 0;
    uint64_t accepted_inside_points = 0;
};
// ---------------- Box-level debug (Table C) ----------------
struct BoxDebug {
    array<long long,5> lo;
    array<long long,5> hi;
    int depth = 0;
    uint64_t point_count = 0;

    // Lower-bound components
    double min_vec_dist = std::numeric_limits<double>::quiet_NaN();
    double q_height_ms  = std::numeric_limits<double>::quiet_NaN();
    double q_adj_ms     = std::numeric_limits<double>::quiet_NaN();
    double box_min_h    = std::numeric_limits<double>::quiet_NaN();
    double box_min_adj  = std::numeric_limits<double>::quiet_NaN();

    // Which nodes contributed the minima (may be different nodes)
    std::string min_h_node;
    std::string min_adj_node;

    // Derived bounds
    double base_lb_ms          = std::numeric_limits<double>::quiet_NaN(); // min_vec_dist + q_height + box_min_h
    double adjusted_lb_ms      = std::numeric_limits<double>::quiet_NaN(); // base_lb + q_adj + box_min_adj  (current pruning bound)
    int serf_apply_adj_bestcase = 0; // 1 if adjusted_lb_ms > 0, else 0 (Serf would ignore adjustment)
    double serf_effective_lb_ms = std::numeric_limits<double>::quiet_NaN(); // adjusted if >0 else base
    int serf_effective_lb_gt_T  = 0; // 1 if serf_effective_lb_ms > T_ms
};

static inline std::string box_key_str(const array<long long,5>& lo, const array<long long,5>& hi) {
    // Stable key for this process (single query run)
    std::string k;
    k.reserve(128);
    for (int d=0; d<5; d++) { if (d) k.push_back(','); k += std::to_string(lo[d]); }
    k.push_back('|');
    for (int d=0; d<5; d++) { if (d) k.push_back(','); k += std::to_string(hi[d]); }
    return k;
}

// NEW helper:
// Find optimistic (best-case) Height/Adjustment within this box by scanning the nodes.
// This is intentionally simple and safe because N is small (tens/hundreds), and it avoids
// building complex per-box summary structures.
static inline void box_min_height_adj(
    const vector<array<PU_int,5>>& pts,
    const vector<double>& heights_ms,
    const vector<double>& adjs_ms,
    const vector<string>& names,
    const array<long long,5>& lo,
    const array<long long,5>& hi,
    double& out_min_h,
    double& out_min_adj,
    size_t& out_min_h_idx,
    size_t& out_min_adj_idx,
    uint64_t& out_point_count)
{
    out_min_h = std::numeric_limits<double>::infinity();
    out_min_adj = std::numeric_limits<double>::infinity(); // "most negative" => minimum value
    out_min_h_idx = (size_t)-1;
    out_min_adj_idx = (size_t)-1;
    out_point_count = 0;
    for (size_t i = 0; i < pts.size(); i++) {
        bool inside = true;
        for (int d = 0; d < 5; d++) {
            long long v = (long long)pts[i][d];
            if (v < lo[d] || v > hi[d]) { inside = false; break; }
        }
        if (!inside) continue;
        out_point_count++;
        if (heights_ms[i] < out_min_h) { out_min_h = heights_ms[i]; out_min_h_idx = i; }
        if (adjs_ms[i] < out_min_adj)  { out_min_adj = adjs_ms[i];  out_min_adj_idx = i; }
    }
}
// depth counts subdivision steps; max depth is ORDER (leaf cell size = 1).
static void coverSphere5D_indexed(
    DBASE* DB,
    int ORDER,
    const array<double,5>& center_ms,   // center in Step-2 shifted ms-space
    int T_ms,                           // radius in ms
    double cell_size_ms,                // for mapping grid box -> ms box
    const unordered_map<string, vector<string>>& point_to_names,
    const string& qname,
    const vector<array<PU_int,5>>& pts, // all node grid coords for lower-bound scan
    const vector<double>& heights_ms,   // node Height in ms
    const vector<double>& adjs_ms,      // node Adjustment in ms
    const vector<string>& names,        // full node names (index-aligned)
    double q_height_ms,                 // query Height in ms
    double q_adj_ms,                    // query Adjustment in ms
    unordered_map<string,int>& hilbert_set,
    vector<string>& hilbert_names,
    const array<long long,5>& lo,
    const array<long long,5>& hi,
    int depth,
    CoverStats& st,
    QueryStats& qs,
    unordered_map<string, BoxDebug>* box_debug_map,
    vector<BoxHit>* box_hits)
{
    st.visited++;
    const double r2_ms = (double)T_ms * (double)T_ms;
    // Map grid box [lo..hi] to Step-2 ms-space box [lo_ms..hi_ms]
    array<double,5> lo_ms, hi_ms;
    for (int d=0; d<5; d++) {
        lo_ms[d] = (double)lo[d] * cell_size_ms;
        hi_ms[d] = (double)(hi[d] + 1) * cell_size_ms;
    }
    // Geometric prune: fully outside ms-space Vec-sphere
    const double minDist2 = minDist2_PointToBox_ms(center_ms, lo_ms, hi_ms);
    if (minDist2 > r2_ms) {
        st.pruned_outside++;
        return;
    }
    // Occupancy prune (skip boxes that contain no dataset points)
    if (!box_has_any_point(DB, lo, hi)) {
        st.pruned_empty++;
        return;
    }
    // ------------------------------------------------------------
    // NEW PRUNE (box-level RTT optimistic lower bound)
    //
    // Question:
    //   "Even in the BEST possible case, could any node in this box have RTT <= T_ms?"
    //
    // Best possible case is intentionally optimistic:
    //   - Use the closest possible Vec distance from query to this box (sqrt(minDist2)).
    //   - Inside this box, pick the node with the smallest Height.
    //   - Inside this box, pick the node with the most negative Adjustment.
    //
    // If even this optimistic lower bound is > T_ms, then NO node in this box can satisfy
    // RTT <= T_ms, so the whole box is guaranteed FP and can be pruned BEFORE calling Lawder.
    // ------------------------------------------------------------
    double box_min_h = 0.0, box_min_adj = 0.0;
    size_t box_min_h_idx = (size_t)-1;
    size_t box_min_adj_idx = (size_t)-1;
    uint64_t box_point_count = 0;
    box_min_height_adj(pts, heights_ms, adjs_ms, names, lo, hi,
                      box_min_h, box_min_adj, box_min_h_idx, box_min_adj_idx, box_point_count);
    // box_has_any_point already says there is at least one point, so these should be finite.
    // Still, guard defensively.
    if (std::isfinite(box_min_h) && std::isfinite(box_min_adj)) {
        const double min_vec_dist = std::sqrt(minDist2);
        // RTT lower bound for any node in this box, using Serf coordinate model terms.
        // This is a LOWER bound: real RTTs for nodes in this box are >= this value.
        const double base_lb = min_vec_dist + q_height_ms + box_min_h;
        const double best_case_rtt_lb = base_lb + q_adj_ms + box_min_adj;
        const int serf_apply_adj = (best_case_rtt_lb > 0.0) ? 1 : 0;
        const double serf_effective_lb = serf_apply_adj ? best_case_rtt_lb : base_lb;
        if (best_case_rtt_lb > (double)T_ms) {
            st.pruned_rtt_lb++;
            st.pruned_rtt_lb_points += box_point_count;
            return;
        }
    }
    const bool fully_inside = (maxDist2_PointToBox_ms(center_ms, lo_ms, hi_ms) <= r2_ms);
    // Leaf
    if (depth >= ORDER || isSingleCell(lo, hi)) {
        st.accepted_leaf++;
        st.accepted_leaf_points += box_point_count;
        if (fully_inside) {
            st.accepted_inside++;
            st.accepted_inside_points += box_point_count;
        }

// Optional Table C debug: record this accepted leaf box and all nodes fetched from it.
const std::string bkey = box_key_str(lo, hi);
if (box_debug_map) {
    BoxDebug bd;
    bd.lo = lo;
    bd.hi = hi;
    bd.depth = depth;
    bd.point_count = box_point_count;
    bd.min_vec_dist = std::sqrt(minDist2);
    bd.q_height_ms = q_height_ms;
    bd.q_adj_ms = q_adj_ms;
    bd.box_min_h = box_min_h;
    bd.box_min_adj = box_min_adj;
    bd.min_h_node = (box_min_h_idx != (size_t)-1 && box_min_h_idx < names.size()) ? names[box_min_h_idx] : "";
    bd.min_adj_node = (box_min_adj_idx != (size_t)-1 && box_min_adj_idx < names.size()) ? names[box_min_adj_idx] : "";
    bd.base_lb_ms = bd.min_vec_dist + q_height_ms + box_min_h;
    bd.adjusted_lb_ms = bd.base_lb_ms + q_adj_ms + box_min_adj;
    bd.serf_apply_adj_bestcase = (bd.adjusted_lb_ms > 0.0) ? 1 : 0;
    bd.serf_effective_lb_ms = bd.serf_apply_adj_bestcase ? bd.adjusted_lb_ms : bd.base_lb_ms;
    bd.serf_effective_lb_gt_T = (bd.serf_effective_lb_ms > (double)T_ms) ? 1 : 0;

    // Insert only once per box key (same box may be reached once in this traversal).
    box_debug_map->emplace(bkey, std::move(bd));
}
        run_one_box_query(DB, lo, hi, point_to_names, qname,
                          bkey, box_hits,
                          hilbert_set, hilbert_names, qs);
        return;
    }
    // Subdivide into 32 children (5D)
    array<long long,5> mid;
    for (int d=0; d<5; d++) mid[d] = (lo[d] + hi[d]) >> 1;
    for (int mask=0; mask<(1<<5); mask++) {
        array<long long,5> clo, chi;
        for (int d=0; d<5; d++) {
            if (mask & (1<<d)) { clo[d] = mid[d] + 1; chi[d] = hi[d]; }
            else               { clo[d] = lo[d];      chi[d] = mid[d]; }
        }
        bool empty=false;
        for (int d=0; d<5; d++) if (clo[d] > chi[d]) { empty=true; break; }
        if (empty) continue;
        auto __tq0 = std::chrono::high_resolution_clock::now();
    coverSphere5D_indexed(DB, ORDER, center_ms, T_ms, cell_size_ms,
                              point_to_names, qname,
                              pts, heights_ms, adjs_ms, names, q_height_ms, q_adj_ms,
                              hilbert_set, hilbert_names,
                              clo, chi, depth+1, st, qs, box_debug_map, box_hits);
    }
}
int main(int argc, char** argv) {
    std::string qname;
    int T_ms = -1;
    int ORDER = -1;
    bool rebuild = false;
    int step = 3;
    bool debug = false;
    bool vec_already_ms = false;
    std::string input_json = "cluster-status-15012026.json"; //"cluster-status-15012026.json";
    std::string fp_counts_json;
    // NEW: output CSVs (append unsorted; sort later in a separate step)
    std::string out_summary_csv = "summary_csv/tableA_summary.csv";
    std::string out_nodes_csv   = "nodes_csv/tableB_nodes.csv";
    std::string out_boxes_csv; // optional Table C box-level debug
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--qnode" && i + 1 < argc) qname = argv[++i];
        else if (a == "--rtt" && i + 1 < argc) T_ms = std::stoi(argv[++i]);
        else if (a == "--horder" && i + 1 < argc) ORDER = std::stoi(argv[++i]);
        else if (a == "--rebuild") rebuild = true;
        else if (a == "--step" && i + 1 < argc) step = std::stoi(argv[++i]);
        else if (a == "--debug") debug = true;
        else if (a == "--vec_already_ms") vec_already_ms = true;
        else if (a == "--json" && i + 1 < argc) input_json = argv[++i];
        else if (a == "--fp_counts_json" && i + 1 < argc) fp_counts_json = argv[++i];
        else if (a == "--out_summary_csv" && i + 1 < argc) out_summary_csv = argv[++i];
        else if (a == "--out_nodes_csv" && i + 1 < argc) out_nodes_csv = argv[++i];
        else if (a == "--out_boxes_csv" && i + 1 < argc) out_boxes_csv = argv[++i];
        else if (a == "--help" || a == "-h") { usage(argv[0]); return 0; }
        else { std::cerr << "Unknown or incomplete arg: " << a << "\n"; usage(argv[0]); return 2; }
    }
    if (qname.empty() || T_ms < 0 || ORDER < 1 || ORDER > 30 || step < 1 || step > 3) {
        usage(argv[0]);
        return 2;
    }
    cout << "Query node: " << qname << "\n";
    cout << "RTT threshold (ms): " << T_ms << "\n";
    // ---------- load JSON ----------
    json root;
    if (!read_json_file(input_json, root)) {
        cerr << "Cannot read/parse JSON file: " << input_json << "\n";
        return 1;
    }
    json arr;
    if (!extract_nodes_array(root, arr)) {
        cerr << "JSON must be an array, or an object with a 'nodes' array: " << input_json << "\n";
        return 1;
    }
    const size_t N = arr.size();
    if (N == 0) { cerr << "JSON array is empty\n"; return 1; }
    // ---------- build name -> index map ----------
    unordered_map<string, int> idx;
    idx.reserve(N);
    for (size_t i = 0; i < N; i++) if (arr[i].contains("name")) idx[arr[i]["name"].get<string>()] = (int)i;
    auto it = idx.find(qname);
    if (it == idx.end()) { cerr << "Query node not found in JSON: " << qname << "\n"; return 1; }
    int qi = it->second;
    // ---------- RTT map for q (Vivaldi predicted truth set from JSON) ----------
    if (!arr[qi].contains("rtts") || !arr[qi]["rtts"].is_object()) {
        cerr << "JSON node has no rtts object for: " << qname << "\n"; return 1;
    }
    const json& rtts_q = arr[qi]["rtts"];
    // ---------- extract all names + 5D Vec + Height + Adjustment ----------
    vector<string> names(N);
    vector<array<double,5>> vecs_ms(N);
    // Height and Adjustment are treated in the same unit conversion as Vec.
    // If JSON is in seconds, multiply by 1000 to get ms.
    vector<double> heights_ms(N, 0.0);
    vector<double> adjs_ms(N, 0.0);
    const double VEC_TO_MS = vec_already_ms ? 1.0 : 1000.0;
    for (size_t i = 0; i < N; i++) {
        names[i] = arr[i]["name"].get<string>();
        const json& v = arr[i]["coordinate"]["Vec"];
        for (int d = 0; d < 5; d++) vecs_ms[i][d] = v[d].get<double>() * VEC_TO_MS;
        if (arr[i]["coordinate"].contains("Height")) heights_ms[i] = arr[i]["coordinate"]["Height"].get<double>() * VEC_TO_MS;
        if (arr[i]["coordinate"].contains("Adjustment")) adjs_ms[i] = arr[i]["coordinate"]["Adjustment"].get<double>() * VEC_TO_MS;
    }
    const double q_height_ms = heights_ms[qi];
    const double q_adj_ms = adjs_ms[qi];

    // ---------- GT (truth set for this query) ----------
    vector<pair<string, double>> gt;
    for (size_t j = 0; j < N; j++) {
        if ((int)j == qi) continue;
        const string& nname = names[j];
        if (!rtts_q.contains(nname)) continue;
        double rtt = rtts_q[nname].get<double>();
        if (rtt <= (double)T_ms) gt.push_back({nname, rtt});
    }

    cout << "GT count: " << gt.size() << "\n";
    for (auto& p : gt) cout << "  " << p.first << "  " << p.second << "\n";

    // ============================================================
    // Step 1: latency_max and cell_size (based on JSON rtts)
    // ============================================================
    double latency_max = 83.283;

    const PU_int GRID_MAX = (PU_int)((1u << ORDER) - 1u);
    const double denom = (double)(1u << ORDER);
    const double cell_size_ms = latency_max / denom;

    std::string dbname = "serfdb_o" + std::to_string(ORDER);

    cout << "DB name: " << dbname << " (ORDER=" << ORDER << ", GRID_MAX=" << ((1u<<ORDER)-1u) << ")\n";
    cout << "latency_max(ms): " << latency_max << "\n";
    cout << "cell_size(ms): " << cell_size_ms << "  (latency_max / 2^ORDER)\n";

    if (step == 1) {
        cout << "[STEP1 DONE]\n";
        return 0;
    }

    // ============================================================
    // Step 2: shifted coordinate system (NO global-span scaling) + quantization
    // ============================================================
    array<double,5> mn, mx;
    for (int d=0; d<5; d++) {
        mn[d] =  numeric_limits<double>::infinity();
        mx[d] = -numeric_limits<double>::infinity();
    }

    for (size_t i=0; i<N; i++) {
        for (int d=0; d<5; d++) {
            mn[d] = std::min(mn[d], vecs_ms[i][d]);
            mx[d] = std::max(mx[d], vecs_ms[i][d]);
        }
    }

    if (debug) {
        cout << "Shift minima mn[d] (ms): ";
        for (int d=0; d<5; d++) cout << mn[d] << (d==4? "\n":" ");
        cout << "Shift maxima mx[d] (ms): ";
        for (int d=0; d<5; d++) cout << mx[d] << (d==4? "\n":" ");
        cout << "Scaling mode: SHIFT-ONLY (no global-span normalization)\n";
    }

    auto quantize_prof = [&](double x_ms, int d) -> PU_int {
        double shifted = x_ms - mn[d];
        if (shifted < 0.0) shifted = 0.0;
        if (shifted > latency_max) shifted = latency_max;
        long long gi = (long long)floor(shifted / cell_size_ms);
        return clamp_pu(gi, (PU_int)0, GRID_MAX);
    };

    vector<array<PU_int,5>> pts(N);
    for (size_t i=0; i<N; i++) {
        for (int d=0; d<5; d++) pts[i][d] = quantize_prof(vecs_ms[i][d], d);
    }

    cout << "Query point grid coords: ";
    for (int d=0; d<5; d++) cout << pts[qi][d] << (d==4? "\n":" ");

    if (step == 2) {
        cout << "[STEP2 DONE]\n";
        return 0;
    }

    // ============================================================
    // Step 3: build DB and do sphere->sub-quadrants->Lawder
    // ============================================================
    unordered_map<string, vector<string>> point_to_names;
    point_to_names.reserve(N*2);
    for (size_t i=0; i<N; i++) point_to_names[point_key(pts[i])].push_back(names[i]);

    const int DIMS = 5;
    const int BT_NODE_ENTRIES = 10;
    const int BUFFER_PAGES = 10;
    const int PAGE_RECORDS = 200;

    bool db_exists =
        file_exists(dbname + ".db") &&
        file_exists(dbname + ".idx") &&
        file_exists(dbname + ".inf");

    if (rebuild || !db_exists) {
        remove((dbname + ".db").c_str());
        remove((dbname + ".idx").c_str());
        remove((dbname + ".inf").c_str());
        remove((dbname + ".fpl").c_str());
    }

    DBASE* DB = new DBASE(dbname, DIMS, BT_NODE_ENTRIES, BUFFER_PAGES, PAGE_RECORDS);

    bool did_build = (rebuild || !db_exists);
    unsigned long long encode_calls_total = 0;
    double t_encode_total_ms = std::numeric_limits<double>::quiet_NaN();
    double t_index_build_ms = std::numeric_limits<double>::quiet_NaN();

    std::chrono::high_resolution_clock::time_point __tb0;
    if (did_build) {
        hilbert_reset_encode_stats();
        __tb0 = std::chrono::high_resolution_clock::now();
    }

    if (did_build) {
        cout << "Creating new DB files...\n";
        if (!DB->db_create()) { cerr << "DB create failed\n"; delete DB; return 1; }
    }

    if (!DB->db_open()) { cerr << "DB open failed\n"; delete DB; return 1; }

    if (debug) {
        // Count total records stored in the Lawder DB by scanning full key space.
        // This is a sanity check to confirm DB cardinality.
        PU_int lb[5] = {0,0,0,0,0};
        PU_int ub[5] = { (PU_int)GRID_MAX, (PU_int)GRID_MAX, (PU_int)GRID_MAX, (PU_int)GRID_MAX, (PU_int)GRID_MAX };

        int set_id = -1;
        if (DB->db_range_open_set(lb, ub, &set_id)) {
            unsigned long long db_rows = 0;
            PU_int rec[5];
            while (DB->db_range_fetch_another(set_id, rec)) {
                db_rows++;
            }
            DB->db_close_set(set_id);
            cerr << "[debug] DB total records (full scan): " << db_rows << "\n";
        } else {
            cerr << "[debug] Could not open full-scan range set\n";
        }
    }


    if (did_build) {
        int inserted = 0;
        for (size_t i=0; i<N; i++) {
            PU_int p[5]; for (int d=0; d<5; d++) p[d] = pts[i][d];
            if (!DB->db_data_insert(p)) {
                cerr << "Insert failed for " << names[i] << "\n";
                DB->db_close(); delete DB; return 1;
            }
            inserted++;
        }
        cout << "Inserted " << inserted << " points into DB\n";
        auto __tb1 = std::chrono::high_resolution_clock::now();
        t_index_build_ms = std::chrono::duration<double, std::milli>(__tb1 - __tb0).count();
        encode_calls_total = hilbert_get_encode_calls();
        t_encode_total_ms = (double)hilbert_get_encode_time_ns() / 1000000.0;
    }

    const double r_cont = (double)T_ms / cell_size_ms;
    const long long r_cells = (long long)ceil(r_cont);

    cout << "radius_in_cells(continuous): " << r_cont << "  (T_ms / cell_size_ms)\n";
    cout << "radius_in_cells(integer): " << r_cells << "\n";

    array<long long,5> center = {
        (long long)pts[qi][0], (long long)pts[qi][1], (long long)pts[qi][2],
        (long long)pts[qi][3], (long long)pts[qi][4]
    };

    array<long long,5> root_lo, root_hi;
    for (int d=0; d<5; d++) {
        root_lo[d] = std::max(0LL, center[d] - r_cells);
        root_hi[d] = std::min((long long)GRID_MAX, center[d] + r_cells);
    }

    cout << "Root AABB around center (clamped):\n";
    cout << "  LB: " << root_lo[0] << " " << root_lo[1] << " " << root_lo[2] << " " << root_lo[3] << " " << root_lo[4] << "\n";
    cout << "  UB: " << root_hi[0] << " " << root_hi[1] << " " << root_hi[2] << " " << root_hi[3] << " " << root_hi[4] << "\n";

    array<double,5> center_ms;
    for (int d=0; d<5; d++) {
        double s = vecs_ms[qi][d] - mn[d];
        if (s < 0.0) s = 0.0;
        if (s > latency_max) s = latency_max;
        center_ms[d] = s;
    }

    vector<string> hilbert_names;
    unordered_map<string,int> hilbert_set;
    hilbert_set.reserve(N*4);

    QueryStats qs_cover;
    CoverStats st;
    double t_hilbert_query_ms = std::numeric_limits<double>::quiet_NaN();

    // Optional Table C (box-level debug): only populated when --out_boxes_csv is provided.
    unordered_map<string, BoxDebug> box_debug_map;
    vector<BoxHit> box_hits;

    auto __tq0 = std::chrono::high_resolution_clock::now();
    coverSphere5D_indexed(DB, ORDER, center_ms, T_ms, cell_size_ms,
                          point_to_names, qname,
                          pts, heights_ms, adjs_ms, names, q_height_ms, q_adj_ms,
                          hilbert_set, hilbert_names,
                          root_lo, root_hi, 0, st, qs_cover,
                          (out_boxes_csv.empty() ? nullptr : &box_debug_map),
                          (out_boxes_csv.empty() ? nullptr : &box_hits));
    auto __tq1 = std::chrono::high_resolution_clock::now();
    t_hilbert_query_ms = std::chrono::duration<double, std::milli>(__tq1 - __tq0).count();

    cout << "Cover stats:\n";
    cout << "  visited=" << st.visited
         << " pruned_outside=" << st.pruned_outside
         << " pruned_empty=" << st.pruned_empty
         << " pruned_rtt_lb=" << st.pruned_rtt_lb
         << " accepted_inside=" << st.accepted_inside
         << " accepted_leaf=" << st.accepted_leaf << "\n";

    cout << "Hilbert result count (excluding self): " << hilbert_names.size() << "\n";
    for (const string& nm : hilbert_names) {
        double rtt = rtts_q.contains(nm) ? rtts_q[nm].get<double>() : -1.0;
        cout << "  " << nm << "  rtt=" << rtt << "\n";
    }

    // ============================================================
    // Professor experiment metrics vs RTT truth set:
    //   Truth set = nodes with RTT(q, node) <= T_ms (from JSON rtts).
    // ============================================================
    unordered_map<string,int> gt_set;
    gt_set.reserve(gt.size()*2 + 10);
    for (auto& p : gt) gt_set.emplace(p.first,1);

    vector<string> fp_nodes;
    fp_nodes.reserve(hilbert_names.size());

    int TP = 0, FP = 0, FN = 0;
    for (const auto& kv : hilbert_set) {
        if (gt_set.find(kv.first) != gt_set.end()) TP++;
        else { FP++; fp_nodes.push_back(kv.first); }
    }
    for (const auto& kv : gt_set) {
        if (hilbert_set.find(kv.first) == hilbert_set.end()) FN++;
    }

    int truth_count     = (int)gt_set.size();
    int retrieved_count = TP + FP;

    double precision = std::numeric_limits<double>::quiet_NaN();
    double recall    = std::numeric_limits<double>::quiet_NaN();
    double jaccard   = std::numeric_limits<double>::quiet_NaN();

    if (retrieved_count > 0) precision = (double)TP / (double)retrieved_count;
    if (truth_count > 0)     recall    = (double)TP / (double)truth_count;
    if (retrieved_count + FN > 0) jaccard = (double)TP / (double)(retrieved_count + FN);

    cout << "Metrics vs RTT truth set:\n";
    cout << "  TP=" << TP << " FP=" << FP << " FN=" << FN << "\n";
    cout << "  precision=" << precision << " recall=" << recall << " jaccard=" << jaccard << "\n";

    if (!fp_nodes.empty()) {
        sort(fp_nodes.begin(), fp_nodes.end());
        cout << "FP nodes (" << fp_nodes.size() << "):";
        for (const auto& n : fp_nodes) cout << " " << n;
        cout << "\n";
    }

    // NEW: Build node lists for Table B (details)
    vector<string> gt_nodes;
    gt_nodes.reserve(gt_set.size());
    for (const auto& kv : gt_set) gt_nodes.push_back(kv.first);

    vector<string> returned_nodes = hilbert_names; // already excludes self
    // Build TP nodes list
    vector<string> tp_nodes;
    tp_nodes.reserve((size_t)TP);
    for (const auto& kv : hilbert_set) {
        if (gt_set.find(kv.first) != gt_set.end()) tp_nodes.push_back(kv.first);
    }
    // Build FN nodes list
    vector<string> fn_nodes;
    fn_nodes.reserve((size_t)FN);
    for (const auto& kv : gt_set) {
        if (hilbert_set.find(kv.first) == hilbert_set.end()) fn_nodes.push_back(kv.first);
    }

    // Accumulate repeat-offender FP counts across multiple runs.
    if (!fp_counts_json.empty()) {
        json counts;
        if (file_exists(fp_counts_json)) {
            if (!read_json_file(fp_counts_json, counts) || !counts.is_object()) counts = json::object();
        } else {
            counts = json::object();
        }

        for (const auto& n : fp_nodes) {
            long long v = 0;
            if (counts.contains(n) && counts[n].is_number()) v = counts[n].get<long long>();
            counts[n] = v + 1;
        }

        ofstream out(fp_counts_json.c_str(), ios::out | ios::binary);
        if (out) {
            out << counts.dump(2);
            out.close();
            cout << "Updated FP counts file: " << fp_counts_json << "\n";
        } else {
            cerr << "Warning: could not write fp counts file: " << fp_counts_json << "\n";
        }
    }

    // NEW: Append Table A summary CSV (updated: appended point-count fields)
    {
        const bool need_header = file_empty_or_missing(out_summary_csv);
        ofstream out(out_summary_csv.c_str(), ios::out | ios::app);
        if (!out) {
            cerr << "Warning: could not append to summary CSV: " << out_summary_csv << "\n";
        } else {
            if (need_header) {
                out << "query_node,rtt_ms,hilbert_order,truth_size,returned,tp,fp,fn,precision,recall,jaccard,latency_max_ms,cell_size_ms,mn0_ms,mn1_ms,mn2_ms,mn3_ms,mn4_ms,mx0_ms,mx1_ms,mx2_ms,mx3_ms,mx4_ms,q_height_ms,q_adj_ms,radius_cells_cont,radius_cells_int,root_lb0,root_lb1,root_lb2,root_lb3,root_lb4,root_ub0,root_ub1,root_ub2,root_ub3,root_ub4,cover_visited,cover_pruned_outside,cover_pruned_empty,cover_pruned_rtt_lb,cover_accepted_inside,cover_accepted_leaf,hilbert_count_excl_self,cover_pruned_rtt_lb_points,cover_accepted_leaf_points,cover_accepted_inside_points,lawder_box_queries,t_hilbert_query_ms,did_build,encode_calls_total,t_encode_total_ms,t_index_build_ms\n";
            }

            out << csv_escape(short_name(qname)) << ","
                << T_ms << ","
                << ORDER << ","
                << truth_count << ","
                << retrieved_count << ","
                << TP << ","
                << FP << ","
                << FN << ","
                << precision << ","
                << recall << ","
                << jaccard << ","
                << latency_max << ","
                << cell_size_ms << ","
                << mn[0] << "," << mn[1] << "," << mn[2] << "," << mn[3] << "," << mn[4] << ","
                << mx[0] << "," << mx[1] << "," << mx[2] << "," << mx[3] << "," << mx[4] << ","
                << q_height_ms << "," << q_adj_ms << ","
                << r_cont << "," << r_cells << ","
                << root_lo[0] << "," << root_lo[1] << "," << root_lo[2] << "," << root_lo[3] << "," << root_lo[4] << ","
                << root_hi[0] << "," << root_hi[1] << "," << root_hi[2] << "," << root_hi[3] << "," << root_hi[4] << ","
                << st.visited << "," << st.pruned_outside << "," << st.pruned_empty << "," << st.pruned_rtt_lb << "," << st.accepted_inside << "," << st.accepted_leaf << ","
                << hilbert_names.size() << "," << st.pruned_rtt_lb_points << "," << st.accepted_leaf_points << "," << st.accepted_inside_points << ","
                << qs_cover.lawder_box_queries << "," << t_hilbert_query_ms << "," << (did_build ? 1 : 0) << "," << encode_calls_total << "," << t_encode_total_ms << "," << t_index_build_ms
                << "\n";
        }
    }

    // NEW: Append Table B nodes CSV (unchanged)
    {
        const bool need_header = file_empty_or_missing(out_nodes_csv);
        ofstream out(out_nodes_csv.c_str(), ios::out | ios::app);
        if (!out) {
            cerr << "Warning: could not append to nodes CSV: " << out_nodes_csv << "\n";
        } else {
            if (need_header) {
                out << "query_node,rtt_ms,hilbert_order,gt_nodes,returned_nodes,tp_nodes,fp_nodes,fn_nodes\n";
            }

            vector<string> gt_nodes_s       = to_short_names(gt_nodes);
            vector<string> returned_nodes_s = to_short_names(returned_nodes);
            vector<string> tp_nodes_s       = to_short_names(tp_nodes);
            vector<string> fp_nodes_s       = to_short_names(fp_nodes);
            vector<string> fn_nodes_s       = to_short_names(fn_nodes);

            out << csv_escape(short_name(qname)) << ","
                << T_ms << ","
                << ORDER << ","
                << csv_escape(join_names_sorted(gt_nodes_s)) << ","
                << csv_escape(join_names_sorted(returned_nodes_s)) << ","
                << csv_escape(join_names_sorted(tp_nodes_s)) << ","
                << csv_escape(join_names_sorted(fp_nodes_s)) << ","
                << csv_escape(join_names_sorted(fn_nodes_s))
                << "\n";
        }
    }


// Optional Table C: box-level debug rows for FPs.
// One row per (accepted leaf box, FP node) occurrence.
if (!out_boxes_csv.empty()) {
    const bool need_header = file_empty_or_missing(out_boxes_csv);
    ofstream out(out_boxes_csv.c_str(), ios::out | ios::app);
    if (!out) {
        cerr << "Warning: could not append to boxes CSV: " << out_boxes_csv << "\n";
    } else {
        if (need_header) {
            out << "query_node,rtt_ms,hilbert_order,box_lo,box_hi,depth,box_point_count,"
                   "min_vec_dist_ms,q_height_ms,q_adj_ms,box_min_h_ms,box_min_adj_ms,"
                   "min_h_node,min_adj_node,base_lb_ms,adjusted_lb_ms,serf_apply_adj_bestcase,"
                   "serf_effective_lb_ms,serf_effective_lb_gt_T,fp_node\n";
        }

        unordered_map<string,int> fp_set_local;
        fp_set_local.reserve(fp_nodes.size()*2 + 8);
        for (const auto& n : fp_nodes) fp_set_local.emplace(n, 1);

        auto arr_to_str = [&](const array<long long,5>& a) -> std::string {
            std::string s;
            s.reserve(64);
            for (int d=0; d<5; d++) {
                if (d) s.push_back(' ');
                s += std::to_string(a[d]);
            }
            return s;
        };

        for (const auto& hit : box_hits) {
            if (fp_set_local.find(hit.node_name) == fp_set_local.end()) continue;

            auto itb = box_debug_map.find(hit.box_key);
            if (itb == box_debug_map.end()) continue;
            const BoxDebug& bd = itb->second;

            out << csv_escape(short_name(qname)) << ","
                << T_ms << ","
                << ORDER << ","
                << csv_escape(arr_to_str(bd.lo)) << ","
                << csv_escape(arr_to_str(bd.hi)) << ","
                << bd.depth << ","
                << bd.point_count << ","
                << bd.min_vec_dist << ","
                << bd.q_height_ms << ","
                << bd.q_adj_ms << ","
                << bd.box_min_h << ","
                << bd.box_min_adj << ","
                << csv_escape(short_name(bd.min_h_node)) << ","
                << csv_escape(short_name(bd.min_adj_node)) << ","
                << bd.base_lb_ms << ","
                << bd.adjusted_lb_ms << ","
                << bd.serf_apply_adj_bestcase << ","
                << bd.serf_effective_lb_ms << ","
                << bd.serf_effective_lb_gt_T << ","
                << csv_escape(short_name(hit.node_name))
                << "\n";
        }
    }
}

    cout << "Appended CSVs:\n";
    cout << "  " << out_summary_csv << "\n";
    cout << "  " << out_nodes_csv << "\n";
    if (!out_boxes_csv.empty()) cout << "  " << out_boxes_csv << "\n";

    DB->db_close();
    delete DB;
    return 0;
}
