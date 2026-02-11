// tests/serf_driver_naive.cc
// Naive baseline range query for Serf Vivaldi coordinates.
// For a given query node q and RTT threshold T_ms:
//   - Compute predicted RTT(q, v) for every node v using Serf's formula (Vec + Height + Adjustment + guard).
//   - Convert the final RTT from seconds to milliseconds (rtt * 1000.0).
//   - Return nodes with RTT <= T_ms.
// Times only the in-memory query work (no CSV writing inside the timer).
// No shifting, no Hilbert, no Lawder, no indexing.

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>

#include "../tests/third-party/json.hpp"
using json = nlohmann::json;
using namespace std;

// ---------------- helpers ----------------
static bool file_empty_or_missing(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return true;
    return (st.st_size == 0);
}

static bool read_all(const std::string& path, std::string& out) {
    ifstream f(path.c_str(), ios::in | ios::binary);
    if (!f) return false;
    ostringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    return true;
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
        out_nodes = root["nodes"];
        return true;
    }
    return false;
}

static std::string short_name(const std::string& full) {
    size_t pos = full.rfind('-');
    if (pos == std::string::npos) return full;
    if (pos + 1 >= full.size()) return full;
    return full.substr(pos + 1);
}

static std::string csv_escape(const std::string& s) {
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

static void usage(const char* prog) {
    cerr
        << "Usage:\n"
        << "  " << prog << " --qnode <name> --rtt <ms>\n"
        << "         [--json <path>]\n"
        << "         [--out_naive_summary_csv <path>] [--out_naive_nodes_csv <path>]\n"
        << "\n"
        << "Notes:\n"
        << "  - Reads Serf coordinates from JSON in seconds (Vec/Height/Adjustment).\n"
        << "  - Computes RTT in seconds, applies Serf adjustment guard, then converts to milliseconds (rtt * 1000.0).\n"
        << "  - Timing includes only: RTT computations over all nodes + threshold filtering + building returned set in memory.\n"
        << "  - Timing excludes any CSV writing.\n";
}

// ---------------- Serf RTT (seconds -> milliseconds) ----------------
// Inputs are in seconds (as in Serf JSON).
// Returns RTT in milliseconds.
static inline double serf_rtt_ms_from_sec(
    const array<double,5>& a_vec_sec, double a_h_sec, double a_adj_sec,
    const array<double,5>& b_vec_sec, double b_h_sec, double b_adj_sec)
{
    double sumsq = 0.0;
    for (int i = 0; i < 5; i++) {
        const double diff = a_vec_sec[i] - b_vec_sec[i];
        sumsq += diff * diff;
    }

    double rtt = std::sqrt(sumsq) + a_h_sec + b_h_sec;

    const double adjusted = rtt + a_adj_sec + b_adj_sec;
    if (adjusted > 0.0) {
        rtt = adjusted;
    }

    // RTT is in seconds by model — convert to milliseconds
    return rtt * 1000.0;
}

struct NodeCoord {
    std::string name;
    array<double,5> vec_sec{};
    double height_sec = 0.0;
    double adj_sec = 0.0;
};

int main(int argc, char** argv) {
    std::string qname;
    double T_ms = -1.0;

    std::string json_path = "cluster-status-2025-07-11-1-rtts_recomputed.json";
    std::string out_naive_summary_csv = "tableA_naive_summary.csv";
    std::string out_naive_nodes_csv; // optional

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];

        auto need = [&](const char* flag) {
            if (i + 1 >= argc) {
                cerr << "Missing value for " << flag << "\n";
                usage(argv[0]);
                exit(1);
            }
        };

        if (a == "--qnode") {
            need("--qnode");
            qname = argv[++i];
        } else if (a == "--rtt") {
            need("--rtt");
            T_ms = atof(argv[++i]);
        } else if (a == "--json") {
            need("--json");
            json_path = argv[++i];
        } else if (a == "--out_naive_summary_csv") {
            need("--out_naive_summary_csv");
            out_naive_summary_csv = argv[++i];
        } else if (a == "--out_naive_nodes_csv") {
            need("--out_naive_nodes_csv");
            out_naive_nodes_csv = argv[++i];
        } else {
            cerr << "Unknown arg: " << a << "\n";
            usage(argv[0]);
            return 1;
        }
    }

    if (qname.empty() || T_ms < 0.0) {
        usage(argv[0]);
        return 1;
    }

    json root;
    if (!read_json_file(json_path, root)) {
        cerr << "Error: could not read JSON: " << json_path << "\n";
        return 1;
    }

    json nodes;
    if (!extract_nodes_array(root, nodes)) {
        cerr << "Error: JSON root must be an array or {\"nodes\":[...]}\n";
        return 1;
    }

    vector<NodeCoord> all;
    all.reserve(nodes.size());

    for (const auto& nd : nodes) {
        if (!nd.is_object()) continue;
        if (!nd.contains("name") || !nd["name"].is_string()) continue;
        if (!nd.contains("coordinate") || !nd["coordinate"].is_object()) continue;

        NodeCoord nc;
        nc.name = nd["name"].get<std::string>();

        const auto& c = nd["coordinate"];
        if (!c.contains("Vec") || !c["Vec"].is_array() || c["Vec"].size() < 5) {
            cerr << "Error: node " << nc.name << " missing coordinate.Vec[5]\n";
            return 1;
        }

        for (int d = 0; d < 5; d++) {
            nc.vec_sec[d] = c["Vec"][d].get<double>(); // seconds
        }

        nc.height_sec = c.contains("Height") ? c["Height"].get<double>() : 0.0;         // seconds
        nc.adj_sec    = c.contains("Adjustment") ? c["Adjustment"].get<double>() : 0.0; // seconds

        all.push_back(nc);
    }

    if (all.empty()) {
        cerr << "Error: no nodes parsed from JSON.\n";
        return 1;
    }

    int qi = -1;
    for (size_t i = 0; i < all.size(); i++) {
        if (all[i].name == qname) { qi = (int)i; break; }
    }
    if (qi < 0) {
        cerr << "Error: query node not found: " << qname << "\n";
        return 1;
    }

    // ---------------- Naive timed region ----------------
    vector<string> returned;
    returned.reserve(all.size());

    const auto t0 = std::chrono::high_resolution_clock::now();

    int naive_rtt_evals = 0;
    const auto& q = all[(size_t)qi];

    for (size_t i = 0; i < all.size(); i++) {
        if ((int)i == qi) continue;

        naive_rtt_evals++;
        const auto& v = all[i];

        const double rtt_ms = serf_rtt_ms_from_sec(
            q.vec_sec, q.height_sec, q.adj_sec,
            v.vec_sec, v.height_sec, v.adj_sec
        );

        if (rtt_ms <= T_ms) {
            returned.push_back(v.name);
        }
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    const double t_naive_query_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    // ---------------- End timed region ----------------

    // Write summary CSV (outside timed region)
    {
        const bool need_header = file_empty_or_missing(out_naive_summary_csv);
        ofstream out(out_naive_summary_csv.c_str(), ios::out | ios::app);
        if (!out) {
            cerr << "Warning: could not append naive summary CSV: " << out_naive_summary_csv << "\n";
        } else {
            if (need_header) {
                out << "query_node,rtt_ms,naive_returned,naive_rtt_evals,t_naive_query_ms\n";
            }
            out << csv_escape(short_name(qname)) << ","
                << T_ms << ","
                << returned.size() << ","
                << naive_rtt_evals << ","
                << t_naive_query_ms
                << "\n";
        }
    }

    // Optional nodes CSV (outside timed region)
    if (!out_naive_nodes_csv.empty()) {
        const bool need_header = file_empty_or_missing(out_naive_nodes_csv);
        ofstream out(out_naive_nodes_csv.c_str(), ios::out | ios::app);
        if (!out) {
            cerr << "Warning: could not append naive nodes CSV: " << out_naive_nodes_csv << "\n";
        } else {
            if (need_header) {
                out << "query_node,rtt_ms,naive_returned,returned_nodes\n";
            }

            vector<string> shorted;
            shorted.reserve(returned.size());
            for (const auto& nm : returned) shorted.push_back(short_name(nm));

            out << csv_escape(short_name(qname)) << ","
                << T_ms << ","
                << returned.size() << ","
                << csv_escape(join_names_sorted(shorted))
                << "\n";
        }
    }

    cout << "Naive query done: qnode=" << qname
         << " T_ms=" << T_ms
         << " evals=" << naive_rtt_evals
         << " returned=" << returned.size()
         << " t_naive_query_ms=" << t_naive_query_ms
         << "\n";

    return 0;
}
