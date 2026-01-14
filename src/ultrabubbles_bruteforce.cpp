#include <bits/stdc++.h>
#include <cassert>
#include <cstdint>
#ifdef __unix__
#include <unistd.h>
#include <limits.h>
#include <sys/wait.h>
#endif

using namespace std;

#ifndef DEBUG_BLOCKS
#define DEBUG_BLOCKS 1
#endif

static inline int sgnIdx(char c) { return (c == '-') ? 1 : 0; }
static inline char oppSign(char c) { return (c == '+') ? '-' : '+'; }

struct IncEdge {
    int v;
    char su; // side at u
    char sv; // side at v
};

struct SnarlKey {
    int a, b;
    char da, db;
    bool operator<(SnarlKey const& o) const {
        return std::tie(a,da,b,db) < std::tie(o.a,o.da,o.b,o.db);
    }
    bool operator==(SnarlKey const& o) const {
        return a==o.a && b==o.b && da==o.da && db==o.db;
    }
};

static inline SnarlKey canon_snarl(int x, char dx, int y, char dy) {
    if (x < y) return {x,y,dx,dy};
    return {y,x,dy,dx};
}

#ifdef __unix__
static std::string get_self_dir() {
    char buf[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (len == -1) return {};
    buf[len] = '\0';
    std::string abs_path(buf);
    auto pos = abs_path.find_last_of('/');
    if (pos == std::string::npos) return ".";
    return abs_path.substr(0, pos);
}

static bool parse_side_token(const std::string& t,
                            const std::unordered_map<std::string,int>& id,
                            int &v, char &s) {
    if (t.size() < 2) return false;
    char c = t.back();
    if (c != '+' && c != '-') return false;
    std::string nm = t.substr(0, t.size()-1);
    auto it = id.find(nm);
    if (it == id.end()) return false;
    v = it->second;
    s = c;
    return true;
}

static bool run_snarls_bf(const std::string& snarls_bf_path,
                          const std::string& gfa_path,
                          const std::unordered_map<std::string,int>& id,
                          std::set<SnarlKey>& out) {
    out.clear();

    if (snarls_bf_path.empty()) return false;
    if (access(snarls_bf_path.c_str(), X_OK) != 0) return false;

    std::string cmd = "'" + snarls_bf_path + "' '" + gfa_path + "'";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return false;

    char buf[1<<15];
    while (fgets(buf, sizeof(buf), pipe)) {
        std::string line(buf);
        std::stringstream ss(line);
        std::string t1, t2;
        ss >> t1 >> t2;
        if (!ss) continue;

        int x,y; char dx,dy;
        if (!parse_side_token(t1, id, x, dx)) continue;
        if (!parse_side_token(t2, id, y, dy)) continue;

        out.insert(canon_snarl(x, dx, y, dy));
    }

    int status = pclose(pipe);
    if (status == -1) return false;
    if (!WIFEXITED(status)) return false;
    if (WEXITSTATUS(status) != 0) return false;
    return true;
}
#endif

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc != 2) {
        cerr << "Usage: ./ultrabubbles_bf input.gfa\n";
        return 1;
    }

#ifdef __unix__
    const std::string self_dir = get_self_dir();
    const std::string snarls_bf_path = self_dir.empty()
                                        ? ""
                                        : (self_dir + "/snarls_bf");
#else
    const std::string snarls_bf_path;
#endif

    ifstream fin(argv[1]);
    if (!fin) {
        cerr << "Cannot open input.\n";
        return 1;
    }

    unordered_map<string,int> id;
    vector<string> nameOf;
    auto getId = [&](const string& s)->int{
        auto it = id.find(s);
        if (it != id.end()) return it->second;
        int nid = (int)nameOf.size();
        id[s] = nid;
        nameOf.push_back(s);
        return nid;
    };

    vector<vector<IncEdge>> adj;
    vector<uint8_t> seenSegment;

    auto ensureN = [&](int n){
        if ((int)adj.size() < n) adj.resize(n);
        if ((int)seenSegment.size() < n) seenSegment.resize(n, 0);
    };

    string line;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        string type;
        ss >> type;
        if (!ss) continue;

        if (type == "S") {
            string seg, seq;
            ss >> seg >> seq;
            if (!ss) continue;
            int u = getId(seg);
            ensureN(u+1);
            seenSegment[u] = 1;
        } else if (type == "L") {
            string a, b, ovl;
            char oa = 0, ob = 0;
            ss >> a >> oa >> b >> ob >> ovl;
            if (!ss) continue;

            int u = getId(a);
            int v = getId(b);
            ensureN(max(u,v)+1);
            seenSegment[u] = seenSegment[v] = 1;

            // bidirected incidence convention:
            // at u we store oa, at v we store flip(ob)
            char su = oa;
            char sv = (ob == '-') ? '+' : '-';
            adj[u].push_back({v, su, sv});
            adj[v].push_back({u, sv, su});
        }
    }

    int N = (int)nameOf.size();
    ensureN(N);

    vector<int> nodes;
    nodes.reserve(N);
    for (int v = 0; v < N; v++) if (seenSegment[v]) nodes.push_back(v);

    sort(nodes.begin(), nodes.end());
    assert(is_sorted(nodes.begin(), nodes.end()));

    cerr << "Vertices: " << nodes.size() << "\n";

    const int REACH_FLAG = 1;
    const int BAD_FLAG   = 2;

    vector<int> seen(N, 0);
    int dfs_it = 1;
    vector<int> snarl_component;

    // DFS used for separability + component construction (net graph semantics: split terminals)
    function<int(int,int,int,char,char,bool)> dfs = [&](int z, int x, int y, char dx, char dy, bool track)->int {
        if (seen[z] == dfs_it) return 0;
        seen[z] = dfs_it;
        if (track) snarl_component.push_back(z);

        int ans = 0;
        if (z == x) {
            for (auto &inc : adj[z]) {
                if (inc.su != dx) continue;
                if (inc.v == y && inc.sv == dy) ans |= REACH_FLAG;
                if (inc.v == y && inc.sv != dy) ans |= BAD_FLAG;
                ans |= dfs(inc.v, x, y, dx, dy, track);
            }
        } else if (z == y) {
            for (auto &inc : adj[z]) {
                if (inc.su != dy) continue;
                if (inc.v == x && inc.sv != dx) ans |= BAD_FLAG;
                ans |= dfs(inc.v, x, y, dx, dy, track);
            }
        } else {
            for (auto &inc : adj[z]) {
                if (inc.v == y && inc.sv == dy) ans |= REACH_FLAG;
                if (inc.v == y && inc.sv != dy) ans |= BAD_FLAG;
                if (inc.v == x && inc.sv != dx) ans |= BAD_FLAG;
                ans |= dfs(inc.v, x, y, dx, dy, track);
            }
        }
        return ans;
    };

    auto test_separability = [&](int x, int y, char dx, char dy)->bool{
        int out = dfs(x, x, y, dx, dy, false);
        dfs_it++;
        return (out & REACH_FLAG) && !(out & BAD_FLAG);
    };

    auto compute_component = [&](int x, int y, char dx, char dy)->void{
        snarl_component.clear();
        dfs(x, x, y, dx, dy, true);
        dfs_it++;
    };

    // separability memo: N*N*4 (dx,dy in {0,1})
    vector<int8_t> sep((size_t)N * (size_t)N * 4, (int8_t)-1);

    auto sepIndex = [&](int x, int y, int dx, int dy)->size_t{
        return ((size_t)x * (size_t)N + (size_t)y) * 4 + (size_t)(dx*2 + dy);
    };

    auto getSep = [&](int x, int y, int dx, int dy)->bool{
        int8_t &cell = sep[sepIndex(x,y,dx,dy)];
        if (cell == (int8_t)-1) {
            char sx = dx ? '-' : '+';
            char sy = dy ? '-' : '+';
            cell = test_separability(x, y, sx, sy) ? (int8_t)1 : (int8_t)0;
        }
        return cell == (int8_t)1;
    };

    // ---- split-terminals filter (ONLY for net-graph checks like tip-free) ----
    int CX=-1, CY=-1;
    char CDX='?', CDY='?';

    auto edge_survives_split = [&](int u, int v, char su, char sv)->bool{
        if (u == CX && su != CDX) return false;
        if (v == CX && sv != CDX) return false;
        if (u == CY && su != CDY) return false;
        if (v == CY && sv != CDY) return false;
        return true;
    };

    vector<uint8_t> in_comp(N, 0);

    // tip-free check in the NET GRAPH (so we keep split filtering)
    auto has_internal_tip_in_component = [&]()->bool{
        vector<array<uint8_t,2>> hasSign(N, {0,0});

        for (int z : snarl_component) {
            for (auto &inc : adj[z]) {
                int w = inc.v;
                if (!in_comp[w]) continue;
                if (!edge_survives_split(z, w, inc.su, inc.sv)) continue;
                hasSign[z][sgnIdx(inc.su)] = 1;
            }
        }

        for (int z : snarl_component) {
            if (z == CX || z == CY) continue;
            if (!hasSign[z][0] || !hasSign[z][1]) return true;
        }
        return false;
    };

    // ---- directed cycle check in the INDUCED SUBGRAPH (NO split filtering) ----
    // This is what removes "weak ultrabubbles" like your {3,5} example.
    vector<array<uint8_t,2>> onstack(N, {0,0});
    vector<array<int,2>> seen2(N, {0,0});
    int dfs_it2 = 1;

    auto state_has_outgoing_induced = [&](int u, char in_sign)->bool{
        char out_needed = oppSign(in_sign);
        for (auto &inc : adj[u]) {
            int v = inc.v;
            if (!in_comp[v]) continue;          // induced subgraph
            if (inc.su == out_needed) return true;
        }
        return false;
    };

    function<bool(int,char,int,char,int)> dfs_cycle_or_opposite_induced =
        [&](int u, char in_sign, int start_u, char start_sign, int depth)->bool {

        int idx = sgnIdx(in_sign);

        if (onstack[u][idx]) return true;

        // optional strengthening: if we can return to the same vertex with opposite sign,
        // it corresponds to a directed cycle in the state graph.
        if (depth > 0 && u == start_u && in_sign == oppSign(start_sign)) return true;

        if (seen2[u][idx] == dfs_it2) return false;
        seen2[u][idx] = dfs_it2;
        onstack[u][idx] = 1;

        char out_needed = oppSign(in_sign);

        for (auto &inc : adj[u]) {
            int v = inc.v;
            if (!in_comp[v]) continue;          // induced subgraph
            if (inc.su != out_needed) continue; // must exit u on opposite side

            if (dfs_cycle_or_opposite_induced(v, inc.sv, start_u, start_sign, depth + 1)) {
                onstack[u][idx] = 0;
                return true;
            }
        }

        onstack[u][idx] = 0;
        return false;
    };

    auto has_directed_cycle_in_component_induced = [&]()->bool{
        for (int u : snarl_component) {
            for (char s : {'+','-'}) {
                if (!state_has_outgoing_induced(u, s)) continue;

                dfs_it2++;
                for (int z : snarl_component) onstack[z][0] = onstack[z][1] = 0;

                if (dfs_cycle_or_opposite_induced(u, s, u, s, 0)) return true;
            }
        }
        return false;
    };

    cerr << "Enumerating ultrabubbles\n";

    const size_t nV = nodes.size();
    const size_t expected_candidates = (nV * (nV - 1) / 2) * 4;
    size_t tested_candidates = 0;

    struct Eval { bool checked_tip=false, checked_cycle=false; };
    std::map<SnarlKey, Eval> snarl_eval;
    std::set<SnarlKey> snarls_ours;

    for (size_t ii = 0; ii < nodes.size(); ii++) {
        int x = nodes[ii];
        for (size_t jj = ii + 1; jj < nodes.size(); jj++) {
            int y = nodes[jj];

            for (int dx = 0; dx < 2; dx++) {
                for (int dy = 0; dy < 2; dy++) {
                    tested_candidates++;

                    bool separable = getSep(x, y, dx, dy);

#if DEBUG_BLOCKS
                    assert(separable == getSep(y, x, dy, dx));
#endif
                    if (!separable) continue;

                    CDX = dx ? '-' : '+';
                    CDY = dy ? '-' : '+';
                    compute_component(x, y, CDX, CDY);

                    fill(in_comp.begin(), in_comp.end(), 0);
                    for (int z : snarl_component) in_comp[z] = 1;

#if DEBUG_BLOCKS
                    bool hasX=false, hasY=false;
                    for (int z : snarl_component) { if (z==x) hasX=true; if (z==y) hasY=true; }
                    assert(hasX && hasY);
#endif

                    // minimality
                    bool not_minimal = false;
                    for (int z : snarl_component) {
                        if (z == x || z == y) continue;
                        if (getSep(x, z, dx, 0) && getSep(z, y, 1, dy)) { not_minimal = true; break; }
                        if (getSep(x, z, dx, 1) && getSep(z, y, 0, dy)) { not_minimal = true; break; }
                    }
                    if (not_minimal) continue;

                    SnarlKey sk = canon_snarl(x, CDX, y, CDY);
                    snarls_ours.insert(sk);

                    // terminals for split-filtered checks (tip-free in net graph)
                    CX = x; CY = y;

                    bool internal_tip = has_internal_tip_in_component();
                    snarl_eval[sk].checked_tip = true;

                    // IMPORTANT FIX: cycle check is now on induced subgraph => rejects weak ultrabubbles
                    bool directed_cycle = has_directed_cycle_in_component_induced();
                    snarl_eval[sk].checked_cycle = true;

                    if (!internal_tip && !directed_cycle) {
                        cout << nameOf[x] << (dx ? "-" : "+") << " "
                             << nameOf[y] << (dy ? "-" : "+") << "\n";
                    }
                }
            }
        }
    }

    assert(tested_candidates == expected_candidates);

    for (auto const& [k, ev] : snarl_eval) {
        (void)k;
        assert(ev.checked_tip);
        assert(ev.checked_cycle);
    }
    for (auto const& k : snarls_ours) {
        auto it = snarl_eval.find(k);
        assert(it != snarl_eval.end());
        assert(it->second.checked_tip && it->second.checked_cycle);
    }

#if DEBUG_BLOCKS
#ifdef __unix__
    std::set<SnarlKey> snarls_ref;
    bool ok = run_snarls_bf(snarls_bf_path, argv[1], id, snarls_ref);
    if (ok) {
        if (snarls_ref != snarls_ours) {
            cerr << "\n\n\n Mismatch with ./snarls_bf on snarl set\n";
            for (auto &s : snarls_ours) if (!snarls_ref.count(s))
                cerr << "Extra:   " << nameOf[s.a] << s.da << " " << nameOf[s.b] << s.db << "\n";
            for (auto &s : snarls_ref) if (!snarls_ours.count(s))
                cerr << "Missing: " << nameOf[s.a] << s.da << " " << nameOf[s.b] << s.db << "\n";
            assert(false && "snarl set mismatch with snarls bf");
        }
    } else {
        cerr << "Could not run ./snarls_bf from same directory. Skipping external snarl cross check.\n";
    }
#endif
#endif

    return 0;
}