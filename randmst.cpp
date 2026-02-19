#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <cstdint>
#include <string>
#include <random>
#include <algorithm>
#include <random>
#include <chrono>


using namespace std;
using u64 = uint64_t;

// splitmix64: fast deterministic mixing function (good for hashing)
// splitmix64: a fast 64-bit mixing function.
// We use this to deterministically convert (seed, u, v) into a
// well-distributed pseudo-random number.
// This avoids storing edge weights while keeping them consistent
// within each trial.
static inline u64 splitmix64(u64 x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// edge_weight01(u, v, seed):
// Returns a deterministic pseudo-random weight in [0,1) for edge (u,v).
// Important properties:
// 1. Symmetric: weight(u,v) == weight(v,u)
// 2. Same edge always gets same weight within a trial
// 3. Different trials use different seeds
// This lets us simulate a random complete graph WITHOUT storing all edges.

static inline double edge_weight01(int u, int v, u64 seed) {
    if (u > v) std::swap(u, v);

    // pack (seed, u, v) into one 64-bit stream and mix it
    u64 x = seed;
    x ^= (u64)u * 0x9e3779b97f4a7c15ULL;
    x ^= (u64)v * 0xbf58476d1ce4e5b9ULL;

    u64 h = splitmix64(x);

    // convert top 53 bits to a double in [0,1)
    const u64 mantissa = h >> 11;               // keep 53 bits
    return (double)mantissa * (1.0 / (1ULL << 53));
}
// hypercube_neighbors(v, n, out):
// Generates all neighbors of vertex v in the assignment's "hypercube" graph.
// Two vertices are connected if |a - b| = 2^i for some i.
// So neighbors are v ± 1, v ± 2, v ± 4, v ± 8, ...
// Degree is O(log n), making this graph sparse.

static inline void hypercube_neighbors(int v, int n, vector<int>& out) {
    out.clear();
    for (int step = 1; step < n; step <<= 1) {
        int a = v - step;
        int b = v + step;
        if (a >= 0) out.push_back(a);
        if (b < n) out.push_back(b);
    }
}
// mst_hypercube(n, seed):
// Computes MST weight for the hypercube graph (dimension 1)
// using Prim's algorithm with a min-heap (priority_queue).
//
// Why heap-based Prim?
// - Hypercube graph is sparse (≈ n log n edges).
// - Heap Prim runs in roughly O(n (log n)^2).
// - Efficient enough for n up to 262144 as required.
//
// We:
// 1. Start from vertex 0.
// 2. Push its outgoing edges into a min-heap.
// 3. Repeatedly take the smallest edge leading to a new vertex.
// 4. Add that vertex and push its edges.
// 5. Stop when all vertices are included.

static double mst_hypercube(int n, u64 seed) {
    if (n <= 1) return 0.0;

    vector<char> in_tree(n, 0);
    vector<int> neigh;
    neigh.reserve(64);

    // min-heap of (weight, vertex)
    using Item = pair<double, int>;
    priority_queue<Item, vector<Item>, greater<Item>> pq;

    in_tree[0] = 1;
    int added = 1;
    double total = 0.0;

    hypercube_neighbors(0, n, neigh);
    for (int nb : neigh) {
        pq.push({edge_weight01(0, nb, seed), nb});
    }

    while (added < n) {
        auto [w, v] = pq.top();
        pq.pop();
        if (in_tree[v]) continue;

        in_tree[v] = 1;
        total += w;
        added++;

        hypercube_neighbors(v, n, neigh);
        for (int nb : neigh) {
            if (!in_tree[nb]) {
                pq.push({edge_weight01(v, nb, seed), nb});
            }
        }
    }

    return total;
}
// mst_complete_dim0(n, seed):
// Naive Prim for complete graph where each edge weight is uniform in [0,1].
// We do NOT store edges; we compute weights on demand using edge_weight01.
// Runtime: O(n^2), Memory: O(n).
static double mst_complete_dim0(int n, u64 seed) {
    if (n <= 1) return 0.0;

    vector<char> in_tree(n, 0);
    vector<double> best(n, 1e100);

    best[0] = 0.0;
    double total = 0.0;

    for (int iter = 0; iter < n; ++iter) {
        int u = -1;
        double u_best = 1e100;

        // pick next vertex with smallest best[] not in tree
        for (int i = 0; i < n; ++i) {
            if (!in_tree[i] && best[i] < u_best) {
                u_best = best[i];
                u = i;
            }
        }

        in_tree[u] = 1;
        total += u_best;

        // relax edges (u, v) for all v not in tree
        for (int v = 0; v < n; ++v) {
            if (!in_tree[v]) {
                double w = edge_weight01(u, v, seed);
                if (w < best[v]) best[v] = w;
            }
        }
    }

    return total;
}
struct Point {
    double x, y, z, w;
};

// Generate n random points in [0,1]^dim using a per-trial RNG
static vector<Point> gen_points(int n, int dim, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    vector<Point> pts(n);
    for (int i = 0; i < n; ++i) {
        Point p{0,0,0,0};
        p.x = dist(rng);
        if (dim >= 2) p.y = dist(rng);
        if (dim >= 3) p.z = dist(rng);
        if (dim >= 4) p.w = dist(rng);
        pts[i] = p;
    }
    return pts;
}
// adding a point0type +point generator
static inline double sqr(double a) { return a * a; }

static inline double euclid(const Point& a, const Point& b, int dim) {
    double s = sqr(a.x - b.x);
    if (dim >= 2) s += sqr(a.y - b.y);
    if (dim >= 3) s += sqr(a.z - b.z);
    if (dim >= 4) s += sqr(a.w - b.w);
    return std::sqrt(s);
}
// mst_geometric(n, dim, rng):
// Complete graph on random points in [0,1]^dim with Euclidean edge weights.
// We generate points once per trial, then run naive Prim in O(n^2).
static double mst_geometric(int n, int dim, std::mt19937_64& rng) {
    if (n <= 1) return 0.0;

    vector<Point> pts = gen_points(n, dim, rng);

    vector<char> in_tree(n, 0);
    vector<double> best(n, 1e100);

    best[0] = 0.0;
    double total = 0.0;

    for (int iter = 0; iter < n; ++iter) {
        int u = -1;
        double u_best = 1e100;

        for (int i = 0; i < n; ++i) {
            if (!in_tree[i] && best[i] < u_best) {
                u_best = best[i];
                u = i;
            }
        }

        in_tree[u] = 1;
        total += u_best;

        for (int v = 0; v < n; ++v) {
            if (!in_tree[v]) {
                double w = euclid(pts[u], pts[v], dim);
                if (w < best[v]) best[v] = w;
            }
        }
    }

    return total;
}


// at top of file make sure you have:
// #include <random>
// #include <chrono>

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc != 5) {
        cerr << "Usage: ./randmst 0 numpoints numtrials dimension\n";
        return 2;
    }

    int flag = stoi(argv[1]); (void)flag;
    int n = stoi(argv[2]);
    int trials = stoi(argv[3]);
    int dim = stoi(argv[4]);

    if (n < 0 || trials <= 0) {
        cerr << "Error: numpoints must be >= 0 and numtrials must be > 0\n";
        return 2;
    }

    // --- RNG: create a master RNG with a nondeterministic-ish seed ---
    std::random_device rd;
    uint64_t time_seed = static_cast<uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()
    );
    // Combine multiple sources of entropy
    uint64_t seed0 = (uint64_t(rd()) << 32) ^ uint64_t(rd()) ^ time_seed;
    std::mt19937_64 master_rng(seed0);

    double sum = 0.0;
    for (int t = 0; t < trials; ++t) {
        // Draw a fresh 64-bit seed for this trial
        uint64_t trial_seed = master_rng();

        if (dim == 0) {
            // For dim 0: use the seed in your deterministic edge-weight hash
            sum += mst_complete_dim0(n, trial_seed);
        } else if (dim == 1) {
            // For dim 1: same idea, pass the trial seed to your hypercube MST
            sum += mst_hypercube(n, trial_seed);
        } else if (dim == 2 || dim == 3 || dim == 4) {
            // Create a trial-local RNG seeded by trial_seed for reproducible point generation
            std::mt19937_64 trial_rng(trial_seed);
            sum += mst_geometric(n, dim, trial_rng);
        } else {
            cerr << "Error: dimension must be 0,1,2,3,4\n";
            return 2;
        }
    }

    double average = sum / static_cast<double>(trials);
    cout << average << " " << n << " " << trials << " " << dim << "\n";
    return 0;
}
