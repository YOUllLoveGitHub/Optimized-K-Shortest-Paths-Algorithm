#include<ext/pb_ds/priority_queue.hpp>
#include <bits/stdc++.h>

using namespace __gnu_pbds;
using namespace std;

typedef long long ll;
typedef pair<int, int> pii;

const int MAXM = 30;
const int MAXN = 1e6 + 5;
const int MAXT = 1e6 + 5;
const int inf = 0x3f3f3f3f;

/** Profit is the value class of the smallest storage unit, currently contains only one integer variable, which can be modified if needed. */
struct Profit {
    int val;

    Profit(int X = 0)
        :val(X) {}

    bool operator != (const Profit& P) const {
        return val != P.val;
    }
    bool operator < (const Profit& P) const {
        return val < P.val;
    }
    bool operator <= (const Profit& P) const {
        return val <= P.val;
    }
    bool operator > (const Profit& P) const {
        return val > P.val;
    }
    Profit operator - () const {
        return Profit{-val};
    }
    Profit operator + (const Profit &a) const {
        return Profit{val + a.val};
    }
    Profit operator - (const Profit &a) const {
        return Profit{val - a.val};
    }
    Profit operator * (const int &f) const {
        return Profit{val * f};
    }
};

/** The global / local profit mentioned in the paper are the input that needs to be provided.
    local_profit[j][i]: the local profit that tile i is stored on server j;
    global_profit[i]: the profit of the first copy of tile i stored in the collaboration domain. */
Profit local_profit[MAXM][MAXN], global_profit[MAXN];

/** storages[j]: Storage capacity of MEC server j;
    N: number of tile files; M number of servers */
int storages[MAXM], N, M;

/** Optimal profit calculated by oksp algorithm */
ll ans = 0;


/** Auxiliary variables needed for OKSP algorithm */
int cnt[MAXN];
pair<Profit, int> PathArray[MAXM][MAXN];
int PathArrayInd[MAXM];
__gnu_pbds::priority_queue<pair<Profit, int>, greater<pair<Profit, int> >,
                           pairing_heap_tag> LossQueue[MAXM][MAXM];
int remain[MAXM];
bool contain[MAXM][MAXN];
bool debug = false;
pair<Profit, int> d[MAXM][MAXM];
Profit pi[MAXM];

/** Returns the profit and determines if global_profit needs to be added based on f */
Profit vals(int ind, int j, int f) {
    return global_profit[ind] * f + local_profit[j][ind];
}

/** The process of constructing G', mentioned in the paper, is the same as described in the pseudocode. */
void set_G() {
    // {Profit{inf}, -1} means unreachable
    d[0][M + 1] = {Profit{inf}, -1};
    for(int di = 1; di <= M; di++) {
        int i = di - 1;
        pair<Profit, int> best_val = {Profit{inf}, -1};
        if(PathArrayInd[i]) {
            best_val = PathArray[i][ PathArrayInd[i] - 1 ];
            best_val = {-best_val.first, best_val.second};
        }

        d[0][di] = best_val;
        d[di][M + 1] = {Profit{inf}, -1};
        if(remain[i]) d[di][M + 1] = {Profit{0}, -1};
        for(int dj = 1; dj <= M; dj++) {
            int j = dj - 1;
            d[di][dj] = {Profit{inf}, -1};
            if(i == j) continue;
            while(!LossQueue[i][j].empty()) {
                pair<Profit, int> pq_top = LossQueue[i][j].top();
                if(!(contain[i][pq_top.second] && !contain[j][pq_top.second])) {
                    LossQueue[i][j].pop();
                    continue;
                }
                d[di][dj] = {pq_top.first, pq_top.second};
                break;
            }
        }
    }
}

/** Reweight the graph, making all edge weights non-negative.
    Then we can apply the Dijkstra algorithm to compute the shortest path trees at each iteration. */
void set_pi(int m) {
    pi[0] = Profit{0};
    for(int k = 0; k < m; k++)
        for(int i = 0; i < m; i++)
            for(int j = 0; j < m; j++)
                if (d[i][j].first > d[i][k].first + d[k][j].first)
                    d[i][j] = {d[i][k].first + d[k][j].first, d[i][k].second};

    for(int i = 1; i < m; i++) pi[i] = d[0][i].first;
}

/** Auxiliary arrays for Dijkstra algorithm */
bool vis[MAXN];
Profit dis[MAXN];
int pre[MAXM];

/** Dijkstra for O(M ^ 2) */
void dijkstra(int m) {
    fill(dis, dis + m, Profit(inf));
    fill(vis, vis + m, false);
    dis[0] = Profit{0};

    while(true) {
        int v = -1;
        for(int i = 0; i < m; i++) {
            if(!vis[i] && (v == -1 || dis[i] < dis[v])) v = i;
        }
        if(v == -1) break;
        vis[v] = true;

        for(int i = 0; i < m; i++) {
            Profit w = d[v][i].first + pi[v] - pi[i];
            if (dis[i] > dis[v] + w) {
                dis[i] = dis[v] + w;
                pre[i] = v;;
            }
        }
    }
}

/** Backtracking shortest path of Dijkstra algorithm */
vector<pii> getPath(int e) {
    vector<pii> path;
    while (e) {
        path.push_back({pre[e] - 1, d[ pre[e] ][e].second});
        e = pre[e];
    }
    reverse(path.begin(), path.end());
    return path;
}

ll oksp() {
    clock_t t1 = clock();
    ans = 0;
    for(int i = 0; i < M + 2; i++)
        for(int j = 0; j < M + 2; j++)
            d[i][j] = {Profit{inf}, -1};

    // Finding the maximum profit path that could exist for each MEC server
    for(int i = 0; i < N; i++) cnt[i] = 0;
    for(int j = 0; j < M; j++) {
        PathArrayInd[j] = 0;
        remain[j] = storages[j];
        for(int i = 0; i < N; i++) {
            // Two kinds of profits of tile i store in MEC server j.
            contain[j][i] = false;
            PathArray[j][ PathArrayInd[j]++ ] = {vals(i, j, 0), i};
            PathArray[j][ PathArrayInd[j]++ ] = {vals(i, j, 1), i};
        }
        // Sort PathArray in descending order
        sort(PathArray[j], PathArray[j] + PathArrayInd[j]);
        for(int k = 0; k < M; k++) LossQueue[j][k].clear();
    }

    // Build G', and reweight the weight of its edges
    set_G();
    set_pi(M + 2);

    // K is the maximum number of paths
    int K = 0;
    for(int i = 0; i < M; i++) K += storages[i];
    K = min(M * N, K);

    // Checking the sorting time
    int totPathLen = 0;
    printf("The Sorting Time: %ld ms\n", clock() - t1);
    t1 = clock();
    printf("K: %d\n", K);

    // Iterate K times to find the K-shortest path
    while(K--) {
        // Rebuild G'
        for(int j = 0; j < M; j++) {
            while(PathArrayInd[j] > 0) {
                pair<Profit, int> init_val = PathArray[j][ PathArrayInd[j] - 1 ];
                int ind = init_val.second;
                if (cnt[ind]) {
                    Profit local_benefit = Profit{local_profit[j][ind]};
                    // If tile ind is already stored on another MEC server
                    // or server j has already stored the tile ind
                    if (init_val.first != local_benefit || contain[j][ind]) {
                        // The path is invalid
                        PathArrayInd[j]--;
                        continue;
                    }
                }
                break;
            }
        }
        set_G();

        //
        dijkstra(M + 2);
        Profit mi = -(dis[M + 1] + pi[M + 1] - pi[0]);
        if(mi < Profit{0}) {
            printf("There is no path exist, K: %d\n", K);
            break;
        }

        for(int i = 0; i < M + 2; i++) {
            pi[i] = pi[i] + dis[i];
            pi[i] = min(pi[i], Profit{inf});
        }

        vector<pii> path = getPath(M + 1);

        totPathLen += path.size();

        // Maintaining LossQueue
        cnt[path[0].second] += 1;
        remain[path[ path.size() - 1 ].first] -= 1;
        PathArrayInd[path[1].first]--;
        for(int i = 1; i < (int)path.size(); i++) {
            int mec_add = path[i].first;
            int item = path[i - 1].second;
            contain[mec_add][item] = true;
            // Adding the new transfer nodes.
            for(int j = 0; j < M; j++) {
                if(j == mec_add) continue;
                if(!contain[j][item]) {
                    Profit val = vals(item, mec_add, 0) - vals(item, j, 0);
                    LossQueue[mec_add][j].push({val, item});
                }
            }
            if(i != 1) {
                int mec_remove = path[i - 1].first;
                contain[mec_remove][item] = false;
                // Remove infeasible transfer nodes.
                for(int j = 0; j < M; j++) {
                    if(j == mec_remove) continue;
                    if(contain[j][item]) {
                        Profit val = vals(item, j, 0) - vals(item, mec_remove, 0);
                        LossQueue[j][mec_remove].push({val, item});
                    }
                }
            }
        }

        ans += mi.val;
    }
    printf("Total Path Length: %d\n", totPathLen);
    printf("Maximum Profit: %lld\n", ans);
    printf("OKSP() Time: %ld ms\n", clock() - t1);
    printf("\n");
    return ans;
}

// Converting problems into linear programming input files for correctness checking
void output_lp_file() {
    FILE* fp = fopen("lp.input", "w");
    fprintf(fp, "Maximize\n");

    fprintf(fp, "obj: ");
    for(int n = 0; n < N; n++) {
        for(int m = 0; m < M; m++) {
            fprintf(fp, "%d f2%d%d + ", local_profit[m][n].val, n, m);
        }
        fprintf(fp, "%d f1s%d0 %c ", global_profit[n].val, n, N - 1 == n ? '\n' : '+');
    };

    fprintf(fp, "\nSubject To\n");
    int cnt = 1;

    string in_flows, out_flows;

    // Add flow conservation constraints and flow variable constraints
    for(int n = 0; n < N; n++) {
        in_flows = "";
        for(int d = 0; d < M; d++) {
            in_flows += "f1s" + to_string(n) + to_string(d);
            if(d != M - 1) in_flows += " + ";
        }
        out_flows = "";
        for(int m = 0; m < M; m++) {
            out_flows += " - f2" + to_string(n) + to_string(m);
        }

        fprintf(fp, "\nc%d: %s\n%s = 0\n", cnt++, in_flows.c_str(), out_flows.c_str());
    }

    out_flows = "";
    for(int n = 0; n < N; n++) {
        for(int d = 0; d < M; d++) {
            out_flows += "f1s" + to_string(n) + to_string(d);
            if(d != M - 1) out_flows += " + ";
        }
    }
    in_flows = "";
    for(int m = 0; m < M; m++) {
        for(int d = 0; d < storages[m]; d++) {
            in_flows += " - f3" + to_string(m) + "t" + to_string(d);
        }
    }
    fprintf(fp, "\nc%d: %s\n%s = 0\n", cnt++, out_flows.c_str(), in_flows.c_str());

    for(int m = 0; m < M; m++) {
        in_flows = "";
        for(int n = 0; n < N; n++) {
            in_flows += "f2" + to_string(n) + to_string(m);
            if(n != N - 1) in_flows += " + ";
        }
        out_flows = "";
        for(int d = 0; d < storages[m]; d++) {
            in_flows += " - f3" + to_string(m) + "t" + to_string(d);
        }
        fprintf(fp, "\nc%d: %s\n%s = 0\n", cnt++, out_flows.c_str(), in_flows.c_str());
    }

    fprintf(fp, "\nBounds\n");

    for(int n = 0; n < N; n++) {
        for(int d = 0; d < M; d++) {
            fprintf(fp, "\n0 <= f1s%d%d <= 1\n", n, d);
        }
        for(int m = 0; m < M; m++) {
            fprintf(fp, "\n0 <= f2%d%d <= 1\n", n, m);
        }
    }

    for(int m = 0; m < M; m++) {
         for(int d = 0; d < storages[m]; d++) {
            fprintf(fp, "\n0 <= f3%dt%d <= 1\n", m, d);
         }
    }

    fprintf(fp, "End\n");

    fclose(fp);
}

int main() {
    srand(time(0));
    N = 5000; M = 10;
    for(int j = 0; j < M; j++) {
        storages[j] = 200;
        for(int i = 0; i < N; i++) {
            local_profit[j][i] = rand() % 1000 + 1;
        }
    }
    for(int i = 0; i < N; i++) global_profit[i] = rand() % 1000 * M  + 1;
    oksp();
    output_lp_file();
}
