#include "dataStruct.h"
#include <stdlib.h>
#include <math.h>

void histogram_single(vector<int32_t> &order, Histogram &Hist, uint16_t *maps, double *G, double *H) {
    uint16_t bin;
    for (auto i : order) {
        bin = maps[i];
        ++Hist.cnt[bin];
        Hist.g[bin] += G[i];
        Hist.h[bin] += H[i];
    }
    // integration
    for (int i = 1; i < Hist.cnt.size(); ++i) {
        Hist.cnt[i] += Hist.cnt[i - 1];
        Hist.g[i] += Hist.g[i - 1];
        Hist.h[i] += Hist.h[i - 1];
    }
}

void histogram_multi(vector<int32_t> &order, Histogram &Hist, uint16_t *maps, double *G, double *H, int out_dim) {
    int j, ind, bin;
    int Xp[out_dim];
    create_vsrp(Xp, out_dim);
    for (auto i : order) {
        bin = maps[i] * out_dim;
        ind = i * out_dim;
        ++Hist.cnt[maps[i]];
        for (j = 0; j < out_dim; ++j) {
            if (Xp[j] == 1) {
                Hist.g[bin] += G[ind];
                Hist.h[bin++] += H[ind++];
            }else if (Xp[j] == -1) {
                Hist.g[bin] -= G[ind];
                Hist.h[bin++] -= H[ind++];
            }else {
                ++bin;
                ++ind;
            }
        }
    }
    // integration
    ind = 0;
    for (int i = 1; i < Hist.cnt.size(); ++i) {
        Hist.cnt[i] += Hist.cnt[i - 1];
        for (j = 0; j < out_dim; ++j) {
            if (Xp[j] == 1) {
                Hist.g[ind + out_dim] += Hist.g[ind];
                Hist.h[ind + out_dim] += Hist.h[ind];
                ++ind;
            } else if (Xp[j] == -1) {
                Hist.g[ind + out_dim] -= Hist.g[ind];
                Hist.h[ind + out_dim] -= Hist.h[ind];
                ++ind;
            } else {
                ++ind;
            }
        }
    }
}

void create_vsrp(int *X, int M)
{
    int j;
    double r;
    double tmp1 = 1.0 / (2.0 * sqrt(M));

    for (j = 0; j < M; j++) {
        r = drand48();
        if (r < tmp1)
            X[j] = -1;
        else if (r >= 1.0 - tmp1)
            X[j] = 1;
        else
            X[j] = 0;
    }
}