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

void histogram_multi(vector<int32_t> &order, Histogram &Hist, uint16_t *maps, const double *G, const double *H, const int *Xp, int out_dim) {
    int j, ind, bin;
    for (auto i : order) {
        bin = maps[i] * out_dim;
        ind = i * out_dim;
        ++Hist.cnt[maps[i]];
        for (j = 0; j < out_dim; j++) {
//            printf("%f\n",H[ind]);
            Hist.g[bin] += G[ind] * Xp[j];
            Hist.h[bin++] += H[ind++];// * Xp[j];
        }
    }


//    exit(0);
    // integration
    ind = 0;
    for (int i = 1; i < Hist.cnt.size(); ++i) {
        Hist.cnt[i] += Hist.cnt[i - 1];
        for (j = 0; j < out_dim; j++) {
//            printf("%f\n",Xp[j]);
            Hist.g[ind + out_dim] += Hist.g[ind];// * Xp[j];
            Hist.h[ind + out_dim] += Hist.h[ind];// * Xp[j];
            ++ind;
        }
    }
//    printf("%d\n" ,Hist.cnt.size());
//
//    int g_size=Hist.g.size();
//    for (int i = 1; i < Hist.g.size(); ++i) {
//        printf("%f\n", Hist.g[j]);
//    }
//    exit(0);
    }



//    for (auto i : order) {
//        int g_size= sizeof(&G)/sizeof(G[0]);
//        int h_size= sizeof(&H)/sizeof(H[0]);
//        printf("%d %d - %d %d\n",g_size, h_size, i,i + out_dim);
//        for (j = 0; j < out_dim; ++j) {
//            if (Xp[j] == 1) {
//                Hist.g[j] += G[i + out_dim];
//                Hist.h[j] += H[i + out_dim];
//            } else if (Xp[j] == -1) {
//                Hist.g[j] -= G[i + out_dim];
//                Hist.h[j] -= H[i + out_dim];
//            } else {
//            }
//        }
//    }


//// OFFICIAL GBDTMO
//    int j, ind, bin;
//    for (auto i : order) {
//        bin = maps[i] * out_dim;
//        ind = i * out_dim;
//        ++Hist.cnt[maps[i]];
//        for (j = 0; j < out_dim; ++j) {
//            Hist.g[bin] += G[ind];
//            Hist.h[bin++] += H[ind++];
//        }
//    }
//    // integration
////    int g_size= sizeof(Hist.g)/sizeof(Hist.g[0]);
////    int h_size= sizeof(Hist.h)/sizeof(Hist.h[0]);
////    printf("%d %d - \n",g_size, h_size);
//    ind = 0;
//    for (int i = 1; i < Hist.cnt.size(); ++i) {
//        Hist.cnt[i] += Hist.cnt[i - 1];
//        for (j = 0; j < out_dim; ++j) {
//            Hist.g[ind + out_dim] += Hist.g[ind];
//            Hist.h[ind + out_dim] += Hist.h[ind];
//            printf("%d %d - %d %d\n",bin, ind, i,ind + out_dim);
//            ++ind;
//        }
//    }
//}

