#include "../TYPES/WMMStructs.h"
#include <vector>
#include <map>
#include <iostream>

using namespace std;

namespace wmm {

    double GetInterpValue(wmm::Wmm_<double> &wave, wmm::NodeD &dp, wmm::NodeD &dd, wmm::NodeD &dn, wmm::NodeD &f0,
                          wmm::NodeD &f1, wmm::NodeD &fn, double epsilon, int side) {

        double ft, value, y0 = wave.v[0], y1 = wave.v[side + 1];

        ft = wave.fm[2 * side + 2] * epsilon * epsilon + wave.fm[2 * side + 1] * epsilon + wave.fm[0];
        double norm = wmm::norm(dn - dp - epsilon * dd);
        value = wave.m[2 * side + 2] * epsilon * epsilon + wave.m[2 * side + 1] * epsilon + wave.m[0] +
                norm * (ft + wmm::norm(fn)) / 2.0;
        if (value < y0) {
            double v0 = (1.0 - epsilon) * y0;
            double v1 = epsilon * y1;
            value = v0 + v1 + norm * ((1.0 - epsilon) * wmm::norm(f0) + epsilon * wmm::norm(f1) + wmm::norm(fn)) / 2.0;
        }

        return value;
    }

    double GetEpsilonGradient(wmm::NodeD &dd, wmm::NodeD &dp, wmm::NodeD &dn, wmm::NodeD &fn) {
        double epsilon;
        double A = -dd.y, B = dd.x, C = dd.y * dp.x - dd.x * dp.y;
        double den = B * fn.y;
        double t = (A * dn.x + B * dn.y + C) / den;

        wmm::NodeD x(dn.y - t * fn.y, dn.x);

        if (fabs(dd.x) > 0.0 && fabs(den) > 0.0) {
            epsilon = (x.x - dp.x) / dd.x;
        } else if (fabs(dd.y) > 0.0 && fabs(den) > 0.0) {
            epsilon = (x.y - dp.y) / dd.y;
        } else if (fabs(den) == 0.0 && wmm::norm(dd) > 0.0) {
            double dist = fabs(A * dn.x + B * dn.y + C) / sqrt(A * A + B * B);
            epsilon = (wmm::norm(dn - dp) - dist) / (fabs(dd.x) + fabs(dd.y));
        } else {
            return 0.0;
        }

        if (epsilon < 0.0)
            epsilon = 0.0;
        else if (epsilon > 1.0)
            epsilon = 1.0;
        return epsilon;
    }

    double GetInterpValueSimple(wmm::Wmm_<double> &wave, wmm::NodeD &dp, wmm::NodeD &dd, wmm::NodeD &dn, double f0,
                                double f1, double fn, double epsilon, int side) {

        double ft, value, y0 = wave.v[0], y1 = wave.v[side + 1];

        ft = wave.fm[2 * side + 2] * epsilon * epsilon + wave.fm[2 * side + 1] * epsilon + wave.fm[0];
        double norm = wmm::norm(dn - dp - epsilon * dd);
        value = wave.m[2 * side + 2] * epsilon * epsilon + wave.m[2 * side + 1] * epsilon + wave.m[0] +
                norm * (ft + fn) / 2.0;
        if (value < y0) {
            double v0 = (1.0 - epsilon) * y0;
            double v1 = epsilon * y1;
            value = v0 + v1 + norm * ((1.0 - epsilon) * f0 + epsilon * f1 + fn) / 2.0;
        }

        return value;
    }

    double GetEpsilonGradientSimple(wmm::NodeD &dd, wmm::NodeD &dp, wmm::NodeD &dn, double fn) {
        double epsilon;
        double A = -dd.y, B = dd.x, C = dd.y * dp.x - dd.x * dp.y;
        double den = A * fn + B * fn;
        double t = (A * dn.x + B * dn.y + C) / den;

        wmm::NodeD x(dn.y - t * fn, dn.x - t * fn);

        if (fabs(dd.x) > 0.0 && fabs(den) > 0.0) {
            epsilon = (x.x - dp.x) / dd.x;
        } else if (fabs(dd.y) > 0.0 && fabs(den) > 0.0) {
            epsilon = (x.y - dp.y) / dd.y;
        } else if (fabs(den) == 0.0 && wmm::norm(dd) > 0.0) {
            double dist = fabs(A * dn.x + B * dn.y + C) / sqrt(A * A + B * B);
            epsilon = (wmm::norm(dn - dp) - dist) / (fabs(dd.x) + fabs(dd.y));
        } else {
            return 0.0;
        }

        if (epsilon < 0.0)
            epsilon = 0.0;
        else if (epsilon > 1.0)
            epsilon = 1.0;
        return epsilon;
    }

    double GetVal2DSimple(wmm::Grid &image, wmm::Grid &u_surface, wmm::Wmm_<double> &wave, wmm::Node &neigh, wmm::NodeD &h) {
        wmm::NodeD f0(image.at(wave.p, 0), 0), fn(image.at(neigh, 0), 0);
        double y0 = wave.v[0];

        if (isinf(wmm::norm(f0)) || isnan(wmm::norm(f0)))
            f0 = fn;

        double val;
        if (wave.dir < 0) {
            wmm::NodeD diff(h.y * (neigh.y - wave.p.y), h.x * (neigh.x - wave.p.x));
            val = y0 + wmm::norm(diff) * (wmm::norm(f0) + wmm::norm(fn)) / 2.0;
        } else {

            wmm::Node p(wave.p.y + wmm::yarray[(wave.dir + 1) % 8] - wmm::yarray[wave.dir],
                        wave.p.x + wmm::xarray[(wave.dir + 1) % 8] - wmm::xarray[wave.dir]);
            double res1 = wmm::MAX_VAL;
            wmm::NodeD dp(h.y * wave.p.y, h.x * wave.p.x), dn(h.y * neigh.y, h.x * neigh.x);

            if (u_surface.contains(p)) {
                wmm::NodeD dd(h.y * (wmm::yarray[(wave.dir + 1) % 8] - wmm::yarray[wave.dir]),
                              h.x * (wmm::xarray[(wave.dir + 1) % 8] - wmm::xarray[wave.dir]));

                wmm::NodeD f1(image.at(p, 0), 0);
                if (isinf(wmm::norm(f1)) || isnan(wmm::norm(f1)))
                    f1 = fn;

                double epsilon = GetEpsilonGradient(dd, dp, dn, fn);

                res1 = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, epsilon, wmm::S_RIGHT);
            }

            p = wmm::Node(wave.p.y + wmm::yarray[(wave.dir + 7) % 8] - wmm::yarray[wave.dir],
                          wave.p.x + wmm::xarray[(wave.dir + 7) % 8] - wmm::xarray[wave.dir]);
            double res2 = wmm::MAX_VAL;

            if (u_surface.contains(p)) {
                wmm::NodeD dd(h.y * (wmm::yarray[(wave.dir + 7) % 8] - wmm::yarray[wave.dir]),
                              h.x * (wmm::xarray[(wave.dir + 7) % 8] - wmm::xarray[wave.dir]));

                wmm::NodeD f1(image.at(p, 0), 0);
                if (isinf(wmm::norm(f1)) || isnan(wmm::norm(f1)))
                    f1 = fn;

                double epsilon = GetEpsilonGradient(dd, dp, dn, fn);

                res2 = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, epsilon, wmm::S_LEFT);

            }

            val = std::min(res1, res2);

        }

        return val;
    }

    void setCoeffs2D(const double *y, double *m, int pos) {
        m[0] = y[pos];
        if (pos % 2 == 0) {
            m[2] = (y[(pos + 2) % 8] + y[pos] - 2.0 * y[(pos + 1) % 8]) / 2.0;
            m[1] = y[(pos + 1) % 8] - y[pos] - m[2];
            m[4] = (y[(pos + 6) % 8] + y[pos] - 2.0 * y[(pos + 7) % 8]) / 2.0;
            m[3] = y[(pos + 7) % 8] - y[pos] - m[4];
        } else {
            m[2] = m[4] = (y[(pos + 1) % 8] + y[(pos + 7) % 8] - 2.0 * y[pos]) / 2.0;
            m[1] = y[pos] - y[(pos + 7) % 8] + m[2];
            m[3] = y[pos] - y[(pos + 1) % 8] + m[4];
        }

    }

    wmm::Grid WmmIsoSurface2D(wmm::Grid &image, std::vector<wmm::Node> &initials, wmm::NodeD &h) {
        bool isnewpos[8];
        double valcenter[8];
        double imcenter[8];

        wmm::Grid u_surface = wmm::Grid(wmm::MAX_VAL, image.rows, image.cols, 1);
        wmm::Grid_<unsigned char> state = wmm::Grid_<unsigned char>(image.rows, image.cols);

        std::multimap<double, wmm::Wmm_<double> > trial_set;
        std::map<int, typename std::multimap<double, wmm::Wmm_<double> >::iterator> mapa_trial;

        typename std::multimap<double, wmm::Wmm_<double> >::iterator trial_set_it;
        typename std::map<int, typename std::multimap<double, wmm::Wmm_<double> >::iterator>::iterator mapa_trial_it;
        std::pair<double, wmm::Wmm_<double> > pr_trial;
        std::pair<int, typename std::multimap<double, wmm::Wmm_<double> >::iterator> pr_mapa;

        int key, i;
        wmm::Wmm_<double> winner, new_w;
        wmm::Node neigh;

        // Initialization
        for (i = 0; i < initials.size(); i++) {
            key = initials[i].y * u_surface.cols + initials[i].x;
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                winner.dir = -1;
                winner.v[0] = 0.0;
                winner.p = initials[i];
                state.at(initials[i]) = wmm::P_TRIAL;
                pr_trial = std::pair<double, wmm::Wmm_<double> >(0.0, winner);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, typename std::multimap<double, wmm::Wmm_<double> >::iterator>(key,
                                                                                                       trial_set_it);
                mapa_trial.insert(pr_mapa);
            }
        }

        while (!trial_set.empty()) {
            trial_set_it = trial_set.begin();
            key = trial_set_it->second.p.y * u_surface.cols + trial_set_it->second.p.x;
            mapa_trial_it = mapa_trial.find(key);

            if (mapa_trial_it == mapa_trial.end()) {
                printf("ERROR: bad map alloc");
                exit(-1);
            }

            if (mapa_trial_it->second != trial_set_it) {
                printf("ERROR: bad trial/map alloc");
                exit(-1);
            }

            winner = trial_set_it->second;

            trial_set.erase(trial_set_it);
            mapa_trial.erase(mapa_trial_it);

            state.at(winner.p) = wmm::P_ALIVE;

            // Neighbour temptative value computation
            for (i = 0; i < 8; i++) {
                neigh = winner.p + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
                isnewpos[i] = false;
                valcenter[i] = u_surface.contains(neigh) ? u_surface.at(neigh) : wmm::MAX_VAL;
                imcenter[i] = u_surface.contains(neigh) ? image.at(neigh, 0) : wmm::MAX_VAL;
                if (u_surface.contains(neigh) && state.at(neigh) != wmm::P_ALIVE) {
                    double val_neigh = GetVal2DSimple(image, u_surface, winner, neigh, h);
                    if (val_neigh < valcenter[i]) {
                        valcenter[i] = val_neigh;
                        isnewpos[i] = true;
                    }
                }
            }

            // Update
            for (i = 0; i < 8; i++) {
                if (isnewpos[i]) {
                    neigh = winner.p + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
                    key = neigh.y * u_surface.cols + neigh.x;
                    if (state.at(neigh) == wmm::P_TRIAL) {
                        mapa_trial_it = mapa_trial.find(key);
                        trial_set.erase(mapa_trial_it->second);
                        mapa_trial.erase(mapa_trial_it);
                    } else {
                        state.at(neigh) = wmm::P_TRIAL;
                    }
                    new_w.p = neigh;
                    new_w.dir = i;

                    new_w.v[0] = valcenter[i];
                    new_w.v[1] = valcenter[(i + 1) % 8];
                    new_w.v[2] = valcenter[(i + 7) % 8];

                    setCoeffs2D(valcenter, new_w.m, i);
                    setCoeffs2D(imcenter, new_w.fm, i);

                    pr_trial = std::pair<double, wmm::Wmm_<double> >(valcenter[i], new_w);
                    trial_set_it = trial_set.insert(pr_trial);
                    pr_mapa = std::pair<int, typename std::multimap<double, wmm::Wmm_<double> >::iterator>(
                            key, trial_set_it);
                    mapa_trial.insert(pr_mapa);

                    u_surface.at(new_w.p) = valcenter[i];
                }
            }
        }

        free(state.data);
        return u_surface;
    }
}

#include <iomanip>

int main() {
    int rows = 8;
    int cols = 8;

    /* Create output arrays */
    double *out = (double *) malloc(sizeof(double) * rows * cols);
    double surface[] = {
            2, 2, 2, 2, 2, 2, 2, 4,
            2, 2, 2, 2, 2, 2, 2, 4,
            2, 2, 2, 2, 2, 2, 2, 4,
            2, 2, 2, 2, 2, 2, 3, 3,
            2, 2, 2, 2, 2, 2, 3, 3,
            2, 2, 2, 2, 2, 2, 3, 3,
            2, 2, 2, 2, 2, 2, 3, 3,
            2, 2, 2, 2, 2, 2, 3, 3,

            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
    };

    wmm::Grid imagen(surface, rows, cols);

    std::vector<wmm::Node> initials;
    //initials.emplace_back(1, 2);
    initials.emplace_back(1, 1);
    //initials.emplace_back(0, 0);
    wmm::NodeD hs = wmm::NodeD(1, 1);
    wmm::Grid out_surface(out, rows, cols);
    out_surface = wmm::WmmIsoSurface2D(imagen, initials, hs);

    for (int i = 0; i < rows * cols; i++) {
        if (i % cols == 0) {
            cout << "\n";
        }
        cout << surface[i] << " ";
    }

    cout << "\n-----------" << std::fixed << std::setprecision(2);

    for (int i = 0; i < rows * cols; i++) {
        if (i % cols == 0) {
            cout << "\n";
        }
        if (out_surface.data[i] > 1000) {
            cout << "XXXX" << " ";
        } else {
            cout << out_surface.data[i] << "   ";
        }
    }
}
