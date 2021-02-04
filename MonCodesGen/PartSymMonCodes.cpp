#include <set>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>

#include "PartSymMonCodes.h"

unsigned long long binomial(unsigned n, unsigned k) {
    if (k > n)
        return 0;
    unsigned long long res = 1;
    if (k > n - k) k = n - k;

    for (unsigned i = 0; i < k; ++i) {
        res = res * (n - i);
        res = res / (i + 1);
    }

    return res;
}

unsigned wt(unsigned i)
{
    unsigned w = 0;
    while (i > 0) {
        if (i % 2) w++;
        i >>= 1;
    }
    return w;
}

unsigned deg(unsigned i, unsigned m)
{
    return m - wt(i);
}


unsigned gcd(unsigned a, unsigned b)
{
    if (a == 0) return b;
    if (b == 0) return a;
    if (a >= b) return gcd(b, a % b);
    return gcd(a, b % a);
}


unsigned lcm(unsigned a, unsigned b)
{
    if (a == 0 && b == 0) return 0;
    return (a * b) / gcd(a, b);
}

std::vector<unsigned> vTmp;
std::vector<unsigned> vHist;
std::vector<std::set<unsigned>> vInfSymbols;
std::vector<std::set<unsigned>> vVars;
bool dfs_greedy2(unsigned target_k, unsigned cur_k, unsigned fin_deg, unsigned target_w, unsigned target_dec)
{
    if (cur_k == target_k)
        return true;
    auto&& worst_ind = *std::max_element(vTmp.begin(), vTmp.end(),
                                     [](unsigned a, unsigned b)
                                     {
                                         return vVars[a].size() < vVars[b].size();
                                     });
    auto&& worst = vVars[worst_ind];
    std::set<std::pair<unsigned, unsigned>> tbl;
    for (auto sym : worst)
    {
        if (wt(sym) != target_w)
            continue;
        bool f = true;
        for (auto&& j : vTmp)
        {
            if (vVars[j].size() == fin_deg && vVars[j].find(sym) != vVars[j].end())
            {
                f = false;
                break;
            }
        }
        if (f)
        {
            unsigned deg = 0;
            for (auto&& elem : vInfSymbols[sym])
                deg += vVars[elem].size();
            tbl.insert({ deg, sym });
       }
    }
    for (auto&& elem = tbl.rbegin(); elem != tbl.rend(); ++elem)
    {
        unsigned sym = elem->second;
        vHist[cur_k] = sym;
        auto&& inf_node = vInfSymbols[sym];

        unsigned dec = 0;
        for (auto&& ind : inf_node)
        {
            if (vVars[ind].find(sym) != vVars[ind].end())
            {
                vVars[ind].erase(sym);
                dec++;
            }
        }
        if (dec == target_dec)
        {
            if (dfs_greedy2(target_k, cur_k + 1, fin_deg, target_w, target_dec))
                return true;
        }
        for (auto&& ind : inf_node)
            vVars[ind].insert(sym);
    }
    return false;
}

unsigned c_gran(unsigned t, unsigned l)
{
    return lcm(t, l) / l;
}
unsigned p_gran(unsigned t, unsigned l)
{
    return lcm(t, l) / t;
}


void CPartSymMonCodeGen::Generate(unsigned m, unsigned k, unsigned t, unsigned r)
{
    unsigned n = 1u << m;
    std::set<unsigned> M;//initial tgenerating set
    for (unsigned j = 0; j < n; ++j)
        M.insert(j);

    //Create bipartite graph with target variables on one side and inf symbols on the other one
    vVars.resize(m);
    vInfSymbols.resize(n);
    for (unsigned i = m - t; i < m; ++i)
    {
        for (unsigned j = 0; j < n; ++j)
        {
            if (((j >> i) & 1) == 0)
            {//index j corresponds to monomial that contains variable x_i
                vInfSymbols[j].insert(i);
                vVars[i].insert(j);
            }
        }
    }
    unsigned cur_k = n;
    unsigned cur_l = t;
    unsigned proj_deg = cur_k / 2;
    unsigned cur_l_tot = (1ull << (m - t)) * binomial(t, cur_l);
    unsigned cur_p_gran = p_gran(t, cur_l);
    unsigned cur_c_gran = c_gran(t, cur_l);
    unsigned k_red = 0;
    for (unsigned i = r + 1; i <= m; ++i)
    {
        k_red += binomial(m, i);
        unsigned proj_red = 0;
        for (unsigned l = 1; l <= t; ++l)
        {
            if (l > i)
                break;
            unsigned l_imp = binomial(t, l) * binomial(m - t, i - l);
            unsigned c_g = c_gran(t, l);
            unsigned p_g = p_gran(t, l);
            proj_red += l_imp / c_g * p_g;
        }
        proj_deg -= proj_red;
    }
    if (cur_k - k_red < k)
        throw std::logic_error("Cannot achieve the target dimension, it can be at most " + std::to_string(cur_k - k_red));
    cur_k -= k_red;
    //removing all monomials with degree less than maxDeg
    for (unsigned j = 0; j < n; ++j)
    {
        if (deg(j, m) > r)
        {
            for (auto&& a : vVars)
                a.erase(j);
            M.erase(j);
        }
    }
    if (cur_l <= r)
    {
        unsigned val = (r + 1 >= cur_l) ? r + 1 - cur_l : 0;
        for (unsigned tmp_w = val; tmp_w <= m - t; ++tmp_w)
            cur_l_tot -= binomial(t, cur_l) * binomial(m - t, tmp_w);
    }
    else
        cur_l_tot = 0; 
    //all l-variate coverings
    while (cur_k - cur_l_tot >= k)
    {
        cur_k -= cur_l_tot;
        proj_deg -= cur_l_tot / cur_c_gran * cur_p_gran;
        for (unsigned j = 0; j < n; ++j)
        {
            if (vInfSymbols[j].size() == cur_l && M.find(j) != M.end())
            {
                M.erase(j);
                for (auto&& a : vVars)
                    a.erase(j);
            }
        }
        cur_l--;
        if (cur_l == 0)
            break;
        cur_l_tot = (1ull << (m - t)) * binomial(t, cur_l);
        if (cur_l <= r) {
            unsigned val = (r + 1 >= cur_l) ? r + 1 - cur_l : 0;
            for (unsigned tmp_w = val; tmp_w <= m - t; ++tmp_w)
                cur_l_tot -= binomial(t, cur_l) * binomial(m - t, tmp_w);
        }
        else
            cur_l_tot = 0;
        cur_p_gran = p_gran(t, cur_l);
        cur_c_gran = c_gran(t, cur_l);
    }
    if (cur_l == 0)
    {//dimensions of target derivatives are zero, so we can remove any monomials

        std::cout << "Removing " << (cur_k - k) << " extra monomials\n";
        for (unsigned j = 0; j < n && M.size() > k; ++j)
        {
            M.erase(M.begin());
            cur_k--;
        }
    } 
    else
    {
        unsigned cur_w = m - t;
        if (cur_l > r)
            cur_w = 0;
        else if (cur_l + cur_w > r)
            cur_w = r - cur_l;
        unsigned cur_w_tot = binomial(m - t, cur_w) * binomial(t, cur_l);
        //all l-variate coverings with monomials of degree (cur_w + cur_l)
        while (cur_k - cur_w_tot >= k) {
            cur_k -= cur_w_tot;
            proj_deg -= cur_w_tot / cur_c_gran * cur_p_gran;
            for (unsigned j = 0; j < n; ++j) {
                if (vInfSymbols[j].size() == cur_l && deg(j, m) == (cur_w + cur_l) && M.find(j) != M.end()) {
                    M.erase(j);
                    for (auto&& a : vVars)
                        a.erase(j);
                }
            }
            cur_w--;
            cur_w_tot = binomial(m - t, cur_w) * binomial(t, cur_l);
        }
        //finding set of monomials to remove while keeping t-symmetry as finding a biregular subgraph of biregular bipartite graph
        if (cur_k > k) {
            unsigned n_steps = (cur_k - k) / cur_c_gran;
            if ((cur_k - k) % cur_c_gran)
                n_steps++;
            unsigned diff = n_steps * cur_c_gran;
            cur_k -= diff;
            vTmp.resize(t);
            std::iota(vTmp.begin(), vTmp.end(), m - t);
            vHist.resize(diff);
            //max-flow algorithm always gives the correct solution, but dfs finds it too and sufficiently fast 
            if (!dfs_greedy2(diff, 0, proj_deg - n_steps * cur_p_gran, m - (cur_w + cur_l), cur_l))
                throw std::logic_error("Cannot construct code..terminating\n");
            for (unsigned i = 0; i < diff; ++i) {
                unsigned sym = vHist[i];
                M.erase(sym);
            }
        }
        //due to t-symmetry some values of k cannot be reached
        if (cur_k != k)
            std::cout << "k cannot be achieved, using closest value " << cur_k << std::endl;
    }
    //printing the code specification
    std::cerr << n << ' ' << cur_k << std::endl;
    for (unsigned j = 0; j < n; ++j)
    {
        if (M.find(j) == M.end())
            std::cerr << 1 << ' ' << j << std::endl;
    }
}


