#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <numeric> 
#include "omp.h"
using namespace std;


void get_data(vector<vector<int>> &arestas, map<int, int> &demandas, ifstream &file) {
    string line;
    int n_nos, n_rotas, id, demanda, origem, destino, peso;
    int i = 0;
    while (getline(file, line)) {
        istringstream iss(line);
        if (i == 0) {
            iss >> n_nos;
            arestas.resize(n_nos, vector<int>(n_nos, -1));
        }
        else if (i < n_nos){
            iss >> id >> demanda;
            demandas[id] = demanda;
        }
        else if (i == n_nos) {
            iss >> n_rotas;
        }
        else {
            iss >> origem >> destino >> peso;
            arestas[origem][destino] = peso;
        }
        i++;
    }
    file.close();
    return;
}


void combine(const vector<int>& locais, int start, int depth, vector<int>& current, vector<vector<int>>& result) {
    if (depth == 0) {
        do {
            vector<int> combination = {0};
            combination.insert(combination.end(), current.begin(), current.end());
            combination.push_back(0);
            result.push_back(combination);
        } while (next_permutation(current.begin(), current.end()));
        return;
    }

    for (int i = start; i < locais.size(); ++i) {
        current.push_back(locais[i]);
        combine(locais, i + 1, depth - 1, current, result);
        current.pop_back();
    }
}

void get_combinations(map<int, int>& demandas, vector<vector<int>>& combinacoes) {
    vector<int> locais;
    for (const auto& pair : demandas) {
        locais.push_back(pair.first);
    }

    for (int len = 1; len <= locais.size(); ++len) {
        vector<int> current;
        combine(locais, 0, len, current, combinacoes);
    }
}

void validate_combinations(map<int, int>& demandas, vector<vector<int>>& arestas, vector<vector<int>>& combinacoes, int max_capacity, vector<vector<int>>& combinacoes_validas) {
    
    #pragma omp parallel for
    for (const auto& combinacao : combinacoes) {
        bool is_valid = true;
        int current_capacity = 0;
        
        for (int i = 0; i < combinacao.size() - 1; ++i) {
            if (arestas[combinacao[i]][combinacao[i + 1]] == -1) {
                is_valid = false;
                break;
            }
            int inc = demandas[combinacao[i + 1]];
            if (current_capacity + inc > max_capacity) {
                is_valid = false;
                break;
            }
            current_capacity += inc;
        }
        
        if (is_valid) {
            #pragma omp critical
            combinacoes_validas.push_back(combinacao);
        }
    }
}


void get_best_route(vector<vector<int>>& combinacoes_validas, vector<vector<int>>& arestas, map<int, int>& demandas, vector<int>& resp, int& best_cost) {
    float best_value_per_cost = 0.0;
    int current_cost, current_value, src, dst;
    float ratio;

    #pragma omp parallel for private(current_cost, current_value, src, dst, ratio)
    for (int i = 0; i < combinacoes_validas.size(); ++i) {
        current_cost = 0;
        current_value = 0;

        for (int j = 0; j < combinacoes_validas[i].size() - 1; ++j) {
            src = combinacoes_validas[i][j];
            dst = combinacoes_validas[i][j + 1];

            current_value += demandas[dst];
            current_cost += arestas[src][dst];
        }

        ratio = static_cast<float>(current_value) / current_cost;

        #pragma omp critical
        {
            if (ratio > best_value_per_cost) {
                best_value_per_cost = ratio;
                best_cost = current_cost;
                resp = combinacoes_validas[i];
            }
        }
    }
}


int main(int argc, char* argv[]) {
    if (argc!= 2) {
        cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }

    ifstream file(argv[1]);
    if (!file.is_open()) {
        cerr << "Unable to open file: " << argv[1] << '\n';
        return 1; 
    }

    vector<vector<int>> arestas, combinacoes, combinacoes_validas;
    map<int, int> demandas;

    get_data(arestas, demandas, file);
    get_combinations(demandas, combinacoes);
    validate_combinations(demandas, arestas, combinacoes, 15, combinacoes_validas);

    /*cout << "Demandas:\n\n";
    for (auto &val : demandas) {
        std::cout << "(" << val.first << ", " << val.second << ") ";
    }
    cout << endl;
    
    cout << "Distancias:\n\n";
    for (const auto &row : arestas) {
        for (const auto &element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }

    cout << "Combinacoes possiveis:\n\n";
    for (const auto &row : combinacoes) {
        for (const auto &element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }

    cout << "Combinacoes validas:\n\n";

    for (const auto &row : combinacoes_validas) {
        for (const auto &element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }*/

    vector<int> resp;
    int best_cost;

    cout << "Melhor rota:\n\n";
    get_best_route(combinacoes_validas, arestas, demandas, resp, best_cost);

    for (const auto &element : resp) {
        std::cout << element << " ";
    }
    std::cout << std::endl;

    cout << "Custo: " << best_cost << endl;

    file.close();

    return 0;
}