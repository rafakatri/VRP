#include <mpi.h>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
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
            combinacoes_validas.push_back(combinacao);
        }
    }
}

void get_best_route(vector<vector<int>>& combinacoes_validas, vector<vector<int>>& arestas, map<int, int>& demandas, vector<int> &resp, int &best_cost, float &best_ratio) {
    best_ratio = 0;
    int current_cost, current_value, src, dst;
    float ratio;

    for (const auto& combinacao : combinacoes_validas) {
        current_cost = 0;
        current_value = 0;
        for (int j = 0; j < combinacao.size() - 1; j++) {
            src = combinacao[j];
            dst = combinacao[j + 1];

            current_value += demandas[dst];
            current_cost += arestas[src][dst];
        }
        ratio = static_cast<float>(current_value) / current_cost;
        if (ratio > best_ratio) {
            best_ratio = ratio;
            best_cost = current_cost;
            resp = combinacao;
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector<vector<int>> arestas, combinacoes, combinacoes_validas;
    map<int, int> demandas;
    vector<int> global_resp;
    int global_best_cost;

    if (rank == 0) {
        if (argc != 2) {
            cerr << "Usage: " << argv[0] << " <filename>\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        ifstream file(argv[1]);
        if (!file.is_open()) {
            cerr << "Unable to open file: " << argv[1] << '\n';
            MPI_Abort(MPI_COMM_WORLD, 1); 
        }

        get_data(arestas, demandas, file);
        get_combinations(demandas, combinacoes);
        validate_combinations(demandas, arestas, combinacoes, 15, combinacoes_validas);

        file.close();
    }

    // Broadcast arestas and demandas to all processes
    int n_nos = arestas.size();
    MPI_Bcast(&n_nos, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        arestas.resize(n_nos, vector<int>(n_nos));
    }

    for (int i = 0; i < n_nos; ++i) {
        MPI_Bcast(arestas[i].data(), n_nos, MPI_INT, 0, MPI_COMM_WORLD);
    }

    int demandas_size = demandas.size();
    MPI_Bcast(&demandas_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        for (int i = 0; i < demandas_size; ++i) {
            int id, demanda;
            MPI_Bcast(&id, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&demanda, 1, MPI_INT, 0, MPI_COMM_WORLD);
            demandas[id] = demanda;
        }
    } else {
        for (const auto& pair : demandas) {
            int id = pair.first;
            int demanda = pair.second;
            MPI_Bcast(&id, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&demanda, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
    }

    // Broadcast the number of valid combinations
    int num_combinacoes_validas = combinacoes_validas.size();
    MPI_Bcast(&num_combinacoes_validas, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Flatten the combinacoes_validas for scattering
    vector<int> flat_combinacoes_validas;
    vector<int> combinacoes_sizes;
    if (rank == 0) {
        for (const auto& combinacao : combinacoes_validas) {
            combinacoes_sizes.push_back(combinacao.size());
            flat_combinacoes_validas.insert(flat_combinacoes_validas.end(), combinacao.begin(), combinacao.end());
        }
    }

    // Broadcast the sizes of each combination
    int sizes_size = combinacoes_sizes.size();
    MPI_Bcast(&sizes_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        combinacoes_sizes.resize(sizes_size);
    }
    MPI_Bcast(combinacoes_sizes.data(), sizes_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast the flattened combinations
    int total_size = flat_combinacoes_validas.size();
    MPI_Bcast(&total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        flat_combinacoes_validas.resize(total_size);
    }
    MPI_Bcast(flat_combinacoes_validas.data(), total_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Determine the portion each process will handle
    int portion_size = (num_combinacoes_validas + size - 1) / size; // Ensure at least one combination per process
    int start = rank * portion_size;
    int end = min(start + portion_size, num_combinacoes_validas);

    // Extract local combinations from the flattened data
    vector<vector<int>> local_combinacoes_validas;
    int index = 0;
    for (int i = 0; i < start; ++i) {
        index += combinacoes_sizes[i];
    }
    for (int i = start; i < end; ++i) {
        vector<int> combinacao(combinacoes_sizes[i]);
        copy(flat_combinacoes_validas.begin() + index, flat_combinacoes_validas.begin() + index + combinacoes_sizes[i], combinacao.begin());
        local_combinacoes_validas.push_back(combinacao);
        index += combinacoes_sizes[i];
    }

    // Find the best local route
    vector<int> local_resp;
    int local_best_cost = numeric_limits<int>::max();
    float local_best_ratio = -1;
    if (local_combinacoes_validas.size() > 0) {
        get_best_route(local_combinacoes_validas, arestas, demandas, local_resp, local_best_cost, local_best_ratio);
    }

    // Gather the best costs and responses
    float global_best_ratio;
    MPI_Allreduce(&local_best_ratio, &global_best_ratio, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

    // Gather the best costs and responses
    if (local_best_ratio == global_best_ratio) {
        global_resp = local_resp;
        global_best_cost = local_best_cost;
    }

    // Broadcast the global best response
    int resp_size = global_resp.size();
    MPI_Bcast(&resp_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        global_resp.resize(resp_size);
    }
    MPI_Bcast(global_resp.data(), resp_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast the global best cost
    MPI_Bcast(&global_best_cost, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Print the result in rank 0
    if (rank == 0) {
        cout << "Best route: ";
        for (const int& node : global_resp) {
            cout << node << " ";
        }
        cout << "\nBest cost: " << global_best_cost << endl;
    }

    MPI_Finalize();
    return 0;
}
