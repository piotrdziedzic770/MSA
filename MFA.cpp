#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>

using namespace std;

const int MATCH = 1; 
const int MISMATCH = -1;
const int GAP_OPEN = -2; 
const int GAP_EXT = -1;

int evaluatePair(const string& s1, const string& s2) {
    int score = 0;
    bool in_gap1 = false;
    bool in_gap2 = false;

    for (size_t i = 0; i < s1.length(); ++i) {
        char c1 = s1[i];
        char c2 = s2[i];

        if (c1 == '-' && c2 == '-') {
            continue;
        }
        else if (c1 == '-') {
            score += in_gap1 ? GAP_EXT : GAP_OPEN;
            in_gap1 = true;
            in_gap2 = false;
        }
        else if (c2 == '-') {
            score += in_gap2 ? GAP_EXT : GAP_OPEN;
            in_gap2 = true;
            in_gap1 = false;
        }
        else {
            score += (c1 == c2) ? MATCH : MISMATCH;
            in_gap1 = false;
            in_gap2 = false;
        }
    }
    return score;
}

int calculateSP(const vector<string>& alignment) {
    int total_score = 0;
    int n = alignment.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            total_score += evaluatePair(alignment[i], alignment[j]);
        }
    }
    return total_score;
}

void mutate(vector<string>& alignment, mt19937& gen) {
    uniform_int_distribution<> row_dist(0, alignment.size() - 1);
    int row_idx = row_dist(gen);
    string& seq = alignment[row_idx];

    vector<int> gap_positions;
    for (int i = 0; i < seq.length(); ++i) {
        if (seq[i] == '-') {
            gap_positions.push_back(i);
        }
    }

    if (gap_positions.empty()) {
        return;
    }

    uniform_int_distribution<> gap_dist(0, gap_positions.size() - 1);
    int gap_to_remove_idx = gap_positions[gap_dist(gen)];

    seq.erase(gap_to_remove_idx, 1);

    uniform_int_distribution<> new_pos_dist(0, seq.length());
    int new_gap_pos = new_pos_dist(gen);

    seq.insert(new_gap_pos, 1, '-');
}

vector<string> initializeAlignment(const vector<string>& sequences, int extra_buffer = 10) {
    size_t max_len = 0;
    for (const auto& seq : sequences) {
        max_len = max(max_len, seq.length());
    }

    vector<string> alignment = sequences;
    size_t target_len = max_len + extra_buffer;

    for (auto& seq : alignment) {
        while (seq.length() < target_len) {
            seq += '-';
        }
    }
    return alignment;
}

void cleanupColumns(vector<string>& alignment) {
    if (alignment.empty()) return;
    int len = alignment[0].length();
    for (int col = len - 1; col >= 0; --col) {
        bool all_gaps = true;
        for (const auto& seq : alignment) {
            if (seq[col] != '-') {
                all_gaps = false;
                break;
            }
        }
        if (all_gaps) {
            for (auto& seq : alignment) {
                seq.erase(col, 1);
            }
        }
    }
}

vector<string> simulatedAnnealing(const vector<string>& sequences) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> prob_dist(0.0, 1.0);

    double temp = 10000.0;
    double min_temp = 0.01;
    double cooling_rate = 0.99;
    int iter_per_temp = 200;

    vector<string> current_solution = initializeAlignment(sequences);
    int current_score = calculateSP(current_solution);

    vector<string> best_solution = current_solution;
    int best_score = current_score;

    while (temp > min_temp) {
        for (int i = 0; i < iter_per_temp; ++i) {
            vector<string> new_solution = current_solution;

            mutate(new_solution, gen);
            int new_score = calculateSP(new_solution);

            double delta = current_score - new_score;

            if (delta < 0 || exp(-delta / temp) > prob_dist(gen)) {
                current_solution = new_solution;
                current_score = new_score;

                if (current_score > best_score) {
                    best_solution = current_solution;
                    best_score = current_score;
                }
            }
        }
        temp *= cooling_rate;
    }

    cleanupColumns(best_solution);
    return best_solution;
}

void printAlignment(const vector<string>& alignment, int score) {
    cout << "Osiagniety wynik sum-of-pairs: " << score << "\n";
    cout << "Dopasowanie:\n";
    for (size_t i = 0; i < alignment.size(); ++i) {
        cout << "w" << (i + 1) << " : " << alignment[i] << "\n";
    }
}

int main() {
    vector<string> input_sequences = {
        "ATGGCCCTGTGGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        "ATGGCCCTGTGGATCGTACGATCGTAGCTACTAGCTAGCTAGCTAGCTAGCTAGC",
        "ATGGCCCTGTGGATCGTACGATCGTAGCTAGCTAGCTTGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        "GGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        "ATGGCCCTGTGGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        "ATGGCCCTGTGGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        "ATGGCCCTGTGGATCGTACGATCGAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        "ATGGCCCTGTGGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
    };

    vector<string> initial_sol = initializeAlignment(input_sequences);
    cleanupColumns(initial_sol);
    cout << "\nStan poczatkowy\n";
    printAlignment(initial_sol, calculateSP(initial_sol));

    vector<string> best_alignment = simulatedAnnealing(input_sequences);

    cout << "\nStan po optymalizacji\n";
    printAlignment(best_alignment, calculateSP(best_alignment));

    return 0;
}