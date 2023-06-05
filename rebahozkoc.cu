
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include <iostream>
#include <string>

#include <fstream>
#include <sstream>
#include <vector>

#include <array>
#include <omp.h>
#include <stack>


using namespace std;


enum direction {
    d_down,
    d_right,
    none
};

#define COORD std::pair<int, int>

//#define DEBUG

int iter = 0;

/// Auxiliary functions

void display_arr(int* arr, int n) {
    cout << "arr: ";
    for (int i = 0; i < n; i++) {
        cout << arr[i] << " ";
    }
    cout << endl;
}

void print_coords(COORD start, COORD end) {

    cout << "Start:" << start.first << "," << start.second << endl;
    cout << "End:" << end.first << "," << end.second << endl;
}

int find_length(COORD start, COORD end, direction dir) {
    if (dir == d_down)
        return end.first - start.first;
    if (dir == d_right)
        return end.second - start.second;
    return -1;
}

void convert_sol(int** mat, int**& sol_mat, int m, int n) {

    sol_mat = new int* [m]; // Rows
    for (int i = 0; i < m; i++) {
        sol_mat[i] = new int[n]; // Cols
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (mat[i][j] == -2)
                sol_mat[i][j] = -2; // Empty value cell
            else
                sol_mat[i][j] = -1; // Hint or empty cell
        }
    }
}

void print_one_matrix(int** matrix, int m, int n) {
    std::cout << "Matrix: " << std::endl;
    for (int i = 0; i < m; i++) { // rows
        for (int j = 0; j < n; j++) { // cols
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << "\n";
    }
}

void sol_to_file(int** mat, int** sol_mat, int m, int n, string fname) {
    ofstream to_write(fname);

    to_write << m << " " << n << "\n";

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (mat[i][j] != -2)
                to_write << mat[i][j] << " ";
            else
                to_write << sol_mat[i][j] << " ";
        }
        to_write << "\n";
    }

    to_write.close();
}

void read_matrix(int**& matrix, std::ifstream& afile, int m, int n) {

    matrix = new int* [m]; // rows

    for (int i = 0; i < m; i++) {
        matrix[i] = new int[n]; // cols
    }

    int val;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            afile >> val;
            matrix[i][j] = val;
        }
    }
}

/// Auxiliary functions

struct sum {
    COORD start;
    COORD end;

    int hint;
    int dir;
    int length;
    int posMin;
    int posMax;

    void print_sum() {
        cout << "############################" << endl;
        cout << "Creating sum with: " << endl;
        print_coords(start, end);
        cout << "Hint: " << hint << endl;
        cout << "Direction: " << dir << endl;
        cout << "Length: " << length << endl;
        cout << endl;
        cout << "############################" << endl;
    }

    sum(COORD _start, COORD _end, int _hint, direction _dir) : start(_start), end(_end), hint(_hint), dir(_dir) {
        length = find_length(_start, _end, _dir);
        // This is equal to hint - sum of numbers 9 + 8 + 7
        posMin = hint - 45 + ((8 - length) * (9 - length)) / 2;

        posMax = hint - (length * (length - 1)) / 2;

#ifdef DEBUG
        cout << "############################" << endl;
        cout << "Creating sum with: " << endl;
        print_coords(start, end);
        cout << "Hint: " << hint << endl;
        cout << "Direction: " << dir << endl;
        cout << "Length: " << length << endl;
        cout << "############################" << endl;
#endif
    }
};

COORD find_end(int** matrix, int m, int n, int i, int j, direction dir) { // 0 down 1 right

    if (dir == d_right) {
        for (int jj = j + 1; jj < n; jj++) {
            if (matrix[i][jj] != -2 || jj == n - 1) {
                if (matrix[i][jj] == -2 && jj == n - 1)
                    jj++;
                COORD END = COORD(i, jj);
                return END;
            }
        }
    }

    if (dir == d_down) {
        for (int ii = i + 1; ii < m; ii++) {
            if (matrix[ii][j] != -2 || ii == m - 1) {
                if (matrix[ii][j] == -2 && ii == m - 1)
                    ii++;
                COORD END = COORD(ii, j);
                return END;
            }
        }
    }

    cout << "ERROR: Find end is called with faulty parameters." << endl;
    return COORD(0, 0);
}

vector<sum> get_sums(int** matrix, int m, int n) {

    vector<sum> sums;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int val = matrix[i][j];
            if (val != -1 && val != -2) {
                int hint = val;
                hint = hint / 10;
                // right sum
                if ((hint % 100) == 0) {
                    hint = (int)(hint / 100);
                    COORD START = COORD(i, j + 1);
                    COORD END = find_end(matrix, m, n, i, j, d_right);
                    sum _sum = sum(START, END, hint, d_right);
                    sums.push_back(_sum);
                }

                else {
                    int div = (int)(hint / 100);
                    int rem = (int)(hint % 100);
                    // down sum
                    if (div == 0 && rem != 0) {
                        COORD START = COORD(i + 1, j);
                        COORD END = find_end(matrix, m, n, i, j, d_down);
                        sum _sum = sum(START, END, rem, d_down);
                        sums.push_back(_sum);
                    }
                    // combined sum
                    if (div != 0 && rem != 0) {
                        COORD START1 = COORD(i + 1, j);
                        COORD START2 = COORD(i, j + 1);
                        COORD END1 = find_end(matrix, m, n, i, j, d_down);
                        COORD END2 = find_end(matrix, m, n, i, j, d_right);
                        sum _sum1 = sum(START1, END1, rem, d_down);
                        sum _sum2 = sum(START2, END2, div, d_right);
                        sums.push_back(_sum1);
                        sums.push_back(_sum2);
                    }
                }
            }
        }
    }
    return sums;
}



void print_flattened(int* h_sum_starts_x, int* h_sum_starts_y, int* h_sum_ends_x, int* h_sum_ends_y, int* h_sum_hints, int* h_sum_lengths, int* h_sum_dirs, int no_sums) {

    cout << "###h_sum_starts_x: " << endl;
    for (int i = 0; i < no_sums; i++) {
        cout << h_sum_starts_x[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_starts_y: " << endl;
    for (int i = 0; i < no_sums; i++) {
        cout << h_sum_starts_y[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_ends_x: " << endl;
    for (int i = 0; i < no_sums; i++) {
        cout << h_sum_ends_x[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_ends_y: " << endl;
    for (int i = 0; i < no_sums; i++) {
        cout << h_sum_ends_y[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_hints: " << endl;
    for (int i = 0; i < no_sums; i++) {
        cout << h_sum_hints[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_lengths: " << endl;
    for (int i = 0; i < no_sums; i++) {
        cout << h_sum_lengths[i] << " ";
    }
    cout << endl;

    cout << "###h_sum_dirs: " << endl;
    for (int i = 0; i < no_sums; i++) {
        cout << h_sum_dirs[i] << " ";
    }
    cout << endl;
}

void flatten_sol_mat(int** sol_mat, int* h_sol_mat, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            h_sol_mat[i * n + j] = sol_mat[i][j];
        }
    }
}

void print_flattened_matrix(int* h_sol_mat, int m, int n) {

    cout << "###Flattened matrix: " << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << h_sol_mat[i * n + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}


__global__
void print_flattened_matrix_device(int* d_sol_mat, int m, int n) {

    printf("###Flattened matrix: %d, %d\n", m, n);
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", d_sol_mat[i * n + j]);
        }
        printf("\n");
    }
    printf("\n");
}

__device__
bool checkSumPartial(int* mat, int m, int n, int sum_start_x, int sum_start_y, int sum_end_x, int sum_end_y, int sum_hint, int sum_length, int sum_dir) {
    int sum = 0;
    if (sum_dir == 1) { // right
        for (int j = sum_start_y; j < sum_end_y; j++) {
            if (mat[sum_start_x * n + j] != -2) {
                sum += mat[sum_start_x *n + j];
            }

        }
    }
    else { // down
        for (int i = sum_start_x; i < sum_end_x; i++) {
            if (mat[i*n + sum_start_y] != -2) {
                sum += mat[i * n + sum_start_y];
            }
        }
    }

    if (sum > sum_hint) {
        return false;
    }
    return true;
}
__device__
bool checkSumComplete(int* mat, int m, int n, int sum_start_x, int sum_start_y, int sum_end_x, int sum_end_y, int sum_hint, int sum_length, int sum_dir) {
    int sum = 0;
    if (sum_dir == 1) { // right
        for (int j = sum_start_y; j < sum_end_y; j++) {
            if (mat[sum_start_x * n + j] != -2) {
                sum += mat[sum_start_x * n + j];
            }

        }
    }
    else { // down
        for (int i = sum_start_x; i < sum_end_x; i++) {
            if (mat[i* n + sum_start_y] != -2) {
                sum += mat[i * n + sum_start_y];
            }
        }
    }
    if (sum != sum_hint) {
        return false;
    }
    return true;
}

__device__
bool areElementsUniqueExceptEmpties(int* mat, int m, int n, int sum_start_x, int sum_start_y, int sum_end_x, int sum_end_y, int sum_hint, int sum_length, int sum_dir) {
    if (sum_dir == 1) { // right
        for (int j = sum_start_y; j < sum_end_y; j++) {
            if (mat[sum_start_x * n + j] != -2) {
                for (int k = j + 1; k < sum_end_y; k++) {
                    if (mat[sum_start_x * n +j] == mat[sum_start_x * n + k] && mat[sum_start_x * n + k] != -2) {
                        // If the same non-empty element is found, it means it's a duplicate, and we can return false
                        return false;
                    }
                }
            }
        }
    }
    else { // down
        for (int i = sum_start_x; i < sum_end_x; i++) {
            if (mat[i* n + sum_start_y] != -2) {
                for (int k = i + 1; k < sum_end_x; k++) {
                    if (mat[i* n + sum_start_y] == mat[k* n + sum_start_y] && mat[k * n + sum_start_y] != -2) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

__device__
bool areElementsUnique(int* mat, int m, int n, int sum_start_x, int sum_start_y, int sum_end_x, int sum_end_y, int sum_hint, int sum_length, int sum_dir) {
    if (sum_dir == 1) { // right
        for (int j = sum_start_y; j < sum_end_y; j++) {
            for (int k = j + 1; k < sum_end_y; k++) {
                if (mat[sum_start_x * n + j] == mat[sum_start_x * n + k]) {
                    // If the same element is found, it means it's a duplicate, and we can return false
                    return false;
                }
            }
        }
    }
    else { // down
        for (int i = sum_start_x; i < sum_end_x; i++) {
            for (int k = i + 1; k < sum_end_x; k++) {
                if (mat[i* n + sum_start_y] == mat[k* n + sum_start_y]) {
                    return false;
                }
            }
        }
    }
    return true;
}
__device__
bool isArrFull(int* mat, int m, int n, int sum_start_x, int sum_start_y, int sum_end_x, int sum_end_y, int sum_hint, int sum_length, int sum_dir) {
    if (sum_dir == 1) { // right
        for (int j = sum_start_y; j < sum_end_y; j++) {
            if (mat[sum_start_x * n + j] == -2)
                return false;
        }
    }
    else { // down
        for (int i = sum_start_x; i < sum_end_x; i++) {
            if (mat[i* n + sum_start_y] == -2)
                return false;
        }
    }
    return true;
}

__device__
bool isACandidate(int val, int posMin, int posMax){
    if (val < posMin) {
        return false;
    }
    if (val > posMax) {
        return false;
    }
    return true;
}


__device__
bool isLastCell(int curr_i, int curr_j, int sum_start_x, int sum_start_y, int sum_end_x, int sum_end_y, int sum_hint, int sum_length, int sum_dir) {
    if (sum_dir == 1) { // to right
        if (sum_start_x == curr_i && sum_end_y - 1 == curr_j) {
            return true;
        }
    }
    else {
        if (sum_start_y == curr_j && sum_end_x - 1 == curr_i) {
            {
                return true;
            }
        }
    }
    return false;
}

__device__
bool fullCheck(int* mat, int curr_i, int curr_j, int val, int m, int n, int sum_start_x, int sum_start_y, int sum_end_x, int sum_end_y, int sum_hint, int sum_length, int sum_dir, int sum_min, int sum_max) {
    int i = curr_i;
    int j = curr_j;
    
    if (sum_dir == 1) { // right
        // current coordinate is not in the sum
        if (!((i == sum_start_x) && (j >= sum_start_x) && (j < sum_end_y))) {
            return true;
        }
    }
    else { // down
        if (!(j == sum_start_y && i >= sum_start_x && i < sum_end_x)) {
            return true;
        }
    }
    if (!isACandidate(val, sum_min, sum_max)) {
        return false;
    }
    if (!areElementsUniqueExceptEmpties(mat, m, n, sum_start_x, sum_start_y, sum_end_x, sum_end_y, sum_hint, sum_length, sum_dir)) {
        return false;
    }
    if (isArrFull(mat, m, n, sum_start_x, sum_start_y, sum_end_x, sum_end_y, sum_hint, sum_length, sum_dir) || isLastCell(curr_i, curr_j, sum_start_x, sum_start_y, sum_end_x, sum_end_y, sum_hint, sum_length, sum_dir)) {
        return checkSumComplete(mat, m, n, sum_start_x, sum_start_y, sum_end_x, sum_end_y, sum_hint, sum_length, sum_dir) && areElementsUnique(mat, m, n, sum_start_x, sum_start_y, sum_end_x, sum_end_y, sum_hint, sum_length, sum_dir);
    }
    else {
        return checkSumPartial(mat, m, n, sum_start_x, sum_start_y, sum_end_x, sum_end_y, sum_hint, sum_length, sum_dir);
    }
}


int mat_iter_get_next_i(int* mat, int m, int n, int curr_i, int curr_j) {
    do {
        if (curr_j + 1 < n) {
            curr_j += 1;
        }
        else {
            curr_i += 1;
            curr_j = 0;
        }
    } while (curr_i != m && (mat[curr_i * n + curr_j] != -2));
    if (curr_i == m) {
        curr_i = -999; // end of the iteration
        curr_j = -999;

    }
    return curr_i;
}


int mat_iter_get_next_j(int* mat, int m, int n, int curr_i, int curr_j) {
    do {
        if (curr_j + 1 < n) {
            curr_j += 1;
        }
        else {
            curr_i += 1;
            curr_j = 0;

        }
    } while (curr_i != m && (mat[curr_i * n + curr_j] != -2));
    if (curr_i == m) {
        curr_i = -999; // end of the iteration
        curr_j = -999;

    }
    return curr_j;
}

int mat_iter_init_i(int* mat, int m, int n) {
    int curr_i = 0;
    int curr_j = 0;
    if (mat[curr_i*n + curr_j] != -2) {
        curr_i = mat_iter_get_next_i(mat, m, n, curr_i, curr_j);
        curr_j = mat_iter_get_next_j(mat, m, n, curr_i, curr_j);
    }
    return curr_i;
}

int mat_iter_init_j(int* mat, int m, int n) {
    int curr_i = 0;
    int curr_j = 0;
    if (mat[curr_i * n + curr_j] != -2) {
        curr_i = mat_iter_get_next_i(mat, m, n, curr_i, curr_j);
        curr_j = mat_iter_get_next_j(mat, m, n, curr_i, curr_j);
    }
    return curr_j;
}

void flatten_sums(vector<sum> sums, int* h_sum_starts_x, int* h_sum_starts_y, int* h_sum_ends_x, int* h_sum_ends_y, int* h_sum_hints, int* h_sum_lengths, int* h_sum_dirs, int* h_sum_pos_mins, int* h_sum_pos_maxs,  int no_sums) {

    for (int i = 0; i < no_sums; i++) {

        h_sum_starts_x[i] = sums[i].start.first;
        h_sum_starts_y[i] = sums[i].start.second;

        h_sum_ends_x[i] = sums[i].end.first;
        h_sum_ends_y[i] = sums[i].end.second;

        h_sum_hints[i] = sums[i].hint;
        h_sum_lengths[i] = sums[i].length;

        h_sum_dirs[i] = sums[i].dir;
        
        h_sum_pos_mins[i] = sums[i].posMin;
        h_sum_pos_maxs[i] = sums[i].posMax;
    }
}

__global__
void full_check_kernel(int* mat, int curr_i, int curr_j, int val, int m, int n, int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x, int* d_sum_ends_y, int* d_sum_hints,
    int* d_sum_lengths, int* d_sum_dirs, int * d_sum_mins, int * d_sum_maxs, int no_sums, volatile bool* partial_correct) {
    volatile __shared__ bool someoneFoundIt;

    int i = (blockDim.x * blockIdx.x) + threadIdx.x;
    if (threadIdx.x == 0) someoneFoundIt = *partial_correct;
    __syncthreads();
    if (someoneFoundIt && i < no_sums) {
        bool iFoundItFalse = !fullCheck(mat, curr_i, curr_j, val, m, n,
            d_sum_starts_x[i],
            d_sum_starts_y[i],
            d_sum_ends_x[i],
            d_sum_ends_y[i],
            d_sum_hints[i],
            d_sum_lengths[i],
            d_sum_dirs[i],
            d_sum_mins[i],
            d_sum_maxs[i]);
        if (iFoundItFalse) { someoneFoundIt = false; *partial_correct = false; }
        if (threadIdx.x == 0 && !(*partial_correct)) someoneFoundIt = false;
    }
}
///////////////////
// CUDA FUNCTIONS //
///////////////////


bool solution(int* h_sol_mat, int* d_sol_mat, int m, int n, stack<int> iter_i_stack, stack<int> iter_j_stack, stack<int> val_stack, int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x, int* d_sum_ends_y, int* d_sum_hints,
    int* d_sum_lengths, int* d_sum_dirs, int* d_sum_mins, int* d_sum_maxs, int no_sums) {

    const int GRIDSIZE = (no_sums + 1023) / 1024;
    const int THREADSIZE = min(1024, no_sums);
    bool* partial_correctness = (bool*)malloc(sizeof(bool));
    
    bool* d_partial_correctness;
    cudaMalloc(&d_partial_correctness, sizeof(bool));
    while (!iter_i_stack.empty()) {
        int iter_i = iter_i_stack.top();
        int iter_j = iter_j_stack.top();
        int curr_val = val_stack.top();

        iter_i_stack.pop();
        iter_j_stack.pop();
        val_stack.pop();

        if (curr_val < 10) {
            h_sol_mat[iter_i *n + iter_j] = curr_val;

            * partial_correctness = true;
            cudaMemcpy(d_partial_correctness, partial_correctness, sizeof(bool), cudaMemcpyHostToDevice);
            
            cudaMemcpy(d_sol_mat, h_sol_mat, (m * n) * sizeof(int), cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();
            
            full_check_kernel << < GRIDSIZE, THREADSIZE >> > (d_sol_mat, iter_i, iter_j, curr_val, m, n, d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints,
                d_sum_lengths, d_sum_dirs, d_sum_mins, d_sum_maxs, no_sums, d_partial_correctness);
            cudaDeviceSynchronize();
            cudaMemcpy(partial_correctness, d_partial_correctness, sizeof(bool), cudaMemcpyDeviceToHost);
            if (*partial_correctness) {
                iter_i_stack.push(iter_i);
                iter_j_stack.push(iter_j);
                val_stack.push(curr_val);

                int iter_i_next = mat_iter_get_next_i(h_sol_mat, m, n, iter_i_stack.top(), iter_j_stack.top());
                int iter_j_next = mat_iter_get_next_j(h_sol_mat, m, n, iter_i_stack.top(), iter_j_stack.top());

                if (iter_i_next == -999 || iter_j_next == -999) {
                    cout << "END INSIDE:" << endl;
                    cudaMemcpy(h_sol_mat, d_sol_mat, (m * n) * sizeof(int), cudaMemcpyDeviceToHost);
                    print_flattened_matrix(h_sol_mat, m, n);
                    return true;
                }
                //state_stack.push(state(iter_i_next, iter_j_next, 1));
                iter_i_stack.push(iter_i_next);
                iter_j_stack.push(iter_j_next);
                val_stack.push(1);
            }
            else {
                curr_val += 1;
                h_sol_mat[iter_i * n + iter_j] = -2;
                cudaMemcpy(d_sol_mat, h_sol_mat, (m * n) * sizeof(int), cudaMemcpyHostToDevice);
      
                iter_i_stack.push(iter_i);
                iter_j_stack.push(iter_j);
                val_stack.push(curr_val);
            }

        }
        else {
            while (!iter_i_stack.empty()) {
                int iter_i = iter_i_stack.top();
                int iter_j = iter_j_stack.top();
                int curr_val = val_stack.top();

                iter_i_stack.pop();
                iter_j_stack.pop();
                val_stack.pop();
                h_sol_mat[iter_i * n+ iter_j] = -2;
                cudaMemcpy(d_sol_mat, h_sol_mat, (m * n) * sizeof(int), cudaMemcpyHostToDevice);
                
                curr_val += 1;
                if (curr_val < 10) {
                    iter_i_stack.push(iter_i);
                    iter_j_stack.push(iter_j);
                    val_stack.push(curr_val);
                    break;
                }
            }
        }
    }
    return false;
}


int main(int argc, char** argv) {

    std::string filename("board5_1.kakuro");
    std::ifstream file;
    file.open(filename.c_str());

    int m, n;

    file >> m;
    file >> n;

    int** mat;
    read_matrix(mat, file, m, n);
    print_one_matrix(mat, m, n);

    int** sol_mat;
    convert_sol(mat, sol_mat, m, n);
    print_one_matrix(sol_mat, m, n);

    double start;
    double end;

    // CUDA
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("==prop== Running on device: %d -- %s \n", 0, prop.name);
    printf("==prop== #of SM -- %d \n", prop.multiProcessorCount);
    printf("==prop== Max Threads Per Block: -- %d \n", prop.maxThreadsPerBlock);

    vector<sum> sums = get_sums(mat, m, n);

    int no_sums = sums.size();

    // Flattening sums and matrix
    int* h_sum_starts_x = new int[no_sums];
    int* h_sum_starts_y = new int[no_sums];
    int* h_sum_ends_x = new int[no_sums];
    int* h_sum_ends_y = new int[no_sums];
    int* h_sum_hints = new int[no_sums];
    int* h_sum_lengths = new int[no_sums];
    int* h_sum_dirs = new int[no_sums];
    int* h_sum_mins = new int[no_sums];
    int* h_sum_maxs = new int[no_sums];

    flatten_sums(sums, h_sum_starts_x, h_sum_starts_y, h_sum_ends_x, h_sum_ends_y, h_sum_hints, h_sum_lengths, h_sum_dirs, h_sum_mins, h_sum_maxs, no_sums);

    print_flattened(h_sum_starts_x, h_sum_starts_y, h_sum_ends_x, h_sum_ends_y, h_sum_hints, h_sum_lengths, h_sum_dirs, no_sums);

    int* h_sol_mat;
    h_sol_mat = new int[m * n];
    flatten_sol_mat(sol_mat, h_sol_mat, m, n);

    print_flattened_matrix(h_sol_mat, m, n);

    // Declare device pointers and copy data into device
    int* d_sum_starts_x, * d_sum_starts_y, * d_sum_ends_x, * d_sum_ends_y, * d_sum_hints, * d_sum_lengths, * d_sum_dirs, * d_sum_mins, * d_sum_maxs, * d_sol_mat;

    cudaMalloc(&d_sum_starts_x, no_sums * sizeof(int));
    cudaMalloc(&d_sum_starts_y, no_sums * sizeof(int));
    cudaMalloc(&d_sum_ends_x, no_sums * sizeof(int));
    cudaMalloc(&d_sum_ends_y, no_sums * sizeof(int));
    cudaMalloc(&d_sum_hints, no_sums * sizeof(int));
    cudaMalloc(&d_sum_lengths, no_sums * sizeof(int));
    cudaMalloc(&d_sum_dirs, no_sums * sizeof(int));
    cudaMalloc(&d_sum_mins, no_sums * sizeof(int));
    cudaMalloc(&d_sum_maxs, no_sums * sizeof(int));

    cudaMalloc(&d_sol_mat, (m * n) * sizeof(int));
    
    cudaMemcpy(d_sum_starts_x, h_sum_starts_x, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_starts_y, h_sum_starts_y, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_ends_x, h_sum_ends_x, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_ends_y, h_sum_ends_y, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_hints, h_sum_hints, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_lengths, h_sum_lengths, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_dirs, h_sum_dirs, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_mins, h_sum_mins, no_sums * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sum_maxs, h_sum_maxs, no_sums * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_sol_mat, h_sol_mat, (m * n) * sizeof(int), cudaMemcpyHostToDevice);

    int iter_i = mat_iter_init_i(h_sol_mat, m, n);
    int iter_j = mat_iter_init_j(h_sol_mat, m, n);
    
    stack<int> iter_i_stack, iter_j_stack, val_stack;
    iter_i_stack.push(iter_i);
    iter_j_stack.push(iter_j);
    val_stack.push(1);

    start = omp_get_wtime();
    bool result = solution(h_sol_mat, d_sol_mat, m, n, iter_i_stack, iter_j_stack, val_stack, d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints,
        d_sum_lengths, d_sum_dirs, d_sum_mins, d_sum_maxs, no_sums);
    end = omp_get_wtime();
    printf("Work took %f seconds\n", end - start);

    //sol_to_file(mat, sol_mat, m, n, "solution.kakuro");
    if (result) {
        cout << "SUCCESS" << endl;
    }
    else {
        cout << "COULD NOT SOLVE" << endl;
    }
    cudaDeviceSynchronize();
    // DELETE PART 
    for (int i = 0; i < m; i++) {
        delete mat[i];
        delete sol_mat[i];
    }

    delete mat;
    delete sol_mat;

    delete h_sum_starts_x;
    delete h_sum_starts_y;
    delete h_sum_ends_x;
    delete h_sum_ends_y;
    delete h_sum_hints;
    delete h_sum_lengths;
    delete h_sum_dirs;
    delete h_sol_mat;

    cudaFree(d_sum_starts_x);
    cudaFree(d_sum_starts_y);
    cudaFree(d_sum_ends_x);
    cudaFree(d_sum_ends_y);
    cudaFree(d_sum_hints);
    cudaFree(d_sum_lengths);
    cudaFree(d_sum_dirs);
    cudaFree(d_sol_mat);
    cudaFree(d_sum_mins);
    cudaFree(d_sum_maxs);

    return 0;
}
