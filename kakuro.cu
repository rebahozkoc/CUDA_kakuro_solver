#include <iostream>
#include <string>

#include <fstream>
#include <sstream>
#include <vector>

#include <bits/stdc++.h>
#include <array>
#include <omp.h>
#include <stack>

using namespace std;

int MAX_THREADS;

enum direction{
    d_down,
    d_right,
    none
};

#define COORD std::pair<int, int>

//#define DEBUG

int iter = 0;

/// Auxiliary functions

void display_arr(int *arr, int n){
    cout << "arr: ";
    for (int i = 0; i < n; i++){
        cout << arr[i] << " ";
    }
    cout << endl;
}

void print_coords(COORD start, COORD end){

    cout << "Start:" << start.first << "," << start.second << endl;
    cout << "End:" << end.first << "," << end.second << endl;
}

int find_length(COORD start, COORD end, direction dir){
    if (dir == d_down)
        return end.first - start.first;
    if (dir == d_right)
        return end.second - start.second;
    return -1;
}

void convert_sol(int **mat, int **&sol_mat, int m, int n){

    sol_mat = new int *[m]; // Rows
    for (int i = 0; i < m; i++){
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

void print_one_matrix(int **matrix, int m, int n){
    std::cout << "Matrix: " << std::endl;
    for (int i = 0; i < m; i++) { // rows
        for (int j = 0; j < n; j++) { // cols
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << "\n";
    }
}

void sol_to_file(int **mat, int **sol_mat, int m, int n, string fname){
    ofstream to_write(fname);

    to_write << m << " " << n << "\n";

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            if (mat[i][j] != -2)
                to_write << mat[i][j] << " ";
            else
                to_write << sol_mat[i][j] << " ";
        }
        to_write << "\n";
    }

    to_write.close();
}

void read_matrix(int **&matrix, std::ifstream &afile, int m, int n){

    matrix = new int *[m]; // rows

    for (int i = 0; i < m; i++){
        matrix[i] = new int[n]; // cols
    }

    int val;
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            afile >> val;
            matrix[i][j] = val;
        }
    }
}

/// Auxiliary functions

struct sum{
    COORD start;
    COORD end;

    int hint;
    int dir;
    int length;
    int posMin;
    int posMax;

    void print_sum(){
        cout << "############################" << endl;
        cout << "Creating sum with: " << endl;
        print_coords(start, end);
        cout << "Hint: " << hint << endl;
        cout << "Direction: " << dir << endl;
        cout << "Length: " << length << endl;
        cout << endl;
        cout << "############################" << endl;
    }

    sum(COORD _start, COORD _end, int _hint, direction _dir) : start(_start), end(_end), hint(_hint), dir(_dir){
        length = find_length(_start, _end, _dir);
        // This is equal to hint - sum of numbers 9 + 8 + 7
        posMin = hint - 45 + ((8-length) * (9-length)) / 2; 

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
    bool checkSumPartial(int** mat) const{
        int sum = 0;
        if (dir == 1){ // right
            for (int j = start.second; j < end.second; j++){
                if (mat[start.first][j] != -2){
                        sum += mat[start.first][j];
                }
                
            }
        }else{ // down
            for (int i = start.first; i < end.first; i++){
                if (mat[i][start.second] != -2){
                    sum += mat[i][start.second];
                }
            }
        }       
        
        if (sum > hint){
            return false;
        }        
        return true;
    }

    bool checkSumComplete(int ** mat) const {
        int sum = 0;
        if (dir == 1){ // right
            for (int j = start.second; j < end.second; j++){
                if (mat[start.first][j] != -2){
                        sum += mat[start.first][j];
                }
                
            }
        }else{ // down
            for (int i = start.first; i < end.first; i++){
                if (mat[i][start.second] != -2){
                    sum += mat[i][start.second];
                }
            }
        }       
        if (sum != hint){
            return false;
        }
        return true;
    }

    bool areElementsUniqueExceptEmpties(int** mat) const {
        std::unordered_set<int> unique_elements;
        if (dir == 1){ // right
            for (int j = start.second; j < end.second; j++){
                if (mat[start.first][j] != -2){
                    if (unique_elements.find(mat[start.first][j]) == unique_elements.end()){
                        unique_elements.insert(mat[start.first][j]);
                    }else{
                        // If the element is found, it means it's a duplicate, and we can return false
                        return false;
                    }  
                }
                
            }
        }
        else{ // down
            for (int i = start.first; i < end.first; i++){
                if (mat[i][start.second] != -2){
                    if (unique_elements.find(mat[i][start.second]) == unique_elements.end()){
                        unique_elements.insert(mat[i][start.second]);
                    }else{
                        return false;
                    }
                }
            }
        }       
        return true;
    }

    bool areElementsUnique(int** mat) const {
        std::unordered_set<int> unique_elements;
        if (dir == 1){ // right
            for (int j = start.second; j < end.second; j++){
                if (unique_elements.find(mat[start.first][j]) == unique_elements.end()){
                    unique_elements.insert(mat[start.first][j]);
                }else{
                    // If the element is found, it means it's a duplicate, and we can return false
                    return false;
                }  
            }
        }else{ // down
            for (int i = start.first; i < end.first; i++){
                if (unique_elements.find(mat[i][start.second]) == unique_elements.end()){
                    unique_elements.insert(mat[i][start.second]);
                }else{
                    return false;
                }
            }
        }       
        return true;
    }

    bool isArrFull(int** mat) const{
        if (dir == 1){ // right
            for (int j = start.second; j < end.second; j++){
                if (mat[start.first][j] == -2)
                    return false;
            }
        }else{ // down
            for (int i = start.first; i < end.first; i++){
                if (mat[i][start.second] == -2)
                    return false;
            }
        }       
        return true;
    }

    bool isACandidate(int val) const{
        if (val < posMin){
            return false;
        }
        if (val > posMax){
            return false;
        }
        return true;
    }

    bool fullCheck(int** mat, COORD curr, int val) const{
        int i = curr.first;
        int j = curr.second;
        if (dir == 1){ // right
            // current coordinate is not in the sum
            if (!((i == start.first) && (j >= start.second) && (j < end.second))){
                return true;
            }
        }else{ // down
            if(!(j == start.second && i >= start.first && i < end.first)){
                return true;
            }
        }
        if (!isACandidate(val)){
            return false;
        }
        if (!areElementsUniqueExceptEmpties(mat)){
            return false;
        }
        if (isArrFull(mat) || isLastCell(curr)){
            return checkSumComplete(mat) && areElementsUnique(mat);
        }else{
            return checkSumPartial(mat);
        }
    }

    bool isLastCell(COORD curr) const{
        if (dir == 1){ // to right
            if (start.first == curr.first && end.second - 1 == curr.second){
                return true;
            }
        }else{
            if (start.second == curr.second && end.first - 1 == curr.first){{
                return true;
            }}
        }
        return false;
    }
};

COORD find_end(int **matrix, int m, int n, int i, int j, direction dir){ // 0 down 1 right

    if (dir == d_right){
        for (int jj = j + 1; jj < n; jj++){
            if (matrix[i][jj] != -2 || jj == n - 1){
                if (matrix[i][jj] == -2 && jj == n - 1)
                    jj++;
                COORD END = COORD(i, jj);
                return END;
            }
        }
    }

    if (dir == d_down){
        for (int ii = i + 1; ii < m; ii++){
            if (matrix[ii][j] != -2 || ii == m - 1){
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

vector<sum> get_sums(int **matrix, int m, int n){

    vector<sum> sums;

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            int val = matrix[i][j];
            if (val != -1 && val != -2){
                int hint = val;
                hint = hint / 10;
                // right sum
                if ((hint % 100) == 0){
                    hint = (int)(hint / 100);
                    COORD START = COORD(i, j + 1);
                    COORD END = find_end(matrix, m, n, i, j, d_right);
                    sum _sum = sum(START, END, hint, d_right);
                    sums.push_back(_sum);
                }

                else{
                    int div = (int)(hint / 100);
                    int rem = (int)(hint % 100);
                    // down sum
                    if (div == 0 && rem != 0){
                        COORD START = COORD(i + 1, j);
                        COORD END = find_end(matrix, m, n, i, j, d_down);
                        sum _sum = sum(START, END, rem, d_down);
                        sums.push_back(_sum);
                    }
                    // combined sum
                    if (div != 0 && rem != 0){
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

struct mat_iter{
    COORD curr;
    int **mat;
    int m, n;

    mat_iter(int **_mat, int _m, int _n) : mat(_mat), m(_m), n(_n){
        curr = COORD(0, 0);
        if (mat[curr.first][curr.second] != -2){
            set_next();
        }
    }

    // Copy constructor
    mat_iter(const mat_iter &other) {
        curr = COORD(other.curr.first, other.curr.second);
        mat = other.mat;
        m = other.m;
        n = other.n;
    }

    void set_next(){
        do {
            if (curr.second + 1 < n){
                curr = COORD(curr.first, curr.second +1);
            }else{
                curr = COORD(curr.first + 1, 0);
            }
        }while (curr.first != m && (mat[curr.first][curr.second] != -2 ));
        if (curr.first  == m ){
                curr.first = -999; // end of the iteration
                curr.second = -999;

        }
    }
};

//////////////////////////////////////////////
// Auxiliary functions for preparing problem //
//////////////////////////////////////////////

///////////////////////////////////////////////////
// Auxiliary functions for preparing CUDA setting //
///////////////////////////////////////////////////

void flatten_sums(vector<sum> sums, int *h_sum_starts_x, int *h_sum_starts_y, int *h_sum_ends_x, int *h_sum_ends_y, int *h_sum_hints, int *h_sum_lengths, int *h_sum_dirs, int no_sums){

	for (int i = 0; i < no_sums; i++){

		h_sum_starts_x[i] = sums[i].start.first;
		h_sum_starts_y[i] = sums[i].start.second;

		h_sum_ends_x[i] = sums[i].end.first;
		h_sum_ends_y[i] = sums[i].end.second;

		h_sum_hints[i] = sums[i].hint;
		h_sum_lengths[i] = sums[i].length;

		h_sum_dirs[i] = sums[i].dir;
	}
}

void print_flattened(int *h_sum_starts_x, int *h_sum_starts_y, int *h_sum_ends_x, int *h_sum_ends_y, int *h_sum_hints, int *h_sum_lengths, int *h_sum_dirs, int no_sums){

	cout << "###h_sum_starts_x: " << endl;
	for (int i = 0; i < no_sums; i++){
		cout << h_sum_starts_x[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_starts_y: " << endl;
	for (int i = 0; i < no_sums; i++){
		cout << h_sum_starts_y[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_ends_x: " << endl;
	for (int i = 0; i < no_sums; i++){
		cout << h_sum_ends_x[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_ends_y: " << endl;
	for (int i = 0; i < no_sums; i++){
		cout << h_sum_ends_y[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_hints: " << endl;
	for (int i = 0; i < no_sums; i++){
		cout << h_sum_hints[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_lengths: " << endl;
	for (int i = 0; i < no_sums; i++){
		cout << h_sum_lengths[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_dirs: " << endl;
	for (int i = 0; i < no_sums; i++){
		cout << h_sum_dirs[i] << " ";
	}
	cout << endl;
}

void flatten_sol_mat(int **sol_mat, int *h_sol_mat, int m, int n){
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			h_sol_mat[i * n + j] = sol_mat[i][j];
		}
	}
}

void print_flattened_matrix(int *h_sol_mat, int m, int n){

	cout << "###Flattened matrix: " << endl;
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			cout << h_sol_mat[i * n + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

///////////////////////////////////////////////////
// Auxiliary functions for preparing CUDA setting //
///////////////////////////////////////////////////
struct state{
    mat_iter iter;
    int val;
    state(mat_iter _iter, int _val): iter(_iter), val(_val){}
};
///////////////////
// CUDA FUNCTIONS //
///////////////////

__global__ 
    void kakuro_kernel(int *d_sum_starts_x, int *d_sum_starts_y, int *d_sum_ends_x, int *d_sum_ends_y, int *d_sum_hints, 
    int *d_sum_lengths, int *d_sum_dirs, int *d_sol_mat, int *d_t_mats, int m, int n, int no_sums,  volatile bool *solved){
        *solved = true;
            

	// TO DO

	// About volatile bool* solved:
	// You can get idea from https://stackoverflow.com/questions/12505750/how-can-a-global-function-return-a-value-or-break-out-like-c-c-does%5B/url%5D for how to break out of a CUDA kernel
	// You may or may not use it
}

///////////////////
// CUDA FUNCTIONS //
///////////////////



bool solution(int **sol_mat, vector<sum> &sums, int m, int n, stack<state> state_stack) {    
    while(!state_stack.empty()){
        mat_iter iter = state_stack.top().iter;
        int curr_val = state_stack.top().val;
        state_stack.pop();
        if (curr_val < 10){
        
            // Update the matrix with the new val for iter.curr
            sol_mat[iter.curr.first][iter.curr.second] = curr_val;
            bool partial_correctness = true;
            for (int i = 0; i< sums.size(); i++){
                if (!sums[i].fullCheck(sol_mat, iter.curr, curr_val)){
                    partial_correctness = false;
                    break;
                }
            }
            if (partial_correctness){
                state_stack.push(state(iter, curr_val));
                mat_iter temp(iter);
                temp.set_next();
                if (temp.curr.first == -999 || temp.curr.second == -999){
                    cout << "END INSIDE:" << endl;
                    print_one_matrix(sol_mat, m, n);
                    return true;
                }
                state_stack.push(state(temp, 1));
            }else{
                curr_val += 1;
                sol_mat[iter.curr.first][iter.curr.second] = -2;
                state_stack.push(state(iter, curr_val));
            }
                
        }else{
            while(!state_stack.empty()){
                mat_iter iter = state_stack.top().iter;
                int curr_val = state_stack.top().val;
                state_stack.pop();
                sol_mat[iter.curr.first][iter.curr.second] = -2;
                curr_val += 1;
                if (curr_val < 10){
                    state_stack.push(state(iter, curr_val));
                    break;
                }
            }
        }
    }
    return false;
}


int main(int argc, char **argv){

    std::string filename(argv[1]);
    std::ifstream file;
    file.open(filename.c_str());

    int m, n;
    
    file >> m;
    file >> n;

    int **mat;
    read_matrix(mat, file, m, n);
    print_one_matrix(mat, m, n);

    int **sol_mat;
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
    
    int grid_dim =  m;		// TO DO
	int block_dim = n;// To DO

    int no_sums = sums.size();

	// Flattening sums and matrix
	int *h_sum_starts_x = new int[no_sums];
	int *h_sum_starts_y = new int[no_sums];
	int *h_sum_ends_x = new int[no_sums];
	int *h_sum_ends_y = new int[no_sums];
	int *h_sum_hints = new int[no_sums];
	int *h_sum_lengths = new int[no_sums];
	int *h_sum_dirs = new int[no_sums];


    flatten_sums(sums, h_sum_starts_x, h_sum_starts_y, h_sum_ends_x, h_sum_ends_y, h_sum_hints, h_sum_lengths, h_sum_dirs, no_sums);

	print_flattened(h_sum_starts_x, h_sum_starts_y, h_sum_ends_x, h_sum_ends_y, h_sum_hints, h_sum_lengths, h_sum_dirs, no_sums);

    int *h_sol_mat;
	h_sol_mat = new int[m * n];
	flatten_sol_mat(sol_mat, h_sol_mat, m, n);

	print_flattened_matrix(h_sol_mat, m, n);

    // Declare device pointers and copy data into device
	int *d_sum_starts_x, *d_sum_starts_y, *d_sum_ends_x, *d_sum_ends_y, *d_sum_hints, *d_sum_lengths, *d_sum_dirs, *d_sol_mat, *d_t_mats;

	cudaMalloc(&d_sum_starts_x, no_sums * sizeof(int));
	cudaMalloc(&d_sum_starts_y, no_sums * sizeof(int));
	cudaMalloc(&d_sum_ends_x, no_sums * sizeof(int));
	cudaMalloc(&d_sum_ends_y, no_sums * sizeof(int));
	cudaMalloc(&d_sum_hints, no_sums * sizeof(int));
	cudaMalloc(&d_sum_lengths, no_sums * sizeof(int));
	cudaMalloc(&d_sum_dirs, no_sums * sizeof(int));
	cudaMalloc(&d_sol_mat, (m * n) * sizeof(int));
	cudaMalloc(&d_t_mats, (m * n * grid_dim * block_dim) * sizeof(int)); // Allocating invidual matrix for each GPU thread
	// You may use this array if you will implement a thread-wise solution

    cudaMemcpy(d_sum_starts_x, h_sum_starts_x, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_starts_y, h_sum_starts_y, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_ends_x, h_sum_ends_x, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_ends_y, h_sum_ends_y, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_hints, h_sum_hints, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_lengths, h_sum_lengths, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_dirs, h_sum_dirs, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sol_mat, h_sol_mat, (m * n) * sizeof(int), cudaMemcpyHostToDevice);

    bool *solved = (bool*)malloc(sizeof(bool));
	*solved = false;
	bool *d_solved;

	cudaMalloc(&d_solved, sizeof(bool));
	cudaMemcpy(d_solved, solved, sizeof(bool), cudaMemcpyHostToDevice);
    
    mat_iter iter = mat_iter(sol_mat,  m, n);
    stack<state> state_stack;
    state_stack.push(state(iter, 1));
    cout << "Number of threads: " << omp_get_num_threads() << endl;
    
    start = omp_get_wtime();
    kakuro_kernel<<<grid_dim, block_dim>>>(d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints,
										   d_sum_lengths, d_sum_dirs, d_sol_mat, d_t_mats, m, n,
										   no_sums, d_solved);
	cudaDeviceSynchronize();
    bool result = solution(sol_mat, sums, m, n, state_stack);
    end = omp_get_wtime(); 
    printf("Work took %f seconds\n", end - start);

    sol_to_file(mat, sol_mat, m, n, "solution.kakuro");
    if (result){
        cout << "SUCCESS" << endl;
    }else{
        cout << "COULD NOT SOLVE" << endl;
    }
    
    cout << "SOLVED BY CUDA before copy: " << *solved << endl; 

    cudaMemcpy(solved, d_solved, sizeof(bool), cudaMemcpyDeviceToHost);
    cout << "SOLVED BY CUDA?: " << *solved << endl; 
    


    // DELETE PART 
    for (int i = 0; i < m; i++){
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

	cudaFree(d_t_mats);
	cudaFree(d_sum_starts_x);
	cudaFree(d_sum_starts_y);
	cudaFree(d_sum_ends_x);
	cudaFree(d_sum_ends_y);
	cudaFree(d_sum_hints);
	cudaFree(d_sum_lengths);
	cudaFree(d_sum_dirs);
	cudaFree(d_sol_mat);

    return 0;
}
