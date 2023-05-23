#include <iostream>
#include <string>

#include <fstream>
#include <sstream>
#include <vector>

#include <bits/stdc++.h>
#include <array>
#include <omp.h>

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


bool solution(int **sol_mat, vector<sum> &sums, int m, int n, mat_iter iter) {    
    if (iter.curr.first == -999 || iter.curr.second == -999){
        cout << "END:" << endl;
        print_one_matrix(sol_mat, m, n);
        return true;
    }else{
        for (int val = 1; val < 10; val++){
            // Update the sums vector elements with the new val for iter.curr
            sol_mat[iter.curr.first][iter.curr.second] = val;
            bool partial_correctness = true;
            
            for (int i = 0; i< sums.size(); i++){
                if (!sums[i].fullCheck(sol_mat, iter.curr, val)){
                    partial_correctness = false;
                    break;
                }
            }
            
            if (partial_correctness){
                bool sol;
                mat_iter temp(iter);
                temp.set_next();
                sol = solution(sol_mat, sums, m, n, temp);                
                if (sol){
                    return true;
                }
            }
            sol_mat[iter.curr.first][iter.curr.second] = -2;
        }
        return false;
    }
    return true;
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


    vector<sum> sums = get_sums(mat, m, n);
    mat_iter iter = mat_iter(sol_mat,  m, n);

    cout << "Number of threads: " << omp_get_num_threads() << endl;
    
    double start;
    double end;
    start = omp_get_wtime();
    bool result = solution(sol_mat, sums, m, n, iter);
    end = omp_get_wtime(); 
    printf("Work took %f seconds\n", end - start);

    sol_to_file(mat, sol_mat, m, n, "solution.kakuro");
    if (result){
        cout << "SUCCESS" << endl;
    }else{
        cout << "COULD NOT SOLVE" << endl;
    }
    
    for (int i = 0; i < m; i++){
        delete mat[i];
        delete sol_mat[i];
    }

    delete mat;
    delete sol_mat;
    return 0;
}
