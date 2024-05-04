#include <iostream>
#include <mpi.h>
#include <vector>
#include <string>
#include <chrono>
using namespace std;
void print_matrix(vector<vector<int>> matrix) {
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[0].size(); ++j) {
            cout << to_string(matrix[i][j]) + " ";
        }
        cout << '\n';
    }
    cout << '\n';
}
vector<vector<int>> initialize_matrix(int sizei, int sizej) {
    vector<vector<int>> matrix(sizei, vector<int>(sizej));
    for (int i = 0; i < sizei; ++i) {
        for (int j = 0; j < sizej; ++j) {
            matrix[i][j] = rand() % 10;
        }
    }
    return matrix;
}
int* flattenVector(vector<vector<int>>& input) {
    int* flattened = new int[input.size()*input[0].size()];
    int index=0;
    for (int i = 0; i < input.size(); ++i) {
        for (int j = 0; j <input[0].size(); ++j) {
            flattened[index] = input[i][j];
            index++;
        }
    }
    return flattened;
}vector<vector<int>> unflattenVector(int* flattened, int rows, int cols) {
    vector<vector<int>> unflattened(rows, vector<int>(cols));
    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            unflattened[i][j] = flattened[index];
            index++;
        }
    }
    return unflattened;
}
vector<vector<int>> split(vector<vector<int>> matrix, int quarterRow, int quarterCol) {
    if (matrix.size() == 1) {
        return matrix;
    }
    int newRows = matrix.size() / 2;
    int newCols = matrix[0].size() / 2;

    int startRow = (quarterRow - 1) * newRows;
    int startCol = (quarterCol - 1) * newCols;

    vector<vector<int>> result(newRows, vector<int>(newCols));
    for (int i = 0; i < newRows; i++) {
        for (int j = 0; j < newCols; j++) {
            result[i][j] = matrix[startRow + i][startCol + j];
        }
    }
    return result;
}
void matrix_multiply(vector<vector<int>>& C, vector<vector<int>>& A, vector<vector<int>>& B) {

    int n = A[0].size();
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
    }
    else {
        vector<vector<int>> T(n, vector<int>(n));

        vector<vector<int>> A11 = split(A, 1, 1);
        vector<vector<int>> A12 = split(A, 1, 2);
        vector<vector<int>> A21 = split(A, 2, 1);
        vector<vector<int>> A22 = split(A, 2, 2);

        vector<vector<int>> B11 = split(B, 1, 1);
        vector<vector<int>> B12 = split(B, 1, 2);
        vector<vector<int>> B21 = split(B, 2, 1);
        vector<vector<int>> B22 = split(B, 2, 2);

        vector<vector<int>> C11 = split(C, 1, 1);
        vector<vector<int>> C12 = split(C, 1, 2);
        vector<vector<int>> C21 = split(C, 2, 1);
        vector<vector<int>> C22 = split(C, 2, 2);

        vector<vector<int>> T11 = split(T, 1, 1);
        vector<vector<int>> T12 = split(T, 1, 2);
        vector<vector<int>> T21 = split(T, 2, 1);
        vector<vector<int>> T22 = split(T, 2, 2);

         matrix_multiply(C11, A11, B11);
         matrix_multiply(C12, A11, B12);
         matrix_multiply(C21, A21, B11);
         matrix_multiply(C22, A21, B12);
         matrix_multiply(C11, A11, B11);
         matrix_multiply(T11, A12, B21);
         matrix_multiply(T12, A12, B22);
         matrix_multiply(T21, A22, B21);
         matrix_multiply(T22, A22, B22);

        for (int i = 0; i < n / 2; i++) {
            for (int j = 0; j < n / 2; j++) {
                C[i][j] = C11[i][j] + T11[i][j];
            }
        }

        for (int i = n / 2; i < n; i++) {
            for (int j = 0; j < n / 2; j++) {
                C[i][j] = C21[i - n / 2][j] + T21[i - n / 2][j];
            }
        }
        for (int i = n / 2; i < n; i++) {
            for (int j = n / 2; j < n; j++) {
                C[i][j] = C22[i - n / 2][j - n / 2] + T22[i - n / 2][j - n / 2];
            }
        }

        for (int i = 0; i < n / 2; i++) {
            for (int j = n / 2; j < n; j++) {
                C[i][j] = C12[i][j - n / 2] + T12[i][j - n / 2];
            }
        }
    }
}
void p_matrix_multiply(vector<vector<int>>& C, vector<vector<int>>& A, vector<vector<int>>& B) {

    int n = A[0].size();
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
    }
    else {
        vector<vector<int>> T(n, vector<int>(n));

        vector<vector<int>> A11;
        vector<vector<int>> A12;
        vector<vector<int>> A21;
        vector<vector<int>> A22;

        vector<vector<int>> B11;
        vector<vector<int>> B12;
        vector<vector<int>> B21;
        vector<vector<int>> B22;

        vector<vector<int>> C11;
        vector<vector<int>> C12;
        vector<vector<int>> C21;
        vector<vector<int>> C22;
        vector<vector<int>> T11;
        vector<vector<int>> T12;
        vector<vector<int>> T21;
        vector<vector<int>> T22;

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            C11 = split(C, 1, 1);
            B11 = split(B, 1, 1);
            A11 = split(A, 1, 1);
            matrix_multiply(C11, A11, B11);
        }
        if (rank == 1) {
            C12 = split(C, 1, 2);
            A11 = split(A, 1, 1);
            B12 = split(B, 1, 2);
            matrix_multiply(C12, A11, B12);
            MPI_Rsend(flattenVector(C12), C12.size() * C12[0].size(), MPI_INT, 0,0, MPI_COMM_WORLD);

        }
        if (rank == 2) {
            C21 = split(C, 2, 1);
            A21 = split(A, 2, 1);
            B11 = split(B, 1, 1);
            matrix_multiply(C21, A21, B11);
            MPI_Rsend(flattenVector(C21), C21.size() * C21[0].size(), MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 3) {
            C22 = split(C, 2, 2);
            A21 = split(A, 2, 1);
            B12 = split(B, 1, 2);
            matrix_multiply(C22, A21, B12);
            MPI_Rsend(flattenVector(C22), C22.size() * C22[0].size(), MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 4) {
            T11 = split(T, 1, 1);
            A12 = split(A, 1, 2);
            B21 = split(B, 2, 1);
            matrix_multiply(T11, A12, B21);
            MPI_Rsend(flattenVector(T11), T11.size() * T11[0].size(), MPI_INT, 0,0, MPI_COMM_WORLD);
        }
        if (rank == 5) {
            T12 = split(T, 1, 2);
            A12 = split(A, 1, 2);
            B22 = split(B, 2, 2);
            matrix_multiply(T12, A12, B22);
            MPI_Rsend(flattenVector(T12), T12.size() * T12[0].size(), MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 6) {
            T21 = split(T, 2, 1);
            A22 = split(A, 2, 2);
            B21 = split(B, 2, 1);
            matrix_multiply(T21, A22, B21);
            MPI_Rsend(flattenVector(T21), T21.size() * T21[0].size(), MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 7) {
            T22 = split(T, 2, 2);
            A22 = split(A, 2, 2);
            B22 = split(B, 2, 2);
            matrix_multiply(T22, A22, B22);
            MPI_Rsend(flattenVector(T22), T22.size() * T22[0].size(), MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 0) {
            MPI_Status status;
            int size = (C.size()/2) * (C[0].size()/2);
            int rc = C.size() / 2;
            int* data = new int[size];
            for (int i = 1; i <= 7; i++) {
                MPI_Recv(data, size, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                if (status.MPI_SOURCE == 1) C12 = unflattenVector(data, rc, rc);
                if (status.MPI_SOURCE == 2) C21 = unflattenVector(data, rc, rc);
                if (status.MPI_SOURCE == 3) C22 = unflattenVector(data, rc, rc);
                if (status.MPI_SOURCE == 4) T11 = unflattenVector(data, rc, rc);
                if (status.MPI_SOURCE == 5) T12 = unflattenVector(data, rc, rc);
                if (status.MPI_SOURCE == 6) T21 = unflattenVector(data, rc, rc);
                if (status.MPI_SOURCE == 7) T22 = unflattenVector(data, rc, rc);
            }

            for (int i = 0; i < n / 2; i++) {
                for (int j = 0; j < n / 2; j++) {
                    C[i][j] = C11[i][j] + T11[i][j];
                }
            }
            for (int i = n / 2; i < n; i++) {
                for (int j = 0; j < n / 2; j++) {
                    C[i][j] = C21[i - n / 2][j] + T21[i - n / 2][j];
                }
            }
            for (int i = n / 2; i < n; i++) {
                for (int j = n / 2; j < n; j++) {
                    C[i][j] = C22[i - n / 2][j - n / 2] + T22[i - n / 2][j - n / 2];
                }
            }
            for (int i = 0; i < n / 2; i++) {
                for (int j = n / 2; j < n; j++) {
                    C[i][j] = C12[i][j - n / 2] + T12[i][j - n / 2];
                }
            }
        }
    }
}
int main() {
    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto A = initialize_matrix(128, 128);
    auto B = initialize_matrix(128, 128);
    /*if (rank == 0) {
        print_matrix(A);
        print_matrix(B);
    }*/
    vector<vector<int>> C(A.size(), vector<int>(A.size()));
    /*p_matrix_multiply(C, A, B);
    if(rank==0) print_matrix(C);*/
    chrono::steady_clock::time_point begin;
    if(rank==0) begin= chrono::steady_clock::now();
    p_matrix_multiply(C, A, B);
    chrono::steady_clock::time_point end;
    if (rank == 0) {
        end = std::chrono::steady_clock::now();
        cout << "parallelism: " + to_string(std::chrono::duration_cast<std::chrono::seconds>(end - begin).count()) + " seconds" << endl;
    }
    MPI_Finalize();
    return 0;
}