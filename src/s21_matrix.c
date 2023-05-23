#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error_code = OK;
  if ((rows < 1) || (columns < 1) || result == NULL) {
    error_code = INCORRECT_MATRIX;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double *));
    if (result->matrix == NULL) {
      error_code = INCORRECT_MATRIX;
    } else {
      for (int i = 0; i < rows; i++)
        result->matrix[i] = calloc(columns, sizeof(double));
    }
  }
  return error_code;
}

void s21_remove_matrix(matrix_t *A) {
  for (int k = 0; k < A->rows; k++) {
    free(A->matrix[k]);
  }
  free(A->matrix);
  A->columns = 0;
  A->rows = 0;
  A->matrix = NULL;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;
  if (s21_it_is_normal_matrix(A), s21_it_is_normal_matrix(B)) {
    if (A->rows == B->rows && A->columns == B->columns) {
      for (int i = 0; i < A->rows && res == SUCCESS; i++) {
        for (int j = 0; j < A->columns && res == SUCCESS; j++) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
            res = FAILURE;
          }
        }
      }
    } else {
      res = FAILURE;
    }
  } else {
    res = FAILURE;
  }
  return res;
}

int s21_sum_or_sub(matrix_t *A, matrix_t *B, matrix_t *result, bool sign) {
  int error_code = CALC_ERROR;
  if (!s21_it_is_normal_matrix(A) || !s21_it_is_normal_matrix(B)) {
    error_code = INCORRECT_MATRIX;
  } else {
    if (s21_it_is_normal_matrix(A) && s21_it_is_normal_matrix(B) &&
        (A->rows == B->rows) && (A->columns == B->columns)) {
      if (s21_create_matrix(A->rows, A->columns, result) == 0) {
        error_code = OK;
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            if (sign == true) {
              result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
            } else {
              result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
            }
          }
        }
      }
    } else {
      error_code = CALC_ERROR;
    }
  }
  return error_code;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  return s21_sum_or_sub(A, B, result, true);
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  return s21_sum_or_sub(A, B, result, false);
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error_code = INCORRECT_MATRIX;
  if (s21_it_is_normal_matrix(A)) {
    if (s21_create_matrix(A->rows, A->columns, result) == 0) {
      error_code = OK;
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    } else
      error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error_code = INCORRECT_MATRIX;
  if (s21_it_is_normal_matrix(A) && s21_it_is_normal_matrix(B)) {
    if (A->columns != B->rows) {
      error_code = CALC_ERROR;
    } else {
      if (s21_create_matrix(A->rows, B->columns, result) == 0) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < B->columns; j++) {
            result->matrix[i][j] = 0;
            for (int k = 0; k < A->columns; k++) {
              result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
            }
          }
        }
        error_code = OK;
      } else {
        error_code = INCORRECT_MATRIX;
      }
    }
  }
  return error_code;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error_code = INCORRECT_MATRIX;
  if (s21_it_is_normal_matrix(A)) {
    if (s21_create_matrix(A->columns, A->rows, result) == 0) {
      error_code = OK;
      for (int i = 0; i < A->columns; i++) {
        for (int j = 0; j < A->rows; j++) {
          result->matrix[i][j] = A->matrix[j][i];
        }
      }
    } else
      error_code = INCORRECT_MATRIX;
  }
  return error_code;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error_code = INCORRECT_MATRIX;
  if (s21_it_is_normal_matrix(A)) {
    error_code = OK;
    if (A->rows != A->columns) {
      error_code = CALC_ERROR;
    } else {
      s21_create_matrix(A->rows, A->columns, result);
      if (A->rows == 1) {
        result->matrix[0][0] = A->matrix[0][0];
      } else {
        int sign = 0;
        matrix_t minor;
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            s21_minor(i, j, A, &minor);
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            double det = 0;
            s21_determinant(&minor, &det);
            result->matrix[i][j] = det * sign;
            s21_remove_matrix(&minor);
          }
        }
      }
    }
  }
  return error_code;
}

int s21_determinant(matrix_t *A, double *result) {
  int error_code = INCORRECT_MATRIX;
  if (s21_it_is_normal_matrix(A)) {
    error_code = OK;
    if (A->rows != A->columns) {
      error_code = CALC_ERROR;
    } else {
      *result = det(A);
    }
  }
  return error_code;
}

double det(matrix_t *A) {
  double result = 0;
  if (A->rows == 1) {
    result = A->matrix[0][0];
  } else {
    if (A->rows == 2) {
      result = (A->matrix[0][0] * A->matrix[1][1]) -
               (A->matrix[0][1] * A->matrix[1][0]);
    } else {
      int sign = 1;
      for (int i = 0; i < A->rows; i++) {
        matrix_t temp;
        s21_minor(0, i, A, &temp);

        result += sign * A->matrix[0][i] * det(&temp);
        sign *= -1;
        s21_remove_matrix(&temp);
      }
    }
  }
  return result;
}

int s21_minor(int row, int column, matrix_t *A, matrix_t *result) {
  int error_code = INCORRECT_MATRIX;
  if (s21_it_is_normal_matrix(A) && A->rows == A->columns) {
    if (s21_create_matrix(A->rows - 1, A->columns - 1, result) == 0 &&
        row < A->rows && column < A->columns) {
      error_code = OK;
      int n = -1;
      for (int i = 0; i < A->rows; i++) {
        if (i != row) {
          n++;
          int m = -1;
          for (int j = 0; j < A->columns; j++) {
            if (j != column) {
              m++;
              result->matrix[n][m] = A->matrix[i][j];
            }
          }
        }
      }
    } else {
      error_code = CALC_ERROR;
    }
  }
  return error_code;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error_code = INCORRECT_MATRIX;
  if (s21_it_is_normal_matrix(A)) {
    error_code = OK;
    double det = 0;
    s21_determinant(A, &det);
    if (A->rows != A->columns || det == 0) {
      error_code = CALC_ERROR;
    } else {
      matrix_t complements;
      matrix_t transponse;
      int compl = s21_calc_complements(A, &complements);
      int create = s21_transpose(&complements, &transponse);
      if (create == 0 && compl == 0) {
        error_code = s21_mult_number(&transponse, 1 / det, result);
        s21_remove_matrix(&complements);
        s21_remove_matrix(&transponse);
      }
    }
  }
  return error_code;
}

int s21_it_is_normal_matrix(matrix_t *A) {
  int error_code = INCORRECT_MATRIX;
  if (A == NULL) {
    error_code = OK;
  } else {
    if (A->matrix == NULL || A->columns <= 0 || A->rows <= 0) {
      error_code = OK;
    }
  }
  return error_code;
}
