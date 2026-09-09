use nalgebra::DMatrix;

fn get_a_submatrix(a: DMatrix<f64>, indices: &[usize]) -> DMatrix<f64> {
    let mut submatrix = DMatrix::zeros(indices.len(), a.ncols());
    for (i, &row) in indices.iter().enumerate() {
        submatrix.row_mut(i).copy_from(&a.row(row));
    }
    submatrix
}

pub fn get_6dof_lateral(a: DMatrix<f64>) -> DMatrix<f64> {
    get_a_submatrix(a, &[0, 1, 4, 7])
}

pub fn get_6dof_longitudal(a: DMatrix<f64>) -> DMatrix<f64> {
    get_a_submatrix(a, &[2, 3, 6, 8])
}
