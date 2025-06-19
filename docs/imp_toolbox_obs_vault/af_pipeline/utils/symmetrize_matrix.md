```python
def symmetrize_matrix(matrix: np.ndarray):
    """Symmetrize a matrix by averaging it with its transpose

    Args:
        matrix (np.ndarray): input matrix

    Returns:
        sym_matrix (np.ndarray): symmetrized matrix
    """

    assert isinstance(matrix, np.ndarray), "Input must be a numpy array"
    assert matrix.ndim == 2, "Input must be a 2D array"

    sym_matrix = (matrix + matrix.T) / 2

    return sym_matrix
```