def csr_add_sparse_vec(sps_mat, sps_vec) :
    """Adds a sparse vector to every row of a sparse matrix"""
    # No checks done, but both arguments should be sparse matrices in CSR
    # format, both should have the same number of columns, and the vector
    # should be a vector and have only one row.

    rows, cols = sps_mat.shape
    nnz_vec = len(sps_vec.data)
    nnz_per_row = np.diff(sps_mat.indptr)
    longest_row = np.max(nnz_per_row)

    old_data = np.zeros((rows * longest_row,), dtype=sps_mat.data.dtype)
    old_cols = np.zeros((rows * longest_row,), dtype=sps_mat.indices.dtype)

    data_idx = np.arange(longest_row) < nnz_per_row[:, None]
    data_idx = data_idx.reshape(-1)
    old_data[data_idx] = sps_mat.data
    old_cols[data_idx] = sps_mat.indices
    old_data = old_data.reshape(rows, -1)
    old_cols = old_cols.reshape(rows, -1)

    new_data = np.zeros((rows, longest_row + nnz_vec,),
                        dtype=sps_mat.data.dtype)
    new_data[:, :longest_row] = old_data
    del old_data
    new_cols = np.zeros((rows, longest_row + nnz_vec,),
                        dtype=sps_mat.indices.dtype)
    new_cols[:, :longest_row] = old_cols
    del old_cols
    new_data[:, longest_row:] = sps_vec.data
    new_cols[:, longest_row:] = sps_vec.indices
    new_data = new_data.reshape(-1)
    new_cols = new_cols.reshape(-1)
    new_pointer = np.arange(0, (rows + 1) * (longest_row + nnz_vec),
                            longest_row + nnz_vec)

    ret = sps.csr_matrix((new_data, new_cols, new_pointer),
                         shape=sps_mat.shape)
    ret.eliminate_zeros()

    return ret
   



#                unions = np.zeros((Nsam,1))
#                for j in range(self.Nsam):
#                    suma = self.X[i,:] + self.X[j,:] 
#                    unions[j] = suma.size
         
         #                print unions.shape

# Anpther way to do it.
#                suma = csr_add_sparse_vec(self.X[:,:],self.X[i,:])
#                unions = np.array(suma.getnnz(axis = 1), dtype= float) - inters
#                suma = csr_add_sparse_vec(self.X[:,:],self.X[i,:])
#                unions = np.array(suma.getnnz(axis = 1), dtype= float) - inters


