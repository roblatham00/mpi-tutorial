/* 
 * A simple program to read a sparse matrix, stored in 
 * compressed-sparse-row (CSR) format.
 *
 * CSR format is defined by 3 arrays and one integer:
 * n - number of rows in the matrix (integer)
 * ia[n+1] - index of data for each row.  Last (n+1st) element gives the 
 * total length of the arrays a and ja (integer)
 * a[nz]   - values of matrix entries (double)
 * ja[nz]  - column number for corresponding a entry
 * All index values are 1-origin (this format was developed for Fortran)
 *
 * For example, the sparse matrix
 *      ( 1  2  0
 *        3  0  0
 *        0  4  5 )
 * is stored as:
 * n = 3
 * ia = ( 1, 3, 4, 6 )
 * a  = ( 1, 2, 3, 4, 5 )
 * ja = ( 1, 2, 1, 2, 3 )
 *
 * The number of elements on the ith row is ia[i+1] i ia[i] (which is why ia
 * and n+1, not just n entries).
 *
 * A file format such as the Boeing-Harwell exchange format will often 
 * include the number of non-zeros separately from the ia[n+1] value.
 *
 * For the purposes of our example, we will store the data in binary (native)
 * representation (later, we can look at external32)
 *
 * title (char, 80 bytes, fixed size)
 * n (int)
 * nz (int)
 * ia[i], i=1,...,n+1 (int array)
 * ja[i], i=1,...,nz  (int array)
 * a[i],  i=1,...,nz  (double array)
 *
 * A simple parallel reader does the following
 *  read_all( title, n, nz )   - All processes read the header
 *  compute location of entries to read:
 *    row range = (n/comm_size) * comm_rank to (n/comm_size) * (comm_rank + 1)-1 *  read row range of ia into a local ia array (ialocal)
 *  compute elements of file to read based on the indices in ialocal:
 *    JAOffset = sizeofHeader + (n+1) * sizeof(int) + (ialocal[o]-1) * sizeof(int)
 *    AOffset = sizeofHeader + (n+1) * sizeof(int) + nz * sizeof(int) + 
 *                (ialocal[o]-1) * sizeof(double)
 *    nelements = ialocal[lastLocalRow+1] - ialocal[0]
 *    read nelements ints at JAOffset
 *    read nelements doubles at AOffset
 */
