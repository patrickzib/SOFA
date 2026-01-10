# cython: language_level=3

from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
import numpy as np
cimport numpy as np

from ._native cimport (
    messi_index,
    messi_index_params,
    messi_index_create,
    messi_index_destroy,
    messi_index_add_file,
    messi_index_search,
)

ctypedef np.float32_t FLOAT32_t
ctypedef np.int64_t INT64_t

cdef int DEFAULT_PAA_SEGMENTS = 16


cdef class Index:
    cdef messi_index* _index
    cdef public int _dim
    cdef bytes _root_dir
    cdef bint _has_data
    cdef bint _is_norm

    def __cinit__(self,
                  int timeseries_size,
                  int sax_bit_cardinality=8,
                  int max_leaf_size=2000,
                  int min_leaf_size=10,
                  int initial_leaf_buffer_size=2000,
                  int max_total_buffer_size=200000,
                  int initial_fbl_buffer_size=100,
                  int total_loaded_leaves=1,
                  int tight_bound=0,
                  int aggressive_check=0,
                  int function_type=3,
                  char simd=1,
                  int sample_size=1000,
                  char is_norm=1,
                  int histogram_type=2,
                  int sample_type=1,
                  int sfa_n_coefficients=32,
                  int filetype_int=0,
                  int max_query_threads=1,
                  root_directory=None):
        cdef messi_index_params params
        cdef bytes root_dir_bytes
        params.timeseries_size = timeseries_size
        params.n_segments = DEFAULT_PAA_SEGMENTS
        params.sax_bit_cardinality = sax_bit_cardinality
        params.max_leaf_size = max_leaf_size
        params.min_leaf_size = min_leaf_size
        params.initial_leaf_buffer_size = initial_leaf_buffer_size
        params.max_total_buffer_size = max_total_buffer_size
        params.initial_fbl_buffer_size = initial_fbl_buffer_size
        params.total_loaded_leaves = total_loaded_leaves
        params.tight_bound = tight_bound
        params.aggressive_check = aggressive_check
        params.function_type = function_type
        params.simd = simd
        params.sample_size = sample_size
        params.is_norm = is_norm
        params.histogram_type = histogram_type
        params.sample_type = sample_type
        params.n_coefficients = sfa_n_coefficients
        params.filetype_int = filetype_int
        params.max_query_threads = max_query_threads
        params.queue_count = max_query_threads
        cdef const char *root_ptr = <const char *> 0

        if root_directory is None:
            root_dir_bytes = b""
        elif isinstance(root_directory, bytes):
            root_dir_bytes = root_directory
        else:
            root_dir_bytes = (<str> root_directory).encode("utf-8")
        self._root_dir = root_dir_bytes
        if self._root_dir:
            root_ptr = <const char *> self._root_dir
        params.root_directory = root_ptr

        self._index = messi_index_create(&params)
        if self._index is NULL:
            raise MemoryError("Failed to create MESSI index")
        self._dim = timeseries_size
        self._has_data = False
        self._is_norm = is_norm

    def __dealloc__(self):
        if self._index is not NULL:
            if self._has_data:
                messi_index_destroy(self._index)
            else:
                # Skip destroying partially initialized native state to avoid crashes.
                pass
            self._index = NULL

    def add(self, filename: str, long ts_num):
        cdef bytes path = filename.encode("utf-8")
        if messi_index_add_file(self._index, path, ts_num) != 0:
            raise RuntimeError("Bulk add failed")
        self._has_data = True

    def search(self, np.ndarray[FLOAT32_t, ndim=2] queries, int k):
        if not self._has_data:
            raise RuntimeError("Index contains no data. Call add() before search().")
        if queries.dtype != np.float32:
            queries = np.asarray(queries, dtype=np.float32)
        if queries.ndim != 2:
            raise ValueError("queries must be 2-D")
        if queries.shape[1] != self._dim:
            raise ValueError("dimension mismatch")
        cdef Py_ssize_t nq = queries.shape[0]
        cdef np.ndarray[FLOAT32_t, ndim=2] distances = np.empty((nq, k), dtype=np.float32)
        cdef np.ndarray[INT64_t, ndim=2] labels = np.empty((nq, k), dtype=np.int64)
        cdef float* q_ptr = <float*> queries.data
        cdef float* d_ptr = <float*> distances.data
        cdef long* l_ptr = <long*> labels.data

        if messi_index_search(self._index,
                               q_ptr,
                               nq,
                               self._dim,
                               k,
                               d_ptr,
                               l_ptr) != 0:
            raise RuntimeError("search failed")

        return distances, labels

    @property
    def is_norm(self):
        return bool(self._is_norm)
