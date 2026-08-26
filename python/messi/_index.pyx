# cython: language_level=3

"""Safe Python-facing wrapper around the file-backed MESSI engine."""

from libc.string cimport memset
import os
import tempfile

import numpy as np
cimport numpy as np

from ._native cimport (
    messi_index, messi_index_params, messi_index_create, messi_index_destroy,
    messi_index_add_file, messi_index_search, messi_index_pca_transform,
)

ctypedef np.float32_t FLOAT32_t
ctypedef np.int64_t INT64_t

cdef int DEFAULT_PAA_SEGMENTS = 16
cdef int MESSI_INDEX_ISAX = 0
cdef int MESSI_INDEX_TRIE = 1

_TRANSFORMS = {"sax": 3, "sfa": 4, "spartan": 5, "pisa": 6}
_LAYOUTS = {"isax": MESSI_INDEX_ISAX, "trie": MESSI_INDEX_TRIE}


cdef class Index:
    cdef messi_index* _index
    cdef public int _dim
    cdef public int _n_segments
    cdef int _transform_dim
    cdef int _function_type
    cdef int _filetype_int
    cdef bytes _root_dir
    cdef bint _has_data
    cdef bint _is_norm
    cdef bint _closed
    cdef object _owned_raw_path
    cdef object _config

    def __cinit__(self,
                  int timeseries_size,
                  int n_segments=DEFAULT_PAA_SEGMENTS,
                  int sax_bit_cardinality=8,
                  int max_leaf_size=2000,
                  int min_leaf_size=10,
                  int initial_leaf_buffer_size=2000,
                  int max_total_buffer_size=200000,
                  int initial_fbl_buffer_size=100,
                  int total_loaded_leaves=1,
                  int tight_bound=1,
                  int aggressive_check=0,
                  int function_type=3,
                  transform=None,
                  layout="isax",
                  char simd=1,
                  int sample_size=1000,
                  char is_norm=1,
                  int histogram_type=2,
                  int sample_type=1,
                  sfa_n_coefficients=None,
                  int filetype_int=0,
                  int max_query_threads=1,
                  int queue_count=0,
                  int sampling_seed=1,
                  int node_split_criterion=1,
                  int trie_mbr_dimensions=0,
                  int trie_record_lb_dimensions=0,
                  int trie_split_dimensions=0,
                  bint trie_record_mbr_suffix_bound=True,
                  int trie_leaf_kmeans=0,
                  int trie_fanout=8,
                  bint trie_dynamic_alphabet=False,
                  int trie_min_fanout=2,
                  int trie_max_fanout=16,
                  int trie_alphabet_budget_bits=3,
                  bint dynamic_root_split_variance=False,
                  root_directory=None):
        cdef messi_index_params params
        cdef bytes root_dir_bytes
        cdef const char *root_ptr = <const char *> 0
        cdef int index_type
        cdef int transform_dim
        cdef int record_lb_dim
        cdef int resolved_sfa_n_coefficients
        cdef object layout_name
        cdef object transform_name

        self._index = NULL
        self._has_data = False
        self._closed = False
        self._owned_raw_path = None
        self._config = None
        if timeseries_size <= 0:
            raise ValueError("timeseries_size must be positive")
        if max_query_threads <= 0 or queue_count < 0:
            raise ValueError("max_query_threads must be positive and queue_count non-negative")
        if sampling_seed < 0 or node_split_criterion < 1 or node_split_criterion > 4:
            raise ValueError("invalid sampling_seed or node_split_criterion")
        if not isinstance(layout, str):
            raise TypeError("layout must be 'isax' or 'trie'")
        layout_name = layout.lower()
        if layout_name not in _LAYOUTS:
            raise ValueError("layout must be 'isax' or 'trie'")
        index_type = _LAYOUTS[layout_name]
        if transform is not None:
            if not isinstance(transform, str) or transform.lower() not in _TRANSFORMS:
                raise ValueError("transform must be sax, sfa, spartan, or pisa")
            transform_name = transform.lower()
            function_type = _TRANSFORMS[transform_name]
        elif function_type not in (3, 4, 5, 6):
            raise ValueError("function_type must be 3 (SAX), 4 (SFA), 5 (SPARTAN), or 6 (PISA)")

        if index_type == MESSI_INDEX_TRIE:
            record_lb_dim = trie_record_lb_dimensions or n_segments
            if record_lb_dim < 16 or record_lb_dim > 64:
                raise ValueError("trie_record_lb_dimensions (or n_segments) must be between 16 and 64")
            transform_dim = trie_mbr_dimensions or min(128, timeseries_size)
            if transform_dim < record_lb_dim or transform_dim > 128 or transform_dim > timeseries_size:
                raise ValueError("trie_mbr_dimensions must be between record bound width and min(128, timeseries_size)")
            if trie_split_dimensions == 0:
                trie_split_dimensions = min(32, transform_dim)
            if trie_split_dimensions < 1 or trie_split_dimensions > transform_dim:
                raise ValueError("trie_split_dimensions must be between 1 and trie_mbr_dimensions")
            if trie_leaf_kmeans and (trie_leaf_kmeans < 2 or trie_leaf_kmeans > 64 or function_type not in (4, 5, 6)):
                raise ValueError("trie_leaf_kmeans requires SFA, SPARTAN, or PISA and a value between 2 and 64")
            if trie_dynamic_alphabet:
                if trie_fanout != 8:
                    raise ValueError("trie_fanout cannot be combined with trie_dynamic_alphabet")
                if trie_min_fanout not in (2, 4, 8, 16, 32, 64, 128, 256) or trie_max_fanout not in (2, 4, 8, 16, 32, 64, 128, 256) or trie_max_fanout < trie_min_fanout:
                    raise ValueError("dynamic trie fanouts must be powers of two between 2 and 256")
            elif trie_fanout not in (2, 4, 8):
                raise ValueError("trie_fanout must be 2, 4, or 8")
        else:
            transform_dim = n_segments
            record_lb_dim = 0
            if n_segments <= 0 or n_segments > timeseries_size:
                raise ValueError("n_segments must be between 1 and timeseries_size")
        if dynamic_root_split_variance:
            if index_type != MESSI_INDEX_ISAX:
                raise ValueError("dynamic_root_split_variance requires layout='isax'")
            if function_type not in (4, 5, 6):
                raise ValueError("dynamic_root_split_variance requires SFA, SPARTAN, or PISA")
        if sfa_n_coefficients is None:
            resolved_sfa_n_coefficients = min(64, timeseries_size)
            if resolved_sfa_n_coefficients % 2 != 0:
                resolved_sfa_n_coefficients -= 1
            if index_type == MESSI_INDEX_TRIE and function_type == 4 and \
                    resolved_sfa_n_coefficients < transform_dim:
                resolved_sfa_n_coefficients = transform_dim
        else:
            resolved_sfa_n_coefficients = int(sfa_n_coefficients)
        if function_type == 4 and (resolved_sfa_n_coefficients <= 0 or
                                   resolved_sfa_n_coefficients % 2 != 0 or
                                   resolved_sfa_n_coefficients < transform_dim or
                                   resolved_sfa_n_coefficients > timeseries_size):
            raise ValueError("sfa_n_coefficients must be even and between the transform width and timeseries_size")

        if root_directory is None:
            root_dir_bytes = b""
        elif isinstance(root_directory, bytes):
            root_dir_bytes = root_directory
        else:
            root_dir_bytes = os.fsencode(os.fspath(root_directory))
        self._root_dir = root_dir_bytes
        if self._root_dir:
            root_ptr = <const char *> self._root_dir

        memset(&params, 0, sizeof(messi_index_params))
        params.root_directory = root_ptr
        params.timeseries_size = timeseries_size
        params.n_segments = transform_dim
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
        params.n_coefficients = resolved_sfa_n_coefficients
        params.filetype_int = filetype_int
        params.max_query_threads = max_query_threads
        params.queue_count = queue_count or max_query_threads
        params.index_type = index_type
        params.sampling_seed = sampling_seed or 1
        params.node_split_criterion = node_split_criterion
        params.trie_bound_dimensions = record_lb_dim
        params.trie_split_dimensions = trie_split_dimensions
        params.trie_record_mbr_suffix_bound = trie_record_mbr_suffix_bound
        params.trie_leaf_kmeans = trie_leaf_kmeans
        params.trie_fanout = trie_fanout
        params.trie_dynamic_alphabet = trie_dynamic_alphabet
        params.trie_min_fanout = trie_min_fanout
        params.trie_max_fanout = trie_max_fanout
        params.trie_alphabet_budget_bits = trie_alphabet_budget_bits
        params.dynamic_root_split_variance = dynamic_root_split_variance
        self._index = messi_index_create(&params)
        if self._index is NULL:
            raise MemoryError("Failed to create MESSI index")
        self._dim = timeseries_size
        self._n_segments = record_lb_dim if index_type == MESSI_INDEX_TRIE else transform_dim
        self._transform_dim = transform_dim
        self._function_type = function_type
        self._filetype_int = filetype_int
        self._is_norm = is_norm
        self._config = {"layout": layout_name, "function_type": function_type,
                        "tight_bound": bool(tight_bound),
                        "sfa_n_coefficients": resolved_sfa_n_coefficients if function_type == 4 else None,
                        "transform_dimensions": transform_dim,
                        "record_lb_dimensions": record_lb_dim if index_type == MESSI_INDEX_TRIE else None,
                        "trie_split_dimensions": trie_split_dimensions if index_type == MESSI_INDEX_TRIE else None,
                        "trie_record_mbr_suffix_bound": bool(trie_record_mbr_suffix_bound),
                        "trie_leaf_kmeans": trie_leaf_kmeans,
                        "dynamic_root_split_variance": bool(dynamic_root_split_variance),
                        "max_query_threads": max_query_threads,
                        "queue_count": queue_count or max_query_threads,
                        "sampling_seed": sampling_seed or 1}

    cdef void _remove_owned_raw(self):
        cdef object path = self._owned_raw_path
        self._owned_raw_path = None
        if path is not None:
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass

    def __dealloc__(self):
        if self._index is not NULL:
            messi_index_destroy(self._index)
            self._index = NULL
        try:
            self._remove_owned_raw()
        except Exception:
            pass

    def close(self):
        """Release native resources and remove an API-owned array snapshot."""
        if self._index is not NULL:
            messi_index_destroy(self._index)
            self._index = NULL
        self._has_data = False
        self._closed = True
        self._remove_owned_raw()

    def __enter__(self):
        if self._closed:
            raise RuntimeError("Index is closed")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    cdef void _ensure_buildable(self):
        if self._closed or self._index is NULL:
            raise RuntimeError("Index is closed")
        if self._has_data:
            raise RuntimeError("Index already contains data; MESSI indexes are single-build")

    cdef void _ensure_searchable(self):
        if self._closed or self._index is NULL:
            raise RuntimeError("Index is closed")
        if not self._has_data:
            raise RuntimeError("Index contains no data. Call add_file() or add_array() before search().")

    def add_file(self, filename, ts_num=None, int dynamic_index=1):
        """Build from a caller-owned raw binary dataset."""
        cdef bytes path
        cdef long count
        cdef long available
        cdef long item_size = 1 if self._filetype_int else 4
        cdef long row_bytes = item_size * self._dim
        cdef object file_path
        self._ensure_buildable()
        file_path = os.fspath(filename)
        try:
            available = os.path.getsize(file_path)
        except OSError as exc:
            raise FileNotFoundError(f"Unable to access dataset file: {file_path}") from exc
        if available <= 0 or available % row_bytes != 0:
            raise ValueError("dataset file size is not a whole number of time-series rows")
        if ts_num is None:
            count = available // row_bytes
        else:
            count = int(ts_num)
            if count <= 0 or count > available // row_bytes:
                raise ValueError("ts_num must be positive and no larger than the number of rows in the dataset file")
        path = os.fsencode(file_path)
        if messi_index_add_file(self._index, path, count, dynamic_index) != 0:
            self.close()
            raise RuntimeError("Bulk add failed; the Index was closed")
        self._has_data = True

    def add(self, filename, ts_num=None, int dynamic_index=1):
        """Compatibility alias for :meth:`add_file`."""
        return self.add_file(filename, ts_num, dynamic_index)

    def add_array(self, data, storage_dir=None, int dynamic_index=1):
        """Build from a finite 2-D array via an owned temporary raw snapshot."""
        cdef np.ndarray[FLOAT32_t, ndim=2, mode="c"] array
        cdef object input_array
        cdef int fd = -1
        cdef object path = None
        cdef object directory = None
        self._ensure_buildable()
        input_array = np.asarray(data)
        if input_array.ndim != 2 or input_array.shape[0] == 0 or input_array.shape[1] != self._dim:
            raise ValueError("data must have shape (n_records, timeseries_size) with at least one row")
        if input_array.dtype.kind not in "fiu":
            raise TypeError("data must have a numeric dtype")
        array = np.ascontiguousarray(input_array, dtype=np.float32)
        if not np.isfinite(array).all():
            raise ValueError("data must contain only finite values")
        if storage_dir is not None:
            directory = os.fspath(storage_dir)
            if not os.path.isdir(directory):
                raise ValueError("storage_dir must be an existing directory")
        try:
            fd, path = tempfile.mkstemp(prefix="messi-raw-", suffix=".f32", dir=directory)
            with os.fdopen(fd, "wb") as raw:
                fd = -1
                array.tofile(raw)
                raw.flush()
                os.fsync(raw.fileno())
            self.add_file(path, array.shape[0], dynamic_index)
            self._owned_raw_path = path
        except Exception:
            if fd >= 0:
                os.close(fd)
            if path is not None:
                try:
                    os.unlink(path)
                except FileNotFoundError:
                    pass
            self.close()
            raise

    def search(self, queries, int k=1, int dynamic_index=1):
        """Return exact 1-NN distances and ``None`` for unavailable indices."""
        cdef np.ndarray[FLOAT32_t, ndim=2, mode="c"] query_array
        cdef np.ndarray[FLOAT32_t, ndim=2] distances
        cdef np.ndarray[INT64_t, ndim=2] labels
        cdef Py_ssize_t nq
        cdef object input_array
        self._ensure_searchable()
        if k != 1:
            raise NotImplementedError("exact top-k search is not implemented; only k=1 is supported")
        input_array = np.asarray(queries)
        if input_array.ndim != 2 or input_array.shape[1] != self._dim or input_array.shape[0] == 0:
            raise ValueError("queries must have shape (n_queries, timeseries_size) with at least one row")
        if input_array.dtype.kind not in "fiu":
            raise TypeError("queries must have a numeric dtype")
        query_array = np.ascontiguousarray(input_array, dtype=np.float32)
        if not np.isfinite(query_array).all():
            raise ValueError("queries must contain only finite values")
        nq = query_array.shape[0]
        distances = np.empty((nq, 1), dtype=np.float32)
        labels = np.empty((nq, 1), dtype=np.int64)
        if messi_index_search(self._index, <float*> query_array.data, nq, self._dim, 1,
                              <float*> distances.data, <long*> labels.data, dynamic_index) != 0:
            raise RuntimeError("search failed")
        return distances, None

    def pca_transform(self, queries):
        """Return the learned SPARTAN PCA projection for query rows."""
        cdef np.ndarray[FLOAT32_t, ndim=2, mode="c"] query_array
        cdef np.ndarray[FLOAT32_t, ndim=2] out
        cdef Py_ssize_t nq
        cdef object input_array
        self._ensure_searchable()
        if self._function_type != 5:
            raise RuntimeError("pca_transform is available only for the SPARTAN transform")
        input_array = np.asarray(queries)
        if input_array.ndim != 2 or input_array.shape[1] != self._dim:
            raise ValueError("queries must have shape (n_queries, timeseries_size)")
        if input_array.dtype.kind not in "fiu":
            raise TypeError("queries must have a numeric dtype")
        query_array = np.ascontiguousarray(input_array, dtype=np.float32)
        if not np.isfinite(query_array).all():
            raise ValueError("queries must contain only finite values")
        nq = query_array.shape[0]
        out = np.empty((nq, self._transform_dim), dtype=np.float32)
        if messi_index_pca_transform(self._index, <float*> query_array.data, nq, self._dim,
                                     <float*> out.data, self._transform_dim) != 0:
            raise RuntimeError("pca_transform failed")
        return out

    @property
    def is_norm(self):
        return bool(self._is_norm)

    @property
    def raw_data_path(self):
        return self._owned_raw_path

    @property
    def config(self):
        return dict(self._config)
