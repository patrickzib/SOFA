# cython: language_level=3

"""Safe Python-facing wrapper around the native S3-Trie/MESSI engine."""

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

cdef int DEFAULT_ISAX_SEGMENTS = 16
cdef int DEFAULT_TRIE_RECORD_LB_SEGMENTS = 64
cdef int MESSI_INDEX_ISAX = 0
cdef int MESSI_INDEX_TRIE = 1

_TRANSFORMS = {"sax": 3, "sfa": 4, "spartan": 5, "pisa": 6}
_LAYOUTS = {"isax": MESSI_INDEX_ISAX, "trie": MESSI_INDEX_TRIE}
_PRESETS = {
    "messi": ("sax", "isax"),
    "sofa": ("sfa", "isax"),
    "s3trie": ("spartan", "trie"),
}


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
                  n_segments=None,
                  int sax_bit_cardinality=8,
                  max_leaf_size=None,
                  min_leaf_size=None,
                  initial_leaf_buffer_size=None,
                  int max_total_buffer_size=200000,
                  int initial_fbl_buffer_size=100,
                  int total_loaded_leaves=1,
                  int tight_bound=1,
                  int aggressive_check=0,
                  function_type=None,
                  index=None,
                  transform=None,
                  layout=None,
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
                  trie_record_mbr_suffix_bound=None,
                  trie_leaf_ivf=None,
                  trie_leaf_ivf_min_size=None,
                  trie_leaf_ivf_raw_ball_bound=None,
                  trie_leaf_ivf_radial_bound=None,
                  bint trie_leaf_ivf_radial_bound_auto=False,
                  int trie_fanout=8,
                  bint trie_dynamic_alphabet=False,
                  int trie_min_fanout=2,
                  int trie_max_fanout=16,
                  int trie_alphabet_budget_bits=3,
                  bint dynamic_root_split_variance=False,
                  root_directory=None,
                  trie_streaming_leaf_scan=None,
                  trie_residual_record_only=None,
                  trie_residual_order=None):
        cdef messi_index_params params
        cdef bytes root_dir_bytes
        cdef const char *root_ptr = <const char *> 0
        cdef int index_type
        cdef int transform_dim
        cdef int record_lb_dim
        cdef int resolved_n_segments
        cdef int resolved_sfa_n_coefficients
        cdef int resolved_max_leaf_size
        cdef int resolved_min_leaf_size
        cdef int resolved_initial_leaf_buffer_size
        cdef int resolved_trie_leaf_ivf
        cdef int resolved_trie_leaf_ivf_min_size
        cdef bint resolved_trie_streaming_leaf_scan
        cdef bint resolved_trie_record_mbr_suffix_bound
        cdef bint resolved_trie_leaf_ivf_raw_ball_bound
        cdef bint resolved_trie_leaf_ivf_radial_bound
        cdef bint resolved_trie_residual_record_only
        cdef int resolved_trie_residual_order
        cdef object preset_name
        cdef object preset_transform
        cdef object preset_layout
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

        # The high-level Python entry point defaults to the paper-oriented
        # S³-Trie preset. Explicit layout/transform/function_type arguments
        # retain the lower-level behavior and bypass the preset.
        if (index is None and layout is None and transform is None and
                function_type is None):
            index = "s3trie"

        preset_name = None
        if index is not None:
            if not isinstance(index, str):
                raise TypeError("index must be 'messi', 'sofa', or 's3trie'")
            preset_name = index.lower()
            if preset_name not in _PRESETS:
                raise ValueError("index must be 'messi', 'sofa', or 's3trie'")
            preset_transform, preset_layout = _PRESETS[preset_name]
            if layout is not None and (not isinstance(layout, str) or layout.lower() != preset_layout):
                raise ValueError(f"index='{preset_name}' conflicts with layout={layout!r}")
            if transform is not None and (not isinstance(transform, str) or transform.lower() != preset_transform):
                raise ValueError(f"index='{preset_name}' conflicts with transform={transform!r}")
            if function_type is not None and int(function_type) != _TRANSFORMS[preset_transform]:
                raise ValueError(f"index='{preset_name}' conflicts with function_type={function_type!r}")
            layout_name = preset_layout
            transform_name = preset_transform
            function_type = _TRANSFORMS[transform_name]
        else:
            layout_name = "isax" if layout is None else layout
            if not isinstance(layout_name, str):
                raise TypeError("layout must be 'isax' or 'trie'")
            layout_name = layout_name.lower()
            if transform is not None:
                if not isinstance(transform, str) or transform.lower() not in _TRANSFORMS:
                    raise ValueError("transform must be sax, sfa, spartan, or pisa")
                transform_name = transform.lower()
                if function_type is not None and int(function_type) != _TRANSFORMS[transform_name]:
                    raise ValueError("transform conflicts with function_type")
                function_type = _TRANSFORMS[transform_name]
            else:
                if function_type is None:
                    function_type = 3
                    if layout is None:
                        preset_name = "messi"
                function_type = int(function_type)
                if function_type not in (3, 4, 5, 6):
                    raise ValueError("function_type must be 3 (SAX), 4 (SFA), 5 (SPARTAN), or 6 (PISA)")
                transform_name = {3: "sax", 4: "sfa", 5: "spartan", 6: "pisa"}[function_type]

        if layout_name not in _LAYOUTS:
            raise ValueError("layout must be 'isax' or 'trie'")
        index_type = _LAYOUTS[layout_name]
        resolved_max_leaf_size = int(20000 if preset_name == "s3trie" and max_leaf_size is None
                                     else 2000 if max_leaf_size is None else max_leaf_size)
        resolved_min_leaf_size = int(20000 if preset_name == "s3trie" and min_leaf_size is None
                                     else 10 if min_leaf_size is None else min_leaf_size)
        resolved_initial_leaf_buffer_size = int(
            20000 if preset_name == "s3trie" and initial_leaf_buffer_size is None
            else 2000 if initial_leaf_buffer_size is None else initial_leaf_buffer_size)
        if (resolved_max_leaf_size <= 0 or resolved_min_leaf_size <= 0 or
                resolved_min_leaf_size > resolved_max_leaf_size or
                resolved_initial_leaf_buffer_size <= 0):
            raise ValueError("leaf sizes must be positive and min_leaf_size cannot exceed max_leaf_size")
        if n_segments is None:
            resolved_n_segments = (
                DEFAULT_TRIE_RECORD_LB_SEGMENTS
                if index_type == MESSI_INDEX_TRIE
                else DEFAULT_ISAX_SEGMENTS
            )
        else:
            resolved_n_segments = n_segments
        resolved_trie_streaming_leaf_scan = (index_type == MESSI_INDEX_TRIE
            if trie_streaming_leaf_scan is None else bool(trie_streaming_leaf_scan))
        resolved_trie_record_mbr_suffix_bound = (index_type == MESSI_INDEX_TRIE
            if trie_record_mbr_suffix_bound is None else bool(trie_record_mbr_suffix_bound))
        resolved_trie_leaf_ivf = int(16 if preset_name == "s3trie" and trie_leaf_ivf is None
                                     else 0 if trie_leaf_ivf is None else trie_leaf_ivf)
        resolved_trie_leaf_ivf_min_size = int(
            4096 if trie_leaf_ivf_min_size is None else trie_leaf_ivf_min_size)
        resolved_trie_leaf_ivf_raw_ball_bound = (
            index_type == MESSI_INDEX_TRIE and resolved_trie_leaf_ivf != 0
            if trie_leaf_ivf_raw_ball_bound is None else bool(trie_leaf_ivf_raw_ball_bound))
        resolved_trie_leaf_ivf_radial_bound = (
            preset_name == "s3trie" and resolved_trie_leaf_ivf != 0
            if trie_leaf_ivf_radial_bound is None else bool(trie_leaf_ivf_radial_bound))
        resolved_trie_residual_record_only = (preset_name == "s3trie"
            if trie_residual_record_only is None else bool(trie_residual_record_only))
        if trie_residual_order is None:
            resolved_trie_residual_order = 0
        elif trie_residual_order == "symbolic-first":
            resolved_trie_residual_order = 0
        elif trie_residual_order == "residual-first":
            resolved_trie_residual_order = 1
        else:
            raise ValueError("trie_residual_order must be 'symbolic-first' or 'residual-first'")

        if index_type == MESSI_INDEX_TRIE:
            record_lb_dim = trie_record_lb_dimensions or resolved_n_segments
            if record_lb_dim < 16 or record_lb_dim > 64:
                raise ValueError("trie_record_lb_dimensions (or n_segments) must be between 16 and 64")
            transform_dim = trie_mbr_dimensions or min(128, timeseries_size)
            if transform_dim < record_lb_dim or transform_dim > 128 or transform_dim > timeseries_size:
                raise ValueError("trie_mbr_dimensions must be between record bound width and min(128, timeseries_size)")
            if trie_split_dimensions == 0:
                trie_split_dimensions = max(record_lb_dim, min(32, transform_dim))
            if trie_split_dimensions < record_lb_dim or trie_split_dimensions > transform_dim:
                raise ValueError("trie_split_dimensions must be between the record bound width and trie_mbr_dimensions")
            if resolved_trie_leaf_ivf and (resolved_trie_leaf_ivf < 2 or resolved_trie_leaf_ivf > 64 or function_type not in (4, 5, 6)):
                raise ValueError("trie_leaf_ivf requires SFA, SPARTAN, or PISA and a value between 2 and 64")
            if resolved_trie_leaf_ivf_min_size <= 0:
                raise ValueError("trie_leaf_ivf_min_size must be positive")
            if (resolved_trie_leaf_ivf_radial_bound or trie_leaf_ivf_radial_bound_auto) and not resolved_trie_leaf_ivf:
                raise ValueError("trie radial bounds require trie_leaf_ivf")
            if resolved_trie_residual_record_only and function_type != 5:
                raise ValueError("trie_residual_record_only requires the SPARTAN transform")
            if trie_dynamic_alphabet:
                if trie_fanout != 8:
                    raise ValueError("trie_fanout cannot be combined with trie_dynamic_alphabet")
                if trie_min_fanout not in (2, 4, 8, 16, 32, 64, 128, 256) or trie_max_fanout not in (2, 4, 8, 16, 32, 64, 128, 256) or trie_max_fanout < trie_min_fanout:
                    raise ValueError("dynamic trie fanouts must be powers of two between 2 and 256")
            elif trie_fanout not in (2, 4, 8):
                raise ValueError("trie_fanout must be 2, 4, or 8")
        else:
            if resolved_trie_streaming_leaf_scan:
                raise ValueError("trie_streaming_leaf_scan requires layout='trie'")
            if resolved_trie_leaf_ivf or resolved_trie_leaf_ivf_radial_bound or trie_leaf_ivf_radial_bound_auto:
                raise ValueError("trie radial bounds require layout='trie'")
            if resolved_trie_residual_record_only:
                raise ValueError("trie_residual_record_only requires layout='trie'")
            transform_dim = resolved_n_segments
            record_lb_dim = 0
            if resolved_n_segments <= 0 or resolved_n_segments > timeseries_size:
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
        params.max_leaf_size = resolved_max_leaf_size
        params.min_leaf_size = resolved_min_leaf_size
        params.initial_leaf_buffer_size = resolved_initial_leaf_buffer_size
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
        params.trie_record_mbr_suffix_bound = resolved_trie_record_mbr_suffix_bound
        params.trie_streaming_leaf_scan = resolved_trie_streaming_leaf_scan
        params.trie_leaf_ivf = resolved_trie_leaf_ivf
        params.trie_leaf_ivf_min_size = resolved_trie_leaf_ivf_min_size
        params.trie_leaf_ivf_raw_ball_bound = resolved_trie_leaf_ivf_raw_ball_bound
        params.trie_leaf_ivf_raw_ball_bound_specified = 1
        params.trie_leaf_ivf_radial_bound = resolved_trie_leaf_ivf_radial_bound or trie_leaf_ivf_radial_bound_auto
        params.trie_leaf_ivf_radial_bound_auto = trie_leaf_ivf_radial_bound_auto
        params.trie_fanout = trie_fanout
        params.trie_dynamic_alphabet = trie_dynamic_alphabet
        params.trie_min_fanout = trie_min_fanout
        params.trie_max_fanout = trie_max_fanout
        params.trie_alphabet_budget_bits = trie_alphabet_budget_bits
        params.dynamic_root_split_variance = dynamic_root_split_variance
        params.trie_residual_norm_bound = resolved_trie_residual_record_only
        params.trie_residual_record_only = resolved_trie_residual_record_only
        params.trie_residual_order = resolved_trie_residual_order
        self._index = messi_index_create(&params)
        if self._index is NULL:
            raise MemoryError("Failed to create native index")
        self._dim = timeseries_size
        self._n_segments = record_lb_dim if index_type == MESSI_INDEX_TRIE else transform_dim
        self._transform_dim = transform_dim
        self._function_type = function_type
        self._filetype_int = filetype_int
        self._is_norm = is_norm
        self._config = {"index": preset_name,
                        "layout": layout_name,
                        "transform": transform_name,
                        "function_type": function_type,
                        "max_leaf_size": resolved_max_leaf_size,
                        "min_leaf_size": resolved_min_leaf_size,
                        "initial_leaf_buffer_size": resolved_initial_leaf_buffer_size,
                        "tight_bound": bool(tight_bound),
                        "histogram_type": histogram_type,
                        "sfa_n_coefficients": resolved_sfa_n_coefficients if function_type == 4 else None,
                        "transform_dimensions": transform_dim,
                        "record_lb_dimensions": record_lb_dim if index_type == MESSI_INDEX_TRIE else None,
                        "trie_split_dimensions": trie_split_dimensions if index_type == MESSI_INDEX_TRIE else None,
                        "trie_record_mbr_suffix_bound": bool(resolved_trie_record_mbr_suffix_bound),
                        "trie_streaming_leaf_scan": bool(resolved_trie_streaming_leaf_scan),
                        "trie_leaf_ivf": resolved_trie_leaf_ivf,
                        "trie_leaf_ivf_min_size": resolved_trie_leaf_ivf_min_size,
                        "trie_leaf_ivf_raw_ball_bound": bool(resolved_trie_leaf_ivf_raw_ball_bound),
                        "trie_leaf_ivf_radial_bound": bool(resolved_trie_leaf_ivf_radial_bound or trie_leaf_ivf_radial_bound_auto),
                        "trie_leaf_ivf_radial_bound_auto": bool(trie_leaf_ivf_radial_bound_auto),
                        "trie_residual_record_only": bool(resolved_trie_residual_record_only),
                        "trie_residual_order": "residual-first" if resolved_trie_residual_order else "symbolic-first",
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
            raise RuntimeError("Index already contains data; indexes are single-build")

    cdef void _ensure_searchable(self):
        if self._closed or self._index is NULL:
            raise RuntimeError("Index is closed")
        if not self._has_data:
            raise RuntimeError("Index contains no data. Call add_file() or add() before search().")

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

    def add(self, data, storage_dir=None, int dynamic_index=1):
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
            fd, path = tempfile.mkstemp(prefix="s3trie-raw-", suffix=".f32", dir=directory)
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
        """Return exact 1-NN distances and zero-based sequential row IDs."""
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
        return distances, labels

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
