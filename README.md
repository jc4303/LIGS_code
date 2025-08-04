# Repository for the paper [Locality-Sensitive Indexing for Graph-Based Approximate Nearest Neighbor Search](https://dl.acm.org/doi/10.1145/3726302.3730028)

## Build Instructions

```bash
mkdir build
cd build
cmake ..
make -j
```

This will build the stand‑alone `.so` library.

## Python Interface Functions

The following functions are provided:

- **initialize_graph_python** – Allocates the initial data structure for the graph.
- **build_initial_python** – Constructs the initial graph by inserting the starting points. A pure construction workload will only involve this function and the one above.
- **delete_nodes_python** – Deletes points from the graph.
- **insert_nodes_python** – Inserts additional points after the initial graph has been constructed.
- **test_vs_recall_python** – Tests recall vs. QPS.

See **example.py** for a minimal usage example.

## Relation to HNSWlib

This repository includes modified or preserved versions of several source files originally from the [HNSWlib](https://github.com/nmslib/hnswlib/tree/master) repository. LIGS was initially intended as a modification to HNSW, and thus some files are either identical or closely follow the structure of HNSWlib.

Files replicated essentially verbatim for the underlying search logic:

- ligslib/bruteforce.h
- ligslib/hnswlib.h
- ligslib/space_ip.h
- ligslib/space_l2.h
- ligslib/visited_list_pool.h

Files with similar structure to HNSWlib counterparts:

- hnswbin_double_graph.h – Similar to hnswalg.h
- bin_methods_double.cpp – Similar to sift_1b.cpp

These files were copied from an older version of HNSWlib to ensure compatibility with our implementation. They may not reflect the current upstream version and are not automatically updated.

We acknowledge the original authors and license of HNSWlib. The project is available at:  
https://github.com/nmslib/hnswlib

## Citation

If you use this work in your research, please cite:

Jun Woo Chung, Huawei Lin, and Weijie Zhao. *Locality-Sensitive Indexing for Graph-Based Approximate Nearest Neighbor Search*. In SIGIR ’25, Padua, Italy, July 2025. https://doi.org/10.1145/3726302.3730028

BibTeX:

@inproceedings{10.1145/3726302.3730028,
  author    = {Chung, Jun Woo and Lin, Huawei and Zhao, Weijie},
  title     = {Locality-Sensitive Indexing for Graph-Based Approximate Nearest Neighbor Search},
  year      = {2025},
  isbn      = {9798400715921},
  publisher = {Association for Computing Machinery},
  address   = {New York, NY, USA},
  url       = {https://doi.org/10.1145/3726302.3730028},
  doi       = {10.1145/3726302.3730028},
  booktitle = {Proceedings of the 48th International ACM SIGIR Conference on Research and Development in Information Retrieval},
  pages     = {2418--2428},
  numpages  = {11},
  keywords  = {approximate nearest neighbor search, databases, information retrieval, recommender systems},
  location  = {Padua, Italy},
  series    = {SIGIR '25}
}
