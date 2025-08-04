import numpy as np
import pandas as pd
import python_interface

input_array = np.random.uniform(-0.1, 0.1, size=(99989, 128)).astype(np.float32)  # Example data
queries = np.random.uniform(-0.1, 0.1, size=(100, 128)).astype(np.float32)  # Example data
insert_array = np.random.uniform(-0.1, 0.1, size=(10, 128)).astype(np.float32)  # Example data

# Initialize graph
graph = python_interface.initialize_graph_python(100000, 128, 4, 4, 32, 10, 0, "euclidean")

# Build initial graph
python_interface.build_initial_python(graph, input_array, list(range(99989)), 100000, 128, "euclidean")

# Insert nodes
python_interface.insert_nodes_python(graph, insert_array, list(range(99989,99999)), 128)

# # Delete nodes
python_interface.delete_nodes_python(graph, list(range(10)))

# Query
python_interface.test_python(graph, queries, False, "test_python_gt", "test_", "LIGS", 1, "euclidean")