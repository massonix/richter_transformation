# This script joins all the summary.csv of all libraries in a single file


# Load modules
import pandas as pd
import numpy as np
import os



#### HASHED LIBRARIES ####

# Load data
summary_path = "projects/{}/jobs/{}/{}/outs/metrics_summary.csv"
iterable_hashed = []
gem_ids_hashed = os.listdir("projects/BCLLATLAS_10/jobs/")
for gem_id in gem_ids_hashed:
    iterable_hashed.append(["BCLLATLAS_10", gem_id])
print(iterable_hashed)
summary_files_paths_hashed = [summary_path.format(subproject, gem_id, gem_id) for subproject, gem_id in iterable_hashed]
summary_files_dfs_hashed = [pd.read_csv(path) for path in summary_files_paths_hashed]


# Merge
summary_df_hashed = summary_files_dfs_hashed[0]
for i in range(1, len(summary_files_dfs_hashed)):
    summary_df_hashed = summary_df_hashed.append(summary_files_dfs_hashed[i])
summary_df_hashed.insert(0, column = "gem_id", value = gem_ids_hashed)


# Save
summary_df_hashed.to_csv("../results/tables/cellranger_mapping/cellranger_mapping_metrics_hashed.csv", header = True, index = None)



#### NOT HASHED LIBRARIES ####

# Load data
iterable_not_hashed = []
gem_ids_not_hashed = os.listdir("projects/BCLLATLAS_29/jobs/")
for gem_id in gem_ids_not_hashed:
    iterable_not_hashed.append(["BCLLATLAS_29", gem_id])
print(iterable_not_hashed)
summary_files_paths_not_hashed = [summary_path.format(subproject, gem_id, gem_id) for subproject, gem_id in iterable_not_hashed]
summary_files_dfs_not_hashed = [pd.read_csv(path) for path in summary_files_paths_not_hashed]


# Merge
summary_df_not_hashed = summary_files_dfs_not_hashed[0]
for i in range(1, len(summary_files_dfs_not_hashed)):
    summary_df_not_hashed = summary_df_not_hashed.append(summary_files_dfs_not_hashed[i])
summary_df_not_hashed.insert(0, column = "gem_id", value = gem_ids_not_hashed)


# Save
summary_df_not_hashed.to_csv("../results/tables/cellranger_mapping/cellranger_mapping_metrics_not_hashed.csv", header = True, index = None)
