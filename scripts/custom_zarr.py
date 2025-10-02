# The following code is written in Python
import zarr
import numpy as np
from zarr.storage import ZipStore
import sys

from pathlib import Path
from typing import Final, NamedTuple

def getArgs():
    zip_file = sys.argv[1]
    cluster_dir = sys.argv[2].split("/")[-1]
    csv_file = sys.argv[3].split("/")[-1]
    analysis_path = sys.argv[4]
    return zip_file, cluster_dir, csv_file, analysis_path

zip_file, cluster_dir, csv_file, analysis_path = getArgs()

print(zip_file)
print(cluster_dir)
print(csv_file)
print(analysis_path)

OUTPUT_ZIP_FILENAME = zip_file
# Filenames
CELL_GROUPS_ZARR_FILENAME = "cell_groups"

# Common array names
INDICES_ARRAY: Final[str] = "indices"
INDPTR_ARRAY: Final[str] = "indptr"

# Clustering names
custom_clusters: Final[str] = cluster_dir
CLUSTERING_NAMES: Final[list[str]] = [
    custom_clusters,
]

# Attribute names
MAJOR_VERSION_ATTR: Final[str] = "major_version"
MINOR_VERSION_ATTR: Final[str] = "minor_version"
NUMBER_GROUPINGS_ATTR: Final[str] = "number_groupings"
GROUPING_NAMES_ATTR: Final[str] = "grouping_names"
GROUP_NAMES_ATTR: Final[str] = "group_names"


class CellGrouping(NamedTuple):
    name: str
    group_names: list[str]
    indices: np.ndarray[int, np.dtype[np.uint32]]
    indptr: np.ndarray[int, np.dtype[np.uint32]]

    def __eq__(self, other):
        return (
            (self.name == other.name)
            & (self.group_names == other.group_names)
            & (self.indices == other.indices).all()
            & (self.indptr == other.indptr).all()
        )


class CellGroupDataset:
    def __init__(self, cell_groupings: list[CellGrouping]):
        self.cell_groupings = cell_groupings

    def __eq__(self, other):
        return self.cell_groupings == other.cell_groupings

    @staticmethod
    def from_analysis_path(analysis_csv: Path) -> "CellGroupDataset":
        cell_groupings = []

        first_cell_ids = None
        for clustering_name in CLUSTERING_NAMES:

            cell_indices, cluster_assignments = CellGroupDataset.read_clustering(
                analysis_csv, clustering_name
            )


            if first_cell_ids is None:
                first_cell_ids = cell_indices

            assert np.all(cell_indices == first_cell_ids)
            num_cells = cell_indices.shape[0]

            max_cluster_ix = np.max(cluster_assignments)
            num_clusters = max_cluster_ix + 1

            indices = np.zeros(num_cells, dtype=np.uint32)
            indptr = np.zeros(num_clusters, dtype=np.uint32)
            curr_indptr = 0
            for cluster_ix in range(num_clusters):
                indptr[cluster_ix] = curr_indptr
                assigned_to_cluster = cell_indices[cluster_assignments == cluster_ix]
                num_cells_assigned = assigned_to_cluster.shape[0]
                next_indptr = curr_indptr + num_cells_assigned
                indices[curr_indptr:next_indptr] = assigned_to_cluster
                curr_indptr = next_indptr

            group_names = [
                f"Cluster {cluster_ix+1}" for cluster_ix in range(num_clusters)
            ]

            cell_groupings.append(
                CellGrouping(
                    name=clustering_name,
                    group_names=group_names,
                    indices=indices,
                    indptr=indptr,
                )
            )

        return CellGroupDataset(cell_groupings=cell_groupings)

    @staticmethod
    def read_clustering(
        analysis_csv: Path, clustering_type: str
    ) -> tuple[
        np.ndarray[int, np.dtype[np.uint32]], np.ndarray[int, np.dtype[np.uint32]]
    ]:

        csv_path = analysis_csv + "/clustering/" + clustering_type + "/" + csv_file
        with open(csv_path) as csv_handle:
            csv_contents = np.loadtxt(
                csv_handle, delimiter=",", skiprows=1, dtype=np.uint32
            )

        return csv_contents[:, 0] - 1, csv_contents[:, 1]


def store_cell_groups_zarr(cell_groupings: list[CellGrouping]) -> None:
    with zarr.storage.ZipStore(OUTPUT_ZIP_FILENAME, mode='w') as store:
        g = zarr.group(store=store)
        group = g.create_group(CELL_GROUPS_ZARR_FILENAME)

        grouping_names = [grouping.name for grouping in cell_groupings]
        group_names = [grouping.group_names for grouping in cell_groupings]
        number_groupings = len(grouping_names)

        group.attrs.update(
            {
                MAJOR_VERSION_ATTR: 1,
                MINOR_VERSION_ATTR: 0,
                GROUPING_NAMES_ATTR: grouping_names,
                GROUP_NAMES_ATTR: group_names,
                NUMBER_GROUPINGS_ATTR: number_groupings,
            }
        )

        for cell_grouping_ix, cell_grouping in enumerate(cell_groupings):

            subgroup = group.create_group(cell_grouping_ix)

            indices_arr = subgroup.zeros(
                INDICES_ARRAY,
                shape=cell_grouping.indices.shape,
                chunks=cell_grouping.indices.shape,
                dtype=np.uint32,
            )
            indices_arr[:] = cell_grouping.indices

            indptr_arr = subgroup.zeros(
                INDPTR_ARRAY,
                shape=cell_grouping.indptr.shape,
                chunks=cell_grouping.indptr.shape,
                dtype=np.uint32,
            )
            indptr_arr[:] = cell_grouping.indptr

        if isinstance(group.store, ZipStore):
            group.store.close()

dataset = CellGroupDataset.from_analysis_path(analysis_path)

store_cell_groups_zarr(dataset.cell_groupings)