"""
utils
===========

- Utility functions for various tasks in IMP_Toolbox.

### Classes:

```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram
    class MatrixPatches {
        - \_\_init__(self, matrix, row_obj, col_obj) None
        + santiy_check_row_col_obj(self)
        + get_patches_from_matrix(self)
        + @staticmethod get_one_sets_from_matrix(matrix, axis)$
        + @classmethod extend_one_sets_by_subsets(cls, one_sets) dict
        + @classmethod split_one_sets(cls, one_sets) dict
        + @staticmethod split_one_set(one_set) list$
        + @staticmethod one_sets_to_df(one_sets, columns) pd.DataFrame$
        + @staticmethod aggregate_df_rows(df, groupby_col, agg_col) pd.DataFrame$
        + @staticmethod combine_dfs(df1, df2, colname_1, colname_2) pd.DataFrame$
        + remove_subset_rows(self, df, colname_1, colname_2)
    }
```
"""