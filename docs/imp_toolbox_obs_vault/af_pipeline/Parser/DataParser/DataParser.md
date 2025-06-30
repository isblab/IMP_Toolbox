```python
class DataParser
    """Class containing methods to parse the AF2/3 data file.

    Attributes:

        data_file_path (str):
            Path to the AF2/3 data file.
    """
```

```mermaid
classDiagram
    class DataParser {
        - __init__(self, data_file_path) None
        + get_data_dict(self) Dict
        + get_token_chain_ids(self, data) list
        + get_token_res_ids(self, data) list
        + get_pae(self, data)
        + get_avg_pae(self, pae)
        + get_contact_probs_mat(self, data)
        + get_avg_contact_probs_mat(self, contact_probs_mat)
    }
```

## Input

- **data_file_path** (`str`)
	- Path to AF3 prediction `json` or `pkl` file

## Attributes

- **data_file_path** (`str`)
	- same as [[#^f5ec42|data_file_path]]

## Methods

- [[get_data_dict]]
- [[get_token_chain_ids]]
- [[get_token_res_ids]]
- [[get_pae]]
- [[get_avg_pae]]
- [[get_contact_probs_mat]]
- [[get_avg_contact_probs_mat]]
- 