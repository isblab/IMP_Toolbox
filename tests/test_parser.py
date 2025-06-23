import pytest
from af_pipeline.parser import AfParser
import numpy as np


@pytest.fixture
def parser0():
    return AfParser(
        data_file_path="./tests/data/af_predictions/dummy_files/dummy_data_afdb.json"
    )

@pytest.fixture
def parser1():
    return AfParser(
        data_file_path="./tests/data/af_predictions/afdb/AF-Q49B96-F1-predicted_aligned_error_v4.json",
        structure_file_path="./tests/data/af_predictions/afdb/AF-Q49B96-F1-model_v4.pdb",
    )

@pytest.fixture
def parser2():
    return AfParser(
        data_file_path="./tests/data/af_predictions/dummy_files/dummy_data_af3.json"
    )


def test_create_mask():
    lengths_dict = {"A": 3, "B": 2, "C": 2, "total": 7}
    hide_interactions_ = ["intrachain", "interchain"]
    masked_value_ = [-100, 1]
    unmasked_value_ = [1, 0]

    expected_result1 = np.array(
        [
            [-100, -100, -100, 1, 1, 1, 1],
            [-100, -100, -100, 1, 1, 1, 1],
            [-100, -100, -100, 1, 1, 1, 1],
            [1, 1, 1, -100, -100, 1, 1],
            [1, 1, 1, -100, -100, 1, 1],
            [1, 1, 1, 1, 1, -100, -100],
            [1, 1, 1, 1, 1, -100, -100],
        ]
    )
    result1 = AfParser.create_mask(
        lengths_dict,
        hide_interactions_[0],
        masked_value_[0],
        unmasked_value_[0],
    )
    np.testing.assert_array_equal(result1, expected_result1)

    expected_result2 = np.array(
        [
            [1, 1, 1, -100, -100, -100, -100],
            [1, 1, 1, -100, -100, -100, -100],
            [1, 1, 1, -100, -100, -100, -100],
            [-100, -100, -100, 1, 1, -100, -100],
            [-100, -100, -100, 1, 1, -100, -100],
            [-100, -100, -100, -100, -100, 1, 1],
            [-100, -100, -100, -100, -100, 1, 1],
        ]
    )
    result2 = AfParser.create_mask(
        lengths_dict,
        hide_interactions_[1],
        masked_value_[0],
        unmasked_value_[0],
    )
    np.testing.assert_array_equal(result2, expected_result2)

    expected_result3 = np.array(
        [
            [1, 1, 1, 0, 0, 0, 0],
            [1, 1, 1, 0, 0, 0, 0],
            [1, 1, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 1, 0, 0],
            [0, 0, 0, 1, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 1],
            [0, 0, 0, 0, 0, 1, 1],
        ]
    )
    result3 = AfParser.create_mask(
        lengths_dict,
        hide_interactions_[0],
        masked_value_[1],
        unmasked_value_[1],
    )
    np.testing.assert_array_equal(result3, expected_result3)

    expected_result4 = np.array([
        [0, 0, 0, 1, 1, 1, 1],
        [0, 0, 0, 1, 1, 1, 1],
        [0, 0, 0, 1, 1, 1, 1],
        [1, 1, 1, 0, 0, 1, 1],
        [1, 1, 1, 0, 0, 1, 1],
        [1, 1, 1, 1, 1, 0, 0],
        [1, 1, 1, 1, 1, 0, 0],
    ])
    result4 = AfParser.create_mask(
        lengths_dict,
        hide_interactions_[1],
        masked_value_[1],
        unmasked_value_[1],
    )
    np.testing.assert_array_equal(result4, expected_result4)

def test_get_min_pae():
    pae_matrix = np.array([
        [0.0, 1.2, 2.3, 8.0,],
        [1.2, 0.0, 1.5, 7.5,],
        [1.3, 1.5, 0.0, 1.0,],
        [8.0, 7.0, 6.0, 0.0,],
    ])
    lengths_dict = {"A": 2, "B": 1, "C": 1, "total": 4}
    along_axis_ = [1, 0]
    hide_interactions_ = ["intrachain", "interchain"]
    return_type_ = ["list", "array", "dict"]

    expected_result1 = [2.3, 1.5, 1.0, 6.0]
    result1 = AfParser.get_min_pae(
        pae_matrix,
        lengths_dict,
        along_axis_[0],
        hide_interactions_[0],
        return_type_[0],
    )
    assert isinstance(result1, list), \
        f"result should be a list, got {type(result1)}"
    np.testing.assert_array_equal(result1, expected_result1)

    expected_result2 = [1.3, 1.5, 1.5, 1.0]
    result2 = AfParser.get_min_pae(
        pae_matrix,
        lengths_dict,
        along_axis_[1],
        hide_interactions_[0],
        return_type_[0],
    )
    assert isinstance(result2, list), \
        f"result should be a list, got {type(result1)}"
    np.testing.assert_array_equal(result2, expected_result2)

    expected_result3 = np.array([0.0, 0.0, 0.0, 0.0])
    result3 = AfParser.get_min_pae(
        pae_matrix,
        lengths_dict,
        along_axis_[0],
        hide_interactions_[1],
        return_type_[1],
    )
    assert isinstance(result3, type(expected_result3)), \
        f"result should be a numpy array, got {type(result3)}"
    np.testing.assert_array_equal(result3, expected_result3)

    expected_result4 = {
        "A": [2.3, 1.5],
        "B": [1.0],
        "C": [6.0],
    }
    result4 = AfParser.get_min_pae(
        pae_matrix,
        lengths_dict,
        along_axis_[0],
        hide_interactions_[0],
        return_type_[2],
    )
    for chain_id, min_pae in result4.items():
        assert isinstance(chain_id, str), "chain_id should be a string"
        assert isinstance(min_pae, list), "min_pae should be a list"
        np.testing.assert_array_equal(
            min_pae, expected_result4[chain_id]
        ), f"min_pae is incorrect for chain {chain_id}"

def test_get_chain_lengths():
    token_chain_ids = ["A", "A", "B", "C"]

    expected_result = {"A": 2, "B": 1, "C": 1, "total": 4}
    results = AfParser.get_chain_lengths(token_chain_ids)
    assert isinstance(results, dict), "results should be a dictionary"
    for chain_id, ch_len in results.items():
        assert isinstance(chain_id, str), "chain_id should be a string"
        assert isinstance(ch_len, int), "chain length should be an integer"
        assert ch_len >= 0, "chain length should be non-negative"
        assert ch_len == expected_result[chain_id], \
            f"Check if length for chain {chain_id} is correct: expected {expected_result[chain_id]}, got {ch_len}"

def test_get_data_dict0(parser0: AfParser):
    result = parser0.dataparser.get_data_dict()
    assert isinstance(result, dict), "data_dict should be a dictionary"
    assert ("predicted_aligned_error" in result), \
        "data_dict should contain 'predicted_aligned_error' key"

def test_get_data_dict1(parser2: AfParser):
    result = parser2.dataparser.get_data_dict()
    assert isinstance(result, dict), "data_dict should be a dictionary"
    assert "contact_probs" in result, \
        "data_dict should contain 'contact_probs' key"
    assert "pae" in result, "data_dict should contain 'pae' key"
    assert "token_chain_ids" in result, \
        "data_dict should contain 'token_chain_ids' key"
    assert "token_res_ids" in result, \
        "data_dict should contain 'token_res_ids' key"

#TODO
# def test_update_token_ids
# def test_update_pae
# def test_update_contact_probs
# def test_extract_perresidue_quantity
# def test renumber_structure
# def test_renumber_chain_res_num
# def test_renumber_region_of_interest
# def test_residue_map

