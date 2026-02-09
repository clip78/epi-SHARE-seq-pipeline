import os
import tempfile
import unittest

import pandas as pd

from src.python.plot_insert_size_hist import get_hist_vals

class TestGetHistVals(unittest.TestCase):

    def test_get_hist_vals(self):
        # Create a sample histogram file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_file:
            temp_file.write("# A sample Picard CollectInsertSizeMetrics histogram file\n")
            temp_file.write("# This file is used for testing the plot_insert_size_hist.py script.\n")
            temp_file.write("\n")
            temp_file.write("insert_size\tAll_Reads.fr_count\n")
            temp_file.write("10\t100\n")
            temp_file.write("20\t200\n")
            temp_file.write("30\t300\n")
            histogram_path = temp_file.name

        try:
            # Call the function with the test file
            df = get_hist_vals(histogram_path)
        finally:
            os.unlink(histogram_path)

        # Check the output dataframe
        expected_df = pd.DataFrame({
            "insert_size": [10, 20, 30],
            "count": [100, 200, 300]
        })
        pd.testing.assert_frame_equal(df, expected_df)

if __name__ == '__main__':
    unittest.main()
