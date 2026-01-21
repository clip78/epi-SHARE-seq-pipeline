import unittest
import pandas as pd
from src.python.plot_insert_size_hist import get_hist_vals

class TestGetHistVals(unittest.TestCase):

    def test_get_hist_vals(self):
        # Create a sample histogram file
        with open("test_histogram.txt", "w") as f:
            f.write("# A sample Picard CollectInsertSizeMetrics histogram file\n")
            f.write("# This file is used for testing the plot_insert_size_hist.py script.\n")
            f.write("\n")
            f.write("insert_size\tAll_Reads.fr_count\n")
            f.write("10\t100\n")
            f.write("20\t200\n")
            f.write("30\t300\n")

        # Call the function with the test file
        df = get_hist_vals("test_histogram.txt")

        # Check the output dataframe
        expected_df = pd.DataFrame({
            "insert_size": [10, 20, 30],
            "count": [100, 200, 300]
        })
        pd.testing.assert_frame_equal(df, expected_df)

if __name__ == '__main__':
    unittest.main()
