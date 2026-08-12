import unittest

from storca.adaptive_xtb_extension import extension_targets


class AdaptiveXTBExtensionTests(unittest.TestCase):
    def test_extension_targets_are_real_cumulative_rounds(self):
        self.assertEqual(extension_targets(40, 100, 20), [60, 80, 100])

    def test_nonmultiple_maximum_is_retained_as_final_round(self):
        self.assertEqual(extension_targets(40, 95, 20), [60, 80, 95])

    def test_invalid_limits_fail(self):
        with self.assertRaises(ValueError):
            extension_targets(40, 20, 20)


if __name__ == "__main__":
    unittest.main()
