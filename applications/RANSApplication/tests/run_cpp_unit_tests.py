from KratosMultiphysics import *
from KratosMultiphysics.RANSApplication import *


def run():
    Tester.SetVerbosity(Tester.Verbosity.FAILED_TESTS_OUTPUTS)  # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosRansFastSuite")


if __name__ == '__main__':
    run()
