import bandsolver
import nmc_model


def test_imports():
    assert nmc_model.__version__
    assert hasattr(bandsolver, "newton")
