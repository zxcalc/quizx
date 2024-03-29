from setuptools import setup  # type: ignore
from setuptools_rust import Binding, RustExtension  # type: ignore

setup(
    name="libquizx",
    version="0.1",
    rust_extensions=[RustExtension("libquizx", binding=Binding.PyO3)],
    zip_safe=False,
)
