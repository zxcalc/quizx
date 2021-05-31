from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="libquizx",
    version="0.1",
    rust_extensions=[RustExtension("libquizx", binding=Binding.PyO3)],
    zip_safe=False,
)
