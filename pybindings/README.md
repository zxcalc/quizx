This folder contains Python bindings for quizx so that it can be used with PyZX.
To install, run the install_package.bat file (on Windows).
Or you can install manually by calling::
    
    python setup.py install
    python setup2.py install

The first installs a user-friendly Python interface over the raw Rust bindings. The second installs the actual quizx Python-Rust interface. For some reason the order in which you install these is important (you get a `ImportError: dynamic module does not define module export function` otherwise).
