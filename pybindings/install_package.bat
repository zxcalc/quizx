:: For some reason this only works when we first install the Python interface "quizx" and only then the Rust bindings "libquizx"
:: If we do it in the other order we get a "ImportError: dynamic module does not define module export function"

python setup.py install
python setup2.py install
pause