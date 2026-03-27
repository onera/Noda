"""
Get correct work directory when running pytest from package root.

With this fixture, test functions will be run with the directory where they are
located as current work directory, allowing Noda to find the input file (Noda
looks for input files based on Path.cwd()).
Fixture obtained from https://stackoverflow.com/a/62055409
"""
import pytest

@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


"""
Explicitly define non-interactive matplotlib backend.

Needed for tests involving matplotlib to pass when called by pytest on Windows
(otherwise tkinter is defined as backend and it is not properly installed).
"""
import matplotlib as mpl

mpl.rcParams['backend'] = 'Agg'
