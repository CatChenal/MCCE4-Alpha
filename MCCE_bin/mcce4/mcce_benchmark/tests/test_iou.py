#!/usr/bin/env python

from mcce4.mcce_benchmark.io_utils import Pathok
from mcce4 import Path
import pytest


class TestPathok:

    def test_return_path_if_exists(self):
        path = "./"
        check_fn = "exists"
        result = Pathok(path, check_fn)

        assert isinstance(result, Path)
        assert result == Path(path)
        assert result.exists()

    def test_return_path_if_directory(self):
        path = "./"
        check_fn = "is_dir"
        result = Pathok(path, check_fn)

        assert isinstance(result, Path)
        assert result == Path(path)
        assert result.is_dir()

    def test_return_path_if_file(self):
        Path("./test_file.txt").touch()
        path = Path("./test_file.txt")
        check_fn = "is_file"
        result = Pathok(path, check_fn)

        assert isinstance(result, Path)
        assert result == Path(path)
        assert result.is_file()
        path.unlink()

    def test_return_path_if_check_fn_invalid_when_path_valid(self):
        path = "./"
        check_fn = "invalid"
        result = Pathok(path, check_fn)

        assert isinstance(result, Path)

    def test_raise_type_error_if_path_invalid(self):
        for path in [123, None, True]:
            with pytest.raises(TypeError):
                Pathok(path)

    def test_raise_file_not_found_error(self):
        with pytest.raises(FileNotFoundError):
            Pathok("nonexistent_path")
