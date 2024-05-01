#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

import sys
import json
from functools import wraps
from pathlib import Path

__module = sys.modules[__name__]
__module.cache_dir = Path.cwd() / ".atomphys"


def disk_cache(func):
    @wraps(func)
    def wrapper(arg, refresh_cache=False):
        __module.cache_dir.mkdir(parents=True, exist_ok=True)
        filename = __module.cache_dir / f"{func.__name__}({arg}).cache"
        if not refresh_cache and filename.exists():
            with open(filename) as fp:
                return json.load(fp)
        data = func(arg, refresh_cache)
        with open(filename, "w+") as fp:
            json.dump(data, fp)
        return data

    return wrapper


def get_cache_dir():
    return __module.cache_dir


def set_cache_dir(path):
    __module.cache_dir = Path(path).resolve()


def clear_cache_dir():
    for _file in __module.cache_dir.iterdir():
        _file.unlink()
    __module.cache_dir.rmdir()
