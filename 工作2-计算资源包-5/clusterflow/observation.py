"""Read runtime observations without letting damaged records hide the audit.

Never use this reader for inputs, authorization, queues or scientific acceptance.
Unreadable observations remain unknown and are explicitly reported.
"""
import json
import time
from pathlib import Path
from .common import readj


def observe(path, errors, *, missing_ok=False):
    path = Path(path)
    for attempt in range(2):
        try:
            value = readj(path)
            if not isinstance(value, dict):
                raise ValueError('Runtime observation must be an object')
            return value
        except FileNotFoundError as exc:
            if missing_ok:
                return None
            error = exc
        except (OSError, ValueError) as exc:
            error = exc
        if attempt == 0:
            time.sleep(.1)
    try:
        size = path.stat().st_size
    except OSError:
        size = None
    errors.append({'path': str(path), 'status': 'UNREADABLE_OBSERVATION',
                   'bytes': size, 'error': type(error).__name__ + ': ' + str(error)})
    return None
