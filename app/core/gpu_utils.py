"""
GPU utility and hardware capability detection.
Detects CUDA availability via CuPy or Numba and provides graceful fallback to CPU.
"""

import logging
from typing import Dict, Tuple

logger = logging.getLogger(__name__)

_GPU_CHECKED = False
_GPU_AVAILABLE = False
_GPU_DEVICE_NAME = None
_GPU_REASON = "Not checked"


def check_gpu_availability() -> Tuple[bool, str, Dict[str, str]]:
    """
    Check if a CUDA-enabled GPU and required acceleration libraries (CuPy) are functional.

    Returns:
        (is_available: bool, device_or_status_message: str, details_dict: dict)
    """
    global _GPU_CHECKED, _GPU_AVAILABLE, _GPU_DEVICE_NAME, _GPU_REASON

    details = {
        "cupy_installed": "no",
        "cuda_runtime": "unavailable",
        "device_count": "0",
        "device_name": "None",
    }

    try:
        import cupy as cp  # type: ignore

        details["cupy_installed"] = "yes"
        count = cp.cuda.runtime.getDeviceCount()
        details["device_count"] = str(count)

        if count > 0:
            dev = cp.cuda.Device(0)
            dev_name = cp.cuda.runtime.getDeviceProperties(0)["name"].decode("utf-8")
            details["cuda_runtime"] = "active"
            details["device_name"] = dev_name

            _GPU_AVAILABLE = True
            _GPU_DEVICE_NAME = dev_name
            _GPU_REASON = f"CUDA GPU available: {dev_name}"
            return True, dev_name, details
        else:
            _GPU_AVAILABLE = False
            _GPU_REASON = "CuPy installed, but 0 CUDA devices found."
            return False, _GPU_REASON, details

    except ImportError:
        _GPU_AVAILABLE = False
        _GPU_REASON = "CuPy not installed; running on optimized CPU mode (NumPy)."
        return False, _GPU_REASON, details
    except Exception as e:
        _GPU_AVAILABLE = False
        _GPU_REASON = f"GPU detection failed ({type(e).__name__}): {str(e)}"
        return False, _GPU_REASON, details


def is_gpu_ready() -> bool:
    """Quick boolean check for GPU readiness."""
    global _GPU_CHECKED, _GPU_AVAILABLE
    if not _GPU_CHECKED:
        is_avail, _, _ = check_gpu_availability()
        _GPU_CHECKED = True
        _GPU_AVAILABLE = is_avail
    return _GPU_AVAILABLE


def get_hardware_status() -> Dict[str, object]:
    """Return comprehensive hardware status for health endpoints and diagnostics."""
    avail, msg, details = check_gpu_availability()
    return {
        "cpu_mode": "Vectorized NumPy SIMD",
        "gpu_available": avail,
        "gpu_device_name": details.get("device_name") if avail else None,
        "status_detail": msg,
        "details": details,
    }

