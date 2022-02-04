from ._core import ffi, lib


__all__ = ('GIT_REVISION', 'TAUOLA_VERSION', 'VERSION')


VERSION = ffi.string(lib.alouette_version()).decode()
TAUOLA_VERSION = ffi.string(lib.tauola_version.version).strip().decode()
GIT_REVISION = ffi.string(lib.alouette_git_revision()).decode()
