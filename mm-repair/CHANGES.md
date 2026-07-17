# Local changes to mm-repair

This directory contains a redistribution of the third-party `mm-repair` project
(Apache License 2.0). The original `LICENSE.md` and `Readme.md` are retained.

As required by Apache-2.0 section 4(b), the files changed with respect to the
upstream project are listed below. All changes are build/portability fixes; none
alter the compression algorithms or affect the numbers reported in the
manuscript.

| File | Change |
|---|---|
| `brepair/makefile` | Drop the x86-only `-m64` flag so the build works on aarch64. |
| `brepair/basics.h`, `makefile` | Build fixes for current toolchains. |
| `tools/xerrors.h` | Guard CPU-affinity calls and use portable `strerror` on non-Linux platforms. |
| `bin2csrv.cpp`, `brepair/irepair0.c`, `matrepair` | Windows compatibility: open files in binary mode; fall back gracefully when `psutil` is unavailable. |
| `.gitignore` | Added; ignores compiled binaries. |
| `Readme.md` | Citation block removed for double-blind review (see below). |

The upstream `Readme.md` ends with a `### Citation` section giving the reference
for the paper that describes this software. That reference has been replaced with
a placeholder because this repository is under double-blind review. It is cited
in the accompanying manuscript and will be restored in the camera-ready artifact.
