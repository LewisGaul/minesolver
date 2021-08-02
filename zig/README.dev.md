## Python wrapper

The Zig executable is wrapped by a minimal Python API in the `zig_minesolver` package. The wrapper assumes the executable is available at the root of the Python package under the name `zig-main` (or `zig-main.exe` on Windows).

### Packaging the Zig code

The `setup.py` script contains an attempt at compiling the Zig code, but the attempt has been abandoned for now because installing the `ziglang` dependency was slow, and it relies on using `distutils` which is in the process of being deprecated.

Instead, the following manual steps can be followed to create wheels and upload to PyPI, which then allows `pip install zig_minesolver` on supported platforms.

- Run the Zig build (from the `zig/` directory):  
  `zig build -Drelease-safe -target x86_64-windows`
- Copy the executable to the Python package:  
  `cp zig-out/bin/zig-main.exe zig_minesolver/`
- Build the Python package:  
  `pip install build && python -m build --wheel`
- Manually change the wheel tag to reflect the supported platform:  
  `mv dist/zig_minesolver-0.1.0-py3-none-any.whl dist/zig_minesolver-0.1.0-py3-none-win_amd64.whl`
- Repeat for other targeted platforms
- Upload all wheels to PyPI:  
  `pip install twine && twine upload dist/*`
