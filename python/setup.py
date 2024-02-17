from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import subprocess
import os

# find python include directory from interpreter path
import sysconfig

class BuildCMakeFirst(build_ext):
    def run(self):
        # Ensure CMake build is done before python extension build
        if not os.path.isfile('lib/nikfemm/build/libnikfemm.so'):
            self.build_cmake()
        super().run()

    def build_cmake(self):
        # Build CMake project
        build_dir = os.path.abspath('../lib/nikfemm/build')
        os.makedirs(build_dir, exist_ok=True)
        debug_flag = '-DCMAKE_BUILD_TYPE=Debug' if self.debug else '-DCMAKE_BUILD_TYPE=Release'
        # specific build commands for windows and linux
        if os.name == 'nt':
            subprocess.check_call(['cmake', debug_flag, '-DNIKFEMM_USE_OPENCV=OFF', '-DNIKFEMM_BUILD_TESTS=OFF', '-G', 'MinGW Makefiles', '..'], cwd=build_dir)
            subprocess.check_call(['mingw32-make'], cwd=build_dir)
        else:
            subprocess.check_call(['cmake', debug_flag, '-DNIKFEMM_USE_OPENCV=OFF', '-DNIKFEMM_BUILD_TESTS=OFF', '..'], cwd=build_dir)
            # adaptively use the number of cores
            subprocess.check_call(['make', '-j'], cwd=build_dir)

# Define Python extension
nikfemm_module = Extension(
    'nikfemm',
    sources=['nikfemm.cpp'],
    include_dirs=['../lib/nikfemm/include', sysconfig.get_paths()['include']],
    extra_objects=['../lib/nikfemm/build/libnikfemm_static.a'],  # static library
    extra_link_args=['-Wl,-Bstatic', '-Wl,-Bdynamic'],
)

# link libnikfemm.so to the extension
setup(
    name='nikfemm',
    version='0.1',
    description='Python bindings for nikfemm',
    ext_modules=[nikfemm_module],
    cmdclass={'build_ext': BuildCMakeFirst},
    include_package_data=True,
)